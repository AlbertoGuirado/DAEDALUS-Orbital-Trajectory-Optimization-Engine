"""Numerical propagation with instantaneous velocity changes, in SI units."""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu, validate_state
from .kepler_propagator import propagate_kepler
from .propagator import PropagationResult, propagate


def _burn(time: ArrayLike, delta_v: ArrayLike) -> tuple[float, NDArray[np.float64]]:
    burn_time = real_array(time, "burn time")
    if burn_time.shape != ():
        raise ValueError("burn time must be a scalar in seconds")
    vector = real_array(delta_v, "delta_v")
    if vector.shape != (3,):
        raise ValueError("delta_v must have shape (3,) in m/s")
    return float(burn_time), vector


def propagate_with_impulses(
    initial_state: ArrayLike,
    impulses: Sequence[tuple[float, ArrayLike]],
    duration: float,
    mu: float = EARTH_MU,
    *,
    method: Literal["RK45", "RK4"] = "RK45",
    max_step: float = 60.0,
    rtol: float = 1e-10,
    atol: ArrayLike | None = None,
) -> NDArray[np.float64]:
    """Return the (6,) state at ``duration`` after instantaneous burns.

    ``impulses`` is a sequence of ``(time_s, delta_v_m_s)`` with times strictly
    increasing in [0,duration]; a burn at ``duration`` is applied before return.
    Each burn adds delta_v to the velocity; the position is unchanged. Between
    burns the state is propagated with ``propagate`` using the same options.
    """
    gm = validate_mu(mu)
    duration_s = real_array(duration, "duration")
    if duration_s.shape != () or duration_s <= 0:
        raise ValueError("duration must be a positive scalar in seconds")
    duration_s = float(duration_s)
    state = validate_state(initial_state)
    clock, last = 0.0, -np.inf
    options = {"method": method, "max_step": max_step, "rtol": rtol, "atol": atol}
    for time, delta_v in impulses:
        burn_time, vector = _burn(time, delta_v)
        if not (last < burn_time and 0.0 <= burn_time <= duration_s):
            raise ValueError("burn times must be strictly increasing within [0,duration]")
        if burn_time > clock:
            state = propagate(state, burn_time - clock, gm, **options).final_state
            clock = burn_time
        state[3:] = state[3:] + vector
        last = burn_time
    if duration_s > clock:
        state = propagate(state, duration_s - clock, gm, **options).final_state
    return state


@dataclass(frozen=True, eq=False)
class ImpulsiveTrajectory:
    """Sampled trajectory of a flight with instantaneous burns, SI units.

    arcs: one PropagationResult per free-flight arc, in order; the times of
    each arc are elapsed seconds from the start of that arc. arc_start_times:
    (K,) absolute start time of every arc, s. burn_times (N,) s and
    burn_vectors (N, 3) m/s: the burns as applied. final_state: (6,) state at
    ``duration``, after a burn applied at that instant. A burn at t = 0 or at
    ``duration`` creates no zero-length arc: arcs only cover intervals of
    positive duration. Arrays are read-only.
    """

    arcs: tuple[PropagationResult, ...]
    arc_start_times: NDArray[np.float64]
    burn_times: NDArray[np.float64]
    burn_vectors: NDArray[np.float64]
    final_state: NDArray[np.float64]
    duration: float
    mu: float

    def __post_init__(self) -> None:
        for name in ("arc_start_times", "burn_times", "burn_vectors", "final_state"):
            array = np.array(getattr(self, name), dtype=float, copy=True)
            array.setflags(write=False)
            object.__setattr__(self, name, array)
        object.__setattr__(self, "arcs", tuple(self.arcs))

    def times(self) -> NDArray[np.float64]:
        """Absolute sample times of all arcs, concatenated; arc junctions appear twice."""
        return np.concatenate([start + arc.times for start, arc in zip(self.arc_start_times, self.arcs)])

    def states(self) -> NDArray[np.float64]:
        """Sampled states (M, 6) of all arcs, concatenated; before and after a burn at junctions."""
        return np.concatenate([arc.states for arc in self.arcs])

    def sampled_radius_extrema(self) -> tuple[float, float]:
        """(min, max) of |r| over the samples, m. Sampled values, not continuous-time bounds."""
        radii = np.linalg.norm(self.states()[:, :3], axis=1)
        return float(radii.min()), float(radii.max())


def propagate_impulsive_arcs(
    initial_state: ArrayLike,
    impulses: Sequence[tuple[float, ArrayLike]],
    duration: float,
    mu: float = EARTH_MU,
    *,
    method: Literal["RK45", "RK4", "Kepler"] = "RK45",
    samples_per_arc: int | None = None,
    max_step: float = 60.0,
    rtol: float = 1e-10,
    atol: ArrayLike | None = None,
) -> ImpulsiveTrajectory:
    """Propagate with instantaneous burns and keep the trajectory of every arc (DT-11).

    Same contract as ``propagate_with_impulses`` for ``initial_state``,
    ``impulses`` and ``duration``, and the same final state for the numerical
    methods. ``method``: "RK45" or "RK4" integrate each arc with ``propagate``
    (``max_step``, ``rtol``, ``atol`` as there); "Kepler" uses the analytic
    ``propagate_kepler`` and ignores those options. ``samples_per_arc`` >= 2
    requests that many equally spaced samples on every arc, end points
    included; None keeps the accepted steps of the integrator (for "Kepler",
    only the two end points). Each arc must be a bound, nonradial ellipse.
    Bad input raises ValueError; integration failure raises RuntimeError.
    """
    gm = validate_mu(mu)
    duration_s = real_array(duration, "duration")
    if duration_s.shape != () or duration_s <= 0:
        raise ValueError("duration must be a positive scalar in seconds")
    duration_s = float(duration_s)
    if samples_per_arc is not None and (isinstance(samples_per_arc, bool)
                                        or not isinstance(samples_per_arc, (int, np.integer)) or samples_per_arc < 2):
        raise ValueError("samples_per_arc must be an integer >= 2 or None")
    if method not in ("RK45", "RK4", "Kepler"):
        raise ValueError("method must be 'RK45', 'RK4' or 'Kepler'")
    state = validate_state(initial_state)
    burns: list[tuple[float, NDArray[np.float64]]] = []
    last = -np.inf
    for time, delta_v in impulses:
        burn_time, vector = _burn(time, delta_v)
        if not (last < burn_time and 0.0 <= burn_time <= duration_s):
            raise ValueError("burn times must be strictly increasing within [0,duration]")
        burns.append((burn_time, vector))
        last = burn_time

    def fly(start_state: NDArray[np.float64], span: float) -> PropagationResult:
        samples = None if samples_per_arc is None else np.linspace(0.0, span, samples_per_arc)
        if method == "Kepler":
            return propagate_kepler(start_state, span, gm, sample_times=samples)
        return propagate(start_state, span, gm, method=method, sample_times=samples,
                         max_step=max_step, rtol=rtol, atol=atol)

    arcs: list[PropagationResult] = []
    starts: list[float] = []
    clock = 0.0
    for burn_time, vector in burns:
        if burn_time > clock:
            arc = fly(state, burn_time - clock)
            arcs.append(arc)
            starts.append(clock)
            state = arc.final_state
            clock = burn_time
        state[3:] = state[3:] + vector
    if duration_s > clock:
        arc = fly(state, duration_s - clock)
        arcs.append(arc)
        starts.append(clock)
        state = arc.final_state
    return ImpulsiveTrajectory(
        arcs=tuple(arcs), arc_start_times=np.array(starts),
        burn_times=np.array([time for time, _ in burns]),
        burn_vectors=np.array([vector for _, vector in burns]).reshape(len(burns), 3),
        final_state=state, duration=duration_s, mu=gm,
    )
