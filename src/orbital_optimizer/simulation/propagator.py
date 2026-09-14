"""Numerical two-body propagation in an inertial frame, using SI units."""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp

from ..dynamics.two_body import EARTH_MU, real_array, two_body_dynamics, validate_mu, validate_state
from ..orbital.elements import state_to_elements


@dataclass(frozen=True)
class PropagationResult:
    """Elapsed seconds and states of shape (N, 6), including both endpoints.

    Arrays are copied and made read-only. ``nfev`` counts derivative evaluations;
    it is not the number of output samples or a wall-clock performance estimate.
    """

    times: NDArray[np.float64]
    states: NDArray[np.float64]
    method: str
    nfev: int
    mu: float

    def __post_init__(self) -> None:
        times = np.array(self.times, dtype=float, copy=True)
        states = np.array(self.states, dtype=float, copy=True)
        if times.ndim != 1 or len(times) < 2:
            raise ValueError("times must be a one-dimensional array with at least two entries")
        if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
            raise ValueError("times must be finite and strictly increasing")
        if states.shape != (len(times), 6) or not np.all(np.isfinite(states)):
            raise ValueError("states must be finite with shape (len(times), 6)")
        times.setflags(write=False)
        states.setflags(write=False)
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "states", states)

    @property
    def final_state(self) -> NDArray[np.float64]:
        """Return an independent copy of the state at the requested duration."""
        return self.states[-1].copy()


def _positive_scalar(value: float, name: str) -> float:
    array = real_array(value, name)
    if array.ndim != 0 or np.iscomplexobj(array):
        raise ValueError(f"{name} must be a positive finite scalar")
    try:
        number = float(array)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive finite scalar") from exc
    if not np.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be a positive finite scalar")
    return number


def _sample_times(values: ArrayLike | None, duration: float) -> NDArray | None:
    if values is None:
        return None
    if np.iscomplexobj(values):
        raise ValueError("sample_times must be real")
    times = real_array(values, "sample_times")
    if times.ndim != 1 or times.size < 1:
        raise ValueError("sample_times must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError("sample_times must be finite and strictly increasing")
    if times[0] < 0.0 or times[-1] > duration:
        raise ValueError("sample_times must lie in [0,duration]")
    return np.unique(np.concatenate(([0.0], times, [duration])))


def _rk4(
    initial_state: NDArray, duration: float, mu: float,
    max_step: float, sample_times: NDArray | None,
) -> PropagationResult:
    state = initial_state.copy()
    time = 0.0
    times = [time]
    states = [state.copy()]
    nfev = 0
    targets = [duration] if sample_times is None else sample_times[1:]
    for target in targets:
        while time < target:
            next_time = min(time + max_step, float(target))
            step = next_time - time
            if step <= 0:
                raise RuntimeError("RK4 step is too small to advance floating-point time")
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                try:
                    k1 = two_body_dynamics(time, state, mu)
                    k2 = two_body_dynamics(time + step / 2, state + step * k1 / 2, mu)
                    k3 = two_body_dynamics(time + step / 2, state + step * k2 / 2, mu)
                    k4 = two_body_dynamics(next_time, state + step * k3, mu)
                    state = state + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                except (ValueError, FloatingPointError) as exc:
                    raise RuntimeError(f"RK4 encountered an invalid state at t={time:g} s") from exc
            nfev += 4
            if not np.all(np.isfinite(state)) or np.linalg.norm(state[:3]) == 0:
                raise RuntimeError(f"RK4 produced an invalid state at t={next_time:g} s")
            time = next_time
            if sample_times is None:
                times.append(time)
                states.append(state.copy())
        if sample_times is not None:
            times.append(time)
            states.append(state.copy())
    return PropagationResult(np.array(times), np.array(states), "RK4", nfev, mu)


def propagate(
    initial_state: ArrayLike,
    duration: float,
    mu: float = EARTH_MU,
    *,
    method: Literal["RK45", "RK4"] = "RK45",
    sample_times: ArrayLike | None = None,
    max_step: float = 60.0,
    rtol: float = 1e-10,
    atol: ArrayLike | None = None,
) -> PropagationResult:
    """Propagate ``[x, y, z, vx, vy, vz]`` for a positive duration in seconds.

    Positions are metres, velocities m/s, and mu m^3/s^2. This autonomous
    point-mass model uses elapsed time from zero and has no body-surface model.
    It must not be used to continue trajectories through a central collision.

    RK45 uses SciPy's adaptive Dormand-Prince 5(4) method. Default absolute
    tolerances are 1e-3 m for position and 1e-6 m/s for velocity; relative
    tolerance is 1e-10. These control local estimates, not global error.
    ``atol`` may be a positive scalar or six positive component tolerances.
    ``max_step`` is a finite positive ceiling in seconds (default 60 s).

    RK4 uses classical fourth-order steps of ``max_step``, shortened to hit
    the end and any requested sample time exactly. It has no error control:
    rtol/atol are validated but only used by RK45. Sample times can therefore
    alter the RK4 integration grid; use sample_times=None for convergence
    studies with prescribed step sizes.

    With sample_times=None, return all accepted step endpoints. Otherwise,
    sample_times must be strictly increasing within [0,duration]; missing
    endpoints are added. RK45 interpolates requested samples on its adaptive grid.
    Bad input raises ValueError; integration failure raises RuntimeError.
    """
    state = validate_state(initial_state)
    mu = validate_mu(mu)
    # The validated Phase 1 domain is bound, nonradial elliptic motion.
    state_to_elements(state, mu)
    duration = _positive_scalar(duration, "duration")
    max_step = _positive_scalar(max_step, "max_step")
    rtol = _positive_scalar(rtol, "rtol")
    # SciPy silently clips smaller values; make that change explicit instead.
    if rtol < 100 * np.finfo(float).eps:
        raise ValueError("rtol must be at least 100 * float64 machine epsilon")
    if atol is None:
        tolerances = np.array([1e-3] * 3 + [1e-6] * 3)
    else:
        if np.iscomplexobj(atol):
            raise ValueError("atol must be real")
        tolerances = real_array(atol, "atol")
    if tolerances.shape not in ((), (6,)):
        raise ValueError("atol must be a scalar or a six-component array")
    if not np.all(np.isfinite(tolerances)) or np.any(tolerances <= 0):
        raise ValueError("atol must contain positive finite values")
    samples = _sample_times(sample_times, duration)
    if method == "RK4":
        return _rk4(state, duration, mu, max_step, samples)
    if method != "RK45":
        raise ValueError("method must be 'RK45' or 'RK4'")

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            solution = solve_ivp(
                two_body_dynamics, (0.0, duration), state, args=(mu,),
                method="RK45", t_eval=samples, max_step=max_step,
                rtol=rtol, atol=tolerances,
            )
    except (ValueError, FloatingPointError) as exc:
        raise RuntimeError("RK45 encountered an invalid state during propagation") from exc
    if not solution.success:
        raise RuntimeError(f"RK45 failed: {solution.message}")
    if not np.all(np.isfinite(solution.y)):
        raise RuntimeError("RK45 returned non-finite states")
    return PropagationResult(solution.t, solution.y.T, "RK45", solution.nfev, mu)
