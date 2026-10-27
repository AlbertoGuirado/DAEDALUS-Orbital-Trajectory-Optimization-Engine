"""Numerical propagation with instantaneous velocity changes, in SI units."""

from collections.abc import Sequence
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu, validate_state
from .propagator import propagate


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
