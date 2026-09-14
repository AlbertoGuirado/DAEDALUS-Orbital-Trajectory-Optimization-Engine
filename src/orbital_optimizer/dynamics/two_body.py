"""Newtonian point-mass gravity: r'' = -mu r / |r|^3.

Equations and scope: docs/phase1.md; inertial Cartesian SI coordinates.
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

# Conventional Earth GM, IERS TN36, Table 1.1; m^3/s^2.
# https://iers-conventions.obspm.fr/content/tn36.pdf
EARTH_MU: float = 3.986004418e14


def real_array(value: ArrayLike, name: str) -> NDArray[np.float64]:
    """Copy finite real numerical data, preserving its caller-defined SI units."""
    try:
        array = np.asarray(value)
        if array.dtype.kind not in "iuf":
            raise ValueError(f"{name} must contain real numbers")
        result = np.array(array, dtype=float, copy=True)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values")
    return result


def validate_mu(mu: float) -> float:
    """Validate a positive finite scalar gravitational parameter, in m^3/s^2."""
    value = real_array(mu, "mu")
    if value.shape != () or value <= 0:
        raise ValueError("mu must be a positive finite scalar in m^3/s^2")
    return float(value)


def validate_state(state: ArrayLike) -> NDArray[np.float64]:
    """Copy a finite real (6,) state [r in m, v in m/s], with nonzero radius."""
    value = real_array(state, "state")
    if value.shape != (6,):
        raise ValueError("state must have shape (6,): [x,y,z,vx,vy,vz]")
    radius = np.hypot.reduce(value[:3])
    if radius == 0 or not np.isfinite(radius):
        raise ValueError("state must have a finite nonzero radius")
    return value


def two_body_dynamics(t: float, state: ArrayLike, mu: float = EARTH_MU) -> NDArray[np.float64]:
    """Return [v, -mu*r/|r|^3], shape (6,), in m/s and m/s^2.

    ``t`` is elapsed seconds (unused in this autonomous Newtonian equation).
    Input state uses metres and m/s; mu uses m^3/s^2. No surface or collision
    continuation is modeled. See docs/phase1.md, Newton's two-body equation.
    """
    value = validate_state(state)
    parameter = validate_mu(mu)
    radius = np.hypot.reduce(value[:3])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        try:
            acceleration = -(parameter / radius / radius) * (value[:3] / radius)
        except FloatingPointError as exc:
            raise RuntimeError("gravity exceeds floating-point range") from exc
    return np.concatenate((value[3:], acceleration))
