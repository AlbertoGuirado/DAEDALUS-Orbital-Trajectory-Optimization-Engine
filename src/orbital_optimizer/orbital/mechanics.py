"""Two-body invariants and impulse costs; see equations in docs/phase1.md."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu, validate_state


def _finite(value: ArrayLike) -> NDArray[np.float64]:
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise RuntimeError("orbital calculation exceeds floating-point range")
    return result


def specific_energy(state: ArrayLike, mu: float = EARTH_MU) -> float:
    """epsilon = v.v/2 - mu/|r|, in J/kg (m^2/s^2); state in m and m/s."""
    y, gm = validate_state(state), validate_mu(mu)
    with np.errstate(over="ignore", invalid="ignore"):
        return float(_finite(np.dot(y[3:], y[3:]) / 2 - gm / np.linalg.norm(y[:3])))


def specific_angular_momentum(state: ArrayLike) -> NDArray[np.float64]:
    """h = r cross v, shape (3,), in m^2/s; state in m and m/s."""
    y = validate_state(state)
    with np.errstate(over="ignore", invalid="ignore"):
        return _finite(np.cross(y[:3], y[3:]))


def eccentricity_vector(state: ArrayLike, mu: float = EARTH_MU) -> NDArray[np.float64]:
    """e = (v cross h)/mu - r/|r|, dimensionless (3,); state in SI."""
    y, gm = validate_state(state), validate_mu(mu)
    with np.errstate(over="ignore", invalid="ignore"):
        return _finite(np.cross(y[3:], specific_angular_momentum(y)) / gm - y[:3] / np.linalg.norm(y[:3]))


def node_vector(state: ArrayLike) -> NDArray[np.float64]:
    """n = z_hat cross h, shape (3,), in m^2/s; zero for equatorial orbits."""
    return np.cross([0., 0., 1.], specific_angular_momentum(state))


def orbital_period(semi_major_axis: float, mu: float = EARTH_MU) -> float:
    """T = 2*pi*sqrt(a^3/mu), seconds; a>0 in metres, mu in m^3/s^2."""
    a, gm = real_array(semi_major_axis, "semi_major_axis"), validate_mu(mu)
    if a.shape != () or a <= 0:
        raise ValueError("semi_major_axis must be a positive scalar in metres")
    with np.errstate(over="ignore", invalid="ignore"):
        return float(_finite(2 * np.pi * a * np.sqrt(a / gm)))


def delta_v(velocity_before: ArrayLike, velocity_after: ArrayLike) -> float:
    """|v_after - v_before| in m/s for two velocities of shape (3,), in m/s."""
    before, after = real_array(velocity_before, "velocity_before"), real_array(velocity_after, "velocity_after")
    if before.shape != (3,) or after.shape != (3,):
        raise ValueError("velocities must have shape (3,)")
    with np.errstate(over="ignore", invalid="ignore"):
        return float(_finite(np.hypot.reduce(after - before)))


def total_delta_v(burn_vectors: ArrayLike) -> float:
    """Sum_k |delta_v_k| in m/s, for (N,3) impulse vectors in m/s; [] costs 0."""
    burns = real_array(burn_vectors, "burn_vectors")
    if burns.shape == (0,):
        return 0.0
    if burns.ndim != 2 or burns.shape[1] != 3:
        raise ValueError("burn_vectors must have shape (N,3)")
    with np.errstate(over="ignore", invalid="ignore"):
        return float(_finite(np.sum(np.hypot.reduce(burns, axis=1))))
