"""Hohmann transfer between coplanar circular orbits; equations in docs/phase2-api.md."""

from dataclasses import dataclass

import numpy as np

from ..dynamics.two_body import EARTH_MU, real_array, validate_mu
from .mechanics import orbital_period


@dataclass(frozen=True)
class HohmannTransfer:
    """Impulse magnitudes and transfer timing; m/s, s and m.

    Both burns are magnitudes, so the result is the same for an ascending or
    descending transfer.
    """

    first_delta_v: float
    second_delta_v: float
    total_delta_v: float
    transfer_time: float
    transfer_semi_major_axis: float


def _radius(value: float, name: str) -> float:
    radius = real_array(value, name)
    if radius.shape != () or radius <= 0:
        raise ValueError(f"{name} must be a positive scalar in metres")
    return float(radius)


def hohmann_transfer(initial_radius: float, final_radius: float, mu: float = EARTH_MU) -> HohmannTransfer:
    """Two-impulse transfer between circular radii in m; mu in m^3/s^2.

    Velocities come from vis-viva at the circular orbits and at the transfer
    apsides. The transfer time is half the period of the transfer ellipse.
    """
    r1 = _radius(initial_radius, "initial_radius")
    r2 = _radius(final_radius, "final_radius")
    gm = validate_mu(mu)
    a = (r1 + r2) / 2
    with np.errstate(over="ignore", invalid="ignore"):
        v_circular_1 = np.sqrt(gm / r1)
        v_circular_2 = np.sqrt(gm / r2)
        v_perigee = np.sqrt(gm * (2 / r1 - 1 / a))
        v_apogee = np.sqrt(gm * (2 / r2 - 1 / a))
    first = float(abs(v_perigee - v_circular_1))
    second = float(abs(v_circular_2 - v_apogee))
    return HohmannTransfer(
        first_delta_v=first,
        second_delta_v=second,
        total_delta_v=first + second,
        transfer_time=orbital_period(a, gm) / 2,
        transfer_semi_major_axis=float(a),
    )
