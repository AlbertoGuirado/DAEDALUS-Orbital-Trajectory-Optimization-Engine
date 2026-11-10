"""Hohmann and bi-elliptic transfers between coplanar circular orbits; equations in docs/phase2-api.md."""

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


@dataclass(frozen=True)
class BiellipticTransfer:
    """Impulse magnitudes and timing of a bi-elliptic transfer; m/s, s and m.

    The three burns are magnitudes. ``transfer_time`` is the sum of the two
    half periods; ``first_semi_major_axis`` and ``second_semi_major_axis``
    are those of the two transfer ellipses.
    """

    first_delta_v: float
    second_delta_v: float
    third_delta_v: float
    total_delta_v: float
    transfer_time: float
    first_semi_major_axis: float
    second_semi_major_axis: float


def bielliptic_transfer(initial_radius: float, final_radius: float, intermediate_radius: float,
                        mu: float = EARTH_MU) -> BiellipticTransfer:
    """Three-impulse transfer between circular radii through an outer apoapsis (RF-2.05).

    Radii in m; mu in m^3/s^2; rb = intermediate_radius >= max(r1, r2). With
    a1 = (r1 + rb)/2 and a2 = (r2 + rb)/2, vis-viva v(r, a) = sqrt(mu (2/r - 1/a))
    and v_c(r) = sqrt(mu/r):

        dv1 = |v(r1, a1) - v_c(r1)|      tangential burn at r1
        dv2 = |v(rb, a2) - v(rb, a1)|    tangential burn at the apoapsis rb
        dv3 = |v_c(r2) - v(r2, a2)|      tangential burn at r2
        t   = pi sqrt(a1^3/mu) + pi sqrt(a2^3/mu)

    (Vallado, Fundamentals of Astrodynamics and Applications, 4th ed.,
    Sec. 6.3.2; Curtis, Orbital Mechanics for Engineering Students, Sec. 6.4.)
    For rb equal to the larger of r1 and r2 the transfer degenerates into the
    Hohmann transfer with one null burn. Coplanar circular orbits, two-body
    motion, instantaneous impulses. Invalid input raises ValueError; a
    non-finite result raises RuntimeError.
    """
    r1 = _radius(initial_radius, "initial_radius")
    r2 = _radius(final_radius, "final_radius")
    rb = _radius(intermediate_radius, "intermediate_radius")
    gm = validate_mu(mu)
    if rb < max(r1, r2):
        raise ValueError("intermediate_radius must not be below the initial and final radii")
    a1, a2 = (r1 + rb) / 2, (r2 + rb) / 2
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        first = float(abs(np.sqrt(gm * (2 / r1 - 1 / a1)) - np.sqrt(gm / r1)))
        second = float(abs(np.sqrt(gm * (2 / rb - 1 / a2)) - np.sqrt(gm * (2 / rb - 1 / a1))))
        third = float(abs(np.sqrt(gm / r2) - np.sqrt(gm * (2 / r2 - 1 / a2))))
        total = first + second + third
    if not np.isfinite(total):
        raise RuntimeError("bi-elliptic transfer exceeds floating-point range")
    return BiellipticTransfer(
        first_delta_v=first,
        second_delta_v=second,
        third_delta_v=third,
        total_delta_v=total,
        transfer_time=orbital_period(a1, gm) / 2 + orbital_period(a2, gm) / 2,
        first_semi_major_axis=float(a1),
        second_semi_major_axis=float(a2),
    )
