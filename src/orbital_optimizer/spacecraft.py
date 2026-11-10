"""Spacecraft mass and propulsion properties (RD-03), SI units.

A data model, not a simulation: it holds the dry mass, the propellant mass,
the maximum thrust and the specific impulse of one engine, and evaluates the
ideal rocket equation with them. With effective exhaust velocity
v_e = Isp * g0 (g0 = constants.STANDARD_GRAVITY, the defining constant of
specific impulse, not local gravity):

    wet mass                 m0 = m_dry + m_prop
    ideal delta-v capacity   dv_cap = v_e ln(m0 / m_dry)
    propellant of a burn     m0 (1 - exp(-dv / v_e))      (Tsiolkovsky)
    burn time at thrust F    t_b = m0 v_e (1 - exp(-dv / v_e)) / F

The burn time follows from the constant mass flow F / v_e. It is a property
of an ideal constant-thrust burn; the impulsive model does not use it, and no
gravity loss is included. A single specific impulse is assumed for every burn.

Distinguish the three masses: ``dry_mass_kg`` excludes propellant,
``propellant_mass_kg`` is what the tanks hold, and ``wet_mass_kg`` is their
sum at the start of the mission.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .constants import STANDARD_GRAVITY
from .dynamics.two_body import real_array


def _scalar(value: Any, name: str) -> float:
    array = real_array(value, name)
    if array.shape != ():
        raise ValueError(f"{name} must be a scalar")
    return float(array)


@dataclass(frozen=True)
class Spacecraft:
    """Mass and propulsion of a vehicle with one engine.

    dry_mass_kg > 0; propellant_mass_kg >= 0; max_thrust_n > 0, the largest
    thrust the engine delivers, in newtons; specific_impulse_s > 0, in
    seconds. Physically impossible values raise ValueError (RD-05).
    """

    dry_mass_kg: float
    propellant_mass_kg: float
    max_thrust_n: float
    specific_impulse_s: float

    def __post_init__(self) -> None:
        for name in ("dry_mass_kg", "propellant_mass_kg", "max_thrust_n", "specific_impulse_s"):
            object.__setattr__(self, name, _scalar(getattr(self, name), name))
        if self.dry_mass_kg <= 0:
            raise ValueError("dry_mass_kg must be positive")
        if self.propellant_mass_kg < 0:
            raise ValueError("propellant_mass_kg must be nonnegative")
        if self.max_thrust_n <= 0:
            raise ValueError("max_thrust_n must be positive, in newtons")
        if self.specific_impulse_s <= 0:
            raise ValueError("specific_impulse_s must be positive, in seconds")
        if not np.isfinite(self.dry_mass_kg + self.propellant_mass_kg):
            raise ValueError("the wet mass must be finite")

    @property
    def wet_mass_kg(self) -> float:
        """Dry mass plus propellant, kg."""
        return self.dry_mass_kg + self.propellant_mass_kg

    @property
    def exhaust_velocity_m_s(self) -> float:
        """Effective exhaust velocity v_e = Isp * g0, m/s."""
        return self.specific_impulse_s * STANDARD_GRAVITY

    @property
    def delta_v_capacity_m_s(self) -> float:
        """Ideal delta-v of the whole propellant load, v_e ln(m_wet / m_dry), m/s."""
        return float(self.exhaust_velocity_m_s * np.log(self.wet_mass_kg / self.dry_mass_kg))

    def propellant_per_burn_kg(self, burn_delta_v_m_s: ArrayLike) -> NDArray[np.float64]:
        """Propellant consumed by each burn of a sequence, kg, shape (N,).

        Burn k starts from the mass left by burns 1..k-1:
        m_k = m_{k-1} exp(-dv_k / v_e), consumption m_{k-1} - m_k, starting
        from the wet mass. The total equals m0 (1 - exp(-sum(dv) / v_e)). The
        propellant actually on board is not checked here: a sequence may ask
        for more than the tanks hold (see constraints). ``burn_delta_v_m_s``
        holds nonnegative magnitudes in m/s, shape (N,).
        """
        burns = real_array(burn_delta_v_m_s, "burn_delta_v_m_s")
        if burns.ndim != 1 or np.any(burns < 0):
            raise ValueError("burn_delta_v_m_s must be a one-dimensional array of nonnegative magnitudes in m/s")
        cumulative = np.concatenate(([0.0], np.cumsum(burns)))
        masses = self.wet_mass_kg * np.exp(-cumulative / self.exhaust_velocity_m_s)
        return -np.diff(masses)

    def propellant_required_kg(self, burn_delta_v_m_s: ArrayLike) -> float:
        """Total propellant of a burn sequence, m0 (1 - exp(-sum(dv) / v_e)), kg."""
        burns = real_array(burn_delta_v_m_s, "burn_delta_v_m_s")
        if burns.ndim != 1 or np.any(burns < 0):
            raise ValueError("burn_delta_v_m_s must be a one-dimensional array of nonnegative magnitudes in m/s")
        return float(self.wet_mass_kg * -np.expm1(-float(np.sum(burns)) / self.exhaust_velocity_m_s))

    def burn_durations_s(self, burn_delta_v_m_s: ArrayLike, thrust_n: float | None = None) -> NDArray[np.float64]:
        """Duration of each burn at constant thrust, s, shape (N,).

        t_k = (propellant of burn k) * v_e / F, from the constant mass flow
        F / v_e. ``thrust_n`` defaults to ``max_thrust_n`` and must be
        positive. Ideal constant-thrust burns without gravity loss; the value
        tells how far a burn is from instantaneous.
        """
        thrust = self.max_thrust_n if thrust_n is None else _scalar(thrust_n, "thrust_n")
        if thrust <= 0:
            raise ValueError("thrust_n must be positive, in newtons")
        return self.propellant_per_burn_kg(burn_delta_v_m_s) * self.exhaust_velocity_m_s / thrust

    def describe(self) -> dict[str, float]:
        """JSON-compatible record of the properties and the derived quantities."""
        return {
            "dry_mass_kg": self.dry_mass_kg, "propellant_mass_kg": self.propellant_mass_kg,
            "wet_mass_kg": self.wet_mass_kg, "max_thrust_n": self.max_thrust_n,
            "specific_impulse_s": self.specific_impulse_s, "exhaust_velocity_m_s": self.exhaust_velocity_m_s,
            "delta_v_capacity_m_s": self.delta_v_capacity_m_s,
        }


def spacecraft_from_mapping(data: Any) -> Spacecraft:
    """Spacecraft from a mapping with exactly the four RD-03 keys; anything else raises ValueError."""
    keys = {"dry_mass_kg", "propellant_mass_kg", "max_thrust_n", "specific_impulse_s"}
    if not hasattr(data, "keys") or set(data.keys()) != keys:
        raise ValueError(f"spacecraft must contain exactly these keys: {sorted(keys)}")
    return Spacecraft(**{key: data[key] for key in keys})


__all__: Sequence[str] = ("Spacecraft", "spacecraft_from_mapping")
