"""Mission limits and their evaluation on a transfer (RF-4.01 to RF-4.06).

Catalogue. Every constraint is g <= 0 (see constraints.model); F is the
thrust level at which the burns are executed, m_wet = m_dry + m_prop,
v_e = Isp g0 and R the reference radius of the central body.

    name          g                                unit  scale        source of the limit
    max_wet_mass  m_wet - m_max                    kg    m_max        mission
    max_thrust    F - min(F_engine, F_max)         N     that limit   vehicle and mission
    max_delta_v   sum_k |dv_k| - dv_max            m/s   dv_max       mission
    max_duration  t_f - t_max                      s     t_max        mission
    min_altitude  h_min - (min_t |r(t)| - R)       m     R + h_min    mission
    propellant    m_wet (1 - exp(-sum dv / v_e))   kg    m_wet        vehicle
                  - m_prop

Notes that matter physically:
  - ``min_altitude`` is a path constraint: it uses the minimum radius over
    the whole trajectory, not the radii at the burns (RF-4.05). The summary
    must therefore carry a trajectory minimum (exact for Keplerian arcs:
    orbital.kepler.radius_extrema_on_arc). Altitude is measured above a
    sphere of radius R; the body has no flattening or terrain here.
  - ``propellant`` and ``max_delta_v`` are different limits. The first
    follows from the vehicle through the rocket equation; the second is an
    independent budget. Neither is inferred from the other.
  - ``max_thrust`` needs a thrust level. The impulsive model has none, so
    without one the constraint is reported as *not applicable*, never as
    satisfied. A thrust level is not inferred from a delta-v.
  - ``max_wet_mass`` and ``max_thrust`` do not depend on the trajectory: they
    reject a vehicle configuration, whatever the transfer.

Feasibility tolerance: tolerance_j = relative_tolerance * scale_j, default
1e-9. The evaluated quantities carry roundoff of at most about 1e-10
relative (conditioning of the analytic propagation), so 1e-9 keeps a
candidate that sits exactly on a limit from being rejected by rounding; it
corresponds to 4e-6 m/s on a 3.9 km/s budget and 7 mm on a LEO radius, far
below the accuracy of the physical model. It is a parameter, recorded with
every report.

``duration`` is the instant of the last burn, with t = 0 at the first event
of the plan; time is elapsed seconds, without an epoch.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from ..constants import EARTH_EQUATORIAL_RADIUS
from ..dynamics.two_body import real_array
from ..spacecraft import Spacecraft
from .model import ConstraintReport, ConstraintValue

CONSTRAINT_NAMES: tuple[str, ...] = (
    "max_wet_mass", "max_thrust", "max_delta_v", "max_duration", "min_altitude", "propellant",
)
DEFAULT_RELATIVE_TOLERANCE: float = 1e-9


def _optional_positive(value: Any, name: str) -> float | None:
    if value is None:
        return None
    array = real_array(value, name)
    if array.shape != () or array <= 0:
        raise ValueError(f"{name} must be a positive scalar or None")
    return float(array)


@dataclass(frozen=True)
class MissionLimits:
    """Limits of a mission; None means that the limit is not imposed.

    max_wet_mass_kg, max_thrust_n, max_delta_v_m_s, max_duration_s: positive.
    min_altitude_m >= 0, above a sphere of body_radius_m > 0 (default: the
    equatorial radius of the Earth). relative_tolerance >= 0: feasibility
    tolerance relative to each constraint's scale. soft: names of the
    constraints that are reported but do not decide feasibility (a subset of
    CONSTRAINT_NAMES; empty by default, so every constraint is hard).
    Invalid values raise ValueError.
    """

    max_wet_mass_kg: float | None = None
    max_thrust_n: float | None = None
    max_delta_v_m_s: float | None = None
    max_duration_s: float | None = None
    min_altitude_m: float | None = None
    body_radius_m: float = EARTH_EQUATORIAL_RADIUS
    relative_tolerance: float = DEFAULT_RELATIVE_TOLERANCE
    soft: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        for name in ("max_wet_mass_kg", "max_thrust_n", "max_delta_v_m_s", "max_duration_s"):
            object.__setattr__(self, name, _optional_positive(getattr(self, name), name))
        if self.min_altitude_m is not None:
            altitude = real_array(self.min_altitude_m, "min_altitude_m")
            if altitude.shape != () or altitude < 0:
                raise ValueError("min_altitude_m must be a nonnegative scalar or None")
            object.__setattr__(self, "min_altitude_m", float(altitude))
        object.__setattr__(self, "body_radius_m", _optional_positive(self.body_radius_m, "body_radius_m"))
        if self.body_radius_m is None:
            raise ValueError("body_radius_m must be a positive scalar")
        tolerance = real_array(self.relative_tolerance, "relative_tolerance")
        if tolerance.shape != () or tolerance < 0:
            raise ValueError("relative_tolerance must be a nonnegative scalar")
        object.__setattr__(self, "relative_tolerance", float(tolerance))
        if isinstance(self.soft, str):
            raise ValueError("soft must be a collection of constraint names, not a string")
        soft = frozenset(self.soft)
        unknown = soft - set(CONSTRAINT_NAMES)
        if unknown:
            raise ValueError(f"unknown constraint names in soft: {sorted(unknown)}")
        object.__setattr__(self, "soft", soft)

    def describe(self) -> dict[str, Any]:
        return {
            "max_wet_mass_kg": self.max_wet_mass_kg, "max_thrust_n": self.max_thrust_n,
            "max_delta_v_m_s": self.max_delta_v_m_s, "max_duration_s": self.max_duration_s,
            "min_altitude_m": self.min_altitude_m, "body_radius_m": self.body_radius_m,
            "relative_tolerance": self.relative_tolerance, "soft": sorted(self.soft),
            "convention": "g <= 0 is feasible; satisfied iff g <= relative_tolerance * scale",
        }


def limits_from_mapping(data: Any) -> MissionLimits:
    """MissionLimits from a mapping of its field names; unknown keys raise ValueError."""
    if not isinstance(data, Mapping):
        raise ValueError("limits must be a mapping")
    known = {"max_wet_mass_kg", "max_thrust_n", "max_delta_v_m_s", "max_duration_s", "min_altitude_m",
             "body_radius_m", "relative_tolerance", "soft"}
    unknown = set(data) - known
    if unknown:
        raise ValueError(f"unknown limit keys {sorted(unknown)}")
    values = dict(data)
    if "soft" in values:
        if isinstance(values["soft"], str):
            raise ValueError("soft must be a list of constraint names")
        values["soft"] = frozenset(values["soft"])
    return MissionLimits(**values)


@dataclass(frozen=True)
class TransferSummary:
    """The quantities of a transfer that the mission constraints depend on.

    burn_delta_v_m_s: magnitudes of the burns in order, m/s, each >= 0.
    duration_s >= 0: instant of the last burn. minimum_radius_m > 0: minimum
    of |r| over the whole trajectory. thrust_n: thrust level of the burns in
    newtons, or None when the model defines none (impulsive burns).
    Non-finite values are allowed in duration and radius only as NaN/inf
    markers of a failed computation, which make the affected constraints
    unsatisfied; negative or non-numeric values raise ValueError.
    """

    burn_delta_v_m_s: tuple[float, ...]
    duration_s: float
    minimum_radius_m: float
    thrust_n: float | None = None

    def __post_init__(self) -> None:
        try:
            burns = tuple(float(value) for value in np.asarray(self.burn_delta_v_m_s, dtype=float).reshape(-1))
        except (TypeError, ValueError) as exc:
            raise ValueError("burn_delta_v_m_s must contain real numbers") from exc
        if any(value < 0 for value in burns):
            raise ValueError("burn magnitudes must be nonnegative")
        object.__setattr__(self, "burn_delta_v_m_s", burns)
        for name in ("duration_s", "minimum_radius_m"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
                raise ValueError(f"{name} must be a real number")
            object.__setattr__(self, name, float(value))
        if self.duration_s < 0 or self.minimum_radius_m <= 0:
            raise ValueError("duration_s must be nonnegative and minimum_radius_m positive")
        object.__setattr__(self, "thrust_n", _optional_positive(self.thrust_n, "thrust_n"))

    @property
    def total_delta_v_m_s(self) -> float:
        """sum_k |dv_k|, m/s."""
        return float(sum(self.burn_delta_v_m_s))

    @classmethod
    def from_metadata(cls, metadata: Mapping[str, Any], thrust_n: float | None = None) -> "TransferSummary":
        """Build from the metadata of a transfer evaluation.

        Required keys: ``burn_delta_v_m_s``, ``transfer_time_s`` and
        ``minimum_radius_m`` (problems.coplanar_multi_impulse). A missing key
        raises ValueError: a constraint is never evaluated on a guess.
        """
        missing = [key for key in ("burn_delta_v_m_s", "transfer_time_s", "minimum_radius_m") if key not in metadata]
        if missing:
            raise ValueError(f"transfer metadata lacks {missing}; the constraints cannot be evaluated")
        return cls(tuple(np.asarray(metadata["burn_delta_v_m_s"], dtype=float)), float(metadata["transfer_time_s"]),
                   float(metadata["minimum_radius_m"]), thrust_n)


def evaluate_constraints(summary: TransferSummary, spacecraft: Spacecraft, limits: MissionLimits) -> ConstraintReport:
    """Evaluate every applicable mission constraint on a transfer (see the module catalogue).

    The propellant constraint is always evaluated; the others when their
    limit is set. ``max_thrust`` is always reported: evaluated when the
    summary has a thrust level, not applicable otherwise.
    """
    if not isinstance(summary, TransferSummary) or not isinstance(spacecraft, Spacecraft) \
            or not isinstance(limits, MissionLimits):
        raise ValueError("evaluate_constraints needs a TransferSummary, a Spacecraft and MissionLimits")
    relative = limits.relative_tolerance
    values: list[ConstraintValue] = []
    not_applicable: dict[str, str] = {}

    def add(name: str, description: str, unit: str, value: float, limit: float, g: float, scale: float) -> None:
        values.append(ConstraintValue(name, description, unit, value, limit, g, scale, relative * scale,
                                      hard=name not in limits.soft))

    if limits.max_wet_mass_kg is not None:
        add("max_wet_mass", "wet mass <= maximum mass", "kg", spacecraft.wet_mass_kg, limits.max_wet_mass_kg,
            spacecraft.wet_mass_kg - limits.max_wet_mass_kg, limits.max_wet_mass_kg)
    thrust_limit = spacecraft.max_thrust_n if limits.max_thrust_n is None \
        else min(spacecraft.max_thrust_n, limits.max_thrust_n)
    if summary.thrust_n is None:
        not_applicable["max_thrust"] = ("the transfer defines no thrust level (impulsive burns); "
                                        "the constraint needs a commanded thrust")
    else:
        add("max_thrust", "burn thrust <= min(engine maximum, mission maximum)", "N", summary.thrust_n,
            thrust_limit, summary.thrust_n - thrust_limit, thrust_limit)
    if limits.max_delta_v_m_s is not None:
        add("max_delta_v", "sum of burn magnitudes <= delta-v budget", "m/s", summary.total_delta_v_m_s,
            limits.max_delta_v_m_s, summary.total_delta_v_m_s - limits.max_delta_v_m_s, limits.max_delta_v_m_s)
    if limits.max_duration_s is not None:
        add("max_duration", "instant of the last burn <= maximum duration", "s", summary.duration_s,
            limits.max_duration_s, summary.duration_s - limits.max_duration_s, limits.max_duration_s)
    if limits.min_altitude_m is not None:
        altitude = summary.minimum_radius_m - limits.body_radius_m
        add("min_altitude", "minimum altitude over the whole trajectory >= minimum altitude", "m", altitude,
            limits.min_altitude_m, limits.min_altitude_m - altitude, limits.body_radius_m + limits.min_altitude_m)
    with np.errstate(over="ignore", invalid="ignore"):
        required = spacecraft.propellant_required_kg(summary.burn_delta_v_m_s) \
            if math.isfinite(summary.total_delta_v_m_s) else math.nan
    add("propellant", "propellant required by the burns (Tsiolkovsky) <= propellant on board", "kg", required,
        spacecraft.propellant_mass_kg, required - spacecraft.propellant_mass_kg, spacecraft.wet_mass_kg)
    return ConstraintReport(tuple(values), not_applicable)
