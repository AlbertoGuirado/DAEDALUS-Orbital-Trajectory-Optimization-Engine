"""Evaluated constraints and feasibility reports (RF-4.07, RF-4.08).

Sign convention, used everywhere: a constraint is written g <= 0. The
residual g is in the unit of the constrained quantity and is positive by the
amount the limit is exceeded. A constraint is *satisfied* iff
g <= tolerance, with an explicit absolute tolerance in the same unit; a
non-finite g is never satisfied.

    violation             = max(g, 0) if not satisfied, else 0     [unit]
    normalized violation  = violation / scale                       [-]
    margin                = -g  (positive: slack left)              [unit]

``scale`` is a positive reference magnitude of the quantity (normally its
limit); it only makes violations of different units commensurable so that
they can be added. The raw residual is always kept: a report never reduces a
violation to a boolean, and it never mixes it into an objective value.

Hard constraints decide feasibility. Soft constraints are evaluated and
reported in the same way but do not make a candidate infeasible.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any


@dataclass(frozen=True)
class ConstraintValue:
    """One constraint evaluated on one candidate.

    name: identifier. description: the inequality in words. unit: unit of
    ``value``, ``limit``, ``g`` and ``tolerance``. value: the constrained
    quantity. limit: its bound. g: residual, feasible side g <= 0. scale > 0
    and tolerance >= 0 in ``unit``. hard: whether it decides feasibility.
    Invalid scale or tolerance raises ValueError.
    """

    name: str
    description: str
    unit: str
    value: float
    limit: float
    g: float
    scale: float
    tolerance: float
    hard: bool = True

    def __post_init__(self) -> None:
        for text in ("name", "description", "unit"):
            if not isinstance(getattr(self, text), str) or not getattr(self, text):
                raise ValueError(f"{text} must be a nonempty string")
        for number in ("value", "limit", "g", "scale", "tolerance"):
            value = getattr(self, number)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{number} must be a real number")
            object.__setattr__(self, number, float(value))
        if not (math.isfinite(self.scale) and self.scale > 0):
            raise ValueError("scale must be positive and finite")
        if not (math.isfinite(self.tolerance) and self.tolerance >= 0):
            raise ValueError("tolerance must be nonnegative and finite")
        if not isinstance(self.hard, bool):
            raise ValueError("hard must be a bool")

    @property
    def satisfied(self) -> bool:
        """g <= tolerance; False for a NaN or +inf residual."""
        return bool(self.g <= self.tolerance)

    @property
    def violation(self) -> float:
        """Amount by which the limit is exceeded, in ``unit``; 0 when satisfied, +inf for a NaN residual."""
        if self.satisfied:
            return 0.0
        return math.inf if math.isnan(self.g) else max(self.g, 0.0)

    @property
    def normalized_violation(self) -> float:
        """violation / scale, dimensionless."""
        return self.violation / self.scale

    @property
    def margin(self) -> float:
        """-g, in ``unit``: positive when there is slack."""
        return -self.g

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name, "description": self.description, "unit": self.unit, "hard": self.hard,
            "value": self.value, "limit": self.limit, "g": self.g, "margin": self.margin,
            "tolerance": self.tolerance, "scale": self.scale, "satisfied": self.satisfied,
            "violation": self.violation, "normalized_violation": self.normalized_violation,
        }


@dataclass(frozen=True)
class ConstraintReport:
    """All constraints of one candidate.

    values: the evaluated constraints, in a fixed order. not_applicable:
    name -> reason, for constraints that the physical model cannot evaluate
    (they are reported, never silently assumed satisfied).
    """

    values: tuple[ConstraintValue, ...]
    not_applicable: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        values = tuple(self.values)
        if not all(isinstance(value, ConstraintValue) for value in values):
            raise ValueError("values must contain ConstraintValue objects")
        names = [value.name for value in values] + list(self.not_applicable)
        if len(set(names)) != len(names):
            raise ValueError("constraint names must be unique")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "not_applicable", MappingProxyType(dict(self.not_applicable)))

    @property
    def feasible(self) -> bool:
        """True iff every hard constraint is satisfied."""
        return all(value.satisfied for value in self.values if value.hard)

    @property
    def total_violation(self) -> float:
        """Sum of the normalized violations of the hard constraints; 0 iff feasible."""
        return float(sum(value.normalized_violation for value in self.values if value.hard))

    @property
    def soft_violation(self) -> float:
        """Sum of the normalized violations of the soft constraints."""
        return float(sum(value.normalized_violation for value in self.values if not value.hard))

    @property
    def violated(self) -> tuple[str, ...]:
        """Names of the hard constraints that are not satisfied."""
        return tuple(value.name for value in self.values if value.hard and not value.satisfied)

    def by_name(self, name: str) -> ConstraintValue:
        for value in self.values:
            if value.name == name:
                return value
        raise KeyError(name)

    def normalized(self) -> dict[str, float]:
        """name -> normalized violation of every hard constraint."""
        return {value.name: value.normalized_violation for value in self.values if value.hard}

    def to_dict(self) -> dict[str, Any]:
        return {
            "feasible": self.feasible, "total_violation": self.total_violation,
            "soft_violation": self.soft_violation, "violated": list(self.violated),
            "constraints": [value.to_dict() for value in self.values],
            "not_applicable": dict(self.not_applicable),
        }
