"""Weighted, nondimensional scalarization of several metrics (RF-3.19 to RF-3.21).

    J(x) = sum_i w_i * m_i(x) / m_ref,i

Each term names a metric m_i that the wrapped problem reports in
``Evaluation.metadata``, the positive reference scale m_ref,i in the same unit
and a dimensionless weight w_i >= 0, so every summand is dimensionless and no
unit conversion hides inside a weight (docs/03 section 10.3). Every term must
state its unit and the rationale of its weight and reference; a term without
them is rejected (RF-3.20).

This module is generic: it knows nothing about orbits. The orbital adapter
declares the delta-v, time, energy and risk terms (problems/objectives.py).

A weighted sum can only reach the convex part of a Pareto front (docs/03
section 14.3). It is a way to state a preference, not a substitute for a
multi-objective analysis.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .problem import Bounds, Evaluation, EvaluationStatus, Problem


def _text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _real(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a real scalar")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


@dataclass(frozen=True)
class ObjectiveTerm:
    """One summand w * metric / reference of a weighted objective.

    name: label of the term. metric: key of the value in Evaluation.metadata.
    weight: dimensionless, >= 0. reference: > 0, in ``unit``, the unit of the
    metric. rationale: why this weight and reference were chosen.
    """

    name: str
    metric: str
    weight: float
    reference: float
    unit: str
    rationale: str

    def __post_init__(self) -> None:
        for field_name in ("name", "metric", "unit", "rationale"):
            _text(getattr(self, field_name), field_name)
        weight = _real(self.weight, "weight")
        reference = _real(self.reference, "reference")
        if weight < 0.0:
            raise ValueError("weight must be nonnegative")
        if reference <= 0.0:
            raise ValueError("reference must be positive")
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "reference", reference)

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "metric": self.metric, "weight": self.weight,
                "reference": self.reference, "unit": self.unit, "rationale": self.rationale}


class WeightedObjective:
    """J = sum_i w_i m_i / m_ref,i over uniquely named terms with at least one positive weight."""

    def __init__(self, terms: Sequence[ObjectiveTerm]) -> None:
        terms = tuple(terms)
        if not terms or not all(isinstance(term, ObjectiveTerm) for term in terms):
            raise ValueError("terms must be a nonempty sequence of ObjectiveTerm")
        if len({term.name for term in terms}) != len(terms):
            raise ValueError("term names must be unique")
        if not any(term.weight > 0.0 for term in terms):
            raise ValueError("at least one weight must be positive")
        self.terms = terms

    def contributions(self, metrics: Mapping[str, Any]) -> dict[str, float]:
        """Dimensionless w_i m_i / m_ref,i per term. A missing or non-scalar metric raises ValueError.

        Terms with zero weight are not read, so a metric that a problem does
        not report can still be declared with weight 0.
        """
        result: dict[str, float] = {}
        for term in self.terms:
            if term.weight == 0.0:
                result[term.name] = 0.0
                continue
            if term.metric not in metrics:
                raise ValueError(f"metric {term.metric!r} of term {term.name!r} is not reported by the problem")
            value = np.asarray(metrics[term.metric])
            if value.shape != () or value.dtype.kind not in "iuf":
                raise ValueError(f"metric {term.metric!r} must be a real scalar")
            result[term.name] = term.weight * float(value) / term.reference
        return result

    def value(self, metrics: Mapping[str, Any]) -> float:
        """J, dimensionless; may be non-finite if a metric is."""
        return float(sum(self.contributions(metrics).values()))

    def describe(self) -> dict[str, Any]:
        return {"form": "sum_i weight_i * metric_i / reference_i",
                "terms": [term.describe() for term in self.terms]}


class ScalarizedProblem:
    """Wrap a Problem so that its objective is a WeightedObjective of its metadata.

    Invalid evaluations of the base problem are returned unchanged. For a
    valid one the new objective is J; the base objective and the per-term
    contributions are added to the metadata as ``base_objective`` and
    ``objective_terms``. A non-finite J gives a NON_FINITE evaluation.
    """

    def __init__(self, problem: Problem, objective: WeightedObjective, name: str | None = None) -> None:
        if not isinstance(getattr(problem, "bounds", None), Bounds) or not callable(getattr(problem, "evaluate", None)):
            raise ValueError("problem must implement the Problem protocol")
        if not isinstance(objective, WeightedObjective):
            raise ValueError("objective must be a WeightedObjective")
        self.problem = problem
        self.objective = objective
        self.bounds = problem.bounds
        self.name = name if name is not None else f"{problem.name}:weighted"
        _text(self.name, "name")

    def evaluate(self, x: ArrayLike) -> Evaluation:
        """Evaluate the base problem at x and scalarize its metadata."""
        base = self.problem.evaluate(x)
        if not base.is_valid:
            return base
        contributions = self.objective.contributions(base.metadata)
        value = float(sum(contributions.values()))
        metadata = dict(base.metadata)
        metadata["base_objective"] = base.objective
        metadata["objective_terms"] = contributions
        if not math.isfinite(value):
            return Evaluation(base.x, value, EvaluationStatus.NON_FINITE,
                              f"weighted objective is {value}", metadata)
        return Evaluation(base.x, value, EvaluationStatus.VALID, base.message, metadata)
