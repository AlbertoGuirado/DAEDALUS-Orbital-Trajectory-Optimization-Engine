"""Termination criteria and immutable optimization results.

Every optimizer shares these structures so that algorithms can be compared on
the same problem with the same evaluation budget (RF-3.10, RF-3.17, RF-3.18).
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .problem import Evaluation, EvaluationStatus

# Fixed order of status codes stored in OptimizationHistory.status (int8).
STATUS_ORDER: tuple[EvaluationStatus, ...] = tuple(EvaluationStatus)


class TerminationReason(Enum):
    """Why a run stopped. Checked in the order of TERMINATION_PRECEDENCE."""

    MAX_EVALUATIONS = "max_evaluations"
    TARGET_REACHED = "target_reached"
    POPULATION_CONVERGED = "population_converged"
    OBJECTIVE_CONVERGED = "objective_converged"
    STAGNATION = "stagnation"
    MAX_GENERATIONS = "max_generations"


# When several criteria hold at the same check, the first listed is reported.
TERMINATION_PRECEDENCE: tuple[TerminationReason, ...] = tuple(TerminationReason)


def _positive_int(value: Any, name: str, *, optional: bool) -> int | None:
    if value is None and optional:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be an integer >= 1" + (" or None" if optional else ""))
    return int(value)


def _optional_float(value: Any, name: str, *, nonnegative: bool) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a real scalar or None")
    number = float(value)
    if not math.isfinite(number) or (nonnegative and number < 0):
        raise ValueError(f"{name} must be finite" + (" and nonnegative" if nonnegative else ""))
    return number


@dataclass(frozen=True)
class TerminationCriteria:
    """Stopping rules shared by all optimizers.

    max_evaluations: hard budget of objective evaluations, counting invalid
        ones; never exceeded (required, >= 1).
    max_generations: iterations after initialization (>= 1), or None.
    target_objective: stop once the best valid objective is <= this value, in
        the objective's unit. Intended for benchmarks with a known optimum.
    x_tol: stop once max_j (max_i u_ij - min_i u_ij) <= x_tol, the largest
        per-variable population range in normalized [0,1] units. Periodic
        variables use the plain range, which is conservative near the seam.
    f_tol: stop once every population member is valid and
        max(f) - min(f) <= f_tol, in the objective's unit.
    stagnation_generations: stop after this many consecutive generations
        without a strict decrease of the best valid objective.

    Optional criteria default to None (disabled): the default run stops only
    when the budget is exhausted. Invalid values raise ValueError.
    """

    max_evaluations: int
    max_generations: int | None = None
    target_objective: float | None = None
    x_tol: float | None = None
    f_tol: float | None = None
    stagnation_generations: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "max_evaluations", _positive_int(self.max_evaluations, "max_evaluations", optional=False))
        object.__setattr__(self, "max_generations", _positive_int(self.max_generations, "max_generations", optional=True))
        object.__setattr__(self, "target_objective", _optional_float(self.target_objective, "target_objective", nonnegative=False))
        object.__setattr__(self, "x_tol", _optional_float(self.x_tol, "x_tol", nonnegative=True))
        object.__setattr__(self, "f_tol", _optional_float(self.f_tol, "f_tol", nonnegative=True))
        object.__setattr__(
            self, "stagnation_generations",
            _positive_int(self.stagnation_generations, "stagnation_generations", optional=True),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "max_evaluations": self.max_evaluations,
            "max_generations": self.max_generations,
            "target_objective": self.target_objective,
            "x_tol": self.x_tol,
            "f_tol": self.f_tol,
            "stagnation_generations": self.stagnation_generations,
        }


def _frozen(values: Any, dtype: Any) -> NDArray:
    array = np.array(values, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True, eq=False)
class OptimizationHistory:
    """Convergence record (RF-3.18). All arrays are read-only copies.

    Per evaluation k = 0..nfev-1:
      best_so_far[k]: best valid objective after evaluation k (+inf before
        the first valid one); status[k]: index into STATUS_ORDER.
    Per generation g = 0..G (g = 0 is the initial population):
      generation_nfev[g]: cumulative evaluations at the end of generation g;
      generation_best[g]: best valid objective at that point;
      generation_valid[g]: valid members in the population;
      generation_spread[g]: largest per-variable normalized range.
    populations: (G+1, NP, D) normalized populations, or None if not recorded.
    """

    best_so_far: NDArray[np.float64]
    status: NDArray[np.int8]
    generation_nfev: NDArray[np.int64]
    generation_best: NDArray[np.float64]
    generation_valid: NDArray[np.int64]
    generation_spread: NDArray[np.float64]
    populations: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "best_so_far", _frozen(self.best_so_far, float))
        object.__setattr__(self, "status", _frozen(self.status, np.int8))
        object.__setattr__(self, "generation_nfev", _frozen(self.generation_nfev, np.int64))
        object.__setattr__(self, "generation_best", _frozen(self.generation_best, float))
        object.__setattr__(self, "generation_valid", _frozen(self.generation_valid, np.int64))
        object.__setattr__(self, "generation_spread", _frozen(self.generation_spread, float))
        if self.populations is not None:
            object.__setattr__(self, "populations", _frozen(self.populations, float))
        if self.best_so_far.shape != self.status.shape or self.best_so_far.ndim != 1:
            raise ValueError("per-evaluation arrays must be one-dimensional with equal length")
        lengths = {len(self.generation_nfev), len(self.generation_best),
                   len(self.generation_valid), len(self.generation_spread)}
        if len(lengths) != 1:
            raise ValueError("per-generation arrays must have equal length")

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "best_so_far": self.best_so_far,
            "status_order": [status.value for status in STATUS_ORDER],
            "status": self.status,
            "generation_nfev": self.generation_nfev,
            "generation_best": self.generation_best,
            "generation_valid": self.generation_valid,
            "generation_spread": self.generation_spread,
        }
        if self.populations is not None:
            data["populations"] = self.populations
        return _jsonable(data)


def _jsonable(value: Any) -> Any:
    """Convert to strict-JSON data; non-finite floats become 'inf'/'-inf'/'nan'."""
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if math.isfinite(number):
            return number
        return "nan" if math.isnan(number) else ("inf" if number > 0 else "-inf")
    if value is None or isinstance(value, str):
        return value
    if hasattr(value, "to_dict"):
        return _jsonable(value.to_dict())
    return repr(value)


def jsonable(value: Any) -> Any:
    """Strict-JSON form of results and reports: arrays to lists, enums to values,

    mappings to dicts with string keys and non-finite floats to the strings
    'inf', '-inf' and 'nan' (strict JSON has no such numbers).
    """
    return _jsonable(value)


@dataclass(frozen=True, eq=False)
class OptimizationResult:
    """Outcome of one optimization run.

    u_best, (D,), is the normalized point the optimizer evaluated, and x_best
    = bounds.from_unit(u_best) is the physical vector that was passed to the
    problem. The map is applied in that direction only: bounds.to_unit(x_best)
    equals u_best only up to floating-point rounding, and for a periodic
    variable it can differ by one period when from_unit wrapped upper to
    lower. Both are None, with f_best = +inf, when no evaluation was valid. ``found_valid`` states only that a valid point was found; it
    makes no claim of optimality. ``status_counts`` maps each status value to
    its number of evaluations. ``wall_time_s`` is informative and excluded from
    reproducibility comparisons.
    """

    x_best: NDArray[np.float64] | None
    u_best: NDArray[np.float64] | None
    f_best: float
    best_evaluation: Evaluation | None
    found_valid: bool
    termination_reason: TerminationReason
    message: str
    nfev: int
    status_counts: Mapping[str, int]
    generations: int
    seed: int
    algorithm: str
    settings: Mapping[str, Any]
    termination: TerminationCriteria
    problem_info: Mapping[str, Any]
    environment: Mapping[str, Any]
    wall_time_s: float
    history: OptimizationHistory

    def __post_init__(self) -> None:
        for name in ("x_best", "u_best"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _frozen(value, float))
        for name in ("status_counts", "settings", "problem_info", "environment"):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))

    @property
    def n_invalid(self) -> int:
        """Number of evaluations whose status is not VALID."""
        return self.nfev - int(self.status_counts.get(EvaluationStatus.VALID.value, 0))

    def to_dict(self, *, include_history: bool = True) -> dict[str, Any]:
        """Strict-JSON-compatible record of the run, configuration and environment."""
        best = None
        if self.best_evaluation is not None:
            best = {
                "x": self.best_evaluation.x,
                "objective": self.best_evaluation.objective,
                "status": self.best_evaluation.status,
                "message": self.best_evaluation.message,
                "metadata": self.best_evaluation.metadata,
            }
        data: dict[str, Any] = {
            "algorithm": self.algorithm,
            "seed": self.seed,
            "settings": self.settings,
            "termination": self.termination,
            "problem": self.problem_info,
            "found_valid": self.found_valid,
            "x_best": self.x_best,
            "u_best": self.u_best,
            "f_best": self.f_best,
            "best_evaluation": best,
            "termination_reason": self.termination_reason,
            "message": self.message,
            "nfev": self.nfev,
            "n_invalid": self.n_invalid,
            "status_counts": self.status_counts,
            "generations": self.generations,
            "environment": self.environment,
            "wall_time_s": self.wall_time_s,
        }
        if include_history:
            data["history"] = self.history
        return _jsonable(data)
