"""Algorithm-independent optimizer interface, seeding and evaluation budget."""

from __future__ import annotations

import copy
import math
import platform
import sys
from typing import Any, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from .problem import Bounds, Evaluation, Problem
from .result import STATUS_ORDER, OptimizationResult, TerminationCriteria

BIT_GENERATOR: str = "PCG64"


@runtime_checkable
class Optimizer(Protocol):
    """Minimize any Problem without modifying it (RF-3.10).

    Implementations must respect ``termination.max_evaluations`` exactly
    (RF-3.17), record convergence history (RF-3.18) and give identical
    results for the same seed, settings, problem and environment (RF-3.16).
    """

    name: str

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult: ...


def validate_seed(seed: Any) -> int:
    """Return seed as int; it must be an integer >= 0 (bool rejected)."""
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be an integer >= 0")
    return int(seed)


def make_rng(seed: int) -> np.random.Generator:
    """Independent generator ``Generator(PCG64(seed))``; no global random state."""
    return np.random.Generator(np.random.PCG64(validate_seed(seed)))


def validate_problem(problem: Any) -> Bounds:
    """Check the Problem protocol structurally and return its bounds."""
    bounds = getattr(problem, "bounds", None)
    if not isinstance(bounds, Bounds):
        raise TypeError("problem.bounds must be a Bounds instance")
    if not callable(getattr(problem, "evaluate", None)):
        raise TypeError("problem must define evaluate(x)")
    if not isinstance(getattr(problem, "name", None), str):
        raise TypeError("problem.name must be a string")
    return bounds


def environment_info() -> dict[str, Any]:
    """Versions and platform needed to interpret reproducibility claims."""
    return {
        "python": sys.version.split()[0],
        "python_implementation": platform.python_implementation(),
        "numpy": np.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "bit_generator": BIT_GENERATOR,
    }


def problem_info(problem: Problem) -> dict[str, Any]:
    """Name, type and bounds description of a problem."""
    return {"name": problem.name, "type": type(problem).__qualname__, "bounds": problem.bounds.describe()}


class EvaluationBudget:
    """Evaluate normalized points under a hard evaluation budget.

    Each call maps u in [0,1]^D to x = lower + u*(upper - lower), validates x,
    calls ``problem.evaluate`` and records status and best-so-far value. The
    best is the first evaluated point with the lowest valid objective.

    The problem receives its own read-only copy of x; the budget keeps a
    separate private copy, so a problem that re-enables writing and mutates
    its argument (now or later) cannot alter ``best_x``. Likewise
    ``best_evaluation`` is not the problem's object but a new Evaluation
    built from the private x with ``copy.deepcopy`` of the metadata.

    Contract checks on the problem's answer:
      - not an Evaluation: TypeError;
      - ``evaluation.x`` without shape (D,), or not exactly equal to the
        evaluated x: ValueError (the result would describe another point).
    Calling ``evaluate`` with the budget exhausted is a programming error
    (RuntimeError). Exceptions raised by the problem propagate. A failed
    contract check raises before anything is recorded.
    """

    def __init__(self, problem: Problem, max_evaluations: int) -> None:
        self.problem = problem
        self.bounds = validate_problem(problem)
        if isinstance(max_evaluations, bool) or not isinstance(max_evaluations, (int, np.integer)) or max_evaluations < 1:
            raise ValueError("max_evaluations must be an integer >= 1")
        self.max_evaluations = int(max_evaluations)
        self._status: list[int] = []
        self._best_so_far: list[float] = []
        self.best_value: float = math.inf
        self.best_x: NDArray[np.float64] | None = None
        self.best_u: NDArray[np.float64] | None = None
        self.best_evaluation: Evaluation | None = None

    @property
    def used(self) -> int:
        return len(self._status)

    @property
    def remaining(self) -> int:
        return self.max_evaluations - self.used

    @property
    def exhausted(self) -> bool:
        return self.used >= self.max_evaluations

    def evaluate(self, u: NDArray[np.float64]) -> Evaluation:
        if self.exhausted:
            raise RuntimeError("evaluation budget exhausted")
        x = self.bounds.validate(self.bounds.from_unit(u))
        given = np.array(x, copy=True)
        given.setflags(write=False)
        evaluation = self.problem.evaluate(given)
        if not isinstance(evaluation, Evaluation):
            raise TypeError("problem.evaluate must return an Evaluation")
        if evaluation.x.shape != x.shape:
            raise ValueError(f"evaluation.x must have shape {x.shape}, got {evaluation.x.shape}")
        if not np.array_equal(evaluation.x, x):
            raise ValueError("evaluation.x differs from the evaluated decision vector")
        self._status.append(STATUS_ORDER.index(evaluation.status))
        if evaluation.is_valid and evaluation.objective < self.best_value:
            self.best_value = evaluation.objective
            self.best_x = x
            self.best_u = np.array(u, dtype=float, copy=True)
            # Budget-owned record: rebuilt from the private x with deep-copied metadata,
            # so later changes to the problem's own Evaluation object cannot reach it.
            self.best_evaluation = Evaluation(
                x, evaluation.objective, evaluation.status, evaluation.message,
                copy.deepcopy(dict(evaluation.metadata)),
            )
        self._best_so_far.append(self.best_value)
        return evaluation

    def status_codes(self) -> NDArray[np.int8]:
        return np.array(self._status, dtype=np.int8)

    def best_so_far(self) -> NDArray[np.float64]:
        return np.array(self._best_so_far, dtype=float)

    def status_counts(self) -> dict[str, int]:
        codes = self.status_codes()
        return {status.value: int(np.count_nonzero(codes == index)) for index, status in enumerate(STATUS_ORDER)}
