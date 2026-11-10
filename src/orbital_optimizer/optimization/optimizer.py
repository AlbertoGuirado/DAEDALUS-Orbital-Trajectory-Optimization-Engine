"""Algorithm-independent optimizer interface, seeding and evaluation budget."""

from __future__ import annotations

import copy
import math
import platform
import sys
from typing import Any, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

from .constraint_handling import ConstraintHandling, FeasibilityRules, Key
from .problem import Bounds, Evaluation, Problem
from .result import (
    STATUS_ORDER, OptimizationHistory, OptimizationResult, TerminationCriteria, TerminationReason,
)

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
    best is the first evaluated point with the lowest *feasible* objective;
    for a problem without constraints every valid point is feasible, so it is
    the lowest valid objective.

    Constraints (Phase 4). ``handling`` ranks candidates during the search
    (``key``); it defaults to ``problem.constraint_handling`` when the
    problem declares one and to Deb's feasibility rules otherwise. Whatever
    the strategy, the reported best is chosen by feasibility: ``best_*`` hold
    the best feasible point and ``least_violation_*`` the valid infeasible
    point with the smallest total violation, kept separately so that an
    infeasible point is never presented as a solution. ``best_value`` -- and
    with it the convergence history, the target criterion and the stagnation
    count -- refers to feasible points only. ``constrained`` becomes True as
    soon as one evaluation carries a constraint violation.

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

    def __init__(self, problem: Problem, max_evaluations: int,
                 handling: ConstraintHandling | None = None) -> None:
        self.problem = problem
        self.bounds = validate_problem(problem)
        if isinstance(max_evaluations, bool) or not isinstance(max_evaluations, (int, np.integer)) or max_evaluations < 1:
            raise ValueError("max_evaluations must be an integer >= 1")
        self.max_evaluations = int(max_evaluations)
        if handling is None:
            handling = getattr(problem, "constraint_handling", None)
        if handling is None:
            handling = FeasibilityRules()
        if not callable(getattr(handling, "key", None)) or not callable(getattr(handling, "describe", None)):
            raise TypeError("constraint handling must define key(evaluation) and describe()")
        self.handling: ConstraintHandling = handling
        self.constrained: bool = False
        self._status: list[int] = []
        self._best_so_far: list[float] = []
        self._violation: list[float] = []
        self.valid_count: int = 0
        self.best_value: float = math.inf
        self.best_x: NDArray[np.float64] | None = None
        self.best_u: NDArray[np.float64] | None = None
        self.best_evaluation: Evaluation | None = None
        self.least_violation: float = math.inf
        self.least_violation_x: NDArray[np.float64] | None = None
        self.least_violation_u: NDArray[np.float64] | None = None
        self.least_violation_evaluation: Evaluation | None = None

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
        if evaluation.is_constrained:
            self.constrained = True
        if evaluation.is_valid:
            self.valid_count += 1
        violation = math.nan if not evaluation.is_valid else (evaluation.constraint_violation or 0.0)
        self._violation.append(violation)
        if evaluation.is_feasible and evaluation.objective < self.best_value:
            self.best_value = evaluation.objective
            self.best_x = x
            self.best_u = np.array(u, dtype=float, copy=True)
            self.best_evaluation = self._own(x, evaluation)
        elif evaluation.is_valid and not evaluation.is_feasible and violation < self.least_violation:
            self.least_violation = violation
            self.least_violation_x = x
            self.least_violation_u = np.array(u, dtype=float, copy=True)
            self.least_violation_evaluation = self._own(x, evaluation)
        self._best_so_far.append(self.best_value)
        return evaluation

    @staticmethod
    def _own(x: NDArray[np.float64], evaluation: Evaluation) -> Evaluation:
        """Budget-owned record: rebuilt from the private x with deep-copied metadata,

        so later changes to the problem's own Evaluation object cannot reach it.
        """
        return Evaluation(
            x, evaluation.objective, evaluation.status, evaluation.message,
            copy.deepcopy(dict(evaluation.metadata)),
            evaluation.constraint_violation, dict(evaluation.constraints),
        )

    def key(self, evaluation: Evaluation) -> Key:
        """Ranking key shared by every optimizer, given by the constraint handling.

        Keys compare lexicographically and lower is better. Without
        constraints it orders valid points by objective and puts invalid
        points last, mutually tied -- the comparison of DE's selection
        (differential_evolution._trial_wins) in key form. With constraints
        the strategy decides (constraint_handling).
        """
        return self.handling.key(evaluation)

    def status_codes(self) -> NDArray[np.int8]:
        return np.array(self._status, dtype=np.int8)

    def best_so_far(self) -> NDArray[np.float64]:
        return np.array(self._best_so_far, dtype=float)

    def violations(self) -> NDArray[np.float64]:
        """Total normalized violation of each evaluation: 0 if feasible, NaN if not valid."""
        return np.array(self._violation, dtype=float)

    def constraint_record(self) -> dict[str, Any]:
        """Keyword arguments of OptimizationResult for the constraint fields; empty without constraints."""
        if not self.constrained:
            return {}
        return {
            "constraint_handling": self.handling.describe(),
            "found_feasible": self.best_evaluation is not None,
            "least_violation": None if self.least_violation_evaluation is None else self.least_violation,
            "x_least_violation": self.least_violation_x,
            "u_least_violation": self.least_violation_u,
            "least_violation_evaluation": self.least_violation_evaluation,
        }

    def status_counts(self) -> dict[str, int]:
        codes = self.status_codes()
        return {status.value: int(np.count_nonzero(codes == index)) for index, status in enumerate(STATUS_ORDER)}


def check_termination(termination: TerminationCriteria, budget: EvaluationBudget, generation: int,
                      spread: float, valid: NDArray, values: NDArray, stagnant: int) -> TerminationReason | None:
    """Return the first satisfied criterion in TERMINATION_PRECEDENCE order, or None.

    ``spread`` is the largest per-variable normalized range of the retained
    population, ``valid``/``values`` its validity flags and objectives (+inf
    for invalid members) and ``stagnant`` the number of consecutive
    generations without a strict decrease of the best valid objective.
    """
    if budget.exhausted:
        return TerminationReason.MAX_EVALUATIONS
    if termination.target_objective is not None and budget.best_value <= termination.target_objective:
        return TerminationReason.TARGET_REACHED
    if termination.x_tol is not None and spread <= termination.x_tol:
        return TerminationReason.POPULATION_CONVERGED
    if (termination.f_tol is not None and np.all(valid)
            and float(np.max(values) - np.min(values)) <= termination.f_tol):
        return TerminationReason.OBJECTIVE_CONVERGED
    if (termination.stagnation_generations is not None and generation > 0
            and stagnant >= termination.stagnation_generations):
        return TerminationReason.STAGNATION
    if termination.max_generations is not None and generation >= termination.max_generations:
        return TerminationReason.MAX_GENERATIONS
    return None


class GenerationLog:
    """Per-generation convergence record shared by the optimizers (RF-3.18).

    ``record`` is called once after the initial population (generation 0) and
    once after every later generation with the retained normalized population,
    shape (N, D), and its validity flags, shape (N,). It returns the spread,
    max_j (max_i u_ij - min_i u_ij).
    """

    def __init__(self, budget: EvaluationBudget, record_population: bool = False) -> None:
        self.budget = budget
        self.record_population = bool(record_population)
        self.nfev: list[int] = []
        self.best: list[float] = []
        self.valid: list[int] = []
        self.spread: list[float] = []
        self.populations: list[NDArray[np.float64]] = []

    def record(self, population: NDArray[np.float64], valid: NDArray[np.bool_]) -> float:
        spread = float(np.max(np.ptp(population, axis=0)))
        self.nfev.append(self.budget.used)
        self.best.append(self.budget.best_value)
        self.valid.append(int(np.count_nonzero(valid)))
        self.spread.append(spread)
        if self.record_population:
            self.populations.append(np.array(population, dtype=float, copy=True))
        return spread

    def history(self) -> OptimizationHistory:
        return OptimizationHistory(
            best_so_far=self.budget.best_so_far(),
            status=self.budget.status_codes(),
            generation_nfev=self.nfev,
            generation_best=self.best,
            generation_valid=self.valid,
            generation_spread=self.spread,
            populations=np.array(self.populations) if self.record_population else None,
            violation=self.budget.violations() if self.budget.constrained else None,
        )


def run_message(budget: EvaluationBudget, reason: TerminationReason, generations: int) -> str:
    """One-line outcome of a run; it names feasibility only when the problem has constraints."""
    head = f"stopped by {reason.value} after {budget.used} evaluations and {generations} generations; "
    if not budget.constrained:
        return head + ("best valid objective recorded" if budget.best_evaluation is not None
                       else "no valid evaluation found")
    if budget.best_evaluation is not None:
        return head + "best feasible objective recorded"
    if budget.least_violation_evaluation is not None:
        return head + f"no feasible point found (least total violation {budget.least_violation:.6g})"
    return head + "no valid evaluation found"


def finish_run(*, problem: Problem, budget: EvaluationBudget, log: GenerationLog,
               termination: TerminationCriteria, reason: TerminationReason, generations: int,
               seed: int, algorithm: str, settings: dict[str, Any], wall_time_s: float) -> OptimizationResult:
    """Assemble the OptimizationResult of a finished run from its budget and log."""
    return OptimizationResult(
        x_best=budget.best_x,
        u_best=budget.best_u,
        f_best=budget.best_value,
        best_evaluation=budget.best_evaluation,
        found_valid=budget.valid_count > 0,
        termination_reason=reason,
        message=run_message(budget, reason, generations),
        nfev=budget.used,
        status_counts=budget.status_counts(),
        generations=generations,
        seed=seed,
        algorithm=algorithm,
        settings=settings,
        termination=termination,
        problem_info=problem_info(problem),
        environment=environment_info(),
        wall_time_s=wall_time_s,
        history=log.history(),
        **budget.constraint_record(),
    )


def wrap_unit(values: NDArray[np.float64], periodic: NDArray[np.bool_]) -> NDArray[np.float64]:
    """Return a copy with periodic components reduced modulo 1 into [0, 1).

    The modulo of a tiny negative number can round to exactly 1.0, which is
    the same point as 0.0 one period later; it is mapped to 0.0.
    """
    result = np.array(values, dtype=float, copy=True)
    if np.any(periodic):
        wrapped = np.mod(result[..., periodic], 1.0)
        wrapped[wrapped >= 1.0] = 0.0
        result[..., periodic] = wrapped
    return result


def signed_unit_difference(target: NDArray[np.float64], origin: NDArray[np.float64],
                           periodic: NDArray[np.bool_]) -> NDArray[np.float64]:
    """target - origin, using the shortest signed arc in [-0.5, 0.5) for periodic components."""
    difference = np.asarray(target, dtype=float) - np.asarray(origin, dtype=float)
    if np.any(periodic):
        difference = np.array(difference, copy=True)
        difference[..., periodic] = np.mod(difference[..., periodic] + 0.5, 1.0) - 0.5
    return difference
