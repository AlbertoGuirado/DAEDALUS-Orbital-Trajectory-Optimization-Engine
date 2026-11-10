"""Canonical Differential Evolution, DE/rand/1/bin with generational selection.

Storn and Price (1997); equations in docs/03 section 11.1. The search runs on
normalized vectors u in [0,1]^D (see problem.Bounds). For generation g:

    mutation:   v_i = u_r1 + F (u_r2 - u_r3),  r1, r2, r3 distinct and != i
    crossover:  t_ij = v_ij if R_ij <= CR or j == jrand_i, else u_ij
    boundary:   periodic j:      t_ij <- t_ij mod 1
                non-periodic j:  t_ij <- (u_ij + 0)/2 if t_ij < 0,
                                 t_ij <- (u_ij + 1)/2 if t_ij > 1
    selection:  u_i <- t_i if key(t_i) <= key(u_i), key = (0, f) if valid
                else (1, 0)

All trials of a generation are built from the population at its start
(generational, not immediate, updating). Invalid evaluations never replace a
valid parent; between two invalid points the trial wins (neutral drift).

Random draws (reproducibility contract, "DE/rand/1/bin@1"), all from one
``Generator(PCG64(seed))``:
  1. initial population: ``rng.random((NP, D))``, rows evaluated in order;
  2. per generation, before any evaluation:
     a. for i = 0..NP-1: ``c = rng.choice(NP - 1, size=3, replace=False)``,
        then c[c >= i] += 1 gives (r1, r2, r3);
     b. ``jrand = rng.integers(0, D, size=NP)``;
     c. ``R = rng.random((NP, D))``.
  Trials are evaluated in index order while budget remains. Changing this
  order is a breaking change of the algorithm version.

Results are bitwise reproducible for the same seed, settings, problem and
environment (Python, NumPy, platform). No global optimality is guaranteed.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .optimizer import (
    EvaluationBudget, check_termination, environment_info, make_rng, problem_info, run_message,
    validate_problem, validate_seed,
)
from .problem import Problem
from .result import OptimizationHistory, OptimizationResult, TerminationCriteria, TerminationReason

logger = logging.getLogger(__name__)

ALGORITHM: str = "DE/rand/1/bin@1"


def _real_scalar(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a real scalar")
    return float(value)


@dataclass(frozen=True)
class DifferentialEvolutionSettings:
    """DE/rand/1/bin hyperparameters (dimensionless).

    population_size: NP >= 4 (rand/1 needs three donors distinct from the
        target); None selects NP = 10*D, the Storn-Price rule of thumb.
    mutation: differential weight F in (0, 2]; default 0.5.
    crossover: binomial crossover probability CR in [0, 1]; default 0.9.
    record_population: store every normalized population in the history.

    Defaults come from the literature (docs/03 section 11.1), not from tuning
    on DAEDALUS benchmarks. Invalid values raise ValueError.
    """

    population_size: int | None = None
    mutation: float = 0.5
    crossover: float = 0.9
    record_population: bool = False

    def __post_init__(self) -> None:
        size = self.population_size
        if size is not None and (isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size < 4):
            raise ValueError("population_size must be an integer >= 4 or None")
        mutation = _real_scalar(self.mutation, "mutation")
        if not 0.0 < mutation <= 2.0:
            raise ValueError("mutation F must lie in (0, 2]")
        crossover = _real_scalar(self.crossover, "crossover")
        if not 0.0 <= crossover <= 1.0:
            raise ValueError("crossover CR must lie in [0, 1]")
        if not isinstance(self.record_population, bool):
            raise ValueError("record_population must be a bool")
        object.__setattr__(self, "population_size", None if size is None else int(size))
        object.__setattr__(self, "mutation", mutation)
        object.__setattr__(self, "crossover", crossover)

    def resolved_population_size(self, dimension: int) -> int:
        """NP for a problem of dimension D."""
        return self.population_size if self.population_size is not None else max(4, 10 * dimension)


def _select_donors(rng: np.random.Generator, population_size: int) -> NDArray[np.int64]:
    """Draw (r1, r2, r3) for every target i: distinct and != i, shape (NP, 3)."""
    donors = np.empty((population_size, 3), dtype=np.int64)
    for i in range(population_size):
        choice = rng.choice(population_size - 1, size=3, replace=False)
        choice[choice >= i] += 1
        donors[i] = choice
    return donors


def _mutate(population: NDArray, donors: NDArray, mutation: float) -> NDArray:
    """v_i = u_r1 + F (u_r2 - u_r3), shape (NP, D)."""
    return population[donors[:, 0]] + mutation * (population[donors[:, 1]] - population[donors[:, 2]])


def _binomial_crossover(population: NDArray, mutants: NDArray, uniforms: NDArray,
                        jrand: NDArray, crossover: float) -> NDArray:
    """t_ij = v_ij if R_ij <= CR or j == jrand_i, else u_ij."""
    take = uniforms <= crossover
    take[np.arange(len(population)), jrand] = True
    return np.where(take, mutants, population)


def _repair(trials: NDArray, parents: NDArray, periodic: NDArray) -> NDArray:
    """Wrap periodic components modulo 1; bounce others midway to the parent's bound."""
    repaired = trials.copy()
    if np.any(periodic):
        wrapped = np.mod(repaired[:, periodic], 1.0)
        # mod of a tiny negative number can round to exactly 1.0, which is 0 + one period.
        wrapped[wrapped >= 1.0] = 0.0
        repaired[:, periodic] = wrapped
    fixed = ~periodic
    if np.any(fixed):
        values, parent = repaired[:, fixed], parents[:, fixed]
        values = np.where(values < 0.0, parent / 2.0, values)
        values = np.where(values > 1.0, (parent + 1.0) / 2.0, values)
        repaired[:, fixed] = values
    return repaired


def _trial_wins(trial_valid: bool, trial_value: float, parent_valid: bool, parent_value: float) -> bool:
    """Selection by key (0, f) for valid points and (1, 0) for invalid ones, with <=.

    It states the selection rule for a problem without constraints. The
    optimizer applies it through ``EvaluationBudget.key``, which gives the
    same decisions there and, for a constrained problem (Phase 4), the
    ranking of its constraint handling instead.
    """
    if trial_valid and parent_valid:
        return trial_value <= parent_value
    if trial_valid != parent_valid:
        return trial_valid
    return True


class DifferentialEvolution:
    """DE/rand/1/bin optimizer implementing the Optimizer protocol."""

    name: str = ALGORITHM

    def __init__(self, settings: DifferentialEvolutionSettings | None = None) -> None:
        if settings is None:
            settings = DifferentialEvolutionSettings()
        if not isinstance(settings, DifferentialEvolutionSettings):
            raise ValueError("settings must be DifferentialEvolutionSettings")
        self.settings = settings

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult:
        """Minimize ``problem`` under ``termination`` with an explicit integer seed >= 0.

        Raises ValueError for an invalid seed or a budget smaller than NP,
        TypeError for an object that is not a Problem. Exceptions raised by
        the problem's evaluate propagate and abort the run.
        """
        bounds = validate_problem(problem)
        if not isinstance(termination, TerminationCriteria):
            raise ValueError("termination must be TerminationCriteria")
        seed = validate_seed(seed)
        dimension = bounds.dimension
        size = self.settings.resolved_population_size(dimension)
        if termination.max_evaluations < size:
            raise ValueError(f"max_evaluations must be at least the population size ({size})")
        mutation, crossover = self.settings.mutation, self.settings.crossover
        periodic = np.array(bounds.periodic, dtype=bool)

        start = time.perf_counter()
        rng = make_rng(seed)
        budget = EvaluationBudget(problem, termination.max_evaluations)
        population = rng.random((size, dimension))
        values = np.full(size, np.inf)
        valid = np.zeros(size, dtype=bool)
        keys = []
        for i in range(size):
            evaluation = budget.evaluate(population[i])
            valid[i] = evaluation.is_valid
            values[i] = evaluation.objective if evaluation.is_valid else np.inf
            keys.append(budget.key(evaluation))

        records: dict[str, list[Any]] = {"nfev": [], "best": [], "valid": [], "spread": [], "populations": []}

        def record() -> float:
            spread = float(np.max(np.ptp(population, axis=0)))
            records["nfev"].append(budget.used)
            records["best"].append(budget.best_value)
            records["valid"].append(int(np.count_nonzero(valid)))
            records["spread"].append(spread)
            if self.settings.record_population:
                records["populations"].append(population.copy())
            return spread

        generation = 0
        stagnant = 0
        previous_best = budget.best_value
        spread = record()
        reason = self._check(termination, budget, generation, spread, valid, values, stagnant)
        while reason is None:
            generation += 1
            donors = _select_donors(rng, size)
            jrand = rng.integers(0, dimension, size=size)
            uniforms = rng.random((size, dimension))
            mutants = _mutate(population, donors, mutation)
            trials = _repair(_binomial_crossover(population, mutants, uniforms, jrand, crossover), population, periodic)
            for i in range(size):
                if budget.exhausted:
                    break
                evaluation = budget.evaluate(trials[i])
                trial_value = evaluation.objective if evaluation.is_valid else np.inf
                # Without constraints this is _trial_wins; with them the handling decides.
                trial_key = budget.key(evaluation)
                if trial_key <= keys[i]:
                    population[i] = trials[i]
                    values[i] = trial_value
                    valid[i] = evaluation.is_valid
                    keys[i] = trial_key
            stagnant = stagnant + 1 if not budget.best_value < previous_best else 0
            previous_best = budget.best_value
            spread = record()
            logger.debug("generation %d: nfev=%d best=%r valid=%d spread=%.3e",
                         generation, budget.used, budget.best_value, records["valid"][-1], spread)
            reason = self._check(termination, budget, generation, spread, valid, values, stagnant)

        wall_time = time.perf_counter() - start
        history = OptimizationHistory(
            best_so_far=budget.best_so_far(),
            status=budget.status_codes(),
            generation_nfev=records["nfev"],
            generation_best=records["best"],
            generation_valid=records["valid"],
            generation_spread=records["spread"],
            populations=np.array(records["populations"]) if self.settings.record_population else None,
            violation=budget.violations() if budget.constrained else None,
        )
        message = run_message(budget, reason, generation)
        logger.info("%s seed=%d: %s", ALGORITHM, seed, message)
        return OptimizationResult(
            x_best=budget.best_x,
            u_best=budget.best_u,
            f_best=budget.best_value,
            best_evaluation=budget.best_evaluation,
            found_valid=budget.valid_count > 0,
            termination_reason=reason,
            message=message,
            nfev=budget.used,
            status_counts=budget.status_counts(),
            generations=generation,
            seed=seed,
            algorithm=ALGORITHM,
            settings={
                "population_size": size,
                "mutation": mutation,
                "crossover": crossover,
                "initialization": "uniform random in [0,1]^D",
                "boundary_handling": "periodic: wrap modulo 1; others: midpoint between parent and violated bound",
                "selection": "generational; trial replaces target if key(trial) <= key(target), "
                             "key=(0,f) valid, (1,0) invalid",
                "record_population": self.settings.record_population,
            },
            termination=termination,
            problem_info=problem_info(problem),
            environment=environment_info(),
            wall_time_s=wall_time,
            history=history,
            **budget.constraint_record(),
        )

    @staticmethod
    def _check(termination: TerminationCriteria, budget: EvaluationBudget, generation: int,
               spread: float, valid: NDArray, values: NDArray, stagnant: int) -> TerminationReason | None:
        """Return the first satisfied criterion in TERMINATION_PRECEDENCE order, or None.

        The rule is shared by every optimizer: optimizer.check_termination.
        """
        return check_termination(termination, budget, generation, spread, valid, values, stagnant)
