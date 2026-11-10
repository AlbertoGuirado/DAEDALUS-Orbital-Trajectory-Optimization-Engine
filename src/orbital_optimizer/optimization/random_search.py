"""Pure uniform random search: the null baseline of every optimizer comparison.

Each candidate is an independent draw u ~ U[0,1)^D mapped to the physical
box by problem.Bounds; nothing is learned from earlier evaluations. A
metaheuristic that does not beat this baseline under the same evaluation
budget has shown no search ability on that problem (docs/03 section 17).

Random draws (reproducibility contract, "RandomSearch/uniform@1"), all from
one ``Generator(PCG64(seed))``: batch g = 0, 1, ... is ``rng.random((B, D))``
and its rows are evaluated in order while budget remains. The sequence of
evaluated points therefore does not depend on the budget (prefix property).

A batch plays the role of a generation in the convergence history: B = 10*D
by default, the population rule of the other optimizers, so that histories
have comparable granularity. The batch size does not change which points are
sampled for a fixed number of evaluations only if B is unchanged.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from .optimizer import (
    EvaluationBudget, GenerationLog, check_termination, finish_run, make_rng, validate_problem,
    validate_seed,
)
from .problem import Problem
from .result import OptimizationResult, TerminationCriteria

ALGORITHM: str = "RandomSearch/uniform@1"


@dataclass(frozen=True)
class RandomSearchSettings:
    """batch_size: samples per recorded iteration (>= 1); None selects max(4, 10*D).

    record_population: store every normalized batch in the history.
    """

    batch_size: int | None = None
    record_population: bool = False

    def __post_init__(self) -> None:
        size = self.batch_size
        if size is not None and (isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size < 1):
            raise ValueError("batch_size must be an integer >= 1 or None")
        if not isinstance(self.record_population, bool):
            raise ValueError("record_population must be a bool")
        object.__setattr__(self, "batch_size", None if size is None else int(size))

    def resolved_batch_size(self, dimension: int) -> int:
        return self.batch_size if self.batch_size is not None else max(4, 10 * dimension)


class RandomSearch:
    """Uniform random sampling implementing the Optimizer protocol."""

    name: str = ALGORITHM

    def __init__(self, settings: RandomSearchSettings | None = None) -> None:
        if settings is None:
            settings = RandomSearchSettings()
        if not isinstance(settings, RandomSearchSettings):
            raise ValueError("settings must be RandomSearchSettings")
        self.settings = settings

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult:
        """Sample ``problem`` uniformly under ``termination`` with an explicit integer seed >= 0."""
        bounds = validate_problem(problem)
        if not isinstance(termination, TerminationCriteria):
            raise ValueError("termination must be TerminationCriteria")
        seed = validate_seed(seed)
        dimension = bounds.dimension
        size = self.settings.resolved_batch_size(dimension)

        start = time.perf_counter()
        rng = make_rng(seed)
        budget = EvaluationBudget(problem, termination.max_evaluations)
        log = GenerationLog(budget, self.settings.record_population)

        generation = -1
        stagnant = 0
        previous_best = np.inf
        reason = None
        while reason is None:
            generation += 1
            batch = rng.random((size, dimension))
            valid = np.zeros(size, dtype=bool)
            values = np.full(size, np.inf)
            for i in range(size):
                if budget.exhausted:
                    break
                evaluation = budget.evaluate(batch[i])
                valid[i] = evaluation.is_valid
                values[i] = evaluation.objective if evaluation.is_valid else np.inf
            if generation > 0:
                stagnant = stagnant + 1 if not budget.best_value < previous_best else 0
            previous_best = budget.best_value
            spread = log.record(batch, valid)
            reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)

        return finish_run(
            problem=problem, budget=budget, log=log, termination=termination, reason=reason,
            generations=generation, seed=seed, algorithm=ALGORITHM,
            settings={"batch_size": size, "sampling": "independent uniform in [0,1)^D",
                      "record_population": self.settings.record_population},
            wall_time_s=time.perf_counter() - start,
        )
