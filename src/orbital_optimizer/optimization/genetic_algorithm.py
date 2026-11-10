"""Real-coded generational genetic algorithm (RF-3.11).

Binary tournament selection, simulated binary crossover (SBX; Deb and
Agrawal, 1995), polynomial mutation (Deb and Goyal, 1996) and elitism;
docs/03 section 11.2. The search runs on normalized vectors u in [0,1]^D
(see problem.Bounds). One generation, for a population of even size NP:

    selection:  parent_i = better of two members drawn with replacement,
                by key = (0, f) if valid else (1, 0); a tie keeps the first.
    crossover:  parents (2k, 2k+1) are recombined with probability p_c.
                For each variable chosen with probability 1/2,
                    beta = (2u)^(1/(eta_c+1))            if u <= 1/2
                           (1/(2(1-u)))^(1/(eta_c+1))    otherwise
                    c1 = ((1+beta) a + (1-beta) b)/2,  c2 = ((1-beta) a + (1+beta) b)/2
                with u ~ U[0,1). For a periodic variable b is first replaced
                by its image nearest to a.
    mutation:   each variable with probability p_m:  x <- x + delta,
                    delta = (2u)^(1/(eta_m+1)) - 1       if u < 1/2
                            1 - (2(1-u))^(1/(eta_m+1))   otherwise
    boundary:   periodic j: wrap modulo 1; non-periodic j: clip to [0,1].
    elitism:    the offspring replace the population; if the best previous
                member is strictly better than the best offspring, it
                replaces the worst offspring (last index among the worst).

Defaults p_c = 0.9, eta_c = 15, p_m = 1/D, eta_m = 20 are the usual values
for these operators in the literature (Deb et al., 2002), not a tuning on
DAEDALUS benchmarks. NP defaults to 10*D, the population rule used for
Differential Evolution.

Random draws (reproducibility contract, "GA/tournament2/SBX/polynomial/elitist@1"),
all from one ``Generator(PCG64(seed))``:
  1. initial population ``rng.random((NP, D))``, rows evaluated in order;
  2. per generation, before any evaluation:
     a. tournaments ``rng.integers(0, NP, size=(NP, 2))``;
     b. ``rng.random(NP // 2)`` (pair recombines if < p_c);
     c. ``rng.random((NP // 2, D))`` (variable takes part if < 1/2);
     d. ``rng.random((NP // 2, D))`` (u of the SBX spread);
     e. ``rng.random((NP, D))`` (variable mutates if < p_m);
     f. ``rng.random((NP, D))`` (u of the mutation step).
Offspring are evaluated in index order while budget remains. No global
optimality is guaranteed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .optimizer import (
    EvaluationBudget, GenerationLog, Key, check_termination, finish_run, make_rng,
    signed_unit_difference, validate_problem, validate_seed, wrap_unit,
)
from .problem import Problem
from .result import OptimizationResult, TerminationCriteria

ALGORITHM: str = "GA/tournament2/SBX/polynomial/elitist@1"


def _real_scalar(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a real scalar")
    return float(value)


@dataclass(frozen=True)
class GeneticAlgorithmSettings:
    """Hyperparameters of the real-coded GA (dimensionless).

    population_size: even NP >= 4; None selects max(4, 10*D).
    crossover_probability: p_c in [0, 1]; default 0.9.
    crossover_index: SBX distribution index eta_c >= 0; default 15.
    mutation_probability: per-variable p_m in [0, 1]; None selects 1/D.
    mutation_index: polynomial distribution index eta_m >= 0; default 20.
    record_population: store every normalized population in the history.
    Invalid values raise ValueError.
    """

    population_size: int | None = None
    crossover_probability: float = 0.9
    crossover_index: float = 15.0
    mutation_probability: float | None = None
    mutation_index: float = 20.0
    record_population: bool = False

    def __post_init__(self) -> None:
        size = self.population_size
        if size is not None and (isinstance(size, bool) or not isinstance(size, (int, np.integer))
                                 or size < 4 or size % 2):
            raise ValueError("population_size must be an even integer >= 4 or None")
        crossover = _real_scalar(self.crossover_probability, "crossover_probability")
        if not 0.0 <= crossover <= 1.0:
            raise ValueError("crossover_probability must lie in [0, 1]")
        mutation = self.mutation_probability
        if mutation is not None:
            mutation = _real_scalar(mutation, "mutation_probability")
            if not 0.0 <= mutation <= 1.0:
                raise ValueError("mutation_probability must lie in [0, 1] or be None")
        indices = []
        for name in ("crossover_index", "mutation_index"):
            index = _real_scalar(getattr(self, name), name)
            if not (index >= 0.0 and np.isfinite(index)):
                raise ValueError(f"{name} must be finite and nonnegative")
            indices.append(index)
        if not isinstance(self.record_population, bool):
            raise ValueError("record_population must be a bool")
        object.__setattr__(self, "population_size", None if size is None else int(size))
        object.__setattr__(self, "crossover_probability", crossover)
        object.__setattr__(self, "mutation_probability", mutation)
        object.__setattr__(self, "crossover_index", indices[0])
        object.__setattr__(self, "mutation_index", indices[1])

    def resolved_population_size(self, dimension: int) -> int:
        return self.population_size if self.population_size is not None else max(4, 10 * dimension)

    def resolved_mutation_probability(self, dimension: int) -> float:
        return self.mutation_probability if self.mutation_probability is not None else 1.0 / dimension


def _sbx_spread(uniforms: NDArray[np.float64], index: float) -> NDArray[np.float64]:
    """SBX spread factor beta(u) for u in [0,1): mean-preserving, beta >= 0."""
    exponent = 1.0 / (index + 1.0)
    lower = (2.0 * uniforms) ** exponent
    with np.errstate(divide="ignore"):
        upper = (1.0 / (2.0 * (1.0 - uniforms))) ** exponent
    return np.where(uniforms <= 0.5, lower, upper)


def _polynomial_step(uniforms: NDArray[np.float64], index: float) -> NDArray[np.float64]:
    """Polynomial mutation step delta(u) in [-1, 1) for u in [0,1)."""
    exponent = 1.0 / (index + 1.0)
    lower = (2.0 * uniforms) ** exponent - 1.0
    upper = 1.0 - (2.0 * (1.0 - uniforms)) ** exponent
    return np.where(uniforms < 0.5, lower, upper)


def _repair(children: NDArray[np.float64], periodic: NDArray[np.bool_]) -> NDArray[np.float64]:
    """Wrap periodic components modulo 1 and clip the others to [0, 1]."""
    repaired = wrap_unit(children, periodic)
    fixed = ~periodic
    if np.any(fixed):
        repaired[:, fixed] = np.clip(repaired[:, fixed], 0.0, 1.0)
    return repaired


class GeneticAlgorithm:
    """Real-coded generational GA implementing the Optimizer protocol."""

    name: str = ALGORITHM

    def __init__(self, settings: GeneticAlgorithmSettings | None = None) -> None:
        if settings is None:
            settings = GeneticAlgorithmSettings()
        if not isinstance(settings, GeneticAlgorithmSettings):
            raise ValueError("settings must be GeneticAlgorithmSettings")
        self.settings = settings

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult:
        """Minimize ``problem`` under ``termination`` with an explicit integer seed >= 0.

        Raises ValueError for an invalid seed or a budget smaller than NP,
        TypeError for an object that is not a Problem.
        """
        bounds = validate_problem(problem)
        if not isinstance(termination, TerminationCriteria):
            raise ValueError("termination must be TerminationCriteria")
        seed = validate_seed(seed)
        dimension = bounds.dimension
        size = self.settings.resolved_population_size(dimension)
        if termination.max_evaluations < size:
            raise ValueError(f"max_evaluations must be at least the population size ({size})")
        p_cross, eta_cross = self.settings.crossover_probability, self.settings.crossover_index
        p_mutation = self.settings.resolved_mutation_probability(dimension)
        eta_mutation = self.settings.mutation_index
        periodic = np.array(bounds.periodic, dtype=bool)
        pairs = size // 2

        start = time.perf_counter()
        rng = make_rng(seed)
        budget = EvaluationBudget(problem, termination.max_evaluations)
        log = GenerationLog(budget, self.settings.record_population)

        population = rng.random((size, dimension))
        keys: list[Key] = []
        valid = np.zeros(size, dtype=bool)
        values = np.full(size, np.inf)
        for i in range(size):
            evaluation = budget.evaluate(population[i])
            keys.append(budget.key(evaluation))
            valid[i] = evaluation.is_valid
            values[i] = evaluation.objective if evaluation.is_valid else np.inf

        generation = 0
        stagnant = 0
        previous_best = budget.best_value
        spread = log.record(population, valid)
        reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)
        while reason is None:
            generation += 1
            tournaments = rng.integers(0, size, size=(size, 2))
            recombine = rng.random(pairs) < p_cross
            take_part = rng.random((pairs, dimension)) < 0.5
            spread_draws = rng.random((pairs, dimension))
            mutate = rng.random((size, dimension)) < p_mutation
            step_draws = rng.random((size, dimension))

            winners = np.array([a if keys[a] <= keys[b] else b for a, b in tournaments.tolist()], dtype=np.int64)
            parents = population[winners]
            first, second = parents[0::2], parents[1::2]
            # Periodic variables: recombine with the image of the second parent nearest to the first.
            second_image = np.where(periodic, first + signed_unit_difference(second, first, periodic), second)
            beta = _sbx_spread(spread_draws, eta_cross)
            active = recombine[:, None] & take_part
            child_a = np.where(active, 0.5 * ((1.0 + beta) * first + (1.0 - beta) * second_image), first)
            child_b = np.where(active, 0.5 * ((1.0 - beta) * first + (1.0 + beta) * second_image), second)
            children = np.empty_like(population)
            children[0::2], children[1::2] = child_a, child_b
            children = np.where(mutate, children + _polynomial_step(step_draws, eta_mutation), children)
            children = _repair(children, periodic)

            elite = min(range(size), key=lambda index: (keys[index], index))
            elite_key, elite_point = keys[elite], population[elite].copy()
            elite_valid, elite_value = bool(valid[elite]), float(values[elite])

            evaluated = 0
            child_keys: list[Key] = []
            child_valid = np.zeros(size, dtype=bool)
            child_values = np.full(size, np.inf)
            for i in range(size):
                if budget.exhausted:
                    break
                evaluation = budget.evaluate(children[i])
                child_keys.append(budget.key(evaluation))
                child_valid[i] = evaluation.is_valid
                child_values[i] = evaluation.objective if evaluation.is_valid else np.inf
                evaluated += 1
            if evaluated == size:
                population, keys, valid, values = children, child_keys, child_valid, child_values
                if elite_key < min(keys):
                    worst = max(range(size), key=lambda index: (keys[index], index))
                    population[worst] = elite_point
                    keys[worst], valid[worst], values[worst] = elite_key, elite_valid, elite_value
            else:
                # Budget ran out mid-generation: evaluated offspring replace members
                # only where they are at least as good, so no evaluated point is lost.
                for i in range(evaluated):
                    if child_keys[i] <= keys[i]:
                        population[i] = children[i]
                        keys[i], valid[i], values[i] = child_keys[i], child_valid[i], child_values[i]
            stagnant = stagnant + 1 if not budget.best_value < previous_best else 0
            previous_best = budget.best_value
            spread = log.record(population, valid)
            reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)

        return finish_run(
            problem=problem, budget=budget, log=log, termination=termination, reason=reason,
            generations=generation, seed=seed, algorithm=ALGORITHM,
            settings={
                "population_size": size, "crossover_probability": p_cross, "crossover_index": eta_cross,
                "mutation_probability": p_mutation, "mutation_index": eta_mutation,
                "selection": "binary tournament with replacement; ties keep the first drawn",
                "crossover": "SBX, each variable with probability 1/2",
                "mutation": "polynomial",
                "boundary_handling": "periodic: wrap modulo 1; others: clip to [0,1]",
                "elitism": "best previous member replaces the worst offspring when strictly better",
                "initialization": "uniform random in [0,1]^D",
                "record_population": self.settings.record_population,
            },
            wall_time_s=time.perf_counter() - start,
        )
