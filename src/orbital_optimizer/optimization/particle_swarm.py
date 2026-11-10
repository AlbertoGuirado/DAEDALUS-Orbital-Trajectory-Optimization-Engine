"""Global-best particle swarm optimization with inertia weight (RF-3.12).

Kennedy and Eberhart (1995); inertia weight of Shi and Eberhart (1998);
docs/03 section 11.3. The search runs on normalized vectors u in [0,1]^D
(see problem.Bounds). For particle i with position x_i, velocity v_i,
personal best p_i and swarm best g:

    v_i <- w v_i + c1 R1_i * d(p_i, x_i) + c2 R2_i * d(g, x_i)
    v_ij <- clip(v_ij, -v_max, v_max)
    x_i <- x_i + v_i

with componentwise products, R1, R2 ~ U[0,1)^D and d(a, b) = a - b, except
for periodic variables, where d is the shortest signed arc in [-0.5, 0.5).

    boundary:   periodic j:      x_ij <- x_ij mod 1 (velocity kept)
                non-periodic j:  x_ij <- 0 or 1 and v_ij <- 0 (absorbing wall)
    memory:     p_i <- x_i if key(x_i) <= key(p_i), key = (0, f) if valid
                else (1, 0)
    swarm best: g <- first p_i with the lowest key, updated once per
                iteration after all particles moved (synchronous update)

Defaults w = 0.7298, c1 = c2 = 1.49618 are the constriction-equivalent
constants of Clerc and Kennedy (2002), a literature choice and not a tuning
on DAEDALUS benchmarks. v_max = 0.5 limits one step to half the normalized
range. The swarm size defaults to 10*D, the population rule used for
Differential Evolution, so both run the same number of iterations under an
equal evaluation budget.

Random draws (reproducibility contract, "PSO/gbest/inertia@1"), all from one
``Generator(PCG64(seed))``:
  1. positions ``rng.random((NP, D))``, rows evaluated in order;
  2. velocities ``rng.random((NP, D)) - 0.5``;
  3. per iteration, before any evaluation: ``R1 = rng.random((NP, D))`` then
     ``R2 = rng.random((NP, D))``.
Particles are evaluated in index order while budget remains. The convergence
history records the personal bests (the retained solutions). No global
optimality is guaranteed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from .optimizer import (
    EvaluationBudget, GenerationLog, Key, check_termination, finish_run, make_rng,
    signed_unit_difference, validate_problem, validate_seed, wrap_unit,
)
from .problem import Problem
from .result import OptimizationResult, TerminationCriteria

ALGORITHM: str = "PSO/gbest/inertia@1"


def _real_scalar(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a real scalar")
    return float(value)


@dataclass(frozen=True)
class ParticleSwarmSettings:
    """Global-best PSO hyperparameters (dimensionless).

    swarm_size: NP >= 2; None selects max(4, 10*D).
    inertia: w in [0, 1]; default 0.7298.
    cognitive, social: c1, c2 >= 0 with c1 + c2 > 0; default 1.49618 each.
    max_velocity: per-variable clamp in (0, 1], fraction of the normalized range.
    record_population: store the personal bests of every iteration.
    Invalid values raise ValueError.
    """

    swarm_size: int | None = None
    inertia: float = 0.7298
    cognitive: float = 1.49618
    social: float = 1.49618
    max_velocity: float = 0.5
    record_population: bool = False

    def __post_init__(self) -> None:
        size = self.swarm_size
        if size is not None and (isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size < 2):
            raise ValueError("swarm_size must be an integer >= 2 or None")
        inertia = _real_scalar(self.inertia, "inertia")
        cognitive = _real_scalar(self.cognitive, "cognitive")
        social = _real_scalar(self.social, "social")
        max_velocity = _real_scalar(self.max_velocity, "max_velocity")
        if not 0.0 <= inertia <= 1.0:
            raise ValueError("inertia w must lie in [0, 1]")
        if not (cognitive >= 0.0 and social >= 0.0 and cognitive + social > 0.0
                and np.isfinite(cognitive + social)):
            raise ValueError("cognitive and social must be finite, nonnegative and not both zero")
        if not 0.0 < max_velocity <= 1.0:
            raise ValueError("max_velocity must lie in (0, 1]")
        if not isinstance(self.record_population, bool):
            raise ValueError("record_population must be a bool")
        object.__setattr__(self, "swarm_size", None if size is None else int(size))
        object.__setattr__(self, "inertia", inertia)
        object.__setattr__(self, "cognitive", cognitive)
        object.__setattr__(self, "social", social)
        object.__setattr__(self, "max_velocity", max_velocity)

    def resolved_swarm_size(self, dimension: int) -> int:
        return self.swarm_size if self.swarm_size is not None else max(4, 10 * dimension)


def _first_best(keys: list[Key]) -> int:
    """Index of the first lowest key."""
    best = 0
    for index in range(1, len(keys)):
        if keys[index] < keys[best]:
            best = index
    return best


class ParticleSwarm:
    """Global-best PSO implementing the Optimizer protocol."""

    name: str = ALGORITHM

    def __init__(self, settings: ParticleSwarmSettings | None = None) -> None:
        if settings is None:
            settings = ParticleSwarmSettings()
        if not isinstance(settings, ParticleSwarmSettings):
            raise ValueError("settings must be ParticleSwarmSettings")
        self.settings = settings

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult:
        """Minimize ``problem`` under ``termination`` with an explicit integer seed >= 0.

        Raises ValueError for an invalid seed or a budget smaller than the
        swarm, TypeError for an object that is not a Problem.
        """
        bounds = validate_problem(problem)
        if not isinstance(termination, TerminationCriteria):
            raise ValueError("termination must be TerminationCriteria")
        seed = validate_seed(seed)
        dimension = bounds.dimension
        size = self.settings.resolved_swarm_size(dimension)
        if termination.max_evaluations < size:
            raise ValueError(f"max_evaluations must be at least the swarm size ({size})")
        inertia, cognitive, social = self.settings.inertia, self.settings.cognitive, self.settings.social
        limit = self.settings.max_velocity
        periodic = np.array(bounds.periodic, dtype=bool)
        fixed = ~periodic

        start = time.perf_counter()
        rng = make_rng(seed)
        budget = EvaluationBudget(problem, termination.max_evaluations)
        log = GenerationLog(budget, self.settings.record_population)

        positions = rng.random((size, dimension))
        velocities = rng.random((size, dimension)) - 0.5
        best_positions = positions.copy()
        best_keys: list[Key] = []
        valid = np.zeros(size, dtype=bool)
        values = np.full(size, np.inf)
        for i in range(size):
            evaluation = budget.evaluate(positions[i])
            best_keys.append(budget.key(evaluation))
            valid[i] = evaluation.is_valid
            values[i] = evaluation.objective if evaluation.is_valid else np.inf
        leader = _first_best(best_keys)

        generation = 0
        stagnant = 0
        previous_best = budget.best_value
        spread = log.record(best_positions, valid)
        reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)
        while reason is None:
            generation += 1
            first = rng.random((size, dimension))
            second = rng.random((size, dimension))
            to_personal = signed_unit_difference(best_positions, positions, periodic)
            to_leader = signed_unit_difference(best_positions[leader], positions, periodic)
            velocities = inertia * velocities + cognitive * first * to_personal + social * second * to_leader
            velocities = np.clip(velocities, -limit, limit)
            positions = positions + velocities
            if np.any(fixed):
                outside = fixed & ((positions < 0.0) | (positions > 1.0))
                positions = np.where(fixed, np.clip(positions, 0.0, 1.0), positions)
                velocities = np.where(outside, 0.0, velocities)
            positions = wrap_unit(positions, periodic)
            for i in range(size):
                if budget.exhausted:
                    break
                evaluation = budget.evaluate(positions[i])
                key = budget.key(evaluation)
                if key <= best_keys[i]:
                    best_keys[i] = key
                    best_positions[i] = positions[i]
                    valid[i] = evaluation.is_valid
                    values[i] = evaluation.objective if evaluation.is_valid else np.inf
            leader = _first_best(best_keys)
            stagnant = stagnant + 1 if not budget.best_value < previous_best else 0
            previous_best = budget.best_value
            spread = log.record(best_positions, valid)
            reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)

        return finish_run(
            problem=problem, budget=budget, log=log, termination=termination, reason=reason,
            generations=generation, seed=seed, algorithm=ALGORITHM,
            settings={
                "swarm_size": size, "inertia": inertia, "cognitive": cognitive, "social": social,
                "max_velocity": limit,
                "initialization": "positions uniform in [0,1)^D; velocities uniform in [-0.5,0.5)^D",
                "boundary_handling": "periodic: wrap modulo 1; others: clip to the bound and zero the velocity",
                "topology": "global best, synchronous update",
                "record_population": self.settings.record_population,
            },
            wall_time_s=time.perf_counter() - start,
        )
