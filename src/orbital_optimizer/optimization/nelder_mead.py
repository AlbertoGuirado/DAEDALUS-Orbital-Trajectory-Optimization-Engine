"""Nelder-Mead simplex search for local refinement (RF-3.13).

Nelder and Mead (1965) with the standard coefficients analysed by Lagarias,
Reeds, Wright and Wright (1998): reflection 1, expansion 2, contraction 1/2,
shrink 1/2; docs/03 section 12.1. It is a derivative-free *local* method: it
only compares candidates, so it uses the ranking key of the other optimizers
(valid points first, then lower objective) and treats invalid points as
worst. It is not a global optimizer and gives no optimality guarantee.

The simplex has D+1 vertices in normalized coordinates. With the vertices
ordered best to worst (stable order) and c the centroid of all but the worst
vertex w:

    reflect   r = c + (c - w)
    if key(best) <= key(r) < key(second worst):  accept r
    elif key(r) < key(best):   expand e = c + 2 (c - w); accept e if
                               key(e) < key(r), else r
    elif key(r) < key(w):      outside contraction o = c + (r - c)/2;
                               accept o if key(o) <= key(r), else shrink
    else:                      inside contraction i = c + (w - c)/2;
                               accept i if key(i) < key(w), else shrink
    shrink    v_k <- best + (v_k - best)/2 for every other vertex

Boundary handling: a non-periodic coordinate of a trial point that leaves
[0,1] is mirrored at the violated face (u < 0 -> -u, u > 1 -> 2 - u) and the
mirrored point is the one stored and evaluated. Mirroring keeps the simplex
full-dimensional; clipping instead makes it collapse onto a face, where it can
converge to a non-stationary point (observed on a convex quadratic with an
interior minimum, 2026-11-10). A periodic coordinate is kept unwrapped in the
simplex and reduced modulo 1 only for evaluation, so the simplex can cross
the seam of an angle.

Start: ``settings.initial_point`` (normalized, shape (D,)) or, when it is
None, ``rng.random(D)`` from ``Generator(PCG64(seed))`` -- the only random
draw of the algorithm ("NelderMead/Lagarias1998@1"). Initial vertices are the
start u0 and u0 + s e_j for each axis j (u0 - s e_j when a non-periodic
coordinate would leave the box), with s = ``initial_step`` = 0.05 by default:
5 % of each variable's range, the conventional relative simplex size.

One iteration is one generation of the convergence history; the recorded
population is the simplex and its spread is the largest per-variable range
of the vertices.

Known limitations. (1) The method can converge to a non-stationary point even
on smooth convex functions (McKinnon, 1998). (2) Invalid evaluations act as a
discontinuous barrier: next to a region of invalid points the simplex can
collapse against it and stall (observed 2026-11-10: linear objective beside an
invalid half-plane, stalled at f = 0.033 with minimum 0). An optimum lying on
such a boundary -- e.g. a tangency limit of an orbital problem -- may
therefore not be refined. (3) Started where every point is invalid, all
vertices tie and the search only shrinks; it reports that no valid point was
found.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .optimizer import (
    EvaluationBudget, GenerationLog, Key, check_termination, finish_run, make_rng, validate_problem,
    validate_seed, wrap_unit,
)
from .problem import Problem
from .result import OptimizationResult, TerminationCriteria

ALGORITHM: str = "NelderMead/Lagarias1998@1"

REFLECTION: float = 1.0
EXPANSION: float = 2.0
CONTRACTION: float = 0.5
SHRINK: float = 0.5


@dataclass(frozen=True)
class NelderMeadSettings:
    """initial_point: normalized start u0 in [0,1]^D, or None for a seeded random start.

    initial_step: normalized simplex edge s in (0, 0.5]; default 0.05.
    record_population: store the simplex of every iteration.
    Invalid values raise ValueError.
    """

    initial_point: tuple[float, ...] | None = None
    initial_step: float = 0.05
    record_population: bool = False

    def __post_init__(self) -> None:
        if self.initial_point is not None:
            try:
                point = np.array(self.initial_point, dtype=float)
            except (TypeError, ValueError) as exc:
                raise ValueError("initial_point must be a sequence of real numbers or None") from exc
            if point.ndim != 1 or point.size < 1 or not np.all(np.isfinite(point)) \
                    or np.any(point < 0.0) or np.any(point > 1.0):
                raise ValueError("initial_point must be a finite vector in [0,1]^D")
            object.__setattr__(self, "initial_point", tuple(float(value) for value in point))
        step = self.initial_step
        if isinstance(step, bool) or not isinstance(step, (int, float, np.integer, np.floating)) \
                or not 0.0 < float(step) <= 0.5:
            raise ValueError("initial_step must be a real scalar in (0, 0.5]")
        if not isinstance(self.record_population, bool):
            raise ValueError("record_population must be a bool")
        object.__setattr__(self, "initial_step", float(step))


def _mirror(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Reflect coordinates at the faces of [0,1]: u < 0 -> -u, u > 1 -> 2 - u; then clip.

    One reflection brings back any point within one range of the box; the
    final clip only acts on points further away than that.
    """
    mirrored = np.where(values < 0.0, -values, values)
    mirrored = np.where(mirrored > 1.0, 2.0 - mirrored, mirrored)
    return np.clip(mirrored, 0.0, 1.0)


class _Simplex:
    """Vertices (stored coordinates), their keys and validity, kept in arbitrary order."""

    def __init__(self, budget: EvaluationBudget, periodic: NDArray[np.bool_]) -> None:
        self.budget = budget
        self.periodic = periodic
        self.fixed = ~periodic
        self.points: list[NDArray[np.float64]] = []
        self.keys: list[Key] = []
        self.valid: list[bool] = []
        self.values: list[float] = []

    def evaluate(self, point: NDArray[np.float64]) -> tuple[NDArray[np.float64], Key, bool, float]:
        """Mirror non-periodic coordinates into the box, evaluate the wrapped point, return its record."""
        stored = np.array(point, dtype=float, copy=True)
        stored[self.fixed] = _mirror(stored[self.fixed])
        evaluation = self.budget.evaluate(wrap_unit(stored, self.periodic))
        value = evaluation.objective if evaluation.is_valid else np.inf
        return stored, self.budget.key(evaluation), evaluation.is_valid, value

    def append(self, record: tuple[NDArray[np.float64], Key, bool, float]) -> None:
        point, key, valid, value = record
        self.points.append(point)
        self.keys.append(key)
        self.valid.append(valid)
        self.values.append(value)

    def replace(self, index: int, record: tuple[NDArray[np.float64], Key, bool, float]) -> None:
        self.points[index], self.keys[index], self.valid[index], self.values[index] = record

    def sort(self) -> None:
        order = sorted(range(len(self.keys)), key=lambda index: (self.keys[index], index))
        self.points = [self.points[index] for index in order]
        self.keys = [self.keys[index] for index in order]
        self.valid = [self.valid[index] for index in order]
        self.values = [self.values[index] for index in order]

    def array(self) -> NDArray[np.float64]:
        return np.array(self.points)


class NelderMead:
    """Nelder-Mead local search implementing the Optimizer protocol."""

    name: str = ALGORITHM

    def __init__(self, settings: NelderMeadSettings | None = None) -> None:
        if settings is None:
            settings = NelderMeadSettings()
        if not isinstance(settings, NelderMeadSettings):
            raise ValueError("settings must be NelderMeadSettings")
        self.settings = settings

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult:
        """Refine from the configured start under ``termination``; seed >= 0 is used only for a random start.

        Raises ValueError for an invalid seed, a start of the wrong dimension
        or a budget smaller than D+1; TypeError for an object that is not a Problem.
        """
        bounds = validate_problem(problem)
        if not isinstance(termination, TerminationCriteria):
            raise ValueError("termination must be TerminationCriteria")
        seed = validate_seed(seed)
        dimension = bounds.dimension
        if termination.max_evaluations < dimension + 1:
            raise ValueError(f"max_evaluations must be at least D + 1 ({dimension + 1})")
        periodic = np.array(bounds.periodic, dtype=bool)
        step = self.settings.initial_step

        start_time = time.perf_counter()
        rng = make_rng(seed)
        if self.settings.initial_point is None:
            start = rng.random(dimension)
            origin = "seeded uniform random point"
        else:
            start = np.array(self.settings.initial_point, dtype=float)
            if start.shape != (dimension,):
                raise ValueError(f"initial_point must have shape ({dimension},)")
            origin = "given initial point"
        budget = EvaluationBudget(problem, termination.max_evaluations)
        log = GenerationLog(budget, self.settings.record_population)
        simplex = _Simplex(budget, periodic)

        simplex.append(simplex.evaluate(start))
        for j in range(dimension):
            vertex = start.copy()
            vertex[j] += step if periodic[j] or start[j] + step <= 1.0 else -step
            simplex.append(simplex.evaluate(vertex))
        simplex.sort()

        generation = 0
        stagnant = 0
        previous_best = budget.best_value
        operations = {"reflection": 0, "expansion": 0, "outside_contraction": 0,
                      "inside_contraction": 0, "shrink": 0}

        def snapshot() -> tuple[float, NDArray[np.bool_], NDArray[np.float64]]:
            valid = np.array(simplex.valid, dtype=bool)
            return log.record(simplex.array(), valid), valid, np.array(simplex.values, dtype=float)

        spread, valid, values = snapshot()
        reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)
        while reason is None:
            generation += 1
            self._iterate(simplex, budget, operations)
            simplex.sort()
            stagnant = stagnant + 1 if not budget.best_value < previous_best else 0
            previous_best = budget.best_value
            spread, valid, values = snapshot()
            reason = check_termination(termination, budget, generation, spread, valid, values, stagnant)

        settings: dict[str, Any] = {
            "initial_point": None if self.settings.initial_point is None else list(self.settings.initial_point),
            "start": origin,
            "initial_step": step,
            "coefficients": {"reflection": REFLECTION, "expansion": EXPANSION,
                             "contraction": CONTRACTION, "shrink": SHRINK},
            "boundary_handling": "periodic: unwrapped simplex, modulo 1 at evaluation; "
                                 "others: mirror at the violated face",
            "operations": operations,
            "record_population": self.settings.record_population,
        }
        return finish_run(
            problem=problem, budget=budget, log=log, termination=termination, reason=reason,
            generations=generation, seed=seed, algorithm=ALGORITHM, settings=settings,
            wall_time_s=time.perf_counter() - start_time,
        )

    @staticmethod
    def _iterate(simplex: _Simplex, budget: EvaluationBudget, operations: dict[str, int]) -> None:
        """One Nelder-Mead iteration on a sorted simplex; stops early if the budget runs out."""
        if budget.exhausted:
            return
        points = simplex.array()
        best_key, second_worst_key, worst_key = simplex.keys[0], simplex.keys[-2], simplex.keys[-1]
        worst = points[-1]
        centroid = np.mean(points[:-1], axis=0)
        reflected = simplex.evaluate(centroid + REFLECTION * (centroid - worst))
        if best_key <= reflected[1] < second_worst_key:
            simplex.replace(-1, reflected)
            operations["reflection"] += 1
            return
        if reflected[1] < best_key:
            if budget.exhausted:
                simplex.replace(-1, reflected)
                operations["reflection"] += 1
                return
            expanded = simplex.evaluate(centroid + EXPANSION * (centroid - worst))
            if expanded[1] < reflected[1]:
                simplex.replace(-1, expanded)
                operations["expansion"] += 1
            else:
                simplex.replace(-1, reflected)
                operations["reflection"] += 1
            return
        outside = reflected[1] < worst_key
        if budget.exhausted:
            if outside:
                simplex.replace(-1, reflected)
                operations["reflection"] += 1
            return
        if outside:
            contracted = simplex.evaluate(centroid + CONTRACTION * (reflected[0] - centroid))
            accepted = contracted[1] <= reflected[1]
            name = "outside_contraction"
        else:
            contracted = simplex.evaluate(centroid + CONTRACTION * (worst - centroid))
            accepted = contracted[1] < worst_key
            name = "inside_contraction"
        if accepted:
            simplex.replace(-1, contracted)
            operations[name] += 1
            return
        operations["shrink"] += 1
        best = simplex.points[0]
        for index in range(1, len(simplex.points)):
            if budget.exhausted:
                return
            simplex.replace(index, simplex.evaluate(best + SHRINK * (simplex.points[index] - best)))
