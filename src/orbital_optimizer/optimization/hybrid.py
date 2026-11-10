"""Global search followed by one final local refinement (RF-3.14).

Memetic scheme with "final refinement" coupling (docs/03 section 12.2): a
global optimizer spends a fixed share of the evaluation budget, then
Nelder-Mead refines its best point with what remains. The local phase never
modifies the global population, so the comparison between metaheuristics
stays clean; Lamarckian and Baldwinian couplings are not implemented.

Budget: with a total of B evaluations the global phase receives
B_g = max(floor(global_fraction * B), minimum of the global optimizer) and
stops at B_g or earlier by any other criterion of ``termination``. The local
phase receives exactly B - (evaluations actually used), so the total never
exceeds B and equals B when both phases stop by budget.

Criteria: ``target_objective`` applies to both phases. ``x_tol``, ``f_tol``,
``stagnation_generations`` and ``max_generations`` apply to the global phase
only; the local phase stops by budget or target. If the global phase reaches
the target, or fewer than D+1 evaluations remain, the local phase is skipped
and that is recorded in the result's settings.

Local start: the global best point. Initial simplex edge: the final
normalized population spread of the global phase (the scale at which the
global search stopped resolving), limited to [1e-6, 0.05]; 0.05 is the
default Nelder-Mead edge and 1e-6 keeps the simplex non-degenerate in
float64 arithmetic on a unit range.

The merged history concatenates both phases: per-evaluation arrays in
evaluation order, per-generation arrays with the local iterations appended
(evaluation counts offset, best values as a running minimum). The seed is
passed unchanged to both phases; Nelder-Mead uses it only if the global phase
found no valid point and a random start is needed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np

from .differential_evolution import DifferentialEvolution
from .nelder_mead import ALGORITHM as LOCAL_ALGORITHM, NelderMead, NelderMeadSettings
from .optimizer import Optimizer, validate_problem, validate_seed
from .problem import Problem
from .result import OptimizationHistory, OptimizationResult, TerminationCriteria, TerminationReason

MINIMUM_LOCAL_STEP: float = 1e-6
MAXIMUM_LOCAL_STEP: float = 0.05


@dataclass(frozen=True)
class HybridSettings:
    """global_fraction: share of the evaluation budget of the global phase, in (0, 1); default 0.8.

    The default leaves one fifth of the budget for refinement; it is a design
    parameter recorded with every result, not a tuned value.
    """

    global_fraction: float = 0.8

    def __post_init__(self) -> None:
        fraction = self.global_fraction
        if isinstance(fraction, bool) or not isinstance(fraction, (int, float, np.integer, np.floating)) \
                or not 0.0 < float(fraction) < 1.0:
            raise ValueError("global_fraction must be a real scalar in (0, 1)")
        object.__setattr__(self, "global_fraction", float(fraction))


def merge_results(first: OptimizationResult, second: OptimizationResult, *, algorithm: str,
                  settings: dict[str, Any], termination: TerminationCriteria,
                  wall_time_s: float) -> OptimizationResult:
    """Concatenate two consecutive runs on the same problem into one result.

    The best point is the first run's unless the second found a strictly
    lower valid objective. Populations are dropped (the phases differ in shape).
    """
    take_second = second.found_valid and (not first.found_valid or second.f_best < first.f_best)
    best = second if take_second else first
    first_history, second_history = first.history, second.history
    history = OptimizationHistory(
        best_so_far=np.concatenate((first_history.best_so_far,
                                    np.minimum(second_history.best_so_far, first.f_best))),
        status=np.concatenate((first_history.status, second_history.status)),
        generation_nfev=np.concatenate((first_history.generation_nfev,
                                        first.nfev + second_history.generation_nfev)),
        generation_best=np.concatenate((first_history.generation_best,
                                        np.minimum(second_history.generation_best, first.f_best))),
        generation_valid=np.concatenate((first_history.generation_valid, second_history.generation_valid)),
        generation_spread=np.concatenate((first_history.generation_spread, second_history.generation_spread)),
    )
    counts = {name: int(first.status_counts.get(name, 0)) + int(second.status_counts.get(name, 0))
              for name in set(first.status_counts) | set(second.status_counts)}
    nfev = first.nfev + second.nfev
    generations = len(history.generation_nfev) - 1
    found = best.found_valid
    message = (f"stopped by {second.termination_reason.value} after {nfev} evaluations "
               f"({first.nfev} global + {second.nfev} local) and {generations} generations; "
               + ("best valid objective recorded" if found else "no valid evaluation found"))
    return OptimizationResult(
        x_best=best.x_best, u_best=best.u_best, f_best=best.f_best, best_evaluation=best.best_evaluation,
        found_valid=found, termination_reason=second.termination_reason, message=message, nfev=nfev,
        status_counts=counts, generations=generations, seed=first.seed, algorithm=algorithm,
        settings=settings, termination=termination, problem_info=first.problem_info,
        environment=first.environment, wall_time_s=wall_time_s, history=history,
    )


class HybridOptimizer:
    """Global optimizer plus final Nelder-Mead refinement, implementing the Optimizer protocol."""

    def __init__(self, global_optimizer: Optimizer | None = None,
                 settings: HybridSettings | None = None) -> None:
        if global_optimizer is None:
            global_optimizer = DifferentialEvolution()
        if not callable(getattr(global_optimizer, "minimize", None)) \
                or not isinstance(getattr(global_optimizer, "name", None), str):
            raise ValueError("global_optimizer must implement the Optimizer protocol")
        if settings is None:
            settings = HybridSettings()
        if not isinstance(settings, HybridSettings):
            raise ValueError("settings must be HybridSettings")
        self.global_optimizer = global_optimizer
        self.settings = settings
        self.name = f"Hybrid({global_optimizer.name} -> {LOCAL_ALGORITHM})@1"

    def minimize(self, problem: Problem, termination: TerminationCriteria, *, seed: int) -> OptimizationResult:
        """Minimize ``problem`` under ``termination`` with an explicit integer seed >= 0."""
        bounds = validate_problem(problem)
        if not isinstance(termination, TerminationCriteria):
            raise ValueError("termination must be TerminationCriteria")
        seed = validate_seed(seed)
        dimension = bounds.dimension
        total = termination.max_evaluations
        fraction = self.settings.global_fraction

        start = time.perf_counter()
        share = int(np.floor(fraction * total))
        global_result = None
        # The global optimizer rejects a budget below its own minimum (e.g. its population size);
        # in that case it receives the whole budget and the local phase is skipped.
        for budget in (share, total):
            if budget < 1:
                continue
            try:
                global_result = self.global_optimizer.minimize(
                    problem,
                    TerminationCriteria(
                        max_evaluations=budget, max_generations=termination.max_generations,
                        target_objective=termination.target_objective, x_tol=termination.x_tol,
                        f_tol=termination.f_tol, stagnation_generations=termination.stagnation_generations,
                    ),
                    seed=seed,
                )
                break
            except ValueError:
                if budget == total:
                    raise
        assert global_result is not None
        remaining = total - global_result.nfev
        settings: dict[str, Any] = {
            "global_fraction": fraction,
            "global_algorithm": global_result.algorithm,
            "global_settings": dict(global_result.settings),
            "global_evaluations": global_result.nfev,
            "global_termination_reason": global_result.termination_reason.value,
            "local_algorithm": LOCAL_ALGORITHM,
        }
        skip = None
        if global_result.termination_reason is TerminationReason.TARGET_REACHED:
            skip = "target reached in the global phase"
        elif remaining < dimension + 1:
            skip = f"fewer than D + 1 = {dimension + 1} evaluations remained"
        if skip is not None:
            settings["local_phase"] = f"skipped: {skip}"
            return OptimizationResult(
                x_best=global_result.x_best, u_best=global_result.u_best, f_best=global_result.f_best,
                best_evaluation=global_result.best_evaluation, found_valid=global_result.found_valid,
                termination_reason=global_result.termination_reason, message=global_result.message,
                nfev=global_result.nfev, status_counts=global_result.status_counts,
                generations=global_result.generations, seed=seed, algorithm=self.name, settings=settings,
                termination=termination, problem_info=global_result.problem_info,
                environment=global_result.environment, wall_time_s=time.perf_counter() - start,
                history=global_result.history,
            )

        spread = float(global_result.history.generation_spread[-1])
        step = float(min(MAXIMUM_LOCAL_STEP, max(MINIMUM_LOCAL_STEP, spread)))
        initial = None if global_result.u_best is None else tuple(float(v) for v in global_result.u_best)
        local = NelderMead(NelderMeadSettings(initial_point=initial, initial_step=step))
        local_result = local.minimize(
            problem, TerminationCriteria(max_evaluations=remaining, target_objective=termination.target_objective),
            seed=seed,
        )
        settings["local_phase"] = "executed"
        settings["local_settings"] = dict(local_result.settings)
        settings["local_evaluations"] = local_result.nfev
        return merge_results(global_result, local_result, algorithm=self.name, settings=settings,
                             termination=termination, wall_time_s=time.perf_counter() - start)
