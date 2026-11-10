"""Run optimizers over fixed scenarios and independent seeds (RF-6A.01 to RF-6A.04).

Benchmark configuration (JSON-compatible, schema ``daedalus.benchmark/1``):

    {"schema": "daedalus.benchmark/1", "name": ..., "description": ...,
     "seeds": {"first": 0, "count": 30},
     "algorithms": [{"id": "de", "type": ..., "settings": {...}}, ...],
     "scenarios": [{"id": ..., "type": ..., "max_evaluations": ..., ...}, ...],
     "acceptance": {...}}        # optional, interpreted by the caller

Rules enforced here:
  - every algorithm runs every scenario with exactly the scenario's
    ``max_evaluations`` and no other stopping criterion, so all algorithms
    receive the same evaluation budget (RF-6A.03);
  - the same seeds are used for every (scenario, algorithm) pair, and each
    run uses its own generator ``Generator(PCG64(seed))`` (RR-01);
  - wall time is recorded and never used as a stopping criterion.

Runs are independent, so they can execute in worker processes (``jobs``).
The records do not depend on ``jobs`` except for ``wall_time_s``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..optimization.optimizer import Optimizer
from ..optimization.result import OptimizationResult, TerminationCriteria
from .registry import Scenario, build_optimizer, build_scenario
from .statistics import evaluations_to_target

SCHEMA: str = "daedalus.benchmark/1"


@dataclass(frozen=True, eq=False)
class RunRecord:
    """One optimizer run on one scenario with one seed.

    f_best is +inf when no valid point was found. gap = f_best - reference
    (None without a reference). success: see registry.Scenario (None when the
    scenario declares no threshold). anomaly: the run is better than the
    reference by more than the scenario's lower threshold.
    evaluations_to_success: evaluations until the best-so-far value first met
    the success threshold (None if never). metrics: scalars copied from the
    best evaluation. verification: record of the scenario's verifier, or
    None. best_so_far: read-only (nfev,) convergence history.
    """

    scenario: str
    algorithm: str
    algorithm_name: str
    seed: int
    f_best: float
    found_valid: bool
    gap: float | None
    success: bool | None
    anomaly: bool
    evaluations_to_success: int | None
    nfev: int
    n_invalid: int
    generations: int
    termination_reason: str
    wall_time_s: float
    x_best: tuple[float, ...] | None
    metrics: Mapping[str, float]
    verification: Mapping[str, Any] | None
    best_so_far: NDArray[np.float64]

    def to_row(self) -> dict[str, Any]:
        """Flat, JSON-compatible record without the convergence history."""
        row: dict[str, Any] = {
            "scenario": self.scenario, "algorithm": self.algorithm, "algorithm_name": self.algorithm_name,
            "seed": self.seed, "f_best": self.f_best, "found_valid": self.found_valid, "gap": self.gap,
            "success": self.success, "anomaly": self.anomaly,
            "evaluations_to_success": self.evaluations_to_success, "nfev": self.nfev,
            "n_invalid": self.n_invalid, "generations": self.generations,
            "termination_reason": self.termination_reason, "wall_time_s": self.wall_time_s,
            "verification_passed": None if self.verification is None else bool(self.verification["passed"]),
            "x_best": None if self.x_best is None else list(self.x_best),
        }
        row.update({f"metric_{name}": value for name, value in self.metrics.items()})
        return row


def record_run(scenario: Scenario, algorithm_id: str, optimizer: Optimizer, seed: int,
               problem: Any, result: OptimizationResult) -> RunRecord:
    """Build the RunRecord of a finished run, applying the scenario's criteria."""
    gap = None if scenario.reference is None else result.f_best - scenario.reference
    verification = None
    if scenario.verifier is not None and result.found_valid:
        verification = scenario.verifier(problem, result)
    success = None
    reached = None
    if scenario.success_threshold is not None and scenario.reference is not None:
        target = scenario.reference + scenario.success_threshold
        success = bool(result.found_valid and result.f_best <= target
                       and (verification is None or verification["passed"]))
        reached = evaluations_to_target(result.history.best_so_far, target)
    anomaly = bool(gap is not None and scenario.lower_threshold is not None and result.found_valid
                   and gap < -scenario.lower_threshold)
    metrics: dict[str, float] = {}
    if result.best_evaluation is not None:
        for name in scenario.metrics:
            value = result.best_evaluation.metadata.get(name)
            if value is not None and np.ndim(value) == 0:
                metrics[name] = float(value)
    history = np.array(result.history.best_so_far, dtype=float, copy=True)
    history.setflags(write=False)
    return RunRecord(
        scenario=scenario.id, algorithm=algorithm_id, algorithm_name=optimizer.name, seed=int(seed),
        f_best=float(result.f_best), found_valid=bool(result.found_valid), gap=gap, success=success,
        anomaly=anomaly, evaluations_to_success=reached, nfev=int(result.nfev), n_invalid=int(result.n_invalid),
        generations=int(result.generations), termination_reason=result.termination_reason.value,
        wall_time_s=float(result.wall_time_s),
        x_best=None if result.x_best is None else tuple(float(value) for value in result.x_best),
        metrics=metrics, verification=verification, best_so_far=history,
    )


def run_once(scenario: Scenario, algorithm_id: str, optimizer: Optimizer, seed: int) -> RunRecord:
    """Run ``optimizer`` on a fresh problem of ``scenario`` with its full budget and one seed."""
    problem = scenario.problem_factory()
    result = optimizer.minimize(problem, TerminationCriteria(max_evaluations=scenario.max_evaluations), seed=seed)
    return record_run(scenario, algorithm_id, optimizer, seed, problem, result)


def _run_from_config(task: tuple[Mapping[str, Any], Mapping[str, Any], int]) -> RunRecord:
    """Worker entry point: rebuild the scenario and the optimizer from their configurations."""
    scenario_config, algorithm_config, seed = task
    optimizer_config = {key: value for key, value in algorithm_config.items() if key != "id"}
    return run_once(build_scenario(scenario_config), algorithm_config["id"], build_optimizer(optimizer_config), seed)


def validate_config(config: Mapping[str, Any]) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]], list[int]]:
    """Validate a benchmark configuration; return (scenario configs, algorithm configs, seeds).

    Every scenario and optimizer is built once, so an invalid entry is
    reported before any run starts. Invalid configuration raises ValueError.
    """
    if not isinstance(config, Mapping):
        raise ValueError("benchmark configuration must be a mapping")
    required = {"schema", "name", "seeds", "algorithms", "scenarios"}
    unknown = set(config) - required - {"description", "acceptance"}
    if required - set(config) or unknown:
        raise ValueError(f"benchmark configuration: missing keys {sorted(required - set(config))}, "
                         f"unknown keys {sorted(unknown)}")
    if config["schema"] != SCHEMA:
        raise ValueError(f"unsupported schema {config['schema']!r}; expected {SCHEMA!r}")
    seeds_config = config["seeds"]
    if not isinstance(seeds_config, Mapping) or set(seeds_config) != {"first", "count"}:
        raise ValueError("seeds must be {'first': int >= 0, 'count': int >= 1}")
    first, count = seeds_config["first"], seeds_config["count"]
    for value, minimum, name in ((first, 0, "seeds.first"), (count, 1, "seeds.count")):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    algorithms, scenarios = config["algorithms"], config["scenarios"]
    for entries, name in ((algorithms, "algorithms"), (scenarios, "scenarios")):
        if not isinstance(entries, Sequence) or isinstance(entries, str) or not entries:
            raise ValueError(f"{name} must be a nonempty list")
    identifiers = []
    for algorithm in algorithms:
        if not isinstance(algorithm, Mapping) or not isinstance(algorithm.get("id"), str) or not algorithm["id"]:
            raise ValueError("every algorithm needs a nonempty string 'id'")
        identifiers.append(algorithm["id"])
        build_optimizer({key: value for key, value in algorithm.items() if key != "id"})
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("algorithm ids must be unique")
    names = [build_scenario(scenario).id for scenario in scenarios]
    if len(set(names)) != len(names):
        raise ValueError("scenario ids must be unique")
    return list(scenarios), list(algorithms), list(range(first, first + count))


def run_benchmark(config: Mapping[str, Any], *, jobs: int = 1,
                  progress: Callable[[int, int, RunRecord], None] | None = None) -> list[RunRecord]:
    """Execute every (scenario, algorithm, seed) of a validated configuration.

    Records are returned in configuration order (scenario, then algorithm,
    then seed), whatever ``jobs`` is. ``jobs`` >= 1 is the number of worker
    processes; 1 runs in this process. ``progress(done, total, record)`` is
    called after each run, in completion order of the ordered iteration.
    """
    if isinstance(jobs, bool) or not isinstance(jobs, int) or jobs < 1:
        raise ValueError("jobs must be an integer >= 1")
    scenarios, algorithms, seeds = validate_config(config)
    tasks = [(scenario, algorithm, seed) for scenario in scenarios for algorithm in algorithms for seed in seeds]
    records: list[RunRecord] = []
    if jobs == 1:
        iterator: Any = map(_run_from_config, tasks)
        for record in iterator:
            records.append(record)
            if progress is not None:
                progress(len(records), len(tasks), record)
        return records
    with ProcessPoolExecutor(max_workers=jobs) as executor:
        for record in executor.map(_run_from_config, tasks):
            records.append(record)
            if progress is not None:
                progress(len(records), len(tasks), record)
    return records
