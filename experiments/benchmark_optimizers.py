"""Reproducible multi-seed benchmark of optimizers on fixed scenarios (RR-06).

    python -m experiments.benchmark_optimizers --config configs/benchmarks/<name>.json \
        --output outputs/<new-folder> [--plots] [--jobs N]

Reads a versioned configuration (schema ``daedalus.benchmark/1``), runs every
(scenario, algorithm, seed) with the scenario's evaluation budget, repeats the
first seed of every pair to check bitwise reproducibility (RV-3.03), verifies
the best plan of every transfer run by numerical propagation (RV-3.04) and
writes, in SI units:

    report.json       configuration, environment, code version, statistical
                      summary, reproducibility check and acceptance result
    runs.csv          one row per run
    convergence.csv   quantiles across seeds of the best-so-far objective
    *.png             with --plots: convergence, final distribution, success
                      rate, cost against wall time and the best transfer of
                      each transfer scenario

Exit code 0 if the acceptance criteria declared in the configuration hold and
every repeated run was identical, 1 if not, 2 on a configuration or runtime
error.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np

from orbital_optimizer.benchmarking import (
    RunRecord, build_optimizer, build_scenario, evaluate_acceptance, records_to_rows, run_benchmark,
    run_once, summarize_benchmark, validate_config,
)
from orbital_optimizer.constants import EARTH_EQUATORIAL_RADIUS, EARTH_MU
from orbital_optimizer.optimization.result import jsonable
from orbital_optimizer.provenance import code_version, environment_record
from orbital_optimizer.simulation.impulsive import propagate_impulsive_arcs

LIMITATIONS = [
    "Results hold for the configured scenarios, budgets and settings; they are not general rankings of algorithms.",
    "Hyperparameters are literature defaults for every algorithm; none was tuned on these scenarios.",
    "Wall time depends on the machine, the interpreter and the number of worker processes; it is context only.",
    "Transfer scenarios use the two-body impulsive model: no perturbations, finite burns or mission constraints.",
    "A success is a run within the declared threshold of the reference; it is not a proof of global optimality.",
]


def _same_run(first: RunRecord, second: RunRecord) -> bool:
    """True if two records of the same (scenario, algorithm, seed) are bitwise identical except wall time."""
    left, right = first.to_row(), second.to_row()
    left.pop("wall_time_s")
    right.pop("wall_time_s")
    return left == right and first.best_so_far.tobytes() == second.best_so_far.tobytes()


def check_reproducibility(config: dict[str, Any], records: list[RunRecord]) -> dict[str, Any]:
    """Repeat the first seed of every (scenario, algorithm) in this process and compare with the recorded run."""
    scenarios, algorithms, seeds = validate_config(config)
    reference = {(record.scenario, record.algorithm, record.seed): record for record in records}
    different = []
    for scenario_config in scenarios:
        scenario = build_scenario(scenario_config)
        for algorithm in algorithms:
            optimizer = build_optimizer({key: value for key, value in algorithm.items() if key != "id"})
            repeated = run_once(scenario, algorithm["id"], optimizer, seeds[0])
            if not _same_run(repeated, reference[(scenario.id, algorithm["id"], seeds[0])]):
                different.append(f"{scenario.id}/{algorithm['id']}")
    checked = len(scenarios) * len(algorithms)
    return {"seed": seeds[0], "checked": checked, "identical": checked - len(different), "different": different,
            "passed": not different,
            "method": "second execution of the same configuration in the parent process; records compared "
                      "field by field and best-so-far histories byte by byte, wall time excluded"}


def scenario_table(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Reference, unit, budget and criteria of every scenario, as resolved from the configuration."""
    table = {}
    for scenario_config in config["scenarios"]:
        scenario = build_scenario(scenario_config)
        table[scenario.id] = {
            "description": scenario.description, "objective_unit": scenario.objective_unit,
            "max_evaluations": scenario.max_evaluations, "reference": scenario.reference,
            "success_threshold": scenario.success_threshold, "lower_threshold": scenario.lower_threshold,
            "verified_by_propagation": scenario.verifier is not None,
        }
    return table


def write_runs(rows: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({name: json.dumps(value) if isinstance(value, list) else value
                             for name, value in row.items()})


def write_convergence(summary: dict[str, Any], scenarios: dict[str, dict[str, Any]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scenario", "algorithm", "objective_unit", "evaluations", "best_so_far_q25",
                         "best_so_far_median", "best_so_far_q75"])
        for scenario, cells in summary["scenarios"].items():
            for algorithm, cell in cells.items():
                band = cell["convergence"]
                if band is None:
                    continue
                for row in zip(band["evaluations"], band["q25"], band["median"], band["q75"]):
                    writer.writerow([scenario, algorithm, scenarios[scenario]["objective_unit"], *row])


def save_best_transfers(config: dict[str, Any], records: list[RunRecord], output: Path) -> list[str]:
    """One figure per transfer scenario: the verified run with the lowest objective."""
    from orbital_optimizer.visualization.transfer_plots import save_transfer_figure

    written = []
    for scenario_config in config["scenarios"]:
        scenario = build_scenario(scenario_config)
        problem = scenario.problem_factory()
        if not hasattr(problem, "plan"):
            continue
        candidates = [record for record in records if record.scenario == scenario.id and record.found_valid
                      and record.verification is not None and record.verification["passed"]]
        if not candidates:
            continue
        best = min(candidates, key=lambda record: record.f_best)
        plan = problem.plan(np.array(best.x_best))
        trajectory = propagate_impulsive_arcs(plan.initial_state, plan.impulses(), plan.transfer_time, plan.mu,
                                              method="Kepler", samples_per_arc=361)
        path = output / f"transfer_{scenario.id}.png"
        save_transfer_figure(
            trajectory, problem.initial_radius, problem.final_radius, path,
            body_radius=EARTH_EQUATORIAL_RADIUS if plan.mu == EARTH_MU else None,
            title=f"{scenario.id}: best verified run ({best.algorithm}, seed {best.seed}), "
                  f"total Δv {best.f_best:.3f} m/s",
        )
        written.append(path.name)
    return written


def run(config: dict[str, Any], output: Path, *, jobs: int = 1, plots: bool = False,
        quiet: bool = False) -> dict[str, Any]:
    """Execute the benchmark and write its artifacts to ``output``; return the report."""
    validate_config(config)
    scenarios = scenario_table(config)
    started = perf_counter()

    def progress(done: int, total: int, record: RunRecord) -> None:
        if not quiet and (done == total or done % 25 == 0):
            print(f"  {done}/{total} runs ({record.scenario}/{record.algorithm})", file=sys.stderr, flush=True)

    records = run_benchmark(config, jobs=jobs, progress=progress)
    reproducibility = check_reproducibility(config, records)
    summary = summarize_benchmark(config, records)
    acceptance = evaluate_acceptance(config.get("acceptance"), summary)
    elapsed = perf_counter() - started
    passed = bool(acceptance["passed"] and reproducibility["passed"])

    output.mkdir(parents=True, exist_ok=True)
    rows = records_to_rows(records)
    write_runs(rows, output / "runs.csv")
    write_convergence(summary, scenarios, output / "convergence.csv")
    figures: list[str] = []
    if plots:
        from orbital_optimizer.benchmarking.plots import (
            save_convergence_figure, save_cost_figure, save_final_distribution_figure, save_success_figure,
        )
        references = {name: info["reference"] for name, info in scenarios.items()}
        units = {name: info["objective_unit"] for name, info in scenarios.items()}
        figures.append(save_convergence_figure(summary, references, units, output / "convergence.png").name)
        figures.append(save_final_distribution_figure(rows, units, output / "final_distribution.png").name)
        success = save_success_figure(summary, output / "success_rate.png")
        if success is not None:
            figures.append(success.name)
        figures.append(save_cost_figure(summary, units, output / "cost.png").name)
        figures.extend(save_best_transfers(config, records, output))

    report = {
        "experiment": "optimizer_benchmark", "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": config,
        "environment": environment_record(),
        "code": code_version(),
        "execution": {"jobs": jobs, "runs": len(records), "total_wall_time_s": elapsed},
        "randomness": "one Generator(PCG64(seed)) per run; seeds listed in the configuration",
        "scenarios": scenarios,
        "summary": summary,
        "reproducibility": reproducibility,
        "acceptance": acceptance,
        "passed": passed,
        "artifacts": ["report.json", "runs.csv", "convergence.csv", *figures],
        "limitations": LIMITATIONS,
    }
    (output / "report.json").write_text(json.dumps(jsonable(report), indent=2, allow_nan=False) + "\n",
                                        encoding="utf-8")
    return report


def headline(report: dict[str, Any]) -> dict[str, Any]:
    """Short console summary: per scenario and algorithm, median gap and success."""
    table: dict[str, Any] = {}
    for scenario, cells in report["summary"]["scenarios"].items():
        table[scenario] = {}
        for algorithm, cell in cells.items():
            statistics = cell["gap"] if cell["gap"] is not None else cell["objective"]
            entry = {"median": statistics["median"], "worst": statistics["worst"]}
            if cell["success"] is not None:
                entry["success"] = f"{cell['success']['successes']}/{cell['success']['runs']}"
            table[scenario][algorithm] = entry
    return table


def main() -> int:
    """CLI: see the module docstring."""
    parser = argparse.ArgumentParser(description="Multi-seed optimizer benchmark on fixed scenarios.")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1, help="worker processes (default 1)")
    parser.add_argument("--plots", action="store_true", help="export PNG figures with Matplotlib")
    parser.add_argument("--quiet", action="store_true", help="do not print progress")
    args = parser.parse_args()
    try:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        report = run(config, args.output, jobs=args.jobs, plots=args.plots, quiet=args.quiet)
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        parser.exit(2, f"Benchmark failed: {exc}\n")
    print(json.dumps(jsonable({
        "passed": report["passed"], "acceptance": report["acceptance"]["checks"],
        "reproducibility": {key: report["reproducibility"][key] for key in ("checked", "identical")},
        "results": headline(report), "report": str((args.output / "report.json").resolve()),
    }), indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
