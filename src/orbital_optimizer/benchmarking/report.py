"""Aggregate benchmark runs into a reproducible report (RF-6A.01, RF-6A.05).

``summarize_benchmark`` turns the RunRecords of ``runner.run_benchmark`` into
a JSON-compatible summary: per (scenario, algorithm) order statistics of the
final objective and gap, success rate with its Wilson interval, evaluations
to success, validity, termination reasons, wall time and the convergence
band; per scenario, pairwise Mann-Whitney tests with Holm-adjusted p-values
and the A12 effect size; across scenarios, a Friedman test on the median
final objectives when at least three algorithms and two scenarios exist.

Nothing here decides that an algorithm "is better": it reports the evidence
and the class of problems it refers to (the configured scenarios).
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from itertools import combinations
from typing import Any

import numpy as np

from .runner import RunRecord, validate_config
from .statistics import (
    convergence_band, friedman, holm_correction, mann_whitney, summarize, wilcoxon_signed_rank,
    wilson_interval,
)

BAND_QUANTILES: tuple[float, float, float] = (0.25, 0.5, 0.75)


def _cell(records: Sequence[RunRecord], band_points: int) -> dict[str, Any]:
    """Summary of the runs of one algorithm on one scenario."""
    finals = [record.f_best for record in records]
    cell: dict[str, Any] = {
        "algorithm_name": records[0].algorithm_name,
        "runs": len(records),
        "seeds": [record.seed for record in records],
        "found_valid": sum(record.found_valid for record in records),
        "objective": summarize(finals),
        "anomalies": sum(record.anomaly for record in records),
        "termination_reasons": dict(Counter(record.termination_reason for record in records)),
        "nfev": summarize([record.nfev for record in records]),
        "invalid_fraction": summarize([record.n_invalid / record.nfev for record in records]),
        "wall_time_s": summarize([record.wall_time_s for record in records]),
    }
    gaps = [record.gap for record in records if record.gap is not None]
    cell["gap"] = summarize(gaps) if len(gaps) == len(records) else None
    if records[0].success is None:
        cell["success"] = None
        cell["evaluations_to_success"] = None
    else:
        successes = sum(bool(record.success) for record in records)
        low, high = wilson_interval(successes, len(records))
        cell["success"] = {"successes": successes, "runs": len(records), "rate": successes / len(records),
                           "wilson_95": [low, high]}
        reached = [record.evaluations_to_success for record in records
                   if record.success and record.evaluations_to_success is not None]
        cell["evaluations_to_success"] = summarize(reached) if reached else None
    verified = [record.verification for record in records if record.verification is not None]
    cell["verification"] = None if not verified else {
        "checked": len(verified), "passed": sum(bool(item["passed"]) for item in verified)}
    names = sorted({name for record in records for name in record.metrics})
    cell["metrics"] = {name: summarize([record.metrics[name] for record in records if name in record.metrics])
                       for name in names}
    budget = len(records[0].best_so_far)
    if all(len(record.best_so_far) == budget for record in records):
        grid = np.unique(np.round(np.linspace(1, budget, min(band_points, budget))).astype(int))
        band = convergence_band([record.best_so_far for record in records], grid, BAND_QUANTILES)
        cell["convergence"] = {"statistic": "inverted-CDF quantiles across seeds of the best-so-far objective",
                               "evaluations": band["evaluations"],
                               "q25": band["quantiles"][0.25], "median": band["quantiles"][0.5],
                               "q75": band["quantiles"][0.75]}
    else:
        cell["convergence"] = None
    return cell


def summarize_benchmark(config: Mapping[str, Any], records: Sequence[RunRecord], *,
                        band_points: int = 200) -> dict[str, Any]:
    """JSON-compatible statistical summary of a finished benchmark (see module docstring).

    ``records`` must contain every (scenario, algorithm, seed) of ``config``
    exactly once. Non-finite values are kept as floats; serialize with
    optimization.result.jsonable.
    """
    scenarios, algorithms, seeds = validate_config(config)
    scenario_ids = [scenario["id"] for scenario in scenarios]
    algorithm_ids = [algorithm["id"] for algorithm in algorithms]
    grouped: dict[tuple[str, str], list[RunRecord]] = {}
    for record in records:
        grouped.setdefault((record.scenario, record.algorithm), []).append(record)
    expected = {(scenario, algorithm) for scenario in scenario_ids for algorithm in algorithm_ids}
    if set(grouped) != expected or any(sorted(r.seed for r in runs) != seeds for runs in grouped.values()):
        raise ValueError("records do not match the scenarios, algorithms and seeds of the configuration")

    summary: dict[str, Any] = {"scenarios": {}, "comparisons": {}, "overall": None}
    for scenario in scenario_ids:
        summary["scenarios"][scenario] = {algorithm: _cell(grouped[(scenario, algorithm)], band_points)
                                          for algorithm in algorithm_ids}
        pairs = list(combinations(algorithm_ids, 2))
        tests = [mann_whitney([r.f_best for r in grouped[(scenario, first)]],
                              [r.f_best for r in grouped[(scenario, second)]]) for first, second in pairs]
        adjusted = holm_correction([test["p_value"] for test in tests]) if tests else []
        summary["comparisons"][scenario] = [
            {"first": first, "second": second, **test, "p_value_holm": p_holm,
             "family": f"{len(pairs)} pairwise comparisons within the scenario"}
            for (first, second), test, p_holm in zip(pairs, tests, adjusted)]

    medians = np.array([[summary["scenarios"][scenario][algorithm]["objective"]["median"]
                         for algorithm in algorithm_ids] for scenario in scenario_ids])
    # Only the order within a scenario matters for rank tests: +inf becomes "worse than every finite median".
    ranked = medians.copy()
    for row in ranked:
        finite = row[np.isfinite(row)]
        row[~np.isfinite(row)] = (finite.max() if finite.size else 0.0) + 1.0
    if len(algorithm_ids) >= 3 and len(scenario_ids) >= 2:
        summary["overall"] = {"scores": "median final objective per scenario", "algorithms": algorithm_ids,
                              "scenarios": scenario_ids, **friedman(ranked)}
    elif len(algorithm_ids) == 2 and len(scenario_ids) >= 2:
        summary["overall"] = {"scores": "median final objective per scenario", "algorithms": algorithm_ids,
                              "scenarios": scenario_ids, **wilcoxon_signed_rank(ranked[:, 0], ranked[:, 1])}
    return summary


def evaluate_acceptance(acceptance: Mapping[str, Any] | None, summary: Mapping[str, Any]) -> dict[str, Any]:
    """Check a benchmark summary against the acceptance criteria declared in its configuration.

    Supported keys (all optional; unknown keys raise ValueError):
      - ``max_anomalies``: largest admitted number of runs better than the
        reference by more than the scenario's lower threshold, over all cells;
      - ``require_verification``: every verified best plan must have passed;
      - ``minimum_success_rate``: {scenario: {algorithm: rate}}.
    Returns {"checks": {name: bool}, "details": {...}, "passed": bool}; with
    no criteria it passes vacuously and says so.
    """
    acceptance = dict(acceptance or {})
    unknown = set(acceptance) - {"max_anomalies", "require_verification", "minimum_success_rate"}
    if unknown:
        raise ValueError(f"unknown acceptance keys {sorted(unknown)}")
    cells = {(scenario, algorithm): cell for scenario, row in summary["scenarios"].items()
             for algorithm, cell in row.items()}
    checks: dict[str, bool] = {}
    details: dict[str, Any] = {}
    if "max_anomalies" in acceptance:
        anomalies = sum(cell["anomalies"] for cell in cells.values())
        checks["anomalies"] = anomalies <= acceptance["max_anomalies"]
        details["anomalies"] = {"observed": anomalies, "admitted": acceptance["max_anomalies"]}
    if acceptance.get("require_verification"):
        failed = {f"{scenario}/{algorithm}": cell["verification"] for (scenario, algorithm), cell in cells.items()
                  if cell["verification"] is not None and cell["verification"]["passed"] != cell["verification"]["checked"]}
        checks["verification"] = not failed
        details["verification"] = {"cells_with_failures": failed}
    for scenario, row in dict(acceptance.get("minimum_success_rate", {})).items():
        for algorithm, minimum in dict(row).items():
            cell = cells.get((scenario, algorithm))
            if cell is None or cell["success"] is None:
                raise ValueError(f"minimum_success_rate refers to {scenario}/{algorithm}, which has no success criterion")
            name = f"success_rate:{scenario}/{algorithm}"
            checks[name] = cell["success"]["rate"] >= minimum
            details[name] = {"observed": cell["success"]["rate"], "minimum": minimum,
                             "wilson_95": cell["success"]["wilson_95"]}
    return {"checks": checks, "details": details, "passed": all(checks.values()),
            "note": "no acceptance criteria declared" if not checks else ""}


def records_to_rows(records: Sequence[RunRecord]) -> list[dict[str, Any]]:
    """Flat rows (one per run) with a stable union of columns, for CSV export."""
    rows = [record.to_row() for record in records]
    columns: list[str] = []
    for row in rows:
        columns.extend(name for name in row if name not in columns)
    return [{name: row.get(name) for name in columns} for row in rows]


def format_number(value: Any) -> str:
    """Compact text for report tables: 4 significant digits, 'inf'/'nan' spelled out."""
    if value is None:
        return "-"
    if isinstance(value, (bool, str)):
        return str(value)
    number = float(value)
    if math.isnan(number):
        return "nan"
    if math.isinf(number):
        return "inf" if number > 0 else "-inf"
    return f"{number:.4g}"
