"""Multi-seed benchmarking of optimizers on fixed scenarios (Phases 3c and 6A).

Layer: it may import the generic optimizers, the problem adapters and the
physics; nothing in those layers imports it. Statistics use scipy.stats.

statistics: order statistics, Wilson interval, convergence bands, A12,
    Mann-Whitney, Wilcoxon, Friedman and Holm correction.
registry: optimizers and scenarios built from versioned configuration data.
runner: executes every (scenario, algorithm, seed) with equal budgets.
report: aggregates the runs into a statistical summary and checks the
    acceptance criteria declared in the configuration.
plots: figures of the summary (optional Matplotlib).

Entry point: ``python -m experiments.benchmark_optimizers --config ...``.
"""

from .registry import Scenario, build_handling, build_optimizer, build_scenario
from .report import evaluate_acceptance, records_to_rows, summarize_benchmark
from .runner import SCHEMA, RunRecord, record_run, run_benchmark, run_configured, run_once, validate_config
from .statistics import (
    convergence_band, evaluations_to_target, friedman, holm_correction, mann_whitney, quantile,
    summarize, vargha_delaney_a12, wilcoxon_signed_rank, wilson_interval,
)

__all__ = [
    "SCHEMA", "Scenario", "build_optimizer", "build_handling", "build_scenario",
    "RunRecord", "record_run", "run_once", "run_configured", "run_benchmark", "validate_config",
    "summarize_benchmark", "evaluate_acceptance", "records_to_rows",
    "quantile", "summarize", "wilson_interval", "evaluations_to_target", "convergence_band",
    "vargha_delaney_a12", "mann_whitney", "wilcoxon_signed_rank", "friedman", "holm_correction",
]
