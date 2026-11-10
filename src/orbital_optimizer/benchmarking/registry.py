"""Build optimizers and benchmark scenarios from versioned configuration data (RR-02).

A benchmark is defined by JSON-compatible dictionaries, never by parameters
embedded in code. This module turns them into objects; unknown keys and
unknown types are errors, so a typo cannot silently fall back to a default.

Optimizer configuration: {"type": <name>, "settings": {...}} with types
``differential_evolution``, ``particle_swarm``, ``genetic_algorithm``,
``random_search``, ``nelder_mead`` and ``hybrid`` (settings: ``global``, an
optimizer configuration, and ``global_fraction``).

Scenario configuration: {"id", "type", "max_evaluations", ...}:
  - ``function``: ``function`` in {sphere, rosenbrock, rastrigin},
    ``dimension``, ``lower``, ``upper`` and ``success_threshold`` (absolute
    gap f - f*; the known minimum is the reference).
  - ``coplanar_transfer``: ``initial_radius_m``, ``final_radius_m``,
    ``impulses`` (default 2), optional ``mu_m3_s2``, ``success_threshold_m_s``
    and ``lower_threshold_m_s`` (default 1e-6). The reference is the Hohmann
    total delta-v; with ``reference`` = {"type": "bielliptic",
    "intermediate_radius_m": rb} it is the bi-elliptic one. Every best plan
    is verified by numerical propagation with the options in
    ``verification`` (radius_tolerance_m, eccentricity_tolerance, rtol, atol,
    max_step_s).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any

from ..constants import EARTH_MU
from ..constraints.limits import evaluate_constraints, limits_from_mapping
from ..constraints.verification import numerical_transfer_summary
from ..optimization.benchmarks import rastrigin_benchmark, rosenbrock_benchmark, sphere_benchmark
from ..optimization.constraint_handling import handling_from_mapping
from ..optimization.differential_evolution import DifferentialEvolution, DifferentialEvolutionSettings
from ..optimization.genetic_algorithm import GeneticAlgorithm, GeneticAlgorithmSettings
from ..optimization.hybrid import HybridOptimizer, HybridSettings
from ..optimization.nelder_mead import NelderMead, NelderMeadSettings
from ..optimization.optimizer import Optimizer
from ..optimization.particle_swarm import ParticleSwarm, ParticleSwarmSettings
from ..optimization.problem import Problem
from ..optimization.random_search import RandomSearch, RandomSearchSettings
from ..optimization.result import OptimizationResult
from ..orbital.transfers import bielliptic_transfer, hohmann_transfer
from ..problems.constrained import ConstrainedTransferProblem
from ..problems.coplanar_multi_impulse import CoplanarMultiImpulseProblem, verify_impulsive_plan
from ..problems.references import time_limited_reference
from ..spacecraft import spacecraft_from_mapping

OPTIMIZERS: dict[str, tuple[type, type]] = {
    "differential_evolution": (DifferentialEvolution, DifferentialEvolutionSettings),
    "particle_swarm": (ParticleSwarm, ParticleSwarmSettings),
    "genetic_algorithm": (GeneticAlgorithm, GeneticAlgorithmSettings),
    "random_search": (RandomSearch, RandomSearchSettings),
    "nelder_mead": (NelderMead, NelderMeadSettings),
}
FUNCTIONS = {"sphere": sphere_benchmark, "rosenbrock": rosenbrock_benchmark, "rastrigin": rastrigin_benchmark}
TRANSFER_METRICS: tuple[str, ...] = (
    "delta_v_total_m_s", "transfer_time_s", "minimum_radius_m", "maximum_radius_m",
    "impulse_kinetic_energy_j_kg", "closest_approach_index",
)


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _keys(config: Mapping[str, Any], name: str, required: set[str], optional: set[str] = frozenset()) -> None:
    missing = required - set(config)
    unknown = set(config) - required - optional
    if missing or unknown:
        raise ValueError(f"{name}: missing keys {sorted(missing)}, unknown keys {sorted(unknown)}")


def _real(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return float(value)


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be an integer >= 1")
    return value


def build_handling(config: Mapping[str, Any]) -> Any:
    """Constraint handling of an algorithm entry, or None when it does not declare one.

    The entry may carry ``constraint_handling`` = {"type": "feasibility_rules" |
    "death_penalty" | "static_penalty", ...}. The runner then sets it on the
    problem of each run, so that strategies are compared on the same problem
    with the same optimizer and seeds (RF-4.10).
    """
    handling = _mapping(config, "algorithm configuration").get("constraint_handling")
    return None if handling is None else handling_from_mapping(handling)


def build_optimizer(config: Mapping[str, Any]) -> Optimizer:
    """Optimizer from {"type": ..., "settings": {...}}; invalid configuration raises ValueError.

    An optional ``constraint_handling`` entry is validated here and applied by the runner.
    """
    config = _mapping(config, "optimizer configuration")
    _keys(config, "optimizer configuration", {"type"}, {"settings", "constraint_handling"})
    build_handling(config)
    kind = config["type"]
    settings = dict(_mapping(config.get("settings", {}), "settings"))
    if kind == "hybrid":
        _keys(settings, "hybrid settings", set(), {"global", "global_fraction"})
        inner = build_optimizer(settings.get("global", {"type": "differential_evolution"}))
        fraction = settings.get("global_fraction", HybridSettings().global_fraction)
        return HybridOptimizer(inner, HybridSettings(global_fraction=fraction))
    if kind not in OPTIMIZERS:
        raise ValueError(f"unknown optimizer type {kind!r}; expected one of {sorted(OPTIMIZERS) + ['hybrid']}")
    optimizer_class, settings_class = OPTIMIZERS[kind]
    if "initial_point" in settings and settings["initial_point"] is not None:
        settings["initial_point"] = tuple(settings["initial_point"])
    try:
        return optimizer_class(settings_class(**settings))
    except TypeError as exc:
        raise ValueError(f"invalid settings for {kind}: {exc}") from exc


@dataclass(frozen=True)
class Scenario:
    """A fixed benchmark problem with its budget, reference and success criterion.

    problem_factory: builds a fresh Problem for every run. reference: known
    optimum or analytic reference of the objective, or None. success_threshold:
    a run succeeds iff it found a reportable point (valid, and feasible for a
    constrained problem) with f_best - reference <= success_threshold (and,
    when ``verifier`` exists, its verification passed). lower_threshold: a
    run with f_best - reference < -lower_threshold is an anomaly (better than
    a proven optimum). objective_unit labels the objective.
    verifier(problem, result) returns a JSON-compatible record with a boolean
    ``passed``. metrics: metadata keys copied from the best evaluation.
    expect_infeasible: the scenario has no feasible solution by construction;
    a run then succeeds iff it reports that it found none.
    """

    id: str
    description: str
    problem_factory: Callable[[], Problem]
    max_evaluations: int
    objective_unit: str
    reference: float | None = None
    success_threshold: float | None = None
    lower_threshold: float | None = None
    verifier: Callable[[Problem, OptimizationResult], dict[str, Any]] | None = None
    metrics: tuple[str, ...] = ()
    expect_infeasible: bool = False


def _function_scenario(config: Mapping[str, Any]) -> Scenario:
    _keys(config, "function scenario", {"id", "type", "function", "dimension", "max_evaluations"},
          {"lower", "upper", "success_threshold", "description"})
    name = config["function"]
    if name not in FUNCTIONS:
        raise ValueError(f"unknown function {name!r}; expected one of {sorted(FUNCTIONS)}")
    arguments = {key: _real(config[key], key) for key in ("lower", "upper") if key in config}
    benchmark = FUNCTIONS[name](_positive_int(config["dimension"], "dimension"), **arguments)
    threshold = config.get("success_threshold")
    return Scenario(
        id=config["id"],
        description=config.get("description", f"{benchmark.name} on [{benchmark.bounds.lower[0]:g}, "
                                              f"{benchmark.bounds.upper[0]:g}]^D, known minimum {benchmark.minimum:g}"),
        problem_factory=benchmark.problem,
        max_evaluations=_positive_int(config["max_evaluations"], "max_evaluations"),
        objective_unit="dimensionless",
        reference=benchmark.minimum,
        success_threshold=None if threshold is None else _real(threshold, "success_threshold"),
        lower_threshold=0.0,
    )


def _transfer_scenario(config: Mapping[str, Any]) -> Scenario:
    _keys(config, "coplanar_transfer scenario",
          {"id", "type", "initial_radius_m", "final_radius_m", "max_evaluations"},
          {"impulses", "mu_m3_s2", "success_threshold_m_s", "lower_threshold_m_s", "reference",
           "verification", "description"})
    r1, r2 = _real(config["initial_radius_m"], "initial_radius_m"), _real(config["final_radius_m"], "final_radius_m")
    mu = _real(config.get("mu_m3_s2", EARTH_MU), "mu_m3_s2")
    impulses = _positive_int(config.get("impulses", 2), "impulses")
    CoplanarMultiImpulseProblem(r1, r2, mu, impulses=impulses)          # validate now, not in a worker
    reference_config = dict(_mapping(config.get("reference", {"type": "hohmann"}), "reference"))
    kind = reference_config.get("type")
    if kind == "hohmann":
        _keys(reference_config, "reference", {"type"})
        reference, label = hohmann_transfer(r1, r2, mu).total_delta_v, "Hohmann"
    elif kind == "bielliptic":
        _keys(reference_config, "reference", {"type", "intermediate_radius_m"})
        rb = _real(reference_config["intermediate_radius_m"], "intermediate_radius_m")
        reference, label = bielliptic_transfer(r1, r2, rb, mu).total_delta_v, f"bi-elliptic through {rb:g} m"
    else:
        raise ValueError(f"unknown reference type {kind!r}; expected 'hohmann' or 'bielliptic'")
    verification = dict(_mapping(config.get("verification", {}), "verification"))
    _keys(verification, "verification", set(),
          {"radius_tolerance_m", "eccentricity_tolerance", "rtol", "atol", "max_step_s"})
    options: dict[str, Any] = {key: _real(verification[key], key)
                               for key in ("radius_tolerance_m", "eccentricity_tolerance", "rtol") if key in verification}
    if "max_step_s" in verification:
        options["max_step"] = _real(verification["max_step_s"], "max_step_s")
    if "atol" in verification:
        options["atol"] = verification["atol"]

    def factory() -> Problem:
        return CoplanarMultiImpulseProblem(r1, r2, mu, impulses=impulses)

    def verifier(problem: Problem, result: OptimizationResult) -> dict[str, Any]:
        outcome = verify_impulsive_plan(problem.plan(result.x_best), **options)   # type: ignore[attr-defined]
        return {"passed": outcome.passed, "message": outcome.message, "checks": dict(outcome.checks),
                "metrics": dict(outcome.metrics), "tolerances": dict(outcome.tolerances),
                "captured_value_errors": problem.captured_value_errors,           # type: ignore[attr-defined]
                "captured_runtime_errors": problem.captured_runtime_errors}       # type: ignore[attr-defined]

    threshold = config.get("success_threshold_m_s")
    return Scenario(
        id=config["id"],
        description=config.get("description", f"coplanar {impulses}-impulse transfer {r1:g} m -> {r2:g} m; "
                                              f"reference: {label}"),
        problem_factory=factory,
        max_evaluations=_positive_int(config["max_evaluations"], "max_evaluations"),
        objective_unit="m/s",
        reference=reference,
        success_threshold=None if threshold is None else _real(threshold, "success_threshold_m_s"),
        lower_threshold=_real(config.get("lower_threshold_m_s", 1e-6), "lower_threshold_m_s"),
        verifier=verifier,
        metrics=TRANSFER_METRICS,
    )


CONSTRAINED_METRICS: tuple[str, ...] = TRANSFER_METRICS + ("propellant_required_kg", "final_mass_kg")
# Tolerance of the numerical re-evaluation of the constraints, relative to each constraint's scale.
# The integrator's global error over a LEO-GEO transfer is about 1e-8 of the radius
# (0.06 m of 4.2e7 m, docs/phase1-phase2-audit.md); 1e-6 leaves two orders of margin.
VERIFICATION_RELATIVE_TOLERANCE: float = 1e-6


def _constrained_scenario(config: Mapping[str, Any]) -> Scenario:
    _keys(config, "constrained_transfer scenario",
          {"id", "type", "initial_radius_m", "final_radius_m", "max_evaluations", "spacecraft", "limits", "reference"},
          {"impulses", "mu_m3_s2", "commanded_thrust_n", "success_threshold_m_s", "lower_threshold_m_s",
           "verification", "description"})
    r1, r2 = _real(config["initial_radius_m"], "initial_radius_m"), _real(config["final_radius_m"], "final_radius_m")
    mu = _real(config.get("mu_m3_s2", EARTH_MU), "mu_m3_s2")
    impulses = _positive_int(config.get("impulses", 2), "impulses")
    spacecraft = spacecraft_from_mapping(config["spacecraft"])
    limits = limits_from_mapping(config["limits"])
    thrust = config.get("commanded_thrust_n")
    thrust = None if thrust is None else _real(thrust, "commanded_thrust_n")
    ConstrainedTransferProblem(CoplanarMultiImpulseProblem(r1, r2, mu, impulses=impulses), spacecraft, limits,
                               commanded_thrust_n=thrust)               # validate now, not in a worker

    reference_config = dict(_mapping(config["reference"], "reference"))
    _keys(reference_config, "reference", {"type"})
    kind = reference_config["type"]
    expect_infeasible = False
    if kind == "hohmann":
        reference, label = hohmann_transfer(r1, r2, mu).total_delta_v, "Hohmann (the constraints do not bind at it)"
    elif kind == "time_limited":
        if limits.max_duration_s is None or impulses != 2 or r2 <= r1:
            raise ValueError("the time_limited reference needs max_duration_s and an ascending two-impulse transfer")
        solution = time_limited_reference(r1, r2, limits.max_duration_s, mu)
        reference = solution.delta_v_total_m_s
        label = f"minimum under the time limit ({solution.method})"
    elif kind == "infeasible":
        reference, label, expect_infeasible = None, "none: no feasible transfer exists by construction", True
    else:
        raise ValueError(f"unknown reference type {kind!r}; expected 'hohmann', 'time_limited' or 'infeasible'")

    verification = dict(_mapping(config.get("verification", {}), "verification"))
    _keys(verification, "verification", set(),
          {"radius_tolerance_m", "eccentricity_tolerance", "rtol", "atol", "max_step_s", "constraint_relative_tolerance"})
    relative = _real(verification.get("constraint_relative_tolerance", VERIFICATION_RELATIVE_TOLERANCE),
                     "constraint_relative_tolerance")
    arrival: dict[str, Any] = {key: _real(verification[key], key)
                               for key in ("radius_tolerance_m", "eccentricity_tolerance", "rtol") if key in verification}
    integrator: dict[str, Any] = {"rtol": arrival["rtol"]} if "rtol" in arrival else {}
    if "max_step_s" in verification:
        arrival["max_step"] = integrator["max_step"] = _real(verification["max_step_s"], "max_step_s")
    if "atol" in verification:
        arrival["atol"] = integrator["atol"] = verification["atol"]
    # Same limits, with the tolerance of the numerical re-evaluation instead of the analytic one.
    numerical_limits = limits_from_mapping({**dict(config["limits"]), "relative_tolerance": relative})

    def factory() -> Problem:
        return ConstrainedTransferProblem(CoplanarMultiImpulseProblem(r1, r2, mu, impulses=impulses),
                                          spacecraft, limits, commanded_thrust_n=thrust)

    def verifier(problem: Problem, result: OptimizationResult) -> dict[str, Any]:
        plan = problem.plan(result.x_best)                                        # type: ignore[attr-defined]
        arrived = verify_impulsive_plan(plan, **arrival)
        summary = numerical_transfer_summary(plan.initial_state, plan.impulses(), plan.mu, thrust_n=thrust,
                                             **integrator)
        recomputed = evaluate_constraints(summary, spacecraft, numerical_limits)
        return {
            "passed": bool(arrived.passed and recomputed.feasible),
            "message": arrived.message if recomputed.feasible
            else f"{arrived.message}; constraints violated on numerical re-evaluation: {list(recomputed.violated)}",
            "arrival_passed": arrived.passed, "checks": dict(arrived.checks), "metrics": dict(arrived.metrics),
            "tolerances": dict(arrived.tolerances),
            "constraints_feasible": recomputed.feasible, "constraints_violated": list(recomputed.violated),
            "constraint_relative_tolerance": relative,
            "numerical_minimum_radius_m": summary.minimum_radius_m,
            "captured_value_errors": problem.captured_value_errors,               # type: ignore[attr-defined]
            "captured_runtime_errors": problem.captured_runtime_errors,           # type: ignore[attr-defined]
        }

    threshold = config.get("success_threshold_m_s")
    return Scenario(
        id=config["id"],
        description=config.get("description", f"constrained coplanar {impulses}-impulse transfer "
                                              f"{r1:g} m -> {r2:g} m; reference: {label}"),
        problem_factory=factory,
        max_evaluations=_positive_int(config["max_evaluations"], "max_evaluations"),
        objective_unit="m/s",
        reference=reference,
        success_threshold=None if threshold is None else _real(threshold, "success_threshold_m_s"),
        lower_threshold=_real(config.get("lower_threshold_m_s", 1e-6), "lower_threshold_m_s"),
        verifier=verifier,
        metrics=CONSTRAINED_METRICS,
        expect_infeasible=expect_infeasible,
    )


SCENARIO_BUILDERS: dict[str, Callable[[Mapping[str, Any]], Scenario]] = {
    "function": _function_scenario,
    "coplanar_transfer": _transfer_scenario,
    "constrained_transfer": _constrained_scenario,
}


def build_scenario(config: Mapping[str, Any]) -> Scenario:
    """Scenario from its configuration; invalid configuration raises ValueError."""
    config = _mapping(config, "scenario configuration")
    kind = config.get("type")
    if kind not in SCENARIO_BUILDERS:
        raise ValueError(f"unknown scenario type {kind!r}; expected one of {sorted(SCENARIO_BUILDERS)}")
    if not isinstance(config.get("id"), str) or not config["id"]:
        raise ValueError("scenario id must be a nonempty string")
    return SCENARIO_BUILDERS[kind](config)
