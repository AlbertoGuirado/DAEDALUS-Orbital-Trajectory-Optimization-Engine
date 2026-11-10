"""Numerical evidence behind the tolerances and remarks of the Phase 3 documentation (RNF-10).

    python -m experiments.phase3_numerical_evidence --output outputs/<new-folder> [--quick]

Every figure quoted in docs/phase3-api.md and docs/phase3-validation.md that
does not come from the benchmark comes from this script, so that it can be
reproduced from the repository. Sections of ``report.json``:

    two_impulse_tangency   reach margin Q of the float Hohmann decision over a
                           grid of mu, scale and radius ratio (tolerance K = 32)
    band_edge              admitted infeasible band of the two-impulse problem,
                           ascending and descending, down to nearly equal radii
    arc_tolerance          reach margin of exact bi-elliptic decisions after a
                           coast arc, on a grid and on a seeded random sample:
                           the constant the arc term needs, and how many
                           decisions the problem accepts, excludes or rejects
    arrival_conditioning   arrival residual against the conditioning of the
                           arcs, with the conditioning cap and the arrival
                           guard lifted: the measurement behind the cap
    near_parabolic         decisions next to the escape bound: no VALID
                           evaluation may miss the target radius
    coincident_split       circularization split at the apoapsis (DT-13)
    seam                   Differential Evolution on a descending transfer with
                           the branch cut of the angle at the optimum and
                           opposite to it, 30 seeds each
    long_transfer          verification of a 13-day bi-elliptic transfer with
                           the default and with a tightened integrator
    rosenbrock             DE on Rosenbrock-5D: final gap, population spread
                           and last improving generation
    evaluation_cost        wall time per evaluation (machine dependent)

Deterministic except for ``evaluation_cost``; random samples use fixed seeds.
``--quick`` shrinks the grids and budgets for a smoke run (not for evidence).
Exit code 0 if every check in ``checks`` holds, 1 otherwise, 2 on error.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np

import orbital_optimizer.problems.coplanar_multi_impulse as multi_impulse
from orbital_optimizer import EARTH_MU, bielliptic_transfer, hohmann_transfer
from orbital_optimizer.benchmarking.statistics import mann_whitney, summarize, wilson_interval
from orbital_optimizer.optimization import (
    Bounds, DifferentialEvolution, EvaluationStatus, TerminationCriteria, rosenbrock_benchmark,
)
from orbital_optimizer.optimization.result import jsonable
from orbital_optimizer.problems import (
    CoplanarMultiImpulseProblem, CoplanarTwoImpulseProblem, reach_margin, verify_impulsive_plan,
)
from orbital_optimizer.problems.coplanar_multi_impulse import (
    ARC_TANGENCY_FACTOR, ARRIVAL_RELATIVE_TOLERANCE, MAX_ARC_CONDITIONING, arc_conditioning,
)
from orbital_optimizer.problems.coplanar_transfer import TANGENCY_FACTOR, UNIT_ROUNDOFF
from orbital_optimizer.provenance import code_version, environment_record

LEO, GEO = 6_678_000.0, 42_164_000.0
SUN_MU = 1.32712440018e20       # heliocentric gravitational parameter, m^3/s^2 (IAU 2015 nominal value)
MOON_MU = 4.9048695e12          # lunar gravitational parameter, m^3/s^2 (DE421 value)
SEAM_SUCCESS_M_S = 1.0          # success threshold of the benchmark scenario geo-leo


def two_impulse_tangency(quick: bool) -> dict[str, Any]:
    """Q of the float Hohmann first impulse: |Q|/(u S) and how often Q is negative."""
    ratios = np.geomspace(1.0001, 100.0, 8 if quick else 40)
    worst, negative, rejected, cases = 0.0, 0, 0, 0
    for mu in (EARTH_MU, SUN_MU, MOON_MU):
        for radius in (7.0e6, 7.0e8, 1.5e11):
            for ratio in ratios:
                final = radius * float(ratio)
                speed = math.sqrt(mu / radius) + hohmann_transfer(radius, final, mu).first_delta_v
                q, s = reach_margin(radius, final, mu, speed, 0.0)
                cases += 1
                worst = max(worst, abs(q) / (UNIT_ROUNDOFF * s))
                negative += q < 0
                rejected += q < -TANGENCY_FACTOR * UNIT_ROUNDOFF * s
    return {"grid": "mu in (Earth, Sun, Moon) x r1 in (7e6, 7e8, 1.5e11) m x ratios geomspace(1.0001, 100)",
            "cases": cases, "max_abs_q_over_u_s": worst, "q_negative": int(negative),
            "hohmann_rejected": int(rejected), "tolerance_factor": TANGENCY_FACTOR}


def band_edge(quick: bool) -> dict[str, Any]:
    """Smallest accepted first impulse below the float Hohmann value, by bisection, tangential burn.

    For each ratio r2/r1 the two-impulse problem is evaluated along the
    Hohmann burn direction with a first impulse slightly too small to reach
    the target. The last accepted value marks the edge of the admitted band;
    its objective minus the Hohmann optimum is the largest advantage the band
    can give, which an anomaly threshold has to exceed.
    """
    excess = (GEO / LEO - 1, 1e-2) if quick else (GEO / LEO - 1, 1e-2, 1e-3, 1e-4)
    rows = []
    for ratio in [1 + value for value in excess] + [1 / (1 + value) for value in excess]:
        final = LEO * ratio
        problem = CoplanarMultiImpulseProblem(LEO, final)
        reference = hohmann_transfer(LEO, final)
        angle = float(problem.hohmann_decision[1])
        accepted, refused = reference.first_delta_v, reference.first_delta_v * (1 - 1e-3)
        if not problem.evaluate([accepted, angle]).is_valid or problem.evaluate([refused, angle]).is_valid:
            raise RuntimeError(f"band_edge: the bisection bracket is not valid for ratio {ratio}")
        while np.nextafter(refused, accepted) < accepted:
            middle = 0.5 * (accepted + refused)
            if problem.evaluate([middle, angle]).is_valid:
                accepted = middle
            else:
                refused = middle
        edge = problem.evaluate([accepted, angle])
        rows.append({
            "ratio": ratio, "direction": "ascending" if ratio > 1 else "descending",
            "band_delta_v1_m_s": reference.first_delta_v - accepted,
            "band_delta_v1_relative": (reference.first_delta_v - accepted) / reference.first_delta_v,
            "band_delta_v1_ulps": (reference.first_delta_v - accepted) / float(np.spacing(reference.first_delta_v)),
            "apsis_deficit_m": abs(edge.metadata["arrival_radius_residual_m"]),
            "apsis_deficit_relative": abs(edge.metadata["arrival_radius_residual_m"]) / final,
            "objective_minus_hohmann_m_s": edge.objective - reference.total_delta_v,
            "advantage_times_ratio_excess_m_s": (reference.total_delta_v - edge.objective) * abs(ratio - 1),
        })
    return {"departure_radius_m": LEO, "problem": "CoplanarMultiImpulseProblem, two impulses", "rows": rows,
            "largest_advantage_m_s": max(-row["objective_minus_hohmann_m_s"] for row in rows)}


def _arc_case(initial: float, final: float, intermediate: float) -> dict[str, Any]:
    problem = CoplanarMultiImpulseProblem(initial, final, impulses=3)
    evaluation = problem.evaluate(problem.bielliptic_decision(intermediate))
    meta = evaluation.metadata
    record: dict[str, Any] = {"ratio": final / initial, "intermediate_over_larger": intermediate / max(initial, final),
                              "status": evaluation.status.value, "message": evaluation.message}
    if "reach_margin_relative" in meta:
        margin, conditioning = meta["reach_margin_relative"], meta["arc_conditioning"]
        # Constant C that the arc term would need for this decision: -Q <= K u S + C u kappa T.
        needed = max(0.0, -margin - TANGENCY_FACTOR * UNIT_ROUNDOFF) / (
            UNIT_ROUNDOFF * conditioning * meta["arc_error_scale_relative"])
        record.update({"reach_margin_relative": margin, "arc_conditioning": conditioning,
                       "needed_constant": needed, "margin_over_tolerance": -margin / meta["tangency_tolerance_relative"]})
    if evaluation.is_valid:
        reference = bielliptic_transfer(initial, final, intermediate)
        record["objective_minus_closed_form_m_s"] = evaluation.objective - reference.total_delta_v
        record["arrival_radius_residual_relative"] = meta["arrival_radius_residual_m"] / final
    return record


def _arc_summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    excluded = [r for r in records if "too close to parabolic" in r["message"]]
    rejected = [r for r in records if "does not reach" in r["message"]]
    other = [r for r in records if r["status"] != "valid" and r not in excluded and r not in rejected]
    measured = [r for r in records if "needed_constant" in r]
    worst = max(measured, key=lambda r: r["needed_constant"])
    valid = [r for r in records if r["status"] == "valid"]
    return {
        "cases": len(records), "valid": len(valid), "excluded_by_conditioning": len(excluded),
        "rejected_by_reach_test": len(rejected), "other_failures": [r["message"] for r in other],
        "reach_margin_negative": sum(r["reach_margin_relative"] < 0 for r in measured),
        "rejected_without_the_arc_term": sum(-r["reach_margin_relative"] > TANGENCY_FACTOR * UNIT_ROUNDOFF
                                             for r in measured),
        "max_arc_conditioning_of_valid": max(r["arc_conditioning"] for r in valid),
        "max_needed_constant": worst["needed_constant"],
        "max_margin_over_tolerance": max(r["margin_over_tolerance"] for r in measured),
        "worst_case": {key: worst[key] for key in ("ratio", "intermediate_over_larger", "arc_conditioning",
                                                   "needed_constant", "margin_over_tolerance")},
        "max_abs_objective_error_m_s": max(abs(r["objective_minus_closed_form_m_s"]) for r in valid),
        "max_abs_arrival_residual_relative": max(abs(r["arrival_radius_residual_relative"]) for r in valid),
    }


def arc_tolerance(quick: bool) -> dict[str, Any]:
    """Exact bi-elliptic decisions (tangent after a coast arc), ascending and descending."""
    ratios = [r for r in np.geomspace(0.02, 50.0, 8 if quick else 31) if abs(math.log(r)) > 0.05]
    factors = np.geomspace(1.001, 600.0, 6 if quick else 25)
    grid = [_arc_case(LEO, LEO * float(ratio), max(LEO, LEO * float(ratio)) * float(factor))
            for ratio in ratios for factor in factors]
    rng = np.random.Generator(np.random.PCG64(20261110))
    sample = []
    for _ in range(200 if quick else 2000):
        ratio = float(np.exp(rng.uniform(math.log(0.02), math.log(50.0))))
        if abs(math.log(ratio)) <= 0.05:
            continue
        factor = float(np.exp(rng.uniform(math.log(1.001), math.log(600.0))))
        sample.append(_arc_case(LEO, LEO * ratio, max(LEO, LEO * ratio) * factor))
    return {
        "constant_in_use": ARC_TANGENCY_FACTOR, "max_arc_conditioning": MAX_ARC_CONDITIONING,
        "domain": "r1 = 6 678 km; r2/r1 in [0.02, 50] excluding |ln ratio| <= 0.05; rb/max(r1, r2) in [1.001, 600]",
        "grid": _arc_summary(grid), "random_sample": {"seed": 20261110, **_arc_summary(sample)},
    }


@contextmanager
def _lifted_limits() -> Iterator[None]:
    """Lift the conditioning cap and the arrival guard of the multi-impulse problem.

    Only for the measurement of ``arrival_conditioning``: it shows what the
    evaluator would return beyond the limits it enforces. Restored on exit.
    """
    saved = multi_impulse.MAX_ARC_CONDITIONING, multi_impulse.ARRIVAL_RELATIVE_TOLERANCE
    multi_impulse.MAX_ARC_CONDITIONING, multi_impulse.ARRIVAL_RELATIVE_TOLERANCE = 1.0e12, 1.0e-1
    try:
        yield
    finally:
        multi_impulse.MAX_ARC_CONDITIONING, multi_impulse.ARRIVAL_RELATIVE_TOLERANCE = saved


def _residual_sample(initial: float, final: float, samples: int, seed: int, mode: str) -> dict[str, Any]:
    problem = CoplanarMultiImpulseProblem(initial, final, impulses=3)
    lower, upper = np.array(problem.bounds.lower), np.array(problem.bounds.upper)
    top = np.nextafter(upper, -np.inf)                    # periodic and open upper ends stay inside the box
    rng = np.random.Generator(np.random.PCG64(seed))
    rows, beyond = [], 0
    for _ in range(samples):
        if mode == "near_escape":    # first burn next to the escape bound, tangential; anything afterwards
            depth = 10.0 ** rng.uniform(-8.0, -1.0)
            x = np.array([upper[0] * (1 - depth), rng.normal(0.0, 1e-3), rng.uniform(lower[2], upper[2]),
                          upper[3] * 10.0 ** rng.uniform(-6.0, 0.0), rng.uniform(lower[4], upper[4])])
        else:                        # uniform in the decision box
            x = lower + rng.random(lower.size) * (upper - lower)
        meta = problem.evaluate(np.clip(x, lower, top)).metadata
        if "arrival_radius_residual_m" not in meta:
            continue                 # rejected by a predicate before the arrival
        if "final_arc_eccentricity" not in meta:
            beyond += 1              # residual above the lifted guard (0.1 r2)
            continue
        rows.append((max(meta["arc_conditioning"], arc_conditioning(meta["final_arc_eccentricity"])),
                     abs(meta["arrival_radius_residual_m"]) / final))
    data = np.array(rows).reshape(-1, 2)
    decades = []
    for exponent in range(0, 12):
        chosen = data[(data[:, 0] >= 10.0 ** exponent) & (data[:, 0] < 10.0 ** (exponent + 1)), 1]
        if chosen.size:
            decades.append({"conditioning_from": 10.0 ** exponent, "conditioning_to": 10.0 ** (exponent + 1),
                            "evaluations": int(chosen.size), "median": float(np.median(chosen)),
                            "max": float(chosen.max()),
                            "fraction_above_tolerance": float(np.mean(chosen > ARRIVAL_RELATIVE_TOLERANCE))})
    inside = data[data[:, 0] <= MAX_ARC_CONDITIONING, 1]
    return {"transfer": f"{initial:g} m -> {final:g} m, three impulses", "sampling": mode, "seed": seed,
            "samples": samples, "evaluations_with_an_arrival": int(data.shape[0]) + beyond,
            "residual_above_0.1_r2": beyond, "by_decade_of_conditioning": decades,
            "evaluations_inside_the_cap": int(inside.size),
            "max_relative_residual_inside_the_cap": float(inside.max()) if inside.size else None}


def arrival_conditioning(quick: bool) -> dict[str, Any]:
    """|arrival radius - r2| / r2 against the largest kappa of the arcs flown, limits lifted."""
    samples = 1500 if quick else 30000
    with _lifted_limits():
        studies = [_residual_sample(LEO, GEO, samples, 11, "near_escape"),
                   _residual_sample(LEO, GEO, samples, 13, "uniform"),
                   _residual_sample(GEO, LEO, samples, 12, "uniform")]
    worst = max(study["max_relative_residual_inside_the_cap"] or 0.0 for study in studies)
    return {"note": "conditioning cap and arrival guard lifted for this measurement only",
            "conditioning": "largest kappa = (3+e)/(1-e) over the coast arcs and the final arc",
            "arrival_relative_tolerance": ARRIVAL_RELATIVE_TOLERANCE, "max_arc_conditioning": MAX_ARC_CONDITIONING,
            "studies": studies, "max_relative_residual_inside_the_cap": worst,
            "tolerance_over_worst_residual_inside_the_cap": ARRIVAL_RELATIVE_TOLERANCE / worst if worst else None}


def near_parabolic(quick: bool) -> dict[str, Any]:
    """Three-impulse decisions next to the escape bound of the first impulse, limits in force."""
    problem = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)
    upper = float(problem.bounds.upper[0])
    rng = np.random.Generator(np.random.PCG64(7))
    counts: dict[str, int] = {}
    messages: dict[str, int] = {}
    worst_residual, longest = 0.0, 0.0
    samples = 2000 if quick else 20000
    for _ in range(samples):
        depth = 10.0 ** rng.uniform(-13.0, -1.0)            # relative distance to the escape bound
        x = np.array([upper * (1 - depth), rng.normal(0.0, 1e-3), rng.uniform(0.0, 2 * np.pi),
                      upper * 10.0 ** rng.uniform(-15.0, 0.0), rng.uniform(-np.pi, np.pi)])
        x[1] = float(np.clip(x[1], -np.pi, np.pi - 1e-12))
        x[3] = min(x[3], upper)
        evaluation = problem.evaluate(x)
        counts[evaluation.status.value] = counts.get(evaluation.status.value, 0) + 1
        if evaluation.is_valid:
            worst_residual = max(worst_residual, abs(evaluation.metadata["arrival_radius_residual_m"]))
            longest = max(longest, evaluation.metadata["transfer_time_s"])
        else:
            reason = evaluation.message.split(" (")[0].split(" by ")[0]
            messages[reason] = messages.get(reason, 0) + 1
    named = {"A": [3200.1474929757806, 3.26e-7, 2.5883, 1.196e-12, -1.8058],
             "B": [upper * (1 - 1e-9), 0.0, 1e-9, 1.60007, math.pi - 1e-12],
             "C": [3200.146, 8.2e-4, 4.7117, 8.64e-3, -0.386]}
    cases = {}
    for label, decision in named.items():
        decision[0] = min(decision[0], upper)
        evaluation = problem.evaluate(decision)
        cases[label] = {"status": evaluation.status.value, "message": evaluation.message}
    return {
        "samples": samples, "status_counts": counts, "rejection_reasons": messages,
        "arrival_tolerance_m": ARRIVAL_RELATIVE_TOLERANCE * GEO,
        "max_abs_arrival_residual_of_valid_m": worst_residual, "longest_valid_transfer_s": longest,
        "audit_cases": cases, "arrival_failures": problem.arrival_failures,
        "captured_value_errors": problem.captured_value_errors,
        "captured_runtime_errors": problem.captured_runtime_errors,
    }


def coincident_split() -> dict[str, Any]:
    """Second burn at the apoapsis of the Hohmann ellipse, i.e. on the target radius."""
    problem = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)
    reference = hohmann_transfer(LEO, GEO)
    rows = []
    for sweep, fraction in ((math.pi, 0.25), (math.pi, 0.5), (math.pi - 1e-3, 0.5), (math.pi + 1e-6, 0.5),
                            (math.pi + 1e-2, 0.5)):
        evaluation = problem.evaluate([reference.first_delta_v, 0.0, sweep, fraction * reference.second_delta_v, 0.0])
        rows.append({"coast_minus_pi_rad": sweep - math.pi, "second_burn_fraction": fraction,
                     "status": evaluation.status.value, "message": evaluation.message,
                     "transfer_time_s": evaluation.metadata.get("transfer_time_s"),
                     "objective_minus_hohmann_m_s": (evaluation.objective - reference.total_delta_v)
                     if evaluation.is_valid else None})
    return {"hohmann_transfer_time_s": reference.transfer_time, "rows": rows,
            "captured_runtime_errors": problem.captured_runtime_errors}


def seam(quick: bool) -> dict[str, Any]:
    """DE on GEO->LEO with the angle cut at the retrograde optimum (alpha = +-pi) and opposite to it."""
    reference = hohmann_transfer(GEO, LEO).total_delta_v
    budget, seeds = (600, range(4)) if quick else (2000, range(30))
    result: dict[str, Any] = {"transfer": "GEO -> LEO, two impulses", "evaluations": budget, "seeds": list(seeds),
                              "success_threshold_m_s": SEAM_SUCCESS_M_S}
    gaps: dict[str, list[float]] = {}
    for label in ("cut_opposite_to_the_optimum", "cut_at_the_optimum"):
        gaps[label] = []
        for seed in seeds:
            problem = CoplanarMultiImpulseProblem(GEO, LEO)
            if label == "cut_at_the_optimum":
                bounds = problem.bounds
                problem.bounds = Bounds([0.0, -np.pi], [float(bounds.upper[0]), np.pi], names=bounds.names,
                                        units=bounds.units, periodic=np.array(bounds.periodic))
            outcome = DifferentialEvolution().minimize(problem, TerminationCriteria(budget), seed=seed)
            gaps[label].append(outcome.f_best - reference)
        successes = sum(gap <= SEAM_SUCCESS_M_S for gap in gaps[label])
        result[label] = {"gaps_m_s": gaps[label], "summary_m_s": summarize(gaps[label]),
                         "successes": successes, "wilson_95": list(wilson_interval(successes, len(gaps[label])))}
    result["mann_whitney_opposite_vs_at"] = mann_whitney(gaps["cut_opposite_to_the_optimum"],
                                                         gaps["cut_at_the_optimum"])
    return result


def long_transfer() -> dict[str, Any]:
    """Bi-elliptic transfer r2 = 20 r1 through rb = 60 r1: default against tightened RK45."""
    problem = CoplanarMultiImpulseProblem(LEO, 20 * LEO, impulses=3)
    plan = problem.plan(problem.bielliptic_decision(60 * LEO))
    rows = {}
    for label, options in (("default", {}), ("tightened", {"rtol": 1e-12, "atol": [1e-6] * 3 + [1e-9] * 3})):
        outcome = verify_impulsive_plan(plan, **options)
        rows[label] = {"passed": outcome.passed, "message": outcome.message, "metrics": dict(outcome.metrics),
                       "options": options or "propagate defaults: rtol 1e-10, atol 1e-3 m and 1e-6 m/s, max_step 60 s"}
    return {"transfer_time_days": plan.transfer_time / 86400.0, "delta_v_total_m_s": plan.total_delta_v,
            "hohmann_delta_v_m_s": hohmann_transfer(LEO, 20 * LEO).total_delta_v, "verification": rows}


def rosenbrock(quick: bool) -> dict[str, Any]:
    """Final gap and normalized population spread of DE on Rosenbrock-5D (benchmark box and budget)."""
    benchmark = rosenbrock_benchmark(5, -2.048, 4.096)
    budget, seeds = (5000, range(3)) if quick else (50000, range(10))
    rows = []
    for seed in seeds:
        outcome = DifferentialEvolution().minimize(benchmark.problem(), TerminationCriteria(budget), seed=seed)
        best = outcome.history.generation_best
        improved = np.nonzero(np.diff(best) < 0)[0]
        rows.append({"seed": seed, "final_gap": outcome.f_best - benchmark.minimum,
                     "final_spread": float(outcome.history.generation_spread[-1]),
                     "last_improving_generation": int(improved[-1]) + 1 if improved.size else 0,
                     "generations": outcome.generations})
    stalled = [row for row in rows if row["final_gap"] > 1e-8]
    return {"evaluations": budget, "seeds": [row["seed"] for row in rows], "rows": rows,
            "solved_to_1e-8": len(rows) - len(stalled),
            "unsolved_still_improving_in_the_last_three_generations":
                sum(row["last_improving_generation"] >= row["generations"] - 2 for row in stalled),
            "max_final_spread_of_unsolved": max((row["final_spread"] for row in stalled), default=None),
            "min_final_spread_of_unsolved": min((row["final_spread"] for row in stalled), default=None)}


def evaluation_cost(quick: bool) -> dict[str, Any]:
    """Best of five timings of one evaluation; single process; depends on the machine."""
    cases = [("reference two-impulse, valid", CoplanarTwoImpulseProblem(LEO, GEO), [2600.0, 0.3]),
             ("general two-impulse, valid", CoplanarMultiImpulseProblem(LEO, GEO), [2600.0, 0.3]),
             ("general three-impulse, valid", CoplanarMultiImpulseProblem(LEO, GEO, impulses=3),
              [2425.8, 0.0, 3.0, 40.0, 0.0]),
             ("general two-impulse, rejected by a predicate", CoplanarMultiImpulseProblem(LEO, GEO), [0.0, 0.0])]
    repeats = 50 if quick else 300
    rows = {}
    for name, problem, decision in cases:
        point = np.array(decision)
        best = math.inf
        for _ in range(5):
            start = perf_counter()
            for _ in range(repeats):
                problem.evaluate(point)
            best = min(best, (perf_counter() - start) / repeats)
        rows[name] = {"seconds_per_evaluation": best, "status": problem.evaluate(point).status.value}
    return {"note": "wall time, best of five repetitions, one process; informative only", "rows": rows}


def run(quick: bool = False) -> dict[str, Any]:
    """Compute every section and the checks that the documentation relies on."""
    sections = {
        "two_impulse_tangency": two_impulse_tangency(quick),
        "band_edge": band_edge(quick),
        "arc_tolerance": arc_tolerance(quick),
        "arrival_conditioning": arrival_conditioning(quick),
        "near_parabolic": near_parabolic(quick),
        "coincident_split": coincident_split(),
        "seam": seam(quick),
        "long_transfer": long_transfer(),
        "rosenbrock": rosenbrock(quick),
        "evaluation_cost": evaluation_cost(quick),
    }
    arcs, parabolic = sections["arc_tolerance"], sections["near_parabolic"]
    worst_inside = sections["arrival_conditioning"]["max_relative_residual_inside_the_cap"]
    checks = {
        "hohmann_never_rejected_by_the_tangency_tolerance": sections["two_impulse_tangency"]["hohmann_rejected"] == 0,
        "no_exact_bielliptic_decision_rejected_in_the_supported_domain":
            arcs["grid"]["rejected_by_reach_test"] == 0 and arcs["random_sample"]["rejected_by_reach_test"] == 0,
        "no_unexplained_failure_of_exact_bielliptic_decisions":
            not arcs["grid"]["other_failures"] and not arcs["random_sample"]["other_failures"],
        "arc_constant_at_least_four_times_the_largest_needed":
            ARC_TANGENCY_FACTOR >= 4 * max(arcs["grid"]["max_needed_constant"],
                                           arcs["random_sample"]["max_needed_constant"]),
        "arrival_residual_inside_the_cap_ten_times_below_the_tolerance":
            worst_inside <= ARRIVAL_RELATIVE_TOLERANCE / 10,
        "every_valid_near_parabolic_evaluation_arrives_on_the_target":
            parabolic["max_abs_arrival_residual_of_valid_m"] <= parabolic["arrival_tolerance_m"],
        "audit_cases_are_not_valid": all(case["status"] != "valid" for case in parabolic["audit_cases"].values()),
        "coincident_final_impulse_is_invalid_domain":
            all(row["status"] == EvaluationStatus.INVALID_DOMAIN.value
                for row in sections["coincident_split"]["rows"] if row["coast_minus_pi_rad"] == 0.0),
        "no_physics_exception_captured":
            parabolic["captured_value_errors"] == 0 and parabolic["captured_runtime_errors"] == 0
            and sections["coincident_split"]["captured_runtime_errors"] == 0,
        "long_transfer_fails_by_default_and_passes_tightened":
            (not sections["long_transfer"]["verification"]["default"]["passed"])
            and sections["long_transfer"]["verification"]["tightened"]["passed"],
    }
    return {
        "experiment": "phase3_numerical_evidence", "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": {"quick": quick, "leo_radius_m": LEO, "geo_radius_m": GEO,
                          "tangency_factor": TANGENCY_FACTOR, "arc_tangency_factor": ARC_TANGENCY_FACTOR,
                          "max_arc_conditioning": MAX_ARC_CONDITIONING,
                          "arrival_relative_tolerance": ARRIVAL_RELATIVE_TOLERANCE},
        "environment": environment_record(), "code": code_version(),
        "sections": sections, "checks": checks, "passed": all(checks.values()),
        "limitations": [
            "The sweeps cover the stated grids and samples; they are measurements, not proofs of a bound.",
            "The seam comparison is one transfer, one optimizer and one budget; it is not a general result.",
            "evaluation_cost depends on the machine and its load.",
            "With --quick the grids are too small to serve as evidence.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Numerical evidence of the Phase 3 tolerances and remarks.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quick", action="store_true", help="small grids and budgets: smoke run, not evidence")
    args = parser.parse_args()
    try:
        report = run(args.quick)
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "report.json").write_text(json.dumps(jsonable(report), indent=2, allow_nan=False) + "\n",
                                                 encoding="utf-8")
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"Evidence run failed: {exc}\n")
    print(json.dumps(jsonable({"passed": report["passed"], "checks": report["checks"],
                               "report": str((args.output / "report.json").resolve())}), indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
