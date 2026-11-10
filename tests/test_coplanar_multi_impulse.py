"""Generalized coplanar transfer problem (Phase 3 increment d): N impulses, both directions.

Independent references:
  - N = 2 ascending: the reference problem problems.coplanar_transfer, itself
    validated against a closed-form oracle (tests/test_coplanar_transfer_problem.py).
  - N = 2 descending: the Hohmann closed form and an RK45 propagation.
  - N = 3: the bi-elliptic closed form (tests/test_arc_extrema_and_bielliptic.py),
    the invariance of the cost under a null intermediate burn, and an RK45
    propagation of the resulting plan.
  - Optimality: below r2/r1 = 11.94 no time-free impulsive transfer between
    coplanar circular orbits is cheaper than Hohmann (Edelbaum, "How many
    impulses?", 1967; docs/03 section 8.3), so no feasible sample may beat it.

Thresholds: eps_low = 1e-6 m/s ("not significantly better than the
reference") as in the reference problem's tests. Optimizer runs with small
budgets assert properties only (never better than the reference, plan
verified by propagation); convergence statistics belong to the multi-seed
benchmark.
"""
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import EARTH_MU, bielliptic_transfer, hohmann_transfer
from orbital_optimizer.optimization import (
    DifferentialEvolution, EvaluationStatus, ScalarizedProblem, TerminationCriteria,
)
from orbital_optimizer.problems import (
    CoplanarMultiImpulseProblem, CoplanarTwoImpulseProblem, ImpulsivePlan, PlanVerification,
    hohmann_two_impulse_plan, transfer_objective, verify_impulsive_plan,
)
from orbital_optimizer.problems.coplanar_multi_impulse import ARC_TANGENCY_FACTOR
from orbital_optimizer.problems.coplanar_transfer import TANGENCY_FACTOR

LEO, GEO = 6_678_000.0, 42_164_000.0
EPS_LOW = 1e-6
UNIT_ROUNDOFF = np.finfo(float).eps / 2


def tangency_time_tolerance(meta, apsis: str, mu: float = EARTH_MU) -> float:
    """A-priori bound (s) on the transfer-time error when the arrival is tangent to the target.

    The crossing anomaly is arccos(c) with c = +-1 at tangency, where
    arccos has square-root conditioning: an error d in c moves the angle by
    sqrt(2 d). The reach test admits |d| of order F u (F = tangency_factor,
    u = unit roundoff); allowing a factor 8 for the division by e in
    c = (p/r2 - 1)/e gives d_nu <= sqrt(2 * 8 F u) = 4 sqrt(F u). The vehicle
    sweeps that angle at the apsis rate nu_dot = n sqrt(1-e^2)/(1 -+ e)^2
    (periapsis/apoapsis), so dt <= 4 sqrt(F u) (1 -+ e)^2 / (n sqrt(1-e^2)).
    """
    a, e = meta["final_arc_semi_major_axis_m"], meta["final_arc_eccentricity"]
    mean_motion = math.sqrt(mu / a**3)
    factor = (1 - e) ** 2 if apsis == "periapsis" else (1 + e) ** 2
    return 4 * math.sqrt(meta["tangency_factor"] * UNIT_ROUNDOFF) * factor / (mean_motion * math.sqrt(1 - e * e))


class ConstructionTests(unittest.TestCase):
    def test_bounds_for_two_and_three_impulses(self) -> None:
        two = CoplanarMultiImpulseProblem(LEO, GEO)
        reference = CoplanarTwoImpulseProblem(LEO, GEO)
        self.assertEqual(two.name, "coplanar-2-impulse")
        self.assertEqual(two.bounds.names, ("delta_v1", "alpha1"))
        assert_allclose(two.bounds.lower, reference.bounds.lower, rtol=0, atol=0)
        assert_allclose(two.bounds.upper, reference.bounds.upper, rtol=0, atol=0)
        three = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)
        self.assertEqual(three.bounds.dimension, 5)
        self.assertEqual(three.bounds.names, ("delta_v1", "alpha1", "coast2", "delta_v2", "alpha2"))
        self.assertEqual(three.bounds.units, ("m/s", "rad", "rad", "m/s", "rad"))
        self.assertEqual(three.bounds.periodic.tolist(), [False, True, True, False, True])
        self.assertEqual(CoplanarMultiImpulseProblem(LEO, GEO, impulses=5).bounds.dimension, 11)
        upper = (math.sqrt(2) - 1) * math.sqrt(EARTH_MU / LEO)
        self.assertAlmostEqual(three.bounds.upper[0], upper, places=9)
        self.assertAlmostEqual(three.bounds.upper[3], upper, places=9)

    def test_descending_angles_are_cut_opposite_to_the_retrograde_direction(self) -> None:
        down = CoplanarMultiImpulseProblem(GEO, LEO)
        self.assertFalse(down.ascending)
        self.assertEqual((down.bounds.lower[1], down.bounds.upper[1]), (0.0, 2 * np.pi))
        self.assertEqual(down.hohmann_decision.tolist(), [hohmann_transfer(GEO, LEO).first_delta_v, np.pi])
        # The bound is the escape excess of the lower orbit and exceeds the Hohmann burn.
        self.assertAlmostEqual(down.bounds.upper[0], (math.sqrt(2) - 1) * math.sqrt(EARTH_MU / LEO), places=9)
        self.assertGreater(down.bounds.upper[0], hohmann_transfer(GEO, LEO).first_delta_v)

    def test_invalid_construction(self) -> None:
        for arguments in ((LEO, LEO), (-1.0, GEO), (LEO, np.nan), (0.0, GEO)):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                CoplanarMultiImpulseProblem(*arguments)
        for impulses in (1, 0, True, 2.0, "2"):
            with self.subTest(impulses=impulses), self.assertRaises(ValueError):
                CoplanarMultiImpulseProblem(LEO, GEO, impulses=impulses)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            CoplanarMultiImpulseProblem(LEO, GEO, mu=0.0)
        with self.assertRaises(RuntimeError):
            CoplanarMultiImpulseProblem(1e-10, 2e-10, mu=1e308)
        with self.assertRaises(ValueError):
            CoplanarMultiImpulseProblem(LEO, GEO, impulses=3).hohmann_decision
        with self.assertRaises(ValueError):
            CoplanarMultiImpulseProblem(LEO, GEO).bielliptic_decision(2 * GEO)


class TwoImpulseTests(unittest.TestCase):
    def test_ascending_case_agrees_with_the_reference_problem(self) -> None:
        general, reference = CoplanarMultiImpulseProblem(LEO, GEO), CoplanarTwoImpulseProblem(LEO, GEO)
        rng = np.random.Generator(np.random.PCG64(0))
        valid = 0
        for u in rng.random((3000, 2)):
            x = reference.bounds.from_unit(u)
            first, second = general.evaluate(x), reference.evaluate(x)
            self.assertIs(first.status, second.status, x)
            if second.is_valid:
                valid += 1
                # Same formulas on the same inputs: 1e-9 m/s and s only leaves room for a
                # different summation order across platforms.
                self.assertAlmostEqual(first.objective, second.objective, delta=1e-9)
                self.assertAlmostEqual(first.metadata["transfer_time_s"], second.metadata["transfer_time_s"], delta=1e-9)
                assert_allclose(first.metadata["burn_vectors_m_s"][0], second.metadata["first_burn_m_s"], atol=1e-9)
                assert_allclose(first.metadata["burn_vectors_m_s"][1], second.metadata["second_burn_m_s"], atol=1e-9)
        self.assertGreater(valid, 100)
        self.assertEqual((general.captured_value_errors, general.captured_runtime_errors), (0, 0))

    def test_hohmann_decision_in_both_directions(self) -> None:
        for r1, r2 in ((LEO, GEO), (GEO, LEO)):
            with self.subTest(r1=r1, r2=r2):
                problem, reference = CoplanarMultiImpulseProblem(r1, r2), hohmann_transfer(r1, r2)
                evaluation = problem.evaluate(problem.hohmann_decision)
                self.assertIs(evaluation.status, EvaluationStatus.VALID)
                meta = evaluation.metadata
                self.assertLessEqual(abs(evaluation.objective - reference.total_delta_v), EPS_LOW)
                self.assertAlmostEqual(meta["transfer_time_s"], reference.transfer_time, places=6)
                assert_allclose(meta["burn_delta_v_m_s"], [reference.first_delta_v, reference.second_delta_v],
                                rtol=0, atol=EPS_LOW)
                assert_allclose(meta["burn_times_s"], [0.0, reference.transfer_time], rtol=0, atol=1e-6)
                # The Hohmann ellipse stays between the two circles.
                self.assertAlmostEqual(meta["minimum_radius_m"], min(r1, r2), delta=1e-3)
                self.assertAlmostEqual(meta["maximum_radius_m"], max(r1, r2), delta=1e-3)
                self.assertAlmostEqual(meta["closest_approach_index"], 1.0, places=9)
                self.assertAlmostEqual(meta["impulse_kinetic_energy_j_kg"],
                                       (reference.first_delta_v**2 + reference.second_delta_v**2) / 2, delta=1e-3)
                self.assertEqual(meta["tangency_factor"], TANGENCY_FACTOR)
                self.assertEqual(meta["impulses"], 2)

    def test_descending_plan_matches_the_oriented_hohmann_plan_and_propagates(self) -> None:
        problem = CoplanarMultiImpulseProblem(GEO, LEO)
        plan = problem.plan(problem.hohmann_decision)
        reference = hohmann_two_impulse_plan(GEO, LEO)
        self.assertIsInstance(plan, ImpulsivePlan)
        # sin(pi) is 1.2e-16, so the "tangential" first burn has a 2e-13 m/s radial component.
        assert_allclose(plan.burn_vectors[0], reference.first_burn, rtol=0, atol=1e-9)
        assert_allclose(plan.burn_vectors[1], reference.second_burn, rtol=0, atol=EPS_LOW)
        self.assertAlmostEqual(plan.transfer_time, reference.transfer_time, places=6)
        verification = verify_impulsive_plan(plan)
        self.assertIsInstance(verification, PlanVerification)
        self.assertTrue(verification.passed, verification.message)
        self.assertLess(abs(verification.metrics["final_semi_major_axis_error_m"]), 1.0)

    def test_minimum_radius_covers_the_whole_arc(self) -> None:
        """An inward first burn dips below the departure radius before rising to the target.

        Both end points of the transfer arc (r1 and r2) are then above the
        periapsis: a check of the end points alone would miss it (RV-4.02).
        """
        problem, reference = CoplanarMultiImpulseProblem(LEO, GEO), CoplanarTwoImpulseProblem(LEO, GEO)
        inward = problem.evaluate([2600.0, -0.3])
        periapsis = reference.evaluate([2600.0, -0.3]).metadata["periapsis_radius_m"]
        self.assertLess(periapsis, LEO - 1e4)
        self.assertAlmostEqual(inward.metadata["minimum_radius_m"], periapsis, delta=1e-3)
        self.assertGreater(inward.metadata["closest_approach_index"], 1.0)
        self.assertAlmostEqual(inward.metadata["maximum_radius_m"], GEO, delta=1e-3)
        outward = problem.evaluate([2600.0, 0.3])
        # Mirror transfer: same ellipse, but the vehicle leaves outward and never visits the
        # periapsis. The start radius is recomputed from the conic equation: a few ulps of a.
        self.assertAlmostEqual(outward.metadata["minimum_radius_m"], LEO, delta=1e-6)
        self.assertAlmostEqual(outward.metadata["closest_approach_index"], 1.0, places=12)
        self.assertAlmostEqual(outward.objective, inward.objective, delta=1e-6)

    def test_domain_predicates(self) -> None:
        problem = CoplanarMultiImpulseProblem(GEO, LEO)
        upper = problem.bounds.upper[0]
        cases = {
            (0.0, np.pi): "does not reach",
            (upper, 0.0): "unbound",                        # prograde burn beyond the local escape excess
            (math.sqrt(EARTH_MU / GEO) * 1.01, np.pi): "radial or retrograde",
        }
        for x, text in cases.items():
            with self.subTest(x=x):
                evaluation = problem.evaluate(x)
                self.assertIs(evaluation.status, EvaluationStatus.INVALID_DOMAIN)
                self.assertIn(text, evaluation.message)
        self.assertEqual((problem.captured_value_errors, problem.captured_runtime_errors), (0, 0))
        with self.assertRaises(ValueError):
            problem.evaluate([100.0, 2 * np.pi])            # periodic upper bound excluded
        with self.assertRaises(ValueError):
            problem.plan([0.0, np.pi])

    def test_no_feasible_descending_sample_beats_hohmann(self) -> None:
        problem = CoplanarMultiImpulseProblem(GEO, LEO)
        reference = hohmann_transfer(GEO, LEO).total_delta_v
        rng = np.random.Generator(np.random.PCG64(2026))
        objectives = [evaluation.objective for evaluation in
                      (problem.evaluate(x) for x in problem.bounds.from_unit(rng.random((6000, 2))))
                      if evaluation.is_valid]
        self.assertGreater(len(objectives), 500)
        self.assertGreaterEqual(min(objectives), reference - EPS_LOW)


class ThreeImpulseTests(unittest.TestCase):
    def test_bielliptic_decision_reproduces_the_closed_form(self) -> None:
        # (r2/r1, rb/max(r1, r2)): ascending and descending, moderate and high eccentricity.
        for ratio, factor in ((5.0, 8.0), (20.0, 3.0), (20.0, 60.0), (0.2, 3.0), (0.05, 2.0)):
            with self.subTest(ratio=ratio, factor=factor):
                r2 = ratio * LEO
                rb = factor * max(LEO, r2)
                problem = CoplanarMultiImpulseProblem(LEO, r2, impulses=3)
                reference = bielliptic_transfer(LEO, r2, rb)
                evaluation = problem.evaluate(problem.bielliptic_decision(rb))
                self.assertIs(evaluation.status, EvaluationStatus.VALID, evaluation.message)
                meta = evaluation.metadata
                self.assertLessEqual(abs(evaluation.objective - reference.total_delta_v), EPS_LOW)
                assert_allclose(meta["burn_delta_v_m_s"],
                                [reference.first_delta_v, reference.second_delta_v, reference.third_delta_v],
                                rtol=0, atol=EPS_LOW)
                # The arrival is tangent at the periapsis of the second ellipse: the crossing time
                # has square-root conditioning there (see tangency_time_tolerance). The propagation
                # itself adds dt/t <= 1.5 * 4 eps (3+e)/(1-e) per arc (tests/test_kepler.py), i.e.
                # at most 3.2e-12 relative for kappa <= 2400, covered by the 1e-10 term.
                time_tolerance = tangency_time_tolerance(meta, "periapsis") + 1e-10 * reference.transfer_time
                self.assertLessEqual(abs(meta["transfer_time_s"] - reference.transfer_time), time_tolerance)
                self.assertLess(time_tolerance, 1e-6 * reference.transfer_time)
                self.assertAlmostEqual(meta["maximum_radius_m"] / rb, 1.0, delta=1e-10)
                self.assertAlmostEqual(meta["minimum_radius_m"] / min(LEO, r2), 1.0, delta=1e-10)
                first_arc_eccentricity = (rb - LEO) / (rb + LEO)
                self.assertAlmostEqual(
                    meta["tangency_factor"],
                    TANGENCY_FACTOR + ARC_TANGENCY_FACTOR * (3 + first_arc_eccentricity) / (1 - first_arc_eccentricity),
                    delta=1e-6 * meta["tangency_factor"])
                self.assertEqual((problem.captured_value_errors, problem.captured_runtime_errors), (0, 0))

    def test_three_impulses_beat_hohmann_above_the_second_threshold(self) -> None:
        ratio, rb = 20.0, 60.0 * 20.0 * LEO
        problem = CoplanarMultiImpulseProblem(LEO, ratio * LEO, impulses=3)
        evaluation = problem.evaluate(problem.bielliptic_decision(rb))
        hohmann = hohmann_transfer(LEO, ratio * LEO).total_delta_v
        self.assertLess(evaluation.objective, hohmann - 10.0)
        self.assertEqual(evaluation.metadata["reference_delta_v_m_s"], hohmann)

    def test_bielliptic_plan_is_verified_by_propagation(self) -> None:
        problem = CoplanarMultiImpulseProblem(LEO, 3.0 * LEO, impulses=3)
        plan = problem.plan(problem.bielliptic_decision(5.0 * LEO))
        self.assertEqual(plan.burn_vectors.shape, (3, 3))
        self.assertTrue(np.all(np.diff(plan.burn_times) > 0))
        verification = verify_impulsive_plan(plan)
        self.assertTrue(verification.passed, verification.message)
        self.assertAlmostEqual(plan.total_delta_v, bielliptic_transfer(LEO, 3.0 * LEO, 5.0 * LEO).total_delta_v, delta=EPS_LOW)

    def test_null_intermediate_burn_reduces_to_hohmann(self) -> None:
        """A zero second burn anywhere on the transfer ellipse leaves the Hohmann transfer unchanged."""
        problem = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)
        reference = hohmann_transfer(LEO, GEO)
        for sweep in (0.3, 1.5, 3.0):          # mean-anomaly sweep before apoapsis (pi)
            for alpha in (0.0, 1.0, -2.0):     # direction of a zero burn is irrelevant
                with self.subTest(sweep=sweep, alpha=alpha):
                    evaluation = problem.evaluate([reference.first_delta_v, 0.0, sweep, 0.0, alpha])
                    self.assertIs(evaluation.status, EvaluationStatus.VALID, evaluation.message)
                    self.assertLessEqual(abs(evaluation.objective - reference.total_delta_v), EPS_LOW)
                    # Tangent arrival at the apoapsis after a coast arc: square-root conditioning.
                    time_tolerance = tangency_time_tolerance(evaluation.metadata, "apoapsis")
                    self.assertLess(time_tolerance, 0.05)
                    self.assertLessEqual(abs(evaluation.metadata["transfer_time_s"] - reference.transfer_time),
                                         time_tolerance)
                    self.assertEqual(evaluation.metadata["burn_delta_v_m_s"][1], 0.0)

    def test_coincident_impulses_are_rejected(self) -> None:
        problem = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)
        evaluation = problem.evaluate([2400.0, 0.0, 0.0, 30.0, 0.0])
        self.assertIs(evaluation.status, EvaluationStatus.INVALID_DOMAIN)
        self.assertIn("coincident impulses", evaluation.message)

    def test_no_feasible_sample_beats_hohmann_below_the_first_threshold(self) -> None:
        problem = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)      # r2/r1 = 6.31 < 11.94
        reference = hohmann_transfer(LEO, GEO).total_delta_v
        rng = np.random.Generator(np.random.PCG64(7))
        objectives = [evaluation.objective for evaluation in
                      (problem.evaluate(x) for x in problem.bounds.from_unit(rng.random((4000, 5))))
                      if evaluation.is_valid]
        self.assertGreater(len(objectives), 200)
        self.assertGreaterEqual(min(objectives), reference - EPS_LOW)
        self.assertEqual((problem.captured_value_errors, problem.captured_runtime_errors), (0, 0))


class PlanContractTests(unittest.TestCase):
    def test_plan_validation_and_immutability(self) -> None:
        state = [LEO, 0.0, 0.0, 0.0, math.sqrt(EARTH_MU / LEO), 0.0]
        plan = ImpulsivePlan(state, [0.0, 100.0], [[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]], GEO)
        self.assertEqual(plan.total_delta_v, 3.0)
        self.assertEqual(plan.transfer_time, 100.0)
        self.assertEqual([time for time, _ in plan.impulses()], [0.0, 100.0])
        with self.assertRaises(ValueError):
            plan.burn_vectors[0, 0] = 1.0
        valid = {"initial_state": state, "burn_times": [0.0, 100.0],
                 "burn_vectors": [[0.0, 1.0, 0.0], [0.0, 2.0, 0.0]], "final_radius": GEO}
        for change in ({"burn_times": [100.0, 100.0]}, {"burn_times": [-1.0, 5.0]}, {"burn_times": []},
                       {"burn_times": [[0.0, 1.0]]}, {"burn_vectors": [[0.0, 1.0, 0.0]]},
                       {"burn_vectors": [[0.0, 1.0], [0.0, 2.0]]}, {"final_radius": 0.0}, {"mu": -1.0},
                       {"initial_state": [0.0] * 6}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                ImpulsivePlan(**(valid | change))

    def test_verification_detects_wrong_plans(self) -> None:
        problem = CoplanarMultiImpulseProblem(LEO, GEO)
        plan = problem.plan(problem.hohmann_decision)
        incomplete = ImpulsivePlan(plan.initial_state, plan.burn_times,
                                   [plan.burn_vectors[0], [0.0, 0.0, 0.0]], GEO)
        result = verify_impulsive_plan(incomplete)
        self.assertFalse(result.passed)
        self.assertFalse(result.checks["final_eccentricity"])
        self.assertGreater(result.metrics["final_eccentricity"], 0.7)
        wrong_target = ImpulsivePlan(plan.initial_state, plan.burn_times, plan.burn_vectors, GEO + 5.0)
        self.assertFalse(verify_impulsive_plan(wrong_target).passed)
        single = ImpulsivePlan(plan.initial_state, [0.0], [[0.0, 10.0, 0.0]], GEO)
        self.assertIn("no transfer arc", verify_impulsive_plan(single).message)
        with self.assertRaises(ValueError):
            verify_impulsive_plan("plan")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            verify_impulsive_plan(plan, radius_tolerance_m=0.0)

    def test_verification_tolerances_are_explicit(self) -> None:
        problem = CoplanarMultiImpulseProblem(LEO, GEO)
        plan = problem.plan(problem.hohmann_decision)
        default = verify_impulsive_plan(plan)
        self.assertEqual(dict(default.tolerances), {"radius_m": 1.0, "eccentricity": 1e-6})
        # The measured semi-major-axis error (0.064 m, docs/phase1-phase2-audit.md) fails a 1 cm bound.
        strict = verify_impulsive_plan(plan, radius_tolerance_m=0.01)
        self.assertFalse(strict.passed)
        self.assertEqual(strict.tolerances["radius_m"], 0.01)


class OptimizerEndToEndTests(unittest.TestCase):
    def test_descending_run_never_beats_hohmann_and_its_plan_propagates(self) -> None:
        reference = hohmann_transfer(GEO, LEO).total_delta_v
        problem = CoplanarMultiImpulseProblem(GEO, LEO)
        result = DifferentialEvolution().minimize(problem, TerminationCriteria(max_evaluations=2000), seed=0)
        self.assertTrue(result.found_valid)
        self.assertGreaterEqual(result.f_best - reference, -EPS_LOW)
        self.assertEqual((problem.captured_value_errors, problem.captured_runtime_errors), (0, 0))
        self.assertTrue(verify_impulsive_plan(problem.plan(result.x_best)).passed)

    def test_three_impulse_run_never_beats_hohmann_and_its_plan_propagates(self) -> None:
        reference = hohmann_transfer(LEO, GEO).total_delta_v
        problem = CoplanarMultiImpulseProblem(LEO, GEO, impulses=3)
        result = DifferentialEvolution().minimize(problem, TerminationCriteria(max_evaluations=1000), seed=0)
        self.assertTrue(result.found_valid)
        self.assertGreaterEqual(result.f_best - reference, -EPS_LOW)
        plan = problem.plan(result.x_best)
        self.assertEqual(len(plan.burn_times), 3)
        self.assertTrue(verify_impulsive_plan(plan).passed)


class TransferObjectiveTests(unittest.TestCase):
    def setUp(self) -> None:
        self.problem = CoplanarMultiImpulseProblem(LEO, GEO)
        self.reference = hohmann_transfer(LEO, GEO)

    def test_every_term_equals_one_at_the_hohmann_transfer(self) -> None:
        objective = transfer_objective(self.reference, delta_v=1.0, time=1.0, energy=1.0, risk=1.0,
                                       rationale="unit weights to expose every term")
        metrics = self.problem.evaluate(self.problem.hohmann_decision).metadata
        contributions = objective.contributions(metrics)
        self.assertEqual(set(contributions), {"delta_v", "time", "energy", "risk"})
        for name, value in contributions.items():
            with self.subTest(term=name):
                self.assertAlmostEqual(value, 1.0, places=9)
        description = objective.describe()["terms"]
        self.assertEqual([term["unit"] for term in description], ["m/s", "s", "J/kg", "dimensionless"])
        self.assertTrue(all("unit weights to expose every term" in term["rationale"] for term in description))

    def test_default_weights_give_the_delta_v_objective(self) -> None:
        objective = transfer_objective(self.reference, rationale="minimum propellant: delta-v only")
        scalarized = ScalarizedProblem(self.problem, objective)
        for x in ([2600.0, 0.3], [2500.0, -0.1], self.problem.hohmann_decision):
            evaluation = scalarized.evaluate(x)
            self.assertAlmostEqual(evaluation.objective * self.reference.total_delta_v,
                                   evaluation.metadata["base_objective"], places=6)

    def test_time_weight_trades_delta_v_for_a_faster_transfer(self) -> None:
        objective = transfer_objective(self.reference, delta_v=1.0, time=1.0,
                                       rationale="equal importance of propellant and duration, for the test")
        result = DifferentialEvolution().minimize(ScalarizedProblem(self.problem, objective),
                                                  TerminationCriteria(max_evaluations=2000), seed=0)
        meta = result.best_evaluation.metadata
        # The Hohmann point is feasible with J = 2, so the minimum cannot be above it.
        self.assertLessEqual(result.f_best, 2.0 + 1e-9)
        self.assertLess(meta["transfer_time_s"], self.reference.transfer_time)
        self.assertGreater(meta["delta_v_total_m_s"], self.reference.total_delta_v)

    def test_validation(self) -> None:
        with self.assertRaises(ValueError):
            transfer_objective("reference", rationale="x")  # type: ignore[arg-type]
        for rationale in ("", "   ", None):
            with self.subTest(rationale=rationale), self.assertRaises(ValueError):
                transfer_objective(self.reference, rationale=rationale)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            transfer_objective(self.reference, delta_v=-1.0, rationale="negative weight")
        with self.assertRaises(ValueError):
            transfer_objective(self.reference, delta_v=0.0, rationale="no positive weight")


if __name__ == "__main__":
    unittest.main()
