"""Phase 4: mission constraints on transfer problems (RV-4.01, RV-4.02) and reference solutions.

Independent references:
  - the Hohmann closed form and Tsiolkovsky's equation, evaluated in the test;
  - for the trajectory minimum, the periapsis of the transfer ellipse from
    vis-viva and angular momentum, and a numerical RK45 propagation;
  - for the time-limited optimum, ``problems.references`` (closed-form
    transfer, boundary search and exhaustive grid), which shares nothing
    with the problem's evaluation chain; it is itself tested here against
    the problem's evaluator on a local grid.

Thresholds: feasibility tolerance 1e-9 of each constraint's scale (the
documented default); 1 m and 1e-6 relative for quantities that come from the
numerical integrator; eps_low = 1e-6 m/s for "not better than the reference".
"""
import math
import unittest

import numpy as np

from orbital_optimizer import EARTH_MU, hohmann_transfer
from orbital_optimizer.constants import EARTH_EQUATORIAL_RADIUS, STANDARD_GRAVITY
from orbital_optimizer.constraints import (
    MissionLimits, TransferSummary, evaluate_constraints, numerical_transfer_summary,
)
from orbital_optimizer.optimization import (
    DeathPenalty, DifferentialEvolution, EvaluationStatus, FeasibilityRules, StaticPenalty, TerminationCriteria,
)
from orbital_optimizer.problems import (
    ConstrainedTransferProblem, CoplanarMultiImpulseProblem, reach_delta_v, time_limited_reference,
    two_impulse_closed_form,
)
from orbital_optimizer.problems.references import CHECK_TOLERANCE_M_S
from orbital_optimizer.spacecraft import Spacecraft

LEO, GEO = 6_678_000.0, 42_164_000.0
EPS_LOW = 1e-6
HOHMANN = hohmann_transfer(LEO, GEO)
VEHICLE = Spacecraft(dry_mass_kg=500.0, propellant_mass_kg=2500.0, max_thrust_n=400.0, specific_impulse_s=320.0)
EXHAUST = 320.0 * STANDARD_GRAVITY
LEO_ALTITUDE = LEO - EARTH_EQUATORIAL_RADIUS


def constrained(limits: MissionLimits, vehicle: Spacecraft = VEHICLE, **kwargs) -> ConstrainedTransferProblem:
    return ConstrainedTransferProblem(CoplanarMultiImpulseProblem(LEO, GEO), vehicle, limits, **kwargs)


class InfeasibleSolutionTests(unittest.TestCase):
    """RV-4.01: a deliberately infeasible solution is detected, quantified and never reported as a solution."""

    def test_hohmann_over_a_delta_v_budget_is_rejected_with_its_magnitude(self) -> None:
        budget = HOHMANN.total_delta_v - 50.0
        problem = constrained(MissionLimits(max_delta_v_m_s=budget))
        evaluation = problem.evaluate(problem.hohmann_decision)
        self.assertIs(evaluation.status, EvaluationStatus.VALID)          # a real transfer...
        self.assertFalse(evaluation.is_feasible)                          # ...that breaks a limit
        self.assertAlmostEqual(evaluation.objective, HOHMANN.total_delta_v, delta=EPS_LOW)   # raw objective
        report = evaluation.metadata["constraint_report"]
        self.assertEqual((report["feasible"], report["violated"]), (False, ["max_delta_v"]))
        violated = problem.report(problem.hohmann_decision).by_name("max_delta_v")
        self.assertAlmostEqual(violated.violation, 50.0, delta=EPS_LOW)   # the amount, in m/s
        self.assertAlmostEqual(evaluation.constraint_violation, 50.0 / budget, delta=1e-12)
        self.assertEqual(set(evaluation.constraints), {"max_delta_v", "propellant"})
        self.assertEqual(evaluation.constraints["propellant"], 0.0)
        self.assertFalse(evaluation.metadata["feasible"])

    def test_optimizer_reports_that_no_feasible_transfer_exists(self) -> None:
        # Hohmann is the cheapest two-impulse transfer: with a smaller budget nothing is feasible.
        budget = HOHMANN.total_delta_v - 50.0
        result = DifferentialEvolution().minimize(constrained(MissionLimits(max_delta_v_m_s=budget)),
                                                  TerminationCriteria(2000), seed=0)
        self.assertTrue(result.found_valid)
        self.assertFalse(result.found_feasible)
        self.assertEqual((result.x_best, result.f_best), (None, math.inf))
        self.assertIn("no feasible point found", result.message)
        # The least violation cannot be below that of the Hohmann transfer itself.
        self.assertGreaterEqual(result.least_violation, 50.0 / budget - 1e-9)
        self.assertLess(result.least_violation, 50.0 / budget + 5e-3)
        self.assertGreater(result.least_violation_evaluation.objective, budget)

    def test_every_strategy_treats_the_infeasible_candidate_explicitly(self) -> None:
        limits = MissionLimits(max_delta_v_m_s=HOHMANN.total_delta_v - 50.0)
        infeasible = constrained(limits).evaluate([HOHMANN.first_delta_v, 0.0])
        feasible_problem = constrained(MissionLimits(max_delta_v_m_s=HOHMANN.total_delta_v + 500.0))
        feasible = feasible_problem.evaluate([2450.0, 0.02])            # 4066.6 m/s: within the larger budget
        self.assertTrue(feasible.is_feasible)
        self.assertGreater(feasible.objective, infeasible.objective)      # the feasible point is more expensive
        penalty = StaticPenalty(1e7, "test: large enough to dominate a 1 % violation")
        for strategy in (FeasibilityRules(), DeathPenalty(), penalty):
            with self.subTest(strategy=strategy.name):
                self.assertLess(strategy.key(feasible), strategy.key(infeasible))
        self.assertAlmostEqual(penalty.penalized(infeasible),
                               infeasible.objective + 1e7 * infeasible.constraint_violation**2, places=6)


class AltitudeConstraintTests(unittest.TestCase):
    """RV-4.02: the minimum altitude is checked over the whole trajectory, not at its ends."""

    def test_periapsis_below_the_threshold_with_valid_end_points(self) -> None:
        limits = MissionLimits(min_altitude_m=280_000.0)
        problem = constrained(limits)
        inward = problem.evaluate([2600.0, -0.3])           # first burn with an inward radial component
        # Periapsis of the transfer ellipse from vis-viva and angular momentum, independently.
        tangential = math.sqrt(EARTH_MU / LEO) + 2600.0 * math.cos(-0.3)
        radial = 2600.0 * math.sin(-0.3)
        a = -EARTH_MU / (tangential**2 + radial**2 - 2 * EARTH_MU / LEO)
        periapsis = a * (1 - math.sqrt(1 - (LEO * tangential) ** 2 / (EARTH_MU * a)))
        periapsis_altitude = periapsis - EARTH_EQUATORIAL_RADIUS
        # Both ends of the transfer arc satisfy the limit (299.9 km and GEO); only the path
        # between them does not (periapsis at 256.2 km).
        self.assertGreater(LEO_ALTITUDE, 280_000.0)
        self.assertGreater(GEO - EARTH_EQUATORIAL_RADIUS, 280_000.0)
        self.assertLess(periapsis_altitude, 280_000.0)
        self.assertFalse(inward.is_feasible)
        value = problem.report([2600.0, -0.3]).by_name("min_altitude")
        self.assertAlmostEqual(value.value, periapsis_altitude, delta=1e-3)
        self.assertAlmostEqual(value.violation, 280_000.0 - periapsis_altitude, delta=1e-3)
        self.assertAlmostEqual(inward.constraint_violation,
                               (280_000.0 - periapsis_altitude) / (EARTH_EQUATORIAL_RADIUS + 280_000.0), delta=1e-9)
        # An end-point-only check would have accepted it.
        end_points_only = TransferSummary(tuple(inward.metadata["burn_delta_v_m_s"]),
                                          inward.metadata["transfer_time_s"], min(LEO, GEO))
        self.assertTrue(evaluate_constraints(end_points_only, VEHICLE, limits).feasible)

    def test_mirror_transfer_that_leaves_outward_is_feasible(self) -> None:
        problem = constrained(MissionLimits(min_altitude_m=280_000.0))
        outward = problem.evaluate([2600.0, 0.3])
        self.assertTrue(outward.is_feasible)
        self.assertAlmostEqual(outward.objective, problem.evaluate([2600.0, -0.3]).objective, delta=1e-6)

    def test_numerical_propagation_confirms_the_violation(self) -> None:
        limits = MissionLimits(min_altitude_m=280_000.0)
        problem = constrained(limits)
        plan = problem.plan([2600.0, -0.3])
        analytic = problem.report([2600.0, -0.3]).by_name("min_altitude")
        numerical = evaluate_constraints(numerical_transfer_summary(plan.initial_state, plan.impulses(), plan.mu),
                                         VEHICLE, limits).by_name("min_altitude")
        self.assertFalse(numerical.satisfied)
        self.assertAlmostEqual(numerical.value, analytic.value, delta=1.0)        # integrator accuracy

    def test_limit_equal_to_the_departure_altitude_is_met_by_hohmann(self) -> None:
        # Boundary case: the Hohmann transfer never goes below the initial orbit.
        problem = constrained(MissionLimits(min_altitude_m=LEO_ALTITUDE))
        evaluation = problem.evaluate(problem.hohmann_decision)
        self.assertTrue(evaluation.is_feasible, evaluation.metadata["constraint_report"])
        self.assertLessEqual(abs(problem.report(problem.hohmann_decision).by_name("min_altitude").g), 1e-2)
        self.assertFalse(problem.evaluate([2600.0, -0.01]).is_feasible)


class CatalogueOnTransfersTests(unittest.TestCase):
    def test_boundary_budget_equal_to_hohmann_is_feasible(self) -> None:
        problem = constrained(MissionLimits(max_delta_v_m_s=HOHMANN.total_delta_v, max_duration_s=HOHMANN.transfer_time))
        evaluation = problem.evaluate(problem.hohmann_decision)
        self.assertTrue(evaluation.is_feasible, evaluation.metadata["constraint_report"])
        self.assertEqual(evaluation.constraint_violation, 0.0)

    def test_propellant_follows_tsiolkovsky(self) -> None:
        problem = constrained(MissionLimits())
        evaluation = problem.evaluate(problem.hohmann_decision)
        required = 3000.0 * (1 - math.exp(-HOHMANN.total_delta_v / EXHAUST))
        self.assertAlmostEqual(evaluation.metadata["propellant_required_kg"], required, delta=1e-6)
        self.assertAlmostEqual(evaluation.metadata["final_mass_kg"], 3000.0 - required, delta=1e-6)
        after_first = 3000.0 * math.exp(-HOHMANN.first_delta_v / EXHAUST)
        np.testing.assert_allclose(evaluation.metadata["propellant_per_burn_kg"],
                                   [3000.0 - after_first, after_first - (3000.0 - required)], rtol=1e-9)
        self.assertTrue(evaluation.is_feasible)                           # 2132 kg needed, 2500 kg on board
        short = Spacecraft(500.0, 2000.0, 400.0, 320.0)                   # wet 2500 kg: needs 1777 kg... of 2000
        self.assertTrue(constrained(MissionLimits(), short).evaluate(problem.hohmann_decision).is_feasible)
        dry = Spacecraft(1500.0, 1500.0, 400.0, 320.0)                    # wet 3000 kg: needs 2132 kg of 1500
        starved = constrained(MissionLimits(), dry).evaluate(problem.hohmann_decision)
        self.assertFalse(starved.is_feasible)
        self.assertAlmostEqual(starved.constraint_violation, (required - 1500.0) / 3000.0, delta=1e-9)

    def test_thrust_and_mass_are_configuration_constraints(self) -> None:
        within = constrained(MissionLimits(max_thrust_n=500.0, max_wet_mass_kg=3500.0), commanded_thrust_n=400.0)
        self.assertTrue(within.configuration_report().feasible)
        evaluation = within.evaluate(within.hohmann_decision)
        self.assertTrue(evaluation.is_feasible)
        self.assertEqual(evaluation.metadata["burn_thrust_n"], 400.0)
        np.testing.assert_allclose(evaluation.metadata["burn_durations_s"],
                                   evaluation.metadata["propellant_per_burn_kg"] * EXHAUST / 400.0, rtol=1e-12)
        over = constrained(MissionLimits(max_thrust_n=300.0, max_wet_mass_kg=2800.0), commanded_thrust_n=350.0)
        report = over.configuration_report()
        self.assertEqual(report.violated, ("max_wet_mass", "max_thrust"))
        self.assertAlmostEqual(report.total_violation, 200.0 / 2800.0 + 50.0 / 300.0, places=12)
        # Every candidate inherits the configuration violation: none can be feasible.
        for x in (over.hohmann_decision, [2600.0, 0.3]):
            self.assertFalse(over.evaluate(x).is_feasible)
        self.assertAlmostEqual(over.evaluate(over.hohmann_decision).constraint_violation,
                               200.0 / 2800.0 + 50.0 / 300.0, places=12)
        impulsive = constrained(MissionLimits(max_thrust_n=300.0))
        self.assertIn("max_thrust", impulsive.configuration_report().not_applicable)
        self.assertIsNone(impulsive.evaluate(impulsive.hohmann_decision).metadata["burn_durations_s"])
        self.assertIn("max_thrust", impulsive.evaluate(impulsive.hohmann_decision).metadata["constraint_report"]["not_applicable"])

    def test_simultaneous_violations_add_up(self) -> None:
        # The candidate costs 4566.9 m/s, takes 13 251.8 s and descends to 256.2 km: three limits broken.
        limits = MissionLimits(max_delta_v_m_s=4500.0, max_duration_s=13_000.0, min_altitude_m=280_000.0)
        problem = constrained(limits)
        evaluation = problem.evaluate([2600.0, -0.3])
        report = problem.report([2600.0, -0.3])
        self.assertEqual(set(report.violated), {"max_delta_v", "max_duration", "min_altitude"})
        self.assertAlmostEqual(evaluation.constraint_violation,
                               sum(report.by_name(name).normalized_violation for name in report.violated), places=12)
        self.assertEqual(set(evaluation.constraints), {"max_delta_v", "max_duration", "min_altitude", "propellant"})

    def test_invalid_candidates_stay_invalid_and_wrapper_contracts(self) -> None:
        problem = constrained(MissionLimits(max_delta_v_m_s=4000.0))
        invalid = problem.evaluate([0.0, 0.0])
        self.assertIs(invalid.status, EvaluationStatus.INVALID_DOMAIN)
        self.assertIsNone(invalid.constraint_violation)
        with self.assertRaises(ValueError):
            problem.report([0.0, 0.0])
        self.assertEqual(problem.name, "coplanar-2-impulse:constrained")
        self.assertIs(problem.bounds, problem.problem.bounds)
        self.assertEqual(problem.reference.total_delta_v, HOHMANN.total_delta_v)         # delegated attribute
        self.assertEqual(len(problem.plan(problem.hohmann_decision).burn_times), 2)
        self.assertIsInstance(problem.constraint_handling, FeasibilityRules)
        description = problem.describe()
        self.assertEqual(description["limits"]["max_delta_v_m_s"], 4000.0)
        self.assertEqual(description["spacecraft"]["wet_mass_kg"], 3000.0)
        self.assertEqual(description["constraint_handling"]["name"], "feasibility-rules")
        with self.assertRaises(AttributeError):
            problem.does_not_exist
        base = CoplanarMultiImpulseProblem(LEO, GEO)
        for arguments in ((object(), VEHICLE, MissionLimits()), (base, "vehicle", MissionLimits()), (base, VEHICLE, {})):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                ConstrainedTransferProblem(*arguments)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ConstrainedTransferProblem(base, VEHICLE, MissionLimits(), handling="deb")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ConstrainedTransferProblem(base, VEHICLE, MissionLimits(), commanded_thrust_n=-5.0)


class ReferenceSolutionTests(unittest.TestCase):
    def test_closed_form_agrees_with_the_problem_evaluator(self) -> None:
        problem = CoplanarMultiImpulseProblem(LEO, GEO)
        rng = np.random.Generator(np.random.PCG64(11))
        points = problem.bounds.from_unit(rng.random((2500, 2)))
        closed = two_impulse_closed_form(LEO, GEO, EARTH_MU, points[:, 0], points[:, 1])
        compared = 0
        for index, point in enumerate(points):
            evaluation = problem.evaluate(point)
            if not evaluation.is_valid or not closed["valid"][index]:
                continue                 # the two domain tests may differ within rounding at the tangency
            compared += 1
            meta = evaluation.metadata
            self.assertAlmostEqual(closed["delta_v_total_m_s"][index], evaluation.objective, delta=1e-6)
            self.assertAlmostEqual(closed["minimum_radius_m"][index], meta["minimum_radius_m"], delta=1e-3)
            if 1 + (meta["final_arc_semi_major_axis_m"] * (1 - meta["final_arc_eccentricity"] ** 2) / GEO - 1) \
                    / meta["final_arc_eccentricity"] >= 1e-2:
                # Away from a tangent arrival, where the crossing time is ill-conditioned.
                self.assertAlmostEqual(closed["transfer_time_s"][index] / meta["transfer_time_s"], 1.0, delta=1e-9)
        self.assertGreater(compared, 80)
        with self.assertRaises(ValueError):
            two_impulse_closed_form(GEO, LEO, EARTH_MU, 100.0, 0.0)

    def test_reach_delta_v_is_the_tangency_boundary(self) -> None:
        self.assertAlmostEqual(float(reach_delta_v(LEO, GEO, EARTH_MU, 0.0)), HOHMANN.first_delta_v, delta=1e-9)
        for alpha in (-0.6, 0.2, 0.9):
            boundary = float(reach_delta_v(LEO, GEO, EARTH_MU, alpha))
            above = two_impulse_closed_form(LEO, GEO, EARTH_MU, boundary * (1 + 1e-9), alpha)
            below = two_impulse_closed_form(LEO, GEO, EARTH_MU, boundary * (1 - 1e-9), alpha)
            self.assertTrue(bool(above["valid"]) and not bool(below["valid"]), alpha)

    def test_reference_without_a_binding_limit_is_hohmann(self) -> None:
        for factor in (1.0, 1.5):
            reference = time_limited_reference(LEO, GEO, factor * HOHMANN.transfer_time)
            self.assertFalse(reference.active)
            self.assertEqual((reference.delta_v_total_m_s, reference.alpha1_rad), (HOHMANN.total_delta_v, 0.0))

    def test_binding_limit_gives_a_feasible_boundary_point(self) -> None:
        previous = HOHMANN.total_delta_v
        for fraction in (0.95, 0.8, 0.6):
            limit = fraction * HOHMANN.transfer_time
            reference = time_limited_reference(LEO, GEO, limit)
            with self.subTest(fraction=fraction):
                self.assertTrue(reference.active)
                self.assertLessEqual(reference.transfer_time_s, limit)                    # feasible
                self.assertAlmostEqual(reference.transfer_time_s / limit, 1.0, delta=1e-9)  # and on the boundary
                self.assertGreater(reference.delta_v_total_m_s, previous)                 # a tighter limit costs more
                self.assertGreaterEqual(reference.check_grid_best_m_s, reference.delta_v_total_m_s - CHECK_TOLERANCE_M_S)
                self.assertEqual(reference.to_dict()["active"], True)
                previous = reference.delta_v_total_m_s

    def test_no_evaluator_point_near_the_reference_is_cheaper(self) -> None:
        """Local exhaustive search with the problem's own evaluator, around the reference.

        Window of +-2 % in dv1 and +-0.1 rad in alpha, 81 x 81 nodes: no node
        that meets the time limit may cost less than the reference by more
        than eps_low, and the best node must be within the cost variation of
        one cell (taken from the grid itself).
        """
        limit = 0.8 * HOHMANN.transfer_time
        reference = time_limited_reference(LEO, GEO, limit)
        problem = CoplanarMultiImpulseProblem(LEO, GEO)
        feasible_costs = []
        for dv in np.linspace(0.98, 1.02, 81) * reference.delta_v1_m_s:
            for alpha in reference.alpha1_rad + np.linspace(-0.1, 0.1, 81):
                evaluation = problem.evaluate([dv, alpha])
                if evaluation.is_valid and evaluation.metadata["transfer_time_s"] <= limit:
                    feasible_costs.append(evaluation.objective)
        self.assertGreater(len(feasible_costs), 500)
        self.assertGreaterEqual(min(feasible_costs), reference.delta_v_total_m_s - EPS_LOW)
        self.assertLess(min(feasible_costs) - reference.delta_v_total_m_s, 5.0)
        evaluation = problem.evaluate([reference.delta_v1_m_s, reference.alpha1_rad])
        self.assertAlmostEqual(evaluation.objective, reference.delta_v_total_m_s, delta=EPS_LOW)

    def test_invalid_and_infeasible_requests(self) -> None:
        with self.assertRaises(ValueError):               # faster than any transfer within the decision bounds
            time_limited_reference(LEO, GEO, 0.3 * HOHMANN.transfer_time)
        for arguments in ((GEO, LEO, 1000.0), (LEO, GEO, 0.0), (LEO, GEO, -5.0), (LEO, GEO, math.nan)):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                time_limited_reference(*arguments)
        with self.assertRaises(ValueError):
            time_limited_reference(LEO, GEO, 0.8 * HOHMANN.transfer_time, directions=2)


class ConstrainedOptimizationTests(unittest.TestCase):
    def test_time_limited_run_is_feasible_and_not_better_than_the_reference(self) -> None:
        limit = 0.8 * HOHMANN.transfer_time
        reference = time_limited_reference(LEO, GEO, limit)
        problem = constrained(MissionLimits(max_duration_s=limit))
        result = DifferentialEvolution().minimize(problem, TerminationCriteria(2000), seed=0)
        self.assertTrue(result.found_feasible)
        meta = result.best_evaluation.metadata
        self.assertLessEqual(meta["transfer_time_s"], limit * (1 + 1e-9))        # within the feasibility tolerance
        # The tolerance admits a duration 1e-9 above the limit; the cost is not measurably lower for it.
        self.assertGreaterEqual(result.f_best, reference.delta_v_total_m_s - 1e-3)
        self.assertGreater(result.f_best, HOHMANN.total_delta_v + 100.0)         # the limit costs about 207 m/s
        self.assertTrue(meta["feasible"])

    def test_unconstrained_limits_reduce_to_the_base_problem(self) -> None:
        base = CoplanarMultiImpulseProblem(LEO, GEO)
        # Delta-v capacity 43 km/s: no candidate of the decision box runs out of propellant.
        tanker = Spacecraft(dry_mass_kg=1.0, propellant_mass_kg=1.0e6, max_thrust_n=400.0, specific_impulse_s=320.0)
        self.assertGreater(tanker.delta_v_capacity_m_s, 40_000.0)
        wrapped = ConstrainedTransferProblem(CoplanarMultiImpulseProblem(LEO, GEO), tanker, MissionLimits())
        first = DifferentialEvolution().minimize(base, TerminationCriteria(600), seed=4)
        second = DifferentialEvolution().minimize(wrapped, TerminationCriteria(600), seed=4)
        # Every valid candidate is feasible, so the search is the same as without constraints.
        self.assertEqual(first.f_best, second.f_best)
        np.testing.assert_array_equal(first.history.best_so_far, second.history.best_so_far)
        self.assertTrue(second.found_feasible)
        self.assertIsNone(second.least_violation)


if __name__ == "__main__":
    unittest.main()
