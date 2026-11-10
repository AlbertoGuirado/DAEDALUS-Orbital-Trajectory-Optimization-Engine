"""Phase 4: spacecraft model and mission constraint catalogue (RF-4.01 to RF-4.08).

Every expected value is computed in the test from the defining expression
(Tsiolkovsky, sums, differences), never read back from the code under test.
The feasibility tolerance is the documented one: 1e-9 of each constraint's
scale. The tolerance of the numerical verification is the accuracy of the
integrator (1 m, the Phase 2 arrival criterion), not of the search.
"""
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose

from orbital_optimizer import EARTH_MU, OrbitalElements, elements_to_state, hohmann_transfer, propellant_mass
from orbital_optimizer.constants import EARTH_EQUATORIAL_RADIUS, STANDARD_GRAVITY
from orbital_optimizer.constraints import (
    CONSTRAINT_NAMES, ConstraintReport, ConstraintValue, MissionLimits, TransferSummary, arc_minimum_radius,
    evaluate_constraints, limits_from_mapping, numerical_transfer_summary, trajectory_minimum_radius,
)
from orbital_optimizer.constraints.verification import _parabola_minimum
from orbital_optimizer.orbital.kepler import radius_extrema_on_arc
from orbital_optimizer.simulation.impulsive import propagate_impulsive_arcs
from orbital_optimizer.simulation.propagator import propagate
from orbital_optimizer.spacecraft import Spacecraft, spacecraft_from_mapping

LEO, GEO = 6_678_000.0, 42_164_000.0
VEHICLE = Spacecraft(dry_mass_kg=500.0, propellant_mass_kg=1500.0, max_thrust_n=400.0, specific_impulse_s=320.0)
EXHAUST = 320.0 * STANDARD_GRAVITY


class SpacecraftTests(unittest.TestCase):
    def test_masses_and_rocket_equation(self) -> None:
        self.assertEqual(VEHICLE.wet_mass_kg, 2000.0)
        self.assertEqual(VEHICLE.exhaust_velocity_m_s, EXHAUST)
        self.assertAlmostEqual(VEHICLE.delta_v_capacity_m_s, EXHAUST * math.log(2000.0 / 500.0), places=9)
        # Burning the whole capacity consumes exactly the propellant on board.
        self.assertAlmostEqual(VEHICLE.propellant_required_kg([VEHICLE.delta_v_capacity_m_s]), 1500.0, places=9)
        # One burn: the Phase 2 function with the wet mass as initial mass.
        self.assertAlmostEqual(VEHICLE.propellant_required_kg([1000.0]), propellant_mass(2000.0, 1000.0, 320.0), places=10)

    def test_sequential_burns(self) -> None:
        burns = [2425.769, 1466.839]
        per_burn = VEHICLE.propellant_per_burn_kg(burns)
        after_first = 2000.0 * math.exp(-burns[0] / EXHAUST)
        after_second = after_first * math.exp(-burns[1] / EXHAUST)
        assert_allclose(per_burn, [2000.0 - after_first, after_first - after_second], rtol=1e-13)
        # The total only depends on the sum of the burns (same specific impulse).
        self.assertAlmostEqual(per_burn.sum(), VEHICLE.propellant_required_kg(burns), places=9)
        self.assertAlmostEqual(VEHICLE.propellant_required_kg(burns), VEHICLE.propellant_required_kg(burns[::-1]), places=9)
        self.assertAlmostEqual(VEHICLE.propellant_required_kg(burns), 2000.0 * (1 - math.exp(-sum(burns) / EXHAUST)), places=9)
        self.assertEqual(VEHICLE.propellant_per_burn_kg([]).shape, (0,))
        self.assertEqual(VEHICLE.propellant_required_kg([]), 0.0)
        # The request is not capped by what the tanks hold: that is a constraint, checked elsewhere.
        self.assertGreater(VEHICLE.propellant_required_kg([6000.0]), VEHICLE.propellant_mass_kg)

    def test_burn_durations_follow_the_mass_flow(self) -> None:
        burns = [300.0, 150.0]
        durations = VEHICLE.burn_durations_s(burns)
        # mass flow = F / v_e, so propellant = F t / v_e.
        assert_allclose(durations * 400.0 / EXHAUST, VEHICLE.propellant_per_burn_kg(burns), rtol=1e-13)
        assert_allclose(VEHICLE.burn_durations_s(burns, 200.0), 2 * durations, rtol=1e-13)
        with self.assertRaises(ValueError):
            VEHICLE.burn_durations_s(burns, 0.0)

    def test_validation_and_mapping(self) -> None:
        valid = {"dry_mass_kg": 500.0, "propellant_mass_kg": 1500.0, "max_thrust_n": 400.0, "specific_impulse_s": 320.0}
        self.assertEqual(spacecraft_from_mapping(valid), VEHICLE)
        self.assertEqual(Spacecraft(**(valid | {"propellant_mass_kg": 0.0})).delta_v_capacity_m_s, 0.0)
        for change in ({"dry_mass_kg": 0.0}, {"dry_mass_kg": -1.0}, {"propellant_mass_kg": -1.0}, {"max_thrust_n": 0.0},
                       {"specific_impulse_s": 0.0}, {"dry_mass_kg": math.nan}, {"dry_mass_kg": math.inf},
                       {"dry_mass_kg": [1.0, 2.0]}, {"specific_impulse_s": "320"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                Spacecraft(**(valid | change))
        for bad in ({**valid, "area_m2": 1.0}, {"dry_mass_kg": 500.0}, [500.0, 1500.0, 400.0, 320.0]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                spacecraft_from_mapping(bad)
        for burns in ([-1.0], [[1.0, 2.0]], [math.nan]):
            with self.subTest(burns=burns), self.assertRaises(ValueError):
                VEHICLE.propellant_per_burn_kg(burns)
        description = VEHICLE.describe()
        self.assertEqual((description["wet_mass_kg"], description["exhaust_velocity_m_s"]), (2000.0, EXHAUST))


class ConstraintModelTests(unittest.TestCase):
    def value(self, g: float, tolerance: float = 0.5, **kwargs) -> ConstraintValue:
        return ConstraintValue("limit", "value <= limit", "m/s", 100.0 + g, 100.0, g, 100.0, tolerance, **kwargs)

    def test_sign_convention_violation_and_margin(self) -> None:
        inside = self.value(-20.0)
        self.assertTrue(inside.satisfied)
        self.assertEqual((inside.violation, inside.normalized_violation, inside.margin), (0.0, 0.0, 20.0))
        outside = self.value(30.0)
        self.assertFalse(outside.satisfied)
        self.assertEqual((outside.violation, outside.normalized_violation, outside.margin), (30.0, 0.3, -30.0))

    def test_tolerance_edges(self) -> None:
        self.assertTrue(self.value(0.0).satisfied)                       # exactly on the limit
        within = self.value(0.5)                                         # exceeds by exactly the tolerance
        self.assertTrue(within.satisfied)
        self.assertEqual(within.violation, 0.0)
        beyond = self.value(0.5000001)
        self.assertFalse(beyond.satisfied)
        self.assertEqual(beyond.violation, 0.5000001)                    # the full residual, not g - tolerance
        self.assertFalse(self.value(1e-12, tolerance=0.0).satisfied)     # zero tolerance is strict

    def test_non_finite_residuals_are_never_satisfied(self) -> None:
        for g in (math.nan, math.inf):
            with self.subTest(g=g):
                value = self.value(g)
                self.assertFalse(value.satisfied)
                self.assertEqual(value.violation, math.inf)
                self.assertEqual(value.normalized_violation, math.inf)
        self.assertTrue(self.value(-math.inf).satisfied)

    def test_report_aggregates(self) -> None:
        hard_ok = ConstraintValue("a", "a", "kg", 1.0, 2.0, -1.0, 2.0, 0.0)
        hard_bad = ConstraintValue("b", "b", "s", 12.0, 10.0, 2.0, 10.0, 0.0)
        hard_worse = ConstraintValue("c", "c", "m", 5.0, 4.0, 1.0, 4.0, 0.0)
        soft_bad = ConstraintValue("d", "d", "N", 9.0, 6.0, 3.0, 6.0, 0.0, hard=False)
        report = ConstraintReport((hard_ok, hard_bad, hard_worse, soft_bad), {"e": "no thrust level"})
        self.assertFalse(report.feasible)
        self.assertAlmostEqual(report.total_violation, 2.0 / 10.0 + 1.0 / 4.0, places=15)   # hard only
        self.assertAlmostEqual(report.soft_violation, 3.0 / 6.0, places=15)
        self.assertEqual(report.violated, ("b", "c"))
        self.assertEqual(report.normalized(), {"a": 0.0, "b": 0.2, "c": 0.25})
        self.assertIs(report.by_name("c"), hard_worse)
        with self.assertRaises(KeyError):
            report.by_name("z")
        record = report.to_dict()
        self.assertEqual((record["feasible"], record["violated"], record["not_applicable"]),
                         (False, ["b", "c"], {"e": "no thrust level"}))
        self.assertEqual(record["constraints"][1]["g"], 2.0)
        # A soft violation alone does not make the candidate infeasible.
        self.assertTrue(ConstraintReport((hard_ok, soft_bad)).feasible)
        self.assertTrue(ConstraintReport(()).feasible)
        self.assertEqual(ConstraintReport(()).total_violation, 0.0)

    def test_validation(self) -> None:
        valid = {"name": "a", "description": "a", "unit": "m", "value": 1.0, "limit": 2.0, "g": -1.0,
                 "scale": 2.0, "tolerance": 0.0}
        for change in ({"name": ""}, {"unit": ""}, {"scale": 0.0}, {"scale": -1.0}, {"scale": math.inf},
                       {"tolerance": -1.0}, {"tolerance": math.nan}, {"g": "1"}, {"value": True}, {"hard": 1}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                ConstraintValue(**(valid | change))
        value = ConstraintValue(**valid)
        with self.assertRaises(ValueError):
            ConstraintReport((value, value))                             # duplicated name
        with self.assertRaises(ValueError):
            ConstraintReport((value,), {"a": "also not applicable"})
        with self.assertRaises(ValueError):
            ConstraintReport(("a",))  # type: ignore[arg-type]


class CatalogueTests(unittest.TestCase):
    """Each residual against its mathematical expression."""

    def setUp(self) -> None:
        self.summary = TransferSummary((2500.0, 1500.0), 20_000.0, 6_500_000.0, thrust_n=350.0)
        self.limits = MissionLimits(max_wet_mass_kg=2200.0, max_thrust_n=380.0, max_delta_v_m_s=3900.0,
                                    max_duration_s=19_000.0, min_altitude_m=200_000.0)

    def test_every_residual_matches_its_expression(self) -> None:
        report = evaluate_constraints(self.summary, VEHICLE, self.limits)
        self.assertEqual([value.name for value in report.values], list(CONSTRAINT_NAMES))
        required = 2000.0 * (1 - math.exp(-4000.0 / EXHAUST))
        altitude = 6_500_000.0 - EARTH_EQUATORIAL_RADIUS
        expected = {
            # name: (value, limit, g, scale, unit)
            "max_wet_mass": (2000.0, 2200.0, -200.0, 2200.0, "kg"),
            "max_thrust": (350.0, 380.0, -30.0, 380.0, "N"),
            "max_delta_v": (4000.0, 3900.0, 100.0, 3900.0, "m/s"),
            "max_duration": (20_000.0, 19_000.0, 1000.0, 19_000.0, "s"),
            "min_altitude": (altitude, 200_000.0, 200_000.0 - altitude, EARTH_EQUATORIAL_RADIUS + 200_000.0, "m"),
            "propellant": (required, 1500.0, required - 1500.0, 2000.0, "kg"),
        }
        for name, (value, limit, g, scale, unit) in expected.items():
            with self.subTest(constraint=name):
                got = report.by_name(name)
                self.assertAlmostEqual(got.value, value, places=9)
                self.assertEqual((got.limit, got.scale, got.unit), (limit, scale, unit))
                self.assertAlmostEqual(got.g, g, places=9)
                self.assertEqual(got.tolerance, 1e-9 * scale)
                self.assertTrue(got.hard)
        # altitude 121 863 m < 200 km and 1441.5 kg of propellant < 1500 kg on board.
        self.assertEqual(report.violated, ("max_delta_v", "max_duration", "min_altitude"))
        self.assertAlmostEqual(report.total_violation,
                               100.0 / 3900.0 + 1000.0 / 19_000.0 + (200_000.0 - altitude) / 6_578_137.0, places=12)
        self.assertFalse(report.feasible)

    def test_thrust_limit_is_the_lower_of_engine_and_mission(self) -> None:
        mission_lower = evaluate_constraints(self.summary, VEHICLE, MissionLimits(max_thrust_n=300.0))
        self.assertEqual((mission_lower.by_name("max_thrust").limit, mission_lower.by_name("max_thrust").g), (300.0, 50.0))
        engine_lower = evaluate_constraints(self.summary, VEHICLE, MissionLimits(max_thrust_n=900.0))
        self.assertEqual(engine_lower.by_name("max_thrust").limit, 400.0)
        no_mission_limit = evaluate_constraints(self.summary, VEHICLE, MissionLimits())
        self.assertEqual(no_mission_limit.by_name("max_thrust").limit, 400.0)
        over = TransferSummary((100.0,), 10.0, 7.0e6, thrust_n=450.0)
        self.assertEqual(evaluate_constraints(over, VEHICLE, MissionLimits()).violated, ("max_thrust",))

    def test_thrust_without_a_thrust_level_is_not_applicable(self) -> None:
        impulsive = TransferSummary((2500.0, 1500.0), 20_000.0, 6_500_000.0)
        report = evaluate_constraints(impulsive, VEHICLE, self.limits)
        self.assertNotIn("max_thrust", [value.name for value in report.values])
        self.assertIn("no thrust level", report.not_applicable["max_thrust"])
        self.assertIn("max_thrust", report.to_dict()["not_applicable"])

    def test_only_imposed_limits_are_evaluated_and_propellant_always_is(self) -> None:
        report = evaluate_constraints(TransferSummary((100.0,), 50.0, 7.0e6), VEHICLE, MissionLimits())
        self.assertEqual([value.name for value in report.values], ["propellant"])
        self.assertTrue(report.feasible)

    def test_exact_boundaries_are_feasible_and_tolerance_edges(self) -> None:
        summary = TransferSummary((3900.0,), 19_000.0, EARTH_EQUATORIAL_RADIUS + 200_000.0, thrust_n=380.0)
        vehicle = Spacecraft(500.0, 1700.0, 400.0, 320.0)
        limits = MissionLimits(max_wet_mass_kg=2200.0, max_thrust_n=380.0, max_delta_v_m_s=3900.0,
                               max_duration_s=19_000.0, min_altitude_m=200_000.0)
        report = evaluate_constraints(summary, vehicle, limits)
        self.assertTrue(report.feasible, report.violated)
        for name in ("max_wet_mass", "max_thrust", "max_delta_v", "max_duration", "min_altitude"):
            self.assertEqual(report.by_name(name).g, 0.0, name)
        self.assertEqual(report.total_violation, 0.0)
        # Half the tolerance above the limit is still feasible; twice the tolerance is not.
        inside = TransferSummary((3900.0 * (1 + 0.5e-9),), 19_000.0, 7.0e6)
        outside = TransferSummary((3900.0 * (1 + 2e-9),), 19_000.0, 7.0e6)
        self.assertTrue(evaluate_constraints(inside, vehicle, limits).by_name("max_delta_v").satisfied)
        value = evaluate_constraints(outside, vehicle, limits).by_name("max_delta_v")
        self.assertFalse(value.satisfied)
        self.assertAlmostEqual(value.violation, 3900.0 * 2e-9, delta=1e-12)
        strict = MissionLimits(max_delta_v_m_s=3900.0, relative_tolerance=0.0)
        self.assertFalse(evaluate_constraints(inside, vehicle, strict).by_name("max_delta_v").satisfied)

    def test_soft_constraints_do_not_decide_feasibility(self) -> None:
        soft = MissionLimits(max_delta_v_m_s=3900.0, max_duration_s=19_000.0, soft=frozenset({"max_duration"}))
        vehicle = Spacecraft(500.0, 5000.0, 400.0, 320.0)
        slow = evaluate_constraints(TransferSummary((3000.0,), 25_000.0, 7.0e6), vehicle, soft)
        self.assertTrue(slow.feasible)
        self.assertFalse(slow.by_name("max_duration").hard)
        self.assertAlmostEqual(slow.soft_violation, 6000.0 / 19_000.0, places=15)
        self.assertEqual(slow.total_violation, 0.0)
        expensive = evaluate_constraints(TransferSummary((4000.0,), 25_000.0, 7.0e6), vehicle, soft)
        self.assertEqual(expensive.violated, ("max_delta_v",))

    def test_non_finite_quantities_are_never_feasible(self) -> None:
        limits = MissionLimits(max_delta_v_m_s=3900.0, max_duration_s=19_000.0, min_altitude_m=200_000.0)
        for summary, names in ((TransferSummary((math.nan,), 100.0, 7.0e6), {"max_delta_v", "propellant"}),
                               (TransferSummary((100.0,), math.nan, 7.0e6), {"max_duration"}),
                               (TransferSummary((100.0,), 100.0, math.nan), {"min_altitude"}),
                               (TransferSummary((math.inf,), 100.0, 7.0e6), {"max_delta_v", "propellant"})):
            with self.subTest(summary=summary):
                report = evaluate_constraints(summary, VEHICLE, limits)
                self.assertFalse(report.feasible)
                self.assertEqual(set(report.violated), names)
                self.assertEqual(report.total_violation, math.inf)

    def test_limits_and_summary_validation(self) -> None:
        for kwargs in ({"max_wet_mass_kg": 0.0}, {"max_thrust_n": -1.0}, {"max_delta_v_m_s": math.nan},
                       {"max_duration_s": math.inf}, {"min_altitude_m": -1.0}, {"body_radius_m": 0.0},
                       {"relative_tolerance": -1e-9}, {"soft": frozenset({"max_speed"})}, {"soft": "max_delta_v"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                MissionLimits(**kwargs)
        self.assertEqual(MissionLimits(min_altitude_m=0.0).min_altitude_m, 0.0)
        limits = limits_from_mapping({"max_delta_v_m_s": 4000.0, "min_altitude_m": 150_000.0, "soft": ["max_delta_v"]})
        self.assertEqual((limits.max_delta_v_m_s, limits.soft, limits.body_radius_m),
                         (4000.0, frozenset({"max_delta_v"}), EARTH_EQUATORIAL_RADIUS))
        self.assertEqual(limits.describe()["relative_tolerance"], 1e-9)
        for bad in ({"max_speed": 1.0}, [1.0], {"soft": "max_delta_v"}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                limits_from_mapping(bad)
        for arguments in (((-1.0,), 10.0, 7e6), ((1.0,), -1.0, 7e6), ((1.0,), 10.0, 0.0), (("a",), 10.0, 7e6),
                          ((1.0,), "10", 7e6)):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                TransferSummary(*arguments)
        with self.assertRaises(ValueError):
            TransferSummary((1.0,), 10.0, 7e6, thrust_n=0.0)
        self.assertEqual(TransferSummary((1.5, 2.5), 10.0, 7e6).total_delta_v_m_s, 4.0)
        metadata = {"burn_delta_v_m_s": np.array([1.0, 2.0]), "transfer_time_s": 5.0, "minimum_radius_m": 7e6}
        self.assertEqual(TransferSummary.from_metadata(metadata, 100.0),
                         TransferSummary((1.0, 2.0), 5.0, 7e6, 100.0))
        with self.assertRaises(ValueError):
            TransferSummary.from_metadata({"transfer_time_s": 5.0})
        with self.assertRaises(ValueError):
            evaluate_constraints("summary", VEHICLE, MissionLimits())  # type: ignore[arg-type]


class NumericalVerificationTests(unittest.TestCase):
    def test_parabola_vertex(self) -> None:
        # y = 3 + 2 (t - 0.3)^2 sampled at t = -1, 0, 1: the vertex value is 3.
        samples = np.array([3 + 2 * (t - 0.3) ** 2 for t in (-1.0, 0.0, 1.0)])
        self.assertAlmostEqual(_parabola_minimum(samples), 3.0, places=13)
        self.assertEqual(_parabola_minimum(np.array([2.0, 1.0, 2.0])), 1.0)          # symmetric
        self.assertEqual(_parabola_minimum(np.array([1.0, 1.0, 1.0])), 1.0)          # flat: no curvature

    def test_arc_minimum_matches_the_exact_periapsis(self) -> None:
        # Inbound start: the periapsis lies inside the arc, between integration steps.
        elements = OrbitalElements(1.6e7, 0.55, 0.4, 1.0, 2.0, 5.0)
        state = elements_to_state(elements)
        duration = 9000.0
        exact, _ = radius_extrema_on_arc(elements.true_anomaly, duration, 1.6e7, 0.55)
        self.assertAlmostEqual(exact, 1.6e7 * 0.45, delta=1e-3)
        arc = propagate(state, duration)
        coarse = float(np.min(np.linalg.norm(arc.states[:, :3], axis=1)))
        refined = arc_minimum_radius(arc)
        # The accepted steps (60 s apart) miss the periapsis by far more than the refined estimate.
        self.assertGreater(coarse - exact, 10.0)
        self.assertAlmostEqual(refined, exact, delta=1.0)
        self.assertLess(abs(refined - exact), 0.01 * (coarse - exact))
        with self.assertRaises(ValueError):
            arc_minimum_radius(arc, refine_samples=2)

    def test_arc_minimum_at_an_end_point(self) -> None:
        elements = OrbitalElements(1.6e7, 0.55, 0.4, 1.0, 2.0, 0.5)       # outbound: radius grows along the arc
        state = elements_to_state(elements)
        arc = propagate(state, 3000.0)
        self.assertAlmostEqual(arc_minimum_radius(arc), float(np.linalg.norm(state[:3])), delta=1e-6)

    def test_numerical_summary_of_a_plan(self) -> None:
        transfer = hohmann_transfer(LEO, GEO)
        state = np.array([LEO, 0.0, 0.0, 0.0, math.sqrt(EARTH_MU / LEO), 0.0])
        # Inward first burn: the vehicle dips to the periapsis before climbing to GEO.
        first = 2600.0 * np.array([math.sin(-0.3), math.cos(-0.3), 0.0])
        impulses = [(0.0, first), (transfer.transfer_time, [0.0, -100.0, 0.0])]
        summary = numerical_transfer_summary(state, impulses, thrust_n=250.0)
        self.assertAlmostEqual(summary.burn_delta_v_m_s[0], 2600.0, places=9)
        self.assertEqual((summary.burn_delta_v_m_s[1], summary.duration_s, summary.thrust_n),
                         (100.0, transfer.transfer_time, 250.0))
        after = state + np.concatenate(([0.0, 0.0, 0.0], first))
        speed2 = float(after[3:] @ after[3:])
        a = -EARTH_MU / (speed2 - 2 * EARTH_MU / LEO)
        momentum = LEO * after[4]
        periapsis = a * (1 - math.sqrt(1 - momentum**2 / (EARTH_MU * a)))
        self.assertLess(periapsis, LEO - 1e4)
        self.assertAlmostEqual(summary.minimum_radius_m, periapsis, delta=1.0)
        trajectory = propagate_impulsive_arcs(state, impulses, transfer.transfer_time)
        self.assertAlmostEqual(trajectory_minimum_radius(trajectory), summary.minimum_radius_m, places=6)
        with self.assertRaises(ValueError):
            numerical_transfer_summary(state, [])


if __name__ == "__main__":
    unittest.main()
