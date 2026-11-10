"""Coplanar two-impulse transfer problem (Phase 3 increment b): RV-3.01 and RV-3.04 on LEO->GEO.

The evaluator (state_to_elements + time of flight + propagate_kepler) is
compared with an independent closed-form oracle written here from vis-viva,
the flight-path geometry and Kepler's equation in eccentric anomaly; the
optimality reference is the Hohmann closed form of docs/03 section 8.2.

Error model for evaluator-vs-oracle tolerances, per point (eps = 2.2e-16):
  - Kepler propagation conditioning (tests/test_kepler.py): delta_a/a <=
    4 eps (3+e)/(1-e), giving an arrival timing error dt = 1.5 tau delta_a/a.
    The second burn changes at most at rate mu/r2^2 + v_c2 v/r2 (gravity plus
    rotation of the target direction), so |dJ| <= (mu/r2^2 + v_c2 v/r2) dt.
  - Arithmetic, including the oracle's cancellation in v_r^2 = v^2 - v_theta^2:
    64 eps (J + v^2/dv2).
  - Time of flight, in the well-conditioned region 1 + cos(nu2) >= 1e-2 only:
    64 eps tau kappa with kappa = sqrt((1+e)/(1-e))/sqrt(1 - cos(nu2)^2)
    (conditioning of nu -> E and of arccos), plus the period term 64 eps T.
Thresholds of the optimizer tests are those of the approved design: eps_low =
1e-6 m/s for "not better than Hohmann" and 1 m/s for practical convergence.
"""
import math
import unittest

import numpy as np

from orbital_optimizer import EARTH_MU, hohmann_transfer
from orbital_optimizer.optimization import (
    DifferentialEvolution, Evaluation, EvaluationStatus, TerminationCriteria, TerminationReason,
)
from orbital_optimizer.problems import (
    TANGENCY_FACTOR, CoplanarTwoImpulseProblem, PlanVerification, TwoImpulsePlan,
    hohmann_two_impulse_plan, reach_margin, verify_two_impulse_plan,
)

LEO, GEO = 6_678_000.0, 42_164_000.0
EPS = np.finfo(float).eps
EPS_LOW = 1e-6            # m/s: "not significantly better than Hohmann" (RV-3.01)
PRACTICAL_GAP = 1.0       # m/s: practical convergence threshold of the design


def oracle(r1: float, r2: float, mu: float, dv: float, alpha: float) -> dict:
    """Independent closed form for the first outward crossing of r2 (no propagation)."""
    vc1, vc2 = math.sqrt(mu / r1), math.sqrt(mu / r2)
    vt, vr = vc1 + dv * math.cos(alpha), dv * math.sin(alpha)
    energy = (vt * vt + vr * vr) / 2 - mu / r1
    h = r1 * vt
    a = -mu / (2 * energy)
    e = math.sqrt(max(0.0, 1 - h * h / (mu * a)))
    v2 = 2 * (energy + mu / r2)
    v_theta = h / r2
    v_radial = math.sqrt(max(0.0, v2 - v_theta * v_theta))
    dv2 = math.hypot(vc2 - v_theta, v_radial)
    n = math.sqrt(mu / a**3)
    # e sin E = (r . v)/sqrt(mu a), e cos E = 1 - r/a.
    e1 = math.atan2(r1 * vr / math.sqrt(mu * a), 1 - r1 / a)
    cos_e2 = max(-1.0, min(1.0, (1 - r2 / a) / e)) if e > 0 else 1.0
    e2 = math.acos(cos_e2)
    sweep = ((e2 - e * math.sin(e2)) - (e1 - e * math.sin(e1))) % (2 * math.pi)
    p = a * (1 - e * e)
    cos_nu2 = (p / r2 - 1) / e if e > 0 else 1.0
    return {"J": dv + dv2, "tau": sweep / n, "a": a, "e": e, "apoapsis": a * (1 + e), "dv2": dv2,
            "speed_r2": math.sqrt(max(v2, 0.0)), "cos_nu2": cos_nu2, "period": 2 * math.pi / n}


def j_tolerance(problem: CoplanarTwoImpulseProblem, ref: dict) -> float:
    mu, r2 = problem.mu, problem.final_radius
    dt = 1.5 * ref["tau"] * 4 * EPS * (3 + ref["e"]) / (1 - ref["e"])
    rate = mu / r2**2 + math.sqrt(mu / r2) * ref["speed_r2"] / r2
    return rate * dt + 64 * EPS * (ref["J"] + ref["speed_r2"] ** 2 / max(ref["dv2"], 1.0))


class ConstructionTests(unittest.TestCase):
    def test_bounds_and_reference(self) -> None:
        problem = CoplanarTwoImpulseProblem(LEO, GEO)
        reference = hohmann_transfer(LEO, GEO)
        self.assertEqual(problem.bounds.names, ("delta_v1", "alpha1"))
        self.assertEqual(problem.bounds.units, ("m/s", "rad"))
        self.assertEqual(problem.bounds.periodic.tolist(), [False, True])
        self.assertAlmostEqual(problem.bounds.upper[0], (math.sqrt(2) - 1) * math.sqrt(EARTH_MU / LEO), places=9)
        self.assertEqual((problem.bounds.lower[1], problem.bounds.upper[1]), (-np.pi, np.pi))
        self.assertEqual(problem.reference.total_delta_v, reference.total_delta_v)
        self.assertEqual(problem.hohmann_decision.tolist(), [reference.first_delta_v, 0.0])

    def test_invalid_construction(self) -> None:
        for r1, r2 in ((LEO, LEO), (GEO, LEO), (-1.0, GEO), (LEO, np.nan), (0.0, GEO)):
            with self.subTest(r1=r1, r2=r2), self.assertRaises(ValueError):
                CoplanarTwoImpulseProblem(r1, r2)
        with self.assertRaises(ValueError):
            CoplanarTwoImpulseProblem(LEO, GEO, mu=0.0)
        # sqrt(mu/r1) overflows: the Hohmann reference is not finite at this scale.
        with self.assertRaises(RuntimeError):
            CoplanarTwoImpulseProblem(1e-10, 2e-10, mu=1e308)


class EvaluationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.problem = CoplanarTwoImpulseProblem(LEO, GEO)
        self.reference = hohmann_transfer(LEO, GEO)

    def test_hohmann_decision_reproduces_hohmann(self) -> None:
        evaluation = self.problem.evaluate(self.problem.hohmann_decision)
        self.assertIs(evaluation.status, EvaluationStatus.VALID)
        self.assertLessEqual(abs(evaluation.objective - self.reference.total_delta_v), EPS_LOW)
        meta = evaluation.metadata
        self.assertAlmostEqual(meta["transfer_time_s"], self.reference.transfer_time, places=6)
        self.assertAlmostEqual(meta["delta_v1_m_s"], self.reference.first_delta_v, places=9)
        self.assertLessEqual(abs(meta["delta_v2_m_s"] - self.reference.second_delta_v), EPS_LOW)
        self.assertAlmostEqual(meta["apoapsis_radius_m"], GEO, delta=1e-3)
        self.assertAlmostEqual(meta["periapsis_radius_m"], LEO, delta=1e-3)
        self.assertEqual((self.problem.captured_value_errors, self.problem.captured_runtime_errors), (0, 0))

    def test_metadata_units_and_values(self) -> None:
        meta = self.problem.evaluate([2600.0, 0.3]).metadata
        for key in ("delta_v_total_m_s", "delta_v1_m_s", "delta_v2_m_s", "transfer_time_s", "semi_major_axis_m",
                    "eccentricity", "periapsis_radius_m", "apoapsis_radius_m", "arrival_radius_residual_m",
                    "departure_true_anomaly_rad", "arrival_true_anomaly_rad", "tangency_clip",
                    "reach_margin_relative", "reference_delta_v_m_s", "reference_transfer_time_s"):
            self.assertIn(key, meta)
        self.assertEqual(meta["first_burn_m_s"].shape, (3,))
        self.assertEqual(meta["second_burn_m_s"].shape, (3,))
        self.assertEqual(meta["tangency_clip"], 0.0)
        self.assertGreater(meta["apoapsis_radius_m"], GEO)
        self.assertLess(meta["periapsis_radius_m"], LEO)
        # Analytic f/g arrival: roundoff of order 1e-12 a (tests/test_kepler.py model), << 1 mm.
        self.assertLess(abs(meta["arrival_radius_residual_m"]), 1e-3)

    def test_evaluator_matches_closed_form_oracle_on_feasible_grid(self) -> None:
        upper = self.problem.bounds.upper[0]
        compared = timed = 0
        for dv in np.linspace(0.0, upper, 61)[1:-1]:
            for alpha in np.linspace(-np.pi, np.pi, 61)[:-1]:
                ref = oracle(LEO, GEO, EARTH_MU, dv, alpha)
                evaluation = self.problem.evaluate([dv, alpha])
                # Classification margin 1e-12: about 15x the admitted tangency band (6.3e-14).
                if ref["apoapsis"] < GEO * (1 - 1e-12):
                    self.assertIs(evaluation.status, EvaluationStatus.INVALID_DOMAIN, (dv, alpha))
                    continue
                if ref["apoapsis"] <= GEO * (1 + 1e-12):
                    continue
                self.assertIs(evaluation.status, EvaluationStatus.VALID, (dv, alpha))
                compared += 1
                with self.subTest(dv=dv, alpha=alpha):
                    self.assertLessEqual(abs(evaluation.objective - ref["J"]), j_tolerance(self.problem, ref))
                    if 1 + ref["cos_nu2"] >= 1e-2:
                        timed += 1
                        kappa = math.sqrt((1 + ref["e"]) / (1 - ref["e"])) / math.sqrt(1 - ref["cos_nu2"] ** 2)
                        tolerance = 64 * EPS * (ref["tau"] * kappa + ref["period"])
                        self.assertLessEqual(abs(evaluation.metadata["transfer_time_s"] - ref["tau"]), tolerance)
        self.assertGreater(compared, 100)
        self.assertGreater(timed, 50)
        self.assertEqual((self.problem.captured_value_errors, self.problem.captured_runtime_errors), (0, 0))

    def test_mirror_symmetry(self) -> None:
        # alpha and -alpha give conics mirrored about the departure radius: equal cost.
        for dv in np.linspace(2500.0, 3150.0, 8):
            for alpha in (0.05, 0.3, 0.7, 1.2):
                plus, minus = self.problem.evaluate([dv, alpha]), self.problem.evaluate([dv, -alpha])
                self.assertIs(plus.status, minus.status)
                if plus.is_valid:
                    with self.subTest(dv=dv, alpha=alpha):
                        tolerance = 2 * j_tolerance(self.problem, oracle(LEO, GEO, EARTH_MU, dv, alpha))
                        self.assertLessEqual(abs(plus.objective - minus.objective), tolerance)

    def test_no_feasible_sample_beats_hohmann(self) -> None:
        # RV-3.01 as a property: the Hohmann optimality theorem for this class is the oracle.
        rng = np.random.Generator(np.random.PCG64(2026))
        candidates = self.problem.bounds.from_unit(rng.random((40_000, 2)))
        feasible = [x for x in candidates
                    if oracle(LEO, GEO, EARTH_MU, x[0], x[1])["apoapsis"] > GEO * (1 + 1e-12)][:1000]
        self.assertEqual(len(feasible), 1000)
        objectives = []
        for x in feasible:
            evaluation = self.problem.evaluate(x)
            self.assertIs(evaluation.status, EvaluationStatus.VALID)
            objectives.append(evaluation.objective)
        self.assertGreaterEqual(min(objectives), self.reference.total_delta_v - EPS_LOW)

    def test_domain_predicates_reject_before_physics(self) -> None:
        upper = self.problem.bounds.upper[0]
        cases = {
            (0.0, 0.0): "does not reach",
            (1000.0, 0.0): "does not reach",
            (upper, 0.0): "unbound",
            (upper, np.pi - 1e-9): "does not reach",
        }
        for x, text in cases.items():
            with self.subTest(x=x):
                evaluation = self.problem.evaluate(x)
                self.assertIs(evaluation.status, EvaluationStatus.INVALID_DOMAIN)
                self.assertIn(text, evaluation.message)
        self.assertEqual((self.problem.captured_value_errors, self.problem.captured_runtime_errors), (0, 0))
        with self.assertRaises(ValueError):
            self.problem.evaluate([upper * 1.01, 0.0])     # outside the bounds: caller error
        with self.assertRaises(ValueError):
            self.problem.evaluate([2600.0, np.pi])         # periodic upper bound excluded

    def test_tangency_tolerance(self) -> None:
        u = EPS / 2
        vc1 = math.sqrt(EARTH_MU / LEO)
        q, s = reach_margin(LEO, GEO, EARTH_MU, vc1 + self.reference.first_delta_v, 0.0)
        self.assertGreaterEqual(q, -TANGENCY_FACTOR * u * s)
        # A deficit far beyond rounding (1e-9 relative in dv1) is rejected.
        evaluation = self.problem.evaluate([self.reference.first_delta_v * (1 - 1e-9), 0.0])
        self.assertIs(evaluation.status, EvaluationStatus.INVALID_DOMAIN)

    def test_plan_of_invalid_decision_raises(self) -> None:
        with self.assertRaises(ValueError):
            self.problem.plan([0.0, 0.0])


class VerificationTests(unittest.TestCase):
    def test_hohmann_plan_from_problem_passes(self) -> None:
        problem = CoplanarTwoImpulseProblem(LEO, GEO)
        plan = problem.plan(problem.hohmann_decision)
        self.assertAlmostEqual(plan.total_delta_v, hohmann_transfer(LEO, GEO).total_delta_v, delta=EPS_LOW)
        result = verify_two_impulse_plan(plan)
        self.assertTrue(result.passed, result.message)
        self.assertTrue(all(result.checks.values()))

    def test_descending_geo_to_leo_plan_passes(self) -> None:
        # Q8: descent is verification-only; this closes the gap of a propagated descending transfer.
        plan = hohmann_two_impulse_plan(GEO, LEO)
        self.assertLess(plan.first_burn[1], 0.0)        # retrograde first burn
        self.assertGreater(plan.second_burn[1], 0.0)    # at perigee, moving along -y, slow down
        result = verify_two_impulse_plan(plan)
        self.assertTrue(result.passed, result.message)
        self.assertLess(abs(result.metrics["final_semi_major_axis_error_m"]), 1.0)

    def test_missing_second_burn_is_detected(self) -> None:
        plan = hohmann_two_impulse_plan(LEO, GEO)
        incomplete = TwoImpulsePlan(plan.initial_state, plan.first_burn, [0.0, 0.0, 0.0], plan.transfer_time, GEO)
        result = verify_two_impulse_plan(incomplete)
        self.assertFalse(result.passed)
        self.assertFalse(result.checks["final_eccentricity"])
        self.assertGreater(result.metrics["final_eccentricity"], 0.7)

    def test_plan_and_verification_contracts(self) -> None:
        plan = hohmann_two_impulse_plan(LEO, GEO)
        with self.assertRaises(ValueError):
            plan.first_burn[0] = 1.0
        for kwargs in ({"transfer_time": 0.0}, {"final_radius": -1.0}, {"first_burn": [1.0, 2.0]}):
            arguments = {"initial_state": plan.initial_state, "first_burn": plan.first_burn,
                         "second_burn": plan.second_burn, "transfer_time": plan.transfer_time,
                         "final_radius": GEO} | kwargs
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                TwoImpulsePlan(**arguments)
        with self.assertRaises(ValueError):
            hohmann_two_impulse_plan(LEO, LEO)
        with self.assertRaises(ValueError):
            verify_two_impulse_plan("plan")  # type: ignore[arg-type]
        result = verify_two_impulse_plan(plan)
        self.assertIsInstance(result, PlanVerification)
        with self.assertRaises(TypeError):
            result.checks["final_radius"] = False  # type: ignore[index]
        self.assertEqual(result.tolerances["radius_m"], 1.0)


class DifferentialEvolutionEndToEndTests(unittest.TestCase):
    def test_leo_to_geo_small_budget(self) -> None:
        # Regression guard, not the 30-seed benchmark (increment c).
        reference = hohmann_transfer(LEO, GEO).total_delta_v
        for seed in (0, 1, 2):
            with self.subTest(seed=seed):
                problem = CoplanarTwoImpulseProblem(LEO, GEO)
                result = DifferentialEvolution().minimize(problem, TerminationCriteria(max_evaluations=2000), seed=seed)
                self.assertEqual(result.nfev, 2000)
                self.assertIs(result.termination_reason, TerminationReason.MAX_EVALUATIONS)
                self.assertTrue(result.found_valid)
                gap = result.f_best - reference
                self.assertGreaterEqual(gap, -EPS_LOW)
                self.assertLessEqual(gap, PRACTICAL_GAP)
                self.assertEqual((problem.captured_value_errors, problem.captured_runtime_errors), (0, 0))
                self.assertIsInstance(result.best_evaluation, Evaluation)
                verification = verify_two_impulse_plan(problem.plan(result.x_best))
                self.assertTrue(verification.passed, verification.message)


if __name__ == "__main__":
    unittest.main()
