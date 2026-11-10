"""RF-3.19 to RF-3.21: weighted, nondimensional objective over reported metrics."""
import math
import unittest

import numpy as np

from orbital_optimizer.optimization import (
    Bounds, DifferentialEvolution, Evaluation, EvaluationStatus, ObjectiveTerm, ScalarizedProblem,
    TerminationCriteria, WeightedObjective,
)


class MetricsProblem:
    """Two competing metrics on [0,1]: cost = 10 + 90 u (m/s), duration = 3600 (2 - u) (s).

    u < 0.1 is outside the domain. The cheapest point is u = 0.1, the fastest u = 1.
    """

    name = "metrics"
    bounds = Bounds([0.0], [1.0], names=("u",))

    def evaluate(self, x) -> Evaluation:
        u = float(x[0])
        if u < 0.1:
            return Evaluation(x, math.nan, EvaluationStatus.INVALID_DOMAIN, "outside domain")
        cost, duration = 10.0 + 90.0 * u, 3600.0 * (2.0 - u)
        return Evaluation(x, cost, EvaluationStatus.VALID, "", {"cost_m_s": cost, "duration_s": duration})


def term(name: str, metric: str, weight: float, reference: float, unit: str) -> ObjectiveTerm:
    return ObjectiveTerm(name, metric, weight, reference, unit, rationale="declared for the test")


class ObjectiveTermTests(unittest.TestCase):
    def test_terms_are_dimensionless_contributions(self) -> None:
        objective = WeightedObjective([term("cost", "cost_m_s", 0.75, 100.0, "m/s"),
                                       term("time", "duration_s", 0.25, 7200.0, "s")])
        contributions = objective.contributions({"cost_m_s": 55.0, "duration_s": 5400.0})
        self.assertEqual(contributions, {"cost": 0.75 * 55.0 / 100.0, "time": 0.25 * 5400.0 / 7200.0})
        self.assertAlmostEqual(objective.value({"cost_m_s": 55.0, "duration_s": 5400.0}), 0.6, places=15)
        # Changing the unit of a metric and of its reference together leaves J unchanged.
        in_km_s = WeightedObjective([term("cost", "cost_km_s", 0.75, 0.1, "km/s"),
                                     term("time", "duration_s", 0.25, 7200.0, "s")])
        self.assertAlmostEqual(in_km_s.value({"cost_km_s": 0.055, "duration_s": 5400.0}), 0.6, places=15)

    def test_description_records_every_parameter(self) -> None:
        objective = WeightedObjective([term("cost", "cost_m_s", 1.0, 100.0, "m/s")])
        description = objective.describe()
        self.assertEqual(description["terms"][0], {
            "name": "cost", "metric": "cost_m_s", "weight": 1.0, "reference": 100.0, "unit": "m/s",
            "rationale": "declared for the test"})
        self.assertIn("weight_i * metric_i / reference_i", description["form"])

    def test_zero_weight_terms_are_not_read(self) -> None:
        objective = WeightedObjective([term("cost", "cost_m_s", 1.0, 100.0, "m/s"),
                                       term("risk", "not_reported", 0.0, 1.0, "dimensionless")])
        self.assertEqual(objective.contributions({"cost_m_s": 50.0}), {"cost": 0.5, "risk": 0.0})

    def test_invalid_terms_and_metrics(self) -> None:
        valid = {"name": "cost", "metric": "cost_m_s", "weight": 1.0, "reference": 100.0, "unit": "m/s",
                 "rationale": "why"}
        for change in ({"name": ""}, {"metric": " "}, {"unit": ""}, {"rationale": ""}, {"rationale": None},
                       {"weight": -1.0}, {"weight": math.nan}, {"weight": True}, {"reference": 0.0},
                       {"reference": -5.0}, {"reference": math.inf}, {"reference": "100"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                ObjectiveTerm(**(valid | change))
        with self.assertRaises(ValueError):
            WeightedObjective([])
        with self.assertRaises(ValueError):
            WeightedObjective(["term"])  # type: ignore[list-item]
        with self.assertRaises(ValueError):       # duplicated names
            WeightedObjective([term("a", "m", 1.0, 1.0, "u"), term("a", "n", 1.0, 1.0, "u")])
        with self.assertRaises(ValueError):       # no positive weight
            WeightedObjective([term("a", "m", 0.0, 1.0, "u")])
        objective = WeightedObjective([term("cost", "cost_m_s", 1.0, 100.0, "m/s")])
        with self.assertRaises(ValueError):
            objective.contributions({})
        with self.assertRaises(ValueError):
            objective.contributions({"cost_m_s": np.array([1.0, 2.0])})
        with self.assertRaises(ValueError):
            objective.contributions({"cost_m_s": "55"})


class ScalarizedProblemTests(unittest.TestCase):
    def setUp(self) -> None:
        self.base = MetricsProblem()

    def objective(self, cost_weight: float, time_weight: float) -> WeightedObjective:
        # References: the values at the cheapest valid point u = 0.1.
        return WeightedObjective([term("cost", "cost_m_s", cost_weight, 19.0, "m/s"),
                                  term("time", "duration_s", time_weight, 6840.0, "s")])

    def test_evaluation_keeps_the_base_objective_and_the_terms(self) -> None:
        problem = ScalarizedProblem(self.base, self.objective(1.0, 1.0))
        self.assertIs(problem.bounds, self.base.bounds)
        self.assertEqual(problem.name, "metrics:weighted")
        evaluation = problem.evaluate([0.1])
        self.assertTrue(evaluation.is_valid)
        self.assertAlmostEqual(evaluation.objective, 2.0, places=15)      # both terms equal 1 at the reference
        self.assertEqual(evaluation.metadata["base_objective"], 19.0)
        self.assertEqual(set(evaluation.metadata["objective_terms"]), {"cost", "time"})
        self.assertEqual(evaluation.metadata["cost_m_s"], 19.0)

    def test_invalid_base_evaluations_pass_through(self) -> None:
        problem = ScalarizedProblem(self.base, self.objective(1.0, 1.0))
        evaluation = problem.evaluate([0.05])
        self.assertIs(evaluation.status, EvaluationStatus.INVALID_DOMAIN)
        self.assertEqual(evaluation.message, "outside domain")

    def test_non_finite_weighted_value(self) -> None:
        class Overflowing(MetricsProblem):
            def evaluate(self, x) -> Evaluation:
                return Evaluation(x, 1.0, EvaluationStatus.VALID, "", {"cost_m_s": math.inf, "duration_s": 1.0})

        evaluation = ScalarizedProblem(Overflowing(), self.objective(1.0, 1.0)).evaluate([0.5])
        self.assertIs(evaluation.status, EvaluationStatus.NON_FINITE)
        self.assertEqual(evaluation.objective, math.inf)

    def test_weights_move_the_optimum_along_the_trade_off(self) -> None:
        """Both metrics are linear in u, so J is linear and its minimum is an end point.

        J(u) = wc (10 + 90 u)/19 + wt 3600 (2 - u)/6840 has slope 90 wc/19 - 3600 wt/6840
        = 4.737 wc - 0.526 wt: u* = 0.1 (cheapest) when wc/wt > 1/9 and u* = 1 (fastest) otherwise.
        """
        for cost_weight, time_weight, expected in ((1.0, 1.0, 0.1), (1.0, 0.0, 0.1), (0.05, 1.0, 1.0), (0.0, 1.0, 1.0)):
            with self.subTest(cost_weight=cost_weight, time_weight=time_weight):
                problem = ScalarizedProblem(self.base, self.objective(cost_weight, time_weight))
                result = DifferentialEvolution().minimize(problem, TerminationCriteria(600), seed=0)
                self.assertTrue(result.found_valid)
                self.assertAlmostEqual(float(result.x_best[0]), expected, delta=1e-6)
                self.assertGreaterEqual(float(result.x_best[0]), 0.1)

    def test_constructor_validation(self) -> None:
        objective = self.objective(1.0, 1.0)
        with self.assertRaises(ValueError):
            ScalarizedProblem(object(), objective)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ScalarizedProblem(self.base, "objective")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            ScalarizedProblem(self.base, objective, name="")
        missing = WeightedObjective([term("risk", "not_reported", 1.0, 1.0, "dimensionless")])
        with self.assertRaises(ValueError):        # contract error of the problem: propagates
            ScalarizedProblem(self.base, missing).evaluate([0.5])


if __name__ == "__main__":
    unittest.main()
