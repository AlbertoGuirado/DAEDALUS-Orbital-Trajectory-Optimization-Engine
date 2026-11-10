"""Phase 4: constraint handling in the generic optimizers (RF-4.07 to RF-4.10).

Reference problem with a known constrained optimum, independent of orbits:

    minimize f(x) = x0 + x1 on [0, 1]^2 subject to x0 x1 >= 1/4.

By the arithmetic-geometric mean inequality x0 + x1 >= 2 sqrt(x0 x1) >= 1 on
the feasible set, with equality only at x0 = x1 = 1/2: f* = 1. Any feasible
point with f < 1 would be an arithmetic contradiction, so "f_best >= 1" is
exact, with no tolerance. The unconstrained minimum is f = 0 at the origin,
which is infeasible: an optimizer that ignores the constraint ends there.
"""
import json
import math
import unittest

import numpy as np

from orbital_optimizer.optimization import (
    Bounds, DeathPenalty, DifferentialEvolution, Evaluation, EvaluationBudget, EvaluationStatus,
    FeasibilityRules, FunctionProblem, GeneticAlgorithm, HybridOptimizer, NelderMead, NelderMeadSettings,
    ObjectiveTerm, ParticleSwarm, RandomSearch, ScalarizedProblem, StaticPenalty, TerminationCriteria,
    WeightedObjective, feasibility_key, handling_from_mapping, sphere,
)
from orbital_optimizer.optimization.differential_evolution import _trial_wins

SCALE = 0.25


class ProductProblem:
    """f = x0 + x1 with the hard constraint 1/4 - x0 x1 <= 0 (normalized by 1/4); x0 < 0.02 is invalid."""

    name = "product"
    bounds = Bounds([0.0, 0.0], [1.0, 1.0])

    def __init__(self, handling=None) -> None:
        if handling is not None:
            self.constraint_handling = handling

    def evaluate(self, x) -> Evaluation:
        if x[0] < 0.02:
            return Evaluation(x, math.nan, EvaluationStatus.INVALID_DOMAIN, "outside the domain")
        violation = max(0.0, SCALE - float(x[0] * x[1])) / SCALE
        return Evaluation(x, float(x[0] + x[1]), EvaluationStatus.VALID, "", {"cost": float(x[0] + x[1])},
                          violation, {"product": violation})


def evaluation(objective: float, violation: float | None, *, valid: bool = True, parts=None) -> Evaluation:
    status = EvaluationStatus.VALID if valid else EvaluationStatus.INVALID_DOMAIN
    constraints = {} if violation is None else (parts if parts is not None else {"c": violation})
    return Evaluation([0.5], objective if valid else math.nan, status, "", {}, violation if valid else None,
                      constraints if valid else {})


class EvaluationFieldTests(unittest.TestCase):
    def test_feasibility_flags(self) -> None:
        unconstrained = evaluation(3.0, None)
        self.assertEqual((unconstrained.is_constrained, unconstrained.is_feasible), (False, True))
        feasible = evaluation(3.0, 0.0)
        self.assertEqual((feasible.is_constrained, feasible.is_feasible), (True, True))
        infeasible = evaluation(3.0, 0.4, parts={"a": 0.1, "b": 0.3})
        self.assertEqual((infeasible.is_constrained, infeasible.is_feasible), (True, False))
        self.assertEqual(dict(infeasible.constraints), {"a": 0.1, "b": 0.3})
        invalid = evaluation(3.0, None, valid=False)
        self.assertFalse(invalid.is_feasible)
        self.assertFalse(evaluation(1.0, math.inf).is_feasible)          # a constraint that could not be evaluated
        with self.assertRaises(TypeError):
            infeasible.constraints["a"] = 0.0  # type: ignore[index]

    def test_invalid_constraint_fields(self) -> None:
        for kwargs in ({"constraint_violation": -0.1}, {"constraint_violation": math.nan},
                       {"constraint_violation": True}, {"constraint_violation": "0"},
                       {"constraint_violation": 0.1, "constraints": {"a": -1.0}},
                       {"constraint_violation": 0.1, "constraints": {"a": math.nan}},
                       {"constraint_violation": 0.1, "constraints": {1: 0.1}},
                       {"constraints": {"a": 0.1}}):                      # per-constraint values without a total
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                Evaluation([0.5], 1.0, EvaluationStatus.VALID, "", {}, **kwargs)


class StrategyTests(unittest.TestCase):
    def test_feasibility_rules_are_deb_three_rules(self) -> None:
        rules = FeasibilityRules()
        cheap_infeasible, costly_feasible = evaluation(1.0, 0.2), evaluation(50.0, 0.0)
        self.assertLess(rules.key(costly_feasible), rules.key(cheap_infeasible))         # feasible beats infeasible
        self.assertLess(rules.key(evaluation(2.0, 0.0)), rules.key(evaluation(3.0, 0.0)))  # both feasible: objective
        self.assertLess(rules.key(evaluation(9.0, 0.1)), rules.key(evaluation(1.0, 0.2)))  # both infeasible: violation
        self.assertLess(rules.key(cheap_infeasible), rules.key(evaluation(0.0, None, valid=False)))
        self.assertEqual(rules.key(evaluation(4.0, None)), (0, 4.0))                     # unconstrained: by objective
        self.assertEqual(rules.key(evaluation(1.0, math.inf)), (1, math.inf))
        self.assertEqual(feasibility_key(costly_feasible), rules.key(costly_feasible))

    def test_death_penalty_ties_everything_infeasible(self) -> None:
        death = DeathPenalty()
        self.assertEqual(death.key(evaluation(2.0, 0.0)), (0, 2.0))
        tied = [death.key(evaluation(1.0, 0.01)), death.key(evaluation(9.0, 5.0)), death.key(evaluation(0.0, None, valid=False))]
        self.assertEqual(tied, [(1, 0.0)] * 3)

    def test_static_penalty_value(self) -> None:
        penalty = StaticPenalty(100.0, "test coefficient")
        both = evaluation(10.0, 0.5, parts={"a": 0.1, "b": 0.4})
        self.assertAlmostEqual(penalty.penalized(both), 10.0 + 100.0 * (0.1**2 + 0.4**2), places=12)
        self.assertEqual(penalty.key(both), (0, penalty.penalized(both)))
        self.assertEqual(penalty.penalized(evaluation(10.0, 0.0)), 10.0)
        self.assertEqual(penalty.penalized(evaluation(10.0, None)), 10.0)
        total_only = Evaluation([0.5], 10.0, EvaluationStatus.VALID, "", {}, 0.3)
        self.assertAlmostEqual(penalty.penalized(total_only), 10.0 + 100.0 * 0.09, places=12)
        self.assertEqual(penalty.key(evaluation(0.0, None, valid=False)), (1, 0.0))
        # A cheap, slightly infeasible point can outrank a feasible one: the penalty is inexact.
        self.assertLess(penalty.key(evaluation(1.0, 0.01)), penalty.key(evaluation(2.0, 0.0)))

    def test_configuration_and_validation(self) -> None:
        self.assertIsInstance(handling_from_mapping({"type": "feasibility_rules"}), FeasibilityRules)
        self.assertIsInstance(handling_from_mapping({"type": "death_penalty"}), DeathPenalty)
        penalty = handling_from_mapping({"type": "static_penalty", "coefficient": 7.0, "rationale": "why"})
        self.assertEqual((penalty.coefficient, penalty.describe()["rationale"], penalty.describe()["name"]),
                         (7.0, "why", "static-penalty"))
        for bad in ({"type": "barrier"}, {"type": "static_penalty"}, {"type": "static_penalty", "coefficient": 1.0},
                    {"type": "feasibility_rules", "coefficient": 1.0}, {"coefficient": 1.0}, "feasibility_rules"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                handling_from_mapping(bad)
        for coefficient, rationale in ((0.0, "x"), (-1.0, "x"), (math.inf, "x"), (True, "x"), (1.0, ""), (1.0, None)):
            with self.subTest(coefficient=coefficient, rationale=rationale), self.assertRaises(ValueError):
                StaticPenalty(coefficient, rationale)  # type: ignore[arg-type]
        for strategy in (FeasibilityRules(), DeathPenalty(), penalty):
            json.dumps(strategy.describe())

    def test_key_selection_equals_the_unconstrained_rule_of_differential_evolution(self) -> None:
        rules = FeasibilityRules()
        cases = [(True, 1.0, True, 1.0), (True, 0.5, True, 1.0), (True, 2.0, True, 1.0),
                 (True, 1e300, False, math.inf), (False, math.inf, True, 1e300), (False, math.inf, False, math.inf)]
        for trial_valid, trial_value, parent_valid, parent_value in cases:
            with self.subTest(case=(trial_valid, trial_value, parent_valid, parent_value)):
                trial = evaluation(trial_value, None, valid=trial_valid)
                parent = evaluation(parent_value, None, valid=parent_valid)
                self.assertEqual(rules.key(trial) <= rules.key(parent),
                                 _trial_wins(trial_valid, trial_value, parent_valid, parent_value))


class BudgetTests(unittest.TestCase):
    def test_best_is_the_best_feasible_point(self) -> None:
        budget = EvaluationBudget(ProductProblem(), 10)
        points = [[0.01, 0.5],   # invalid
                  [0.1, 0.1],    # valid, f = 0.2, x0 x1 = 0.01: violation 0.96
                  [0.4, 0.5],    # valid, f = 0.9, x0 x1 = 0.20: violation 0.20
                  [0.9, 0.9],    # feasible, f = 1.8
                  [0.3, 0.3],    # valid, f = 0.6, violation 0.64
                  [0.6, 0.5]]    # feasible, f = 1.1
        for point in points:
            budget.evaluate(np.array(point))
        self.assertTrue(budget.constrained)
        self.assertEqual((budget.best_value, budget.best_x.tolist()), (1.1, [0.6, 0.5]))
        self.assertTrue(budget.best_evaluation.is_feasible)
        self.assertAlmostEqual(budget.least_violation, 0.2, places=12)
        self.assertEqual(budget.least_violation_x.tolist(), [0.4, 0.5])
        self.assertEqual(budget.least_violation_evaluation.objective, 0.9)
        self.assertEqual(budget.valid_count, 5)
        # The history counts feasible points only: +inf until the first one.
        self.assertEqual(budget.best_so_far().tolist(), [math.inf, math.inf, math.inf, 1.8, 1.8, 1.1])
        violations = budget.violations()
        self.assertTrue(math.isnan(violations[0]))
        np.testing.assert_allclose(violations[1:], [0.96, 0.2, 0.0, 0.64, 0.0], atol=1e-12)
        record = budget.constraint_record()
        self.assertEqual((record["found_feasible"], record["constraint_handling"]["name"]), (True, "feasibility-rules"))
        self.assertAlmostEqual(record["least_violation"], 0.2, places=12)

    def test_unconstrained_problem_keeps_the_previous_behaviour(self) -> None:
        budget = EvaluationBudget(FunctionProblem(sphere, Bounds([-1.0, -1.0], [1.0, 1.0])), 5)
        for point in ([0.9, 0.9], [0.5, 0.6], [0.1, 0.2]):
            budget.evaluate(np.array(point))
        self.assertFalse(budget.constrained)
        self.assertEqual(budget.constraint_record(), {})
        self.assertIsNone(budget.least_violation_evaluation)
        self.assertEqual(budget.violations().tolist(), [0.0, 0.0, 0.0])
        self.assertAlmostEqual(budget.best_value, 0.04, places=12)       # (0.5 -> 0, 0.6 -> 0.2): 0^2 + 0.2^2

    def test_handling_resolution_order(self) -> None:
        penalty, death = StaticPenalty(5.0, "test"), DeathPenalty()
        self.assertIsInstance(EvaluationBudget(ProductProblem(), 5).handling, FeasibilityRules)
        self.assertIs(EvaluationBudget(ProductProblem(penalty), 5).handling, penalty)
        self.assertIs(EvaluationBudget(ProductProblem(penalty), 5, handling=death).handling, death)
        with self.assertRaises(TypeError):
            EvaluationBudget(ProductProblem(), 5, handling="deb")  # type: ignore[arg-type]
        budget = EvaluationBudget(ProductProblem(penalty), 5)
        infeasible = budget.evaluate(np.array([0.4, 0.5]))
        self.assertEqual(budget.key(infeasible), penalty.key(infeasible))


class OptimizerTests(unittest.TestCase):
    """Every optimizer, with the default feasibility rules, on the product problem."""

    def test_feasible_optimum_is_found_and_never_undercut(self) -> None:
        optimizers = [DifferentialEvolution(), ParticleSwarm(), GeneticAlgorithm(), RandomSearch(), HybridOptimizer(),
                      NelderMead(NelderMeadSettings(initial_point=(0.8, 0.7)))]
        for optimizer in optimizers:
            with self.subTest(optimizer=optimizer.name):
                result = optimizer.minimize(ProductProblem(), TerminationCriteria(1500), seed=3)
                self.assertTrue(result.is_constrained)
                self.assertTrue(result.found_feasible)
                self.assertGreaterEqual(result.f_best, 1.0)                 # exact: AM-GM on the feasible set
                self.assertGreaterEqual(float(result.x_best[0] * result.x_best[1]), SCALE)
                self.assertEqual(result.best_evaluation.constraint_violation, 0.0)
                self.assertEqual(result.constraint_handling["name"], "feasibility-rules")
                self.assertEqual(result.f_best, result.history.best_so_far[-1])
                self.assertEqual(len(result.history.violation), 1500)
                self.assertIn("feasible", result.message)
        # The population-based methods converge to the constrained optimum f* = 1 (declared loose bound).
        for optimizer in (DifferentialEvolution(), ParticleSwarm(), HybridOptimizer()):
            with self.subTest(optimizer=optimizer.name, check="convergence"):
                result = optimizer.minimize(ProductProblem(), TerminationCriteria(1500), seed=3)
                self.assertLessEqual(result.f_best, 1.01)

    def test_infeasible_points_are_cheaper_but_never_reported_as_the_solution(self) -> None:
        result = DifferentialEvolution().minimize(ProductProblem(), TerminationCriteria(1500), seed=3)
        # The search evaluated infeasible points below f* = 1; they are kept apart as diagnostics.
        self.assertIsNotNone(result.least_violation_evaluation)
        self.assertGreater(result.least_violation, 0.0)
        self.assertFalse(result.least_violation_evaluation.is_feasible)
        record = result.to_dict()
        self.assertEqual(record["constraints"]["found_feasible"], True)
        self.assertEqual(record["best_evaluation"]["constraint_violation"], 0.0)
        self.assertEqual(record["constraints"]["least_violation_evaluation"]["constraint_violation"],
                         result.least_violation)
        self.assertEqual(len(record["history"]["violation"]), 1500)
        json.dumps(record, allow_nan=False)

    def test_static_penalty_searches_infeasible_points_but_the_result_stays_feasible(self) -> None:
        """Quadratic penalty with rho = 1: on the diagonal x0 = x1 = s the penalized function is

        2 s + ((1/4 - s^2) / (1/4))^2, whose derivative at s = 1/2 is 2 > 0, so its minimizer lies at
        s < 1/2, inside the infeasible region. The reported best must still be a feasible point.
        """
        penalty = StaticPenalty(1.0, "deliberately weak, to expose the inexactness of the penalty")

        class Recording(ProductProblem):
            def __init__(self) -> None:
                super().__init__(penalty)
                self.seen: list[Evaluation] = []

            def evaluate(self, x) -> Evaluation:
                self.seen.append(super().evaluate(x))
                return self.seen[-1]

        problem = Recording()
        result = DifferentialEvolution().minimize(problem, TerminationCriteria(3000), seed=3)
        self.assertEqual(result.constraint_handling["name"], "static-penalty")
        self.assertEqual(result.constraint_handling["coefficient"], 1.0)
        # Analytic minimum of the penalized function, on the diagonal by symmetry:
        # phi(s) = 2 s + 16 (1/4 - s^2)^2, phi'(s) = 2 - 64 s (1/4 - s^2) = 0.
        low, high = 0.40, 0.45
        for _ in range(60):
            middle = (low + high) / 2
            low, high = (middle, high) if 2 - 64 * middle * (0.25 - middle**2) < 0 else (low, middle)
        expected = 2 * low + 16 * (0.25 - low**2) ** 2
        self.assertAlmostEqual(expected, 0.927, places=3)
        self.assertLess(expected, 1.0)                                   # cheaper than every feasible point
        searched = min((e for e in problem.seen if e.is_valid), key=penalty.penalized)
        self.assertGreaterEqual(penalty.penalized(searched), expected - 1e-12)
        self.assertLessEqual(penalty.penalized(searched), expected + 1e-2)
        self.assertFalse(searched.is_feasible)                           # the search converged to an infeasible point
        self.assertAlmostEqual(float(searched.x[0]), low, delta=0.05)
        if result.found_feasible:
            self.assertGreaterEqual(result.f_best, 1.0)
            self.assertEqual(result.best_evaluation.constraint_violation, 0.0)
        else:
            self.assertEqual((result.f_best, result.x_best), (math.inf, None))

    def test_no_feasible_point_is_reported_as_such(self) -> None:
        class Impossible(ProductProblem):
            def evaluate(self, x) -> Evaluation:
                base = super().evaluate(x)
                if not base.is_valid:
                    return base
                violation = 1.0 + float(x[0])                            # never zero; smallest at x0 -> 0.02
                return Evaluation(x, base.objective, EvaluationStatus.VALID, "", {}, violation, {"never": violation})

        for optimizer in (DifferentialEvolution(), ParticleSwarm(), GeneticAlgorithm(), RandomSearch(), HybridOptimizer()):
            with self.subTest(optimizer=optimizer.name):
                result = optimizer.minimize(Impossible(), TerminationCriteria(400), seed=1)
                self.assertTrue(result.found_valid)
                self.assertFalse(result.found_feasible)
                self.assertEqual((result.x_best, result.u_best, result.best_evaluation), (None, None, None))
                self.assertEqual(result.f_best, math.inf)
                self.assertTrue(np.all(np.isinf(result.history.best_so_far)))
                self.assertGreaterEqual(result.least_violation, 1.02)
                self.assertIn("no feasible point found", result.message)
                json.dumps(result.to_dict(), allow_nan=False)
        # The feasibility rules descend on the violation: DE approaches its infimum 1.02.
        result = DifferentialEvolution().minimize(Impossible(), TerminationCriteria(1500), seed=1)
        self.assertLess(result.least_violation, 1.03)

    def test_unconstrained_record_has_no_constraint_section(self) -> None:
        result = DifferentialEvolution().minimize(FunctionProblem(sphere, Bounds([-1.0, -1.0], [1.0, 1.0])),
                                                  TerminationCriteria(200), seed=0)
        self.assertFalse(result.is_constrained)
        self.assertEqual((result.found_feasible, result.constraint_handling, result.least_violation), (None, None, None))
        record = result.to_dict()
        self.assertNotIn("constraints", record)
        self.assertNotIn("violation", record["history"])
        self.assertNotIn("constraint_violation", record["best_evaluation"])
        self.assertIsNone(result.history.violation)

    def test_scalarized_problem_keeps_the_constraints(self) -> None:
        penalty = StaticPenalty(3.0, "test")
        base = ProductProblem(penalty)
        objective = WeightedObjective([ObjectiveTerm("cost", "cost", 2.0, 1.0, "dimensionless", "doubles the cost")])
        scalarized = ScalarizedProblem(base, objective)
        self.assertIs(scalarized.constraint_handling, penalty)
        evaluated = scalarized.evaluate(np.array([0.4, 0.5]))
        self.assertAlmostEqual(evaluated.objective, 1.8, places=12)
        self.assertAlmostEqual(evaluated.constraint_violation, 0.2, places=12)
        self.assertEqual(set(evaluated.constraints), {"product"})
        result = DifferentialEvolution().minimize(ScalarizedProblem(ProductProblem(), objective),
                                                  TerminationCriteria(1500), seed=3)
        self.assertTrue(result.found_feasible)
        self.assertGreaterEqual(result.f_best, 2.0)                       # twice the constrained optimum


if __name__ == "__main__":
    unittest.main()
