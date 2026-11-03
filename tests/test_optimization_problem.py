"""Phase 3 problem contracts: bounds, normalization, evaluations and FunctionProblem."""
import math
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from orbital_optimizer.optimization import Bounds, Evaluation, EvaluationStatus, FunctionProblem, Problem


class BoundsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.bounds = Bounds([0.0, -np.pi], [3200.0, np.pi], names=("dv", "angle"),
                             units=("m/s", "rad"), periodic=np.array([False, True]))

    def test_valid_bounds_are_copied_read_only_with_metadata(self) -> None:
        lower, upper, periodic = np.array([0.0, -1.0]), np.array([1.0, 1.0]), np.array([False, True])
        bounds = Bounds(lower, upper, periodic=periodic)
        lower[0], upper[0], periodic[0] = 99.0, 99.0, True
        self.assertEqual((bounds.lower[0], bounds.upper[0], bool(bounds.periodic[0])), (0.0, 1.0, False))
        for array, source in ((bounds.lower, lower), (bounds.upper, upper), (bounds.periodic, periodic)):
            self.assertFalse(np.shares_memory(array, source))
        bounds = Bounds(np.array([0.0, -1.0]), [1.0, 1.0])
        self.assertEqual(bounds.dimension, 2)
        self.assertEqual(bounds.names, ("x0", "x1"))
        self.assertEqual(bounds.units, ("dimensionless", "dimensionless"))
        assert_array_equal(bounds.periodic, [False, False])
        for array in (bounds.lower, bounds.upper, bounds.periodic):
            with self.assertRaises(ValueError):
                array[0] = 0
        self.assertEqual(self.bounds.describe()["units"], ["m/s", "rad"])

    def test_invalid_bounds_are_rejected(self) -> None:
        cases = [
            ([], []),                         # empty
            ([0.0, 0.0], [1.0]),              # shape mismatch
            ([[0.0]], [[1.0]]),               # not one-dimensional
            ([0.0, 1.0], [1.0, 1.0]),         # lower == upper
            ([2.0], [1.0]),                   # lower > upper
            ([np.nan], [1.0]),                # non-finite
            ([0.0], [np.inf]),                # non-finite
            (["a"], ["b"]),                   # not numeric
            ([True], [True]),                 # boolean dtype
            ([-1e308], [1e308]),              # width overflows to inf
        ]
        for lower, upper in cases:
            with self.subTest(lower=lower, upper=upper), self.assertRaises(ValueError):
                Bounds(lower, upper)

    def test_invalid_labels_and_periodic_flags_are_rejected(self) -> None:
        for kwargs in (
            {"names": ("a",)}, {"names": ("a", "a")}, {"names": "ab"}, {"names": ("a", "")},
            {"units": ("m",)}, {"periodic": np.array([1, 0])}, {"periodic": np.array([True])},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                Bounds([0.0, 0.0], [1.0, 1.0], **kwargs)

    def test_unit_map_round_trip_and_endpoints(self) -> None:
        rng = np.random.Generator(np.random.PCG64(7))
        u = rng.random((100, 2))
        x = self.bounds.from_unit(u)
        assert_allclose(self.bounds.to_unit(x), u, rtol=0, atol=1e-15)
        assert_array_equal(self.bounds.from_unit([0.0, 0.0]), self.bounds.lower)
        # u=1 reaches upper for the non-periodic speed, but the periodic angle wraps to lower.
        assert_array_equal(self.bounds.from_unit([1.0, 1.0]), [3200.0, -np.pi])
        self.assertTrue(np.all(x >= self.bounds.lower) and np.all(x[:, 0] <= 3200.0) and np.all(x[:, 1] < np.pi))

    def test_periodic_variables_stay_in_half_open_period(self) -> None:
        # In [1,2) the float64 spacing is 2**-52, so 1 + (1 - 2**-53) rounds up to exactly 2.0.
        bounds = Bounds([1.0], [2.0], periodic=np.array([True]))
        self.assertEqual(1.0 + (1.0 - 2.0**-53), 2.0)
        self.assertEqual(bounds.from_unit([1.0 - 2.0**-53])[0], 1.0)
        self.assertEqual(bounds.from_unit([1.0])[0], 1.0)
        self.assertEqual(bounds.from_unit([1.0 - 2.0**-52])[0], 2.0 - 2.0**-52)
        near_one = 1.0 - np.arange(0, 65)[:, None] * 2.0**-53
        for lower, upper in ((1.0, 2.0), (-np.pi, np.pi), (0.0, 2 * np.pi), (-1e-300, 3.0), (1e6, 1e6 + 0.1)):
            with self.subTest(lower=lower, upper=upper):
                periodic = Bounds([lower], [upper], periodic=np.array([True]))
                x = periodic.from_unit(near_one)
                self.assertTrue(np.all(x >= lower) and np.all(x < upper))
                for value in x[:, 0]:
                    self.assertTrue(periodic.contains([value]))
                    periodic.validate([value])
                self.assertFalse(periodic.contains([upper]))
                with self.assertRaises(ValueError):
                    periodic.validate([upper])
        # A non-periodic variable keeps the closed interval.
        closed = Bounds([1.0], [2.0])
        self.assertEqual(closed.from_unit([1.0])[0], 2.0)
        self.assertTrue(closed.contains([2.0]))

    def test_unit_map_rejects_points_outside_the_hypercube_or_bad_shape(self) -> None:
        for u in ([-1e-12, 0.5], [0.5, 1.0 + 1e-12], [np.nan, 0.5], [0.5], [[0.5, 0.5, 0.5]]):
            with self.subTest(u=u), self.assertRaises(ValueError):
                self.bounds.from_unit(u)
        with self.assertRaises(ValueError):
            self.bounds.to_unit([1.0, 2.0, 3.0])

    def test_validate_returns_read_only_copy_and_rejects_bad_points(self) -> None:
        point = np.array([100.0, 0.0])
        validated = self.bounds.validate(point)
        point[0] = -5.0
        self.assertEqual(validated[0], 100.0)
        with self.assertRaises(ValueError):
            validated[0] = 1.0
        for x in ([-1.0, 0.0], [100.0, 4.0], [np.nan, 0.0], [1.0], [[1.0, 0.0]]):
            with self.subTest(x=x), self.assertRaises(ValueError):
                self.bounds.validate(x)
        self.assertTrue(self.bounds.contains([3200.0, -np.pi]))
        self.assertFalse(self.bounds.contains([3200.0, np.pi]))   # periodic: [lower, upper)
        self.assertFalse(self.bounds.contains([3200.1, 0.0]))


class EvaluationTests(unittest.TestCase):
    def test_valid_evaluation_requires_finite_objective(self) -> None:
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Evaluation(np.zeros(2), value, EvaluationStatus.VALID)
        evaluation = Evaluation(np.zeros(2), math.nan, EvaluationStatus.INVALID_DOMAIN, "outside")
        self.assertFalse(evaluation.is_valid)
        self.assertTrue(math.isnan(evaluation.objective))

    def test_evaluation_is_immutable(self) -> None:
        x = np.array([1.0, 2.0])
        evaluation = Evaluation(x, 3.0, EvaluationStatus.VALID, metadata={"tau_s": 1.0})
        x[0] = 9.0
        self.assertEqual(evaluation.x[0], 1.0)
        with self.assertRaises(ValueError):
            evaluation.x[0] = 0.0
        with self.assertRaises(TypeError):
            evaluation.metadata["tau_s"] = 2.0  # type: ignore[index]

    def test_invalid_fields_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            Evaluation(np.zeros(2), 1.0, "valid")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            Evaluation(np.zeros(2), [1.0, 2.0], EvaluationStatus.VALID)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            Evaluation(np.zeros((2, 2)), 1.0, EvaluationStatus.VALID)


class FunctionProblemTests(unittest.TestCase):
    def setUp(self) -> None:
        self.bounds = Bounds([-1.0, -1.0], [1.0, 1.0])

    def test_satisfies_problem_protocol_and_evaluates(self) -> None:
        problem = FunctionProblem(lambda x: float(np.sum(x**2)), self.bounds, "sq")
        self.assertIsInstance(problem, Problem)
        evaluation = problem.evaluate([0.5, -0.5])
        self.assertIs(evaluation.status, EvaluationStatus.VALID)
        self.assertEqual(evaluation.objective, 0.5)

    def test_non_finite_values_are_reported_with_raw_value(self) -> None:
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value):
                evaluation = FunctionProblem(lambda x, v=value: v, self.bounds).evaluate([0.0, 0.0])
                self.assertIs(evaluation.status, EvaluationStatus.NON_FINITE)
                self.assertFalse(evaluation.is_valid)
                self.assertEqual(repr(evaluation.objective), repr(value))

    def test_exceptions_propagate_unless_declared(self) -> None:
        def failing(x: np.ndarray) -> float:
            raise ArithmeticError("diverged")

        with self.assertRaises(ArithmeticError):
            FunctionProblem(failing, self.bounds).evaluate([0.0, 0.0])
        problem = FunctionProblem(failing, self.bounds, invalid_exceptions=(ArithmeticError,))
        evaluation = problem.evaluate([0.0, 0.0])
        self.assertIs(evaluation.status, EvaluationStatus.NUMERICAL_FAILURE)
        self.assertIn("diverged", evaluation.message)

    def test_function_receives_read_only_vector(self) -> None:
        def mutating(x: np.ndarray) -> float:
            x[0] = 0.0
            return 0.0

        with self.assertRaises(ValueError):
            FunctionProblem(mutating, self.bounds).evaluate([0.5, 0.5])

    def test_bad_outputs_inputs_and_configuration_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            FunctionProblem(lambda x: x, self.bounds).evaluate([0.0, 0.0])
        with self.assertRaises(ValueError):
            FunctionProblem(lambda x: "1", self.bounds).evaluate([0.0, 0.0])
        with self.assertRaises(ValueError):
            FunctionProblem(lambda x: 0.0, self.bounds).evaluate([2.0, 0.0])
        with self.assertRaises(ValueError):
            FunctionProblem(lambda x: 0.0, self.bounds, invalid_exceptions=(KeyboardInterrupt,))  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            FunctionProblem(lambda x: 0.0, self.bounds, invalid_exceptions=ValueError)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            FunctionProblem(0.0, self.bounds)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            FunctionProblem(lambda x: 0.0, ([-1.0], [1.0]))  # type: ignore[arg-type]


if __name__ == "__main__":
    unittest.main()
