"""RV-3.02: standard test functions with known optima, and DE regression guards on them.

The function checks compare against hand-computed values of the textbook
formulas. The DE guards use fixed seeds (so they are deterministic) and
asymmetric boxes so that the optimum is not the box centre. Their thresholds
were declared before the first run and are deliberately loose: they detect a
broken optimizer, they are not statistical performance evidence (that belongs
to a multi-seed benchmark report). Exception: the Rosenbrock guard was revised
after an independent verification measured its per-seed failure rate (see that
test's docstring); the change is recorded there, not hidden.
"""
import math
import unittest

import numpy as np

from orbital_optimizer.optimization import (
    DifferentialEvolution, TerminationCriteria, rastrigin, rastrigin_benchmark, rosenbrock,
    rosenbrock_benchmark, sphere, sphere_benchmark,
)

GUARD_SEEDS = (0, 1, 2, 3, 4)


class FunctionValueTests(unittest.TestCase):
    def test_known_optima(self) -> None:
        for dimension in (1, 2, 5, 10):
            self.assertEqual(sphere(np.zeros(dimension)), 0.0)
            self.assertEqual(rastrigin(np.zeros(dimension)), 0.0)
        for dimension in (2, 5, 10):
            self.assertEqual(rosenbrock(np.ones(dimension)), 0.0)

    def test_hand_computed_values(self) -> None:
        self.assertEqual(sphere([1.0, 2.0]), 5.0)
        self.assertEqual(rosenbrock([0.0, 0.0]), 1.0)
        # 100*(1 - 1.44)^2 + (1 + 1.2)^2 = 19.36 + 4.84
        self.assertAlmostEqual(rosenbrock([-1.2, 1.0]), 24.2, places=12)
        # 20 + 2*(1 - 10*cos(2*pi)) = 2
        self.assertAlmostEqual(rastrigin([1.0, 1.0]), 2.0, places=12)
        # x = 0.5: 10 + 0.25 - 10*cos(pi) = 20.25
        self.assertAlmostEqual(rastrigin([0.5]), 20.25, places=12)

    def test_invalid_inputs(self) -> None:
        for function, x in ((sphere, []), (sphere, [[1.0]]), (rastrigin, ["a"]), (rosenbrock, [1.0])):
            with self.subTest(function=function.__name__, x=x), self.assertRaises(ValueError):
                function(x)

    def test_benchmark_metadata(self) -> None:
        for benchmark in (sphere_benchmark(3), rosenbrock_benchmark(3), rastrigin_benchmark(3),
                          sphere_benchmark(2, -2.56, 5.12)):
            with self.subTest(name=benchmark.name):
                self.assertTrue(benchmark.bounds.contains(benchmark.minimizer))
                self.assertEqual(benchmark.function(benchmark.minimizer), benchmark.minimum)
                self.assertEqual(benchmark.problem().evaluate(benchmark.minimizer).objective, benchmark.minimum)
        self.assertEqual(rosenbrock_benchmark(2).bounds.lower.tolist(), [-5.0, -5.0])
        with self.assertRaises(ValueError):
            sphere_benchmark(2, lower=1.0, upper=2.0)   # excludes the minimizer
        with self.assertRaises(ValueError):
            rosenbrock_benchmark(1)
        with self.assertRaises(ValueError):
            sphere_benchmark(0)


class DifferentialEvolutionGuardTests(unittest.TestCase):
    def finals(self, benchmark, budget: int) -> list[float]:
        values = []
        for seed in GUARD_SEEDS:
            result = DifferentialEvolution().minimize(
                benchmark.problem(), TerminationCriteria(max_evaluations=budget), seed=seed)
            self.assertEqual(result.nfev, budget)
            self.assertTrue(result.found_valid)
            self.assertGreaterEqual(result.f_best, benchmark.minimum)
            values.append(result.f_best - benchmark.minimum)
        return values

    def test_sphere_2d_asymmetric(self) -> None:
        gaps = self.finals(sphere_benchmark(2, -2.56, 5.12), 2_000)
        self.assertLessEqual(max(gaps), 1e-8, gaps)

    def test_sphere_5d_asymmetric(self) -> None:
        gaps = self.finals(sphere_benchmark(5, -2.56, 5.12), 5_000)
        self.assertLessEqual(max(gaps), 1e-3, gaps)

    def test_rosenbrock_2d_asymmetric(self) -> None:
        """At least 4 of 5 seeds reach gap <= 1e-4; every gap is finite.

        Per-seed failure rate at this budget and box: an independent
        verification measured about 9.5 % of 200 seeds with gap > 1e-4. It
        was reproduced on 2026-11-03 (seeds 0-199: 19/200 failures, largest
        gap 0.45; Python 3.12.1, NumPy 2.3.3). Requiring all 5 seeds would make
        the guard depend on the seed choice; seeds 0-4 happen to pass.
        """
        gaps = self.finals(rosenbrock_benchmark(2, -2.048, 4.096), 4_000)
        self.assertGreaterEqual(sum(gap <= 1e-4 for gap in gaps), 4, gaps)
        self.assertTrue(all(math.isfinite(gap) for gap in gaps), gaps)

    def test_rastrigin_2d_asymmetric(self) -> None:
        # Multimodal: allow one seed in five to end in a non-global basin.
        gaps = self.finals(rastrigin_benchmark(2, -2.56, 5.12), 4_000)
        self.assertGreaterEqual(sum(gap <= 1e-4 for gap in gaps), 4, gaps)
        self.assertLessEqual(max(gaps), 2.0, gaps)
        self.assertTrue(all(math.isfinite(gap) for gap in gaps))


if __name__ == "__main__":
    unittest.main()
