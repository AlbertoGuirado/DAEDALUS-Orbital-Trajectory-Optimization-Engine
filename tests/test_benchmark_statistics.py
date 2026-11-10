"""Benchmark statistics against hand-computed values, brute force and exact enumeration.

No test compares a statistic with a value produced earlier by the same code:
quantiles and Holm adjustments are worked by hand, A12 and U are recounted
pair by pair, the Mann-Whitney p-value is enumerated over all rank
arrangements and the Wilson interval is checked through its defining
quadratic equation.
"""
import itertools
import math
import unittest

import numpy as np
from scipy import stats

from orbital_optimizer.benchmarking import (
    convergence_band, evaluations_to_target, friedman, holm_correction, mann_whitney, quantile,
    summarize, vargha_delaney_a12, wilcoxon_signed_rank, wilson_interval,
)

INF = math.inf


class OrderStatisticsTests(unittest.TestCase):
    def test_inverted_cdf_quantiles(self) -> None:
        values = [7.0, 1.0, 5.0, 3.0]                       # sorted: 1, 3, 5, 7
        self.assertEqual([quantile(values, q) for q in (0.0, 0.25, 0.5, 0.75, 1.0)], [1.0, 1.0, 3.0, 5.0, 7.0])
        self.assertEqual(quantile(values, 0.26), 3.0)       # ceil(0.26 * 4) = 2
        self.assertEqual(quantile([2.0], 0.5), 2.0)
        self.assertEqual(quantile([1.0, INF, 2.0], 0.5), 2.0)
        self.assertEqual(quantile([1.0, INF, INF], 0.5), INF)

    def test_summary(self) -> None:
        summary = summarize([4.0, 1.0, 3.0, 2.0, 10.0])
        self.assertEqual((summary["n"], summary["n_finite"]), (5, 5))
        self.assertEqual((summary["best"], summary["q1"], summary["median"], summary["q3"], summary["worst"]),
                         (1.0, 2.0, 3.0, 4.0, 10.0))
        self.assertEqual(summary["iqr"], 2.0)
        self.assertAlmostEqual(summary["mean"], 4.0, places=15)
        self.assertAlmostEqual(summary["std"], math.sqrt(50.0 / 4), places=14)   # sum of squares 50, ddof 1

    def test_summary_with_failed_runs(self) -> None:
        summary = summarize([1.0, INF, 3.0, INF])
        self.assertEqual((summary["n"], summary["n_finite"], summary["worst"]), (4, 2, INF))
        self.assertEqual((summary["median"], summary["q3"]), (3.0, INF))
        self.assertTrue(math.isnan(summary["iqr"]))
        self.assertEqual(summary["mean"], 2.0)                # finite values only
        all_failed = summarize([INF, INF])
        self.assertTrue(math.isnan(all_failed["mean"]) and math.isnan(all_failed["std"]))
        self.assertTrue(math.isnan(summarize([5.0])["std"]))

    def test_invalid_samples(self) -> None:
        for values in ([], [[1.0, 2.0]], [math.nan], [-INF], ["a"]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                summarize(values)
        for q in (-0.1, 1.1, True, "0.5"):
            with self.subTest(q=q), self.assertRaises(ValueError):
                quantile([1.0, 2.0], q)  # type: ignore[arg-type]


class WilsonIntervalTests(unittest.TestCase):
    def test_bounds_solve_the_score_equation(self) -> None:
        """The Wilson limits are the roots p of (p_hat - p)^2 = z^2 p (1 - p) / n."""
        z = stats.norm.ppf(0.975)
        for successes, trials in ((0, 30), (1, 30), (15, 30), (27, 30), (30, 30), (3, 7)):
            low, high = wilson_interval(successes, trials)
            observed = successes / trials
            with self.subTest(successes=successes, trials=trials):
                self.assertTrue(0.0 <= low <= observed <= high <= 1.0)
                for bound in (low, high):
                    if 0.0 < bound < 1.0:
                        self.assertAlmostEqual((observed - bound) ** 2, z * z * bound * (1 - bound) / trials, places=12)
        self.assertEqual(wilson_interval(0, 30)[0], 0.0)
        self.assertEqual(wilson_interval(30, 30)[1], 1.0)
        # Rule of three: with 0/n successes the 95 % upper limit is about 3.84/(n + 3.84).
        self.assertAlmostEqual(wilson_interval(0, 30)[1], z * z / (30 + z * z), places=12)

    def test_confidence_and_validation(self) -> None:
        narrow, wide = wilson_interval(10, 20, 0.80), wilson_interval(10, 20, 0.99)
        self.assertLess(wide[0], narrow[0])
        self.assertGreater(wide[1], narrow[1])
        for arguments in ((-1, 10), (11, 10), (1, 0), (1.5, 10), (True, 10)):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                wilson_interval(*arguments)  # type: ignore[arg-type]
        for confidence in (0.0, 1.0, 1.5):
            with self.subTest(confidence=confidence), self.assertRaises(ValueError):
                wilson_interval(1, 10, confidence)


class ConvergenceTests(unittest.TestCase):
    def test_evaluations_to_target(self) -> None:
        history = [INF, 9.0, 9.0, 4.0, 4.0, 1.0]
        self.assertEqual(evaluations_to_target(history, 9.0), 2)
        self.assertEqual(evaluations_to_target(history, 4.5), 4)
        self.assertEqual(evaluations_to_target(history, 1.0), 6)
        self.assertIsNone(evaluations_to_target(history, 0.5))
        self.assertEqual(evaluations_to_target([3.0], 3.0), 1)

    def test_band_reports_observed_values_per_evaluation(self) -> None:
        histories = [[INF, 5.0, 2.0, 2.0], [9.0, 9.0, 3.0, 1.0], [8.0, 4.0, 4.0, 4.0]]
        band = convergence_band(histories)
        self.assertEqual(band["evaluations"], [1, 2, 3, 4])
        self.assertEqual(band["quantiles"][0.5], [9.0, 5.0, 3.0, 2.0])
        self.assertEqual(band["quantiles"][0.25], [8.0, 4.0, 2.0, 1.0])
        self.assertEqual(band["quantiles"][0.75], [INF, 9.0, 4.0, 4.0])
        subset = convergence_band(histories, np.array([1, 4]), quantiles=(0.5,))
        self.assertEqual(subset, {"evaluations": [1, 4], "quantiles": {0.5: [9.0, 2.0]}})

    def test_band_validation(self) -> None:
        with self.assertRaises(ValueError):
            convergence_band([])
        with self.assertRaises(ValueError):
            convergence_band([[1.0, 2.0], [1.0]])
        for evaluations in (np.array([0]), np.array([3]), np.array([1.5]), np.array([[1]])):
            with self.subTest(evaluations=evaluations), self.assertRaises(ValueError):
                convergence_band([[1.0, 2.0]], evaluations)
        with self.assertRaises(ValueError):
            convergence_band([[1.0, 2.0]], quantiles=(1.5,))


def exact_two_sided_p(first: list[float], second: list[float]) -> float:
    """Exact Mann-Whitney p-value by enumerating every assignment of the pooled values (no ties)."""
    pooled = first + second
    size = len(first)

    def statistic(sample: tuple[float, ...], rest: list[float]) -> int:
        return sum(a > b for a in sample for b in rest)

    centre = len(first) * len(second) / 2
    observed = abs(statistic(tuple(first), second) - centre)
    extreme = total = 0
    for indices in itertools.combinations(range(len(pooled)), size):
        sample = tuple(pooled[i] for i in indices)
        rest = [pooled[i] for i in range(len(pooled)) if i not in indices]
        total += 1
        extreme += abs(statistic(sample, rest) - centre) >= observed - 1e-12
    return extreme / total


class ComparisonTests(unittest.TestCase):
    def test_a12_by_counting_pairs(self) -> None:
        first, second = [1.0, 2.0, 4.0], [3.0, 4.0, 5.0, 6.0]
        lower = sum(a < b for a in first for b in second)       # 1: 4 pairs, 2: 4 pairs, 4: 2 pairs
        equal = sum(a == b for a in first for b in second)      # (4, 4)
        self.assertEqual((lower, equal), (10, 1))
        self.assertAlmostEqual(vargha_delaney_a12(first, second), 10.5 / 12, places=15)
        self.assertEqual(vargha_delaney_a12([1.0, 2.0], [3.0, 4.0]), 1.0)     # first always better
        self.assertEqual(vargha_delaney_a12([3.0, 4.0], [1.0, 2.0]), 0.0)
        self.assertEqual(vargha_delaney_a12([1.0, 2.0], [1.0, 2.0]), 0.5)
        self.assertEqual(vargha_delaney_a12([INF, INF], [INF]), 0.5)          # failed runs tie
        self.assertEqual(vargha_delaney_a12([1.0], [INF]), 1.0)
        # Complementarity: A12(x, y) + A12(y, x) = 1.
        self.assertAlmostEqual(vargha_delaney_a12(first, second) + vargha_delaney_a12(second, first), 1.0, places=15)

    def test_mann_whitney_against_exact_enumeration(self) -> None:
        first, second = [1.1, 2.3, 3.7, 9.0], [4.2, 5.5, 6.1, 7.4, 8.8]
        result = mann_whitney(first, second)
        greater = sum(a > b for a in first for b in second)     # U of the first sample
        self.assertEqual(result["statistic"], float(greater))
        self.assertAlmostEqual(result["p_value"], exact_two_sided_p(first, second), places=12)
        self.assertAlmostEqual(result["a12"], 1.0 - greater / 20, places=15)
        self.assertEqual((result["n_first"], result["n_second"]), (4, 5))
        separated = mann_whitney([1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0])
        self.assertAlmostEqual(separated["p_value"], 2 / math.comb(10, 5), places=12)   # the two extreme arrangements

    def test_mann_whitney_with_failed_runs_and_ties(self) -> None:
        result = mann_whitney([1.0, 2.0, INF], [3.0, INF, INF])
        self.assertTrue(0.0 < result["p_value"] <= 1.0)
        # 1 and 2 beat all three values of the second sample (6 pairs); inf ties with its two infs.
        self.assertAlmostEqual(result["a12"], (6 + 0.5 * 2) / 9, places=15)
        tied = mann_whitney([2.0, 2.0], [2.0, 2.0, 2.0])
        self.assertEqual((tied["p_value"], tied["statistic"], tied["a12"]), (1.0, 3.0, 0.5))
        self.assertEqual(mann_whitney([INF], [INF, INF])["p_value"], 1.0)

    def test_wilcoxon_signed_rank(self) -> None:
        first = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        second = [value + shift for value, shift in zip(first, (0.1, 0.2, 0.3, 0.4, 0.5, 0.6))]
        result = wilcoxon_signed_rank(first, second)
        # All six differences are negative: exact two-sided p = 2 / 2^6.
        self.assertAlmostEqual(result["p_value"], 2 / 64, places=12)
        self.assertEqual((result["statistic"], result["n_pairs"], result["n_nonzero"]), (0.0, 6, 6))
        same = wilcoxon_signed_rank(first, first)
        self.assertEqual((same["p_value"], same["n_nonzero"]), (1.0, 0))
        with self.assertRaises(ValueError):
            wilcoxon_signed_rank([1.0, 2.0], [1.0])
        with self.assertRaises(ValueError):
            wilcoxon_signed_rank([1.0, INF], [1.0, 2.0])

    def test_friedman_ranks_and_statistic(self) -> None:
        # Algorithm 0 is best on every problem, 2 is worst: ranks (1, 2, 3) in each of 4 rows.
        table = [[1.0, 2.0, 3.0], [10.0, 20.0, 30.0], [0.1, 0.5, 0.9], [5.0, 6.0, 7.0]]
        result = friedman(table)
        self.assertEqual(result["average_ranks"], [1.0, 2.0, 3.0])
        # chi2 = 12 n /(k (k+1)) * sum_j (R_j - (k+1)/2)^2 = 12*4/12 * (1 + 0 + 1) = 8.
        self.assertAlmostEqual(result["statistic"], 8.0, places=12)
        self.assertAlmostEqual(result["p_value"], math.exp(-4.0), places=12)    # chi2 with 2 dof: exp(-x/2)
        self.assertEqual((result["n_problems"], result["n_algorithms"]), (4, 3))
        tied = friedman([[1.0, 1.0, 3.0], [2.0, 2.0, 5.0]])
        self.assertEqual(tied["average_ranks"], [1.5, 1.5, 3.0])
        for table in ([[1.0, 2.0, 3.0]], [[1.0, 2.0], [3.0, 4.0]], [[1.0, 2.0, INF], [1.0, 2.0, 3.0]]):
            with self.subTest(table=table), self.assertRaises(ValueError):
                friedman(table)

    def test_holm_correction(self) -> None:
        # Sorted: 0.01, 0.03, 0.04 -> 3*0.01, max(0.03, 2*0.03), max(0.06, 1*0.04).
        adjusted = holm_correction([0.01, 0.04, 0.03])
        for got, expected in zip(adjusted, [0.03, 0.06, 0.06]):
            self.assertAlmostEqual(got, expected, places=15)
        self.assertEqual(holm_correction([0.5]), [0.5])
        self.assertEqual(holm_correction([0.6, 0.7]), [1.0, 1.0])       # capped at 1
        adjusted = holm_correction([0.02, 0.001, 0.04, 0.03, 0.2])
        self.assertTrue(all(a >= p for a, p in zip(adjusted, [0.02, 0.001, 0.04, 0.03, 0.2])))
        order = np.argsort([0.02, 0.001, 0.04, 0.03, 0.2])
        self.assertTrue(np.all(np.diff(np.array(adjusted)[order]) >= 0))   # monotone in the raw p-values
        for values in ([-0.1], [1.2], [math.nan], []):
            with self.subTest(values=values), self.assertRaises(ValueError):
                holm_correction(values)


if __name__ == "__main__":
    unittest.main()
