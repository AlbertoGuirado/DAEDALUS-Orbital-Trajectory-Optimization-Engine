"""Statistics for multi-seed optimizer benchmarks (docs/03 section 17).

Final objectives of stochastic optimizers are skewed and may contain +inf
(a run that found no valid point), so the summaries are order statistics and
the comparisons are non-parametric (RF-6A.05). Order statistics use the
inverted empirical CDF (no interpolation): the q-quantile of n sorted values
is the value at position ceil(q n), so every reported quantile is an observed
value and +inf needs no arithmetic.

Significance tests come from scipy.stats; the multiple-comparison correction
(Holm, 1979) and the effect size (Vargha and Delaney, 2000) are implemented
here. A p-value is not an effect size: report both (docs/03 section 17.4).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import stats


def _sample(values: ArrayLike, name: str, *, allow_infinite: bool = True) -> NDArray[np.float64]:
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain real numbers") from exc
    if array.ndim != 1 or array.size < 1:
        raise ValueError(f"{name} must be a nonempty one-dimensional sample")
    if np.any(np.isnan(array)) or np.any(array == -np.inf):
        raise ValueError(f"{name} must not contain NaN or -inf")
    if not allow_infinite and not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def quantile(values: ArrayLike, q: float) -> float:
    """q-quantile, 0 <= q <= 1, by the inverted empirical CDF: sorted value at position ceil(q n).

    q = 0 returns the minimum. +inf is allowed and sorts last.
    """
    array = np.sort(_sample(values, "values"))
    if isinstance(q, bool) or not isinstance(q, (int, float, np.integer, np.floating)) or not 0.0 <= q <= 1.0:
        raise ValueError("q must be a real number in [0, 1]")
    position = max(1, math.ceil(float(q) * array.size))
    return float(array[position - 1])


def summarize(values: ArrayLike) -> dict[str, Any]:
    """Order statistics of a sample of final values (lower is better).

    Keys: n, n_finite, best (minimum), q1, median, q3, worst (maximum), iqr
    (q3 - q1, NaN if q3 is infinite), and mean and std (sample standard
    deviation, ddof = 1) over the finite values only (NaN when undefined).
    """
    array = _sample(values, "values")
    finite = array[np.isfinite(array)]
    q1, median, q3 = (quantile(array, q) for q in (0.25, 0.5, 0.75))
    return {
        "n": int(array.size),
        "n_finite": int(finite.size),
        "best": float(np.min(array)),
        "q1": q1,
        "median": median,
        "q3": q3,
        "worst": float(np.max(array)),
        "iqr": q3 - q1 if math.isfinite(q3) else math.nan,
        "mean": float(np.mean(finite)) if finite.size else math.nan,
        "std": float(np.std(finite, ddof=1)) if finite.size > 1 else math.nan,
    }


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion (Wilson, 1927).

    With p = successes/trials and z the standard normal quantile of
    (1 + confidence)/2:

        centre = (p + z^2/(2n)) / (1 + z^2/n)
        half   = z sqrt(p(1-p)/n + z^2/(4n^2)) / (1 + z^2/n)

    Unlike the normal approximation it stays inside [0, 1] and is usable for
    0 or n successes. Returns (low, high).
    """
    for value, name in ((successes, "successes"), (trials, "trials")):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{name} must be an integer")
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("trials must be >= 1 and 0 <= successes <= trials")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie in (0, 1)")
    z = float(stats.norm.ppf((1.0 + confidence) / 2.0))
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    centre = (proportion + z * z / (2 * trials)) / denominator
    half = z * math.sqrt(proportion * (1 - proportion) / trials + z * z / (4 * trials * trials)) / denominator
    # The exact interval always contains the observed proportion (its limits are 0 and 1 for 0
    # and n successes); the outer min/max keep that property under floating-point rounding.
    return min(max(0.0, centre - half), proportion), max(min(1.0, centre + half), proportion)


def evaluations_to_target(best_so_far: ArrayLike, target: float) -> int | None:
    """Number of evaluations after which the best-so-far value first is <= target, or None.

    ``best_so_far[k]`` is the best value after evaluation k (0-based), so the
    result is k + 1 for the first such k.
    """
    history = _sample(best_so_far, "best_so_far")
    reached = np.nonzero(history <= target)[0]
    return int(reached[0]) + 1 if reached.size else None


def convergence_band(histories: Sequence[ArrayLike], evaluations: ArrayLike | None = None,
                     quantiles: Sequence[float] = (0.25, 0.5, 0.75)) -> dict[str, Any]:
    """Quantiles across runs of the best-so-far value at given evaluation counts.

    ``histories``: one best-so-far array per run, all of the same length B
    (equal budgets). ``evaluations``: 1-based evaluation counts in [1, B] at
    which to report; None reports every evaluation. Returns
    {"evaluations": [...], "quantiles": {q: [...]}} with inverted-CDF
    quantiles, so +inf (no valid point yet) is reported as such.
    """
    if len(histories) < 1:
        raise ValueError("at least one history is required")
    stacked = np.array([_sample(history, "history") for history in histories])
    if stacked.ndim != 2:
        raise ValueError("all histories must have the same length")
    budget = stacked.shape[1]
    if evaluations is None:
        counts = np.arange(1, budget + 1)
    else:
        counts = np.asarray(evaluations)
        if counts.ndim != 1 or counts.size < 1 or counts.dtype.kind not in "iu" \
                or np.any(counts < 1) or np.any(counts > budget):
            raise ValueError("evaluations must be integer counts in [1, budget]")
    columns = np.sort(stacked[:, counts - 1], axis=0)
    runs = stacked.shape[0]
    result: dict[str, Any] = {"evaluations": counts.astype(int).tolist(), "quantiles": {}}
    for q in quantiles:
        if not 0.0 <= q <= 1.0:
            raise ValueError("quantiles must lie in [0, 1]")
        position = max(1, math.ceil(q * runs))
        result["quantiles"][float(q)] = columns[position - 1].tolist()
    return result


def vargha_delaney_a12(first: ArrayLike, second: ArrayLike) -> float:
    """A12 = P(X < Y) + P(X = Y)/2 for X from ``first`` and Y from ``second`` (lower is better).

    Probability that a run of the first algorithm ends with a better (lower)
    value than a run of the second; 0.5 means no difference. Conventional
    magnitudes (Vargha and Delaney, 2000): |A12 - 0.5| >= 0.06 small, 0.14
    medium, 0.21 large. +inf is allowed (two +inf values tie).
    """
    x, y = _sample(first, "first"), _sample(second, "second")
    lower = np.count_nonzero(x[:, None] < y[None, :])
    equal = np.count_nonzero(x[:, None] == y[None, :])
    return float((lower + 0.5 * equal) / (x.size * y.size))


def mann_whitney(first: ArrayLike, second: ArrayLike) -> dict[str, Any]:
    """Two-sided Mann-Whitney U test for two independent samples, with A12.

    Null hypothesis: P(X < Y) = P(X > Y). ``statistic`` is U of the first
    sample (scipy.stats.mannwhitneyu, exact for small samples without ties,
    otherwise the tie-corrected normal approximation). +inf values are ranked
    as the largest values.
    """
    x, y = _sample(first, "first"), _sample(second, "second")
    # Ranks only depend on the order: replace +inf by a value above every finite one.
    finite = np.concatenate((x[np.isfinite(x)], y[np.isfinite(y)]))
    ceiling = (float(np.max(finite)) if finite.size else 0.0) + 1.0
    x_ranked, y_ranked = np.where(np.isfinite(x), x, ceiling), np.where(np.isfinite(y), y, ceiling)
    if np.all(x_ranked == x_ranked[0]) and np.all(y_ranked == x_ranked[0]):
        statistic, p_value = x.size * y.size / 2.0, 1.0      # all values tied: no evidence of a difference
    else:
        outcome = stats.mannwhitneyu(x_ranked, y_ranked, alternative="two-sided")
        statistic, p_value = float(outcome.statistic), float(outcome.pvalue)
    return {"test": "mann-whitney-u", "statistic": statistic, "p_value": p_value,
            "a12": vargha_delaney_a12(x, y), "n_first": int(x.size), "n_second": int(y.size)}


def wilcoxon_signed_rank(first: ArrayLike, second: ArrayLike) -> dict[str, Any]:
    """Two-sided Wilcoxon signed-rank test for paired finite samples (e.g. per-scenario medians).

    Pairs with zero difference are discarded (Wilcoxon's method). With no
    nonzero difference the p-value is 1.
    """
    x = _sample(first, "first", allow_infinite=False)
    y = _sample(second, "second", allow_infinite=False)
    if x.shape != y.shape:
        raise ValueError("paired samples must have the same length")
    differences = x - y
    if not np.any(differences != 0):
        return {"test": "wilcoxon-signed-rank", "statistic": 0.0, "p_value": 1.0, "n_pairs": int(x.size),
                "n_nonzero": 0}
    outcome = stats.wilcoxon(x, y, zero_method="wilcox", alternative="two-sided")
    return {"test": "wilcoxon-signed-rank", "statistic": float(outcome.statistic), "p_value": float(outcome.pvalue),
            "n_pairs": int(x.size), "n_nonzero": int(np.count_nonzero(differences))}


def friedman(matrix: ArrayLike) -> dict[str, Any]:
    """Friedman test for k >= 3 algorithms (columns) on n >= 2 problems (rows), with average ranks.

    Entries are finite scores, lower is better; rank 1 is the best in its
    row and ties share the average rank.
    """
    try:
        table = np.asarray(matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("matrix must contain real numbers") from exc
    if table.ndim != 2 or table.shape[0] < 2 or table.shape[1] < 3 or not np.all(np.isfinite(table)):
        raise ValueError("matrix must be finite with at least 2 rows (problems) and 3 columns (algorithms)")
    ranks = np.array([stats.rankdata(row, method="average") for row in table])
    outcome = stats.friedmanchisquare(*table.T)
    return {"test": "friedman", "statistic": float(outcome.statistic), "p_value": float(outcome.pvalue),
            "average_ranks": ranks.mean(axis=0).tolist(), "n_problems": int(table.shape[0]),
            "n_algorithms": int(table.shape[1])}


def holm_correction(p_values: ArrayLike) -> list[float]:
    """Holm-Bonferroni adjusted p-values (Holm, 1979), in the input order.

    With the m p-values sorted ascending, p_(1) <= ... <= p_(m):
        adjusted_(i) = max_{j <= i} min(1, (m - j + 1) p_(j)).
    Rejecting when adjusted <= alpha controls the family-wise error rate at
    alpha without assuming independence.
    """
    p = _sample(p_values, "p_values", allow_infinite=False)
    if np.any(p < 0.0) or np.any(p > 1.0):
        raise ValueError("p-values must lie in [0, 1]")
    order = np.argsort(p, kind="stable")
    m = p.size
    adjusted = np.empty(m)
    running = 0.0
    for position, index in enumerate(order):
        running = max(running, min(1.0, (m - position) * float(p[index])))
        adjusted[index] = running
    return adjusted.tolist()
