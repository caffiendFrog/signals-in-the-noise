import itertools
import logging
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import stats

logger = logging.getLogger(__name__)

Alternative = Literal["two-sided", "greater", "less"]
TwoSampleStatistic = Callable[[np.ndarray, np.ndarray], float]


def fdr_to_stars(fdr: float) -> str:
    """Convert an FDR q-value to a significance-star annotation string.

    Thresholds follow the conventional three-tier scheme:

    - ``< 0.01`` → ``' ***'``
    - ``< 0.05`` → ``' **'``
    - ``< 0.1``  → ``' *'``
    - ``>= 0.1`` → ``''`` (empty string, no annotation)

    Args:
        fdr: FDR-corrected q-value in the range ``[0, 1]``.

    Returns:
        A star-annotation string (with a leading space when non-empty) suitable
        for appending to axis tick labels or table cells.
    """
    if fdr < 0.01:
        return " ***"
    if fdr < 0.05:
        return " **"
    if fdr < 0.1:
        return " *"
    return ""


def difference_in_means(a: np.ndarray, b: np.ndarray) -> float:
    """Return ``mean(a) - mean(b)``."""
    return float(np.mean(a) - np.mean(b))


def difference_in_medians(a: np.ndarray, b: np.ndarray) -> float:
    """Return ``median(a) - median(b)``."""
    return float(np.median(a) - np.median(b))


@dataclass(frozen=True)
class PermutationTestResult:
    """Outcome of :func:`exact_permutation_test`.

    Attributes:
        observed: Statistic for the observed group labels.
        p_value: Fraction of relabelings at least as extreme as ``observed``
            (the observed labeling is one of them, so ``p_value >= 1 / n_permutations``).
        n_permutations: Number of distinct relabelings enumerated.
        alternative: Direction of the test.
    """

    observed: float
    p_value: float
    n_permutations: int
    alternative: Alternative

    @property
    def min_achievable_p_value(self) -> float:
        """Smallest p-value this design can produce."""
        return 1.0 / self.n_permutations


def exact_permutation_test(
    group_a: Sequence[float],
    group_b: Sequence[float],
    *,
    statistic: TwoSampleStatistic = difference_in_means,
    alternative: Alternative = "two-sided",
    max_permutations: int = 1_000_000,
) -> PermutationTestResult:
    """Exact two-sample permutation test enumerating every relabeling of the pooled values.

    Intended for small donor-level designs, e.g. 4 BRCA1 vs 8 WT donors gives
    ``C(12, 4) = 495`` relabelings and a minimum p-value of ``1/495``.

    Args:
        group_a: Values for the group of interest (e.g. BRCA1 donors).
        group_b: Values for the reference group (e.g. WT donors).
        statistic: Function of ``(a, b)`` returning a scalar. Defaults to
            :func:`difference_in_means`.
        alternative: ``"greater"`` tests ``statistic > 0``, ``"less"`` tests
            ``statistic < 0``, ``"two-sided"`` compares absolute values.
        max_permutations: Refuse to enumerate more relabelings than this.

    Returns:
        A :class:`PermutationTestResult`.

    Raises:
        ValueError: If either group is empty, ``alternative`` is unknown, or the
            number of relabelings exceeds ``max_permutations``.
    """
    a = np.asarray(group_a, dtype=float)
    b = np.asarray(group_b, dtype=float)
    if a.size == 0 or b.size == 0:
        raise ValueError("Both groups must contain at least one value.")
    if alternative not in ("two-sided", "greater", "less"):
        raise ValueError(f"Unknown alternative {alternative!r}.")

    pooled = np.concatenate([a, b])
    n_total, n_a = pooled.size, a.size
    n_permutations = math.comb(n_total, n_a)
    if n_permutations > max_permutations:
        raise ValueError(
            f"{n_permutations} relabelings exceed max_permutations={max_permutations}."
        )

    observed = statistic(a, b)
    null = np.empty(n_permutations)
    for i, members in enumerate(itertools.combinations(range(n_total), n_a)):
        in_a = np.zeros(n_total, dtype=bool)
        in_a[list(members)] = True
        null[i] = statistic(pooled[in_a], pooled[~in_a])

    # Tolerance so relabelings tied with the observed statistic count as extreme
    # despite floating-point rounding.
    tol = 1e-12 * max(1.0, abs(observed))
    if alternative == "greater":
        extreme = null >= observed - tol
    elif alternative == "less":
        extreme = null <= observed + tol
    else:
        extreme = np.abs(null) >= abs(observed) - tol

    result = PermutationTestResult(
        observed=observed,
        p_value=float(extreme.mean()),
        n_permutations=n_permutations,
        alternative=alternative,
    )
    logger.debug("exact permutation test: %s", result)
    return result


def bootstrap_mean_difference_ci(
    group_a: Sequence[float],
    group_b: Sequence[float],
    *,
    confidence: float = 0.95,
    n_resamples: int = 10_000,
    seed: int | None = None,
) -> tuple[float, float]:
    """Percentile bootstrap CI for ``mean(a) - mean(b)``, resampling donors within each group."""
    a = np.asarray(group_a, dtype=float)
    b = np.asarray(group_b, dtype=float)
    rng = np.random.default_rng(seed)
    a_means = a[rng.integers(0, a.size, size=(n_resamples, a.size))].mean(axis=1)
    b_means = b[rng.integers(0, b.size, size=(n_resamples, b.size))].mean(axis=1)
    tail = (1.0 - confidence) / 2.0
    low, high = np.quantile(a_means - b_means, [tail, 1.0 - tail])
    return float(low), float(high)


def pooled_t_interval(
    group_a: Sequence[float], group_b: Sequence[float], *, confidence: float = 0.95
) -> tuple[float, float]:
    """Pooled-variance t interval for ``mean(a) - mean(b)``."""
    a = np.asarray(group_a, dtype=float)
    b = np.asarray(group_b, dtype=float)
    df = a.size + b.size - 2
    pooled_var = ((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1)) / df
    half_width = stats.t.ppf(0.5 + confidence / 2.0, df) * math.sqrt(
        pooled_var * (1.0 / a.size + 1.0 / b.size)
    )
    difference = float(a.mean() - b.mean())
    return difference - half_width, difference + half_width
