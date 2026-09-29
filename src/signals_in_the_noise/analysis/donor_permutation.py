"""Exact permutation tests with the donor as the unit of replication.

Group sizes stay fixed. For 4 versus 8 donors every assignment of the labels
is equally likely, so the p-value is a multiple of 1 over the number of
assignments, and the observed assignment is counted.
"""

import functools
import math
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.stats import false_discovery_control, norm

_ALTERNATIVES = ("two-sided", "greater", "less")


@dataclass(frozen=True)
class PermutationResult:
    """One exact permutation test."""

    statistic: float
    p_value: float
    n_permutations: int
    n_extreme: int
    alpha: float
    alternative: str
    significant: bool


@dataclass(frozen=True)
class ShiftInterval:
    """Confidence interval from inverting a shift permutation test."""

    estimate: float
    low: float
    high: float
    level: float
    alternative: str
    step: float
    bootstrap_low: float | None = None
    bootstrap_high: float | None = None


@dataclass(frozen=True)
class MDEResult:
    """Minimum effect with the requested power, on the log-odds scale and in points."""

    delta_log_odds: float
    points: float
    normal_points: float
    power: float
    power_at_zero: float
    attained: bool
    pooled_log_odds: float


def mean_diff(values, labels, positive: str) -> float:
    """Unweighted mean of ``positive`` minus the unweighted mean of the other group."""
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    in_group = labels == positive
    if in_group.sum() == 0 or (~in_group).sum() == 0:
        raise ValueError("Both groups must be present.")
    return float(values[in_group].mean() - values[~in_group].mean())


def attainable_alpha(n_perm: int, alpha: float) -> float:
    """Largest multiple of ``1/n_perm`` that does not exceed ``alpha``."""
    if n_perm < 1:
        raise ValueError(f"n_perm must be positive, got {n_perm}.")
    return math.floor(alpha * n_perm + 1e-12) / n_perm


def bh_adjust(p_values) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, in the input order."""
    return np.asarray(false_discovery_control(np.asarray(p_values, dtype=float)), dtype=float)


def empirical_log_odds(removed, total) -> np.ndarray:
    """Log-odds of ``(removed + 0.5) / (total + 1)`` so a zero count stays finite."""
    removed = np.asarray(removed, dtype=float)
    total = np.asarray(total, dtype=float)
    proportion = (removed + 0.5) / (total + 1.0)
    return np.log(proportion / (1.0 - proportion))


def expit(values) -> np.ndarray:
    """Inverse logit."""
    values = np.asarray(values, dtype=float)
    return 1.0 / (1.0 + np.exp(-values))


def exact_permutation_test(
    values,
    labels,
    *,
    positive: str,
    alternative: str = "two-sided",
    alpha: float = 0.01,
) -> PermutationResult:
    """Exact test of an unweighted mean difference, conditional on the group sizes.

    ``alternative`` is ``"two-sided"``, ``"greater"``, or ``"less"``. The
    two-sided p-value counts assignments whose absolute mean difference is at
    least as large as the observed one. The observed assignment is included, so
    the minimum p-value is ``1/n_permutations``.
    """
    _check_alternative(alternative)
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    if len(values) != len(labels):
        raise ValueError("values and labels must have the same length.")
    positive_index = tuple(sorted(np.flatnonzero(labels == positive)))
    if len(positive_index) == 0 or len(positive_index) == len(labels):
        raise ValueError("Both groups must be present.")
    combos = _combos(len(values), len(positive_index))
    stats = _stats_for_combos(values, combos)
    match = np.all(combos == positive_index, axis=1)
    observed = float(stats[match][0])
    n_extreme, p_value = _tail_count(stats, observed, alternative)
    return PermutationResult(
        statistic=observed,
        p_value=p_value,
        n_permutations=int(combos.shape[0]),
        n_extreme=n_extreme,
        alpha=alpha,
        alternative=alternative,
        significant=bool(p_value <= alpha),
    )


def permutation_ci(
    values,
    labels,
    *,
    positive: str,
    level: float = 0.99,
    alternative: str = "two-sided",
    step: float = 0.001,
    n_bootstrap: int = 2000,
    rng: np.random.Generator | None = None,
) -> ShiftInterval:
    """Confidence interval for the mean difference by inverting the shift test.

    A hypothesized difference is inside the interval when shifting the
    ``positive`` group by that amount is not rejected at ``1 - level``. A
    one-sided alternative leaves the other end unbounded. A percentile
    bootstrap interval, resampling donors inside each group, is attached when
    ``n_bootstrap`` is positive. With four donors that bootstrap under-covers,
    so it is a secondary check.
    """
    _check_alternative(alternative)
    if not 0 < level < 1:
        raise ValueError(f"level must be in (0, 1), got {level}.")
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    alpha = 1.0 - level
    estimate = mean_diff(values, labels, positive)
    in_group = labels == positive

    def accepts(delta: float) -> bool:
        shifted = values.copy()
        shifted[in_group] = shifted[in_group] - delta
        result = exact_permutation_test(
            shifted, labels, positive=positive, alternative=alternative, alpha=alpha
        )
        return not result.significant

    if alternative == "two-sided":
        low = _walk_until_rejected(estimate, -step, accepts, -1.0)
        high = _walk_until_rejected(estimate, step, accepts, 1.0)
    elif alternative == "greater":
        low = _walk_until_rejected(estimate, -step, accepts, -1.0)
        high = math.inf
    else:
        low = -math.inf
        high = _walk_until_rejected(estimate, step, accepts, 1.0)

    bootstrap_low = bootstrap_high = None
    if n_bootstrap > 0:
        bootstrap_low, bootstrap_high = _bootstrap_interval(
            values, labels, positive, level, alternative, n_bootstrap, rng
        )
    return ShiftInterval(
        estimate=estimate,
        low=low,
        high=high,
        level=level,
        alternative=alternative,
        step=step,
        bootstrap_low=bootstrap_low,
        bootstrap_high=bootstrap_high,
    )


def leave_one_donor_out(
    values,
    labels,
    donor_ids,
    *,
    positive: str,
    level: float = 0.98,
    step: float = 0.001,
) -> list[dict]:
    """Estimate and a two-sided interval after dropping each donor in turn.

    The interval uses ``level`` for every dropped donor, including when the
    full-data test uses a stricter level. Dropping a member of the smaller
    group leaves fewer assignments than dropping a member of the larger group.
    """
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    donor_ids = np.asarray(donor_ids)
    rows = []
    for index, donor in enumerate(donor_ids):
        keep = np.ones(len(values), dtype=bool)
        keep[index] = False
        kept_values = values[keep]
        kept_labels = labels[keep]
        interval = permutation_ci(
            kept_values,
            kept_labels,
            positive=positive,
            level=level,
            alternative="two-sided",
            step=step,
            n_bootstrap=0,
        )
        n_positive = int((kept_labels == positive).sum())
        rows.append(
            {
                "dropped": donor,
                "dropped_group": labels[index],
                "estimate": interval.estimate,
                "ci_low": interval.low,
                "ci_high": interval.high,
                "n_permutations": math.comb(int(keep.sum()), n_positive),
                "n_positive": n_positive,
                "n_negative": int(keep.sum()) - n_positive,
            }
        )
    return rows


def minimum_detectable_effect(
    removed,
    total,
    labels,
    *,
    positive: str,
    alpha: float,
    alternative: str = "two-sided",
    power: float = 0.8,
    delta_step: float = 0.05,
    delta_max: float = 5.0,
) -> MDEResult:
    """Smallest log-odds shift with at least ``power`` under exact enumeration.

    Donor counts are converted to log-odds, each genotype is centered on its
    own mean, and the pooled mean is added back. That removes the observed
    genotype difference and keeps the donor spread. Power at a shift is the
    share of label assignments that the same exact test rejects after the
    shift is added to the positive group. The result in percentage points is
    ``100 * (expit(m + shift) - expit(m))`` at the pooled log-odds ``m``.

    The normal approximation uses the pooled within-genotype standard deviation
    on the log-odds scale. It is a cross-check, not the reported effect.
    """
    _check_alternative(alternative)
    log_odds = empirical_log_odds(removed, total)
    labels = np.asarray(labels)
    centered = _center_groups(log_odds, labels)
    n_positive = int((labels == positive).sum())
    if n_positive == 0 or n_positive == len(labels):
        raise ValueError("Both groups must be present.")
    power_at_zero = _randomization_power(centered, n_positive, 0.0, alternative, alpha)
    delta, attained, achieved = _smallest_shift(
        centered, n_positive, alternative, alpha, power, delta_step, delta_max
    )
    pooled = float(np.mean(log_odds))
    return MDEResult(
        delta_log_odds=delta,
        points=_to_points(pooled, delta),
        normal_points=_normal_points(log_odds, labels, positive, alpha, alternative, power, pooled),
        power=achieved,
        power_at_zero=power_at_zero,
        attained=attained,
        pooled_log_odds=pooled,
    )


def mde_resample_spread(
    removed,
    total,
    labels,
    *,
    positive: str,
    alpha: float,
    alternative: str = "two-sided",
    power: float = 0.8,
    n_resamples: int = 20,
    rng: np.random.Generator | None = None,
    delta_step: float = 0.1,
) -> dict[str, float]:
    """Min, median, and max MDE from resampling donors inside each genotype."""
    rng = np.random.default_rng() if rng is None else rng
    removed = np.asarray(removed, dtype=float)
    total = np.asarray(total, dtype=float)
    labels = np.asarray(labels)
    positive_index = np.flatnonzero(labels == positive)
    negative_index = np.flatnonzero(labels != positive)
    negative_label = labels[negative_index[0]]
    points = []
    for _ in range(n_resamples):
        take_positive = rng.choice(positive_index, size=len(positive_index), replace=True)
        take_negative = rng.choice(negative_index, size=len(negative_index), replace=True)
        take = np.concatenate([take_positive, take_negative])
        resampled_labels = np.array(
            [positive] * len(take_positive) + [negative_label] * len(take_negative)
        )
        result = minimum_detectable_effect(
            removed[take],
            total[take],
            resampled_labels,
            positive=positive,
            alpha=alpha,
            alternative=alternative,
            power=power,
            delta_step=delta_step,
        )
        points.append(result.points)
    array = np.asarray(points, dtype=float)
    return {"min": float(array.min()), "median": float(np.median(array)), "max": float(array.max())}


def log_odds_shift_points(pooled_log_odds: float, delta: float) -> float:
    """Percentage points for a log-odds shift away from ``pooled_log_odds``."""
    return _to_points(pooled_log_odds, delta)


def _to_points(pooled_log_odds: float, delta: float) -> float:
    return float(100.0 * (expit(pooled_log_odds + delta) - expit(pooled_log_odds)))


def _check_alternative(alternative: str) -> None:
    if alternative not in _ALTERNATIVES:
        raise ValueError(f"alternative must be one of {_ALTERNATIVES}, got {alternative!r}.")


@functools.lru_cache(maxsize=8)
def _combos(n: int, k: int) -> np.ndarray:
    return np.asarray(list(combinations(range(n), k)), dtype=int)


def _stats_for_combos(values: np.ndarray, combos: np.ndarray) -> np.ndarray:
    k = combos.shape[1]
    positive_sum = values[combos].sum(axis=1)
    negative_n = len(values) - k
    return positive_sum / k - (values.sum() - positive_sum) / negative_n


def _tail_count(stats: np.ndarray, observed: float, alternative: str) -> tuple[int, float]:
    if alternative == "two-sided":
        extreme = np.abs(stats) >= abs(observed)
    elif alternative == "greater":
        extreme = stats >= observed
    else:
        extreme = stats <= observed
    n_extreme = int(extreme.sum())
    return n_extreme, n_extreme / len(stats)


def _walk_until_rejected(start: float, step: float, accepts, limit: float) -> float:
    """Last grid point, walking from ``start`` by ``step``, that is still accepted."""
    current = start
    if not accepts(current):
        return current
    nxt = current + step
    while (step < 0 and nxt >= limit) or (step > 0 and nxt <= limit):
        if not accepts(nxt):
            return current
        current = nxt
        nxt = current + step
    return current


def _bootstrap_interval(values, labels, positive, level, alternative, n_bootstrap, rng):
    rng = np.random.default_rng() if rng is None else rng
    positive_values = values[labels == positive]
    negative_values = values[labels != positive]
    stats = np.empty(n_bootstrap)
    for index in range(n_bootstrap):
        drawn_positive = rng.choice(positive_values, size=len(positive_values), replace=True)
        drawn_negative = rng.choice(negative_values, size=len(negative_values), replace=True)
        stats[index] = drawn_positive.mean() - drawn_negative.mean()
    alpha = 1.0 - level
    if alternative == "two-sided":
        low, high = np.quantile(stats, [alpha / 2, 1 - alpha / 2])
        return float(low), float(high)
    if alternative == "greater":
        return float(np.quantile(stats, alpha)), math.inf
    return -math.inf, float(np.quantile(stats, 1 - alpha))


def _center_groups(values: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Subtract each group's mean and add back the pooled mean."""
    centered = np.asarray(values, dtype=float).copy()
    pooled = float(centered.mean())
    for group in np.unique(labels):
        mask = labels == group
        centered[mask] = values[mask] - values[mask].mean() + pooled
    return centered


def _randomization_power(
    centered: np.ndarray, n_positive: int, delta: float, alternative: str, alpha: float
) -> float:
    """Share of label assignments rejected after shifting the positive group by ``delta``."""
    n = len(centered)
    combos = _combos(n, n_positive)
    n_assignments = combos.shape[0]
    shifted = np.repeat(centered[None, :], n_assignments, axis=0)
    rows = np.repeat(np.arange(n_assignments), n_positive)
    shifted[rows, combos.ravel()] += delta
    positive_sum = shifted[:, combos].sum(axis=2)
    total = shifted.sum(axis=1)
    stats = positive_sum / n_positive - (total[:, None] - positive_sum) / (n - n_positive)
    observed = np.diag(stats)
    if alternative == "two-sided":
        extreme = np.abs(stats) >= np.abs(observed)[:, None]
    elif alternative == "greater":
        extreme = stats >= observed[:, None]
    else:
        extreme = stats <= observed[:, None]
    p_values = extreme.mean(axis=1)
    return float((p_values <= alpha).mean())


def _smallest_shift(centered, n_positive, alternative, alpha, power, delta_step, delta_max):
    delta = 0.0
    achieved = _randomization_power(centered, n_positive, delta, alternative, alpha)
    if achieved >= power:
        return 0.0, True, achieved
    while delta < delta_max - 1e-12:
        delta = min(delta_max, delta + delta_step)
        achieved = _randomization_power(centered, n_positive, delta, alternative, alpha)
        if achieved >= power:
            fine_step = delta_step / 10
            fine = np.arange(delta - delta_step + fine_step, delta + fine_step / 2, fine_step)
            for candidate in fine:
                candidate_power = _randomization_power(
                    centered, n_positive, float(candidate), alternative, alpha
                )
                if candidate_power >= power:
                    return float(candidate), True, candidate_power
            return delta, True, achieved
    return delta, False, achieved


def _normal_points(log_odds, labels, positive, alpha, alternative, power, pooled) -> float:
    parts = []
    degrees = 0
    for group in np.unique(labels):
        group_values = log_odds[labels == group]
        if len(group_values) < 2:
            return math.nan
        parts.append((len(group_values) - 1) * float(group_values.var(ddof=1)))
        degrees += len(group_values) - 1
    scale = math.sqrt(sum(parts) / degrees)
    n_positive = int((labels == positive).sum())
    n_negative = len(labels) - n_positive
    z_alpha = norm.ppf(1 - alpha / 2) if alternative == "two-sided" else norm.ppf(1 - alpha)
    delta = (z_alpha + norm.ppf(power)) * scale * math.sqrt(1 / n_positive + 1 / n_negative)
    return _to_points(pooled, float(delta))
