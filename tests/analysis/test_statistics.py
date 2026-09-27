"""Tests for signals_in_the_noise.analysis.statistics."""

import pytest

from signals_in_the_noise.analysis.statistics import (
    bootstrap_mean_difference_ci,
    difference_in_medians,
    exact_permutation_test,
    fdr_to_stars,
    pooled_t_interval,
)


# ---------------------------------------------------------------------------
# Return-type contract
# ---------------------------------------------------------------------------


def test_fdr_to_stars_returns_string():
    assert isinstance(fdr_to_stars(0.05), str)


# ---------------------------------------------------------------------------
# Threshold boundary tests
# ---------------------------------------------------------------------------


def test_fdr_to_stars_below_0_01_returns_three_stars():
    assert fdr_to_stars(0.001) == " ***"


def test_fdr_to_stars_exactly_0_01_returns_two_stars():
    assert fdr_to_stars(0.01) == " **"


def test_fdr_to_stars_between_0_01_and_0_05_returns_two_stars():
    assert fdr_to_stars(0.03) == " **"


def test_fdr_to_stars_exactly_0_05_returns_one_star():
    assert fdr_to_stars(0.05) == " *"


def test_fdr_to_stars_between_0_05_and_0_1_returns_one_star():
    assert fdr_to_stars(0.07) == " *"


def test_fdr_to_stars_exactly_0_1_returns_empty_string():
    assert fdr_to_stars(0.1) == ""


def test_fdr_to_stars_above_0_1_returns_empty_string():
    assert fdr_to_stars(0.5) == ""


def test_fdr_to_stars_at_zero_returns_three_stars():
    assert fdr_to_stars(0.0) == " ***"


def test_fdr_to_stars_at_one_returns_empty_string():
    assert fdr_to_stars(1.0) == ""


# ---------------------------------------------------------------------------
# exact_permutation_test
# ---------------------------------------------------------------------------

SEPARATED_A = [10.0, 11.0, 12.0, 13.0]
SEPARATED_B = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]


def test_permutation_4_vs_8_enumerates_495_relabelings():
    result = exact_permutation_test(SEPARATED_A, SEPARATED_B)
    assert result.n_permutations == 495
    assert result.min_achievable_p_value == pytest.approx(1 / 495)


def test_permutation_complete_separation_reaches_minimum_p_value():
    result = exact_permutation_test(SEPARATED_A, SEPARATED_B, alternative="greater")
    assert result.p_value == pytest.approx(1 / 495)


def test_permutation_less_is_not_significant_when_a_is_larger():
    result = exact_permutation_test(SEPARATED_A, SEPARATED_B, alternative="less")
    assert result.p_value == pytest.approx(1.0)


def test_permutation_identical_values_give_p_one():
    result = exact_permutation_test([5.0] * 4, [5.0] * 8)
    assert result.observed == 0.0
    assert result.p_value == pytest.approx(1.0)


def test_permutation_reports_observed_statistic():
    result = exact_permutation_test([3.0, 5.0], [1.0, 1.0, 2.0], statistic=difference_in_medians)
    assert result.observed == pytest.approx(3.0)


def test_permutation_rejects_empty_group():
    with pytest.raises(ValueError, match="at least one"):
        exact_permutation_test([], [1.0])


def test_permutation_rejects_unknown_alternative():
    with pytest.raises(ValueError, match="alternative"):
        exact_permutation_test([1.0], [2.0], alternative="bigger")


def test_permutation_refuses_intractable_designs():
    with pytest.raises(ValueError, match="max_permutations"):
        exact_permutation_test(range(10), range(10), max_permutations=100)


# ---------------------------------------------------------------------------
# Confidence intervals for a difference in means
# ---------------------------------------------------------------------------


def test_bootstrap_ci_brackets_observed_difference():
    low, high = bootstrap_mean_difference_ci(SEPARATED_A, SEPARATED_B, seed=0)
    assert low < 11.5 - 4.5 < high


def test_bootstrap_ci_is_reproducible_with_seed():
    first = bootstrap_mean_difference_ci(SEPARATED_A, SEPARATED_B, seed=7)
    assert bootstrap_mean_difference_ci(SEPARATED_A, SEPARATED_B, seed=7) == first


def test_bootstrap_ci_collapses_for_constant_groups():
    assert bootstrap_mean_difference_ci([2.0] * 4, [1.0] * 8, seed=0) == pytest.approx((1.0, 1.0))


def test_pooled_t_interval_is_centred_on_difference():
    low, high = pooled_t_interval(SEPARATED_A, SEPARATED_B)
    assert (low + high) / 2 == pytest.approx(7.0)


def test_pooled_t_interval_matches_textbook_value():
    # sd_pooled^2 = (3 * 1.667 + 7 * 6.0) / 10 = 4.7; se = sqrt(4.7 * 3 / 8); t(0.975, 10) = 2.228
    low, high = pooled_t_interval(SEPARATED_A, SEPARATED_B)
    assert (high - low) / 2 == pytest.approx(2.2281 * (4.7 * 3 / 8) ** 0.5, rel=1e-3)
