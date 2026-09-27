"""Tests for signals_in_the_noise.analysis.statistics."""

import pytest

from signals_in_the_noise.analysis.statistics import (
    difference_in_medians,
    exact_permutation_test,
    fdr_to_stars,
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
