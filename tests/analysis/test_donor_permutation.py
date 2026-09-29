"""Exact donor permutation tests on synthetic vectors."""

import math

import numpy as np
import pytest

from signals_in_the_noise.analysis.donor_permutation import (
    attainable_alpha,
    bh_adjust,
    exact_permutation_test,
    leave_one_donor_out,
    minimum_detectable_effect,
    permutation_ci,
)

_BRCA1_WT = np.array(["BRCA1"] * 4 + ["WT"] * 8)


def test_tied_values_have_a_zero_statistic_and_a_point_interval():
    values = np.full(12, 0.5)
    result = exact_permutation_test(values, _BRCA1_WT, positive="BRCA1", alpha=0.01)
    assert result.statistic == 0
    assert result.p_value == 1
    assert result.n_permutations == 495
    assert not result.significant
    interval = permutation_ci(
        values, _BRCA1_WT, positive="BRCA1", level=0.99, step=0.01, n_bootstrap=0
    )
    assert interval.low == 0
    assert interval.high == 0


def test_separated_groups_reach_the_minimum_p_value_and_exclude_zero():
    values = np.array([1.0] * 4 + [0.0] * 8)
    result = exact_permutation_test(values, _BRCA1_WT, positive="BRCA1", alpha=0.01)
    assert result.p_value == pytest.approx(1 / 495)
    assert result.n_extreme == 1
    assert result.statistic == pytest.approx(1)
    assert result.significant
    greater = exact_permutation_test(
        values, _BRCA1_WT, positive="BRCA1", alternative="greater", alpha=0.01
    )
    lesser = exact_permutation_test(
        values, _BRCA1_WT, positive="BRCA1", alternative="less", alpha=0.01
    )
    assert greater.p_value == pytest.approx(1 / 495)
    assert greater.significant
    assert lesser.p_value == 1
    assert not lesser.significant
    interval = permutation_ci(
        values, _BRCA1_WT, positive="BRCA1", level=0.99, step=0.01, n_bootstrap=0
    )
    assert interval.low > 0
    one_sided = permutation_ci(
        values,
        _BRCA1_WT,
        positive="BRCA1",
        level=0.99,
        alternative="greater",
        step=0.01,
        n_bootstrap=0,
    )
    assert one_sided.low > 0
    assert math.isinf(one_sided.high)


def test_swapping_which_group_is_high_negates_the_statistic_and_keeps_the_two_sided_p():
    high_brca1 = np.array([1.0] * 4 + [0.0] * 8)
    high_wt = np.array([0.0] * 4 + [1.0] * 8)
    left = exact_permutation_test(high_brca1, _BRCA1_WT, positive="BRCA1")
    right = exact_permutation_test(high_wt, _BRCA1_WT, positive="BRCA1")
    assert right.statistic == pytest.approx(-left.statistic)
    assert right.p_value == pytest.approx(left.p_value)
    assert left.statistic > 0
    assert right.statistic < 0


def test_leave_one_out_permutations_depend_on_which_group_is_dropped():
    values = np.array([1.0] * 4 + [0.0] * 8)
    donors = [f"B{i}" for i in range(4)] + [f"W{i}" for i in range(8)]
    rows = leave_one_donor_out(values, _BRCA1_WT, donors, positive="BRCA1", level=0.98, step=0.05)
    by_donor = {row["dropped"]: row for row in rows}
    assert by_donor["B0"]["n_permutations"] == math.comb(11, 3)
    assert by_donor["W0"]["n_permutations"] == math.comb(11, 4)


def test_zero_within_group_spread_makes_the_exact_mde_one_grid_step_and_the_normal_mde_zero():
    labels = np.array(["BRCA1"] * 4 + ["WT"] * 4)
    total = np.full(8, 100.0)
    flat = np.full(8, 10.0)
    step = 0.2
    flat_mde = minimum_detectable_effect(
        flat, total, labels, positive="BRCA1", alpha=0.05, delta_step=step
    )
    assert flat_mde.power_at_zero == 0
    assert flat_mde.normal_points == pytest.approx(0)
    assert flat_mde.attained
    assert 0 < flat_mde.delta_log_odds <= step
    assert flat_mde.points > 0
    noisy = np.array([0, 100, 0, 100, 0, 100, 0, 100], dtype=float)
    noisy_mde = minimum_detectable_effect(
        noisy, total, labels, positive="BRCA1", alpha=0.05, delta_step=step
    )
    assert noisy_mde.power_at_zero <= 0.05
    assert noisy_mde.points > flat_mde.points


def test_attainable_alpha_is_the_largest_multiple_that_does_not_exceed_alpha():
    assert attainable_alpha(495, 0.01) == pytest.approx(4 / 495)
    assert attainable_alpha(495, 0.02) == pytest.approx(9 / 495)


def test_bh_adjust_matches_the_hand_computed_four_test_adjustment():
    # p = 0.01, 0.04, 0.03, 0.20. Rank 1 adjusts to 0.04. The rank-2 value
    # 0.06 is pulled down to the rank-3 value 0.04 * 4 / 3. The largest stays 0.20.
    adjusted = bh_adjust([0.01, 0.04, 0.03, 0.20])
    assert adjusted[0] == pytest.approx(0.04)
    assert adjusted[1] == pytest.approx(0.04 * 4 / 3)
    assert adjusted[2] == pytest.approx(0.04 * 4 / 3)
    assert adjusted[3] == pytest.approx(0.20)
