"""Tests for signals_in_the_noise.analysis.design_sensitivity."""

import math

import numpy as np
import pytest

from signals_in_the_noise.analysis.design_sensitivity import (
    DesignDecision,
    ceiling_effect,
    fold_enrichment_effect,
    minimum_detectable_effect,
    normal_errors,
    permutation_p_values,
    relabeling_masks,
    resampled_errors,
    sd_with_ci,
    simulate_ci_coverage,
    simulate_power,
    t_interval_half_width,
)
from signals_in_the_noise.analysis.statistics import exact_permutation_test

# ---------------------------------------------------------------------------
# SD estimate
# ---------------------------------------------------------------------------


def test_sd_with_ci_brackets_sample_sd():
    estimate = sd_with_ci([0.01, 0.03, -0.02, 0.0, 0.02, -0.01, 0.04, 0.01])
    assert estimate.low < estimate.sd < estimate.high
    assert estimate.n == 8


def test_sd_with_ci_matches_chi_square_bounds_for_eight_donors():
    estimate = sd_with_ci([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    assert estimate.low / estimate.sd == pytest.approx(0.661, abs=1e-3)
    assert estimate.high / estimate.sd == pytest.approx(2.035, abs=1e-3)


def test_sd_with_ci_rejects_single_value():
    with pytest.raises(ValueError, match="two values"):
        sd_with_ci([1.0])


# ---------------------------------------------------------------------------
# Vectorised permutation test
# ---------------------------------------------------------------------------


def test_relabeling_masks_enumerate_495_designs_with_observed_first():
    masks = relabeling_masks(4, 8)
    assert masks.shape == (495, 12)
    assert masks.sum(axis=1).tolist() == [4] * 495
    assert masks[0, :4].all() and not masks[0, 4:].any()


@pytest.mark.parametrize("alternative", ["two-sided", "greater", "less"])
def test_permutation_p_values_match_exact_permutation_test(alternative):
    data = np.random.default_rng(1).normal(size=(5, 12))
    vectorised = permutation_p_values(data, 4, alternative=alternative)
    expected = [
        exact_permutation_test(row[:4], row[4:], alternative=alternative).p_value for row in data
    ]
    assert vectorised == pytest.approx(expected)


def test_permutation_p_values_reject_unknown_alternative():
    with pytest.raises(ValueError, match="alternative"):
        permutation_p_values(np.zeros((1, 12)), 4, alternative="bigger")


# ---------------------------------------------------------------------------
# Power and MDE
# ---------------------------------------------------------------------------


def test_power_under_null_does_not_exceed_alpha():
    errors = normal_errors(2000, 12, seed=0)
    power = simulate_power([0.0], errors, 4, alpha=0.05)
    assert power[0] <= 0.05 + 0.015


def test_power_increases_with_effect_and_reaches_one():
    errors = normal_errors(500, 12, seed=0)
    power = simulate_power([0.0, 1.0, 2.0, 6.0], errors, 4)
    assert np.all(np.diff(power) >= 0)
    assert power[-1] == 1.0


def test_mde_for_4_vs_8_is_close_to_t_approximation():
    errors = normal_errors(4000, 12, seed=0)
    effects = np.linspace(0.0, 4.0, 81)
    mde = minimum_detectable_effect(effects, simulate_power(effects, errors, 4), target_power=0.8)
    # Noncentral-t power for a two-sided pooled t test with n = 4 and 8 is 0.8 at about 1.93 SD.
    assert mde == pytest.approx(1.93, abs=0.15)


def test_mde_scales_with_sd():
    effects = np.linspace(0.0, 8.0, 161)
    unit = minimum_detectable_effect(effects, simulate_power(effects, normal_errors(1000, 12, seed=3), 4))
    doubled = minimum_detectable_effect(
        effects, simulate_power(effects, normal_errors(1000, 12, sd=2.0, seed=3), 4)
    )
    assert doubled == pytest.approx(2 * unit, rel=1e-6)


def test_mde_interpolates_between_grid_points():
    assert minimum_detectable_effect([0.0, 1.0, 2.0], [0.05, 0.6, 1.0]) == pytest.approx(1.5)


def test_mde_is_nan_when_target_power_never_reached():
    assert math.isnan(minimum_detectable_effect([0.0, 1.0], [0.05, 0.5]))


def test_resampled_errors_have_sample_sd_and_zero_centre():
    residuals = [0.3, -0.1, 0.2, 0.0, -0.4, 0.1, 0.5, -0.2]
    errors = resampled_errors(residuals, 20000, 12, seed=0)
    assert errors.mean() == pytest.approx(0.0, abs=0.01)
    assert errors.std() == pytest.approx(np.std(residuals, ddof=1), rel=0.02)


# ---------------------------------------------------------------------------
# Intervals
# ---------------------------------------------------------------------------


def test_t_interval_half_width_for_4_vs_8():
    assert t_interval_half_width(1.0, 4, 8) == pytest.approx(2.228 * math.sqrt(3 / 8), rel=1e-3)


def test_ci_coverage_t_is_nominal_and_bootstrap_undercovers():
    coverage = simulate_ci_coverage(4, 8, n_simulations=400, n_resamples=500, seed=0)
    coverage = coverage.set_index("method")["coverage"]
    assert coverage["pooled-variance t"] == pytest.approx(0.95, abs=0.04)
    assert coverage["donor bootstrap (percentile)"] < coverage["pooled-variance t"]


# ---------------------------------------------------------------------------
# Effect benchmarks
# ---------------------------------------------------------------------------


def test_ceiling_effect_all_rescued_epithelium_lp():
    # w = 0.1, p = 0.4: Δ could reach 0.06; WT Δ averages 0.01.
    assert ceiling_effect([0.1, 0.1], [0.4, 0.4], [0.0, 0.02]) == pytest.approx(0.05)


def test_fold_enrichment_effect_caps_rescued_lp_fraction_at_one():
    assert fold_enrichment_effect([0.1, 0.2], [0.3, 0.6], fold=2.0) == pytest.approx(
        (0.1 * 0.3 + 0.2 * 0.4) / 2
    )


def test_fold_enrichment_effect_ignores_donors_without_rescued_epithelium():
    assert fold_enrichment_effect([0.0, 0.1], [float("nan"), 0.2], fold=2.0) == pytest.approx(0.01)


# ---------------------------------------------------------------------------
# Decision
# ---------------------------------------------------------------------------


def _decision(**overrides) -> DesignDecision:
    fields = dict(
        n_test_donors=4,
        n_reference_donors=8,
        delta_sd=0.01,
        delta_sd_low=0.0066,
        delta_sd_high=0.02,
        alpha=0.05,
        alternative="two-sided",
        target_power=0.8,
        mde_sd_units=1.95,
        mde=0.0195,
        mde_low=0.0129,
        mde_high=0.0397,
        mde_resampled=0.02,
        ceiling_effect=0.03,
        plausible_effect=0.01,
        plausible_fold=2.0,
        exclusion_bound=0.0136,
        identity_method="test",
    )
    return DesignDecision(**{**fields, **overrides})


def test_decision_is_estimation_when_mde_exceeds_plausible_effect():
    decision = _decision()
    assert decision.framing == "estimation"
    assert decision.ceiling_detectable
    assert "excludes effects larger than about 0.0136" in decision.summary()


def test_decision_is_hypothesis_test_when_plausible_effect_detectable():
    decision = _decision(plausible_effect=0.025)
    assert decision.framing == "hypothesis_test"
    assert decision.summary().startswith("HYPOTHESIS TEST")


def test_decision_flags_undetectable_ceiling():
    decision = _decision(ceiling_effect=0.015)
    assert not decision.ceiling_detectable
    assert "even the largest possible one" in decision.summary()


def test_decision_with_undefined_mde_is_estimation():
    assert _decision(mde=float("nan")).framing == "estimation"


def test_decision_json_round_trip(tmp_path):
    decision = _decision()
    path = tmp_path / "design.json"
    decision.to_json(path)
    assert DesignDecision.from_json(path) == decision
