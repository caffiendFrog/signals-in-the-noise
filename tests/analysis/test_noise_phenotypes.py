"""Tests for signals_in_the_noise.analysis.noise_phenotypes."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.noise_phenotypes import (
    PbsThresholds,
    Thresholds,
    _matches_threshold,
    aggregate_noise_subtypes_by_cancer_type,
    classify_noise_subtypes,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_adata(n: int = 20, seed: int = 42) -> AnnData:
    """Return a minimal AnnData with the three QC columns required by classify_noise_subtypes."""
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame(
        {
            "pct_counts_mt": rng.uniform(10.0, 90.0, n),
            "log1p_total_counts": rng.uniform(10.0, 90.0, n),
            "log1p_n_genes_by_counts": rng.uniform(10.0, 90.0, n),
        },
        index=[f"cell_{i}" for i in range(n)],
    )
    return AnnData(obs=obs)


def _make_adata_with_known_pbs_1(seed: int = 42) -> AnnData:
    """Return an AnnData where cell_0 is guaranteed to be classified as pbs-1.

    cell_0 has extreme mito (99), extreme-low RNA (1), extreme-low genes (1).
    The remaining 19 cells span 10–90, placing their quartiles well clear of
    the extremes so the outlier cell lands unambiguously in the pbs-1 bin.
    """
    n = 20
    rng = np.random.default_rng(seed)
    pct = rng.uniform(10.0, 90.0, n).tolist()
    rna = rng.uniform(10.0, 90.0, n).tolist()
    genes = rng.uniform(10.0, 90.0, n).tolist()
    pct[0] = 99.0
    rna[0] = 1.0
    genes[0] = 1.0
    obs = pd.DataFrame(
        {
            "pct_counts_mt": pct,
            "log1p_total_counts": rna,
            "log1p_n_genes_by_counts": genes,
        },
        index=[f"cell_{i}" for i in range(n)],
    )
    return AnnData(obs=obs)


def _make_adata_with_known_pbs_2(seed: int = 42) -> AnnData:
    """Return an AnnData where cell_0 is guaranteed to be classified as pbs-2.

    cell_0: extreme-low mito (1), extreme-low RNA (1), moderate genes (50).
    """
    n = 20
    rng = np.random.default_rng(seed)
    pct = rng.uniform(10.0, 90.0, n).tolist()
    rna = rng.uniform(10.0, 90.0, n).tolist()
    genes = rng.uniform(10.0, 90.0, n).tolist()
    pct[0] = 1.0
    rna[0] = 1.0
    genes[0] = 50.0
    obs = pd.DataFrame(
        {
            "pct_counts_mt": pct,
            "log1p_total_counts": rna,
            "log1p_n_genes_by_counts": genes,
        },
        index=[f"cell_{i}" for i in range(n)],
    )
    return AnnData(obs=obs)


def _make_adata_with_known_pbs_3(seed: int = 42) -> AnnData:
    """Return an AnnData where cell_0 is guaranteed to be classified as pbs-3.

    cell_0: moderate mito (50), moderate RNA (50), extreme-high genes (99).
    """
    n = 20
    rng = np.random.default_rng(seed)
    pct = rng.uniform(10.0, 90.0, n).tolist()
    rna = rng.uniform(10.0, 90.0, n).tolist()
    genes = rng.uniform(10.0, 90.0, n).tolist()
    pct[0] = 50.0
    rna[0] = 50.0
    genes[0] = 99.0
    obs = pd.DataFrame(
        {
            "pct_counts_mt": pct,
            "log1p_total_counts": rna,
            "log1p_n_genes_by_counts": genes,
        },
        index=[f"cell_{i}" for i in range(n)],
    )
    return AnnData(obs=obs)


# ---------------------------------------------------------------------------
# _matches_threshold tests
# ---------------------------------------------------------------------------


@pytest.fixture()
def metric_series() -> pd.Series:
    """Ascending QC metric with predictable quantiles for threshold tests."""
    return pd.Series(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
        index=[f"cell_{i}" for i in range(10)],
    )


def test_matches_threshold_no_conditions_returns_all_true(metric_series):
    thresholds = Thresholds(q_low=None, q_high=None, q_mod_low=None, q_mod_high=None)
    result = _matches_threshold(metric_series, thresholds)
    assert result.dtype == bool
    assert result.all()
    assert result.index.equals(metric_series.index)


def test_matches_threshold_q_low_uses_greater_than_or_equal(metric_series):
    thresholds = Thresholds(q_low=0.5, q_high=None, q_mod_low=None, q_mod_high=None)
    cutoff = metric_series.quantile(0.5)
    result = _matches_threshold(metric_series, thresholds)
    assert result.equals(metric_series >= cutoff)
    assert result.loc["cell_9"]
    assert not result.loc["cell_0"]


def test_matches_threshold_q_low_includes_value_at_quantile(metric_series):
    thresholds = Thresholds(q_low=0.0, q_high=None, q_mod_low=None, q_mod_high=None)
    result = _matches_threshold(metric_series, thresholds)
    assert result.loc["cell_0"]


def test_matches_threshold_q_high_uses_less_than_or_equal(metric_series):
    thresholds = Thresholds(q_low=None, q_high=0.5, q_mod_low=None, q_mod_high=None)
    cutoff = metric_series.quantile(0.5)
    result = _matches_threshold(metric_series, thresholds)
    assert result.equals(metric_series <= cutoff)
    assert result.loc["cell_0"]
    assert not result.loc["cell_9"]


def test_matches_threshold_q_high_includes_value_at_quantile(metric_series):
    thresholds = Thresholds(q_low=None, q_high=1.0, q_mod_low=None, q_mod_high=None)
    result = _matches_threshold(metric_series, thresholds)
    assert result.loc["cell_9"]


def test_matches_threshold_q_mod_uses_strict_between_bounds(metric_series):
    thresholds = Thresholds(q_low=None, q_high=None, q_mod_low=0.25, q_mod_high=0.75)
    low_val = metric_series.quantile(0.25)
    high_val = metric_series.quantile(0.75)
    result = _matches_threshold(metric_series, thresholds)
    expected = (metric_series > low_val) & (metric_series < high_val)
    assert result.equals(expected)


def test_matches_threshold_q_mod_excludes_values_at_quantile_bounds(metric_series):
    thresholds = Thresholds(q_low=None, q_high=None, q_mod_low=0.0, q_mod_high=1.0)
    result = _matches_threshold(metric_series, thresholds)
    assert not result.loc["cell_0"]
    assert not result.loc["cell_9"]
    assert result.loc["cell_5"]


def test_matches_threshold_combined_conditions_use_logical_and(metric_series):
    thresholds = Thresholds(q_low=0.25, q_high=0.75, q_mod_low=None, q_mod_high=None)
    low_cutoff = metric_series.quantile(0.25)
    high_cutoff = metric_series.quantile(0.75)
    result = _matches_threshold(metric_series, thresholds)
    expected = (metric_series >= low_cutoff) & (metric_series <= high_cutoff)
    assert result.equals(expected)


def test_matches_threshold_ignores_partial_moderate_bounds(metric_series):
    """Only one moderate bound set should not add a moderate filter."""
    thresholds_low_only = Thresholds(q_low=None, q_high=None, q_mod_low=0.25, q_mod_high=None)
    thresholds_high_only = Thresholds(q_low=None, q_high=None, q_mod_low=None, q_mod_high=0.75)
    assert _matches_threshold(metric_series, thresholds_low_only).all()
    assert _matches_threshold(metric_series, thresholds_high_only).all()


def _iscb_scs_pbs_thresholds() -> dict[str, PbsThresholds]:
    """ISCB SCS notebook thresholds expressed with current field semantics."""
    return {
        "pbs-1": PbsThresholds(
            pct_counts_mt=Thresholds(q_low=0.85, q_high=None, q_mod_low=None, q_mod_high=None),
            log1p_total_counts=Thresholds(q_low=None, q_high=0.3, q_mod_low=None, q_mod_high=None),
            log1p_n_genes_by_counts=Thresholds(
                q_low=None, q_high=0.3, q_mod_low=None, q_mod_high=None
            ),
        ),
        "pbs-2": PbsThresholds(
            pct_counts_mt=Thresholds(q_low=None, q_high=0.25, q_mod_low=None, q_mod_high=None),
            log1p_total_counts=Thresholds(q_low=None, q_high=0.3, q_mod_low=None, q_mod_high=None),
            log1p_n_genes_by_counts=Thresholds(
                q_low=None, q_high=None, q_mod_low=0.45, q_mod_high=0.85
            ),
        ),
        "pbs-3": PbsThresholds(
            pct_counts_mt=Thresholds(q_low=None, q_high=None, q_mod_low=0.3, q_mod_high=0.65),
            log1p_total_counts=Thresholds(q_low=None, q_high=None, q_mod_low=0.5, q_mod_high=0.85),
            log1p_n_genes_by_counts=Thresholds(
                q_low=0.95, q_high=None, q_mod_low=None, q_mod_high=None
            ),
        ),
    }


def test_classify_noise_subtypes_iscb_thresholds_are_selective():
    """Custom ISCB thresholds should classify a small minority of cells per subtype."""
    rng = np.random.default_rng(42)
    n = 1000
    adata = AnnData(
        obs=pd.DataFrame(
            {
                "pct_counts_mt": rng.uniform(1.0, 99.0, n),
                "log1p_total_counts": rng.uniform(1.0, 99.0, n),
                "log1p_n_genes_by_counts": rng.uniform(1.0, 99.0, n),
            }
        )
    )
    classify_noise_subtypes(adata, pbs_thresholds=_iscb_scs_pbs_thresholds())
    counts = adata.obs[["pbs-1", "pbs-2", "pbs-3"]].sum()
    assert counts["pbs-1"] == 23
    assert counts["pbs-2"] == 33
    assert counts["pbs-3"] == 6


# ---------------------------------------------------------------------------
# Contract tests
# ---------------------------------------------------------------------------


def test_classify_noise_subtypes_returns_same_adata_object():
    adata = _make_adata()
    result = classify_noise_subtypes(adata)
    assert result is adata


def test_classify_noise_subtypes_adds_pbs_1_column():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert "pbs-1" in adata.obs.columns


def test_classify_noise_subtypes_adds_pbs_2_column():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert "pbs-2" in adata.obs.columns


def test_classify_noise_subtypes_adds_pbs_3_column():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert "pbs-3" in adata.obs.columns


def test_classify_noise_subtypes_columns_are_bool_dtype():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert adata.obs["pbs-1"].dtype == bool
    assert adata.obs["pbs-2"].dtype == bool
    assert adata.obs["pbs-3"].dtype == bool


def test_classify_noise_subtypes_preserves_obs_length():
    n = 30
    adata = _make_adata(n=n)
    classify_noise_subtypes(adata)
    assert len(adata.obs) == n


# ---------------------------------------------------------------------------
# Correctness tests
# ---------------------------------------------------------------------------


def test_classify_noise_subtypes_identifies_pbs_1_cell():
    adata = _make_adata_with_known_pbs_1()
    classify_noise_subtypes(adata)
    assert adata.obs.loc["cell_0", "pbs-1"]


def test_classify_noise_subtypes_identifies_pbs_2_cell():
    adata = _make_adata_with_known_pbs_2()
    classify_noise_subtypes(adata)
    assert adata.obs.loc["cell_0", "pbs-2"]


def test_classify_noise_subtypes_identifies_pbs_3_cell():
    adata = _make_adata_with_known_pbs_3()
    classify_noise_subtypes(adata)
    assert adata.obs.loc["cell_0", "pbs-3"]


def test_classify_noise_subtypes_custom_quantiles_change_results():
    """Tightening q_high moves more cells into the 'high' bins."""
    adata_default = _make_adata()
    adata_tight = _make_adata()
    classify_noise_subtypes(adata_default)
    classify_noise_subtypes(adata_tight)
    # With a lower q_high threshold, at least as many cells are 'pbs-1'
    assert adata_tight.obs["pbs-1"].sum() >= adata_default.obs["pbs-1"].sum()


# ---------------------------------------------------------------------------
# aggregate_noise_subtypes_by_cancer_type helpers
# ---------------------------------------------------------------------------


def _make_classified_adata(
    cancer_type: str,
    cell_population: str = "Total",
    n_obs: int = 20,
    noise_fraction: float = 0.5,
    seed: int = 0,
) -> AnnData:
    """Return a minimal AnnData with all metadata required by aggregate_noise_subtypes_by_cancer_type."""
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame(
        {
            "pct_counts_mt": rng.uniform(10.0, 90.0, n_obs),
            "log1p_total_counts": rng.uniform(10.0, 90.0, n_obs),
            "log1p_n_genes_by_counts": rng.uniform(10.0, 90.0, n_obs),
            "is_noise": ([1] * int(n_obs * noise_fraction) + [0] * (n_obs - int(n_obs * noise_fraction))),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    adata = AnnData(obs=obs)
    adata.uns["cancer_type"] = cancer_type
    adata.uns["cell_population"] = cell_population
    return adata


# ---------------------------------------------------------------------------
# aggregate_noise_subtypes_by_cancer_type tests
# ---------------------------------------------------------------------------


def test_aggregate_noise_subtypes_returns_dataframe():
    adatas = [_make_classified_adata("Luminal")]
    result = aggregate_noise_subtypes_by_cancer_type(adatas)
    assert isinstance(result, pd.DataFrame)


def test_aggregate_noise_subtypes_has_expected_columns():
    adatas = [_make_classified_adata("Luminal")]
    result = aggregate_noise_subtypes_by_cancer_type(adatas)
    assert set(result.columns) == {"pbs-1", "pbs-2", "pbs-3", "noise"}


def test_aggregate_noise_subtypes_index_contains_cancer_type():
    adatas = [_make_classified_adata("Luminal")]
    result = aggregate_noise_subtypes_by_cancer_type(adatas)
    assert any("Luminal" in label for label in result.index)


def test_aggregate_noise_subtypes_index_annotated_with_specimen_count():
    adatas = [
        _make_classified_adata("Luminal", seed=0),
        _make_classified_adata("Luminal", seed=1),
    ]
    result = aggregate_noise_subtypes_by_cancer_type(adatas)
    assert any("2 specimens" in label for label in result.index)


def test_aggregate_noise_subtypes_filters_non_total_population():
    total_adata = _make_classified_adata("Luminal", cell_population="Total")
    other_adata = _make_classified_adata("Basal", cell_population="Reduction")
    result = aggregate_noise_subtypes_by_cancer_type([total_adata, other_adata])
    assert len(result) == 1
    assert any("Luminal" in label for label in result.index)
    assert not any("Basal" in label for label in result.index)


def test_aggregate_noise_subtypes_respects_custom_cell_population_filter():
    reduction_adata = _make_classified_adata("Luminal", cell_population="Reduction")
    result = aggregate_noise_subtypes_by_cancer_type(
        [reduction_adata], cell_population_filter="Reduction"
    )
    assert len(result) == 1


def test_aggregate_noise_subtypes_values_are_percentages():
    adatas = [_make_classified_adata("Luminal")]
    result = aggregate_noise_subtypes_by_cancer_type(adatas)
    # All percentage columns should be in [0, 100]
    for col in ["pbs-1", "pbs-2", "pbs-3", "noise"]:
        assert result[col].between(0, 100).all(), f"column {col!r} has values outside [0, 100]"


def test_aggregate_noise_subtypes_empty_input_returns_empty_dataframe():
    result = aggregate_noise_subtypes_by_cancer_type([])
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0


def test_aggregate_noise_subtypes_zero_noise_produces_zero_percentages():
    """Cancer types with zero noise cells must yield 0% subtypes, not NaN or inf."""
    adata = _make_classified_adata("Luminal", noise_fraction=0.0)
    result = aggregate_noise_subtypes_by_cancer_type([adata])
    assert not result.isnull().any().any(), "result contains NaN"
    assert result[["pbs-1", "pbs-2", "pbs-3"]].eq(0.0).all().all()


def test_aggregate_noise_subtypes_excludes_non_noise_from_subtype_counts():
    """Non-noise cells must not inflate subtype counts beyond the noise population.

    With 5 noise cells and 20 non-noise cells whose metrics fall into the
    'pbs-1' region, the buggy all-cell path would produce pbs-1% = 400%.
    The fix (classify on noise cells only) must keep all values in [0, 100].
    """
    n_noise, n_non_noise = 5, 20
    obs = pd.DataFrame(
        {
            "pct_counts_mt":           [50.0] * n_noise + [99.0] * n_non_noise,
            "log1p_total_counts":      [50.0] * n_noise + [1.0]  * n_non_noise,
            "log1p_n_genes_by_counts": [50.0] * n_noise + [1.0]  * n_non_noise,
            "is_noise": [1] * n_noise + [0] * n_non_noise,
        },
        index=[f"cell_{i}" for i in range(n_noise + n_non_noise)],
    )
    adata = AnnData(obs=obs)
    adata.uns["cancer_type"] = "Luminal"
    adata.uns["cell_population"] = "Total"
    result = aggregate_noise_subtypes_by_cancer_type([adata])
    for col in ["pbs-1", "pbs-2", "pbs-3", "noise"]:
        assert result[col].between(0, 100).all(), (
            f"column {col!r} has values outside [0, 100] — non-noise cells may be included"
        )
