"""Tests for signals_in_the_noise.analysis.noise_phenotypes."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.noise_phenotypes import (
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


def _make_adata_with_known_damaged(seed: int = 42) -> AnnData:
    """Return an AnnData where cell_0 is guaranteed to be classified as damaged.

    cell_0 has extreme mito (99), extreme-low RNA (1), extreme-low genes (1).
    The remaining 19 cells span 10–90, placing their quartiles well clear of
    the extremes so the outlier cell lands unambiguously in the damaged bin.
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


def _make_adata_with_known_dormant(seed: int = 42) -> AnnData:
    """Return an AnnData where cell_0 is guaranteed to be classified as dormant.

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


def _make_adata_with_known_multifunction(seed: int = 42) -> AnnData:
    """Return an AnnData where cell_0 is guaranteed to be classified as multifunction.

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
# Contract tests
# ---------------------------------------------------------------------------


def test_classify_noise_subtypes_returns_same_adata_object():
    adata = _make_adata()
    result = classify_noise_subtypes(adata)
    assert result is adata


def test_classify_noise_subtypes_adds_damaged_column():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert "damaged" in adata.obs.columns


def test_classify_noise_subtypes_adds_dormant_column():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert "dormant" in adata.obs.columns


def test_classify_noise_subtypes_adds_multifunction_column():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert "multifunction" in adata.obs.columns


def test_classify_noise_subtypes_columns_are_bool_dtype():
    adata = _make_adata()
    classify_noise_subtypes(adata)
    assert adata.obs["damaged"].dtype == bool
    assert adata.obs["dormant"].dtype == bool
    assert adata.obs["multifunction"].dtype == bool


def test_classify_noise_subtypes_preserves_obs_length():
    n = 30
    adata = _make_adata(n=n)
    classify_noise_subtypes(adata)
    assert len(adata.obs) == n


# ---------------------------------------------------------------------------
# Correctness tests
# ---------------------------------------------------------------------------


def test_classify_noise_subtypes_identifies_damaged_cell():
    adata = _make_adata_with_known_damaged()
    classify_noise_subtypes(adata)
    assert adata.obs.loc["cell_0", "damaged"]


def test_classify_noise_subtypes_identifies_dormant_cell():
    adata = _make_adata_with_known_dormant()
    classify_noise_subtypes(adata)
    assert adata.obs.loc["cell_0", "dormant"]


def test_classify_noise_subtypes_identifies_multifunction_cell():
    adata = _make_adata_with_known_multifunction()
    classify_noise_subtypes(adata)
    assert adata.obs.loc["cell_0", "multifunction"]


def test_classify_noise_subtypes_custom_quantiles_change_results():
    """Tightening q_high moves more cells into the 'high' bins."""
    adata_default = _make_adata()
    adata_tight = _make_adata()
    classify_noise_subtypes(adata_default, q_low=0.25, q_high=0.75)
    classify_noise_subtypes(adata_tight, q_low=0.25, q_high=0.50)
    # With a lower q_high threshold, at least as many cells are 'damaged'
    assert adata_tight.obs["damaged"].sum() >= adata_default.obs["damaged"].sum()


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
    assert set(result.columns) == {"damaged", "dormant", "multifunction", "noise"}


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
    for col in ["damaged", "dormant", "multifunction", "noise"]:
        assert result[col].between(0, 100).all(), f"column {col!r} has values outside [0, 100]"


def test_aggregate_noise_subtypes_empty_input_returns_empty_dataframe():
    result = aggregate_noise_subtypes_by_cancer_type([])
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0


def test_aggregate_noise_subtypes_excludes_non_noise_from_subtype_counts():
    """Non-noise cells must not inflate subtype counts beyond the noise population.

    With 5 noise cells and 20 non-noise cells whose metrics fall into the
    'damaged' region, the buggy all-cell path would produce damaged% = 400%.
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
    for col in ["damaged", "dormant", "multifunction", "noise"]:
        assert result[col].between(0, 100).all(), (
            f"column {col!r} has values outside [0, 100] — non-noise cells may be included"
        )
