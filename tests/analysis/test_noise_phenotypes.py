"""Tests for signals_in_the_noise.analysis.noise_phenotypes."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.noise_phenotypes import classify_noise_subtypes


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
