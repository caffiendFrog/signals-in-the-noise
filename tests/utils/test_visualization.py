"""Tests for signals_in_the_noise.utils.visualization."""

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from signals_in_the_noise.utils.visualization import plot_pathway_heatmap, plot_score_heatmap


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_adata(n_obs: int = 6, n_vars: int = 10, seed: int = 0) -> AnnData:
    """Return a minimal AnnData with a sparse count matrix and specimen_id obs."""
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.integers(0, 10, size=(n_obs, n_vars)).astype(float))
    adata = AnnData(X)
    adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.obs["specimen_id"] = ["specA", "specA", "specB", "specB", "specC", "specC"]
    return adata


# ---------------------------------------------------------------------------
# Contract tests
# ---------------------------------------------------------------------------


def test_plot_pathway_heatmap_returns_list_of_axes():
    adata = _make_adata()
    genes = [f"GENE_{i}" for i in range(6)]
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes)
    assert isinstance(result, list)
    assert all(isinstance(ax, matplotlib.axes.Axes) for ax in result)
    plt.close("all")


def test_plot_pathway_heatmap_returns_n_panels_axes():
    adata = _make_adata()
    genes = [f"GENE_{i}" for i in range(6)]
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes, n_panels=2)
    assert len(result) == 2
    plt.close("all")


def test_plot_pathway_heatmap_uses_provided_axes():
    adata = _make_adata()
    genes = [f"GENE_{i}" for i in range(4)]
    fig, provided = plt.subplots(1, 2)
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes, n_panels=2, axes=list(provided))
    assert result[0] is provided[0]
    assert result[1] is provided[1]
    plt.close("all")


def test_plot_pathway_heatmap_drops_absent_genes():
    """Genes not in adata.var_names should be silently excluded."""
    adata = _make_adata()
    genes = ["GENE_0", "GENE_1", "NOT_IN_DATASET"]
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes)
    assert result is not None
    plt.close("all")


def test_plot_pathway_heatmap_does_not_raise_with_single_panel():
    adata = _make_adata()
    genes = [f"GENE_{i}" for i in range(4)]
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes, n_panels=1)
    assert len(result) == 1
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_score_heatmap helpers
# ---------------------------------------------------------------------------


def _make_score_adata(n_obs: int = 9, seed: int = 0) -> AnnData:
    """Return a minimal AnnData with tp_score, cr_score, and specimen_id obs."""
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame(
        {
            "tp_score": rng.uniform(-1.0, 1.0, n_obs),
            "cr_score": rng.uniform(-1.0, 1.0, n_obs),
            "specimen_id": ["specA", "specB", "specC"] * (n_obs // 3),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    return AnnData(obs=obs)


# ---------------------------------------------------------------------------
# plot_score_heatmap tests
# ---------------------------------------------------------------------------


def test_plot_score_heatmap_returns_axes():
    adata = _make_score_adata()
    result = plot_score_heatmap(adata, score_columns=["tp_score", "cr_score"])
    assert isinstance(result, matplotlib.axes.Axes)
    plt.close("all")


def test_plot_score_heatmap_uses_provided_axes():
    adata = _make_score_adata()
    fig, ax_in = plt.subplots()
    ax_out = plot_score_heatmap(adata, score_columns=["tp_score", "cr_score"], ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


def test_plot_score_heatmap_accepts_single_score_column():
    adata = _make_score_adata()
    result = plot_score_heatmap(adata, score_columns=["tp_score"])
    assert result is not None
    plt.close("all")


def test_plot_score_heatmap_does_not_raise_on_valid_input():
    adata = _make_score_adata()
    plot_score_heatmap(adata, score_columns=["tp_score", "cr_score"])
    plt.close("all")
