"""Tests for signals_in_the_noise.utils.visualization."""

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from signals_in_the_noise.utils.visualization import (
    plot_empty_cell_violin_comparison,
    plot_gene_signature_score_distribution,
    plot_gsea_nes_heatmap,
    plot_noise_subtype_comparison,
    plot_pathway_heatmap,
    plot_score_heatmap,
)


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


def test_plot_pathway_heatmap_dense_matrix():
    """adata.X as a dense NumPy array must not raise AttributeError."""
    adata = _make_adata()
    adata.X = np.asarray(adata.X.toarray())
    genes = [f"GENE_{i}" for i in range(6)]
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes)
    assert isinstance(result, list)
    assert all(isinstance(ax, matplotlib.axes.Axes) for ax in result)
    plt.close("all")


def test_plot_pathway_heatmap_all_genes_absent():
    """When no requested genes exist in adata the function must not raise."""
    adata = _make_adata()
    genes = ["NOT_A", "NOT_B", "NOT_C"]
    result = plot_pathway_heatmap(adata, "TEST_PATHWAY", genes)
    assert isinstance(result, list)
    plt.close("all")


def test_plot_pathway_heatmap_raises_on_axes_length_mismatch():
    """Providing the wrong number of axes must raise ValueError immediately."""
    adata = _make_adata()
    genes = [f"GENE_{i}" for i in range(6)]
    fig, wrong_ax = plt.subplots(1, 1)  # 1 axis, but n_panels defaults to 2
    with pytest.raises(ValueError, match="n_panels"):
        plot_pathway_heatmap(adata, "TEST_PATHWAY", genes, axes=[wrong_ax])
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


# ---------------------------------------------------------------------------
# plot_gsea_nes_heatmap helpers
# ---------------------------------------------------------------------------


def _make_gsea_df() -> pd.DataFrame:
    """Return a minimal GSEA results DataFrame."""
    return pd.DataFrame(
        {
            "Term": ["HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_DNA_REPAIR"],
            "NES": [1.8, -1.2, 0.5],
            "FDR q-val": [0.005, 0.04, 0.3],
        }
    )


# ---------------------------------------------------------------------------
# plot_gsea_nes_heatmap tests
# ---------------------------------------------------------------------------


def test_plot_gsea_nes_heatmap_returns_axes():
    df = _make_gsea_df()
    result = plot_gsea_nes_heatmap(df)
    assert isinstance(result, matplotlib.axes.Axes)
    plt.close("all")


def test_plot_gsea_nes_heatmap_uses_provided_axes():
    df = _make_gsea_df()
    fig, ax_in = plt.subplots()
    ax_out = plot_gsea_nes_heatmap(df, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


def test_plot_gsea_nes_heatmap_does_not_raise_on_valid_input():
    df = _make_gsea_df()
    plot_gsea_nes_heatmap(df)
    plt.close("all")


def test_plot_gsea_nes_heatmap_ytick_labels_contain_stars_for_significant_terms():
    """Significant terms (FDR < 0.1) should have star suffixes on the y-axis."""
    df = _make_gsea_df()
    fig, ax = plt.subplots()
    plot_gsea_nes_heatmap(df, ax=ax)
    tick_texts = [t.get_text() for t in ax.get_yticklabels()]
    # HALLMARK_G2M_CHECKPOINT has FDR 0.005 → " ***"
    assert any("***" in t for t in tick_texts)
    # HALLMARK_DNA_REPAIR has FDR 0.3 → no stars
    assert any("DNA_REPAIR" in t and "*" not in t for t in tick_texts)
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_noise_subtype_comparison helpers
# ---------------------------------------------------------------------------


def _make_norm_df() -> pd.DataFrame:
    """Return a minimal normalised noise-subtype DataFrame."""
    import pandas as pd
    return pd.DataFrame(
        {
            "damaged": [15.0, 20.0],
            "dormant": [10.0, 8.0],
            "multifunction": [25.0, 18.0],
            "noise": [30.0, 35.0],
        },
        index=["Luminal (2 specimens)", "Basal (3 specimens)"],
    )


# ---------------------------------------------------------------------------
# plot_noise_subtype_comparison tests
# ---------------------------------------------------------------------------


def test_plot_noise_subtype_comparison_returns_two_axes():
    df = _make_norm_df()
    result = plot_noise_subtype_comparison(df)
    assert len(result) == 2
    assert all(isinstance(ax, matplotlib.axes.Axes) for ax in result)
    plt.close("all")


def test_plot_noise_subtype_comparison_uses_provided_axes():
    df = _make_norm_df()
    fig, (ax1_in, ax2_in) = plt.subplots(1, 2)
    ax1_out, ax2_out = plot_noise_subtype_comparison(df, axes=(ax1_in, ax2_in))
    assert ax1_out is ax1_in
    assert ax2_out is ax2_in
    plt.close("all")


def test_plot_noise_subtype_comparison_does_not_raise_on_valid_input():
    df = _make_norm_df()
    plot_noise_subtype_comparison(df)
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_empty_cell_violin_comparison helpers
# ---------------------------------------------------------------------------


def _make_violin_adata(n_obs: int = 12, seed: int = 0) -> AnnData:
    """Return a minimal AnnData suitable for sc.pl.violin with zero_mito groupby specimen_id."""
    import scipy.sparse as sp_scipy
    rng = np.random.default_rng(seed)
    X = sp_scipy.csr_matrix(rng.integers(0, 5, size=(n_obs, 5)).astype(float))
    obs = pd.DataFrame(
        {
            "zero_mito": rng.integers(0, 2, n_obs).astype(float),
            "specimen_id": ["specA", "specB", "specC"] * (n_obs // 3),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    adata = AnnData(X=X, obs=obs)
    adata.var_names = [f"GENE_{i}" for i in range(5)]
    return adata


# ---------------------------------------------------------------------------
# plot_empty_cell_violin_comparison tests
# ---------------------------------------------------------------------------


def test_plot_empty_cell_violin_comparison_returns_list_of_axes():
    noise = _make_violin_adata(seed=0)
    combined = _make_violin_adata(seed=1)
    result = plot_empty_cell_violin_comparison(noise, combined, groupby="specimen_id")
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(isinstance(ax, matplotlib.axes.Axes) for ax in result)
    plt.close("all")


def test_plot_empty_cell_violin_comparison_uses_provided_axes():
    noise = _make_violin_adata(seed=0)
    combined = _make_violin_adata(seed=1)
    fig, (ax1_in, ax2_in) = plt.subplots(1, 2)
    result = plot_empty_cell_violin_comparison(
        noise, combined, groupby="specimen_id", axes=[ax1_in, ax2_in]
    )
    assert result[0] is ax1_in
    assert result[1] is ax2_in
    plt.close("all")


def test_plot_empty_cell_violin_comparison_does_not_raise_on_valid_input():
    noise = _make_violin_adata(seed=0)
    combined = _make_violin_adata(seed=1)
    plot_empty_cell_violin_comparison(noise, combined, groupby="specimen_id")
    plt.close("all")


def test_plot_empty_cell_violin_comparison_raises_on_axes_length_mismatch():
    """Passing the wrong number of axes must raise ValueError immediately."""
    noise = _make_violin_adata(seed=0)
    combined = _make_violin_adata(seed=1)
    fig, ax = plt.subplots(1, 1)  # 1 axis instead of required 2
    with pytest.raises(ValueError, match="2 axes"):
        plot_empty_cell_violin_comparison(noise, combined, groupby="specimen_id", axes=[ax])
    plt.close("all")


# ---------------------------------------------------------------------------
# plot_gene_signature_score_distribution tests
# ---------------------------------------------------------------------------


def test_plot_gene_signature_score_distribution_returns_axes():
    scores = pd.Series([0.1, -0.3, 0.5, 0.2, -0.1, 0.4])
    result = plot_gene_signature_score_distribution(scores, "basal")
    assert isinstance(result, matplotlib.axes.Axes)
    plt.close("all")


def test_plot_gene_signature_score_distribution_uses_provided_axes():
    scores = pd.Series([0.1, -0.3, 0.5, 0.2, -0.1, 0.4])
    fig, ax_in = plt.subplots()
    ax_out = plot_gene_signature_score_distribution(scores, "basal", ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


def test_plot_gene_signature_score_distribution_xlabel_contains_signature_name():
    scores = pd.Series([0.1, 0.2, 0.3])
    fig, ax = plt.subplots()
    plot_gene_signature_score_distribution(scores, "luminal_progenitor", ax=ax)
    assert "luminal_progenitor" in ax.get_xlabel()
    plt.close("all")


def test_plot_gene_signature_score_distribution_does_not_raise_on_valid_input():
    scores = pd.Series([0.1, -0.2, 0.3, -0.4, 0.5])
    plot_gene_signature_score_distribution(scores, "basal")
    plt.close("all")
