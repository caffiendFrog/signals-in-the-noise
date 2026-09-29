"""Tests for signals_in_the_noise.preprocessing.qc."""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import median_abs_deviation

from signals_in_the_noise.preprocessing.qc import (
    NOISE_FLAG_COLUMNS,
    QC_METRICS,
    QcThresholds,
    annotate_qc_metrics,
    depth_matched_qc_obs,
    exclusive_reason,
    flag_noise,
    library_depth,
    mad_thresholds,
    modal_thresholds,
    relabel_noise,
    select_population,
    thin_counts,
)

THRESHOLDS = QcThresholds(mito_upper=0.2, genes_lower=200, genes_upper=5000, total_upper=30000)


def _obs(genes, mito_pct, totals) -> pd.DataFrame:
    return pd.DataFrame(
        {"n_genes_by_counts": genes, "pct_counts_mt": mito_pct, "total_counts": totals},
        index=[f"cell_{i}" for i in range(len(genes))],
    )


def _count_adata(n_obs: int = 200, n_vars: int = 60, lam: float = 5.0, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(lam, size=(n_obs, n_vars)).astype(np.float32))
    adata = AnnData(X)
    adata.var_names = ["MT-CO1", "MT-ND1"] + [f"GENE_{i}" for i in range(n_vars - 2)]
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


# ---------------------------------------------------------------------------
# QcThresholds / flag_noise
# ---------------------------------------------------------------------------


def test_from_uns_reads_published_keys_as_floats():
    uns = {
        "qc_mito_upper": 0.2,
        "qc_genes_lower": 500,
        "qc_genes_upper": 6000,
        "qc_total_upper": 40000,
    }
    thresholds = QcThresholds.from_uns(uns)
    assert thresholds == QcThresholds(
        mito_upper=0.2, genes_lower=500.0, genes_upper=6000.0, total_upper=40000.0
    )
    assert all(isinstance(v, float) for v in vars(thresholds).values())


def test_flag_noise_returns_int_flags_indexed_like_obs():
    obs = _obs([300], [5.0], [1000])
    flags = flag_noise(obs, THRESHOLDS)
    assert list(flags.columns) == [*NOISE_FLAG_COLUMNS, "is_noise"]
    assert flags.index.equals(obs.index)
    assert all(dtype.kind == "i" for dtype in flags.dtypes)


@pytest.mark.parametrize(
    ("genes", "mito_pct", "totals", "flag"),
    [
        (200, 5.0, 1000, "is_low_num_genes"),  # genes <= lower
        (5001, 5.0, 1000, "is_high_num_genes"),  # genes > upper
        (300, 20.1, 1000, "is_high_mito"),  # fraction > upper
        (300, 5.0, 30000, "is_high_total_count"),  # total >= upper
    ],
)
def test_flag_noise_sets_expected_flag_at_boundary(genes, mito_pct, totals, flag):
    flags = flag_noise(_obs([genes], [mito_pct], [totals]), THRESHOLDS)
    assert flags.loc["cell_0", flag] == 1
    assert flags.loc["cell_0", "is_noise"] == 1
    others = [c for c in NOISE_FLAG_COLUMNS if c != flag]
    assert flags.loc["cell_0", others].sum() == 0


def test_flag_noise_does_not_flag_values_just_inside_bounds():
    flags = flag_noise(_obs([201, 5000], [20.0, 20.0], [29999, 29999]), THRESHOLDS)
    assert flags["is_noise"].sum() == 0


def test_gse161529_apply_one_keeps_pal_comparators():
    """Preprocessing must flag cells exactly as Pal et al.'s rule, written out independently."""
    from signals_in_the_noise.preprocessing.gse161529 import GSE161529

    adata = _count_adata(n_obs=50, n_vars=500, lam=0.5)
    adata.obs["adata-filename"] = GSE161529.EXPECTED_MISMATCHES[0]
    uns = {"qc_genes_lower": 190, "qc_genes_upper": 205, "qc_mito_upper": 0.006}
    adata.uns.update({**uns, "qc_total_upper": 265, "num_cells_after": 0})
    GSE161529._apply_one(adata)

    obs = adata.obs
    expected = pd.DataFrame(
        {
            "is_low_num_genes": obs["n_genes_by_counts"] <= 190,
            "is_high_num_genes": obs["n_genes_by_counts"] > 205,
            "is_high_mito": obs["pct_counts_mt"] / 100 > 0.006,
            "is_high_total_count": obs["total_counts"] >= 265,
        }
    )
    expected["is_noise"] = expected.any(axis=1)
    pd.testing.assert_frame_equal(obs[expected.columns], expected.astype(int))
    for column in expected.columns:
        assert 0 < obs[column].sum() < adata.n_obs, f"{column} does not discriminate"


# ---------------------------------------------------------------------------
# modal_thresholds / mad_thresholds / relabel_noise
# ---------------------------------------------------------------------------


def test_exclusive_reason_uses_fixed_priority_and_leaves_clean_cells_missing():
    flags = pd.DataFrame(
        {
            "is_high_mito": [1, 0, 0, 0],
            "is_low_num_genes": [1, 1, 0, 0],
            "is_high_num_genes": [0, 0, 1, 0],
            "is_high_total_count": [1, 0, 0, 0],
        }
    )
    reason = exclusive_reason(flags)
    assert reason.tolist()[:3] == ["high_mito", "low_genes", "high_genes"]
    assert pd.isna(reason.iloc[3])


def test_exclusive_reason_rejects_a_missing_flag_column():
    with pytest.raises(ValueError, match="missing"):
        exclusive_reason(pd.DataFrame({"is_high_mito": [1]}))


def test_modal_thresholds_uses_mode_of_each_field():
    thresholds = [
        QcThresholds(0.2, 500, 6000, 40000),
        QcThresholds(0.2, 500, 7000, 40000),
        QcThresholds(0.3, 300, 7000, 50000),
    ]
    assert modal_thresholds(thresholds) == QcThresholds(0.2, 500, 7000, 40000)


def test_modal_thresholds_breaks_ties_with_smallest_value():
    thresholds = [QcThresholds(0.2, 500, 6000, 40000), QcThresholds(0.3, 300, 7000, 50000)]
    assert modal_thresholds(thresholds) == QcThresholds(0.2, 300, 6000, 40000)


def test_modal_thresholds_rejects_empty_input():
    with pytest.raises(ValueError):
        modal_thresholds([])


def test_mad_thresholds_matches_manual_computation():
    rng = np.random.default_rng(1)
    obs = pd.DataFrame(
        {
            "pct_counts_mt": rng.uniform(1, 30, 500),
            "log1p_total_counts": rng.normal(8, 0.5, 500),
            "log1p_n_genes_by_counts": rng.normal(7, 0.4, 500),
        }
    )

    def bound(values, sign):
        return np.median(values) + sign * 3 * median_abs_deviation(values, scale="normal")

    thresholds = mad_thresholds(obs, n_mads=3)
    assert thresholds.mito_upper == pytest.approx(bound(obs["pct_counts_mt"] / 100, +1))
    assert thresholds.genes_lower == pytest.approx(
        np.expm1(bound(obs["log1p_n_genes_by_counts"], -1))
    )
    assert thresholds.genes_upper == pytest.approx(
        np.expm1(bound(obs["log1p_n_genes_by_counts"], +1))
    )
    assert thresholds.total_upper == pytest.approx(np.expm1(bound(obs["log1p_total_counts"], +1)))


def test_mad_thresholds_zero_mad_flags_nothing_for_that_metric():
    obs = pd.DataFrame(
        {
            "pct_counts_mt": [5.0] * 10,
            "log1p_total_counts": np.linspace(7, 8, 10),
            "log1p_n_genes_by_counts": np.linspace(6, 7, 10),
        }
    )
    thresholds = mad_thresholds(obs)
    assert thresholds.mito_upper == np.inf


def test_mad_rule_flags_an_extreme_low_gene_cell():
    rng = np.random.default_rng(2)
    log_genes = rng.normal(7, 0.2, 200)
    log_genes[0] = 3.0
    log_total = rng.normal(8, 0.2, 200)
    obs = pd.DataFrame(
        {
            "pct_counts_mt": rng.uniform(4, 6, 200),
            "log1p_total_counts": log_total,
            "log1p_n_genes_by_counts": log_genes,
            "total_counts": np.expm1(log_total),
            "n_genes_by_counts": np.expm1(log_genes),
        }
    )
    flags = flag_noise(obs, mad_thresholds(obs))
    assert flags.loc[0, "is_low_num_genes"] == 1
    assert flags["is_noise"].mean() < 0.1


def test_relabel_noise_returns_copies_with_rule_applied():
    original = {
        "a": _obs([100, 300], [5.0, 5.0], [1000, 1000]),
        "b": _obs([300, 300], [5.0, 50.0], [1000, 1000]),
    }
    relabelled = relabel_noise(original, lambda obs: THRESHOLDS, "is_noise_alt")
    assert relabelled["a"]["is_noise_alt"].tolist() == [1, 0]
    assert relabelled["b"]["is_noise_alt"].tolist() == [0, 1]
    assert "is_noise_alt" not in original["a"].columns


def test_select_population():
    obs = pd.DataFrame({"is_noise": [1, 0, 0]})
    assert len(select_population(obs, "all")) == 3
    assert len(select_population(obs, "retained")) == 2
    assert len(select_population(obs, "noise")) == 1
    with pytest.raises(ValueError):
        select_population(obs, "bogus")


# ---------------------------------------------------------------------------
# thin_counts / depth matching
# ---------------------------------------------------------------------------


def test_thin_counts_fraction_one_returns_equal_copy():
    X = sp.csr_matrix(np.array([[1, 0, 3], [2, 5, 0]], dtype=np.float32))
    thinned = thin_counts(X, 1.0, np.random.default_rng(0))
    assert (thinned != X).nnz == 0
    assert thinned is not X


def test_thin_counts_scales_totals_and_never_increases_entries():
    X = sp.csr_matrix(np.full((200, 50), 10, dtype=np.float32))
    thinned = thin_counts(X, 0.3, np.random.default_rng(0))
    assert thinned.sum() / X.sum() == pytest.approx(0.3, rel=0.02)
    assert (thinned > X).nnz == 0
    assert thinned.dtype == np.float32


def test_thin_counts_is_deterministic_given_seed():
    X = sp.csr_matrix(np.random.default_rng(0).poisson(4, (30, 20)).astype(float))
    a = thin_counts(X, 0.5, np.random.default_rng(7))
    b = thin_counts(X, 0.5, np.random.default_rng(7))
    assert (a != b).nnz == 0


def test_thin_counts_accepts_dense_input_and_drops_zeros():
    thinned = thin_counts(np.array([[1.0, 1.0], [1.0, 1.0]]), 1e-9, np.random.default_rng(0))
    assert sp.isspmatrix_csr(thinned)
    assert thinned.nnz == 0


@pytest.mark.parametrize("fraction", [0.0, -0.1, 1.5])
def test_thin_counts_rejects_invalid_fraction(fraction):
    with pytest.raises(ValueError, match="fraction"):
        thin_counts(sp.csr_matrix(np.ones((2, 2))), fraction, np.random.default_rng(0))


def test_thin_counts_rejects_non_integer_counts():
    with pytest.raises(ValueError, match="integer"):
        thin_counts(sp.csr_matrix(np.array([[0.5, 1.0]])), 0.5, np.random.default_rng(0))


def test_library_depth_uses_mask():
    adata = AnnData(sp.csr_matrix(np.array([[1.0, 1.0], [5.0, 5.0], [10.0, 10.0]])))
    assert library_depth(adata) == 10.0
    assert library_depth(adata, np.array([True, False, True])) == 11.0


def test_annotate_qc_metrics_supports_small_gene_sets_without_percent_top():
    adata = _count_adata(n_obs=10, n_vars=20)
    annotate_qc_metrics(adata, percent_top=None)
    assert set(QC_METRICS) <= set(adata.obs.columns)
    assert adata.var["mt"].tolist()[:3] == [True, True, False]


def test_depth_matched_qc_obs_hits_target_depth():
    adata = _count_adata()
    depth = library_depth(adata)
    obs, fraction = depth_matched_qc_obs(adata, depth / 2, np.random.default_rng(0))
    assert fraction == pytest.approx(0.5)
    assert obs.index.equals(adata.obs_names)
    assert set(QC_METRICS) <= set(obs.columns)
    assert obs["total_counts"].median() == pytest.approx(depth / 2, rel=0.1)


def test_depth_matched_qc_obs_leaves_shallow_samples_unthinned():
    adata = _count_adata()
    obs, fraction = depth_matched_qc_obs(adata, library_depth(adata) * 10, np.random.default_rng(0))
    assert fraction == 1.0
    np.testing.assert_allclose(obs["total_counts"], np.asarray(adata.X.sum(axis=1)).ravel())


def test_depth_matched_qc_obs_does_not_modify_input():
    adata = _count_adata()
    before = adata.X.copy()
    depth_matched_qc_obs(adata, library_depth(adata) / 3, np.random.default_rng(0))
    assert (adata.X != before).nnz == 0
    assert "total_counts" not in adata.obs.columns
