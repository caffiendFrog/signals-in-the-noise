"""Tests for signals_in_the_noise.preprocessing.gse161529."""

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

from signals_in_the_noise.preprocessing.gse161529 import GSE161529


# ---------------------------------------------------------------------------
# Class-level constants — no I/O required
# ---------------------------------------------------------------------------


def test_study_id_constant():
    assert GSE161529.STUDY_ID == "GSE161529"


def test_raw_data_directory_derived_from_study_id():
    assert GSE161529.RAW_DATA_DIRECTORY == "GSE161529_RAW"


def test_features_filename_derived_from_study_id():
    assert GSE161529.FEATURES_FILENAME == "GSE161529_features.tsv.gz"


def test_expected_mismatches_is_list_of_strings():
    assert isinstance(GSE161529.EXPECTED_MISMATCHES, list)
    assert all(isinstance(s, str) for s in GSE161529.EXPECTED_MISMATCHES)


# ---------------------------------------------------------------------------
# GSE161529._apply_one — static method, no disk I/O
# ---------------------------------------------------------------------------


def _make_mock_adata(n_obs: int = 10, n_vars: int = 500) -> AnnData:
    """Create a minimal AnnData with QC thresholds in uns.

    n_vars must be >= 500 because sc.pp.calculate_qc_metrics internally
    computes top-segment proportions at [50, 100, 200, 500].
    """
    X = sp.csr_matrix(np.random.default_rng(0).integers(0, 10, size=(n_obs, n_vars)).astype(float))
    adata = AnnData(X)
    adata.var_names = [f"GENE_{i}" for i in range(n_vars)]
    adata.obs_names = [f"CELL_{i}" for i in range(n_obs)]
    adata.obs["adata-filename"] = "test_sample.h5ad"
    adata.uns["qc_genes_lower"] = 0
    adata.uns["qc_genes_upper"] = 10000
    adata.uns["qc_mito_upper"] = 1.0  # 100 % — nothing is high-mito
    adata.uns["qc_total_upper"] = 10_000_000  # nothing is high-total
    adata.uns["num_cells_after"] = n_obs  # all cells are "real"
    return adata


def test_apply_one_adds_is_noise_column():
    adata = _make_mock_adata()
    GSE161529._apply_one(adata)
    assert "is_noise" in adata.obs.columns


def test_apply_one_adds_all_qc_flag_columns():
    adata = _make_mock_adata()
    GSE161529._apply_one(adata)
    expected_cols = {
        "is_low_num_genes",
        "is_high_num_genes",
        "is_high_mito",
        "is_high_total_count",
        "is_noise",
        "zero_genes",
        "zero_mito",
        "zero_count",
    }
    assert expected_cols.issubset(set(adata.obs.columns))


def test_apply_one_raises_on_count_mismatch():
    adata = _make_mock_adata(n_obs=10)
    # Set expected count to something that will never match
    adata.uns["num_cells_after"] = 999
    with pytest.raises(ValueError, match="Check failed"):
        GSE161529._apply_one(adata)


def test_apply_one_skips_count_check_for_whitelisted_file():
    adata = _make_mock_adata(n_obs=10)
    adata.uns["num_cells_after"] = 999
    adata.obs["adata-filename"] = GSE161529.EXPECTED_MISMATCHES[0]
    # Should not raise even though counts do not match
    GSE161529._apply_one(adata)


def test_apply_one_is_noise_is_binary():
    adata = _make_mock_adata()
    GSE161529._apply_one(adata)
    unique_values = set(adata.obs["is_noise"].unique())
    assert unique_values.issubset({0, 1})


def test_apply_one_annotates_mitochondrial_genes():
    adata = _make_mock_adata()
    # Replace first two var names to include one mito and one non-mito gene.
    new_names = ["MT-CO1", "BRCA1"] + [f"GENE_{i}" for i in range(2, adata.n_vars)]
    adata.var_names = new_names
    GSE161529._apply_one(adata)
    assert bool(adata.var["mt"].iloc[0]) is True
    assert bool(adata.var["mt"].iloc[1]) is False
