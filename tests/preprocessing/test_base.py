"""Tests for signals_in_the_noise.preprocessing.base."""

import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from anndata import AnnData

from signals_in_the_noise.preprocessing.base import PreprocessorConfig, Preprocessor


# ---------------------------------------------------------------------------
# PreprocessorConfig — serialisation roundtrip
# ---------------------------------------------------------------------------


def test_preprocessor_config_to_json_creates_file(tmp_path):
    config = PreprocessorConfig(data_loaded=True, annotations_loaded=False, annotations_applied=False)
    dest = tmp_path / "config.json"
    config.to_json(dest)
    assert dest.exists()


def test_preprocessor_config_to_json_roundtrip(tmp_path):
    original = PreprocessorConfig(
        data_loaded=True,
        annotations_loaded=True,
        annotations_applied=False,
        custom=["step_one", "step_two"],
    )
    dest = tmp_path / "config.json"
    original.to_json(dest)

    restored = PreprocessorConfig.from_json(dest)
    assert restored.data_loaded is True
    assert restored.annotations_loaded is True
    assert restored.annotations_applied is False
    assert restored.custom == ["step_one", "step_two"]


def test_preprocessor_config_to_json_creates_parent_dirs(tmp_path):
    config = PreprocessorConfig(False, False, False)
    nested = tmp_path / "a" / "b" / "c" / "config.json"
    config.to_json(nested)
    assert nested.exists()


def test_preprocessor_config_from_json_raises_on_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        PreprocessorConfig.from_json(tmp_path / "nonexistent.json")


# ---------------------------------------------------------------------------
# Preprocessor — initialisation
# ---------------------------------------------------------------------------


def _make_preprocessor(tmp_path: Path, study_id: str = "TEST_STUDY") -> Preprocessor:
    """Construct a Preprocessor whose config file lives in tmp_path."""
    with patch("signals_in_the_noise.preprocessing.base.get_data_path") as mock:
        mock.return_value = tmp_path / f"{study_id}.json"
        p = Preprocessor(study_id)
    return p


def test_preprocessor_creates_fresh_config_when_none_exists(tmp_path):
    p = _make_preprocessor(tmp_path)
    assert p.config.data_loaded is False
    assert p.config.annotations_loaded is False
    assert p.config.annotations_applied is False


def test_preprocessor_loads_existing_config(tmp_path):
    config_data = {
        "data_loaded": True,
        "annotations_loaded": True,
        "annotations_applied": False,
        "custom": ["foo"],
    }
    config_file = tmp_path / "TEST_STUDY.json"
    config_file.write_text(json.dumps(config_data), encoding="utf-8")

    with patch("signals_in_the_noise.preprocessing.base.get_data_path") as mock:
        mock.return_value = config_file
        p = Preprocessor("TEST_STUDY")

    assert p.config.data_loaded is True
    assert p.config.annotations_loaded is True
    assert p.config.custom == ["foo"]


# ---------------------------------------------------------------------------
# Preprocessor — state-flag methods
# ---------------------------------------------------------------------------


def test_data_loaded_sets_flag_and_persists(tmp_path):
    p = _make_preprocessor(tmp_path)
    assert p.is_data_loaded is False
    p.data_loaded()
    assert p.is_data_loaded is True
    restored = PreprocessorConfig.from_json(p.config_path)
    assert restored.data_loaded is True


def test_annotations_loaded_sets_flag_and_persists(tmp_path):
    p = _make_preprocessor(tmp_path)
    p.annotations_loaded()
    assert p.is_annotations_loaded is True


def test_annotations_applied_sets_flag_and_persists(tmp_path):
    p = _make_preprocessor(tmp_path)
    p.annotations_applied()
    assert p.is_annotations_applied is True


def test_has_custom_returns_false_when_absent(tmp_path):
    p = _make_preprocessor(tmp_path)
    assert p.has_custom("missing_step") is False


def test_add_custom_and_has_custom(tmp_path):
    p = _make_preprocessor(tmp_path)
    p.add_custom("my_step")
    assert p.has_custom("my_step") is True


def test_add_custom_persists_to_disk(tmp_path):
    p = _make_preprocessor(tmp_path)
    p.add_custom("persisted_step")
    restored = PreprocessorConfig.from_json(p.config_path)
    assert "persisted_step" in restored.custom


# ---------------------------------------------------------------------------
# Preprocessor — get_dataset
# ---------------------------------------------------------------------------


def test_get_dataset_returns_empty_adata_for_unknown_key(tmp_path):
    p = _make_preprocessor(tmp_path)
    result = p.get_dataset("nonexistent.h5ad")
    assert isinstance(result, AnnData)


def test_get_dataset_returns_copy_when_found(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)
    stored = AnnData(sp.csr_matrix(np.ones((3, 4))))
    p.objects["sample.h5ad"] = stored

    result = p.get_dataset("sample.h5ad")
    assert result is not stored
    assert result.shape == stored.shape


# ---------------------------------------------------------------------------
# Preprocessor — random_seed
# ---------------------------------------------------------------------------


def test_random_seed_returns_integer(tmp_path):
    p = _make_preprocessor(tmp_path)
    seed = p.random_seed
    assert isinstance(seed, int)
    assert 1 <= seed <= 100


def test_random_seed_is_stable_across_calls(tmp_path):
    p = _make_preprocessor(tmp_path)
    first = p.random_seed
    second = p.random_seed
    assert first == second


# ---------------------------------------------------------------------------
# Preprocessor — cache_raw_gene_expression
# ---------------------------------------------------------------------------


def test_cache_raw_gene_expression_adds_obs_columns(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)

    adata = AnnData(sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]])))
    adata.var_names = ["GENE_A", "GENE_B"]

    genes_of_interest = {"gene_a_expr": "GENE_A"}
    result = p.cache_raw_gene_expression(adata, genes_of_interest)

    assert "gene_a_expr" in result.obs.columns
    np.testing.assert_array_almost_equal(result.obs["gene_a_expr"].values, [1.0, 0.0])


def test_cache_raw_gene_expression_zero_fills_missing_gene(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)

    adata = AnnData(sp.csr_matrix(np.array([[1.0], [2.0]])))
    adata.var_names = ["GENE_A"]

    result = p.cache_raw_gene_expression(adata, {"missing_expr": "MISSING_GENE"})
    np.testing.assert_array_equal(result.obs["missing_expr"].values, [0.0, 0.0])


def test_cache_raw_gene_expression_does_not_modify_original_by_default(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)

    adata = AnnData(sp.csr_matrix(np.ones((3, 2))))
    adata.var_names = ["A", "B"]
    original_obs_cols = list(adata.obs.columns)

    p.cache_raw_gene_expression(adata, {"a_expr": "A"})
    assert list(adata.obs.columns) == original_obs_cols


# ---------------------------------------------------------------------------
# Preprocessor — check_adata_for_genes
# ---------------------------------------------------------------------------


def test_check_adata_for_genes_returns_missing_genes(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)

    adata = AnnData(sp.csr_matrix(np.ones((2, 2))))
    adata.var_names = ["BRCA1", "TP53"]

    missing = p.check_adata_for_genes(adata, ["BRCA1", "NONEXISTENT"])
    assert missing == ["NONEXISTENT"]


def test_check_adata_for_genes_is_case_insensitive(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)

    adata = AnnData(sp.csr_matrix(np.ones((2, 2))))
    adata.var_names = ["BRCA1", "TP53"]

    missing = p.check_adata_for_genes(adata, ["brca1", "tp53"])
    assert missing == []


def test_check_adata_for_genes_returns_empty_when_all_present(tmp_path):
    import scipy.sparse as sp

    p = _make_preprocessor(tmp_path)

    adata = AnnData(sp.csr_matrix(np.ones((2, 1))))
    adata.var_names = ["GENE_X"]

    missing = p.check_adata_for_genes(adata, ["GENE_X"])
    assert missing == []
