"""Functional tests for signals_in_the_noise.io.tenx.

These tests exercise the full TenX loading pipeline against small, committed
fixture files in data/fixtures/tenx/.  No mocking of scanpy or file I/O
occurs — the fixtures are read exactly as real data would be.

The only patching done is to redirect ``get_data_path`` so that TenX writes
its reorganised sample directories and h5ad cache into pytest's tmp_path
instead of the project's data/ directory.

To regenerate the fixture files:

    python scripts/generate_fixtures.py
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from signals_in_the_noise.io.tenx import TenX, DirectoryType
from tests.functional.conftest import (
    EXPECTED_N_CELLS,
    EXPECTED_N_GENES,
    EXPECTED_GENE_NAMES,
    EXPECTED_SAMPLES,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_tenx(raw_dir: Path, features_file: Path, tmp_path: Path) -> TenX:
    """Construct a TenX instance whose output directories land in *tmp_path*."""

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        return TenX(
            str(raw_dir),
            DirectoryType.MULTIPLE,
            features_filename=str(features_file),
        )


# ---------------------------------------------------------------------------
# Fixture loading — no cache
# ---------------------------------------------------------------------------


def test_load_data_returns_one_adata_per_sample(tmp_path, tenx_raw_dir, tenx_features_file):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=False)

    assert len(tenx.multiple_adata) == len(EXPECTED_SAMPLES)


def test_load_data_adata_shape_matches_fixture(tmp_path, tenx_raw_dir, tenx_features_file):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=False)

    for adata in tenx.multiple_adata:
        assert adata.n_obs == EXPECTED_N_CELLS, (
            f"Expected {EXPECTED_N_CELLS} cells, got {adata.n_obs}"
        )
        assert adata.n_vars == EXPECTED_N_GENES, (
            f"Expected {EXPECTED_N_GENES} genes, got {adata.n_vars}"
        )


def test_load_data_gene_names_match_fixture(tmp_path, tenx_raw_dir, tenx_features_file):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=False)

    for adata in tenx.multiple_adata:
        assert list(adata.var_names) == EXPECTED_GENE_NAMES


def test_load_data_barcodes_match_fixture(tmp_path, tenx_raw_dir, tenx_features_file):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=False)

    # Collect all barcode lists across loaded adatas
    loaded_barcodes = {
        frozenset(adata.obs_names.tolist()) for adata in tenx.multiple_adata
    }
    expected_barcodes = {
        frozenset(barcodes) for barcodes in EXPECTED_SAMPLES.values()
    }
    assert loaded_barcodes == expected_barcodes


def test_load_data_sets_adata_filename_obs_column(tmp_path, tenx_raw_dir, tenx_features_file):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=False)

    for adata in tenx.multiple_adata:
        assert "adata-filename" in adata.obs.columns
        filename = adata.obs["adata-filename"].iloc[0]
        assert filename.endswith(".h5ad")
        # filename should be one of the expected sample IDs
        sample_id = filename.replace(".h5ad", "")
        assert sample_id in EXPECTED_SAMPLES


def test_load_data_count_matrix_is_non_negative(tmp_path, tenx_raw_dir, tenx_features_file):
    import numpy as np

    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=False)

    for adata in tenx.multiple_adata:
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
        assert np.all(X >= 0), "Count matrix contains negative values"


# ---------------------------------------------------------------------------
# Fixture loading — with cache
# ---------------------------------------------------------------------------


def test_load_data_with_cache_creates_h5ad_files(tmp_path, tenx_raw_dir, tenx_features_file):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        tenx = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                    features_filename=str(tenx_features_file))
        tenx.load_data(cache=True)

        cache_dir = tenx.cache_directory_name
        h5ad_files = list(cache_dir.glob("*.h5ad"))

    assert len(h5ad_files) == len(EXPECTED_SAMPLES)
    expected_names = {f"{sid}.h5ad" for sid in EXPECTED_SAMPLES}
    actual_names = {f.name for f in h5ad_files}
    assert actual_names == expected_names


def test_load_data_with_cache_then_load_adata_gives_same_shapes(
    tmp_path, tenx_raw_dir, tenx_features_file
):
    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        # First pass: load from raw files and cache
        tenx_first = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                          features_filename=str(tenx_features_file))
        tenx_first.load_data(cache=True)

        # Second pass: load from cache via load_adata()
        tenx_second = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                           features_filename=str(tenx_features_file))
        tenx_second.load_adata()

    assert len(tenx_second.multiple_adata) == len(EXPECTED_SAMPLES)
    for adata in tenx_second.multiple_adata:
        assert adata.n_obs == EXPECTED_N_CELLS
        assert adata.n_vars == EXPECTED_N_GENES


def test_load_data_with_cache_skips_already_cached_samples(
    tmp_path, tenx_raw_dir, tenx_features_file
):
    """Second call with cache=True must skip reconstitution for cached samples."""

    def fake_gdp(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_gdp):
        # First load populates the cache
        tenx1 = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                     features_filename=str(tenx_features_file))
        tenx1.load_data(cache=True)

        # Second load with cache=True; should skip and load from skipped_files
        tenx2 = TenX(str(tenx_raw_dir), DirectoryType.MULTIPLE,
                     features_filename=str(tenx_features_file))
        tenx2.load_data(cache=True)

    assert len(tenx2.multiple_adata) == len(EXPECTED_SAMPLES)
    for adata in tenx2.multiple_adata:
        assert adata.n_obs == EXPECTED_N_CELLS
        assert adata.n_vars == EXPECTED_N_GENES


# ---------------------------------------------------------------------------
# Fixture file integrity
# ---------------------------------------------------------------------------


def test_fixture_files_exist(tenx_fixture_dir, tenx_raw_dir, tenx_features_file):
    """Guard against accidentally deleted or missing fixtures."""
    assert tenx_features_file.exists(), f"Missing: {tenx_features_file}"
    for sample_id in EXPECTED_SAMPLES:
        assert (tenx_raw_dir / f"{sample_id}-barcodes.tsv.gz").exists()
        assert (tenx_raw_dir / f"{sample_id}-matrix.mtx.gz").exists()


def test_fixture_features_file_contains_expected_genes(tenx_features_file):
    import gzip

    with gzip.open(tenx_features_file, "rt", encoding="utf-8") as fh:
        lines = [line.strip() for line in fh if line.strip()]

    assert len(lines) == EXPECTED_N_GENES
    gene_names_in_file = [line.split("\t")[1] for line in lines]
    assert gene_names_in_file == EXPECTED_GENE_NAMES


def test_fixture_barcodes_files_contain_expected_cells(tenx_raw_dir):
    import gzip

    for sample_id, expected_barcodes in EXPECTED_SAMPLES.items():
        path = tenx_raw_dir / f"{sample_id}-barcodes.tsv.gz"
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            loaded = [line.strip() for line in fh if line.strip()]
        assert loaded == expected_barcodes, f"Barcodes mismatch for {sample_id}"
