"""Tests for signals_in_the_noise.io.tenx."""

from pathlib import Path
from unittest.mock import patch

import pytest

from signals_in_the_noise.io.tenx import TenX


# ---------------------------------------------------------------------------
# TenX.__init__ — validation errors
# ---------------------------------------------------------------------------


def test_tenx_init_raises_when_features_file_missing(tmp_path):
    with pytest.raises(FileNotFoundError, match="Required file not found"):
        TenX(str(tmp_path), features_filename="nonexistent_features.tsv.gz")


def test_tenx_init_raises_on_wrong_features_format(tmp_path):
    bad_features = tmp_path / "GSE161529_wrong.tsv.gz"
    bad_features.touch()
    with pytest.raises(ValueError, match="expected format"):
        TenX(str(tmp_path), features_filename=str(bad_features))


def test_tenx_init_requires_study_id_when_no_features_and_directory_not_raw(tmp_path):
    with pytest.raises(ValueError, match="study_id must be provided"):
        TenX(str(tmp_path))


def test_tenx_init_per_sample_features_derives_study_id_from_raw_directory(tmp_path):
    raw_dir = tmp_path / "GSE154932_RAW"
    raw_dir.mkdir()

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(raw_dir))

    assert tenx.study_id == "GSE154932"
    assert tenx.features_path is None
    assert tenx.directory == raw_dir


def test_tenx_init_per_sample_features_with_explicit_study_id(tmp_path):
    raw_dir = tmp_path / "raw_data"
    raw_dir.mkdir()

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(raw_dir), study_id="MYSTUDY")

    assert tenx.study_id == "MYSTUDY"
    assert tenx.features_path is None


# ---------------------------------------------------------------------------
# TenX.__init__ — happy path
# ---------------------------------------------------------------------------


def test_tenx_init_happy_path(tmp_path):
    features_file = tmp_path / "GSE161529_features.tsv.gz"
    features_file.touch()

    with patch("signals_in_the_noise.io.tenx.get_data_path") as mock_get_data_path:
        study_dir = tmp_path / "GSE161529"
        mock_get_data_path.return_value = study_dir

        tenx = TenX(str(tmp_path), features_filename=str(features_file))

    assert tenx.study_id == "GSE161529"
    assert tenx.features_path == features_file
    assert tenx.directory == tmp_path
    assert tenx.multiple_adata == []


# ---------------------------------------------------------------------------
# TenX.load_adata — missing cache directory
# ---------------------------------------------------------------------------


def test_load_adata_does_nothing_when_cache_directory_missing(tmp_path):
    features_file = tmp_path / "STUDY_features.tsv.gz"
    features_file.touch()

    def fake_get_data_path(subpath=None):
        return tmp_path / subpath if subpath else tmp_path

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(tmp_path), features_filename=str(features_file))
        # cache directory is not created — load_adata should be a no-op
        tenx.load_adata()

    assert tenx.multiple_adata == []


# ---------------------------------------------------------------------------
# TenX.cache_directory_name
# ---------------------------------------------------------------------------


def test_cache_directory_name_uses_study_id(tmp_path):
    features_file = tmp_path / "MYSTUDY_features.tsv.gz"
    features_file.touch()

    expected_cache = tmp_path / "MYSTUDY_adata_cache"

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(tmp_path), features_filename=str(features_file))
        result = tenx.cache_directory_name

    assert result == expected_cache


# ---------------------------------------------------------------------------
# TenX._samples_to_file_dictionary
# ---------------------------------------------------------------------------


def test_samples_to_file_dictionary_matches_barcodes_and_matrix(tmp_path):
    features_file = tmp_path / "STUDY_features.tsv.gz"
    features_file.touch()

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    (raw_dir / "SAMPLE1-barcodes.tsv.gz").touch()
    (raw_dir / "SAMPLE1-matrix.mtx.gz").touch()
    (raw_dir / "unrelated_file.txt").touch()

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(raw_dir), features_filename=str(features_file))
        result = tenx._samples_to_file_dictionary()

    assert "SAMPLE1" in result
    filenames = result["SAMPLE1"]
    assert any("barcodes" in f for f in filenames)
    assert any("matrix" in f for f in filenames)
    assert "unrelated_file.txt" not in filenames


def test_samples_to_file_dictionary_matches_per_sample_features(tmp_path):
    raw_dir = tmp_path / "GSE154932_RAW"
    raw_dir.mkdir()
    (raw_dir / "SAMPLE1-barcodes.tsv.gz").touch()
    (raw_dir / "SAMPLE1-matrix.mtx.gz").touch()
    (raw_dir / "SAMPLE1-features.tsv.gz").touch()

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(raw_dir))
        result = tenx._samples_to_file_dictionary()

    assert "SAMPLE1" in result
    filenames = result["SAMPLE1"]
    assert any("barcodes" in f for f in filenames)
    assert any("matrix" in f for f in filenames)
    assert any("features" in f for f in filenames)


def test_samples_to_file_dictionary_ignores_non_matching_files(tmp_path):
    features_file = tmp_path / "STUDY_features.tsv.gz"
    features_file.touch()

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    (raw_dir / "readme.txt").touch()
    (raw_dir / "STUDY_features.tsv.gz").touch()

    def fake_get_data_path(subpath=None):
        if subpath is None:
            return tmp_path
        return tmp_path / subpath

    with patch("signals_in_the_noise.io.tenx.get_data_path", side_effect=fake_get_data_path):
        tenx = TenX(str(raw_dir), features_filename=str(features_file))
        result = tenx._samples_to_file_dictionary()

    assert len(result) == 0
