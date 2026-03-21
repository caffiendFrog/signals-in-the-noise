"""Tests for signals_in_the_noise.config."""

from pathlib import Path

from signals_in_the_noise.config import (
    PROJECT_ROOT,
    DATA_DIRECTORY,
    RESOURCES_DIRECTORY,
    get_data_path,
    get_resources_path,
)


def test_project_root_is_absolute():
    assert PROJECT_ROOT.is_absolute()


def test_project_root_contains_pyproject_toml():
    """Verify the root resolves to the actual repository root."""
    assert (PROJECT_ROOT / "pyproject.toml").exists()


def test_data_directory_is_under_project_root():
    assert DATA_DIRECTORY == PROJECT_ROOT / "data"


def test_resources_directory_is_under_project_root():
    assert RESOURCES_DIRECTORY == PROJECT_ROOT / "resources"


def test_get_data_path_no_arg_returns_data_directory():
    result = get_data_path()
    assert result == DATA_DIRECTORY


def test_get_data_path_with_arg_appends_to_data_directory():
    result = get_data_path("raw/some_file.tar")
    assert result == DATA_DIRECTORY / "raw" / "some_file.tar"


def test_get_data_path_returns_path_object():
    assert isinstance(get_data_path(), Path)
    assert isinstance(get_data_path("subdir"), Path)


def test_get_resources_path_no_arg_returns_resources_directory():
    result = get_resources_path()
    assert result == RESOURCES_DIRECTORY


def test_get_resources_path_with_arg_appends_to_resources_directory():
    result = get_resources_path("GSE161529/gene_sets.gmt")
    assert result == RESOURCES_DIRECTORY / "GSE161529" / "gene_sets.gmt"


def test_get_resources_path_returns_path_object():
    assert isinstance(get_resources_path(), Path)
    assert isinstance(get_resources_path("subdir"), Path)
