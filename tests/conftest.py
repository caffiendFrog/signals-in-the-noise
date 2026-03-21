"""Shared pytest fixtures."""

import json
from pathlib import Path

import pytest


FIXTURE_DIR = Path(__file__).parent.parent / "data" / "fixtures"


@pytest.fixture()
def fixture_dir() -> Path:
    """Return the path to the small test-fixture data directory."""
    return FIXTURE_DIR


@pytest.fixture()
def preprocessor_config_dict() -> dict:
    """Return a minimal valid PreprocessorConfig as a plain dict."""
    return {
        "data_loaded": False,
        "annotations_loaded": False,
        "annotations_applied": False,
        "custom": [],
    }


@pytest.fixture()
def config_json_file(tmp_path, preprocessor_config_dict) -> Path:
    """Write a PreprocessorConfig JSON file to a temp directory and return its path."""
    config_file = tmp_path / "test_study.json"
    config_file.write_text(json.dumps(preprocessor_config_dict), encoding="utf-8")
    return config_file
