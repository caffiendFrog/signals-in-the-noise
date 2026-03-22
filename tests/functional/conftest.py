"""Fixtures shared across functional tests."""

from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Values that mirror scripts/generate_fixtures.py exactly.
# If you change the generator, update these constants too.
# ---------------------------------------------------------------------------

FIXTURE_TENX_DIR = Path(__file__).resolve().parents[2] / "data" / "fixtures" / "tenx"

EXPECTED_N_GENES = 30
EXPECTED_N_CELLS = 5
EXPECTED_GENE_NAMES = [f"GENE{i:04d}" for i in range(EXPECTED_N_GENES)]
EXPECTED_SAMPLES = {
    "SAMPLEONE": [
        "AAACATACAACCAC-1",
        "AAACATTGAGCTAC-1",
        "AAACATTGATCAGC-1",
        "AAACCGTGCTTCCG-1",
        "AAACCGTGTATGCG-1",
    ],
    "SAMPLETWO": [
        "AAACATACAAGGGT-1",
        "AAACATACAGAGCC-1",
        "AAACATACAGCATA-1",
        "AAACATACAGCCTG-1",
        "AAACATACAGTAAG-1",
    ],
}


@pytest.fixture(scope="session")
def tenx_fixture_dir() -> Path:
    """Return the tenx fixture directory and assert it exists."""
    if not FIXTURE_TENX_DIR.exists():
        pytest.fail(
            f"Fixture directory not found: {FIXTURE_TENX_DIR}\n"
            "Run 'python scripts/generate_fixtures.py' to create it."
        )
    return FIXTURE_TENX_DIR


@pytest.fixture(scope="session")
def tenx_raw_dir(tenx_fixture_dir: Path) -> Path:
    """Return the directory containing the per-sample raw files."""
    return tenx_fixture_dir / "raw"


@pytest.fixture(scope="session")
def tenx_features_file(tenx_fixture_dir: Path) -> Path:
    """Return the path to the shared features file."""
    return tenx_fixture_dir / "FIXTURE_features.tsv.gz"
