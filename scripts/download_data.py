"""Download and extract the GSE161529 raw data files into data/raw/."""

import logging
import sys
import tarfile
from pathlib import Path
from urllib import request

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from signals_in_the_noise.utils.logging_config import setup_logging

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = PROJECT_ROOT / "data" / "raw"

RAW_DATA_FILE = DATA_DIRECTORY / "GSE161529_RAW.tar"
FEATURES_FILE = DATA_DIRECTORY / "GSE161529_features.tsv.gz"

RAW_DATA_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161529&format=file"
FEATURES_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSE161529&format=file&file=GSE161529%5Ffeatures%2Etsv%2Egz"
)


def download_file(url: str, dest: Path) -> None:
    """Download a file from ``url`` to ``dest`` if it does not already exist.

    Args:
        url: Source URL.
        dest: Destination path (parent directories must exist).
    """
    if dest.exists():
        logger.info(f"File already exists, skipping download: {dest}")
        return
    logger.info(f"Downloading {url} → {dest}")
    request.urlretrieve(url, dest)
    logger.info("Download complete.")


def extract_tar(archive: Path, extract_to: Path) -> None:
    """Extract a tar archive to the given directory.

    Args:
        archive: Path to the tar archive.
        extract_to: Destination directory (created if needed).
    """
    extract_to.mkdir(parents=True, exist_ok=True)
    mode = "r:gz" if archive.suffix == ".gz" else "r:"
    logger.info(f"Extracting {archive} → {extract_to} (mode={mode})")
    with tarfile.open(archive, mode) as tar:
        tar.extractall(extract_to)
    logger.info("Extraction complete.")


def main() -> None:
    """Download raw data and features file; extract the tar archive."""
    setup_logging()

    DATA_DIRECTORY.mkdir(parents=True, exist_ok=True)

    download_file(RAW_DATA_URL, RAW_DATA_FILE)
    extract_to = DATA_DIRECTORY / "GSE161529_RAW"
    if not extract_to.exists():
        extract_tar(RAW_DATA_FILE, extract_to)

    download_file(FEATURES_URL, FEATURES_FILE)


if __name__ == "__main__":
    main()
