from pathlib import Path
from typing import Optional

# Resolves to the repository root (two levels above this file: signals_in_the_noise/ → src/ → root)
PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]
DATA_DIRECTORY: Path = PROJECT_ROOT / "data"
RESOURCES_DIRECTORY: Path = PROJECT_ROOT / "resources"


def get_data_path(file_or_directory: Optional[str] = None) -> Path:
    """Return the data directory or a path within it.

    Args:
        file_or_directory: Optional relative path to append to the data directory.

    Returns:
        Absolute Path to the data directory, or to the given sub-path within it.
    """
    if file_or_directory is None:
        return DATA_DIRECTORY
    return DATA_DIRECTORY / file_or_directory


def get_resources_path(file_or_directory: Optional[str] = None) -> Path:
    """Return the resources directory or a path within it.

    Args:
        file_or_directory: Optional relative path to append to the resources directory.

    Returns:
        Absolute Path to the resources directory, or to the given sub-path within it.
    """
    if file_or_directory is None:
        return RESOURCES_DIRECTORY
    return RESOURCES_DIRECTORY / file_or_directory
