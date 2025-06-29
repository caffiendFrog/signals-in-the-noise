from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIRECTORY = PROJECT_ROOT / "data"

def get_data_path(file_or_directory: str=None):
    if file_or_directory is None:
        return DATA_DIRECTORY
    else:
        return DATA_DIRECTORY / file_or_directory
