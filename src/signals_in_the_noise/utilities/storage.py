from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIRECTORY = PROJECT_ROOT / "data"
RESOURCES_DIRECTORY = PROJECT_ROOT / "resources"

def get_data_path(file_or_directory: str=None):
    if file_or_directory is None:
        return DATA_DIRECTORY
    else:
        return DATA_DIRECTORY / file_or_directory

def get_resources_path(file_or_directory: str=None):
    if file_or_directory is None:
        return RESOURCES_DIRECTORY
    else:
        return RESOURCES_DIRECTORY / file_or_directory
