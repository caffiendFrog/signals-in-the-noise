from pathlib import Path

from signals_in_the_noise.config import get_data_path
from signals_in_the_noise.io.tenx import TenX
from signals_in_the_noise.preprocessing.base import Preprocessor


class GSE154932(Preprocessor):
    """Preprocesses dataset GSE154932 for analysis.

    Pre-requisites:
    - Raw compressed data from study has been expanded into ``RAW_DATA_DIRECTORY``
      inside the data directory.
    """

    STUDY_ID = "GSE154932"
    RAW_DATA_DIRECTORY = f"{STUDY_ID}_RAW"

    def __init__(self):
        super().__init__(self.STUDY_ID)
        raw_data_directory = get_data_path("raw") / self.RAW_DATA_DIRECTORY
        raw_data = TenX(str(raw_data_directory), study_id=self.STUDY_ID)
        self.cache_directory_path = Path(raw_data.cache_directory_name)
        if not self.is_data_loaded:
            raw_data.load_data()
            self.data_loaded()
        else:
            raw_data.load_adata()

        for adata in raw_data.multiple_adata:
            self.objects[adata.obs["adata-filename"].iloc[0]] = adata
