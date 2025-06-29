from signals_in_the_noise.utilities.storage import get_data_path
from signals_in_the_noise.utilities.tenx_genomics import TenX, DirectoryType


class GSE161529:
    """
    Preprocesses the dataset GSE161529 for analysis.

    Pre-requisites for running this analysis:
    - Raw compressed data from study has been expanded into `RAW_DATA_DIRECTORY` in the data directory
    - Raw compressed features file has been downloaded into the data directory
    - Supplementary excel files exist in the resources directory
    """

    RAW_DATA_DIRECTORY = "GSE161529_RAW"
    FEATURES_FILENAME = "GSE161529_features.tsv.gz"
    RESOURCES_DIRECTORY = "resources/GSE161529"

    def __init__(self):
        self.resources_directory = get_data_path(self.RESOURCES_DIRECTORY)

        raw_data_directory = get_data_path(self.RAW_DATA_DIRECTORY)
        features_filename = get_data_path(self.FEATURES_FILENAME)
        self.raw_data = TenX(str(raw_data_directory), DirectoryType.MULTIPLE, features_filename=str(features_filename))

        # load the raw data, if it's not been loaded yet
        self.raw_data.load_data()
