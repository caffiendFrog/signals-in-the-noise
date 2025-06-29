import pandas as pd

from signals_in_the_noise.utilities.storage import get_data_path, get_resources_path
from signals_in_the_noise.utilities.tenx_genomics import TenX, DirectoryType


class GSE161529:
    """
    Preprocesses the dataset GSE161529 for analysis.

    Pre-requisites for running this analysis:
    - Raw compressed data from study has been expanded into `RAW_DATA_DIRECTORY` in the data directory
    - Raw compressed features file has been downloaded into the data directory
    - Supplementary excel files exist in the resources directory
    """

    STUDY_ID = "GSE161529"
    RAW_DATA_DIRECTORY = f"{STUDY_ID}_RAW"
    FEATURES_FILENAME = f"{STUDY_ID}_features.tsv.gz"
    RESOURCE_FILENAMES = {
        'metadata': f"{STUDY_ID}/table_ev_4.xlsx",
        'qc': f"{STUDY_ID}/table_supplementary_2.xlsx",
    }

    def __init__(self):
        raw_data_directory = get_data_path(self.RAW_DATA_DIRECTORY)
        features_filename = get_data_path(self.FEATURES_FILENAME)
        self.raw_data = TenX(str(raw_data_directory), DirectoryType.MULTIPLE, features_filename=str(features_filename))

        # load the raw data, if it's not been loaded yet
        self.raw_data.load_data()

        # load the metadata information into a dataframe
        metadata_path = get_resources_path(self.RESOURCE_FILENAMES['metadata'])
        # -- header is located on row 3
        self.df_metadata = pd.read_excel(metadata_path, header=2)

        # load the qc information into a dataframe
        qc_path = get_resources_path(self.RESOURCE_FILENAMES['qc'])
        self.df_qc = pd.read_excel(qc_path)
