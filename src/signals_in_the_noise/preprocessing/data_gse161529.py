from functools import reduce

import pandas as pd
import scanpy as sc
from slugify import slugify

from signals_in_the_noise.preprocessing.prep_config import Prep
from signals_in_the_noise.utilities.logging import get_logger
from signals_in_the_noise.utilities.storage import get_data_path, get_resources_path
from signals_in_the_noise.utilities.tenx_genomics import TenX, DirectoryType


L = get_logger(__name__)


class GSE161529(Prep):
    """
    Preprocesses the dataset GSE161529 for analysis.

    Pre-requisites for running this analysis:
    - Raw compressed data from study has been expanded into `RAW_DATA_DIRECTORY` in the data directory
    - Raw compressed features file has been downloaded into the data directory
    - Supplementary Excel files exist in the resources directory
    """

    STUDY_ID = "GSE161529"
    RAW_DATA_DIRECTORY = f"{STUDY_ID}_RAW"
    FEATURES_FILENAME = f"{STUDY_ID}_features.tsv.gz"

    def __init__(self):
        super().__init__(self.STUDY_ID)
        raw_data_directory = get_data_path(self.RAW_DATA_DIRECTORY)
        features_filename = get_data_path(self.FEATURES_FILENAME)
        raw_data = TenX(str(raw_data_directory), DirectoryType.MULTIPLE, features_filename=str(features_filename))

        if not self.is_data_loaded:
            raw_data.load_data()
            self.data_loaded()
        else:
            raw_data.load_adata()

        self.objects = raw_data.multiple_adata
        #
        # self.add_annotations()

    def add_annotations(self):
        """
        Adds annotations from resource tables to the anndata objects for the raw data.
        :return:
        """
        resource_df = None
        resource_df_filename = get_data_path(f"{self.STUDY_ID}_annotations_df.csv")
        if not self.is_annotations_loaded:
            annotation_resources = {
                'metadata': (f"{self.STUDY_ID}/table_supplementary_1.xlsx", 0),
                'phenotype': (f"{self.STUDY_ID}/table_ev_4.xlsx", 2),
                'qc': (f"{self.STUDY_ID}/table_supplementary_2.xlsx", 0),
            }
            resource_df = self._prepare_resources_for_annotation(annotation_resources)
            resource_df.to_csv(resource_df_filename, index=False)
            self.annotations_loaded()
        else:
            resource_df = pd.read_csv(resource_df_filename, header=0)

        if not self.is_annotations_added:
            self.raw_data.multiple_adata = []
            annotation_column_names = {
                'title': "title",
                'menopause_status': slugify("menopause status"),
                'cancer_type': slugify("cancer type"),
                'cell_population': slugify("cell population"),
                'num_cells_before': slugify("number of cells"),
                'num_cells_after': slugify("number of cells after filtering"),
                'num_genes_before': slugify("number of genes detected"),
                'num_genes_after': slugify("# genes detected after filtering"),
                'qc_mito_upper': slugify("mito - upper"),
                'qc_genes_lower': slugify("# genes - lower"),
                'qc_genes_upper': slugify("# genes - upper"),
                'qc_size_upper': slugify("library size - upper"),
                # some columns are duplicated in the resource files
                # when they were joined, pandas appended suffices to distinguish, e.g. _x, _y
                'gender': "gender_x",
                'parity': "parity_x",
            }
            for index in range(len(resource_df)):
                filename = get_data_path(f"{self.STUDY_ID}_adata_cache/{resource_df.loc[index, 'adata-filename']}")
                adata = sc.read_h5ad(filename)
                for uns_name, column_name in annotation_column_names.items():
                    adata.uns[uns_name] = resource_df.loc[index, column_name]
                self.raw_data.multiple_adata.append(adata)
                adata.write_h5ad(filename)

            self.annotations_loaded()

    def _prepare_resources_for_annotation(self, resources):
        """
        Loads the resources and prepares them for use to annotate the anndata objects

        :param resources: dictionary in the format of {name: (file, header_row)}
        :return: prepared resources dataframe
        """

        # Read the files into dataframes
        resource_dfs = []
        for resource_name, (resource, header) in resources.items():
            L.info(f"Reading {resource_name} from {resource}...")
            resource_path = get_resources_path(resource)
            resource_df = pd.read_excel(resource_path, header=header)
            # lowercase the column names for consistency
            resource_df.columns = resource_df.columns.str.lower()
            resource_dfs.append(resource_df)

        # Join all the dataframes on the sample name
        join_column = "sample-name"
        resource_df = reduce(lambda left, right: pd.merge(left, right, on=join_column), resource_dfs)
        # slugify the column names for easier use downstream
        resource_df.columns = [slugify(column) for column in resource_df.columns]

        # Add column that is the expected `.h5ad` file name
        resource_df['sample-suffix'] = resource_df['barcodes-file'].str.replace('-barcodes.tsv.gz','')
        resource_df['adata-filename'] = (
            resource_df['geo-id'].astype(str) + '_' +
            resource_df['sample-suffix'].astype(str) + '.h5ad'
        )

        return resource_df

