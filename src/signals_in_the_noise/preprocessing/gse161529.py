from collections import defaultdict
from functools import reduce
from pathlib import Path

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

    EXPECTED_MISMATCHES = [
        'GSM4909296_ER-MH0001.h5ad',
        'GSM4909313_ER-MH0064-T.h5ad',
        'GSM4909319_mER-PM0178.h5ad',
    ]

    def __init__(self):
        super().__init__(self.STUDY_ID)
        raw_data_directory = get_data_path(self.RAW_DATA_DIRECTORY)
        features_filename = get_data_path(self.FEATURES_FILENAME)
        raw_data = TenX(str(raw_data_directory), DirectoryType.MULTIPLE, features_filename=str(features_filename))
        self.cache_directory_path = Path(raw_data.cache_directory_name)
        if not self.is_data_loaded:
            raw_data.load_data()
            self.data_loaded()
        else:
            raw_data.load_adata()

        self.objects = defaultdict()
        for index, adata in enumerate(raw_data.multiple_adata):
            self.objects[adata.uns['adata-filename']] = adata

        if not self.is_annotations_loaded:
            self._load_annotations()

        if not self.is_annotations_applied:
            self._apply_annotations()

    def _load_annotations(self):
        """
        Adds annotations from resource tables to the anndata objects for the raw data.
        :return:
        """
        resource_df_filename = get_data_path(f"{self.STUDY_ID}_annotations_df.csv")

        annotation_resources = {
            'metadata': (f"{self.STUDY_ID}/table_supplementary_1.xlsx", 0),
            'phenotype': (f"{self.STUDY_ID}/table_ev_4.xlsx", 2),
            'qc': (f"{self.STUDY_ID}/table_supplementary_2.xlsx", 0),
        }
        resource_df = self._prepare_resources_for_annotation(annotation_resources)
        resource_df.to_csv(resource_df_filename, index=False)
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
            'qc_total_upper': slugify("library size - upper"),
            # some columns are duplicated in the resource files
            # when they were joined, pandas appended suffices to distinguish, e.g. _x, _y
            'gender': slugify("gender_x"),
            'parity': slugify("parity_x"),
        }
        for index in range(len(resource_df)):
            filename = str(resource_df.loc[index, 'adata-filename'])
            adata = self.objects[filename]
            for uns_name, column_name in annotation_column_names.items():
                adata.uns[uns_name] = resource_df.loc[index, column_name]
            # update the h5ad file with annotations
            self.cache_adata_object(adata, filename)
        self.annotations_loaded()

    def _apply_annotations(self):
        """
        Iterates over the objects (adatas) and uses the annotations to engineer 5 new binary features:

        1. `is_low_num_genes`
            * 1 if the number of genes detected is lower than or equal to a threshold (`qc_genes_lower`)
            * 0 otherwise
        2. `is_high_num_genes`
            * 1 if the number of genes detected is higher than to a threshold (`qc_genes_upper`)
            * 0 otherwise
        3. `is_high_mito`
            * 1 if the percentage (0 to 1) is higher than a threshold (`qc_mito_upper`)
            * 0 otherwise
        4. `is_high_total_count`
            * 1 if the library size is higher than or equal to a threshold (`qc_total_upper`)
            * 0 otherwise
        5. `is_noise`
            * 1 if any of the above are 1
            * 0 otherwise

        :return:
        """
        failed = {}
        success = {}
        for index, adata in enumerate(self.objects.values()):
            filename = adata.uns['adata-filename']
            try:
                L.info(f"Applying annotations for {filename}")
                self._apply_one(adata)
                success[filename] = adata
            except ValueError as value_error:
                L.warning(f"Value error for {filename}: {value_error}")
                failed[index] = value_error

        if len(failed) == 0:
            for filename, adata in success.items():
                self.cache_adata_object(adata, filename)
            self.annotations_applied()
        else:
            L.error(f"Failed to apply annotations at indices {failed}, objects not updated on disk.")

    @staticmethod
    def _apply_one(adata):
        # Annotate mitochondrial genes before getting QC metrics
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')

        # Use scanpy to calculate the QC metrics
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

        # Identify observations that meet the threshold
        adata.obs['is_low_num_genes'] = (adata.obs['n_genes_by_counts'] <= adata.uns['qc_genes_lower']).astype(int)
        adata.obs['is_high_num_genes'] = (adata.obs['n_genes_by_counts'] > adata.uns['qc_genes_upper']).astype(int)
        adata.obs['is_high_mito'] = (adata.obs['pct_counts_mt'] / 100 > adata.uns['qc_mito_upper']).astype(int)
        adata.obs['is_high_total_count'] = (adata.obs['total_counts'] >= adata.uns['qc_total_upper']).astype(int)

        # Additional features while we're here
        adata.obs['zero_genes'] = (adata.obs['n_genes_by_counts'] == 0).count()
        adata.obs['zero_mito'] = (adata.obs['pct_counts_mt'] == 0).count()
        adata.obs['zero_count'] = (adata.obs['total_counts'] == 0).count()

        # Identify observation as noise
        adata.obs['is_noise'] = (
                adata.obs['is_low_num_genes'] |
                adata.obs['is_high_num_genes'] |
                adata.obs['is_high_mito'] |
                adata.obs['is_high_total_count']
        ).astype(int)

        # The number of observations that are not noise should match with the published "after" count
        actual_count = adata[adata.obs['is_noise'] == 0, :].shape[0]
        expected_count = adata.uns['num_cells_after']
        if not bool(actual_count == expected_count) and not (adata.uns['adata-filename'] in GSE161529.EXPECTED_MISMATCHES):
            raise ValueError(f"Check failed! Expected {expected_count} but got {actual_count}.")


    @staticmethod
    def _prepare_resources_for_annotation(resources):
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
        join_column = "sample name"
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
