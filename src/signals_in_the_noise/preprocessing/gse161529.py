from collections import defaultdict
from functools import reduce
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from openTSNE import TSNE
from scipy import sparse
from slugify import slugify

from signals_in_the_noise.preprocessing.prep_config import Preprocessor
from signals_in_the_noise.utilities.logging import get_logger
from signals_in_the_noise.utilities.storage import get_data_path, get_resources_path
from signals_in_the_noise.utilities.tenx_genomics import TenX, DirectoryType


L = get_logger(__name__)


class GSE161529(Preprocessor):
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

    EPI_CELL_TYPING_FILENAMES = {
        "N-1105-epi": "GSM4909260_N-N1105-Epi.h5ad",
        "N-0280-epi": "GSM4909255_N-N280-Epi.h5ad",
        "N-0230.16-epi": "GSM4909264_N-N1B-Epi.h5ad",
        "N-0408-epi": "GSM4909259_N-NE-Epi.h5ad",
        "N-1469-epi": "GSM4909258_N-NF-Epi.h5ad",
        "N-0123-epi": "GSM4909267_N-MH0023-Epi.h5ad",
        "N-0064-epi": "GSM4909262_N-MH0064-Epi.h5ad",
        "N-0093-epi": "GSM4909256_N-PM0095-Epi.h5ad",
        "N-0342-epi": "GSM4909269_N-PM0342-Epi.h5ad",
        "N-0372-epi": "GSM4909275_N-PM0372-Epi.h5ad",
        "N-0275-epi": "GSM4909273_N-MH275-Epi.h5ad",
    }

    EPI_CELL_TYPING_GENES = {
        "EPCAM_expr": "EPCAM",
        "KRT5_expr": "KRT5",
        "ACTA2_expr": "ACTA2",
        "MYLK_expr": "MYLK",
        "SNAI2_expr": "SNAI2",
        "NOTCH4_expr": "NOTCH4",
        "DKK3_expr": "DKK3",
        "CD49f_expr": "ITGA6",
        "ESR1_expr": "ESR1",
        "PGR_expr": "PGR",
        "FOXA1_expr": "FOXA1",
        "TNFRSF11A_expr": "TNFRSF11A",
        "KIT_expr": "KIT",
        "SOX10_expr": "SOX10",
    }

    G2M_CHECKPOINT_GENES = {
        "CDC25C_expr": "CDC25C",
        "CDK1_expr": "CDK1",
        "PLK1_expr": "PLK1",
        "AURKA_expr": "AURKA",
        "WEE1_expr": "WEE1",
    }

    E2F_REGULATION_GENES = {
        "CDK1_expr": "CDK1",
        "CCNA2_expr": "CCNA2",
        "E2F1_expr": "E2F1",
        "E2F2_expr": "E2F2",
    }

    DDR_GENES = {
        "BRCA1_expr": "BRCA1",
        "BRCA2_expr": "BRCA2",
        "RAD51_expr": "RAD51",
        "ATM_expr": "ATM",
        "CHEK2_expr": "CHEK2",
        "XRCC1_expr": "XRCC1",
    }

    UPR_GENES = {
        "ATF6_expr": "ATF6",
        "XBP1_expr": "XBP1",
        "CHOP_expr": "DDIT3",
        "GRP78_expr": "HSPA5",
        "PERK_expr": "EIF2AK3",
    }

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

        self.random_kwargs = {
            'use_rep': 'X_pca',
            'random_state': 43,
        }

        for index, adata in enumerate(raw_data.multiple_adata):
            self.objects[adata.obs['adata-filename'].iloc[0]] = adata

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
            filename = adata.obs['adata-filename'].iloc[0]
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
        """
        Applies the following annotations as observations (obs) to the AnnData object (IN PLACE modifications):

        1. 'mt' = list of mitochondrial genes present (strings)
        2. 'is_low_num_genes' = True (1) if 'n_genes_by_counts' is within the given threshold (int)
        3. 'is_high_num_genes' = True (1) if 'n_genes_by_counts' is above the given threshold (int)
        4. 'is_high_mito' = True (1) if 'pct_count_mt' is above the given threshold (int)
        5. 'is_high_total_count' = True (1) if 'total_counts' is above the given threshold (int)
        6. 'is_noise' = True (1) if any of annotatinos 2 - 5 are True (int)
        7. 'zero_genes' = True (1) if 'n_genes_by_count' is 0 (int)
        8. 'zero_mito' = True (1) if 'pct_counts_mt' is 0 (int)
        9. 'zero_count' = True (1) if 'total_counts' is 0 (int)
        10. raw gene expression of selected genes

        The given threshold is provided by a spreadsheet and loaded into the uns field of the AnnData object.

        Note: There are 3 files whose counts for "real" cells did not line up with the material provided by
        the authors. These are whitelisted from the value check at the end.

        :param adata: AnnData object o annotate
        :return: None
        """
        # Annotate mitochondrial genes before getting QC metrics
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
        # adata.obs['mt'] = list(adata.var['mt'])

        # Use scanpy to calculate the QC metrics
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

        # Identify observations that meet the threshold
        adata.obs['is_low_num_genes'] = (adata.obs['n_genes_by_counts'] <= adata.uns['qc_genes_lower']).astype(int)
        adata.obs['is_high_num_genes'] = (adata.obs['n_genes_by_counts'] > adata.uns['qc_genes_upper']).astype(int)
        adata.obs['is_high_mito'] = (adata.obs['pct_counts_mt'] / 100 > adata.uns['qc_mito_upper']).astype(int)
        adata.obs['is_high_total_count'] = (adata.obs['total_counts'] >= adata.uns['qc_total_upper']).astype(int)

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
        if not bool(actual_count == expected_count) and not (adata.obs['adata-filename'].iloc[0] in GSE161529.EXPECTED_MISMATCHES):
            raise ValueError(f"Check failed! Expected {expected_count} but got {actual_count}.")

        # Additional features while we're here
        adata.obs['zero_genes'] = (adata.obs['n_genes_by_counts'] == 0).astype(int)
        adata.obs['zero_mito'] = (adata.obs['pct_counts_mt'] == 0).astype(int)
        adata.obs['zero_count'] = (adata.obs['total_counts'] == 0).astype(int)

        # # Cache raw gene expression for selected genes
        # genes_to_cache = []
        # for selected_genes in (
        #     self.G2M_CHECKPOINT_GENES,
        #     self.E2F_REGULATION_GENES,
        #     self.DDR_GENES,
        #     self.UPR_GENES,
        # ):
        #     genes_to_cache.extend(list(selected_genes.values()))
        # self.cache_raw_gene_expression(adata, genes_to_cache, in_place=True)

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


    def annotate_epithial_cell_typing(self, adata, *, hvg_only=True):
        """
        Annotates the cells for epitihial cell types:
            - basal
            - luminal progenitor
            - mature luminal
            - other (none of the above)
        :param adata: dataset to annotate
        :return: dataset annotated with 4 new observations:
            - score_basal
                positive basal gene signature expression
            - score_lp
                positive luminal progenitor gene signature expression
            - score_ml
                postive mature luminal gene signature expression
            - score_other
                negative basal/lp/ml gene signature expression (may contain nan)
        """
        gene_signature_filenames = {
            'basal': 'epithial_cell_typing/41591_2009_BFnm2000_MOESM13_ESM.xls',
            'lp': 'epithial_cell_typing/41591_2009_BFnm2000_MOESM14_ESM.xls',
            'ml': 'epithial_cell_typing/41591_2009_BFnm2000_MOESM15_ESM.xls',
            'stromal': 'epithial_cell_typing/41591_2009_BFnm2000_MOESM16_ESM.xls',
        }

        # score the dataset for expression of gene signatures
        adata = self.score_gene_signature_expression(
            adata=adata,
            gene_signature_filenames=gene_signature_filenames,
            log_normalize=True,
            hvg_only=hvg_only,
            # reference article specifically mentions seurat
            hvg_flavor='seurat'
        )

        # initial classification using the gene signature with the highest score
        score_cols = [f"score_{k}" for k in gene_signature_filenames.keys()]
        adata.obs['predicted_type'] = adata.obs[score_cols].idxmax(axis=1).str.replace('score_', '')
        adata.obs['predicted_type_score'] = adata.obs[score_cols].max(axis=1)

        # replace any negative scores with `other`, we only care about upregulation
        adata.obs.loc[adata.obs['predicted_type_score'] <= 0, 'predicted_type'] = "other"
        adata.obs['score_other'] = np.where(
            adata.obs['predicted_type'] == 'other',
            adata.obs['predicted_type_score'],
            np.nan
        )

        return adata.copy()

    def get_combined_epithilial_dataset(
            self,
            *,
            real_filename: str = "combined_epi_normal_real.h5ad",
            noise_filename: str = "combined_epi_normal_noise.h5ad",
            use_leiden: bool = True,
            remove_stromal: bool = True,
            hvg_post_stromal: bool = False,
            genes_to_check: list = None,
            real_pca_kwargs: dict = None,
            noise_pca_kwargs: dict = None,
    ):
        """
        Recreates the combined epitihilial dataset used by the authors for their Figure 1 visualizations.

        Two separate paramters provided to customize the clustering done on the "real" and "noise" cells given
        that they have very different data characteristics (number of cells, number of genes, etc.).

        :param use_leiden: True to use the leiden algorithm for clustering. Currently the only option installed.
            To use louvain algorithm, first install it and its dependencies (including cmake). Then modify
            `find_clusters` to remove the forcing to True.
        :param remove_stromal: True to remove stromal typed cells from the combined dataset.
        :param hvg_post_stromal: True to re-filter down to highly variable genes after filtering out stromal cells
            Does nothing if remove_stromal is False.
        :param genes_to_check: Optional list of genes to check for presence of, produces statements in conosle
        :param real_pca_kwargs: kwargs to use for clustering on the "real" cells
        :param noise_pca_kwargs: kwargs to use for clustering on the "noise" cells
        :return:
        """
        all_real_filename = get_data_path(real_filename)
        all_noise_filename = get_data_path(noise_filename)

        all_real_path = Path(all_real_filename)
        all_noise_path = Path(all_noise_filename)

        if all([all_real_path.exists(), all_noise_path.exists()]):
            L.info("Loading combined datasets...")
            adatas_all_real = ad.read_h5ad(all_real_filename)
            adatas_all_noise = ad.read_h5ad(all_noise_filename)
        else:
            L.info("Combining datasets...")
            adatas_real = []
            adatas_noise = []
            for idx, (specimen_id, filename) in enumerate(self.EPI_CELL_TYPING_FILENAMES.items()):
                adata = self.get_dataset(filename)
                adata.obs_names = [f"{filename}_{i}" for i in range(adata.n_obs)]

                for is_noise in (0, 1):
                    adata_subset = adata[adata.obs['is_noise'] == is_noise].copy()
                    if genes_to_check:
                        missing = self.check_adata_for_genes(adata_subset, genes_to_check)
                        L.info(f"[before annotation] Sanity check genes missing ({len(missing)}) {', '.join(missing)}")
                    adata_subset = self.cache_raw_gene_expression(adata_subset, self.EPI_CELL_TYPING_GENES)
                    adata_subset = self.annotate_epithial_cell_typing(adata_subset, hvg_only=False)
                    if remove_stromal:
                        # remove stromal cells - "...removed the stromal subset..."
                        mask = ~adata_subset.obs['predicted_type'].str.lower().str.contains('stromal')
                        adata_subset = adata_subset[mask].copy()
                        if hvg_post_stromal:
                            sc.pp.highly_variable_genes(adata_subset)
                            adata_subset.raw = adata_subset.copy()
                            adata_subset = adata_subset[:, adata_subset.var['highly_variable']].copy()

                        if genes_to_check:
                            missing = self.check_adata_for_genes(adata_subset, genes_to_check)
                            L.info(f"[after stromal filter], sanity check genes missing ({len(missing)}) {', '.join(missing)}")
                    else:
                        adata_subset.raw = adata_subset.copy()

                    # additional features for visualizations
                    adata_subset.obs['specimen_id'] = specimen_id
                    adata_subset.obs['hormonal_status'] = adata_subset.uns['menopause_status']
                    adata_subset.obs['cancer_type'] = adata_subset.uns['cancer_type']
                    adata_subset.obs['cell_population'] = adata_subset.uns['cell_population']
                    adata_subset.obs['gender'] = adata_subset.uns['gender']
                    adata_subset.obs['parity'] = adata_subset.uns['parity']
                    # just in case
                    if sparse.issparse(adata_subset.X):
                        adata_subset.X = adata_subset.X.toarray()

                    if not is_noise:
                        adatas_real.append(adata_subset)
                    else:
                        adatas_noise.append(adata_subset)

            adatas_all_real = ad.concat(adatas_real, join='inner')
            adatas_all_noise = ad.concat(adatas_noise, join='inner')

            if all([
                real_pca_kwargs is not None,
                noise_pca_kwargs is not None,
            ]):
                for adata, pca_kwargs in (
                    (adatas_all_real, real_pca_kwargs),
                    (adatas_all_noise, noise_pca_kwargs),
                ):
                    self.find_clusters(adata, use_leiden=use_leiden, pca_kwargs=pca_kwargs)
                    self.calculate_tsne(adata)

            # write out the combined datasets
            adatas_all_real.write(all_real_filename)
            adatas_all_noise.write(all_noise_filename)

        return adatas_all_real, adatas_all_noise

    def find_clusters(self, adata, *, use_leiden=True, pca_kwargs=None, neighbors_kwargs=None, cluster_kwargs=None):
        """
        Find clusters in the dataset.

        The default parameter values for the clustering algorithm are:
        * resolution = 0.015 (from article)
        * n_neighbors = 15
        * n_comps = 50

        Paramters can be customized by providing one of the kwargs

        Note: n_comps & n_pcs have been explicitly set for determinism.

        :param adata: AnnData object (IN PLACE modification)
        :param use_leiden: True to use Leiden algorithm, False to use Louvain algorithm
        :param pca_kwargs: dictionary of kwargs to pass to the PCA method
        :param neighbors_kwargs: dictionary of kwargs to pass to the neighbors method
        :param cluster_kwargs: dictionary of kwargs to pass to the clustering algorithm
        :return:
        """
        if not use_leiden:
            L.warning('Packages for louvain have not been configured or installed. Switching to leiden.')
            use_leiden = True

        sc.pp.scale(adata)

        default_pca_parameters = {
            'n_comps': 50,
            'random_state': self.random_seed,
        }
        if pca_kwargs:
            default_pca_parameters.update(pca_kwargs)
        sc.pp.pca(adata, **default_pca_parameters)

        default_neighbors_parameters = {
            'n_neighbors': 15,
            'n_pcs': None,
        }
        default_neighbors_parameters.update(self.random_kwargs)
        if neighbors_kwargs:
            default_neighbors_parameters.update(neighbors_kwargs)
        sc.pp.neighbors(adata, **default_neighbors_parameters)

        default_cluster_parameters = {
            'resolution': 0.015,
            'random_state': self.random_seed,
        }
        if cluster_kwargs:
            default_cluster_parameters.update(cluster_kwargs)
        if use_leiden:
            sc.tl.leiden(adata, **default_cluster_parameters)
        else:
            sc.tl.louvain(adata, **default_cluster_parameters)

    def calculate_tsne(self, adata, *, tsne_kwargs=None):
        """
        Calculates low dimensional embeddings of the data.

        Intentionally configuring TSNE manually rather than using scanpy for finer controls.

        :param adata: AnnData object (IN PLACE modification)
        :param tsne_kwargs: dictionary of kwargs to pass to the TSNE
        :return:
        """
        default_tsne_kwargs = {
            'random_state': self.random_seed,
            'initialization': 'pca',
            'n_jobs': 1,  # single thread for determinism
            'n_iter': 1000,  # default value scanpy uses
            'learning_rate': 200,  # default value scanpy uses
        }
        if tsne_kwargs:
            default_tsne_kwargs.update(tsne_kwargs)
        tsne = TSNE(**default_tsne_kwargs)
        # -- for determinism, round the value to guard against floating point noise
        X_pca = adata.obsm['X_pca']
        X_embedding = tsne.fit(np.round(X_pca, decimals=10))
        adata.obsm['X_tsne'] = np.asarray(X_embedding)

    def visualize_tsne(self, adata, color, *, use_raw=False, plot_kwargs=None):
        default_plot_kwargs = {
            'color': color,
            'use_raw': use_raw,
        }
        if plot_kwargs:
            default_plot_kwargs.update(plot_kwargs)
            kwargs = plot_kwargs.copy()

        sc.pl.tsne(adata, **default_plot_kwargs)
