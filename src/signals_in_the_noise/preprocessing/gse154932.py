from pathlib import Path

import scanpy as sc

from signals_in_the_noise.config import get_data_path
from signals_in_the_noise.io.tenx import TenX
from signals_in_the_noise.preprocessing.base import Preprocessor
from signals_in_the_noise.utils.log import get_logger

logger = get_logger(__name__)


class GSE154932(Preprocessor):
    """Preprocesses dataset GSE154932 for analysis.

    Pre-requisites:
    - Raw compressed data from study has been expanded into ``RAW_DATA_DIRECTORY``
      inside the data directory.
    """

    STUDY_ID = "GSE154932"
    RAW_DATA_DIRECTORY = f"{STUDY_ID}_RAW"

    ## Thresholds from the paper
    # (a) a minimum of 10,000 and maximum of 80,000 transcript counts
    QC_TOTAL_COUNTS_LOWER = 10000
    QC_TOTAL_COUNTS_UPPER = 80000
    # (b) a maximum of 20% of counts attributed to mitochondrial genes
    QC_MITO_UPPER_PCT = 20
    # (c) a minimum of 3,000 genes detected.
    QC_GENES_LOWER = 3000
    # (d) a minimum gene detection in at least 20 cells
    QC_GENE_CELLS_LOWER = 20

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

        if not self.is_annotations_loaded:
            self.annotations_loaded()

        if not self.is_annotations_applied:
            self._apply_annotations()

    def _apply_annotations(self) -> None:
        """Compute and attach QC flag columns to each AnnData object.

        Cell-level thresholds from the paper:
        (a) 10,000--80,000 transcript counts
        (b) maximum 20% mitochondrial counts
        (c) minimum 3,000 genes detected

        ``obs['is_noise']`` is 1 when a cell fails any of (a)--(c).
        ``var['is_noise']`` is 1 when a gene fails (d).
        """
        for adata in self.objects.values():
            filename = adata.obs["adata-filename"].iloc[0]
            logger.info(f"Applying annotations for {filename}")
            self._apply_one(adata)
            # All samples are cancerous
            adata.obs["has_cancer"] = 1
            self.cache_adata_object(adata, filename)

        self.annotations_applied()

    @staticmethod
    def _apply_one(adata) -> None:
        """Annotate a single AnnData object with QC flags (in-place).

        Adds the following ``var`` columns:
        1. ``mt`` — mitochondrial gene indicator
        2. ``is_noise`` — 1 if gene is detected in fewer than 20 cells

        Adds the following ``obs`` columns:
        1. ``is_low_total_count`` — 1 if ``total_counts`` < 10,000
        2. ``is_high_total_count`` — 1 if ``total_counts`` > 80,000
        3. ``is_high_mito`` — 1 if ``pct_counts_mt`` > 20
        4. ``is_low_num_genes`` — 1 if ``n_genes_by_counts`` < 3,000
        5. ``is_noise`` — 1 if any of the above obs flags are set
        """
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        adata.uns["qc_total_counts_lower"] = GSE154932.QC_TOTAL_COUNTS_LOWER
        adata.uns["qc_total_counts_upper"] = GSE154932.QC_TOTAL_COUNTS_UPPER
        adata.uns["qc_mito_upper_pct"] = GSE154932.QC_MITO_UPPER_PCT
        adata.uns["qc_genes_lower"] = GSE154932.QC_GENES_LOWER
        adata.uns["qc_gene_cells_lower"] = GSE154932.QC_GENE_CELLS_LOWER

        adata.obs["is_low_total_count"] = (
            adata.obs["total_counts"] < GSE154932.QC_TOTAL_COUNTS_LOWER
        ).astype(int)
        adata.obs["is_high_total_count"] = (
            adata.obs["total_counts"] > GSE154932.QC_TOTAL_COUNTS_UPPER
        ).astype(int)
        adata.obs["is_high_mito"] = (
            adata.obs["pct_counts_mt"] > GSE154932.QC_MITO_UPPER_PCT
        ).astype(int)
        adata.obs["is_low_num_genes"] = (
            adata.obs["n_genes_by_counts"] < GSE154932.QC_GENES_LOWER
        ).astype(int)

        adata.obs["is_noise"] = (
            adata.obs["is_low_total_count"]
            | adata.obs["is_high_total_count"]
            | adata.obs["is_high_mito"]
            | adata.obs["is_low_num_genes"]
        ).astype(int)

        adata.var["is_noise"] = (
            adata.var["n_cells_by_counts"] < GSE154932.QC_GENE_CELLS_LOWER
        ).astype(int)
