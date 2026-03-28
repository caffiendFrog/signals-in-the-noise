import logging

from anndata import AnnData

logger = logging.getLogger(__name__)


def classify_noise_subtypes(
    adata: AnnData,
    q_low: float = 0.25,
    q_high: float = 0.75,
) -> AnnData:
    """Classify noise cells into damaged, dormant, and multifunction subtypes.

    Uses per-metric quartile thresholds on mitochondrial read fraction,
    total RNA count, and gene count to assign each cell to a noise-cell
    phenotype.  Three boolean columns are written to ``adata.obs`` in-place
    and the object is returned.

    Subtype definitions (all thresholds derived from the input population):

    - **damaged**: high mitochondrial fraction, low total RNA, low gene count.
    - **dormant**: low mitochondrial fraction, low total RNA, moderate gene count.
    - **multifunction**: moderate mitochondrial fraction, moderate total RNA,
      high gene count.

    Args:
        adata: AnnData object containing noise cells.  Must have the following
            columns in ``adata.obs``: ``pct_counts_mt``,
            ``log1p_total_counts``, ``log1p_n_genes_by_counts``.
        q_low: Lower quantile used as the "low" boundary for each metric.
            Defaults to 0.25 (25th percentile).
        q_high: Upper quantile used as the "high" boundary for each metric.
            Defaults to 0.75 (75th percentile).

    Returns:
        The same ``adata`` object with three new boolean columns appended to
        ``adata.obs``: ``damaged``, ``dormant``, and ``multifunction``.
    """
    obs = adata.obs

    mito_lo = obs["pct_counts_mt"].quantile(q_low)
    mito_high = obs["pct_counts_mt"].quantile(q_high)

    rna_lo = obs["log1p_total_counts"].quantile(q_low)
    rna_high = obs["log1p_total_counts"].quantile(q_high)

    gene_lo = obs["log1p_n_genes_by_counts"].quantile(q_low)
    gene_high = obs["log1p_n_genes_by_counts"].quantile(q_high)

    mask_mito_high = obs["pct_counts_mt"] >= mito_high
    mask_mito_low = obs["pct_counts_mt"] <= mito_lo
    mask_mito_moderate = (obs["pct_counts_mt"] > mito_lo) & (obs["pct_counts_mt"] < mito_high)

    mask_rna_low = obs["log1p_total_counts"] <= rna_lo
    mask_rna_moderate = (obs["log1p_total_counts"] > rna_lo) & (obs["log1p_total_counts"] < rna_high)

    mask_genes_low = obs["log1p_n_genes_by_counts"] <= gene_lo
    mask_genes_high = obs["log1p_n_genes_by_counts"] >= gene_high
    mask_genes_moderate = (obs["log1p_n_genes_by_counts"] > gene_lo) & (obs["log1p_n_genes_by_counts"] < gene_high)

    adata.obs["damaged"] = (mask_mito_high & mask_rna_low & mask_genes_low).astype(bool)
    adata.obs["dormant"] = (mask_mito_low & mask_rna_low & mask_genes_moderate).astype(bool)
    adata.obs["multifunction"] = (mask_mito_moderate & mask_rna_moderate & mask_genes_high).astype(bool)

    logger.debug(
        "classified %d noise cells: %d damaged, %d dormant, %d multifunction",
        adata.n_obs,
        adata.obs["damaged"].sum(),
        adata.obs["dormant"].sum(),
        adata.obs["multifunction"].sum(),
    )

    return adata
