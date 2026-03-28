import logging
from collections import defaultdict

import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)


def classify_noise_subtypes(
    adata: AnnData,
    q_low: float = 0.25,
    q_high: float = 0.75,
) -> AnnData:
    """Classify noise cells into pbs-1, pbs-2, and pbs-3 subtypes.

    Uses per-metric quartile thresholds on mitochondrial read fraction,
    total RNA count, and gene count to assign each cell to a noise-cell
    phenotype.  Three boolean columns are written to ``adata.obs`` in-place
    and the object is returned.

    Subtype definitions (all thresholds derived from the input population):

    - **pbs-1**: high mitochondrial fraction, low total RNA, low gene count.
    - **pbs-2**: low mitochondrial fraction, low total RNA, moderate gene count.
    - **pbs-3**: moderate mitochondrial fraction, moderate total RNA,
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
        ``adata.obs``: ``pbs-1``, ``pbs-2``, and ``pbs-3``.
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

    adata.obs["pbs-1"] = (mask_mito_high & mask_rna_low & mask_genes_low).astype(bool)
    adata.obs["pbs-2"] = (mask_mito_low & mask_rna_low & mask_genes_moderate).astype(bool)
    adata.obs["pbs-3"] = (mask_mito_moderate & mask_rna_moderate & mask_genes_high).astype(bool)

    logger.debug(
        "classified %d noise cells: %d pbs-1, %d pbs-2, %d pbs-3",
        adata.n_obs,
        adata.obs["pbs-1"].sum(),
        adata.obs["pbs-2"].sum(),
        adata.obs["pbs-3"].sum(),
    )

    return adata


def aggregate_noise_subtypes_by_cancer_type(
    adatas: list[AnnData],
    cell_population_filter: str = "Total",
) -> pd.DataFrame:
    """Aggregate noise-subtype counts across specimens and normalise by cancer type.

    For each AnnData in ``adatas`` whose ``uns['cell_population']`` matches
    ``cell_population_filter``:

    1. Subsets to noise cells (``is_noise == 1``) and applies
       :func:`classify_noise_subtypes` so thresholds and labels are derived
       exclusively from the noise population.
    2. Groups specimens by ``uns['cancer_type']``.
    3. Sums pbs-1, pbs-2, pbs-3, noise, and total cell counts.
    4. Normalises biological-signal counts to the total noise count and the
       noise count to the total cell count (expressed as percentages).

    The returned DataFrame index is annotated with specimen counts, e.g.
    ``"Luminal (3 specimens)"``.

    Args:
        adatas: Collection of AnnData objects to process.  Each object must
            have ``uns['cell_population']`` and ``uns['cancer_type']`` set,
            and the ``obs`` columns required by
            :func:`classify_noise_subtypes` (``pct_counts_mt``,
            ``log1p_total_counts``, ``log1p_n_genes_by_counts``) plus
            ``is_noise``.
        cell_population_filter: Value of ``uns['cell_population']`` used to
            select relevant specimens.  Defaults to ``'Total'``.

    Returns:
        DataFrame indexed by annotated cancer-type labels with four columns:
        ``pbs-1``, ``pbs-2``, ``pbs-3`` (each as a percentage of
        total noise cells), and ``noise`` (as a percentage of all cells).
    """
    grouped: dict[str, list[AnnData]] = defaultdict(list)
    for adata in adatas:
        if adata.uns.get("cell_population") != cell_population_filter:
            continue
        noise_adata = adata[adata.obs["is_noise"] == 1].copy()
        classify_noise_subtypes(noise_adata)
        noise_adata.uns["_total_cells"] = adata.shape[0]
        grouped[adata.uns["cancer_type"]].append(noise_adata)

    annotated: dict[str, list[AnnData]] = {
        f"{cancer_type} ({len(specimens)} specimens)": specimens
        for cancer_type, specimens in grouped.items()
    }

    counts: dict[str, tuple] = {}
    for label, specimens in annotated.items():
        counts[label] = (
            sum(a.obs["pbs-1"].sum() for a in specimens),
            sum(a.obs["pbs-2"].sum() for a in specimens),
            sum(a.obs["pbs-3"].sum() for a in specimens),
            sum(a.shape[0] for a in specimens),
            sum(a.uns["_total_cells"] for a in specimens),
        )

    raw_df = pd.DataFrame(counts, index=["pbs-1", "pbs-2", "pbs-3", "noise", "total"]).T

    norm_df = raw_df[["pbs-1", "pbs-2", "pbs-3"]].div(raw_df["noise"], axis=0)
    norm_df["noise"] = raw_df["noise"].div(raw_df["total"], axis=0)
    norm_df = norm_df * 100

    zero_noise_labels = list(raw_df.index[raw_df["noise"] == 0])
    if zero_noise_labels:
        logger.warning(
            "cancer type(s) %s have zero noise cells; subtype percentages set to 0",
            zero_noise_labels,
        )
    norm_df = norm_df.fillna(0.0)

    logger.debug(
        "aggregated noise subtypes for %d cancer type(s) from %d specimens",
        len(annotated),
        sum(len(v) for v in annotated.values()),
    )

    return norm_df
