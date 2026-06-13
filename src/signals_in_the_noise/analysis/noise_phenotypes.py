import logging
from collections import defaultdict
from dataclasses import dataclass

import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Thresholds:
    """Quantile boundaries for low, high, and moderate classification of one metric."""

    q_low: float
    q_high: float
    q_mod_low: float
    q_mod_high: float


@dataclass(frozen=True)
class PbsThresholds:
    """Per-metric quantile thresholds used to classify one PBS subtype."""

    pct_counts_mt: Thresholds
    log1p_total_counts: Thresholds
    log1p_n_genes_by_counts: Thresholds


def _thresh(q_low: float, q_high: float) -> Thresholds:
    return Thresholds(q_low, q_high, q_low, q_high)


_DEFAULT = _thresh(0.25, 0.75)


DEFAULT_PBS_THRESHOLDS: dict[str, PbsThresholds] = {
    "pbs-1": PbsThresholds(_DEFAULT, _DEFAULT, _DEFAULT),
    "pbs-2": PbsThresholds(_DEFAULT, _DEFAULT, _DEFAULT),
    "pbs-3": PbsThresholds(_DEFAULT, _DEFAULT, _DEFAULT),
}
"""Default per-PBS, per-metric quantile thresholds for noise subtype classification."""


def classify_noise_subtypes(
    adata: AnnData, *,
    pbs_thresholds: dict[str, PbsThresholds] | None = None,
) -> AnnData:
    """Classify noise cells into pbs-1, pbs-2, and pbs-3 subtypes.

    Assigns each cell to a noise-cell phenotype using per-metric quantile
    thresholds on mitochondrial read fraction, total RNA count, and gene
    count. Three boolean columns are written to ``adata.obs`` in place and
    the object is returned.

    Subtype definitions (quantile cutoffs are derived from the input
    population):

    - **pbs-1**: high mitochondrial fraction, low total RNA, low gene count.
    - **pbs-2**: low mitochondrial fraction, low total RNA, moderate gene count.
    - **pbs-3**: moderate mitochondrial fraction, moderate total RNA,
      high gene count.

    Default per-PBS, per-metric quantiles are defined in
    :data:`DEFAULT_PBS_THRESHOLDS`.

    Args:
        adata: AnnData object containing noise cells. Must include
            ``pct_counts_mt``, ``log1p_total_counts``, and
            ``log1p_n_genes_by_counts`` in ``adata.obs``.
        pbs_thresholds: Optional mapping from PBS label (``"pbs-1"``,
            ``"pbs-2"``, or ``"pbs-3"``) to :class:`PbsThresholds`. When
            omitted, :data:`DEFAULT_PBS_THRESHOLDS` is used.

    Returns:
        The same ``adata`` object with ``pbs-1``, ``pbs-2``, and ``pbs-3``
        boolean columns added to ``adata.obs``.
    """
    obs = adata.obs

    thresholds = DEFAULT_PBS_THRESHOLDS if pbs_thresholds is None else pbs_thresholds

    t1 = thresholds["pbs-1"]
    t2 = thresholds["pbs-2"]
    t3 = thresholds["pbs-3"]

    mito = obs["pct_counts_mt"]
    rna = obs["log1p_total_counts"]
    genes = obs["log1p_n_genes_by_counts"]

    adata.obs["pbs-1"] = (
        (mito >= mito.quantile(t1.pct_counts_mt.q_high))
        & (rna <= rna.quantile(t1.log1p_total_counts.q_low))
        & (genes <= genes.quantile(t1.log1p_n_genes_by_counts.q_low))
    ).astype(bool)

    adata.obs["pbs-2"] = (
        (mito <= mito.quantile(t2.pct_counts_mt.q_low))
        & (rna <= rna.quantile(t2.log1p_total_counts.q_low))
        & (genes > genes.quantile(t2.log1p_n_genes_by_counts.q_mod_low))
        & (genes < genes.quantile(t2.log1p_n_genes_by_counts.q_mod_high))
    ).astype(bool)

    adata.obs["pbs-3"] = (
        (mito > mito.quantile(t3.pct_counts_mt.q_mod_low))
        & (mito < mito.quantile(t3.pct_counts_mt.q_mod_high))
        & (rna > rna.quantile(t3.log1p_total_counts.q_mod_low))
        & (rna < rna.quantile(t3.log1p_total_counts.q_mod_high))
        & (genes >= genes.quantile(t3.log1p_n_genes_by_counts.q_high))
    ).astype(bool)

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
    cell_population_filter: str = "Total", *,
    pbs_thresholds: dict[str, PbsThresholds] | None = None,
) -> pd.DataFrame:
    """Aggregate noise-subtype counts across specimens and normalize by cancer type.

    For each AnnData in ``adatas`` whose ``uns['cell_population']`` matches
    ``cell_population_filter``:

    1. Subsets to noise cells (``is_noise == 1``) and applies
       :func:`classify_noise_subtypes` so thresholds and labels are derived
       exclusively from the noise population.
    2. Groups specimens by ``uns['cancer_type']``.
    3. Sums pbs-1, pbs-2, pbs-3, noise, and total cell counts.
    4. Normalizes biological-signal counts to the total noise count and the
       noise count to the total cell count (expressed as percentages).

    The returned DataFrame index is annotated with specimen counts, for
    example ``"Luminal (3 specimens)"``.

    Args:
        adatas: AnnData objects to process. Each object must set
            ``uns['cell_population']`` and ``uns['cancer_type']``, and its
            ``obs`` must include the columns required by
            :func:`classify_noise_subtypes` plus ``is_noise``.
        cell_population_filter: Value of ``uns['cell_population']`` used to
            select specimens. Defaults to ``"Total"``.
        pbs_thresholds: Optional mapping from PBS label to
            :class:`PbsThresholds`, passed through to
            :func:`classify_noise_subtypes`. When omitted,
            :data:`DEFAULT_PBS_THRESHOLDS` is used.

    Returns:
        DataFrame indexed by annotated cancer-type labels with columns
        ``pbs-1``, ``pbs-2``, ``pbs-3`` (each as a percentage of total noise
        cells), and ``noise`` (as a percentage of all cells).
    """
    thresholds = DEFAULT_PBS_THRESHOLDS if pbs_thresholds is None else pbs_thresholds

    grouped: dict[str, list[AnnData]] = defaultdict(list)
    for adata in adatas:
        if adata.uns.get("cell_population") != cell_population_filter:
            continue
        noise_adata = adata[adata.obs["is_noise"] == 1].copy()
        classify_noise_subtypes(noise_adata, pbs_thresholds=thresholds)
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
