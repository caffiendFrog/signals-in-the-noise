"""Donor registry and per-cell QC table for the BRCA1 pre-neoplastic comparison.

The cohort is the one analysed in Pal et al. (2021) Fig 4A-C / Appendix S1A-C
(author script ``NormBRCA1.R``): 8 normal premenopausal "Total" samples (WT) and
4 BRCA1 pre-neoplastic "Total" samples. Donor order matches the author script.
"""

import logging
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)


class Genotype(StrEnum):
    """Germline BRCA1 status of a donor."""

    WT = "WT"
    BRCA1 = "BRCA1"


AUTHOR_GENOTYPE_LABELS: dict[Genotype, str] = {Genotype.WT: "Normal", Genotype.BRCA1: "BRCA1"}
"""Labels the authors use for each genotype in ``NormBRCA1.R`` figures."""

GENOTYPE_PALETTE: dict[str, str] = {Genotype.WT: "#00008B", Genotype.BRCA1: "#EEC900"}
"""Author colours (R ``darkblue`` / ``gold2``) so figures read the same as Appendix S1A."""


@dataclass(frozen=True)
class Donor:
    """One sample of the cohort.

    Attributes:
        author_id: Sample name used by the authors and in the annotation table.
        geo_stem: GEO file stem (``<GSM>_<name>``) of the raw 10X files.
        genotype: Germline BRCA1 status.
    """

    author_id: str
    geo_stem: str
    genotype: Genotype

    @property
    def h5ad_filename(self) -> str:
        """Filename of the cached AnnData produced by the GSE161529 preprocessor."""
        return f"{self.geo_stem}.h5ad"


DONORS: tuple[Donor, ...] = (
    Donor("N-0019-total", "GSM4909254_N-PM0019-Total", Genotype.WT),
    Donor("N-0233-total", "GSM4909265_N-PM0233-Total", Genotype.WT),
    Donor("N-0092-total", "GSM4909253_N-PM0092-Total", Genotype.WT),
    Donor("N-0230.17-total", "GSM4909261_N-PM0230-Total", Genotype.WT),
    Donor("N-0093-total", "GSM4909257_N-PM0095-Total", Genotype.WT),
    Donor("N-0123-total", "GSM4909268_N-MH0023-Total", Genotype.WT),
    Donor("N-0064-total", "GSM4909263_N-MH0064-Total", Genotype.WT),
    Donor("N-0169-total", "GSM4909266_N-MH0169-Total", Genotype.WT),
    Donor("B1-0894", "GSM4909277_B1-KCF0894", Genotype.BRCA1),
    Donor("B1-0033", "GSM4909278_B1-MH0033", Genotype.BRCA1),
    Donor("B1-0023", "GSM4909279_B1-MH0023", Genotype.BRCA1),
    Donor("B1-0090", "GSM4909280_B1-MH0090", Genotype.BRCA1),
)

PAPER_QC_PASS_TOTALS: dict[Genotype, int] = {Genotype.WT: 36_526, Genotype.BRCA1: 23_240}
"""Cells after QC reported by the authors for this cohort (59,766 in total)."""

QC_METRIC_COLUMNS: tuple[str, ...] = ("total_counts", "n_genes_by_counts", "pct_counts_mt")
"""Per-cell QC metrics computed by ``scanpy.pp.calculate_qc_metrics`` on raw counts."""

CELL_TABLE_COLUMNS: tuple[str, ...] = (
    "donor", "genotype", "barcode", *QC_METRIC_COLUMNS, "is_noise"
)


def build_cell_table(
    get_dataset: Callable[[str], AnnData],
    donors: Iterable[Donor] = DONORS,
) -> pd.DataFrame:
    """Collect per-cell QC metrics and the paper QC flag for every barcode of every donor.

    Args:
        get_dataset: Returns the annotated AnnData for an h5ad filename, e.g.
            ``GSE161529().get_dataset``. ``obs`` must contain
            :data:`QC_METRIC_COLUMNS` and ``is_noise``.
        donors: Donors to include. Defaults to :data:`DONORS`.

    Returns:
        One row per barcode with columns :data:`CELL_TABLE_COLUMNS`. ``donor``
        and ``genotype`` are ordered categoricals.

    Raises:
        ValueError: If a donor's AnnData is empty.
        KeyError: If a required ``obs`` column is missing.
    """
    donors = tuple(donors)
    required = [*QC_METRIC_COLUMNS, "is_noise"]
    frames = []
    for donor in donors:
        adata = get_dataset(donor.h5ad_filename)
        if adata.n_obs == 0:
            raise ValueError(f"Empty AnnData for {donor.author_id} ({donor.h5ad_filename})")
        missing = [column for column in required if column not in adata.obs]
        if missing:
            raise KeyError(f"{donor.author_id} is missing obs columns {missing}")

        frame = adata.obs.loc[:, required].reset_index(drop=True)
        frame.insert(0, "barcode", adata.obs_names.astype(str))
        frame.insert(0, "genotype", donor.genotype.value)
        frame.insert(0, "donor", donor.author_id)
        frames.append(frame)
        logger.info("collected %d barcodes for %s", adata.n_obs, donor.author_id)

    return _with_cohort_dtypes(pd.concat(frames, ignore_index=True), donors)


def load_or_build_cell_table(
    path: Path,
    build: Callable[[], pd.DataFrame],
    *,
    donors: Iterable[Donor] = DONORS,
    force: bool = False,
) -> pd.DataFrame:
    """Return the cached cell table at ``path``, building and caching it when absent.

    Args:
        path: CSV cache location (``.csv.gz`` is compressed automatically).
        build: Zero-argument callable producing the table; only called on a cache miss.
        donors: Donors defining the categorical order of ``donor``.
        force: Rebuild even when the cache exists.
    """
    if path.exists() and not force:
        logger.info("loading cached cell table from %s", path)
        return _with_cohort_dtypes(pd.read_csv(path), tuple(donors))

    table = build()
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False)
    logger.info("wrote cell table (%d rows) to %s", len(table), path)
    return table


def compare_qc_pass_counts_to_paper(
    cell_table: pd.DataFrame, annotations: pd.DataFrame
) -> pd.DataFrame:
    """Compare QC-pass cells per donor with ``number-of-cells-after-filtering`` from the paper.

    Args:
        cell_table: Output of :func:`build_cell_table`.
        annotations: ``GSE161529_annotations_df.csv`` (columns ``sample-name`` and
            ``number-of-cells-after-filtering``).

    Returns:
        One row per donor with ``genotype``, ``qc_pass_observed``, ``qc_pass_paper``
        and ``matches``.
    """
    observed = (
        cell_table.loc[cell_table["is_noise"] == 0]
        .groupby(["donor", "genotype"], observed=True)
        .size()
        .rename("qc_pass_observed")
        .reset_index()
    )
    paper = annotations.set_index("sample-name")["number-of-cells-after-filtering"]
    observed["qc_pass_paper"] = observed["donor"].astype(str).map(paper).astype(int)
    observed["matches"] = observed["qc_pass_observed"] == observed["qc_pass_paper"]
    return observed


def compare_qc_pass_totals_to_paper(cell_table: pd.DataFrame) -> pd.DataFrame:
    """Compare QC-pass cells per genotype with the totals reported in the paper text."""
    observed = (
        cell_table.loc[cell_table["is_noise"] == 0].groupby("genotype", observed=True).size()
    )
    return pd.DataFrame(
        {
            "qc_pass_observed": observed,
            "qc_pass_paper": pd.Series({g.value: n for g, n in PAPER_QC_PASS_TOTALS.items()}),
        }
    ).assign(matches=lambda df: df["qc_pass_observed"] == df["qc_pass_paper"])


def _with_cohort_dtypes(table: pd.DataFrame, donors: tuple[Donor, ...]) -> pd.DataFrame:
    table = table.loc[:, list(CELL_TABLE_COLUMNS)].copy()
    table["donor"] = pd.Categorical(
        table["donor"], categories=[d.author_id for d in donors], ordered=True
    )
    table["genotype"] = pd.Categorical(
        table["genotype"], categories=[g.value for g in Genotype], ordered=True
    )
    table["barcode"] = table["barcode"].astype(str)
    table["is_noise"] = table["is_noise"].astype(int)
    return table
