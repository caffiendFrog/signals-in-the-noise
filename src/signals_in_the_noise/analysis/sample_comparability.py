"""Step 0 of the BRCA1 comparative analysis: are WT and BRCA1 samples comparable?

Every comparison treats the donor, not the cell, as the unit of analysis: each
donor is reduced to the median of each QC metric, and genotypes are compared on
those donor medians. This keeps donors with many cells from dominating and
matches the donor-level tests used later in the analysis.
"""

import json
import logging
from collections.abc import Callable, Iterable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

from signals_in_the_noise.analysis.brca1_cohort import Genotype
from signals_in_the_noise.analysis.statistics import difference_in_medians, exact_permutation_test

logger = logging.getLogger(__name__)

QC_METRICS: dict[str, str] = {
    "total_counts": "UMI counts per cell",
    "n_genes_by_counts": "Genes detected per cell",
    "pct_counts_mt": "Mitochondrial reads (%)",
}
"""QC metrics compared in Step 0, mapped to display labels."""

DEPTH_METRIC = "total_counts"
DEFAULT_DEPTH_TOLERANCE = 0.20
"""Relative BRCA1-vs-WT difference in median depth above which depth matching is mandatory."""

CellMask = Callable[[pd.DataFrame], pd.Series]

POPULATIONS: dict[str, CellMask] = {
    "qc_pass": lambda cells: cells["is_noise"] == 0,
    "noise": lambda cells: cells["is_noise"] == 1,
    "all_barcodes": lambda cells: pd.Series(True, index=cells.index),
}
"""Cell populations compared in Step 0.

- ``qc_pass``: cells kept by the paper QC; the basis of Arm B and of clustering.
- ``noise``: cells discarded by the paper QC; the pool PBS cells are rescued from.
- ``all_barcodes``: everything the authors deposited in GEO.
"""


def summarize_donors(
    cells: pd.DataFrame,
    *,
    populations: Mapping[str, CellMask] = POPULATIONS,
    metrics: Iterable[str] = QC_METRICS,
) -> pd.DataFrame:
    """Summarise each QC metric per donor within each cell population.

    Args:
        cells: Cell table from
            :func:`~signals_in_the_noise.analysis.brca1_cohort.build_cell_table`.
        populations: Mapping of population name to a boolean cell mask.
        metrics: ``cells`` columns to summarise.

    Returns:
        Tidy frame with one row per ``(population, donor, metric)`` and columns
        ``population``, ``donor``, ``genotype``, ``metric``, ``n_cells``,
        ``median``, ``q25``, ``q75``.
    """
    metrics = list(metrics)
    rows = []
    for population, mask in populations.items():
        subset = cells.loc[mask(cells)]
        for (donor, genotype), group in subset.groupby(["donor", "genotype"], observed=True):
            for metric in metrics:
                values = group[metric]
                rows.append(
                    {
                        "population": population,
                        "donor": donor,
                        "genotype": genotype,
                        "metric": metric,
                        "n_cells": len(values),
                        "median": values.median(),
                        "q25": values.quantile(0.25),
                        "q75": values.quantile(0.75),
                    }
                )
    return pd.DataFrame(rows)


def compare_genotypes(
    donor_summary: pd.DataFrame,
    *,
    reference: Genotype = Genotype.WT,
    test: Genotype = Genotype.BRCA1,
) -> pd.DataFrame:
    """Compare donor medians between genotypes for every population and metric.

    The genotype-level value is the median of donor medians. Significance comes
    from an exact two-sided permutation test on donor medians.

    Args:
        donor_summary: Output of :func:`summarize_donors`.
        reference: Genotype in the denominator of ``ratio`` (WT).
        test: Genotype in the numerator of ``ratio`` (BRCA1).

    Returns:
        One row per ``(population, metric)`` with the number of donors per genotype,
        ``<reference>_median``, ``<test>_median``, ``ratio`` (test / reference),
        ``relative_difference`` (``ratio - 1``) and ``permutation_p``.
    """
    rows = []
    for (population, metric), group in donor_summary.groupby(
        ["population", "metric"], sort=False
    ):
        ref = group.loc[group["genotype"] == reference, "median"].to_numpy()
        tst = group.loc[group["genotype"] == test, "median"].to_numpy()
        ref_median = float(pd.Series(ref).median())
        tst_median = float(pd.Series(tst).median())
        ratio = tst_median / ref_median if ref_median else float("nan")
        permutation = exact_permutation_test(tst, ref, statistic=difference_in_medians)
        rows.append(
            {
                "population": population,
                "metric": metric,
                f"n_{reference}": len(ref),
                f"n_{test}": len(tst),
                f"{reference}_median": ref_median,
                f"{test}_median": tst_median,
                "ratio": ratio,
                "relative_difference": ratio - 1.0,
                "permutation_p": permutation.p_value,
            }
        )
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class DepthMatchingDecision:
    """Whether the depth-matched control in Step 5 is mandatory.

    Attributes:
        tolerance: Maximum allowed ``|BRCA1 / WT - 1|`` in median depth.
        relative_differences: Population name to BRCA1-vs-WT relative
            difference in median ``total_counts`` (donor medians).
    """

    tolerance: float
    relative_differences: dict[str, float]

    @property
    def exceeding(self) -> list[str]:
        """Populations whose depth difference exceeds the tolerance."""
        return [
            population
            for population, difference in self.relative_differences.items()
            if abs(difference) > self.tolerance
        ]

    @property
    def mandatory(self) -> bool:
        """True when any compared population exceeds the tolerance."""
        return bool(self.exceeding)

    def summary(self) -> str:
        """One-line human-readable explanation of the decision."""
        details = ", ".join(
            f"{population} {difference:+.1%}"
            for population, difference in self.relative_differences.items()
        )
        verdict = "MANDATORY" if self.mandatory else "optional"
        return (
            f"Depth-matched run in Step 5 is {verdict} "
            f"(BRCA1 vs WT median depth: {details}; tolerance ±{self.tolerance:.0%})."
        )

    def to_json(self, path: Path) -> None:
        """Persist the decision so later steps can read it without recomputing Step 0."""
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {**asdict(self), "exceeding": self.exceeding, "mandatory": self.mandatory}
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @classmethod
    def from_json(cls, path: Path) -> "DepthMatchingDecision":
        """Load a decision written by :meth:`to_json`."""
        payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            tolerance=payload["tolerance"],
            relative_differences=payload["relative_differences"],
        )


def assess_depth_matching(
    genotype_comparison: pd.DataFrame,
    *,
    populations: Iterable[str] = ("qc_pass", "noise"),
    tolerance: float = DEFAULT_DEPTH_TOLERANCE,
) -> DepthMatchingDecision:
    """Decide whether depth matching is mandatory from :func:`compare_genotypes` output.

    Both the QC-pass population (which defines clusters) and the noise
    population (which PBS cells are rescued from) are checked by default, since
    a depth imbalance in either could create a genotype difference on its own.

    Args:
        genotype_comparison: Output of :func:`compare_genotypes`.
        populations: Populations to consider.
        tolerance: Maximum allowed relative difference in median depth.
    """
    depth = genotype_comparison.loc[
        genotype_comparison["metric"] == DEPTH_METRIC
    ].set_index("population")
    populations = list(populations)
    missing = [population for population in populations if population not in depth.index]
    if missing:
        raise KeyError(f"No {DEPTH_METRIC} comparison for populations {missing}")

    decision = DepthMatchingDecision(
        tolerance=tolerance,
        relative_differences={
            population: float(depth.loc[population, "relative_difference"])
            for population in populations
        },
    )
    logger.info(decision.summary())
    return decision
