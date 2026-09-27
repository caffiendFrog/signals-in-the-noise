"""Cell classes, analysis arms and per-donor LP fractions for the BRCA1 comparison.

Every barcode gets one ``cell_class``: ``qc_pass`` for cells kept by the paper
QC, a PBS label for noise cells matching a PBS class, or ``unclassified_noise``.
PBS labels use quantile cutoffs frozen on pooled WT noise cells, so all donors
are classified on the same scale and BRCA1 cells cannot move the cutoffs.

An arm is the set of cell classes it keeps. LP fraction is computed per donor
within epithelium, LP / (LP + basal + ML), so the donor is the unit of analysis.
"""

import logging

import numpy as np
import pandas as pd

from signals_in_the_noise.analysis.brca1_cohort import Genotype
from signals_in_the_noise.analysis.noise_phenotypes import (
    PbsThresholds,
    classify_noise_subtypes_frozen,
)

logger = logging.getLogger(__name__)

QC_PASS = "qc_pass"
UNCLASSIFIED_NOISE = "unclassified_noise"

ARMS: dict[str, tuple[str, ...]] = {
    "B": (QC_PASS,),
    "D": (QC_PASS, "pbs-2", "pbs-3"),
}
"""Cell classes kept by each arm.

- ``B``: standard QC (the paper's cells).
- ``D``: PBS-1 excluded, PBS-2 and PBS-3 retained.

Arm C (ambient- and doublet-corrected) changes counts rather than selecting
cell classes, so it is not expressible here.
"""

EPITHELIAL_TYPES: tuple[str, ...] = ("basal", "lp", "ml")
LP_TYPE = "lp"


def add_log1p_qc_metrics(cells: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``cells`` with the log1p QC columns the PBS thresholds use."""
    cells = cells.copy()
    cells["log1p_total_counts"] = np.log1p(cells["total_counts"])
    cells["log1p_n_genes_by_counts"] = np.log1p(cells["n_genes_by_counts"])
    return cells


def assign_cell_classes(
    cells: pd.DataFrame,
    *,
    reference_genotype: Genotype = Genotype.WT,
    pbs_thresholds: dict[str, PbsThresholds] | None = None,
) -> pd.DataFrame:
    """Label every barcode with its cell class, using PBS cutoffs frozen on reference noise cells.

    Args:
        cells: Cell table from
            :func:`~signals_in_the_noise.analysis.brca1_cohort.build_cell_table`.
        reference_genotype: Genotype whose pooled noise cells define the PBS cutoffs.
        pbs_thresholds: Passed to
            :func:`~signals_in_the_noise.analysis.noise_phenotypes.classify_noise_subtypes_frozen`.

    Returns:
        Copy of ``cells`` with log1p QC columns and a ``cell_class`` column.

    Raises:
        ValueError: If there are no reference noise cells, or a cell matches
            more than one PBS class.
    """
    cells = add_log1p_qc_metrics(cells)
    noise = cells["is_noise"] == 1
    reference = cells.loc[noise & (cells["genotype"] == reference_genotype)]
    if reference.empty:
        raise ValueError(f"No {reference_genotype} noise cells to freeze PBS cutoffs on.")

    labels = classify_noise_subtypes_frozen(
        cells.loc[noise], reference, pbs_thresholds=pbs_thresholds
    )
    overlapping = int((labels.sum(axis=1) > 1).sum())
    if overlapping:
        raise ValueError(f"{overlapping} cells match more than one PBS class.")

    cell_class = pd.Series(UNCLASSIFIED_NOISE, index=cells.index, dtype=object)
    cell_class[~noise] = QC_PASS
    for label in labels.columns:
        cell_class[labels.index[labels[label]]] = label
    cells["cell_class"] = cell_class
    logger.info("cell classes: %s", cell_class.value_counts().to_dict())
    return cells


def rescued_classes(test_arm: str, reference_arm: str) -> tuple[str, ...]:
    """Cell classes ``test_arm`` adds on top of ``reference_arm``.

    Raises:
        ValueError: If ``reference_arm`` is not contained in ``test_arm``.
    """
    reference, test = ARMS[reference_arm], ARMS[test_arm]
    if not set(reference) <= set(test):
        raise ValueError(f"Arm {reference_arm} is not a subset of arm {test_arm}.")
    return tuple(c for c in test if c not in reference)


def delta_by_donor(
    cells: pd.DataFrame,
    *,
    reference_arm: str = "B",
    test_arm: str = "D",
    cell_type_column: str = "cell_type",
) -> pd.DataFrame:
    """Per-donor LP fraction within epithelium under two arms, and their difference Δ.

    With ``N`` epithelial and ``L`` LP cells in the reference arm, and ``M``
    epithelial and ``K`` LP cells among the cells the test arm adds::

        Δ = (L + K) / (N + M) - L / N = w * (r - p)

    where ``p = L / N``, ``r = K / M`` and ``w = M / (N + M)`` is the rescued
    share of epithelium. Whatever the rescued cells are, ``|Δ| <= w``.

    Args:
        cells: Output of :func:`assign_cell_classes` restricted to cells with an
            identity call, plus ``cell_type_column``.
        reference_arm: Arm in :data:`ARMS` Δ is measured from.
        test_arm: Arm in :data:`ARMS` that adds rescued cells.
        cell_type_column: Column with ``basal``/``lp``/``ml`` calls; anything
            else counts as non-epithelial.

    Returns:
        One row per donor with ``genotype``; ``n_epithelial_<arm>`` and
        ``lp_fraction_<arm>`` for both arms; ``n_epithelial_rescued``,
        ``rescued_lp_fraction`` (``r``) and ``rescued_share`` (``w``);
        ``share_<class>`` for each rescued class; and ``delta``.
    """
    added = rescued_classes(test_arm, reference_arm)
    cell_class = cells["cell_class"]
    epithelial = cells[cell_type_column].isin(EPITHELIAL_TYPES)
    lp = cells[cell_type_column] == LP_TYPE
    in_reference = cell_class.isin(ARMS[reference_arm])
    rescued = cell_class.isin(added)

    indicators = pd.DataFrame(
        {
            "epi_reference": epithelial & in_reference,
            "lp_reference": lp & in_reference,
            "epi_rescued": epithelial & rescued,
            "lp_rescued": lp & rescued,
            **{f"epi_{c}": epithelial & (cell_class == c) for c in added},
        }
    )
    indicators["donor"] = cells["donor"]
    indicators["genotype"] = cells["genotype"]
    sums = indicators.groupby(["donor", "genotype"], observed=True).sum()

    n_reference = sums["epi_reference"]
    n_rescued = sums["epi_rescued"]
    n_test = n_reference + n_rescued
    p = sums["lp_reference"] / n_reference.replace(0, np.nan)
    f_test = (sums["lp_reference"] + sums["lp_rescued"]) / n_test.replace(0, np.nan)

    out = pd.DataFrame(
        {
            f"n_epithelial_{reference_arm}": n_reference,
            f"lp_fraction_{reference_arm}": p,
            "n_epithelial_rescued": n_rescued,
            "rescued_lp_fraction": sums["lp_rescued"] / n_rescued.replace(0, np.nan),
            "rescued_share": n_rescued / n_test.replace(0, np.nan),
            **{f"share_{c}": sums[f"epi_{c}"] / n_test.replace(0, np.nan) for c in added},
            f"n_epithelial_{test_arm}": n_test,
            f"lp_fraction_{test_arm}": f_test,
            "delta": f_test - p,
        }
    )
    return out.reset_index()


def rescued_share_by_donor(
    cells: pd.DataFrame, *, reference_arm: str = "B", test_arm: str = "D"
) -> pd.DataFrame:
    """Identity-free rescued share per donor: rescued cells / all cells in the test arm.

    Equals the rescued share of epithelium when rescued cells are epithelial as
    often as reference-arm cells, so it is a check on :func:`delta_by_donor`
    that does not depend on cell identity.

    Returns:
        One row per donor with ``genotype``, ``n_<class>`` for every class in
        the test arm, ``share_<class>`` for each rescued class and ``rescued_share``.
    """
    added = rescued_classes(test_arm, reference_arm)
    counts = (
        cells.loc[cells["cell_class"].isin(ARMS[test_arm])]
        .groupby(["donor", "genotype", "cell_class"], observed=True)
        .size()
        .unstack("cell_class", fill_value=0)
        .reindex(columns=list(ARMS[test_arm]), fill_value=0)
    )
    total = counts.sum(axis=1)
    out = counts.add_prefix("n_")
    for c in added:
        out[f"share_{c}"] = counts[c] / total
    out["rescued_share"] = counts[list(added)].sum(axis=1) / total
    return out.reset_index()
