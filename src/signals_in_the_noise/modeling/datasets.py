"""Cell-level modeling tables built from per-specimen QC frames.

An :class:`Arm` fixes which cells enter the model (discarded, retained, or all)
and how they are featurized (PBS one-hot, raw QC metrics, or within-specimen
QC ranks). Running the same classifier over several arms separates "signal
specific to discarded cells" from "sample-level QC differences".
"""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from signals_in_the_noise.analysis.noise_phenotypes import (
    NO_PBS_LABEL,
    PBS_LABELS,
    PbsThresholds,
    pbs_label,
    pbs_masks,
)
from signals_in_the_noise.preprocessing.qc import QC_METRICS, select_population
from signals_in_the_noise.preprocessing.specimens import (
    CONDITION_COLUMN,
    PATIENT_COLUMN,
    SPECIMEN_COLUMN,
)

CELL_COLUMN = "cell_id"
"""Column holding each cell's original ``obs`` name."""

META_COLUMNS: tuple[str, ...] = (SPECIMEN_COLUMN, PATIENT_COLUMN, CONDITION_COLUMN, CELL_COLUMN)
"""Non-feature columns of a cell table (``patient_id`` is present only when supplied)."""

FeatureBuilder = Callable[[pd.DataFrame], pd.DataFrame]
"""Maps one specimen's selected cells to a feature frame (rows may be dropped)."""


def pbs_features(pbs_thresholds: dict[str, PbsThresholds] | None = None) -> FeatureBuilder:
    """Return a builder of PBS one-hot features, dropping cells with no PBS label.

    Quantile cutoffs are recomputed from whichever cells the builder receives,
    exactly as the published classifier does for each specimen's noise cells.
    """

    def build(cells: pd.DataFrame) -> pd.DataFrame:
        labels = pbs_label(pbs_masks(cells, pbs_thresholds))
        labelled = labels[labels != NO_PBS_LABEL]
        return pd.DataFrame(
            {label: (labelled == label).astype(int) for label in PBS_LABELS}, index=labelled.index
        )

    return build


def qc_features(cells: pd.DataFrame) -> pd.DataFrame:
    """Raw QC metrics; keeps any between-specimen shift in their distributions."""
    return cells[list(QC_METRICS)].astype(float)


def rank_qc_features(cells: pd.DataFrame) -> pd.DataFrame:
    """Within-specimen percentile ranks of the QC metrics.

    Removes every between-specimen difference in the marginal distributions
    while keeping how the metrics co-vary inside the specimen, which is all
    that per-specimen PBS quantiles can see.
    """
    return cells[list(QC_METRICS)].rank(pct=True)


@dataclass(frozen=True)
class Arm:
    """One configuration of which cells to model and how to featurize them.

    Attributes:
        name: Short identifier used in result tables.
        population: ``"noise"``, ``"retained"`` or ``"all"``.
        features: Feature builder applied per specimen.
        size_match: Subsample the population to the specimen's noise-cell count
            before featurizing, so arms see equally many cells.
        description: Human-readable summary for reports.
    """

    name: str
    population: str
    features: FeatureBuilder
    size_match: bool = False
    description: str = ""


CORE_ARM_NAMES: tuple[str, ...] = (
    "noise_pbs",
    "retained_pbs_matched",
    "noise_qc_rank",
    "retained_qc",
)
"""Smallest arm set that still separates the competing explanations."""


def standard_arms(
    pbs_thresholds: dict[str, PbsThresholds] | None = None, names: Iterable[str] | None = None
) -> list[Arm]:
    """The arms used to test whether the PBS signal is specific to discarded cells.

    The first arm reproduces the published classifier and is the reference.

    Args:
        pbs_thresholds: PBS thresholds for the PBS arms.
        names: Optional subset of arm names to keep (standard order is preserved).

    Raises:
        ValueError: If ``names`` contains an unknown arm.
    """
    arms = _all_standard_arms(pbs_thresholds)
    if names is None:
        return arms
    wanted = set(names)
    unknown = wanted - {arm.name for arm in arms}
    if unknown:
        raise ValueError(f"Unknown arm name(s): {sorted(unknown)}.")
    return [arm for arm in arms if arm.name in wanted]


def _all_standard_arms(pbs_thresholds: dict[str, PbsThresholds] | None) -> list[Arm]:
    pbs = pbs_features(pbs_thresholds)
    return [
        Arm("noise_pbs", "noise", pbs, description="PBS one-hot of discarded cells (published)"),
        Arm("retained_pbs", "retained", pbs, description="PBS one-hot of retained cells"),
        Arm(
            "retained_pbs_matched",
            "retained",
            pbs,
            size_match=True,
            description="PBS one-hot of retained cells, subsampled to the noise-cell count",
        ),
        Arm("all_pbs", "all", pbs, description="PBS one-hot of all cells"),
        Arm("noise_qc", "noise", qc_features, description="Raw QC metrics of discarded cells"),
        Arm(
            "noise_qc_rank",
            "noise",
            rank_qc_features,
            description="Within-specimen QC ranks of discarded cells",
        ),
        Arm("retained_qc", "retained", qc_features, description="Raw QC metrics of retained cells"),
    ]


def build_cell_table(
    specimen_obs: Mapping[str, pd.DataFrame],
    conditions: Mapping[str, str],
    arm: Arm,
    *,
    patients: Mapping[str, str] | None = None,
    noise_column: str = "is_noise",
    seed: int = 0,
) -> pd.DataFrame:
    """Featurize every specimen under ``arm`` and stack the results.

    Args:
        specimen_obs: Specimen identifier to cell-level QC frame.
        conditions: Specimen identifier to condition label.
        arm: Which cells to use and how to featurize them.
        patients: Optional specimen identifier to patient identifier; adds a
            ``patient_id`` column so folds can hold out whole patients.
        noise_column: ``obs`` column defining noise vs retained cells.
        seed: Seed for size-matched subsampling.

    Returns:
        DataFrame with the metadata columns of :data:`META_COLUMNS` followed
        by the arm's feature columns, one row per cell, with a fresh ``RangeIndex``.
    """
    rng = np.random.default_rng(seed)
    frames = []
    for specimen_id in sorted(specimen_obs):
        obs = specimen_obs[specimen_id]
        cells = select_population(obs, arm.population, noise_column)
        if arm.size_match:
            n_noise = int((obs[noise_column] == 1).sum())
            if len(cells) > n_noise:
                keep = np.sort(rng.choice(len(cells), size=n_noise, replace=False))
                cells = cells.iloc[keep]
        features = arm.features(cells)
        meta = {SPECIMEN_COLUMN: specimen_id}
        if patients is not None:
            meta[PATIENT_COLUMN] = patients[specimen_id]
        meta[CONDITION_COLUMN] = conditions[specimen_id]
        meta[CELL_COLUMN] = features.index.astype(str)
        frames.append(pd.concat([pd.DataFrame(meta, index=features.index), features], axis=1))
    if not frames:
        return pd.DataFrame(
            columns=[c for c in META_COLUMNS if patients is not None or c != PATIENT_COLUMN]
        )
    return pd.concat(frames, ignore_index=True)


def feature_columns(table: pd.DataFrame) -> list[str]:
    """Columns of a cell or specimen table that are not metadata."""
    return [c for c in table.columns if c not in META_COLUMNS]


def specimen_sizes(table: pd.DataFrame) -> pd.Series:
    """Number of rows (cells) contributed by each specimen."""
    return table.groupby(SPECIMEN_COLUMN).size()


def specimen_composition(table: pd.DataFrame) -> pd.DataFrame:
    """Per-specimen mean of each feature column, plus ``condition`` (and ``patient_id`` if present).

    For PBS one-hot features this is the fraction of labelled cells in each
    subtype, which fully determines what a cell-level classifier on those
    features predicts for the specimen.
    """
    grouped = table.groupby(SPECIMEN_COLUMN)
    composition = grouped[feature_columns(table)].mean()
    meta = [c for c in (PATIENT_COLUMN, CONDITION_COLUMN) if c in table.columns]
    return pd.concat([grouped[meta].first(), composition], axis=1)
