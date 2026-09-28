"""Select specimens from a study and identify the patients they came from."""

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from anndata import AnnData

CONDITION_COLUMN = "condition"
"""Column holding each specimen's condition label (``uns['cancer_type']``)."""

SPECIMEN_COLUMN = "specimen_id"
"""Column holding each specimen's unique identifier."""

PATIENT_COLUMN = "patient_id"
"""Column holding the patient a specimen came from."""

_PATIENT_PATTERN = re.compile(r"Patient\s+(\S+)\s*$")


def filename_specimen_id(adata: AnnData) -> str:
    """Unique specimen identifier from the source file, e.g. ``GSM4909313_ER-MH0064-T``.

    GSE161529 titles are not unique (patient 0029 has two ER+ Total samples
    with the same title), so the GEO file name is used instead.
    """
    return Path(str(adata.obs["adata-filename"].iloc[0])).stem


def patient_id(adata: AnnData) -> str:
    """Patient identifier parsed from a GSE161529 title.

    For example ``"ER+ tumour Total cells from Patient 0029"`` gives ``"0029"``.

    Raises:
        ValueError: If the title does not end with ``Patient <id>``.
    """
    title = str(adata.uns["title"])
    match = _PATIENT_PATTERN.search(title)
    if match is None:
        raise ValueError(f"No patient identifier in title {title!r}.")
    return match.group(1)


def select_specimens(
    adatas: Iterable[AnnData],
    conditions: Iterable[str],
    *,
    cell_population: str = "Total",
) -> dict[str, AnnData]:
    """Return specimens of the requested conditions keyed by :func:`filename_specimen_id`.

    Objects are returned by reference, not copied.

    Raises:
        ValueError: If two selected specimens share an identifier.
    """
    wanted = set(conditions)
    selected: dict[str, AnnData] = {}
    for adata in adatas:
        if adata.uns.get("cell_population") != cell_population:
            continue
        if str(adata.uns.get("cancer_type")) not in wanted:
            continue
        specimen_id = filename_specimen_id(adata)
        if specimen_id in selected:
            raise ValueError(f"Duplicate specimen identifier {specimen_id!r}.")
        selected[specimen_id] = adata
    return dict(sorted(selected.items()))


def specimen_conditions(adatas: Mapping[str, AnnData]) -> dict[str, str]:
    """Map each specimen identifier to its ``uns['cancer_type']``."""
    return {specimen_id: str(adata.uns["cancer_type"]) for specimen_id, adata in adatas.items()}


def specimen_patients(adatas: Mapping[str, AnnData]) -> dict[str, str]:
    """Map each specimen identifier to its :func:`patient_id`."""
    return {specimen_id: patient_id(adata) for specimen_id, adata in adatas.items()}


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _menopause_status(adata: AnnData) -> str:
    """Return ``uns['menopause_status']``, or raise if it is missing."""
    status = adata.uns.get("menopause_status")
    if _is_missing(status):
        raise ValueError(f"Specimen {filename_specimen_id(adata)} has no menopause_status.")
    return str(status).strip()


def _keep_for_normal_menopause(adata: AnnData, normal_menopause: str) -> bool:
    """Keep every non-Normal specimen; keep a Normal specimen only at this menopause status."""
    if str(adata.uns.get("cancer_type")) != "Normal":
        return True
    return _menopause_status(adata) == normal_menopause


@dataclass(frozen=True)
class Cohort:
    """Specimens selected for a comparison, with their conditions and patients."""

    specimens: dict[str, AnnData]
    conditions: dict[str, str]
    patients: dict[str, str]

    @classmethod
    def from_objects(
        cls,
        adatas: Iterable[AnnData],
        conditions: Iterable[str],
        *,
        cell_population: str = "Total",
        normal_menopause: str | None = None,
    ) -> "Cohort":
        """Select specimens with :func:`select_specimens` and record their metadata.

        Args:
            adatas: Candidate specimens.
            conditions: ``uns['cancer_type']`` values to keep.
            cell_population: ``uns['cell_population']`` value to keep.
            normal_menopause: If set, a Normal specimen is kept only when its
                ``menopause_status`` equals this value. Other conditions are
                not filtered. A Normal specimen with no menopause status raises.
        """
        specimens = select_specimens(adatas, conditions, cell_population=cell_population)
        if normal_menopause is not None:
            specimens = {
                specimen_id: adata
                for specimen_id, adata in specimens.items()
                if _keep_for_normal_menopause(adata, normal_menopause)
            }
            specimens = dict(sorted(specimens.items()))
        return cls(specimens, specimen_conditions(specimens), specimen_patients(specimens))

    @property
    def obs(self) -> dict[str, pd.DataFrame]:
        """Cell-level ``obs`` per specimen (by reference)."""
        return {specimen_id: adata.obs for specimen_id, adata in self.specimens.items()}

    def overview(self) -> pd.DataFrame:
        """One row per specimen: condition, patient, and cell count."""
        frame = pd.DataFrame(
            {
                CONDITION_COLUMN: self.conditions,
                PATIENT_COLUMN: self.patients,
                "n_cells": {sid: adata.n_obs for sid, adata in self.specimens.items()},
            }
        )
        frame.index.name = SPECIMEN_COLUMN
        return frame

    def shared_patients(self) -> pd.DataFrame:
        """Specimens whose patient contributed more than one specimen."""
        overview = self.overview()
        shared = overview[overview.duplicated(PATIENT_COLUMN, keep=False)]
        return shared.sort_values([PATIENT_COLUMN, CONDITION_COLUMN])
