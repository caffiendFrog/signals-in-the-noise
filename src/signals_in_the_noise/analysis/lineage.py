"""Lim 2009 lineage labels for the retained cells notebook 03 scores.

Fig 4B cluster ids are not luminal progenitors. This calls the existing
signature scorer on each donor's retained cells and reports the progenitor share.
"""

from collections.abc import Callable, Mapping

import numpy as np
import pandas as pd
from anndata import AnnData

from signals_in_the_noise.analysis.pal_reference import normalize_barcode
from signals_in_the_noise.preprocessing.gse161529 import GSE161529

LP_TYPE = "lp"


def retained_lineage_labels(
    cells: pd.DataFrame,
    adatas: Mapping[str, AnnData],
    donors: pd.DataFrame,
    *,
    scorer: Callable[[AnnData], AnnData] | None = None,
) -> dict:
    """Label retained cells and return each donor's luminal-progenitor share.

    The scorer is ``annotate_epithelial_cell_typing`` with highly variable
    gene filtering left off. Removed cells are not scored. ``lp_share`` is the
    fraction of that donor's retained cells whose winning type is ``lp``.
    """
    if scorer is None:
        scorer = _lim_signature_scorer()
    retained = _retained_cells(cells)
    specimen_of = _specimen_by_donor(donors)
    labels = []
    for donor, barcodes in retained.groupby("donor")["barcode"]:
        specimen_id = specimen_of.get(donor)
        if specimen_id is None or specimen_id not in adatas:
            raise ValueError(f"No count matrix for donor {donor}.")
        scored = scorer(_retained_matrix(adatas[specimen_id], set(barcodes)))
        labels.append(_label_frame(donor, scored, set(barcodes)))
    labeled = pd.concat(labels, ignore_index=True)
    return {"labels": labeled, "shares": _lp_shares(labeled)}


def _lim_signature_scorer() -> Callable[[AnnData], AnnData]:
    """Existing Lim 2009 scorer. Constructing ``GSE161529()`` would load the study."""
    tool = GSE161529.__new__(GSE161529)
    tool.STUDY_ID = GSE161529.STUDY_ID
    tool.random_kwargs = {"random_state": 0}

    def score(adata: AnnData) -> AnnData:
        return tool.annotate_epithelial_cell_typing(adata, hvg_only=False)

    return score


def _retained_cells(cells: pd.DataFrame) -> pd.DataFrame:
    required = {"donor", "barcode", "is_noise"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"Cell table is missing {sorted(missing)}.")
    retained = cells.loc[cells["is_noise"].astype(int) == 0, ["donor", "barcode"]].copy()
    retained["donor"] = retained["donor"].astype(str)
    retained["barcode"] = retained["barcode"].map(normalize_barcode)
    if retained.empty:
        raise ValueError("The cell table has no retained cells.")
    if retained.duplicated().any():
        raise ValueError("A retained donor-barcode pair is repeated.")
    return retained


def _specimen_by_donor(donors: pd.DataFrame) -> dict[str, str]:
    frame = donors.copy()
    if "specimen_id" not in frame.columns:
        frame = frame.reset_index(names="specimen_id")
    if "donor" not in frame.columns:
        raise ValueError("Donor table is missing a donor column.")
    if frame["donor"].duplicated().any():
        raise ValueError("A donor id is repeated.")
    return {
        str(donor): str(specimen)
        for specimen, donor in zip(frame["specimen_id"], frame["donor"], strict=True)
    }


def _retained_matrix(adata: AnnData, barcodes: set[str]) -> AnnData:
    names = pd.Index([normalize_barcode(name) for name in adata.obs_names])
    if names.duplicated().any():
        raise ValueError("A specimen has repeated barcodes.")
    keep = np.asarray(names.isin(barcodes))
    if int(keep.sum()) != len(barcodes):
        raise ValueError("A retained barcode is missing from the count matrix.")
    subset = adata[keep].copy()
    subset.obs_names = pd.Index(names[keep])
    return subset


def _label_frame(donor: str, scored: AnnData, barcodes: set[str]) -> pd.DataFrame:
    if "predicted_type" not in scored.obs.columns:
        raise ValueError("The lineage scorer did not return predicted_type.")
    if scored.obs["predicted_type"].isna().any():
        raise ValueError(f"Donor {donor} has a retained cell with no lineage label.")
    names = pd.Index([normalize_barcode(name) for name in scored.obs_names])
    if set(names) != barcodes or names.duplicated().any():
        raise ValueError(f"Lineage labels for donor {donor} do not match the retained cells.")
    types = scored.obs["predicted_type"].astype(str)
    return pd.DataFrame(
        {"donor": donor, "barcode": names.to_numpy(), "predicted_type": types.to_numpy()}
    )


def _lp_shares(labels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor, group in labels.groupby("donor", sort=True):
        n_retained = len(group)
        n_lp = int((group["predicted_type"] == LP_TYPE).sum())
        rows.append(
            {
                "donor": donor,
                "n_retained": n_retained,
                "n_lp": n_lp,
                "lp_share": n_lp / n_retained,
            }
        )
    return pd.DataFrame(rows).set_index("donor")
