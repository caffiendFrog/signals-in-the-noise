"""Tests for the BRCA1 QC-filtering cohort tables and frozen parameters."""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.qc_filtering import (
    EXPECTED_DONORS,
    assert_published_cell_counts,
    brca1_wt_cohort,
    cell_table,
    donor_label,
    inferred_source_prefix,
    validate_cohort,
    write_phase1_outputs,
)
from signals_in_the_noise.preprocessing.qc import NOISE_FLAG_COLUMNS

ANNOTATIONS = Path(__file__).parents[2] / "data" / "GSE161529_annotations_df.csv"

# Removed-cell flags by file stem. Several flags are set together where the
# exclusive-reason priority has to choose.
_REMOVED_FLAGS = {
    "GSM4909253_N-PM0092-Total": {"is_low_num_genes": 1},
    "GSM4909254_N-PM0019-Total": {"is_high_num_genes": 1, "is_high_total_count": 1},
    "GSM4909257_N-PM0095-Total": {"is_high_total_count": 1},
    "GSM4909277_B1-KCF0894": {"is_high_mito": 1, "is_low_num_genes": 1},
}
_DEFAULT_REMOVED_FLAGS = {"is_high_mito": 1}
# Even count, so the median is the mean of the two middle values (250.5),
# which is neither the minimum nor the mean of all four.
_SHALLOW_COUNTS = [100, 200, 301, 1000]

# stem, condition, patient, menopause, parity, median UMI, mito cap, gene floor
_COHORT = [
    ("GSM4909253_N-PM0092-Total", "Normal", "0092", "Pre", "Nulliparous", 5000, 0.2, 500),
    ("GSM4909254_N-PM0019-Total", "Normal", "0019", "Pre", "Nulliparous", 5000, 0.3, 500),
    ("GSM4909257_N-PM0095-Total", "Normal", "0093", "Pre", "Nulliparous", 5000, 0.2, 500),
    ("GSM4909261_N-PM0230-Total", "Normal", "0230.17", "Pre", "Nulliparous", 5000, 0.3, 500),
    ("GSM4909263_N-MH0064-Total", "Normal", "0064", "Pre", "Nulliparous", 5000, 0.3, 500),
    ("GSM4909265_N-PM0233-Total", "Normal", "0233", "Pre", "Parous", 5000, 0.3, 500),
    ("GSM4909266_N-MH0169-Total", "Normal", "0169", "Pre", "Parous", 5000, 0.3, 500),
    ("GSM4909268_N-MH0023-Total", "Normal", "0123", "Pre", "Parous", 5000, 0.3, 500),
    ("GSM4909277_B1-KCF0894", "BRCA1 pre-neoplastic", "0894", "Pre", "Nulliparous", 5000, 0.2, 500),
    ("GSM4909278_B1-MH0033", "BRCA1 pre-neoplastic", "0033", "Pre", "Nulliparous", 5000, 0.2, 300),
    (
        "GSM4909279_B1-MH0023",
        "BRCA1 pre-neoplastic",
        "0023",
        "Post (oophorectomy)",
        "Parous",
        5000,
        0.2,
        500,
    ),
    (
        "GSM4909280_B1-MH0090",
        "BRCA1 pre-neoplastic",
        "0090",
        "Post (oophorectomy)",
        "Parous",
        5000,
        0.2,
        500,
    ),
]


def _specimen(
    stem: str,
    condition: str,
    patient: str,
    menopause: str,
    parity: str,
    total_counts: float,
    mito_upper: float,
    genes_lower: float,
    *,
    cell_population: str = "Total",
) -> AnnData:
    n_retained = 3
    n = n_retained + 1
    counts = _SHALLOW_COUNTS if stem.endswith("N-PM0095-Total") else [total_counts] * n
    removed_flags = _REMOVED_FLAGS.get(stem, _DEFAULT_REMOVED_FLAGS)
    flag_columns = {
        column: [0] * (n - 1) + [int(removed_flags.get(column, 0))] for column in NOISE_FLAG_COLUMNS
    }
    obs = pd.DataFrame(
        {
            "adata-filename": f"{stem}.h5ad",
            "total_counts": counts,
            "n_genes_by_counts": [1000] * n,
            "pct_counts_mt": [5.0] * n,
            "log1p_total_counts": [float(np.log1p(count)) for count in counts],
            "log1p_n_genes_by_counts": [float(np.log1p(1000))] * n,
            "is_noise": [0] * (n - 1) + [1],
            **flag_columns,
        },
        index=[f"{stem}_{i}" for i in range(n)],
    )
    adata = AnnData(obs=obs)
    adata.uns.update(
        {
            "title": f"{condition} Total cells from Patient {patient}",
            "cancer_type": condition,
            "cell_population": cell_population,
            "menopause_status": menopause,
            "parity": parity,
            "qc_mito_upper": mito_upper,
            "qc_genes_lower": genes_lower,
            "qc_genes_upper": 6500,
            "qc_total_upper": 70000,
            # The cached objects store these counts as numpy floats.
            "num_cells_before": np.float64(n),
            "num_cells_after": np.float64(n_retained),
            "num_genes_before": 100,
        }
    )
    return adata


def _adatas() -> list[AnnData]:
    adatas = [_specimen(*row) for row in _COHORT]
    adatas.append(
        _specimen(
            "GSM4909270_N-PM0342-Total", "Normal", "0342", "Post", "Nulliparous", 9000, 0.2, 500
        )
    )
    epithelial = _specimen(
        "GSM4909256_N-PM0095-Epi", "Normal", "0093", "Pre", "Nulliparous", 5000, 0.2, 500
    )
    epithelial.uns["cell_population"] = "Epithelial"
    adatas.append(epithelial)
    return adatas


def test_inferred_source_prefix_reads_the_sample_token():
    assert inferred_source_prefix("GSM4909277_B1-KCF0894") == "KCF"
    assert inferred_source_prefix("GSM4909257_N-PM0095-Total") == "PM"
    with pytest.raises(ValueError, match="source prefix"):
        inferred_source_prefix("no-prefix")


def test_published_count_mismatch_raises():
    adata = _specimen(*_COHORT[0])
    adata.uns["num_cells_after"] = 99
    with pytest.raises(ValueError, match="GEO"):
        assert_published_cell_counts(adata)


def _annotation_donor_id(sample_name: str) -> str:
    """Donor id from the annotation sample name, not from ``donor_label``."""
    suffix = "-total"
    if sample_name.endswith(suffix):
        return sample_name[: -len(suffix)]
    return sample_name


def test_expected_donors_match_the_annotation_table():
    """The locked donor set is the annotation table's own cohort, not the test fixture."""
    annotations = pd.read_csv(ANNOTATIONS)
    total = annotations[annotations["cell-population"] == "Total"]
    brca1 = total[total["cancer-type"] == "BRCA1 pre-neoplastic"]
    premenopausal_wt = total[
        (total["cancer-type"] == "Normal") & (total["menopause-status"] == "Pre")
    ]
    ids = set(brca1["sample-name"].map(_annotation_donor_id)) | set(
        premenopausal_wt["sample-name"].map(_annotation_donor_id)
    )
    assert ids == set(EXPECTED_DONORS)

    post_brca1 = brca1[brca1["menopause-status"] != "Pre"]
    assert set(post_brca1["sample-name"]) == {"B1-0023", "B1-0090"}
    assert set(post_brca1["menopause-status"]) == {"Post (oophorectomy)"}


def test_cohort_keeps_post_oophorectomy_brca1_and_drops_post_normal():
    cohort = brca1_wt_cohort(_adatas())
    labels = {donor_label(adata) for adata in cohort.specimens.values()}
    assert labels == set(EXPECTED_DONORS)
    assert "N-0342" not in labels
    assert all(adata.uns["cell_population"] == "Total" for adata in cohort.specimens.values())


def test_cell_table_applies_exclusive_reason_priority():
    cells = cell_table(brca1_wt_cohort(_adatas()))
    removed = cells[cells["qc_status"] == "removed"].set_index("donor")
    # Low genes alone, high genes beating high library size, high library size
    # alone, and mitochondrial fraction beating low genes.
    assert removed.loc["N-0092", "exclusive_reason"] == "low_genes"
    assert removed.loc["N-0019", "exclusive_reason"] == "high_genes"
    assert removed.loc["N-0093", "exclusive_reason"] == "high_total"
    assert removed.loc["B1-0894", "exclusive_reason"] == "high_mito"

    retained = cells[cells["qc_status"] == "retained"]
    assert retained["exclusive_reason"].isna().all()
    # Patient 0093 is stored under the PM0095 file name.
    assert set(cells.loc[cells["donor"] == "N-0093", "specimen_id"]) == {
        "GSM4909257_N-PM0095-Total"
    }


def test_cell_table_rejects_flags_that_disagree_with_is_noise():
    adata = _specimen(*_COHORT[0])
    adata.obs.loc[adata.obs.index[0], "is_high_mito"] = 1
    with pytest.raises(ValueError, match="disagree"):
        cell_table(brca1_wt_cohort([adata]))


def _donor_frame() -> pd.DataFrame:
    rows = [
        {
            "donor": donor,
            "genotype": "BRCA1" if donor.startswith("B1-") else "WT",
            "menopause_status": "Pre",
        }
        for donor in EXPECTED_DONORS
    ]
    return pd.DataFrame(rows)


def test_validate_cohort_rejects_a_missing_donor_a_duplicate_and_a_post_wt():
    missing = _donor_frame().iloc[1:]
    with pytest.raises(ValueError, match="Missing"):
        validate_cohort(missing)

    duplicated = pd.concat([_donor_frame(), _donor_frame().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="repeated"):
        validate_cohort(duplicated)

    post_wt = _donor_frame()
    post_wt.loc[post_wt["genotype"] == "WT", "menopause_status"] = "Post"
    with pytest.raises(ValueError, match="premenopausal"):
        validate_cohort(post_wt)


def test_write_phase1_outputs_freezes_parameters_and_roundtrips_cells(tmp_path):
    params_path = tmp_path / "analysis_params.json"
    params_path.write_text(json.dumps({"phase_2": {"kept": True}}), encoding="utf-8")
    handoff = write_phase1_outputs(brca1_wt_cohort(_adatas()), tmp_path)

    donors = pd.read_csv(tmp_path / "donors.csv")
    cells = pd.read_csv(tmp_path / "cells.csv")
    assert not cells.duplicated(["donor", "barcode"]).any()
    assert (cells["qc_status"].eq("removed") == cells["is_noise"].eq(1)).all()
    retained = cells["qc_status"] == "retained"
    assert cells.loc[retained, "exclusive_reason"].isna().all()
    assert not cells.loc[~retained, "exclusive_reason"].isin(["nan", "<NA>", "None"]).any()
    written_rate = cells.groupby("donor")["is_noise"].mean()
    recorded_rate = donors.set_index("donor")["removed_fraction"]
    pd.testing.assert_series_equal(
        written_rate.sort_index(), recorded_rate.sort_index(), check_names=False
    )

    params = json.loads(params_path.read_text(encoding="utf-8"))
    assert params["phase_2"] == {"kept": True}
    phase1 = params["phase_1"]
    assert phase1["alpha"] == {"primary": 0.01, "secondary": 0.02}
    assert phase1["tests"]["overall_removed_fraction"]["family"] is None
    n_brca1 = int((donors["genotype"] == "BRCA1").sum())
    n_wt = int((donors["genotype"] == "WT").sum())
    assert phase1["mde"]["n_assignments"] == math.comb(n_brca1 + n_wt, n_brca1)
    assert phase1["near_threshold_band"]["c"] is None

    # The cohort mode ties. The stored cap is that tie's smaller value, and it
    # is not the mode of the WT donors alone.
    stored_mito = phase1["sensitivity"]["uniform_thresholds"]["mito_upper"]
    cohort_mito = donors["qc_mito_upper"].mode().min()
    wt_mito = donors.loc[donors["genotype"] == "WT", "qc_mito_upper"].mode().min()
    assert stored_mito == pytest.approx(cohort_mito)
    assert stored_mito != pytest.approx(wt_mito)
    genes_lower = phase1["sensitivity"]["uniform_thresholds"]["genes_lower"]
    assert genes_lower == donors["qc_genes_lower"].mode().min()
    assert genes_lower > donors["qc_genes_lower"].min()

    depth = cells.groupby("donor")["total_counts"].median().min()
    assert phase1["sensitivity"]["depth_equalized"]["target_median_umi"] == depth
    assert phase1["sensitivity"]["depth_equalized"]["qc_rule_after_thinning"] == "uniform"
    assert phase1["sensitivity"]["mad"]["n_mads"] == 3.0
    assert len(phase1["pal_proportion_check"]["variants"]) == 4
    assert phase1["batch"] is None

    assert f"{depth:.1f}" in handoff
    assert (tmp_path / "00_handoff.md").read_text(encoding="utf-8") == handoff
