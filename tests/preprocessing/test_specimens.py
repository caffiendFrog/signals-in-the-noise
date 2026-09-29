"""Tests for signals_in_the_noise.preprocessing.specimens."""

import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.preprocessing.specimens import (
    Cohort,
    filename_specimen_id,
    patient_id,
    select_specimens,
    specimen_conditions,
    specimen_patients,
)


def _specimen(
    filename: str, cancer_type: str, patient: str = "0001", cell_population: str = "Total"
) -> AnnData:
    adata = AnnData(obs=pd.DataFrame({"adata-filename": [filename]}, index=["c0"]))
    adata.uns.update(
        {
            "title": f"{cancer_type} {cell_population} cells from Patient {patient}",
            "cancer_type": cancer_type,
            "cell_population": cell_population,
        }
    )
    return adata


def test_filename_specimen_id_strips_extension():
    assert (
        filename_specimen_id(_specimen("GSM4909313_ER-MH0064-T.h5ad", "ER+ tumour"))
        == "GSM4909313_ER-MH0064-T"
    )


@pytest.mark.parametrize(
    ("title", "expected"),
    [
        ("ER+ tumour Total cells from Patient 0029", "0029"),
        ("Normal Total cells from Patient 0230.17", "0230.17"),
        ("ER+ tumour Lymph-node cells from Patient 0040", "0040"),
    ],
)
def test_patient_id_parses_title(title, expected):
    adata = AnnData()
    adata.uns["title"] = title
    assert patient_id(adata) == expected


def test_patient_id_rejects_unexpected_title():
    adata = AnnData()
    adata.uns["title"] = "no patient here"
    with pytest.raises(ValueError, match="patient"):
        patient_id(adata)


def test_select_specimens_filters_condition_and_population_and_sorts_by_id():
    adatas = [
        _specimen("b.h5ad", "ER+ tumour"),
        _specimen("a.h5ad", "Normal"),
        _specimen("c.h5ad", "HER2+ tumour"),
        _specimen("d.h5ad", "Normal", cell_population="Epithelial"),
    ]
    selected = select_specimens(adatas, ["ER+ tumour", "Normal"])
    assert list(selected) == ["a", "b"]
    assert selected["b"] is adatas[0]


def test_select_specimens_keeps_specimens_that_share_a_title():
    """Patient 0029 has two ER+ Total samples with identical titles; both must survive."""
    adatas = [
        _specimen("GSM1_ER-0029-7C.h5ad", "ER+ tumour", patient="0029"),
        _specimen("GSM2_ER-0029-9C.h5ad", "ER+ tumour", patient="0029"),
    ]
    selected = select_specimens(adatas, ["ER+ tumour"])
    assert len(selected) == 2
    assert set(specimen_patients(selected).values()) == {"0029"}


def test_select_specimens_rejects_duplicate_ids():
    with pytest.raises(ValueError, match="Duplicate"):
        select_specimens([_specimen("a.h5ad", "Normal"), _specimen("a.h5ad", "Normal")], ["Normal"])


def test_cohort_collects_metadata_and_reports_shared_patients():
    adatas = [
        _specimen("GSM1_N-0064.h5ad", "Normal", patient="0064"),
        _specimen("GSM2_ER-0064-T.h5ad", "ER+ tumour", patient="0064"),
        _specimen("GSM3_ER-0001.h5ad", "ER+ tumour", patient="0001"),
        _specimen("GSM4_HER2.h5ad", "HER2+ tumour", patient="0308"),
    ]
    cohort = Cohort.from_objects(adatas, ["ER+ tumour", "Normal"])

    assert list(cohort.specimens) == ["GSM1_N-0064", "GSM2_ER-0064-T", "GSM3_ER-0001"]
    assert cohort.obs["GSM3_ER-0001"] is adatas[2].obs
    overview = cohort.overview()
    assert overview.index.name == "specimen_id"
    assert overview.loc["GSM2_ER-0064-T", "patient_id"] == "0064"
    assert overview["n_cells"].tolist() == [1, 1, 1]
    shared = cohort.shared_patients()
    assert set(shared.index) == {"GSM1_N-0064", "GSM2_ER-0064-T"}


def _with_menopause(adata: AnnData, status: str | None) -> AnnData:
    if status is not None:
        adata.uns["menopause_status"] = status
    return adata


def test_normal_menopause_filter_drops_post_normals_and_keeps_post_brca1():
    adatas = [
        _with_menopause(_specimen("GSM_N-PM0092-Total.h5ad", "Normal", "0092"), "Pre"),
        _with_menopause(_specimen("GSM_N-PM0342-Total.h5ad", "Normal", "0342"), "Post"),
        _with_menopause(
            _specimen("GSM_B1-MH0023.h5ad", "BRCA1 pre-neoplastic", "0023"),
            "Post (oophorectomy)",
        ),
        _with_menopause(
            _specimen("GSM_N-PM0092-Epi.h5ad", "Normal", "0092", cell_population="Epithelial"),
            "Pre",
        ),
    ]
    cohort = Cohort.from_objects(adatas, ["Normal", "BRCA1 pre-neoplastic"], normal_menopause="Pre")
    assert list(cohort.specimens) == ["GSM_B1-MH0023", "GSM_N-PM0092-Total"]


def test_normal_menopause_filter_rejects_a_normal_with_no_status():
    adatas = [_specimen("GSM_N-PM0092-Total.h5ad", "Normal", "0092")]
    with pytest.raises(ValueError, match="menopause"):
        Cohort.from_objects(adatas, ["Normal"], normal_menopause="Pre")


def test_specimen_conditions_and_patients():
    selected = select_specimens([_specimen("a.h5ad", "Normal", patient="0064")], ["Normal"])
    assert specimen_conditions(selected) == {"a": "Normal"}
    assert specimen_patients(selected) == {"a": "0064"}
