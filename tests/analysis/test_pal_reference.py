"""Pal reference checks with a known answer that the code has to recover."""

import json

import numpy as np
import pandas as pd
import pytest

from signals_in_the_noise.analysis.pal_reference import (
    donor_cluster_counts,
    label_counts_match,
    largest_enrichment,
    load_pal_labels,
    map_pal_samples,
    near_threshold_enrichment,
    normalize_barcode,
    published_p_in_range,
    quasi_poisson_interaction_p,
    reference_passes,
    run_pal_reference,
    write_pal_reference_outputs,
)

_BARCODE = "AAACCTGAGCAAATCA-1"


def test_barcode_keeps_the_10x_suffix_and_drops_a_sample_prefix():
    assert normalize_barcode(f"N-PM0095_{_BARCODE}") == _BARCODE
    assert normalize_barcode(_BARCODE.lower()) == _BARCODE
    with pytest.raises(ValueError):
        normalize_barcode("not-a-barcode")


def test_library_names_do_not_cross_the_brca1_and_wt_0023_samples():
    donors = pd.DataFrame(
        {
            "specimen_id": ["GSM1_B1-MH0023", "GSM2_N-MH0023-Total", "GSM3_N-PM0095-Total"],
            "donor": ["B1-0023", "N-0123", "N-0093"],
            "genotype": ["BRCA1", "WT", "WT"],
        }
    ).set_index("specimen_id")
    mapped = map_pal_samples(["B1-MH0023", "N-MH0023-total", "N-0093", "N-PM0095-total"], donors)
    assert list(mapped) == ["B1-0023", "N-0123", "N-0093", "N-0093"]


def test_prefixed_export_matches_only_when_the_retained_set_is_exact():
    cells = pd.DataFrame(
        {
            "donor": ["B1", "B1", "W1", "W1"],
            "barcode": [_BARCODE, "AAACCTGAGCAAATCA-2", "AAACCTGAGCAAATCA-3", "AAACCTGAGCAAATCA-4"],
            "is_noise": [0, 0, 0, 1],
        }
    )
    labels = pd.DataFrame(
        {
            "donor": ["B1", "B1", "W1"],
            "barcode": [f"lib_{_BARCODE}", "AAACCTGAGCAAATCA-2", "pref_AAACCTGAGCAAATCA-3"],
        }
    )
    donors = pd.DataFrame(
        {"donor": ["B1", "W1"], "n_retained": [2, 1], "specimen_id": ["s-b", "s-w"]}
    ).set_index("specimen_id")
    same = label_counts_match(labels, cells, donors, expected_total=3)
    assert same["match"]
    short = labels.iloc[:2]
    assert not label_counts_match(short, cells, donors, expected_total=3)["match"]
    published_total = label_counts_match(labels, cells, donors)
    assert published_total["expected_total"] == 59766
    assert not published_total["match"]
    duplicated = pd.concat([labels, labels.iloc[:1]], ignore_index=True)
    with pytest.raises(ValueError, match="repeated"):
        label_counts_match(duplicated, cells, donors, expected_total=3)


def test_missing_clusters_are_explicit_zeros():
    labels = pd.DataFrame(
        {
            "donor": ["B1", "B1", "W1"],
            "genotype": ["BRCA1", "BRCA1", "WT"],
            "cluster": ["LP", "LP", "basal"],
        }
    )
    counts = donor_cluster_counts(labels).set_index(["donor", "cluster"])["count"]
    assert counts.loc[("W1", "LP")] == 0
    assert counts.loc[("B1", "basal")] == 0
    assert counts.loc[("B1", "LP")] == 2


@pytest.mark.filterwarnings(
    "ignore:Perfect separation:statsmodels.tools.sm_exceptions.PerfectSeparationWarning"
)
def test_interaction_p_moves_when_cluster_composition_depends_on_genotype():
    null = _six_donor_counts(separated=False)
    separated = _six_donor_counts(separated=True)
    assert quasi_poisson_interaction_p(null, dispersion="pearson", test="f") == pytest.approx(1)
    assert quasi_poisson_interaction_p(separated, dispersion="deviance", test="chi_square") < 0.01
    skewed = _six_donor_counts(separated=False, skewed=True)
    pearson = quasi_poisson_interaction_p(skewed, dispersion="pearson", test="f")
    deviance = quasi_poisson_interaction_p(skewed, dispersion="deviance", test="f")
    assert pearson != pytest.approx(deviance)


def test_reproduction_requires_exact_counts_and_the_published_p_inside_the_range():
    assert reference_passes(counts_match=True, p_values=[0.10, 0.14, 0.20])
    assert not reference_passes(counts_match=True, p_values=[0.15, 0.30])
    assert not reference_passes(counts_match=False, p_values=[0.10, 0.20])
    assert published_p_in_range([0.10, 0.20], published=0.5) is False
    assert not reference_passes(counts_match=True, p_values=[np.nan], published=0.14)


def test_enrichment_uses_each_donors_own_threshold_and_keeps_ties():
    cells, labels, thresholds = _enrichment_cohort()
    enrichment = near_threshold_enrichment(cells, labels, thresholds, fraction=0.05)
    genes = enrichment.loc[enrichment["threshold"] == "genes_lower"].set_index("cell_type")
    assert genes.loc["LP", "n_band_total"] == 10
    assert genes.loc["LP", "enrichment"] == pytest.approx(20 / 7)
    assert largest_enrichment(enrichment)["c"] == pytest.approx(20 / 7)

    shared = thresholds.copy()
    shared["qc_genes_lower"] = 100
    shared_enrichment = near_threshold_enrichment(cells, labels, shared, fraction=0.05)
    shared_genes = shared_enrichment.loc[shared_enrichment["threshold"] == "genes_lower"]
    shared_lp = shared_genes.loc[shared_genes["cell_type"] == "LP", "enrichment"].iloc[0]
    assert shared_lp == pytest.approx(40 / 7)

    tied = _tied_band()
    tied_enrichment = near_threshold_enrichment(*tied, fraction=0.05)
    tied_genes = tied_enrichment.loc[tied_enrichment["threshold"] == "genes_lower"]
    assert tied_genes["n_band_total"].iloc[0] == 6


def test_frozen_c_is_the_pooled_fold_and_phase_2_stays(tmp_path):
    cells, labels, donors = _enrichment_inputs()
    params = {
        "phase_1": {
            "pal_proportion_check": {
                "published_p": 0.14,
                "variants": [
                    {"dispersion": "pearson", "test": "f"},
                    {"dispersion": "deviance", "test": "f"},
                    {"dispersion": "pearson", "test": "chi_square"},
                    {"dispersion": "deviance", "test": "chi_square"},
                ],
            },
            "near_threshold_band": {"fraction": 0.05, "c": None},
        },
        "phase_2": {"keep": True},
    }
    (tmp_path / "analysis_params.json").write_text(json.dumps(params), encoding="utf-8")
    export = labels.copy()
    export["cell_type_source"] = "cell_type"
    result = run_pal_reference(cells, donors, export, params["phase_1"])
    assert result["variants"] is None
    assert result["c"] == pytest.approx(20 / 7)
    handoff = write_pal_reference_outputs(result, tmp_path)
    stored = json.loads((tmp_path / "analysis_params.json").read_text(encoding="utf-8"))
    assert stored["phase_2"] == {"keep": True}
    assert stored["phase_1"]["near_threshold_band"]["c"] == pytest.approx(20 / 7)
    assert f"{result['c']:.4f}" in handoff


def test_missing_export_names_the_r_script(tmp_path):
    with pytest.raises(FileNotFoundError, match="export_pal_norm_b1_labels.R"):
        load_pal_labels(tmp_path / "missing.csv")


def _six_donor_counts(*, separated: bool, skewed: bool = False) -> pd.DataFrame:
    rows = []
    for index in range(6):
        genotype = "BRCA1" if index < 3 else "WT"
        donor = f"D{index}"
        if separated:
            counts = {
                0: (14, 10),
                1: (16, 12),
                2: (15, 9),
                3: (10, 14),
                4: (12, 16),
                5: (9, 15),
            }[index]
        elif skewed:
            counts = (1, 40) if index % 2 == 0 else (12, 3)
        else:
            counts = (10, 10)
        for cluster, count in zip(("L", "B"), counts, strict=True):
            rows.append({"donor": donor, "genotype": genotype, "cluster": cluster, "count": count})
    return pd.DataFrame(rows)


def _enrichment_cohort():
    cells, labels, donors = _enrichment_inputs()
    thresholds = donors.reset_index()[
        ["donor", "qc_genes_lower", "qc_genes_upper", "qc_mito_upper", "qc_total_upper"]
    ]
    return cells, labels, thresholds


def _enrichment_inputs():
    """Two donors whose nearest cells differ because their gene floors differ.

    Donor A floor 100: the 5 closest cells are LP at 110 genes.
    Donor B floor 4000: the 5 closest cells are basal at 4100 genes.
    Pooled, LP is 5 of 10 band cells and 35 of 200 retained cells, a fold of 20/7.
    """
    rows = []
    labels = []
    spec = [
        ("A", "BRCA1", 100, [("LP", 110, 5), ("LP", 2000, 15), ("basal", 3000, 80)]),
        (
            "B",
            "WT",
            4000,
            [("LP", 200, 5), ("basal", 4100, 5), ("LP", 5000, 10), ("basal", 8000, 80)],
        ),
    ]
    barcode = 0
    for donor, genotype, _floor, groups in spec:
        for cell_type, genes, count in groups:
            for _ in range(count):
                barcode += 1
                name = f"AAACCTGAGCAAATCA-{barcode}"
                rows.append(
                    {
                        "donor": donor,
                        "barcode": name,
                        "is_noise": 0,
                        "n_genes_by_counts": genes,
                        "pct_counts_mt": 1.0,
                        "total_counts": 500,
                    }
                )
                labels.append(
                    {
                        "donor": donor,
                        "sample": donor,
                        "barcode": name,
                        "cluster": "0",
                        "cell_type": cell_type,
                        "genotype": genotype,
                    }
                )
    cells = pd.DataFrame(rows)
    label_frame = pd.DataFrame(labels)
    donors = pd.DataFrame(
        {
            "specimen_id": ["s-a", "s-b"],
            "donor": ["A", "B"],
            "genotype": ["BRCA1", "WT"],
            "n_retained": [100, 100],
            "qc_genes_lower": [100, 4000],
            "qc_genes_upper": [100000, 100000],
            "qc_mito_upper": [0.2, 0.2],
            "qc_total_upper": [30000, 30000],
        }
    ).set_index("specimen_id")
    return cells, label_frame, donors


def _tied_band():
    rows = []
    labels = []
    for index in range(20):
        name = f"AAACCTGAGCAAATCA-{index + 1}"
        genes = 101 if index < 6 else 5000
        rows.append(
            {
                "donor": "A",
                "barcode": name,
                "is_noise": 0,
                "n_genes_by_counts": genes,
                "pct_counts_mt": 1.0,
                "total_counts": 500,
            }
        )
        labels.append(
            {
                "donor": "A",
                "barcode": name,
                "cell_type": "LP" if index < 6 else "basal",
                "cluster": "0",
            }
        )
    thresholds = pd.DataFrame(
        {
            "donor": ["A"],
            "qc_genes_lower": [100],
            "qc_genes_upper": [100000],
            "qc_mito_upper": [0.2],
            "qc_total_upper": [30000],
        }
    )
    return pd.DataFrame(rows), pd.DataFrame(labels), thresholds
