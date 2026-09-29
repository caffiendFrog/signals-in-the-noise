"""Diagnostics for the tipping-point accounting, not echoes of a stored result."""

import json
import math

import numpy as np
import pandas as pd
import pytest

from signals_in_the_noise.analysis.tipping_point import (
    doublet_only,
    full_lp_shares,
    k_max_value,
    removed_lp_share,
    run_tipping_point,
    signature_ceiling,
    smallest_rejecting_k,
    tipping_zone,
    write_tipping_point_outputs,
)


def test_full_share_equals_retained_at_k_1_and_caps_brca1_at_1():
    lp_share = np.array([0.2, 0.4])
    genotype = np.array(["BRCA1", "WT"])
    n_retained = np.array([100, 100])
    n_lp = np.array([20, 40])
    n_removed = np.array([100, 100])

    at_one = removed_lp_share(lp_share, genotype, 1)
    assert at_one == pytest.approx(lp_share)
    assert full_lp_shares(n_retained, n_lp, n_removed, at_one) == pytest.approx(lp_share)

    capped = removed_lp_share(lp_share, genotype, 100)
    assert capped == pytest.approx([1.0, 0.4])
    assert full_lp_shares(n_retained, n_lp, n_removed, capped) == pytest.approx([0.6, 0.4])

    at_three = removed_lp_share(lp_share, genotype, 3)
    assert full_lp_shares(n_retained, n_lp, n_removed, at_three) == pytest.approx([0.4, 0.4])

    extreme = removed_lp_share(lp_share, genotype, 1, extreme=True)
    assert extreme == pytest.approx([1.0, 0.0])
    assert k_max_value(lp_share, genotype) == pytest.approx(5)
    assert math.isinf(k_max_value(np.array([0.0, 0.4]), genotype))


def test_zone_boundaries_and_the_smallest_rejecting_k():
    assert tipping_zone(2, 2, 4) == "plausible"
    assert tipping_zone(2.1, 2, 4) == "possible_but_implausible"
    assert tipping_zone(4, 2, 4) == "possible_but_implausible"
    assert tipping_zone(4.1, 2, 4) == "impossible"
    assert tipping_zone(None, 2, 4) == "impossible"
    assert tipping_zone(3, 2, math.inf) == "possible_but_implausible"

    def p_at(k):
        return 0.5 if k < 2 else 0.001

    found = smallest_rejecting_k(p_at, 5, 0.01)
    assert p_at(found) <= 0.01
    assert p_at(found - 0.01) > 0.01
    assert smallest_rejecting_k(lambda k: 0.5, 5, 0.01) is None
    assert smallest_rejecting_k(lambda k: 0.001, 5, 0.01) == 0.0


def test_doublet_only_keeps_a_cell_that_also_fails_mito():
    cells = pd.DataFrame(
        {
            "is_noise": [1, 1, 1, 0],
            "is_low_num_genes": [0, 0, 0, 0],
            "is_high_mito": [0, 1, 0, 0],
            "is_high_num_genes": [1, 1, 0, 0],
            "is_high_total_count": [0, 0, 1, 0],
        }
    )
    assert doublet_only(cells).tolist() == [True, False, True, False]


def test_signature_ceiling_uses_predicted_type():
    cells, labels, donors = _two_donor_types()
    enrichment, ceiling = signature_ceiling(cells, labels, donors, 0.05)
    genes = enrichment.loc[enrichment["threshold"] == "genes_lower"].set_index("cell_type")
    assert genes.loc["lp", "enrichment"] == pytest.approx(20 / 7)
    assert ceiling["c"] == pytest.approx(20 / 7)
    assert set(ceiling["winners"]["cell_type"]) == {"lp"}


def test_k_star_is_the_hand_solved_boundary_and_cluster_c_stays(tmp_path):
    cells, labels, shares, donors = _twelve_donor_cohort()
    phase1 = _phase1()
    result = run_tipping_point(
        cells,
        donors,
        phase1,
        labels,
        shares,
        rng=np.random.default_rng(0),
    )
    # Equal within genotype, 100 retained and 100 removed.
    # BRCA1 full share exceeds the WT share of 0.4 once k > 3.
    assert result["k_star"] > 3
    assert result["baseline"]["test"].p_value == pytest.approx(1)
    assert result["k_star"] == pytest.approx(3, abs=0.002)
    assert result["k_max"] == pytest.approx(5)
    assert result["zone"] == "possible_but_implausible"
    assert result["signature_c"] == pytest.approx(1)
    assert result["extreme_statistic"] > result["baseline"]["test"].statistic

    # Dropping 25 doublet-only cells leaves 75 removed. The boundary is then 10/3.
    assert result["n_doublet_only"] == 300
    assert result["sensitivity_k_star"] > 10 / 3
    assert result["sensitivity_k_star"] == pytest.approx(10 / 3, abs=0.002)
    assert result["sensitivity_zone"] == "possible_but_implausible"

    params = tmp_path / "analysis_params.json"
    params.write_text(
        json.dumps(
            {
                "phase_1": {
                    "near_threshold_band": {"c": 3.959, "fraction": 0.05},
                    "tipping_point": {"zones": [{"name": "plausible"}]},
                },
                "phase_2": {"keep": True},
            }
        ),
        encoding="utf-8",
    )
    write_tipping_point_outputs(result, tmp_path)
    stored = json.loads(params.read_text(encoding="utf-8"))
    assert stored["phase_1"]["near_threshold_band"]["c"] == 3.959
    assert stored["phase_1"]["tipping_point"]["zones"] == [{"name": "plausible"}]
    assert stored["phase_1"]["tipping_point"]["signature_c"] == pytest.approx(1)
    assert stored["phase_1"]["tipping_point"]["zone"] == "possible_but_implausible"
    assert stored["phase_2"]["keep"] is True
    handoff = (tmp_path / "03_handoff.md").read_text(encoding="utf-8")
    assert "Post (oophorectomy)" in handoff
    assert "3.9590" in handoff


def _phase1() -> dict:
    return {
        "alpha": {"primary": 0.01},
        "confidence": {"primary": 0.99},
        "tests": {
            "lp_baseline_and_tipping_point": {
                "alternative": "greater",
                "contrast": "BRCA1 > WT",
                "alpha": 0.01,
            }
        },
        "near_threshold_band": {"fraction": 0.05, "c": 3.959},
        "mde": {"power": 0.8},
        "loo": {"level": 0.98},
    }


def _twelve_donor_cohort():
    donors = [f"B{i}" for i in range(4)] + [f"W{i}" for i in range(8)]
    genotype = {donor: "BRCA1" if donor.startswith("B") else "WT" for donor in donors}
    rows = []
    labels = []
    barcode = 0
    for donor in donors:
        n_lp = 20 if genotype[donor] == "BRCA1" else 40
        for index in range(100):
            barcode += 1
            name = f"AAACCTGAGCAAATCA-{barcode}"
            rows.append(_cell(donor, name, noise=0, high_genes=0, high_mito=0))
            labels.append(
                {
                    "donor": donor,
                    "barcode": name,
                    "predicted_type": "lp" if index < n_lp else "other",
                }
            )
        for _ in range(25):
            barcode += 1
            rows.append(
                _cell(donor, f"AAACCTGAGCAAATCA-{barcode}", noise=1, high_genes=1, high_mito=0)
            )
        barcode += 1
        rows.append(_cell(donor, f"AAACCTGAGCAAATCA-{barcode}", noise=1, high_genes=1, high_mito=1))
        for _ in range(74):
            barcode += 1
            rows.append(
                _cell(donor, f"AAACCTGAGCAAATCA-{barcode}", noise=1, high_genes=0, high_mito=1)
            )
    cells = pd.DataFrame(rows)
    label_frame = pd.DataFrame(labels)
    shares = pd.DataFrame(
        {
            "n_retained": 100,
            "n_lp": [20 if genotype[donor] == "BRCA1" else 40 for donor in donors],
            "lp_share": [0.2 if genotype[donor] == "BRCA1" else 0.4 for donor in donors],
        },
        index=donors,
    )
    donor_table = pd.DataFrame(
        {
            "specimen_id": [f"s-{donor}" for donor in donors],
            "donor": donors,
            "genotype": [genotype[donor] for donor in donors],
            "menopause_status": [
                "Post (oophorectomy)" if donor == "B0" else "Pre" for donor in donors
            ],
            "qc_genes_lower": 500,
            "qc_genes_upper": 100000,
            "qc_mito_upper": 0.2,
            "qc_total_upper": 30000,
        }
    ).set_index("specimen_id")
    return cells, label_frame, shares, donor_table


def _cell(donor, barcode, *, noise, high_genes, high_mito):
    return {
        "donor": donor,
        "barcode": barcode,
        "is_noise": noise,
        "is_low_num_genes": 0,
        "is_high_num_genes": high_genes,
        "is_high_mito": high_mito,
        "is_high_total_count": 0,
        "n_genes_by_counts": 3000,
        "pct_counts_mt": 1.0,
        "total_counts": 500,
    }


def _two_donor_types():
    rows = []
    labels = []
    spec = [
        ("A", [("lp", 110, 5), ("lp", 2000, 15), ("basal", 3000, 80)]),
        ("B", [("lp", 200, 5), ("basal", 4100, 5), ("lp", 5000, 10), ("basal", 8000, 80)]),
    ]
    barcode = 0
    for donor, groups in spec:
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
                labels.append({"donor": donor, "barcode": name, "predicted_type": cell_type})
    donors = pd.DataFrame(
        {
            "donor": ["A", "B"],
            "qc_genes_lower": [100, 4000],
            "qc_genes_upper": [100000, 100000],
            "qc_mito_upper": [0.2, 0.2],
            "qc_total_upper": [30000, 30000],
        }
    )
    return pd.DataFrame(rows), pd.DataFrame(labels), donors
