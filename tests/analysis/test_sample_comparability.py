"""Tests for signals_in_the_noise.analysis.sample_comparability."""

import pandas as pd
import pytest

from signals_in_the_noise.analysis.sample_comparability import (
    DepthMatchingDecision,
    assess_depth_matching,
    compare_genotypes,
    summarize_donors,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cells(depth_by_donor: dict[tuple[str, str], list[float]], noise_depth: float = 100.0):
    """Build a cell table: listed depths are QC-pass cells, plus one noise cell per donor."""
    rows = []
    for (donor, genotype), depths in depth_by_donor.items():
        for depth in depths:
            rows.append((donor, genotype, depth, 0))
        rows.append((donor, genotype, noise_depth, 1))
    cells = pd.DataFrame(rows, columns=["donor", "genotype", "total_counts", "is_noise"])
    cells["n_genes_by_counts"] = cells["total_counts"] / 2
    cells["pct_counts_mt"] = 5.0
    return cells


WT_DONORS = {(f"N-{i}", "WT"): [1000.0, 1000.0, 1000.0] for i in range(8)}


@pytest.fixture()
def balanced_cells() -> pd.DataFrame:
    brca1 = {(f"B1-{i}", "BRCA1"): [1000.0, 1000.0, 1000.0] for i in range(4)}
    return _cells({**WT_DONORS, **brca1})


@pytest.fixture()
def shallow_brca1_cells() -> pd.DataFrame:
    brca1 = {(f"B1-{i}", "BRCA1"): [500.0, 500.0, 500.0] for i in range(4)}
    return _cells({**WT_DONORS, **brca1})


# ---------------------------------------------------------------------------
# summarize_donors
# ---------------------------------------------------------------------------


def test_summarize_donors_has_row_per_population_donor_metric(balanced_cells):
    summary = summarize_donors(balanced_cells)
    assert len(summary) == 3 * 12 * 3


def test_summarize_donors_counts_cells_per_population(balanced_cells):
    summary = summarize_donors(balanced_cells).set_index(["population", "donor", "metric"])
    assert summary.loc[("qc_pass", "N-0", "total_counts"), "n_cells"] == 3
    assert summary.loc[("noise", "N-0", "total_counts"), "n_cells"] == 1
    assert summary.loc[("all_barcodes", "N-0", "total_counts"), "n_cells"] == 4


def test_summarize_donors_reports_median_and_iqr():
    cells = _cells({("N-0", "WT"): [100.0, 200.0, 300.0, 400.0, 500.0]})
    row = summarize_donors(cells, populations={"qc_pass": lambda c: c["is_noise"] == 0})
    row = row.set_index("metric").loc["total_counts"]
    assert row["median"] == 300.0
    assert row["q25"] == 200.0
    assert row["q75"] == 400.0


# ---------------------------------------------------------------------------
# compare_genotypes
# ---------------------------------------------------------------------------


def test_compare_genotypes_balanced_has_ratio_one(balanced_cells):
    comparison = compare_genotypes(summarize_donors(balanced_cells))
    row = comparison.set_index(["population", "metric"]).loc[("qc_pass", "total_counts")]
    assert row["ratio"] == pytest.approx(1.0)
    assert row["n_WT"] == 8
    assert row["n_BRCA1"] == 4
    assert row["permutation_p"] == pytest.approx(1.0)


def test_compare_genotypes_detects_shallow_brca1(shallow_brca1_cells):
    comparison = compare_genotypes(summarize_donors(shallow_brca1_cells))
    row = comparison.set_index(["population", "metric"]).loc[("qc_pass", "total_counts")]
    assert row["BRCA1_median"] == 500.0
    assert row["WT_median"] == 1000.0
    assert row["relative_difference"] == pytest.approx(-0.5)
    assert row["permutation_p"] < 0.1


# ---------------------------------------------------------------------------
# Depth-matching decision
# ---------------------------------------------------------------------------


def test_assess_depth_matching_not_mandatory_when_balanced(balanced_cells):
    decision = assess_depth_matching(compare_genotypes(summarize_donors(balanced_cells)))
    assert not decision.mandatory
    assert "optional" in decision.summary()


def test_assess_depth_matching_mandatory_when_brca1_shallow(shallow_brca1_cells):
    decision = assess_depth_matching(compare_genotypes(summarize_donors(shallow_brca1_cells)))
    assert decision.mandatory
    assert decision.exceeding == ["qc_pass"]
    assert "MANDATORY" in decision.summary()


def test_decision_tolerance_is_strict_inequality():
    decision = DepthMatchingDecision(tolerance=0.2, relative_differences={"qc_pass": 0.2})
    assert not decision.mandatory


def test_decision_flags_positive_and_negative_differences():
    decision = DepthMatchingDecision(
        tolerance=0.2, relative_differences={"qc_pass": 0.25, "noise": -0.21}
    )
    assert decision.exceeding == ["qc_pass", "noise"]


def test_assess_depth_matching_rejects_unknown_population(balanced_cells):
    comparison = compare_genotypes(summarize_donors(balanced_cells))
    with pytest.raises(KeyError, match="missing_population"):
        assess_depth_matching(comparison, populations=["missing_population"])


def test_decision_json_round_trip(tmp_path):
    decision = DepthMatchingDecision(tolerance=0.2, relative_differences={"qc_pass": -0.3})
    path = tmp_path / "decision.json"
    decision.to_json(path)
    assert DepthMatchingDecision.from_json(path) == decision
