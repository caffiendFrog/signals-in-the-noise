"""Tests for signals_in_the_noise.analysis.rescue_arms."""

import numpy as np
import pandas as pd
import pytest

from signals_in_the_noise.analysis.rescue_arms import (
    QC_PASS,
    UNCLASSIFIED_NOISE,
    assign_cell_classes,
    delta_by_donor,
    rescued_classes,
    rescued_share_by_donor,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _noise_cells(donor: str, genotype: str, n: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "donor": donor,
            "genotype": genotype,
            "total_counts": rng.integers(100, 5000, n).astype(float),
            "n_genes_by_counts": rng.integers(50, 2000, n).astype(float),
            "pct_counts_mt": rng.uniform(0, 60, n),
            "is_noise": 1,
        }
    )


def _classified(rows: list[tuple[str, str, str, str]]) -> pd.DataFrame:
    """Cells from ``(donor, genotype, cell_class, cell_type)`` tuples."""
    return pd.DataFrame(rows, columns=["donor", "genotype", "cell_class", "cell_type"])


def _donor(donor, genotype, *, basal, lp, ml, rescued_lp=0, rescued_ml=0, rescued_class="pbs-3"):
    return (
        [(donor, genotype, QC_PASS, "basal")] * basal
        + [(donor, genotype, QC_PASS, "lp")] * lp
        + [(donor, genotype, QC_PASS, "ml")] * ml
        + [(donor, genotype, QC_PASS, "other")] * 5
        + [(donor, genotype, rescued_class, "lp")] * rescued_lp
        + [(donor, genotype, rescued_class, "ml")] * rescued_ml
        + [(donor, genotype, rescued_class, "other")] * 3
        + [(donor, genotype, "pbs-1", "lp")] * 7
    )


# ---------------------------------------------------------------------------
# assign_cell_classes
# ---------------------------------------------------------------------------


def test_assign_cell_classes_labels_every_barcode():
    qc = _noise_cells("N-0", "WT", 20, seed=0).assign(is_noise=0)
    cells = pd.concat([qc, _noise_cells("N-0", "WT", 500, seed=1)], ignore_index=True)
    classes = assign_cell_classes(cells)["cell_class"]
    assert (classes.iloc[:20] == QC_PASS).all()
    assert set(classes.iloc[20:]) <= {"pbs-1", "pbs-2", "pbs-3", UNCLASSIFIED_NOISE}
    assert (classes.iloc[20:] == "pbs-3").any()


def test_assign_cell_classes_cutoffs_ignore_non_reference_donors():
    wt = _noise_cells("N-0", "WT", 500, seed=1)
    brca1 = _noise_cells("B1-0", "BRCA1", 500, seed=2)
    shifted = brca1.assign(n_genes_by_counts=brca1["n_genes_by_counts"] * 10)

    wt_only = assign_cell_classes(pd.concat([wt, brca1], ignore_index=True))
    with_shift = assign_cell_classes(pd.concat([wt, shifted], ignore_index=True))
    pd.testing.assert_series_equal(
        wt_only["cell_class"].iloc[:500], with_shift["cell_class"].iloc[:500]
    )


def test_assign_cell_classes_requires_reference_noise():
    cells = _noise_cells("B1-0", "BRCA1", 10, seed=0)
    with pytest.raises(ValueError, match="No WT noise cells"):
        assign_cell_classes(cells)


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------


def test_arm_d_rescues_pbs2_and_pbs3():
    assert rescued_classes("D", "B") == ("pbs-2", "pbs-3")


def test_rescued_classes_requires_nested_arms():
    with pytest.raises(ValueError, match="not a subset"):
        rescued_classes("B", "D")


# ---------------------------------------------------------------------------
# delta_by_donor
# ---------------------------------------------------------------------------


def test_delta_by_donor_lp_fractions_within_epithelium():
    cells = _classified(_donor("N-0", "WT", basal=30, lp=40, ml=30, rescued_lp=10, rescued_ml=0))
    row = delta_by_donor(cells).iloc[0]
    assert row["n_epithelial_B"] == 100
    assert row["lp_fraction_B"] == pytest.approx(0.4)
    assert row["n_epithelial_D"] == 110
    assert row["lp_fraction_D"] == pytest.approx(50 / 110)
    assert row["delta"] == pytest.approx(50 / 110 - 0.4)


def test_delta_equals_rescued_share_times_lp_excess():
    cells = _classified(
        _donor("N-0", "WT", basal=30, lp=40, ml=30, rescued_lp=6, rescued_ml=4)
        + _donor("N-1", "WT", basal=10, lp=50, ml=40, rescued_lp=1, rescued_ml=9, rescued_class="pbs-2")
    )
    table = delta_by_donor(cells)
    decomposed = table["rescued_share"] * (table["rescued_lp_fraction"] - table["lp_fraction_B"])
    assert table["delta"].to_numpy() == pytest.approx(decomposed.to_numpy())


def test_delta_by_donor_reports_share_per_rescued_class():
    cells = _classified(
        _donor("N-0", "WT", basal=30, lp=40, ml=30, rescued_lp=5)
        + [("N-0", "WT", "pbs-2", "ml")] * 5
    )
    row = delta_by_donor(cells).iloc[0]
    assert row["share_pbs-3"] == pytest.approx(5 / 110)
    assert row["share_pbs-2"] == pytest.approx(5 / 110)
    assert row["rescued_share"] == pytest.approx(10 / 110)


def test_delta_by_donor_keeps_categorical_donor_order():
    cells = _classified(
        _donor("N-9", "WT", basal=10, lp=10, ml=10) + _donor("N-1", "WT", basal=10, lp=10, ml=10)
    )
    cells["donor"] = pd.Categorical(cells["donor"], categories=["N-9", "N-1"], ordered=True)
    assert delta_by_donor(cells)["donor"].tolist() == ["N-9", "N-1"]


def test_delta_is_zero_without_rescued_epithelium():
    cells = _classified(_donor("N-0", "WT", basal=30, lp=40, ml=30))
    row = delta_by_donor(cells).iloc[0]
    assert row["delta"] == 0.0
    assert row["rescued_share"] == 0.0
    assert np.isnan(row["rescued_lp_fraction"])


# ---------------------------------------------------------------------------
# rescued_share_by_donor
# ---------------------------------------------------------------------------


def test_rescued_share_by_donor_counts_arm_cells_only():
    cells = _classified(_donor("N-0", "WT", basal=30, lp=40, ml=30, rescued_lp=6, rescued_ml=4))
    row = rescued_share_by_donor(cells).iloc[0]
    assert row["n_qc_pass"] == 105
    assert row["n_pbs-3"] == 13
    assert row["n_pbs-2"] == 0
    assert row["rescued_share"] == pytest.approx(13 / 118)
