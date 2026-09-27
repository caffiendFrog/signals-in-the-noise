"""Tests for signals_in_the_noise.analysis.brca1_cohort."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.brca1_cohort import (
    CELL_TABLE_COLUMNS,
    DONORS,
    PAPER_QC_PASS_TOTALS,
    Donor,
    Genotype,
    build_cell_table,
    compare_qc_pass_counts_to_paper,
    compare_qc_pass_totals_to_paper,
    load_or_build_cell_table,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

TWO_DONORS = (
    Donor("N-a", "GSM1_N-a", Genotype.WT),
    Donor("B1-b", "GSM2_B1-b", Genotype.BRCA1),
)


def _make_adata(n: int, n_noise: int, seed: int) -> AnnData:
    rng = np.random.default_rng(seed)
    obs = pd.DataFrame(
        {
            "total_counts": rng.uniform(500, 5000, n),
            "n_genes_by_counts": rng.integers(200, 3000, n),
            "pct_counts_mt": rng.uniform(0, 30, n),
            "is_noise": [1] * n_noise + [0] * (n - n_noise),
            "unused": 0,
        },
        index=[f"BC{i}-1" for i in range(n)],
    )
    return AnnData(obs=obs)


def _get_dataset_factory(sizes: dict[str, tuple[int, int]]):
    def get_dataset(filename: str) -> AnnData:
        n, n_noise = sizes[filename]
        return _make_adata(n, n_noise, seed=len(filename))

    return get_dataset


@pytest.fixture()
def cell_table() -> pd.DataFrame:
    get_dataset = _get_dataset_factory({"GSM1_N-a.h5ad": (10, 3), "GSM2_B1-b.h5ad": (6, 1)})
    return build_cell_table(get_dataset, TWO_DONORS)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


def test_donors_match_author_cohort_of_8_wt_and_4_brca1():
    genotypes = [d.genotype for d in DONORS]
    assert genotypes.count(Genotype.WT) == 8
    assert genotypes.count(Genotype.BRCA1) == 4


def test_donor_ids_are_unique():
    assert len({d.author_id for d in DONORS}) == len(DONORS)


def test_h5ad_filename_appends_extension_to_geo_stem():
    assert DONORS[0].h5ad_filename == "GSM4909254_N-PM0019-Total.h5ad"


def test_paper_totals_sum_to_reported_cohort_size():
    assert sum(PAPER_QC_PASS_TOTALS.values()) == 59_766


# ---------------------------------------------------------------------------
# build_cell_table
# ---------------------------------------------------------------------------


def test_build_cell_table_has_one_row_per_barcode(cell_table):
    assert len(cell_table) == 16


def test_build_cell_table_keeps_only_documented_columns(cell_table):
    assert tuple(cell_table.columns) == CELL_TABLE_COLUMNS


def test_build_cell_table_assigns_donor_genotype(cell_table):
    genotype_by_donor = cell_table.groupby("donor", observed=True)["genotype"].first()
    assert genotype_by_donor["N-a"] == "WT"
    assert genotype_by_donor["B1-b"] == "BRCA1"


def test_build_cell_table_donor_order_follows_registry(cell_table):
    assert list(cell_table["donor"].cat.categories) == ["N-a", "B1-b"]


def test_build_cell_table_keeps_barcodes(cell_table):
    assert cell_table["barcode"].iloc[0] == "BC0-1"


def test_build_cell_table_rejects_empty_dataset():
    with pytest.raises(ValueError, match="Empty AnnData"):
        build_cell_table(lambda _: AnnData(), TWO_DONORS)


def test_build_cell_table_rejects_missing_qc_columns():
    adata = _make_adata(5, 1, seed=0)
    adata.obs = adata.obs.drop(columns=["pct_counts_mt"])
    with pytest.raises(KeyError, match="pct_counts_mt"):
        build_cell_table(lambda _: adata, TWO_DONORS)


# ---------------------------------------------------------------------------
# load_or_build_cell_table
# ---------------------------------------------------------------------------


def test_load_or_build_writes_cache_then_reuses_it(tmp_path, cell_table):
    path = tmp_path / "cells.csv.gz"
    calls = []

    def build():
        calls.append(1)
        return cell_table

    first = load_or_build_cell_table(path, build, donors=TWO_DONORS)
    second = load_or_build_cell_table(path, build, donors=TWO_DONORS)

    assert len(calls) == 1
    pd.testing.assert_frame_equal(first, second)


def test_load_or_build_force_rebuilds(tmp_path, cell_table):
    path = tmp_path / "cells.csv.gz"
    calls = []

    def build():
        calls.append(1)
        return cell_table

    load_or_build_cell_table(path, build, donors=TWO_DONORS)
    load_or_build_cell_table(path, build, donors=TWO_DONORS, force=True)
    assert len(calls) == 2


# ---------------------------------------------------------------------------
# Paper sanity checks
# ---------------------------------------------------------------------------


def test_compare_qc_pass_counts_to_paper_flags_mismatch(cell_table):
    annotations = pd.DataFrame(
        {"sample-name": ["N-a", "B1-b"], "number-of-cells-after-filtering": [7, 99]}
    )
    result = compare_qc_pass_counts_to_paper(cell_table, annotations).set_index("donor")
    assert result.loc["N-a", "qc_pass_observed"] == 7
    assert bool(result.loc["N-a", "matches"])
    assert not bool(result.loc["B1-b", "matches"])


def test_compare_qc_pass_totals_to_paper_reports_both_genotypes(cell_table):
    result = compare_qc_pass_totals_to_paper(cell_table)
    assert result.loc["WT", "qc_pass_observed"] == 7
    assert result.loc["BRCA1", "qc_pass_paper"] == PAPER_QC_PASS_TOTALS[Genotype.BRCA1]
    assert not result["matches"].any()
