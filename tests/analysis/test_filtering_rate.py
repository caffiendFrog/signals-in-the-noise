"""Filtering-rate tests on synthetic donors. Nothing here echoes the 12-donor cohort."""

import math

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.stats import pearsonr

from signals_in_the_noise.analysis.donor_permutation import bh_adjust, exact_permutation_test
from signals_in_the_noise.analysis.filtering_rate import (
    depth_removed_counts,
    donor_measures,
    judge_rerun,
    load_specimen_adatas,
    overall_robustness,
    removed_counts_under_rule,
    run_filtering_rate,
    write_filtering_rate_outputs,
)
from signals_in_the_noise.preprocessing.qc import QcThresholds, annotate_qc_metrics

# donor, genotype, specimen, menopause, source, n_cells, n_low_genes, n_high_mito, n_uniform
_COHORT = [
    ("B1", "BRCA1", "s-b1", "Pre", "KCF", 10, 1, 0, 5),
    ("B2", "BRCA1", "s-b2", "Post (oophorectomy)", "MH", 10, 1, 1, 4),
    ("B3", "BRCA1", "s-b3", "Pre", "MH", 10, 2, 1, 3),
    ("W1", "WT", "s-w1", "Pre", "PM", 10, 2, 2, 2),
    ("W2", "WT", "s-w2", "Pre", "PM", 10, 3, 2, 1),
    ("W3", "WT", "s-w3", "Pre", "PM", 40, 2, 2, 8),
]


def _phase1() -> dict:
    return {
        "alpha": {"primary": 0.01, "secondary": 0.02},
        "confidence": {"primary": 0.8, "secondary": 0.8, "loo": 0.8},
        "sensitivity": {
            "uniform_thresholds": {
                "mito_upper": 0.2,
                "genes_lower": 100,
                "genes_upper": 5000,
                "total_upper": 30000,
            },
            "mad": {"n_mads": 3},
            "depth_equalized": {
                "target_median_umi": 1_000_000,
                "seeds": [0, 1],
                "qc_rule_after_thinning": "uniform",
            },
        },
    }


def _adata(n_cells: int, n_flagged: int) -> AnnData:
    matrix = np.ones((n_cells, 4), dtype=np.int64)
    adata = AnnData(matrix, var=pd.DataFrame(index=["G0", "G1", "G2", "MT-CO1"]))
    annotate_qc_metrics(adata, percent_top=None)
    genes = np.full(n_cells, 1000.0)
    genes[:n_flagged] = 10.0
    adata.obs["n_genes_by_counts"] = genes
    adata.obs["log1p_n_genes_by_counts"] = np.log1p(genes)
    adata.obs["total_counts"] = 500.0
    adata.obs["log1p_total_counts"] = np.log1p(500.0)
    adata.obs["pct_counts_mt"] = 1.0
    return adata


def _cohort() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    cell_rows = []
    donor_rows = []
    adatas = {}
    for donor, genotype, specimen, menopause, source, n_cells, n_low, n_mito, n_uniform in _COHORT:
        for index in range(n_cells):
            low = index < n_low
            mito = n_low <= index < n_low + n_mito
            if mito:
                reason = "high_mito"
            elif low:
                reason = "low_genes"
            else:
                reason = pd.NA
            cell_rows.append(
                {
                    "donor": donor,
                    "genotype": genotype,
                    "is_noise": int(low or mito),
                    "is_low_num_genes": int(low),
                    "is_high_num_genes": 0,
                    "is_high_mito": int(mito),
                    "is_high_total_count": 0,
                    "exclusive_reason": reason,
                }
            )
        donor_rows.append(
            {
                "specimen_id": specimen,
                "donor": donor,
                "genotype": genotype,
                "menopause_status": menopause,
                "source_prefix": source,
                "n_cells": n_cells,
            }
        )
        adatas[specimen] = _adata(n_cells, n_uniform)
    cells = pd.DataFrame(cell_rows)
    donors = pd.DataFrame(donor_rows).set_index("specimen_id")
    return cells, donors, adatas


def test_frozen_thresholds_change_the_removed_count():
    obs = pd.DataFrame(
        {
            "n_genes_by_counts": [100, 100, 100],
            "pct_counts_mt": [0.0, 0.0, 0.0],
            "total_counts": [500, 500, 500],
        }
    )
    specimen_obs = {"s": obs}
    donor_of = pd.Series({"s": "d"})
    strict = QcThresholds(0.2, 500, 5000, 30000)
    loose = QcThresholds(0.2, -1, 1e9, 1e9)
    strict_count = removed_counts_under_rule(specimen_obs, donor_of, lambda _obs: strict)
    loose_count = removed_counts_under_rule(specimen_obs, donor_of, lambda _obs: loose)
    assert int(strict_count["d"]) == 3
    assert int(loose_count["d"]) == 0


def test_same_seed_reproduces_thinned_removal_and_a_different_seed_does_not():
    matrix = np.zeros((50, 3), dtype=np.int64)
    matrix[:, 0] = 200
    matrix[:, 1] = 1
    adata = AnnData(matrix, var=pd.DataFrame(index=["G0", "MT-CO1", "G2"]))
    adatas = {"s1": adata, "s2": adata.copy()}
    donor_of = pd.Series({"s1": "d1", "s2": "d2"})
    # Median depth is 201. Thinning to 50 leaves totals around 50, so a cap of 45
    # flags a large share of cells and that share depends on the seed.
    thresholds = QcThresholds(0.9, -1, 1e9, 45)
    kwargs = {
        "adatas": adatas,
        "specimen_to_donor": donor_of,
        "thresholds": thresholds,
        "target": 50,
    }
    first = depth_removed_counts(**kwargs, seeds=[0, 1])
    second = depth_removed_counts(**kwargs, seeds=[0, 1])
    pd.testing.assert_frame_equal(first, second)
    assert not first.loc[0].equals(first.loc[1])


def test_cache_loader_reads_the_specimen_stem(tmp_path):
    kept = AnnData(np.ones((5, 2), dtype=np.int64))
    decoy = AnnData(np.ones((8, 2), dtype=np.int64))
    kept.write_h5ad(tmp_path / "sample.h5ad")
    decoy.write_h5ad(tmp_path / "sample.h5ad.h5ad")
    loaded = load_specimen_adatas(["sample"], tmp_path)
    assert loaded["sample"].n_obs == 5
    with pytest.raises(FileNotFoundError):
        load_specimen_adatas(["missing"], tmp_path)


def test_robustness_verdicts():
    shared = {"distance_points": 1.0, "mde_points": 10.0}
    assert (
        judge_rerun(
            main_significant=False,
            main_statistic=0.1,
            rerun_significant=False,
            rerun_statistic=0.2,
            **shared,
        )
        == "pass"
    )
    assert (
        judge_rerun(
            main_significant=False,
            main_statistic=0.1,
            rerun_significant=True,
            rerun_statistic=0.2,
            **shared,
        )
        == "fail"
    )
    assert (
        judge_rerun(
            main_significant=True,
            main_statistic=0.2,
            rerun_significant=True,
            rerun_statistic=-0.2,
            distance_points=0.0,
            mde_points=10.0,
        )
        == "fail"
    )
    assert (
        judge_rerun(
            main_significant=False,
            main_statistic=0.1,
            rerun_significant=False,
            rerun_statistic=0.1,
            distance_points=9.0,
            mde_points=10.0,
        )
        == "fail"
    )
    assert (
        judge_rerun(
            main_significant=False,
            main_statistic=0.1,
            rerun_significant=False,
            rerun_statistic=0.1,
            distance_points=1.0,
            mde_points=10.0,
            seed_spread_points=6.0,
        )
        == "cannot_be_judged"
    )
    assert (
        judge_rerun(
            main_significant=False,
            main_statistic=0.1,
            rerun_significant=False,
            rerun_statistic=0.1,
            distance_points=1.0,
            mde_points=10.0,
            seed_spread_points=4.0,
        )
        == "pass"
    )
    assert overall_robustness("pass", "pass") == "robust"
    assert overall_robustness("pass", "fail") == "not_robust"
    assert overall_robustness("fail", "cannot_be_judged") == "not_robust"
    assert overall_robustness("pass", "cannot_be_judged") == "cannot_be_judged"


def test_missing_specimen_raises_before_a_genotype_test():
    cells, donors, _adatas = _cohort()
    with pytest.raises(ValueError, match="No AnnData"):
        run_filtering_rate(cells, donors, {}, _phase1(), n_bootstrap=0, n_resamples=0)


def test_secondary_adjustment_uniform_rule_and_handoff(tmp_path):
    cells, donors, adatas = _cohort()
    phase1 = _phase1()
    result = run_filtering_rate(
        cells,
        donors,
        adatas,
        phase1,
        ci_step=0.05,
        mde_step=0.5,
        n_bootstrap=20,
        n_resamples=0,
        rng=np.random.default_rng(0),
    )
    assert "q_value" not in result["primary"].columns
    assert "p_value" not in result["descriptive"].columns
    assert list(result["secondary"]["measure"]) == [
        "low_genes",
        "high_genes",
        "high_mito",
        "high_total",
    ]
    assert result["secondary"]["q_value"].to_numpy() == pytest.approx(
        bh_adjust(result["secondary"]["p_value"])
    )

    measures = donor_measures(cells)
    low = measures["n_low_genes"].to_numpy() / measures["n_cells"].to_numpy()
    direct = exact_permutation_test(
        low, measures["genotype"].to_numpy(), positive="BRCA1", alpha=0.02
    )
    stored = result["secondary"].set_index("measure").loc["low_genes"]
    assert stored["p_value"] == pytest.approx(direct.p_value)
    assert stored["estimate"] == pytest.approx(direct.statistic)
    assert stored["n_permutations"] == math.comb(6, 3)
    assert result["primary"].iloc[0]["n_permutations"] == math.comb(6, 3)

    rates = result["sensitivity_rates"]
    expected_uniform = {
        donor: n_uniform / n_cells for donor, _, _, _, _, n_cells, _, _, n_uniform in _COHORT
    }
    assert rates["uniform_fraction"].to_dict() == pytest.approx(expected_uniform)
    assert not np.allclose(rates["uniform_fraction"], rates["published_fraction"])

    rho = result["spearman"].set_index("rule").loc["uniform", "rho"]
    pearson = pearsonr(rates["published_fraction"], rates["uniform_fraction"]).statistic
    assert rho < -0.5
    assert rho != pytest.approx(pearson)

    handoff = write_filtering_rate_outputs(result, cells, donors, phase1, tmp_path)
    assert handoff == (tmp_path / "01_handoff.md").read_text(encoding="utf-8")
    assert "percentage points of cells" in handoff
    assert "expit percentage points" in handoff
    assert "Inconclusive" in handoff
    assert "N-0093 is not in this cohort" in handoff
    assert "B2 (BRCA1): Post (oophorectomy)" in handoff
    assert (tmp_path / "01-removal-fractions.png").stat().st_size > 0
    for name in (
        "01_primary.csv",
        "01_secondary.csv",
        "01_descriptive.csv",
        "01_sensitivity_rates.csv",
        "01_robustness.csv",
        "01_loo.csv",
        "01_spearman.csv",
    ):
        assert (tmp_path / name).is_file()
