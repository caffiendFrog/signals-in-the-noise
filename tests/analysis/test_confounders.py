"""Tests for signals_in_the_noise.analysis.confounders."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.confounders import (
    PRESPECIFIED_TECHNICAL_COVARIATES,
    association_table,
    covariates_from_obs,
    feature_covariate_correlations,
    specimen_covariates,
    specimen_qc_summary,
)


def _obs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pct_counts_mt": [10.0, 20.0, 1.0, 3.0],
            "log1p_total_counts": [4.0, 6.0, 8.0, 10.0],
            "log1p_n_genes_by_counts": [2.0, 3.0, 6.0, 5.0],
            "is_noise": [1, 1, 0, 0],
        }
    )


def test_specimen_qc_summary_values():
    summary = specimen_qc_summary(_obs())
    assert summary["n_cells"] == 4
    assert summary["noise_fraction"] == 0.5
    assert summary["median_log1p_total_counts_all"] == 7.0
    assert summary["median_log1p_total_counts_retained"] == 9.0
    assert summary["median_log1p_total_counts_noise"] == 5.0
    assert summary["median_pct_counts_mt_noise"] == 15.0
    assert summary["median_complexity_retained"] == pytest.approx(np.median([6 / 8, 5 / 10]))


def test_specimen_qc_summary_respects_alternative_noise_column():
    obs = _obs().assign(is_noise_alt=[0, 0, 0, 1])
    assert specimen_qc_summary(obs, "is_noise_alt")["noise_fraction"] == 0.25


def test_covariates_from_obs_builds_specimen_indexed_table():
    table = covariates_from_obs({"s1": _obs(), "s2": _obs()}, {"s1": "Normal", "s2": "ER+ tumour"})
    assert table.index.name == "specimen_id"
    assert table.loc["s2", "condition"] == "ER+ tumour"
    assert "median_complexity_retained" in table.columns


def _adata(cancer_type: str, genes_lower: float) -> AnnData:
    adata = AnnData(obs=_obs())
    adata.uns.update(
        {
            "title": f"{cancer_type} Total cells from Patient {int(genes_lower):04d}",
            "cancer_type": cancer_type,
            "qc_mito_upper": 0.2,
            "qc_genes_lower": genes_lower,
            "qc_genes_upper": 6000,
            "qc_total_upper": 40000,
            "num_genes_before": 20000,
            "menopause_status": "Pre",
        }
    )
    return adata


def test_specimen_covariates_adds_thresholds_and_lowered_flag():
    table = specimen_covariates(
        {"a": _adata("Normal", 500), "b": _adata("Normal", 500), "c": _adata("ER+ tumour", 300)}
    )
    assert table["qc_genes_lower"].tolist() == [500, 500, 300]
    assert table["lowered_genes_lower"].tolist() == [0, 0, 1]
    assert table.loc["a", "num_genes_detected"] == 20000
    assert table.loc["c", "menopause_status"] == "Pre"
    assert table.loc["c", "patient_id"] == "0300"
    assert set(PRESPECIFIED_TECHNICAL_COVARIATES) <= set(table.columns)


def _association_input() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "condition": ["ER"] * 6 + ["N"] * 6,
            "shifted": list(range(10, 16)) + list(range(6)),
            "unrelated": [1, 5, 3, 2, 6, 4] * 2,
            "constant": [7.0] * 12,
            "binary": [1, 1, 1, 1, 1, 0] + [0, 0, 0, 0, 0, 1],
        }
    )


def test_association_table_detects_separated_covariate():
    result = association_table(_association_input(), positive="ER", negative="N")
    assert result.loc["shifted", "auc"] == 1.0
    assert result.loc["shifted", "p_value"] < 0.01
    assert result.loc["shifted", "ER"] == 12.5
    assert result.loc["shifted", "N"] == 2.5
    assert result.loc["unrelated", "auc"] == 0.5


def test_association_table_uses_fisher_and_rates_for_binary_columns():
    result = association_table(_association_input(), positive="ER", negative="N")
    assert result.loc["binary", "test"] == "fisher"
    assert result.loc["binary", "center"] == "rate"
    assert result.loc["binary", "ER"] == pytest.approx(5 / 6)
    assert result.loc["shifted", "test"] == "mann-whitney"


def test_association_table_constant_column_is_uninformative():
    result = association_table(_association_input(), positive="ER", negative="N")
    assert result.loc["constant", "p_value"] == 1.0
    assert result.loc["constant", "auc"] == 0.5


def test_association_table_q_values_are_adjusted_and_sorted():
    result = association_table(_association_input(), positive="ER", negative="N")
    assert (result["q_value"] >= result["p_value"] - 1e-12).all()
    assert result["p_value"].is_monotonic_increasing


def test_feature_covariate_correlations_overall_and_within_condition():
    table = pd.DataFrame(
        {
            "condition": ["ER"] * 5 + ["N"] * 5,
            "depth": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
            "pbs-1": [1, 2, 3, 4, 5, 5, 4, 3, 2, 1],
            "flat": [1.0] * 10,
        }
    )
    result = feature_covariate_correlations(table, ["pbs-1", "flat"], ["depth"]).set_index(
        ["feature", "group"]
    )
    assert result.loc[("pbs-1", "ER"), "rho"] == pytest.approx(1.0)
    assert result.loc[("pbs-1", "N"), "rho"] == pytest.approx(-1.0)
    assert result.loc[("pbs-1", "all"), "rho"] == pytest.approx(0.0)
    assert result.loc[("pbs-1", "ER"), "n"] == 5
    assert np.isnan(result.loc[("flat", "all"), "rho"])


def test_feature_covariate_correlations_without_groups():
    table = pd.DataFrame({"a": [1, 2, 3, 4], "b": [2, 4, 6, 8], "condition": ["x"] * 4})
    result = feature_covariate_correlations(table, ["a"], ["b"], group_column=None)
    assert result["group"].tolist() == ["all"]


def test_association_table_respects_explicit_columns():
    result = association_table(
        _association_input(), positive="ER", negative="N", columns=["shifted"]
    )
    assert list(result.index) == ["shifted"]
