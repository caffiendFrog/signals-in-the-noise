"""Specimen-level technical covariates and their association with condition.

These quantify the sample-level differences (depth, library complexity, QC
thresholds, noise burden) that could let a classifier separate conditions
without any cell-level biology.
"""

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import false_discovery_control, fisher_exact, mannwhitneyu, spearmanr

from signals_in_the_noise.preprocessing.qc import (
    POPULATIONS,
    QC_METRICS,
    QcThresholds,
    select_population,
)
from signals_in_the_noise.preprocessing.specimens import (
    CONDITION_COLUMN,
    PATIENT_COLUMN,
    SPECIMEN_COLUMN,
    patient_id,
    specimen_conditions,
)

PRESPECIFIED_TECHNICAL_COVARIATES: tuple[str, ...] = (
    "median_log1p_total_counts_retained",
    "median_complexity_retained",
    "median_pct_counts_mt_all",
    "qc_genes_lower",
)
"""Covariates fixed before seeing results: sequencing depth, library complexity,
a dissociation-stress proxy, and Pal's per-sample gene threshold."""


def specimen_qc_summary(obs: pd.DataFrame, noise_column: str = "is_noise") -> dict[str, float]:
    """Summarize one specimen's QC distributions.

    Returns cell count, noise fraction, per-population medians of each metric
    in :data:`~signals_in_the_noise.preprocessing.qc.QC_METRICS`, and the
    median library complexity (``log1p genes / log1p counts``) of retained cells.
    """
    summary: dict[str, float] = {
        "n_cells": float(len(obs)),
        "noise_fraction": float(obs[noise_column].mean()) if len(obs) else np.nan,
    }
    for population in POPULATIONS:
        cells = select_population(obs, population, noise_column)
        for metric in QC_METRICS:
            summary[f"median_{metric}_{population}"] = float(cells[metric].median())
    retained = select_population(obs, "retained", noise_column)
    complexity = retained["log1p_n_genes_by_counts"] / retained["log1p_total_counts"]
    summary["median_complexity_retained"] = float(complexity.median())
    return summary


def covariates_from_obs(
    specimen_obs: Mapping[str, pd.DataFrame],
    conditions: Mapping[str, str],
    noise_column: str = "is_noise",
) -> pd.DataFrame:
    """Build a specimen-indexed table of QC summaries plus a ``condition`` column."""
    rows = {
        specimen_id: {
            CONDITION_COLUMN: conditions[specimen_id],
            **specimen_qc_summary(obs, noise_column),
        }
        for specimen_id, obs in specimen_obs.items()
    }
    table = pd.DataFrame.from_dict(rows, orient="index")
    table.index.name = SPECIMEN_COLUMN
    return table


def specimen_covariates(
    adatas: Mapping[str, AnnData], noise_column: str = "is_noise"
) -> pd.DataFrame:
    """Build the full technical-covariate table for GSE161529-style specimens.

    Adds the patient identifier, the published per-sample QC thresholds, the
    number of genes detected (``uns['num_genes_before']``), menopause status, and
    ``lowered_genes_lower`` — whether the specimen's gene threshold is below
    the study-wide mode (Pal's low-coverage adjustment).
    """
    table = covariates_from_obs(
        {sid: adata.obs for sid, adata in adatas.items()}, specimen_conditions(adatas), noise_column
    )
    for sid, adata in adatas.items():
        table.loc[sid, PATIENT_COLUMN] = patient_id(adata)
        for key in QcThresholds.UNS_KEYS.values():
            table.loc[sid, key] = float(adata.uns[key])
        table.loc[sid, "num_genes_detected"] = float(adata.uns["num_genes_before"])
        table.loc[sid, "menopause_status"] = str(adata.uns["menopause_status"])
    modal_genes_lower = table["qc_genes_lower"].mode().min()
    table["lowered_genes_lower"] = (table["qc_genes_lower"] < modal_genes_lower).astype(int)
    return table


def feature_covariate_correlations(
    table: pd.DataFrame,
    features: Iterable[str],
    covariates: Iterable[str],
    *,
    group_column: str | None = CONDITION_COLUMN,
) -> pd.DataFrame:
    """Spearman correlation of each feature with each covariate, overall and per group.

    A feature that tracks a technical covariate *within* each condition is
    driven by that covariate regardless of any condition difference.

    Returns:
        Long table with ``feature``, ``covariate``, ``group`` (``"all"`` or a
        condition), ``n``, ``rho`` and ``p_value``.
    """
    subsets = [("all", table)]
    if group_column is not None:
        subsets += list(table.groupby(group_column))
    rows = []
    for group, subset in subsets:
        for feature in features:
            for covariate in covariates:
                pair = subset[[feature, covariate]].dropna()
                if len(pair) < 3 or pair.nunique().min() < 2:
                    rho, p_value = np.nan, np.nan
                else:
                    result = spearmanr(pair[feature], pair[covariate])
                    rho, p_value = float(result.statistic), float(result.pvalue)
                rows.append(
                    {
                        "feature": feature,
                        "covariate": covariate,
                        "group": group,
                        "n": len(pair),
                        "rho": rho,
                        "p_value": p_value,
                    }
                )
    return pd.DataFrame(rows)


def _is_binary(values: pd.Series) -> bool:
    return set(pd.unique(values.dropna())) <= {0, 1, True, False}


def association_table(
    table: pd.DataFrame,
    *,
    positive: str,
    negative: str,
    columns: Iterable[str] | None = None,
    group_column: str = CONDITION_COLUMN,
) -> pd.DataFrame:
    """Test each specimen covariate for a difference between two conditions.

    Binary covariates use Fisher's exact test and report the rate per group;
    continuous ones use a two-sided Mann–Whitney U test and report the median
    per group (columns named after ``positive`` and ``negative``). ``auc`` is
    ``P(positive > negative)`` from the U statistic (0.5 = no separation,
    0 or 1 = perfect separation). ``q_value``
    is Benjamini–Hochberg adjusted across the tested covariates.

    Returns:
        DataFrame indexed by covariate, sorted by ``p_value``.
    """
    if columns is None:
        columns = [c for c in table.select_dtypes("number").columns if c != group_column]
    pos = table[table[group_column] == positive]
    neg = table[table[group_column] == negative]

    rows = {}
    for column in columns:
        a = pos[column].dropna().astype(float)
        b = neg[column].dropna().astype(float)
        binary = _is_binary(table[column])
        if len(a) == 0 or len(b) == 0 or pd.concat([a, b]).nunique() <= 1:
            auc, p_value = 0.5, 1.0
        else:
            u_stat, mwu_p = mannwhitneyu(a, b, alternative="two-sided")
            auc = float(u_stat / (len(a) * len(b)))
            if binary:
                contingency = [[a.sum(), len(a) - a.sum()], [b.sum(), len(b) - b.sum()]]
                p_value = float(fisher_exact(contingency).pvalue)
            else:
                p_value = float(mwu_p)
        rows[column] = {
            "center": "rate" if binary else "median",
            positive: float(a.mean() if binary else a.median()),
            negative: float(b.mean() if binary else b.median()),
            "auc": auc,
            "test": "fisher" if binary else "mann-whitney",
            "p_value": p_value,
        }

    result = pd.DataFrame.from_dict(rows, orient="index")
    if not result.empty:
        result["q_value"] = false_discovery_control(result["p_value"].to_numpy())
        result = result.sort_values("p_value")
    return result
