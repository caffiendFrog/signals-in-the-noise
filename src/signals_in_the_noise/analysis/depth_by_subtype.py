"""Notebook 04's depth check, repeated for every cancer type against Normal.

Each contrast thins to its own shallowest included library and refits the PBS
classifier under the MAD noise rule, before and after thinning. A type with
fewer than two specimens is recorded and skipped.
"""

import logging
from collections.abc import Mapping

import numpy as np
import pandas as pd
from anndata import AnnData

from signals_in_the_noise.analysis.confounders import association_table, covariates_from_obs
from signals_in_the_noise.analysis.noise_phenotypes import ISCB_SCS_PBS_THRESHOLDS
from signals_in_the_noise.modeling.datasets import CORE_ARM_NAMES, build_cell_table, standard_arms
from signals_in_the_noise.modeling.experiments import (
    PUBLISHED_PROTOCOL,
    InsufficientSpecimensError,
    reference_specimens,
    run_arm_suite,
)
from signals_in_the_noise.preprocessing.qc import (
    depth_matched_qc_obs,
    library_depth,
    mad_thresholds,
    relabel_noise,
)

logger = logging.getLogger(__name__)

REFERENCE_CONDITION = "Normal"
NOISE_COLUMN = "is_noise_mad"
PBS_ARM = "noise_pbs"
DEPTH_COVARIATE = "median_log1p_total_counts_retained"
MIN_PER_CONDITION = 2
READING_ALPHA = 0.05


def other_conditions(
    conditions: Mapping[str, str], reference: str = REFERENCE_CONDITION
) -> list[str]:
    """Cancer types other than ``reference``, in alphabetical order."""
    labels = set(conditions.values())
    if reference not in labels:
        raise ValueError(f"No specimens are labelled {reference!r}.")
    return sorted(labels - {reference})


def depth_reading(p_original: float, p_matched: float, alpha: float = READING_ALPHA) -> str:
    """How the PBS specimen-level signal responds to equalizing depth.

    ``alpha`` labels the existing permutation p-values. It is not a new test.
    """
    if not np.isfinite(p_original) or not np.isfinite(p_matched):
        return "not compared"
    if p_original > alpha:
        return "no specimen-level PBS separation"
    if p_matched > alpha:
        return "collapses after thinning"
    return "remains after thinning"


def run_depth_by_subtype(
    specimens: Mapping[str, AnnData],
    conditions: Mapping[str, str],
    patients: Mapping[str, str],
    *,
    reference: str = REFERENCE_CONDITION,
    n_mads: float = 3.0,
    n_permutations: int = 200,
    seed: int = 0,
    min_per_condition: int = MIN_PER_CONDITION,
) -> pd.DataFrame:
    """One row per cancer type: PBS separation before and after depth matching."""
    rows = []
    for positive in other_conditions(conditions, reference):
        kept = [sid for sid, label in conditions.items() if label in {positive, reference}]
        subset = {sid: specimens[sid] for sid in kept}
        labels = {sid: conditions[sid] for sid in kept}
        groups = {sid: patients[sid] for sid in kept}
        logger.info("Depth contrast %s vs %s (%s specimens).", positive, reference, len(subset))
        rows.append(
            _one_contrast(
                subset,
                labels,
                groups,
                positive=positive,
                reference=reference,
                n_mads=n_mads,
                n_permutations=n_permutations,
                seed=seed,
                min_per_condition=min_per_condition,
            )
        )
    return pd.DataFrame(rows)


def depth_by_subtype_handoff(result: pd.DataFrame) -> str:
    """Markdown summary of the per-type depth check."""
    lines = [
        "# Depth check across cancer types",
        "",
        "Each row repeats notebook 04: MAD noise labels, binomial thinning to the",
        "shallowest included library in that contrast, and the PBS classifier of",
        f"discarded cells against {REFERENCE_CONDITION}. The reading uses the specimen",
        f"permutation p-value at {READING_ALPHA:g} only as a label.",
        "",
        _markdown_table(_handoff_table(result)),
        "",
        "A collapse after thinning means the PBS separation in that contrast moves",
        "into the permutation null once library depth is equalized. A contrast that",
        "remains is not explained by depth in the way ER+ was.",
        "",
    ]
    return "\n".join(lines)


def _one_contrast(
    specimens,
    conditions,
    patients,
    *,
    positive,
    reference,
    n_mads,
    n_permutations,
    seed,
    min_per_condition,
) -> dict:
    n_positive = sum(label == positive for label in conditions.values())
    n_reference = sum(label == reference for label in conditions.values())
    base = {
        "condition": positive,
        "reference": reference,
        "n_positive": n_positive,
        "n_reference": n_reference,
    }
    if n_positive < min_per_condition or n_reference < min_per_condition:
        return _skipped(base, f"fewer than {min_per_condition} specimens in a condition")

    obs = {sid: adata.obs for sid, adata in specimens.items()}
    arms = standard_arms(ISCB_SCS_PBS_THRESHOLDS, names=["noise_pbs"])
    reference_table = build_cell_table(obs, conditions, arms[0], patients=patients)
    included = reference_specimens(reference_table, PUBLISHED_PROTOCOL["min_reference_cells"])
    included_labels = [conditions[sid] for sid in included]
    n_positive_included = included_labels.count(positive)
    n_reference_included = included_labels.count(reference)
    base["n_positive"] = n_positive_included
    base["n_reference"] = n_reference_included
    if n_positive_included < min_per_condition or n_reference_included < min_per_condition:
        return _skipped(
            base, f"fewer than {min_per_condition} specimens with enough PBS-labelled cells"
        )

    retained = {sid: specimens[sid].obs["is_noise"].to_numpy() == 0 for sid in specimens}
    depths = pd.Series(
        {sid: library_depth(specimens[sid], retained[sid]) for sid in specimens}, dtype=float
    )
    target = float(depths.loc[included].quantile(0.0))
    rng = np.random.default_rng(seed)
    thinned = {}
    for sid in sorted(specimens):
        thinned[sid], _fraction = depth_matched_qc_obs(
            specimens[sid], target, rng, reference_mask=retained[sid]
        )

    def mad_rule(frame):
        return mad_thresholds(frame, n_mads=n_mads)

    versions = {
        "original": relabel_noise(obs, mad_rule, NOISE_COLUMN),
        "depth-matched": relabel_noise(thinned, mad_rule, NOISE_COLUMN),
    }
    fitted = {}
    depth_auc = {}
    for version, version_obs in versions.items():
        try:
            fitted[version] = run_arm_suite(
                version_obs,
                conditions,
                standard_arms(ISCB_SCS_PBS_THRESHOLDS, names=CORE_ARM_NAMES),
                positive_label=positive,
                negative_label=reference,
                patients=patients,
                noise_column=NOISE_COLUMN,
                n_permutations=n_permutations,
                seed=seed,
                **PUBLISHED_PROTOCOL,
            )
        except InsufficientSpecimensError as error:
            return _skipped(base, str(error))
        covs = covariates_from_obs(
            {sid: version_obs[sid] for sid in included}, conditions, NOISE_COLUMN
        )
        checked = association_table(
            covs, positive=positive, negative=reference, columns=[DEPTH_COVARIATE]
        )
        depth_auc[version] = float(checked.loc[DEPTH_COVARIATE, "auc"])

    original = _arm_metrics(fitted["original"], PBS_ARM)
    matched = _arm_metrics(fitted["depth-matched"], PBS_ARM)
    return {
        **base,
        "target_depth": target,
        "status": "ran",
        "skip_reason": "",
        "noise_pbs_auc_original": original["auc"],
        "noise_pbs_p_original": original["p"],
        "noise_pbs_auc_depth_matched": matched["auc"],
        "noise_pbs_p_depth_matched": matched["p"],
        "depth_auc_original": depth_auc["original"],
        "depth_auc_depth_matched": depth_auc["depth-matched"],
        "reading": depth_reading(original["p"], matched["p"]),
    }


def _arm_metrics(result, arm: str) -> dict:
    summary = result.summary()
    return {
        "auc": float(summary.loc[arm, "specimen_auc"]),
        "p": float(summary.loc[arm, "permutation_p"]),
    }


def _skipped(base: dict, reason: str) -> dict:
    return {
        **base,
        "target_depth": np.nan,
        "status": "skipped",
        "skip_reason": reason,
        "noise_pbs_auc_original": np.nan,
        "noise_pbs_p_original": np.nan,
        "noise_pbs_auc_depth_matched": np.nan,
        "noise_pbs_p_depth_matched": np.nan,
        "depth_auc_original": np.nan,
        "depth_auc_depth_matched": np.nan,
        "reading": "not compared",
    }


def _handoff_table(result: pd.DataFrame) -> pd.DataFrame:
    shown = result[
        [
            "condition",
            "n_positive",
            "n_reference",
            "status",
            "noise_pbs_p_original",
            "noise_pbs_p_depth_matched",
            "reading",
        ]
    ].copy()
    for column in ("noise_pbs_p_original", "noise_pbs_p_depth_matched"):
        shown[column] = shown[column].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    return shown


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(str(row[column]) for column in frame.columns) + " |"
        for _, row in frame.iterrows()
    ]
    return "\n".join([header, separator, *body])
