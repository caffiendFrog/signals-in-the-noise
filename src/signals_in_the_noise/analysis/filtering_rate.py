"""Step 1 of the BRCA1 QC-filtering analysis: does removal differ by genotype?

The donor is the unit. The primary test is the overall removed fraction. The
four non-exclusive reasons are a secondary family. Exclusive reasons and cells
that fail several filters are described, not tested. Uniform, MAD, and
depth-equalized re-runs are stability checks of the overall result.
"""

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from signals_in_the_noise.analysis.donor_permutation import (
    attainable_alpha,
    bh_adjust,
    empirical_log_odds,
    exact_permutation_test,
    leave_one_donor_out,
    log_odds_shift_points,
    mde_resample_spread,
    mean_diff,
    minimum_detectable_effect,
    permutation_ci,
)
from signals_in_the_noise.analysis.qc_filtering import (
    PARAMS_FILENAME,
    output_directory,
)
from signals_in_the_noise.preprocessing.qc import (
    EXCLUSIVE_REASON_PRIORITY,
    NOISE_FLAG_COLUMNS,
    QcThresholds,
    depth_matched_qc_obs,
    flag_noise,
    mad_thresholds,
    relabel_noise,
)
from signals_in_the_noise.utils.log import get_logger
from signals_in_the_noise.utils.visualization import get_figure_axes

logger = get_logger(__name__)

POSITIVE = "BRCA1"
DISPLAY_REASONS = ("low_genes", "high_genes", "high_mito", "high_total")
FLAG_OF_REASON = {label: column for column, label in EXCLUSIVE_REASON_PRIORITY}
REASON_LABELS = {
    "low_genes": "Low gene count",
    "high_genes": "High gene count",
    "high_mito": "High mitochondrial fraction",
    "high_total": "High library size",
}
MENOPAUSE_COLORS = {"Pre": "#0072B2", "Post (oophorectomy)": "#D55E00"}
SOURCE_MARKERS = {"PM": "o", "MH": "s", "KCF": "D"}

PRIMARY_FILENAME = "01_primary.csv"
SECONDARY_FILENAME = "01_secondary.csv"
DESCRIPTIVE_FILENAME = "01_descriptive.csv"
SENSITIVITY_FILENAME = "01_sensitivity_rates.csv"
ROBUSTNESS_FILENAME = "01_robustness.csv"
LOO_FILENAME = "01_loo.csv"
SPEARMAN_FILENAME = "01_spearman.csv"
HANDOFF_FILENAME = "01_handoff.md"
FIGURE_FILENAME = "01-removal-fractions.png"


def donor_measures(cells: pd.DataFrame) -> pd.DataFrame:
    """Per-donor counts for the overall rate, each reason, and cells failing several.

    Fractions use all cells of the donor as the denominator. Non-exclusive
    counts are the four flags. Exclusive counts use the single assigned reason.
    """
    required = {"donor", "genotype", "is_noise", "exclusive_reason", *NOISE_FLAG_COLUMNS}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"Cell table is missing columns {sorted(missing)}.")
    if cells.groupby("donor")["genotype"].nunique().max() != 1:
        raise ValueError("A donor has more than one genotype.")
    rows = []
    for donor, group in cells.groupby("donor", sort=True):
        n_cells = len(group)
        flags = group[list(NOISE_FLAG_COLUMNS)].astype(int)
        row = {
            "donor": donor,
            "genotype": group["genotype"].iloc[0],
            "n_cells": n_cells,
            "n_removed": int(group["is_noise"].astype(int).sum()),
            "n_multi": int((flags.sum(axis=1) >= 2).sum()),
        }
        for label, column in FLAG_OF_REASON.items():
            row[f"n_{label}"] = int(flags[column].sum())
            exclusive = group["exclusive_reason"].astype("string").eq(label)
            row[f"n_exclusive_{label}"] = int(exclusive.sum())
        rows.append(row)
    measures = pd.DataFrame(rows).set_index("donor")
    measures["removed_fraction"] = measures["n_removed"] / measures["n_cells"]
    measures["fraction_multi"] = measures["n_multi"] / measures["n_cells"]
    for label in FLAG_OF_REASON:
        measures[f"fraction_{label}"] = measures[f"n_{label}"] / measures["n_cells"]
        measures[f"fraction_exclusive_{label}"] = (
            measures[f"n_exclusive_{label}"] / measures["n_cells"]
        )
    return measures


def load_specimen_adatas(specimen_ids, cache_dir: Path) -> dict[str, AnnData]:
    """Load the cohort cache files. Keys are specimen ids, without ``.h5ad``."""
    cache_dir = Path(cache_dir)
    loaded = {}
    for specimen_id in specimen_ids:
        path = cache_dir / f"{specimen_id}.h5ad"
        if not path.is_file():
            raise FileNotFoundError(f"Missing cache file {path}.")
        loaded[str(specimen_id)] = sc.read_h5ad(path)
    return loaded


def frozen_uniform_thresholds(phase1: dict) -> QcThresholds:
    """Uniform rule frozen by notebook 00. Not recomputed from the samples."""
    raw = phase1["sensitivity"]["uniform_thresholds"]
    return QcThresholds(
        mito_upper=float(raw["mito_upper"]),
        genes_lower=float(raw["genes_lower"]),
        genes_upper=float(raw["genes_upper"]),
        total_upper=float(raw["total_upper"]),
    )


def removed_counts_under_rule(specimen_obs, specimen_to_donor, rule) -> pd.Series:
    """Removed-cell counts after applying ``rule`` to each specimen's QC table."""
    relabelled = relabel_noise(specimen_obs, rule, "is_noise_rerun")
    counts = {}
    for specimen_id, obs in relabelled.items():
        donor = specimen_to_donor[specimen_id]
        if donor in counts:
            raise ValueError(f"Donor {donor} is attached to more than one specimen.")
        counts[donor] = int(obs["is_noise_rerun"].sum())
    return pd.Series(counts, dtype=int)


def depth_removed_counts(adatas, specimen_to_donor, thresholds, target, seeds) -> pd.DataFrame:
    """Removed counts after thinning each specimen to ``target`` under ``thresholds``.

    Rows are seeds and columns are donors. Each seed uses ``default_rng(seed)``,
    so the same seed repeats. Specimens already at or below the target are not
    thinned. The QC rule is the frozen uniform rule, not a new mode.
    """
    rows = []
    for seed in seeds:
        counts = {}
        for specimen_id, donor in specimen_to_donor.items():
            obs, _fraction = depth_matched_qc_obs(
                adatas[specimen_id], target, np.random.default_rng(int(seed))
            )
            counts[donor] = int(flag_noise(obs, thresholds)["is_noise"].sum())
        counts["seed"] = int(seed)
        rows.append(counts)
        logger.info("Depth-equalized seed %s finished.", seed)
    return pd.DataFrame(rows).set_index("seed")


def judge_rerun(
    *,
    main_significant: bool,
    main_statistic: float,
    rerun_significant: bool,
    rerun_statistic: float,
    distance_points: float,
    mde_points: float,
    seed_spread_points: float | None = None,
) -> str:
    """``pass``, ``fail``, or ``cannot_be_judged`` for one stability re-run.

    Pass requires the same conclusion at the primary alpha, including the
    direction when the main result is significant, and a log-odds distance of
    at most half the primary MDE. A depth re-run whose seeds spread by more
    than that half is ``cannot_be_judged``.
    """
    if not np.isfinite(mde_points):
        return "cannot_be_judged"
    half = mde_points / 2
    if seed_spread_points is not None and seed_spread_points > half:
        return "cannot_be_judged"
    same_conclusion = main_significant == rerun_significant and (
        not main_significant or _sign(main_statistic) == _sign(rerun_statistic)
    )
    if same_conclusion and abs(distance_points) <= half:
        return "pass"
    return "fail"


def overall_robustness(uniform_status: str, depth_status: str) -> str:
    """The overall result is robust only when both re-runs pass."""
    if uniform_status == "pass" and depth_status == "pass":
        return "robust"
    if "fail" in (uniform_status, depth_status):
        return "not_robust"
    return "cannot_be_judged"


def run_filtering_rate(
    cells: pd.DataFrame,
    donors: pd.DataFrame,
    adatas: dict[str, AnnData],
    phase1: dict,
    *,
    ci_step: float = 0.001,
    mde_step: float = 0.05,
    n_bootstrap: int = 2000,
    n_resamples: int = 20,
    rng: np.random.Generator | None = None,
) -> dict[str, pd.DataFrame]:
    """Primary, secondary, descriptive, sensitivity, and robustness tables.

    Confidence intervals and the minimum detectable effect use the alpha,
    confidence levels, thresholds, depth target, and seeds stored in ``phase1``.
    """
    rng = np.random.default_rng() if rng is None else rng
    measures = donor_measures(cells)
    specimen = _specimen_table(donors)
    _check_alignment(measures, specimen, adatas)
    labels = measures["genotype"].to_numpy()
    totals = measures["n_cells"].to_numpy(dtype=float)
    alpha = phase1["alpha"]
    level = phase1["confidence"]
    positive = POSITIVE

    logger.info("Primary removed-fraction test.")
    primary = _tested_row(
        "overall_removed_fraction",
        measures["n_removed"].to_numpy(dtype=float),
        totals,
        labels,
        alpha=alpha["primary"],
        alternative="two-sided",
        level=level["primary"],
        ci_step=ci_step,
        n_bootstrap=n_bootstrap,
        rng=rng,
        mde_step=mde_step,
    )
    if n_resamples:
        logger.info("Resampling donors within genotype for the primary MDE.")
        spread = mde_resample_spread(
            measures["n_removed"].to_numpy(dtype=float),
            totals,
            labels,
            positive=positive,
            alpha=alpha["primary"],
            alternative="two-sided",
            n_resamples=n_resamples,
            rng=rng,
            delta_step=max(mde_step, 0.1),
        )
        primary.update({f"mde_resample_{name}": value for name, value in spread.items()})

    logger.info("Secondary reason tests.")
    secondary = [
        _tested_row(
            reason,
            measures[f"n_{reason}"].to_numpy(dtype=float),
            totals,
            labels,
            alpha=alpha["secondary"],
            alternative="two-sided",
            level=level["secondary"],
            ci_step=ci_step,
            n_bootstrap=n_bootstrap,
            rng=rng,
            mde_step=mde_step,
        )
        for reason in DISPLAY_REASONS
    ]
    secondary = pd.DataFrame(secondary)
    secondary.insert(1, "q_value", bh_adjust(secondary["p_value"]))

    descriptive_specs = [
        (f"exclusive_{reason}", f"n_exclusive_{reason}") for reason in DISPLAY_REASONS
    ]
    descriptive_specs.append(("multi_reason", "n_multi"))
    descriptive = pd.DataFrame(
        [
            _interval_row(
                name,
                measures[column].to_numpy(dtype=float),
                totals,
                labels,
                level=level["secondary"],
                ci_step=ci_step,
            )
            for name, column in descriptive_specs
        ]
    )

    thresholds = frozen_uniform_thresholds(phase1)
    donor_of = specimen.set_index("specimen_id")["donor"]
    specimen_obs = {sid: adatas[sid].obs for sid in donor_of.index}
    uniform_removed = removed_counts_under_rule(
        specimen_obs, donor_of, lambda _obs: thresholds
    ).reindex(measures.index)
    mad_n = float(phase1["sensitivity"]["mad"]["n_mads"])
    mad_removed = removed_counts_under_rule(
        specimen_obs, donor_of, lambda obs: mad_thresholds(obs, n_mads=mad_n)
    ).reindex(measures.index)
    depth = phase1["sensitivity"]["depth_equalized"]
    if depth["qc_rule_after_thinning"] != "uniform":
        raise ValueError("Depth-equalized QC must use the frozen uniform rule.")
    seed_counts = depth_removed_counts(
        adatas, donor_of, thresholds, float(depth["target_median_umi"]), depth["seeds"]
    )
    seed_counts = seed_counts.reindex(columns=measures.index)
    depth_removed = seed_counts.mean(axis=0)

    sensitivity = measures[["genotype", "n_cells", "n_removed", "removed_fraction"]].copy()
    sensitivity = sensitivity.rename(
        columns={"n_removed": "published_removed", "removed_fraction": "published_fraction"}
    )
    sensitivity["uniform_removed"] = uniform_removed
    sensitivity["uniform_fraction"] = uniform_removed / sensitivity["n_cells"]
    sensitivity["mad_removed"] = mad_removed
    sensitivity["mad_fraction"] = mad_removed / sensitivity["n_cells"]
    sensitivity["depth_removed"] = depth_removed
    sensitivity["depth_fraction"] = depth_removed / sensitivity["n_cells"]
    for seed, row in seed_counts.iterrows():
        sensitivity[f"depth_seed_{seed}"] = row / sensitivity["n_cells"]

    pooled = float(primary["mde_pooled_log_odds"])
    mde_points = float(primary["mde_points"]) if primary["mde_attained"] else math.nan
    robustness_rows = []
    reruns = {
        "uniform": uniform_removed.to_numpy(dtype=float),
        "mad": mad_removed.to_numpy(dtype=float),
        "depth_equalized": depth_removed.to_numpy(dtype=float),
    }
    for rule, removed in reruns.items():
        test = exact_permutation_test(
            removed / totals,
            labels,
            positive=positive,
            alternative="two-sided",
            alpha=alpha["primary"],
        )
        distance = log_odds_shift_points(
            pooled,
            mean_diff(empirical_log_odds(removed, totals), labels, positive)
            - mean_diff(empirical_log_odds(measures["n_removed"], totals), labels, positive),
        )
        spread_points = None
        if rule == "depth_equalized":
            seed_effects = [
                mean_diff(
                    empirical_log_odds(seed_counts.loc[seed].to_numpy(dtype=float), totals),
                    labels,
                    positive,
                )
                for seed in seed_counts.index
            ]
            spread_points = abs(
                log_odds_shift_points(pooled, max(seed_effects) - min(seed_effects))
            )
        status = judge_rerun(
            main_significant=bool(primary["significant"]),
            main_statistic=float(primary["estimate"]),
            rerun_significant=test.significant,
            rerun_statistic=test.statistic,
            distance_points=distance,
            mde_points=mde_points,
            seed_spread_points=spread_points,
        )
        published_fraction = measures["removed_fraction"]
        rerun_fraction = pd.Series(removed / totals, index=measures.index)
        robustness_rows.append(
            {
                "rule": rule,
                "in_verdict": rule != "mad",
                "status": status,
                "estimate": test.statistic,
                "significant": test.significant,
                "distance_points": distance,
                "half_mde_points": mde_points / 2 if np.isfinite(mde_points) else np.nan,
                "seed_spread_points": spread_points if spread_points is not None else np.nan,
                "drivers": _driver_text(published_fraction, rerun_fraction)
                if status == "fail"
                else "",
            }
        )
    by_rule = {row["rule"]: row["status"] for row in robustness_rows}
    robustness_rows.append(
        {
            "rule": "overall",
            "in_verdict": True,
            "status": overall_robustness(by_rule["uniform"], by_rule["depth_equalized"]),
            "estimate": np.nan,
            "significant": pd.NA,
            "distance_points": np.nan,
            "half_mde_points": mde_points / 2 if np.isfinite(mde_points) else np.nan,
            "seed_spread_points": np.nan,
            "drivers": "",
        }
    )

    loo_rows = []
    loo_measures = [("overall_removed_fraction", measures["n_removed"].to_numpy(dtype=float))]
    loo_measures.extend(
        (reason, measures[f"n_{reason}"].to_numpy(dtype=float)) for reason in DISPLAY_REASONS
    )
    for name, removed in loo_measures:
        for row in leave_one_donor_out(
            removed / totals,
            labels,
            measures.index.to_numpy(),
            positive=positive,
            level=level["loo"],
            step=ci_step,
        ):
            row["measure"] = name
            row["estimate_points"] = 100 * row["estimate"]
            row["ci_width_points"] = 100 * (row["ci_high"] - row["ci_low"])
            loo_rows.append(row)

    spearman_rows = []
    published = measures["removed_fraction"]
    for rule, column in (
        ("uniform", "uniform_fraction"),
        ("mad", "mad_fraction"),
        ("depth_equalized", "depth_fraction"),
    ):
        other = sensitivity[column]
        if published.nunique() < 2 or other.nunique() < 2:
            rho = math.nan
        else:
            rho = float(spearmanr(published, other).statistic)
        spearman_rows.append({"rule": rule, "rho": rho})

    return {
        "primary": pd.DataFrame([primary]),
        "secondary": secondary,
        "descriptive": descriptive,
        "sensitivity_rates": sensitivity,
        "robustness": pd.DataFrame(robustness_rows),
        "loo": pd.DataFrame(loo_rows),
        "spearman": pd.DataFrame(spearman_rows),
    }


def write_filtering_rate_outputs(
    result: dict[str, pd.DataFrame],
    cells: pd.DataFrame,
    donors: pd.DataFrame,
    phase1: dict,
    output_dir: Path | None = None,
) -> str:
    """Write the Step 1 tables, the strip plot, and the handoff.

    Returns:
        The handoff text.
    """
    output_dir = output_directory() if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result["primary"].to_csv(output_dir / PRIMARY_FILENAME, index=False)
    result["secondary"].to_csv(output_dir / SECONDARY_FILENAME, index=False)
    result["descriptive"].to_csv(output_dir / DESCRIPTIVE_FILENAME, index=False)
    result["sensitivity_rates"].to_csv(output_dir / SENSITIVITY_FILENAME)
    result["robustness"].to_csv(output_dir / ROBUSTNESS_FILENAME, index=False)
    result["loo"].to_csv(output_dir / LOO_FILENAME, index=False)
    result["spearman"].to_csv(output_dir / SPEARMAN_FILENAME, index=False)
    plot_removal_fractions(donor_measures(cells), donors, output_dir / FIGURE_FILENAME)
    handoff = filtering_rate_handoff(result, donors, phase1, output_dir)
    (output_dir / HANDOFF_FILENAME).write_text(handoff, encoding="utf-8")
    return handoff


def plot_removal_fractions(measures: pd.DataFrame, donors: pd.DataFrame, path: Path) -> None:
    """Strip plot of per-donor fractions. Color is menopause and marker is source."""
    meta = _specimen_table(donors).set_index("donor")
    panels = [("removed_fraction", "Removed, any reason")]
    panels.extend((f"fraction_{reason}", REASON_LABELS[reason]) for reason in DISPLAY_REASONS)
    fig, axes = get_figure_axes(len(panels), num_cols=3, subplot_size=(4, 4), share_y=True)
    for ax, (column, title) in zip(axes, panels, strict=True):
        _draw_strip(ax, measures, meta, column)
        ax.set_title(title)
        ax.set_ylabel("Fraction of cells")
        ax.set_xticks([0, 1], ["BRCA1", "WT"])
        ax.set_xlim(-0.55, 1.55)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="black",
            markersize=8,
            label=status,
        )
        for status, color in MENOPAUSE_COLORS.items()
    ]
    handles.extend(
        Line2D(
            [0],
            [0],
            marker=marker,
            color="none",
            markerfacecolor="0.3",
            markeredgecolor="black",
            markersize=8,
            label=source,
        )
        for source, marker in SOURCE_MARKERS.items()
    )
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        title="Menopause and source",
    )
    fig.suptitle("Removed fraction by donor", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def filtering_rate_handoff(
    result: dict[str, pd.DataFrame], donors: pd.DataFrame, phase1: dict, output_dir: Path
) -> str:
    """Markdown summary to paste back after the notebook runs."""
    primary = result["primary"].iloc[0]
    alpha = phase1["alpha"]
    attainable = attainable_alpha(int(primary["n_permutations"]), alpha["primary"])
    lines = [
        "# Notebook 01 handoff",
        "",
        "Paste this file back into the coding session. The donor is the unit of replication.",
        "A non-significant result is inconclusive at this sample size. The minimum detectable",
        "effect is the shift with 80% power, not a bound on the true effect.",
        "",
        f"Output directory: `{output_dir}`",
        "",
        "## Primary test",
        "",
        "Overall removed fraction, BRCA1 minus WT, two-sided. This test is not part of the",
        "reason family.",
        "",
        f"- Alpha: {alpha['primary']}",
        f"- Estimate: {_points(primary['estimate_points'])} percentage points of cells",
        f"- {primary['level']:.0%} permutation interval: "
        f"[{_bound(primary['ci_low_points'])}, {_bound(primary['ci_high_points'])}] "
        "percentage points of cells",
        f"- Bootstrap interval, secondary check: "
        f"[{_bound(primary['bootstrap_low_points'])}, {_bound(primary['bootstrap_high_points'])}]",
        (
            f"- p = {int(primary['n_extreme'])}/{int(primary['n_permutations'])}"
            f" = {primary['p_value']:.6f}"
        ),
        f"- Attainable alpha at {alpha['primary']}: {attainable:.6f}",
        f"- Significant: {bool(primary['significant'])}. {_conclusion(primary)}",
        f"- MDE: {_points(primary['mde_points'])} expit percentage points "
        f"(log-odds shift {primary['mde_delta_log_odds']:.3f}, "
        f"attained: {bool(primary['mde_attained'])})",
        "- Normal approximation cross-check: "
        f"{_points(primary['mde_normal_points'])} expit percentage points",
        "",
        "The estimate and its interval are percentage points of the fraction difference",
        "(100 times BRCA1's mean fraction minus WT's). The MDE is expit percentage points",
        "around the pooled log-odds. Those two scales are not interchangeable.",
        "",
    ]
    if "mde_resample_median" in primary.index:
        lines.extend(
            [
                "Donor resampling within genotype, primary MDE only: "
                f"min {_points(primary['mde_resample_min'])}, "
                f"median {_points(primary['mde_resample_median'])}, "
                f"max {_points(primary['mde_resample_max'])} expit percentage points.",
                "",
            ]
        )
    lines.extend(
        [
            "## Secondary tests",
            "",
            "Four non-exclusive reasons, two-sided. Benjamini-Hochberg adjustment is across",
            (
                "these four p-values at q = "
                f"{alpha['secondary']}. Each MDE uses alpha = {alpha['secondary']},"
            ),
            "which is looser than that adjusted threshold, so the MDE understates the effect",
            "the adjusted test can reliably detect.",
            "",
            _reason_table(result["secondary"]),
            "",
            "## Robustness",
            "",
            "The overall result is robust only when the uniform re-run and the depth-equalized",
            "re-run both keep the conclusion and stay within half the primary MDE. The MAD",
            "re-run is reported and is not part of that verdict. Spearman correlations are a",
            "diagnostic of donor rank, not part of the verdict.",
            "",
            _robustness_text(result["robustness"]),
            "",
            "## Leave-one-donor-out",
            "",
            _n0093_text(result["loo"], float(primary["estimate"]), phase1["confidence"]["loo"]),
            "",
            "## Spearman diagnostic",
            "",
            _spearman_text(result["spearman"]),
            "",
            "## Pooled removal rates",
            "",
            "Sum of removed cells divided by sum of cells. This is descriptive. It sits next",
            "to the published per-sample rule so a stricter or looser re-run is visible.",
            "It is not a genotype test.",
            "",
            _pooled_text(result["sensitivity_rates"]),
            "",
            "## Menopause",
            "",
            _menopause_text(donors),
            "",
        ]
    )
    return "\n".join(lines)


def load_phase1(output_dir: Path) -> dict:
    """Phase 1 block written by notebook 00."""
    params = json.loads((Path(output_dir) / PARAMS_FILENAME).read_text(encoding="utf-8"))
    if "phase_1" not in params:
        raise ValueError(f"{PARAMS_FILENAME} has no phase_1 section.")
    return params["phase_1"]


def _tested_row(
    name, removed, totals, labels, *, alpha, alternative, level, ci_step, n_bootstrap, rng, mde_step
):
    values = np.asarray(removed, dtype=float) / np.asarray(totals, dtype=float)
    test = exact_permutation_test(
        values, labels, positive=POSITIVE, alternative=alternative, alpha=alpha
    )
    interval = permutation_ci(
        values,
        labels,
        positive=POSITIVE,
        level=level,
        alternative=alternative,
        step=ci_step,
        n_bootstrap=n_bootstrap,
        rng=rng,
    )
    logger.info("MDE for %s.", name)
    mde = minimum_detectable_effect(
        removed,
        totals,
        labels,
        positive=POSITIVE,
        alpha=alpha,
        alternative=alternative,
        delta_step=mde_step,
    )
    return {
        "measure": name,
        "estimate": test.statistic,
        "estimate_points": 100 * test.statistic,
        "ci_low": interval.low,
        "ci_high": interval.high,
        "ci_low_points": 100 * interval.low,
        "ci_high_points": 100 * interval.high,
        "level": level,
        "bootstrap_low_points": None
        if interval.bootstrap_low is None
        else 100 * interval.bootstrap_low,
        "bootstrap_high_points": None
        if interval.bootstrap_high is None
        else 100 * interval.bootstrap_high,
        "p_value": test.p_value,
        "n_extreme": test.n_extreme,
        "n_permutations": test.n_permutations,
        "alpha": alpha,
        "alternative": alternative,
        "significant": test.significant,
        "mde_points": mde.points,
        "mde_normal_points": mde.normal_points,
        "mde_delta_log_odds": mde.delta_log_odds,
        "mde_attained": mde.attained,
        "mde_power": mde.power,
        "mde_power_at_zero": mde.power_at_zero,
        "mde_pooled_log_odds": mde.pooled_log_odds,
    }


def _interval_row(name, removed, totals, labels, *, level, ci_step):
    values = np.asarray(removed, dtype=float) / np.asarray(totals, dtype=float)
    interval = permutation_ci(
        values,
        labels,
        positive=POSITIVE,
        level=level,
        alternative="two-sided",
        step=ci_step,
        n_bootstrap=0,
    )
    return {
        "measure": name,
        "estimate": interval.estimate,
        "estimate_points": 100 * interval.estimate,
        "ci_low": interval.low,
        "ci_high": interval.high,
        "ci_low_points": 100 * interval.low,
        "ci_high_points": 100 * interval.high,
        "level": level,
    }


def _specimen_table(donors: pd.DataFrame) -> pd.DataFrame:
    frame = donors.copy()
    if "specimen_id" not in frame.columns:
        frame = frame.reset_index(names="specimen_id")
    if "donor" not in frame.columns:
        raise ValueError("Donor table is missing a donor column.")
    return frame


def _check_alignment(measures, specimen, adatas) -> None:
    if specimen["specimen_id"].duplicated().any() or specimen["donor"].duplicated().any():
        raise ValueError("A donor or specimen id is repeated.")
    donor_ids = set(measures.index)
    if set(specimen["donor"]) != donor_ids:
        raise ValueError("Donor ids in the cell table and the donor table differ.")
    missing = set(specimen["specimen_id"]) - set(adatas)
    if missing:
        raise ValueError(f"No AnnData for specimens {sorted(missing)}.")
    n_cells = specimen.set_index("donor")
    for donor, row in measures.iterrows():
        if "n_cells" in n_cells.columns and int(n_cells.loc[donor, "n_cells"]) != int(
            row["n_cells"]
        ):
            raise ValueError(f"Cell count for {donor} disagrees between the two tables.")
        specimen_id = n_cells.loc[donor, "specimen_id"]
        if adatas[specimen_id].n_obs != int(row["n_cells"]):
            raise ValueError(f"{specimen_id} has a different cell count than the cell table.")


def _draw_strip(ax, measures, meta, column) -> None:
    for x, genotype in ((0, "BRCA1"), (1, "WT")):
        donors = measures.index[measures["genotype"] == genotype]
        order = measures.loc[donors, column].sort_values().index
        offsets = np.linspace(-0.15, 0.15, len(order)) if len(order) else []
        for donor, offset in zip(order, offsets, strict=True):
            status = str(meta.loc[donor, "menopause_status"])
            source = str(meta.loc[donor, "source_prefix"])
            ax.scatter(
                x + offset,
                measures.loc[donor, column],
                color=MENOPAUSE_COLORS.get(status, "0.5"),
                marker=SOURCE_MARKERS.get(source, "x"),
                edgecolor="black",
                linewidth=0.4,
                s=40,
                zorder=3,
            )


def _driver_text(published: pd.Series, rerun: pd.Series) -> str:
    delta = (100 * (rerun - published)).sort_values(
        key=lambda values: values.abs(), ascending=False
    )
    return "; ".join(f"{donor} {change:+.1f} pp" for donor, change in delta.head(3).items())


def _sign(statistic: float) -> int:
    if statistic > 0:
        return 1
    if statistic < 0:
        return -1
    return 0


def _points(value) -> str:
    if pd.isna(value):
        return "not attained"
    return f"{float(value):.2f}"


def _bound(value) -> str:
    if value is None or pd.isna(value):
        return "missing"
    if not np.isfinite(value):
        return "unbounded"
    return f"{float(value):.2f}"


def _conclusion(primary: pd.Series) -> str:
    if not bool(primary["significant"]):
        return "Inconclusive. Compare the estimate with the MDE."
    direction = "higher" if primary["estimate"] > 0 else "lower"
    return f"BRCA1 donors have a {direction} removed fraction than WT donors."


def _reason_table(secondary: pd.DataFrame) -> str:
    shown = pd.DataFrame(
        {
            "reason": secondary["measure"],
            "estimate_pp": secondary["estimate_points"].map(lambda value: f"{value:.2f}"),
            "ci_low_pp": secondary["ci_low_points"].map(lambda value: f"{value:.2f}"),
            "ci_high_pp": secondary["ci_high_points"].map(lambda value: f"{value:.2f}"),
            "p": secondary.apply(
                lambda row: f"{int(row['n_extreme'])}/{int(row['n_permutations'])}", axis=1
            ),
            "q": secondary["q_value"].map(lambda value: f"{value:.4f}"),
            "mde_pp": secondary["mde_points"].map(lambda value: f"{value:.2f}"),
            "normal_pp": secondary["mde_normal_points"].map(lambda value: f"{value:.2f}"),
        }
    )
    return _markdown_table(shown)


def _robustness_text(robustness: pd.DataFrame) -> str:
    lines = []
    for row in robustness.itertuples(index=False):
        extra = ""
        if row.rule == "depth_equalized" and pd.notna(row.seed_spread_points):
            extra = f" Seed spread {row.seed_spread_points:.2f} expit percentage points."
        if row.status == "fail" and row.drivers:
            extra = f" Largest donor changes: {row.drivers}."
        verdict = "in the verdict" if row.in_verdict else "not in the verdict"
        lines.append(f"- {row.rule} ({verdict}): {row.status}.{extra}")
    return "\n".join(lines)


def _n0093_text(loo: pd.DataFrame, full_estimate: float, level: float) -> str:
    overall = loo.loc[loo["measure"] == "overall_removed_fraction"]
    dropped = overall.loc[overall["dropped"] == "N-0093"]
    if dropped.empty:
        return "N-0093 is not in this cohort, so its leave-one-out contrast was not computed."
    row = dropped.iloc[0]
    change = float(row["estimate"]) - full_estimate
    flipped = float(row["estimate"]) * full_estimate < 0
    width = float(row["ci_width_points"])
    larger = abs(100 * change) > width
    low = row["ci_low"] * 100
    high = row["ci_high"] * 100
    return "\n".join(
        [
            f"Every leave-one-out run uses a two-sided {level:.0%} interval. The comparison",
            "below drops N-0093 from the primary overall test. The noise is that interval's width.",
            "",
            f"- Full-data estimate: {100 * full_estimate:.2f} percentage points",
            f"- Without N-0093: {row['estimate_points']:.2f} percentage points",
            f"- {level:.0%} interval: [{low:.2f}, {high:.2f}], width {width:.2f}",
            f"- Sign flip: {'yes' if flipped else 'no'}",
            f"- Absolute change exceeds the interval width: {'yes' if larger else 'no'}",
        ]
    )


def _spearman_text(spearman: pd.DataFrame) -> str:
    lines = []
    for row in spearman.itertuples(index=False):
        rho = "undefined" if pd.isna(row.rho) else f"{row.rho:.3f}"
        lines.append(f"- {row.rule}: {rho}")
    return "\n".join(lines)


def _pooled_text(rates: pd.DataFrame) -> str:
    lines = []
    for rule, removed_column in (
        ("published per-sample", "published_removed"),
        ("uniform", "uniform_removed"),
        ("MAD", "mad_removed"),
        ("depth-equalized", "depth_removed"),
    ):
        pooled = rates[removed_column].sum() / rates["n_cells"].sum()
        lines.append(f"- {rule}: {100 * pooled:.2f}% of cells removed")
    return "\n".join(lines)


def _menopause_text(donors: pd.DataFrame) -> str:
    frame = _specimen_table(donors)
    post = frame.loc[frame["menopause_status"] != "Pre", ["donor", "genotype", "menopause_status"]]
    if post.empty:
        return "Every donor in this run is premenopausal."
    lines = [
        "Menopause is not stratified. Two BRCA1 donors can be post-oophorectomy while every",
        "WT donor is premenopausal, so a genotype contrast is also a menopause contrast.",
        "",
    ]
    lines.extend(
        f"- {row.donor} ({row.genotype}): {row.menopause_status}" for row in post.itertuples()
    )
    return "\n".join(lines)


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(str(row[column]) for column in frame.columns) + " |"
        for _, row in frame.iterrows()
    ]
    return "\n".join([header, separator, *body])
