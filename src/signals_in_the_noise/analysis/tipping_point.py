"""Tipping point for a luminal-progenitor expansion hidden by QC removal.

Removed cells keep no lineage label. Each scenario fills that share in from
the donor's retained share. The test is the same one-sided donor permutation
as the retained-cell baseline.
"""

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from signals_in_the_noise.analysis.donor_permutation import (
    exact_permutation_test,
    leave_one_donor_out,
    minimum_detectable_effect,
    permutation_ci,
)
from signals_in_the_noise.analysis.lineage import LP_TYPE
from signals_in_the_noise.analysis.pal_reference import (
    largest_enrichment,
    near_threshold_enrichment,
)
from signals_in_the_noise.analysis.qc_filtering import PARAMS_FILENAME, output_directory

HAND_OFF = "03_handoff.md"
FIGURE = "03-tipping-point.png"
_K_TOLERANCE = 0.001
_CURVE_POINTS = 201
_ZONE_PHASE2 = {
    "plausible": "Phase 2 goes ahead, and its LP enrichment is compared with k*.",
    "possible_but_implausible": (
        "Phase 2 is a confirmation. The required enrichment is above the signature ceiling."
    ),
    "impossible": "Phase 2 is optional for the BRCA1 question.",
}


def removed_lp_share(lp_share, genotype, k: float, *, extreme: bool = False) -> np.ndarray:
    """Hypothetical LP share among each donor's removed cells.

    At the extreme, every BRCA1 removed cell is LP and no WT removed cell is.
    Otherwise WT removed cells keep that donor's retained share, and BRCA1
    removed cells take ``k`` times theirs, capped at 1.
    """
    share = np.asarray(lp_share, dtype=float)
    geno = np.asarray(genotype)
    if extreme:
        return np.where(geno == "BRCA1", 1.0, 0.0).astype(float)
    out = share.copy()
    brca1 = geno == "BRCA1"
    out[brca1] = np.minimum(1.0, float(k) * share[brca1])
    return out


def full_lp_shares(n_retained, n_lp, n_removed, removed_share) -> np.ndarray:
    """LP share of retained plus removed cells under a removed-cell scenario."""
    n_retained = np.asarray(n_retained, dtype=float)
    n_lp = np.asarray(n_lp, dtype=float)
    n_removed = np.asarray(n_removed, dtype=float)
    removed_share = np.asarray(removed_share, dtype=float)
    n_cells = n_retained + n_removed
    if np.any(n_cells <= 0):
        raise ValueError("A donor has no cells in the scenario.")
    return (n_lp + n_removed * removed_share) / n_cells


def k_max_value(lp_share, genotype) -> float:
    """``1 / smallest BRCA1 retained LP share``. Infinite when that share is 0."""
    brca1 = np.asarray(lp_share, dtype=float)[np.asarray(genotype) == "BRCA1"]
    if len(brca1) == 0:
        raise ValueError("No BRCA1 donor.")
    smallest = float(np.min(brca1))
    if smallest <= 0:
        return math.inf
    return 1.0 / smallest


def search_cap(lp_share, genotype) -> float:
    """Largest k that still changes a BRCA1 removed share. Zero when none can."""
    brca1 = np.asarray(lp_share, dtype=float)[np.asarray(genotype) == "BRCA1"]
    positive = brca1[brca1 > 0]
    if len(positive) == 0:
        return 0.0
    return float(np.max(1.0 / positive))


def smallest_rejecting_k(
    p_at, k_cap: float, alpha: float, tol: float = _K_TOLERANCE
) -> float | None:
    """Smallest k at which ``p_at(k)`` is at most ``alpha``, or None.

    The returned k rejects. It sits within ``tol`` of the boundary.
    """
    if k_cap < 0:
        raise ValueError(f"k_cap must be non-negative, got {k_cap}.")
    if p_at(0.0) <= alpha:
        return 0.0
    if k_cap == 0 or p_at(k_cap) > alpha:
        return None
    low = 0.0
    high = float(k_cap)
    while high - low > tol:
        mid = (low + high) / 2.0
        if p_at(mid) <= alpha:
            high = mid
        else:
            low = mid
    return high


def tipping_zone(k_star: float | None, signature_c: float, k_max: float) -> str:
    """Zone from the frozen rule. A missing k* is above every attainable k."""
    if k_star is None:
        return "impossible"
    if k_star <= signature_c:
        return "plausible"
    if math.isinf(k_max) or k_star <= k_max:
        return "possible_but_implausible"
    return "impossible"


def doublet_only(cells: pd.DataFrame) -> pd.Series:
    """Removed cells whose only failures are the high-gene or high-count flags.

    A cell that also fails the mito or low-gene flag stays in the sensitivity
    count. Those two high flags are Pal's doublet filter.
    """
    low_quality = cells["is_low_num_genes"].astype(int).eq(1) | cells["is_high_mito"].astype(
        int
    ).eq(1)
    doublet = cells["is_high_num_genes"].astype(int).eq(1) | cells["is_high_total_count"].astype(
        int
    ).eq(1)
    return doublet & ~low_quality & cells["is_noise"].astype(int).eq(1)


def run_tipping_point(
    cells: pd.DataFrame,
    donors: pd.DataFrame,
    phase1: dict,
    labels: pd.DataFrame,
    shares: pd.DataFrame,
    *,
    rng: np.random.Generator | None = None,
) -> dict:
    """Baseline, k*, zone, extreme bound, and the doublet-only sensitivity."""
    spec = phase1["tests"]["lp_baseline_and_tipping_point"]
    positive, alternative, alpha = _test_spec(spec)
    level = float(phase1["confidence"]["primary"])
    fraction = float(phase1["near_threshold_band"]["fraction"])
    power = float(phase1["mde"]["power"])
    loo_level = float(phase1["loo"]["level"])
    frame = _analysis_frame(cells, donors, labels, shares, drop_doublet_only=False)
    enrichment, ceiling = signature_ceiling(cells, labels, donors, fraction)
    baseline = _baseline(frame, positive, alternative, alpha, level, power, loo_level, rng)
    primary = _k_search(frame, positive, alternative, alpha)
    extreme = _scenario_test(
        frame, 0.0, extreme=True, positive=positive, alternative=alternative, alpha=alpha
    )
    sensitive = _analysis_frame(cells, donors, labels, shares, drop_doublet_only=True)
    sensitivity = _k_search(sensitive, positive, alternative, alpha)
    zone = tipping_zone(primary["k_star"], ceiling["c"], primary["k_max"])
    sensitivity_zone = tipping_zone(sensitivity["k_star"], ceiling["c"], sensitivity["k_max"])
    return {
        "donor_table": donors,
        "donors": frame,
        "baseline": baseline,
        "curve": primary["curve"],
        "k_star": primary["k_star"],
        "k_max": primary["k_max"],
        "zone": zone,
        "signature_c": ceiling["c"],
        "winners": ceiling["winners"],
        "enrichment": enrichment,
        "cluster_c": phase1["near_threshold_band"].get("c"),
        "extreme_p": extreme.p_value,
        "extreme_statistic": extreme.statistic,
        "n_doublet_only": int(doublet_only(cells).sum()),
        "sensitivity_k_star": sensitivity["k_star"],
        "sensitivity_zone": sensitivity_zone,
        "sensitivity_curve": sensitivity["curve"],
        "alpha": alpha,
        "positive": positive,
        "alternative": alternative,
    }


def signature_ceiling(
    cells: pd.DataFrame,
    labels: pd.DataFrame,
    donors: pd.DataFrame,
    fraction: float,
) -> tuple[pd.DataFrame, dict]:
    """Near-threshold enrichment of the Lim 2009 types. Genotype is not used."""
    if "predicted_type" not in labels.columns:
        raise ValueError("Labels are missing predicted_type.")
    typed = labels.rename(columns={"predicted_type": "cell_type"})
    enrichment = near_threshold_enrichment(
        cells, typed, _with_donor_column(donors), fraction=fraction
    )
    return enrichment, largest_enrichment(enrichment)


def write_tipping_point_outputs(result: dict, output_dir: Path | None = None) -> str:
    """Write the tables, the curve, and the handoff. The cluster enrichment stays put."""
    output_dir = output_directory() if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result["donors"].to_csv(output_dir / "03_donor_shares.csv", index=False)
    result["curve"].to_csv(output_dir / "03_curve.csv", index=False)
    result["baseline"]["loo"].to_csv(output_dir / "03_loo.csv", index=False)
    result["enrichment"].to_csv(output_dir / "03_enrichment.csv", index=False)
    _baseline_row(result).to_csv(output_dir / "03_baseline.csv", index=False)
    plot_tipping_curve(
        result["curve"],
        result["k_star"],
        result["signature_c"],
        result["alpha"],
        output_dir / FIGURE,
    )
    _record_decision(output_dir / PARAMS_FILENAME, result)
    handoff = tipping_point_handoff(result, output_dir)
    (output_dir / HAND_OFF).write_text(handoff, encoding="utf-8")
    return handoff


def plot_tipping_curve(
    curve: pd.DataFrame,
    k_star: float | None,
    signature_c: float,
    alpha: float,
    path: Path,
) -> None:
    """p-value against k, with k* and the signature ceiling marked."""
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.plot(curve["k"], curve["p_value"], color="#0072B2", label="one-sided p-value")
    ax.axhline(alpha, color="#666666", linestyle="--", linewidth=1, label=f"α = {alpha:g}")
    if k_star is not None:
        ax.axvline(k_star, color="#D55E00", linewidth=1, label=f"k* = {k_star:.3f}")
    ax.axvline(
        signature_c,
        color="#009E73",
        linestyle=":",
        linewidth=1,
        label=f"signature c = {signature_c:.3f}",
    )
    ax.set_xlabel("k (BRCA1 removed LP share / retained LP share)")
    ax.set_ylabel("one-sided p-value")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def tipping_point_handoff(result: dict, output_dir: Path) -> str:
    """Markdown summary to paste back after the notebook runs."""
    baseline = result["baseline"]
    test = baseline["test"]
    interval = baseline["interval"]
    mde = baseline["mde"]
    n_brca1 = int((result["donors"]["genotype"] == "BRCA1").sum())
    n_wt = int((result["donors"]["genotype"] == "WT").sum())
    lines = [
        "# Notebook 03 handoff",
        "",
        "Paste this file back into the coding session. Removed cells have no lineage",
        "label. Each k is a hypothetical LP share among BRCA1 removed cells.",
        "",
        f"Output directory: `{output_dir}`",
        "",
        "## Retained LP share",
        "",
        "Labels are Lim 2009 signatures scored with scanpy on retained cells only,",
        "highly variable gene filtering off, random_state 0. The test is one-sided",
        f"({result['positive']} greater than WT) at α = {result['alpha']:g}.",
        "The donor is the unit.",
        "",
        (
            f"BRCA1 minus WT = {100 * test.statistic:.2f} percentage points. "
            f"99% lower bound {100 * interval.low:.2f}. "
            f"p = {test.n_extreme}/{test.n_permutations} = {test.p_value:.6f}."
        ),
        (
            f"MDE at 80% power: {mde.points:.2f} percentage points "
            f"(normal cross-check {mde.normal_points:.2f})."
        ),
        _baseline_verdict(test, n_brca1, n_wt),
        "",
        "At k = 1 the all-cell LP share equals the retained LP share, because removed",
        "cells are given the same share as the cells that were kept.",
        "",
        _markdown_table(_share_table(result["donors"])),
        "",
        "## Tipping point",
        "",
        "WT removed cells take that donor's retained LP share. BRCA1 removed cells",
        "take k times theirs, capped at 1. k* is the smallest k at which the same",
        "one-sided test on the all-cell LP share reaches p ≤ α. The search stops",
        "once every BRCA1 removed share that can rise has reached 1.",
        "",
        _k_sentence(result),
        _zone_sentence(result["zone"], result["signature_c"], result["k_max"]),
        _ZONE_PHASE2[result["zone"]],
        "",
        (
            f"Signature ceiling c = {result['signature_c']:.4f} "
            f"({_winner_text(result['winners'])}). "
            f"This is the largest near-threshold enrichment of a Lim 2009 type."
        ),
        _cluster_sentence(result["cluster_c"]),
        "",
        (
            "Extreme bound, every BRCA1 removed cell LP and no WT removed cell LP: "
            f"p = {result['extreme_p']:.6f}."
        ),
        "",
        "## Doublet-only sensitivity",
        "",
        "Removed cells that fail only the high-gene or high-count flag are left out",
        "of the removed count. A cell that also fails the mito or low-gene flag stays in.",
        f"Doublet-only removed cells: {result['n_doublet_only']}.",
        _sensitivity_sentence(result),
        "",
        "## Leave-one-donor-out",
        "",
        "The estimate is recomputed after each donor is dropped. This is a stability",
        "check. It is not a new test at α = 0.01.",
        "",
        _markdown_table(_loo_table(baseline["loo"], test.statistic)),
        "",
        "## Confounding",
        "",
        "The contrast is genotype-associated. A non-significant result is inconclusive",
        f"at this sample size (n = {n_brca1} versus {n_wt}).",
    ]
    menopause = _menopause_sentence(result)
    if menopause:
        lines.append(menopause)
    return "\n".join(lines) + "\n"


def _test_spec(spec: dict) -> tuple[str, str, float]:
    positive = str(spec["contrast"]).split(">")[0].strip()
    return positive, str(spec["alternative"]), float(spec["alpha"])


def _analysis_frame(cells, donors, labels, shares, *, drop_doublet_only: bool) -> pd.DataFrame:
    frame = _shares_with_donor_column(shares)
    frame["donor"] = frame["donor"].astype(str)
    counted = labels.assign(donor=labels["donor"].astype(str))
    n_lp = counted.groupby("donor")["predicted_type"].apply(lambda s: int((s == LP_TYPE).sum()))
    n_labeled = counted.groupby("donor").size()
    if not frame["donor"].isin(n_lp.index).all():
        raise ValueError("A donor in the share table has no labels.")
    mismatch = frame["n_lp"].to_numpy() != frame["donor"].map(n_lp).to_numpy()
    if mismatch.any():
        raise ValueError("n_lp does not match the lp labels.")
    retained = cells.loc[cells["is_noise"].astype(int) == 0].copy()
    retained["donor"] = retained["donor"].astype(str)
    n_retained = retained.groupby("donor").size()
    if not np.array_equal(
        frame["n_retained"].to_numpy(), frame["donor"].map(n_retained).to_numpy()
    ):
        raise ValueError("n_retained does not match the cell table.")
    if not np.array_equal(frame["n_retained"].to_numpy(), frame["donor"].map(n_labeled).to_numpy()):
        raise ValueError("Label counts do not match retained cells.")
    donor_rows = _with_donor_column(donors)
    if donor_rows["donor"].astype(str).duplicated().any():
        raise ValueError("A donor id is repeated.")
    geno = donor_rows.set_index(donor_rows["donor"].astype(str))["genotype"]
    frame["genotype"] = frame["donor"].map(geno)
    if frame["genotype"].isna().any():
        missing = sorted(frame.loc[frame["genotype"].isna(), "donor"])
        raise ValueError(f"No genotype for donors {missing}.")
    removed = _removed_counts(cells, drop_doublet_only)
    frame["n_removed"] = frame["donor"].map(removed).fillna(0).astype(int)
    frame["lp_share"] = frame["n_lp"] / frame["n_retained"]
    return frame.sort_values("donor").reset_index(drop=True)


def _removed_counts(cells: pd.DataFrame, drop_doublet_only: bool) -> pd.Series:
    noise = cells.loc[cells["is_noise"].astype(int) == 1].copy()
    noise["donor"] = noise["donor"].astype(str)
    if drop_doublet_only:
        noise = noise.loc[~doublet_only(noise).to_numpy()]
    if noise.empty:
        return pd.Series(dtype=int)
    return noise.groupby("donor").size()


def _baseline(frame, positive, alternative, alpha, level, power, loo_level, rng):
    values = frame["lp_share"].to_numpy()
    labels = frame["genotype"].to_numpy()
    test = exact_permutation_test(
        values, labels, positive=positive, alternative=alternative, alpha=alpha
    )
    interval = permutation_ci(
        values,
        labels,
        positive=positive,
        level=level,
        alternative=alternative,
        n_bootstrap=2000,
        rng=rng,
    )
    mde = minimum_detectable_effect(
        frame["n_lp"].to_numpy(),
        frame["n_retained"].to_numpy(),
        labels,
        positive=positive,
        alpha=alpha,
        alternative=alternative,
        power=power,
    )
    loo = pd.DataFrame(
        leave_one_donor_out(
            values,
            labels,
            frame["donor"].to_numpy(),
            positive=positive,
            level=loo_level,
        )
    )
    return {"test": test, "interval": interval, "mde": mde, "loo": loo}


def _k_search(frame, positive, alternative, alpha) -> dict:
    k_max = k_max_value(frame["lp_share"], frame["genotype"])
    k_cap = search_cap(frame["lp_share"], frame["genotype"])

    def p_at(k: float) -> float:
        return _scenario_test(
            frame, k, extreme=False, positive=positive, alternative=alternative, alpha=alpha
        ).p_value

    k_star = smallest_rejecting_k(p_at, k_cap, alpha)
    grid = _k_grid(k_cap, k_star)
    rows = []
    for k in grid:
        result = _scenario_test(
            frame, float(k), extreme=False, positive=positive, alternative=alternative, alpha=alpha
        )
        rows.append({"k": float(k), "statistic": result.statistic, "p_value": result.p_value})
    return {"k_star": k_star, "k_max": k_max, "curve": pd.DataFrame(rows)}


def _scenario_test(frame, k, *, extreme, positive, alternative, alpha):
    removed = removed_lp_share(frame["lp_share"], frame["genotype"], k, extreme=extreme)
    shares = full_lp_shares(frame["n_retained"], frame["n_lp"], frame["n_removed"], removed)
    return exact_permutation_test(
        shares,
        frame["genotype"].to_numpy(),
        positive=positive,
        alternative=alternative,
        alpha=alpha,
    )


def _k_grid(k_cap: float, k_star: float | None) -> np.ndarray:
    if k_cap == 0:
        grid = np.array([0.0])
    else:
        grid = np.linspace(0.0, k_cap, _CURVE_POINTS)
    extra = [1.0] if k_cap >= 1 else []
    if k_star is not None:
        extra.append(k_star)
    values = np.unique(np.concatenate([grid, extra]))
    return values[values <= k_cap + 1e-12]


def _shares_with_donor_column(shares: pd.DataFrame) -> pd.DataFrame:
    """Donor id as a column. The lineage table uses that id as the index name too."""
    frame = shares.copy()
    if "donor" not in frame.columns:
        return frame.rename_axis("donor").reset_index()
    if frame.index.name == "donor":
        return frame.rename_axis(None)
    return frame


def _with_donor_column(donors: pd.DataFrame) -> pd.DataFrame:
    frame = donors.copy()
    if "donor" not in frame.columns:
        frame = frame.reset_index(names="donor")
    return frame


def _baseline_row(result: dict) -> pd.DataFrame:
    test = result["baseline"]["test"]
    interval = result["baseline"]["interval"]
    mde = result["baseline"]["mde"]
    return pd.DataFrame(
        [
            {
                "statistic": test.statistic,
                "p_value": test.p_value,
                "n_extreme": test.n_extreme,
                "n_permutations": test.n_permutations,
                "significant": test.significant,
                "ci_low": interval.low,
                "mde_points": mde.points,
                "mde_normal_points": mde.normal_points,
                "k_star": result["k_star"],
                "k_max": result["k_max"],
                "zone": result["zone"],
                "signature_c": result["signature_c"],
                "cluster_c": result["cluster_c"],
                "extreme_p": result["extreme_p"],
                "sensitivity_k_star": result["sensitivity_k_star"],
                "sensitivity_zone": result["sensitivity_zone"],
            }
        ]
    )


def _record_decision(path: Path, result: dict) -> None:
    existing = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    phase1 = existing.setdefault("phase_1", {})
    block = phase1.setdefault("tipping_point", {})
    block["signature_c"] = result["signature_c"]
    block["signature_c_winners"] = result["winners"].to_dict(orient="records")
    block["k_star"] = result["k_star"]
    block["k_max"] = None if math.isinf(result["k_max"]) else result["k_max"]
    block["k_max_unbounded"] = bool(math.isinf(result["k_max"]))
    block["zone"] = result["zone"]
    block["extreme_p"] = result["extreme_p"]
    block["sensitivity_k_star"] = result["sensitivity_k_star"]
    block["sensitivity_zone"] = result["sensitivity_zone"]
    path.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")


def _share_table(frame: pd.DataFrame) -> pd.DataFrame:
    shown = frame[["donor", "genotype", "n_retained", "n_lp", "n_removed", "lp_share"]].copy()
    shown["lp_share"] = shown["lp_share"].map(lambda value: f"{value:.4f}")
    return shown


def _loo_table(loo: pd.DataFrame, full_estimate: float) -> pd.DataFrame:
    shown = loo[["dropped", "dropped_group", "estimate"]].copy()
    shown["sign_flip"] = np.sign(shown["estimate"]) != np.sign(full_estimate)
    shown.loc[np.isclose(full_estimate, 0.0) & np.isclose(shown["estimate"], 0.0), "sign_flip"] = (
        False
    )
    shown["estimate"] = shown["estimate"].map(lambda value: f"{100 * value:.2f}")
    return shown


def _winner_text(winners: pd.DataFrame) -> str:
    return ", ".join(
        f"{row.cell_type} at {row.threshold}" for row in winners.itertuples(index=False)
    )


def _k_sentence(result: dict) -> str:
    k_max = result["k_max"]
    cap = (
        "unbounded, because a BRCA1 retained LP share is 0" if math.isinf(k_max) else f"{k_max:.4f}"
    )
    if result["k_star"] is None:
        return f"k* is above the search cap. k_max = {cap}. The test does not reach p ≤ α."
    return f"k* = {result['k_star']:.3f}. k_max = {cap}."


def _zone_sentence(zone: str, signature_c: float, k_max: float) -> str:
    return f"Zone: {zone} (signature c = {signature_c:.4f}, k_max = {_format_k(k_max)})."


def _cluster_sentence(cluster_c) -> str:
    if cluster_c is None:
        return "Notebook 02 has not stored a Fig 4B cluster enrichment."
    return (
        f"Notebook 02 stored the Fig 4B cluster enrichment c = {float(cluster_c):.4f} "
        "on its own. The zone uses the signature ceiling."
    )


def _sensitivity_sentence(result: dict) -> str:
    k_star = result["sensitivity_k_star"]
    shown = "above the search cap" if k_star is None else f"{k_star:.3f}"
    change = (
        "The zone is unchanged."
        if result["sensitivity_zone"] == result["zone"]
        else (f"The zone changes from {result['zone']} to {result['sensitivity_zone']}.")
    )
    return f"Sensitivity k* = {shown}. Zone: {result['sensitivity_zone']}. {change}"


def _baseline_verdict(test, n_brca1: int, n_wt: int) -> str:
    if test.significant:
        return f"The retained-cell test reaches α = {test.alpha:g}."
    return (
        f"The retained-cell difference stays inside α = {test.alpha:g}. "
        f"Inconclusive given n = {n_brca1} versus {n_wt}."
    )


def _menopause_sentence(result: dict) -> str:
    donors = result.get("donor_table")
    if donors is None or "menopause_status" not in donors.columns:
        return ""
    post = donors.loc[
        donors["menopause_status"] != "Pre", ["donor", "genotype", "menopause_status"]
    ]
    if post.empty:
        return "Every donor in the table is premenopausal."
    listed = ", ".join(
        f"{row.donor} ({row.genotype}, {row.menopause_status})"
        for row in post.itertuples(index=False)
    )
    return (
        f"Menopause is not stratified. Post-menopausal donors: {listed}. "
        "The genotype contrast is also a menopause contrast."
    )


def _format_k(value: float) -> str:
    if math.isinf(value):
        return "unbounded"
    return f"{value:.4f}"


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(str(row[column]) for column in frame.columns) + " |"
        for _, row in frame.iterrows()
    ]
    return "\n".join([header, separator, *body])
