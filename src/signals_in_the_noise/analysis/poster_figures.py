"""Poster figures in the ISCB Michigan palette.

The hex values match ``notebooks/GSE161529/iscb-scs/12-a-iscb-scs-volcano.ipynb``.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

MI_PALETTE = (
    "#000000",
    "#0072B2",
    "#009E73",
    "#56B4E9",
    "#CC79A7",
    "#E69F00",
    "#D55E00",
)
MI_FONT_COLOR = "#00274C"
MI_BLUE = MI_PALETTE[1]
MI_GREEN = MI_PALETTE[2]
MI_SKY = MI_PALETTE[3]
MI_ORANGE = MI_PALETTE[5]
MI_VERMILLION = MI_PALETTE[6]

_ER_PLUS = "ER+ tumour"
_NORMAL = "Normal"
_POST = "Post (oophorectomy)"


def plot_depth_collapse(summary: pd.DataFrame, depths: pd.DataFrame, path: Path) -> None:
    """Specimen AUC before and after thinning, beside library depth by condition.

    ``summary`` has one row per version and arm, with ``specimen_auc``.
    ``depths`` has ``condition``, ``depth``, and ``included`` for each specimen.
    """
    _save(depth_collapse_figure(summary, depths), path)


def depth_collapse_figure(summary: pd.DataFrame, depths: pd.DataFrame):
    """The depth figure, before it is written to disk.

    Re-plot: the PBS slope is the claim. Draw it last, heavier than the
    kept-cell line, and print the specimen AUC at each end. Those two numbers
    are what the title is asserting. The cohort is ER+ tumour versus Normal.
    """
    scores = _noise_and_qc(summary)
    included = depths.loc[_as_bool(depths["included"])].copy()
    before, after = _pbs_auc(scores)
    fig = _figure(figsize=(10.2, 4.6))
    axes = fig.subplots(1, 2)
    _auc_slopes(axes[0], scores)
    _depth_strip(axes[1], included)
    fig.suptitle(
        "Equalizing depth removes the ER+ signal in discarded cells\n"
        f"PBS specimen AUC {before:.2f} before matching, {after:.2f} after",
        color=MI_FONT_COLOR,
        fontsize=15,
        x=0.02,
        ha="left",
    )
    return fig


def plot_lp_expansion(
    curve: pd.DataFrame,
    shares: pd.DataFrame,
    path: Path,
    *,
    alpha: float,
    signature_c: float,
    extreme_p: float,
) -> None:
    """Tipping-point curve and the retained LP share of each donor.

    ``shares`` includes ``donor``, ``genotype``, ``lp_share``, and, when known,
    ``menopause_status``.
    """
    _save(lp_expansion_figure(curve, shares, alpha, signature_c, extreme_p), path)


def lp_expansion_figure(curve, shares, alpha: float, signature_c: float, extreme_p: float):
    """The progenitor figure, before it is written to disk."""
    fig = _figure(figsize=(10.6, 4.4))
    axes = fig.subplots(1, 2, width_ratios=[1.35, 1])
    _tipping_curve(axes[0], curve, alpha=alpha, signature_c=signature_c, extreme_p=extreme_p)
    _lp_strip(axes[1], shares)
    fig.suptitle(
        "Removed cells cannot create a BRCA1 LP expansion\n"
        "12 donors, including the two post-oophorectomy BRCA1 samples",
        color=MI_FONT_COLOR,
        fontsize=15,
        x=0.02,
        ha="left",
    )
    return fig


def _noise_and_qc(summary: pd.DataFrame) -> pd.DataFrame:
    frame = summary.copy()
    if "version" not in frame.columns or "arm" not in frame.columns:
        frame = frame.reset_index()
    needed = {"version", "arm", "specimen_auc", "permutation_p"}
    missing = needed - set(frame.columns)
    if missing:
        raise ValueError(f"Depth summary is missing {sorted(missing)}.")
    keep = frame["arm"].isin(["noise_pbs", "retained_qc"])
    chosen = frame.loc[keep, ["version", "arm", "specimen_auc", "permutation_p"]].copy()
    versions = set(chosen.loc[chosen["arm"] == "noise_pbs", "version"])
    if versions != {"original", "depth-matched"}:
        raise ValueError(
            f"noise_pbs must have original and depth-matched rows, found {sorted(versions)}."
        )
    return chosen


def _pbs_auc(scores: pd.DataFrame) -> tuple[float, float]:
    rows = scores.loc[scores["arm"] == "noise_pbs"].set_index("version")
    return (
        float(rows.loc["original", "specimen_auc"]),
        float(rows.loc["depth-matched", "specimen_auc"]),
    )


def _auc_slopes(ax, scores: pd.DataFrame) -> None:
    order = ["original", "depth-matched"]
    # Kept cells first, then PBS on top. PBS carries the AUC labels.
    series = (
        ("retained_qc", MI_SKY, 1.6, 6, 2),
        ("noise_pbs", MI_VERMILLION, 3.2, 9, 4),
    )
    labels = {"noise_pbs": "Discarded cells, PBS", "retained_qc": "Kept cells, raw QC"}
    for arm, color, width, marker, zorder in series:
        rows = scores.loc[scores["arm"] == arm].set_index("version")
        if not set(order).issubset(rows.index):
            continue
        values = [float(rows.loc[version, "specimen_auc"]) for version in order]
        ax.plot(
            [0, 1],
            values,
            color=color,
            linewidth=width,
            marker="o",
            markersize=marker,
            zorder=zorder,
        )
        if arm != "noise_pbs":
            continue
        for x, version, value, align, shift in (
            (0, "original", values[0], "left", 8),
            (1, "depth-matched", values[1], "right", -8),
        ):
            ax.annotate(
                f"AUC {value:.2f}",
                (x, value),
                textcoords="offset points",
                xytext=(shift, 8),
                ha=align,
                color=color,
                fontsize=12,
            )
            ax.annotate(
                f"p = {float(rows.loc[version, 'permutation_p']):.2f}",
                (x, value),
                textcoords="offset points",
                xytext=(shift, -14),
                ha=align,
                color=color,
                fontsize=10,
            )
    ax.axhline(0.5, color=MI_FONT_COLOR, linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)
    ax.set_xticks([0, 1], ["Original depth", "Depth-matched"])
    ax.set_ylabel("Specimen AUC, ER+ vs Normal")
    ax.set_ylim(0.0, 1.05)
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=MI_VERMILLION,
                marker="o",
                linewidth=3.2,
                label=labels["noise_pbs"],
            ),
            Line2D([0], [0], color=MI_SKY, marker="o", linewidth=1.6, label=labels["retained_qc"]),
        ],
        frameon=False,
        loc="lower left",
    )
    _style(ax)


def _depth_strip(ax, included: pd.DataFrame) -> None:
    colors = {_ER_PLUS: MI_GREEN, _NORMAL: MI_BLUE}
    rng = np.random.default_rng(0)
    for index, (condition, color) in enumerate(colors.items()):
        values = included.loc[included["condition"] == condition, "depth"].to_numpy(dtype=float)
        if len(values) == 0:
            continue
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s=36,
            color=color,
            edgecolors=MI_FONT_COLOR,
            linewidths=0.4,
            zorder=3,
        )
    ax.set_xticks([0, 1], ["ER+", "Normal"])
    ax.set_ylabel("Median UMIs per retained cell")
    _style(ax)


def _tipping_curve(ax, curve, *, alpha, signature_c, extreme_p) -> None:
    ordered = curve.sort_values("k")
    ax.plot(ordered["k"], ordered["p_value"], color=MI_BLUE, linewidth=2.2)
    ax.axhline(alpha, color=MI_VERMILLION, linestyle="--", linewidth=1.0)
    ax.axvline(signature_c, color=MI_ORANGE, linestyle=":", linewidth=1.2)
    cap = ordered.iloc[-1]
    cap_p = float(cap["p_value"])
    ax.scatter([cap["k"]], [cap_p], s=42, color=MI_VERMILLION, zorder=4)
    ax.annotate(
        "Cap: all removed BRCA1 cells are LP\n"
        "WT removed cells keep their retained share\n"
        f"p = {cap_p:.3f}",
        (cap["k"], cap_p),
        textcoords="offset points",
        xytext=(-8, 12),
        ha="right",
        color=MI_VERMILLION,
        fontsize=10,
    )
    y_top = 0.97
    ax.text(
        signature_c,
        y_top,
        "plausible\nceiling",
        color=MI_ORANGE,
        fontsize=10,
        ha="left",
        va="top",
    )
    ax.text(
        0.02,
        alpha + 0.03,
        f"α = {alpha:g}",
        color=MI_VERMILLION,
        fontsize=11,
        transform=ax.get_yaxis_transform(),
    )
    ax.text(
        0.98,
        0.08,
        f"Separate bound, no LP in removed WT cells: p = {extreme_p:.3f}",
        transform=ax.transAxes,
        ha="right",
        color=MI_FONT_COLOR,
        fontsize=10,
    )
    ax.set_xlim(0, float(ordered["k"].max()) * 1.04)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Removed-cell LP enrichment, k")
    ax.set_ylabel("One-sided p-value")
    _style(ax)


def _lp_strip(ax, shares: pd.DataFrame) -> None:
    colors = {"BRCA1": MI_VERMILLION, "WT": MI_BLUE}
    rng = np.random.default_rng(0)
    menopause = shares["menopause_status"] if "menopause_status" in shares.columns else None
    for index, (genotype, color) in enumerate(colors.items()):
        group = shares.loc[shares["genotype"] == genotype]
        values = group["lp_share"].to_numpy(dtype=float) * 100
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        post = np.zeros(len(group), dtype=bool)
        if menopause is not None:
            post = group["menopause_status"].to_numpy() == _POST
        pre = ~post
        ax.scatter(
            np.full(pre.sum(), index) + jitter[pre],
            values[pre],
            s=42,
            color=color,
            edgecolors=MI_FONT_COLOR,
            linewidths=0.4,
            zorder=3,
        )
        if post.any():
            ax.scatter(
                np.full(post.sum(), index) + jitter[post],
                values[post],
                s=48,
                marker="D",
                color=color,
                edgecolors=MI_FONT_COLOR,
                linewidths=0.4,
                zorder=4,
            )
    ax.set_xticks([0, 1], ["BRCA1", "WT"])
    ax.set_ylabel("Retained LP share (%)")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=MI_VERMILLION,
            markeredgecolor=MI_FONT_COLOR,
            label="BRCA1, pre",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="none",
            markerfacecolor=MI_VERMILLION,
            markeredgecolor=MI_FONT_COLOR,
            label="BRCA1, post-oophorectomy",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=MI_BLUE,
            markeredgecolor=MI_FONT_COLOR,
            label="WT, pre",
        ),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=9)
    _style(ax)


def _figure(*, figsize: tuple[float, float]) -> Figure:
    figure = Figure(figsize=figsize)
    FigureCanvasAgg(figure)
    return figure


def _as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.lower().isin(["true", "1", "yes"])


def _style(ax) -> None:
    ax.set_facecolor("none")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color(MI_FONT_COLOR)
    ax.spines["left"].set_color(MI_FONT_COLOR)
    ax.tick_params(colors=MI_FONT_COLOR, labelsize=12)
    ax.xaxis.label.set_color(MI_FONT_COLOR)
    ax.yaxis.label.set_color(MI_FONT_COLOR)
    ax.title.set_color(MI_FONT_COLOR)


def _save(fig, path: Path) -> None:
    fig.patch.set_facecolor("none")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", transparent=True)
    if path.suffix.lower() == ".png":
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
