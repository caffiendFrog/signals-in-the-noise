import logging
import math

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from anndata import AnnData

from signals_in_the_noise.analysis.statistics import fdr_to_stars

logger = logging.getLogger(__name__)


def get_figure_axes(
    n_plots: int,
    num_cols: int = 2,
    subplot_size: tuple[int, int] = (6, 4),
    share_x: bool = False,
    share_y: bool = False,
    super_title: str = None,
) -> tuple:
    """Return a figure and flat axes array sized for the requested number of plots.

    Args:
        n_plots: Number of plots to display.
        num_cols: Number of columns in the subplot grid.
        subplot_size: Width and height of each individual subplot in inches.
        share_x: Whether subplots share the x-axis.
        share_y: Whether subplots share the y-axis.
        super_title: Optional figure-level title.

    Returns:
        Tuple of (Figure, ndarray of Axes) with unused axes turned off.
    """
    rows = math.ceil(n_plots / num_cols)
    width = num_cols * subplot_size[0]
    height = rows * subplot_size[1]

    fig, axes = plt.subplots(rows, num_cols, figsize=(width, height), sharex=share_x, sharey=share_y)

    # massage axes to guarantee it is a flat array regardless of grid shape
    axes = axes.flatten() if isinstance(axes, np.ndarray) else np.array([axes])

    # turn "off" any plots that are not used, subplot returns a grid
    for i in range(n_plots, len(axes)):
        axes[i].axis("off")

    if super_title:
        fig.suptitle(super_title, fontsize=16)
        fig.subplots_adjust(top=0.92)

    return fig, axes[:n_plots]


def plot_pathway_heatmap(
    adata: AnnData,
    pathway_name: str,
    pathway_genes: list[str],
    n_panels: int = 2,
    cmap: str = "rocket",
    subplot_size: tuple[int, int] = (60, 5),
    axes: list[matplotlib.axes.Axes] | None = None,
) -> list[matplotlib.axes.Axes]:
    """Plot mean raw gene expression per specimen for a single pathway as a heatmap.

    Genes are split evenly across ``n_panels`` vertically stacked sub-panels so
    that wide gene sets remain legible.  Only genes present in ``adata.var_names``
    are included; absent genes are silently dropped.

    Args:
        adata: AnnData object with raw count matrix and an ``obs`` column named
            ``specimen_id``.
        pathway_name: Display name for the pathway used in the figure super-title.
        pathway_genes: Ordered list of gene identifiers belonging to the pathway.
        n_panels: Number of vertical sub-panels to split the gene set across.
            Defaults to 2.
        cmap: Seaborn / Matplotlib colormap name. Defaults to ``'rocket'``.
        subplot_size: ``(width, height)`` in inches for each individual sub-panel.
            Defaults to ``(60, 5)``.
        axes: Pre-existing list of ``n_panels`` axes to draw on.  When ``None``
            a new figure is created internally.

    Returns:
        List of axes, one per panel, with the heatmaps drawn.
    """
    actual_genes = [g for g in pathway_genes if g in adata.var_names]
    gene_splits = np.array_split(actual_genes, n_panels)

    panel_dfs: list[pd.DataFrame] = []
    for genes in gene_splits:
        if len(genes) == 0:
            continue
        X = adata[:, genes].X
        values = X.toarray() if sp.issparse(X) else np.asarray(X)
        df = pd.DataFrame(values, index=adata.obs_names, columns=genes)
        df["specimen_id"] = adata.obs["specimen_id"].values
        panel_dfs.append(df.groupby("specimen_id").mean())

    if axes is None:
        _, axes = get_figure_axes(
            n_panels,
            num_cols=1,
            subplot_size=subplot_size,
            share_y=True,
            super_title=f"Mean Raw Expressions of Genes Grouped by Specimen ID\n{pathway_name}",
        )

    if len(axes) != n_panels:
        raise ValueError(
            f"Expected {n_panels} axes to match n_panels, got {len(axes)}."
        )

    for ax, heatmap_df in zip(axes, panel_dfs):
        sns.heatmap(heatmap_df, cmap=cmap, ax=ax, cbar_kws={"pad": 0.02})

    axes[-1].set_xlabel("Gene Set")
    axes[-1].set_ylabel("Specimen")

    return list(axes)


def plot_score_heatmap(
    adata: AnnData,
    score_columns: list[str],
    title: str = "Raw Gene Set Scores Grouped by Specimen ID",
    cmap: str = "rocket",
    figsize: tuple[float, float] = (60, 4),
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.axes.Axes:
    """Plot a transposed per-cell gene-set score heatmap sorted by specimen ID.

    Cells are sorted by ``specimen_id`` so that cells from the same specimen
    appear together on the x-axis.  The score matrix is transposed so that
    score names appear on the y-axis and individual cells on the x-axis.

    Args:
        adata: AnnData object whose ``obs`` must contain all columns listed in
            ``score_columns`` plus a ``specimen_id`` column.
        score_columns: Names of the score columns in ``adata.obs`` to display.
        title: Figure title string. Defaults to a generic scores-by-specimen label.
        cmap: Colormap name passed to :func:`seaborn.heatmap`.
            Defaults to ``'rocket'``.
        figsize: ``(width, height)`` in inches used when creating a new figure.
            Defaults to ``(60, 4)``.
        ax: Existing axes to draw on.  When ``None`` a new figure is created.

    Returns:
        The axes with the heatmap drawn.
    """
    scores = adata.obs[score_columns]
    specimen_ids = adata.obs["specimen_id"]
    sorted_idx = specimen_ids.sort_values().index
    sorted_scores = scores.loc[sorted_idx]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    sns.heatmap(sorted_scores.T, cmap=cmap, square=True, xticklabels=False, ax=ax)
    ax.set_aspect("auto")
    ax.set_title(title)
    ax.set_xlabel("Cells")
    ax.set_ylabel("Gene Set")

    return ax


def plot_gsea_nes_heatmap(
    gsea_results_df: pd.DataFrame,
    cmap: str = "coolwarm",
    label_fontsize: int = 10,
    fig_width: float = 6.0,
    bar_height: float = 0.5,
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.axes.Axes:
    """Plot a GSEA NES heatmap with FDR significance annotations on the y-axis.

    Gene set terms are shown on the y-axis; the single-column heatmap encodes
    the Normalized Enrichment Score (NES) by colour.  Significance stars from
    :func:`~signals_in_the_noise.analysis.statistics.fdr_to_stars` are
    appended to each term label.

    Args:
        gsea_results_df: DataFrame with at least three columns: ``'Term'``
            (gene set name), ``'NES'`` (normalized enrichment score), and
            ``'FDR q-val'`` (FDR-corrected q-value).  Typically the
            ``res2d`` attribute of a :mod:`gseapy` result object.
        cmap: Diverging colormap name centered at zero.  Defaults to
            ``'coolwarm'``.
        label_fontsize: Font size for y-axis tick labels.  Defaults to ``10``.
        fig_width: Width of the figure in inches when ``ax`` is ``None``.
            Defaults to ``6.0``.
        bar_height: Height per gene-set row in inches, used to compute the
            figure height when ``ax`` is ``None``.  Defaults to ``0.5``.
        ax: Existing axes to draw on.  When ``None`` a new figure is created
            with height ``len(terms) * bar_height``.

    Returns:
        The axes with the heatmap drawn.
    """
    heatmap_data = gsea_results_df.pivot_table(index="Term", values="NES")
    heatmap_data = heatmap_data.apply(pd.to_numeric, errors="coerce")

    fdr_map = gsea_results_df.set_index("Term")["FDR q-val"].apply(fdr_to_stars).to_dict()
    annotated_labels = [term + fdr_map.get(term, "") for term in heatmap_data.index]

    if ax is None:
        _, ax = plt.subplots(figsize=(fig_width, len(heatmap_data) * bar_height))

    sns.heatmap(heatmap_data, annot=True, cmap=cmap, center=0, ax=ax)
    ax.set_yticklabels(annotated_labels, rotation=0, fontsize=label_fontsize)
    ax.set_title("GSEA Normalized Enrichment Score (NES) with Significance")

    return ax


def plot_noise_subtype_comparison(
    norm_df: pd.DataFrame,
    figsize: tuple[float, float] = (13, 6),
    signal_panel_ratio: float = 3.0,
    axes: tuple[matplotlib.axes.Axes, matplotlib.axes.Axes] | None = None,
) -> tuple[matplotlib.axes.Axes, matplotlib.axes.Axes]:
    """Plot a paired bar chart comparing noise-subtype proportions across cancer types.

    The left panel shows the percentage of noise cells that fall into each
    biological-signal subtype (damaged, dormant, multifunction) for every
    cancer type.  The right panel shows the percentage of all cells that are
    noise cells.  A shared legend is placed to the right of the figure.

    Args:
        norm_df: Normalised summary DataFrame as returned by
            :func:`~signals_in_the_noise.analysis.noise_phenotypes.aggregate_noise_subtypes_by_cancer_type`.
            Must have an index of annotated cancer-type labels and columns
            ``'damaged'``, ``'dormant'``, ``'multifunction'``, and ``'noise'``.
        figsize: ``(width, height)`` in inches used when creating a new figure.
            Defaults to ``(13, 6)``.
        signal_panel_ratio: Width ratio of the signal panel relative to the
            noise panel.  Defaults to ``3.0``.
        axes: Pre-existing ``(ax1, ax2)`` tuple to draw on.  When ``None`` a
            new figure is created.

    Returns:
        Tuple of ``(ax1, ax2)`` — the signal panel and the noise panel.
    """
    signal_cols = ["damaged", "dormant", "multifunction"]

    long_signals = (
        norm_df[signal_cols]
        .reset_index()
        .rename(columns={"index": "cancer type"})
        .melt(id_vars="cancer type", var_name="potential signal", value_name="percentage of noise")
    )
    long_noise = (
        norm_df[["noise"]]
        .reset_index()
        .rename(columns={"index": "cancer type"})
        .melt(id_vars="cancer type", var_name="potential signal", value_name="percentage of cells")
    )

    if axes is None:
        _, (ax1, ax2) = plt.subplots(
            1, 2,
            figsize=figsize,
            gridspec_kw={"width_ratios": [signal_panel_ratio, 1]},
            sharey=True,
        )
    else:
        ax1, ax2 = axes

    sns.barplot(data=long_signals, x="potential signal", y="percentage of noise", hue="cancer type", ax=ax1)
    ax1.set_title('Potential Biological Signals In "Noise" Cells')
    ax1.set_xlabel("")
    ax1.set_ylabel("percentage")
    ax1.get_legend().remove()

    sns.barplot(data=long_noise, x="potential signal", y="percentage of cells", hue="cancer type", ax=ax2)
    ax2.set_title('"Noise" Cells')
    ax2.set_xlabel("")
    ax2.set_ylabel("percentage")
    ax2.get_legend().remove()

    handles, labels = ax1.get_legend_handles_labels()
    ax1.figure.legend(
        handles, labels,
        title="Breast Cancer Type",
        loc="center right",
        bbox_to_anchor=(1.05, 0.70),
    )
    ax1.figure.suptitle(
        'Comparative Analysis of "Noise" cells Among Breast Cancer Subtypes',
        fontsize=18,
    )

    return ax1, ax2


def plot_empty_cell_violin_comparison(
    noise_adata: AnnData,
    all_adata: AnnData,
    groupby: str,
    metric: str = "zero_mito",
    ylabel: str = "Empty State (1 = empty, 0 = not empty)",
    subplot_size: tuple[int, int] = (7, 8),
    axes: list[matplotlib.axes.Axes] | None = None,
) -> list[matplotlib.axes.Axes]:
    """Plot side-by-side violin distributions of an empty-cell metric by group.

    Renders two panels: the left shows the distribution for noise cells only;
    the right shows the distribution across all cells (noise + real combined).
    Both panels are grouped by ``groupby`` and share the same y-axis metric.

    Args:
        noise_adata: AnnData containing only noise cells.
        all_adata: AnnData containing all cells (noise and real combined).
        groupby: ``obs`` column name to use as the grouping variable on the
            x-axis (e.g. ``'hormonal_status'``, ``'cancer_type'``).
        metric: ``obs`` column name of the binary empty-cell indicator to
            visualise.  Defaults to ``'zero_mito'``.
        ylabel: Y-axis label shared across both panels.  Defaults to
            ``'Empty State (1 = empty, 0 = not empty)'``.
        subplot_size: ``(width, height)`` in inches for each sub-panel.
            Defaults to ``(7, 8)``.
        axes: Pre-existing list of two axes to draw on.  When ``None`` a new
            figure is created via :func:`get_figure_axes`.

    Returns:
        List of the two axes ``[noise_ax, all_ax]``.
    """
    violin_kwargs = {
        "show": False,
        "keys": metric,
        "rotation": 90,
        "ylabel": ylabel,
    }

    if axes is None:
        _, axes = get_figure_axes(2, num_cols=2, subplot_size=subplot_size)
        axes[0].figure.suptitle("Distribution of Zero Mitochondria")

    if len(axes) != 2:
        raise ValueError(
            f"Expected exactly 2 axes (noise panel and all-cells panel), got {len(axes)}."
        )

    axes[0].set_title('"Noise" cells')
    axes[1].set_title("All cells")

    sc.pl.violin(noise_adata, groupby=groupby, ax=axes[0], **violin_kwargs)
    sc.pl.violin(all_adata, groupby=groupby, ax=axes[1], **violin_kwargs)

    return list(axes)


def plot_gene_signature_score_distribution(
    scores: pd.Series,
    signature_name: str,
    bins: int = 50,
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.axes.Axes:
    """Plot a histogram with KDE and mean line for a single gene signature score.

    Args:
        scores: Series of per-cell gene signature score values, typically
            extracted from ``adata.obs[f'score_{signature_name}']``.
        signature_name: Human-readable name of the gene signature used to
            generate axis labels and the plot title (e.g. ``'basal'``).
        bins: Number of histogram bins.  Defaults to ``50``.
        ax: Existing axes to draw on.  When ``None`` a new figure is created.

    Returns:
        The axes with the histogram, KDE, and mean line drawn.
    """
    if ax is None:
        _, ax = plt.subplots()

    sns.histplot(scores, kde=True, bins=bins, ax=ax)
    ax.axvline(x=scores.mean(), color="red", linestyle="--", label="Mean")
    ax.set_xlabel(f"{signature_name} score")
    ax.set_ylabel("Cell count")
    ax.set_title(f"Distribution of {signature_name} gene signature score")
    ax.legend()

    return ax
