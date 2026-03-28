import logging
import math

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData

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
        df = pd.DataFrame(
            adata[:, genes].X.toarray(),
            index=adata.obs_names,
            columns=genes,
        )
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

    for ax, heatmap_df in zip(axes, panel_dfs):
        sns.heatmap(heatmap_df, cmap=cmap, ax=ax, cbar_kws={"pad": 0.02})

    axes[-1].set_xlabel("Gene Set")
    axes[-1].set_ylabel("Specimen")

    return list(axes)
