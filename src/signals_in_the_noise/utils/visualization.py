import math

import matplotlib.pyplot as plt
import numpy as np


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

    # massage axes to guarantee it is a flat list
    axes = axes.flatten() if isinstance(axes, (list, np.ndarray)) or axes.ndim > 0 else [axes]

    # turn "off" any plots that are not used, subplot returns a grid
    for i in range(n_plots, len(axes)):
        axes[i].axis("off")

    if super_title:
        fig.suptitle(super_title, fontsize=16)
        fig.subplots_adjust(top=0.92)

    return fig, axes[:n_plots]
