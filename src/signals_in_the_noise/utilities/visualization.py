import math
import matplotlib.pyplot as plt
import numpy as np

def get_figure_axes(n_plots, num_cols=2, subplot_size=(6, 4), share_x=False, share_y=False, super_title=None):
    """
    Given the number of plots to display, returns a figure and corresponding axes to use for visualizations.

    :param n_plots:
    :param num_cols:
    :param subplot_size:
    :param share_x:
    :param share_y:
    :param super_title:
    :return:
    """

    # calculate rows and dimensions
    rows = math.ceil(n_plots / num_cols)
    width = num_cols * subplot_size[0]
    height = rows * subplot_size[1]

    # create figure and axes
    fig, axes = plt.subplots(rows, num_cols, figsize=(width, height), sharex=share_x, sharey=share_y)

    # massage axes to guarantee it is a flat list
    axes = axes.flatten() if isinstance(axes, (list, np.ndarray)) or axes.ndim > 0 else [axes]

    # turn "off" any plots that are not used, subplot returns a grid
    for i in range(n_plots, len(axes)):
        axes[i].axis('off')

    # add super title if one exists
    if super_title:
        fig.suptitle(super_title, fontsize=16)
        fig.subplots_adjust(top=0.92)

    return fig, axes[:n_plots]
