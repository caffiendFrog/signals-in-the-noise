"""Poster figures use the ISCB palette and keep the claim visible on the axes."""

from pathlib import Path

import pandas as pd

from signals_in_the_noise.analysis.poster_figures import (
    depth_collapse_figure,
    lp_expansion_figure,
    plot_depth_collapse,
    plot_lp_expansion,
)


def test_depth_figure_draws_the_pbs_drop_and_requires_both_versions(tmp_path):
    summary = pd.DataFrame(
        {
            "version": ["original", "depth-matched", "original", "depth-matched"],
            "arm": ["noise_pbs", "noise_pbs", "retained_qc", "retained_qc"],
            "specimen_auc": [0.82, 0.51, 0.74, 0.53],
            "permutation_p": [0.01, 0.48, 0.02, 0.40],
        }
    )
    depths = pd.DataFrame(
        {
            "condition": ["ER+ tumour", "ER+ tumour", "Normal", "Normal"],
            "depth": [8000, 4200, 3000, 2800],
            "included": [True, True, True, False],
        }
    )
    path = tmp_path / "depth.png"
    plot_depth_collapse(summary, depths, path)
    assert path.is_file()
    assert path.with_suffix(".pdf").is_file()
    figure = depth_collapse_figure(summary, depths)
    pbs = [line for line in figure.axes[0].lines if line.get_color() == "#D55E00"]
    assert len(pbs) == 1
    assert list(pbs[0].get_xdata()) == [0, 1]
    assert list(pbs[0].get_ydata()) == [0.82, 0.51]
    assert pbs[0].get_linewidth() > 2.5
    labels = " ".join(text.get_text() for text in figure.axes[0].texts)
    assert "0.82" in labels
    assert "0.51" in labels
    assert figure._suptitle is None

    summary.loc[summary["version"] == "depth-matched", "version"] = "thinned"
    try:
        plot_depth_collapse(summary, depths, tmp_path / "missing.png")
    except ValueError as error:
        assert "depth-matched" in str(error)
    else:
        raise AssertionError("A summary without the depth-matched PBS row should fail.")


def test_depth_matched_labels_do_not_overlap():
    summary = pd.DataFrame(
        {
            "version": ["original", "depth-matched", "original", "depth-matched"],
            "arm": ["retained_qc", "retained_qc", "noise_pbs", "noise_pbs"],
            "specimen_auc": [0.70, 0.00, 0.00, 0.17],
            "permutation_p": [0.04, 1.00, 1.00, 0.71],
        }
    )
    depths = pd.DataFrame(
        {
            "condition": ["ER+ tumour", "Normal"],
            "depth": [4000, 2500],
            "included": [True, True],
        }
    )
    figure = depth_collapse_figure(summary, depths)
    figure.canvas.draw()
    labels = [text for text in figure.axes[0].texts if "p =" in text.get_text()]
    assert len(labels) == 2
    first, second = (text.get_window_extent().expanded(2, 2) for text in labels)
    assert not first.overlaps(second)


def test_low_auc_labels_stay_inside_the_panel():
    summary = pd.DataFrame(
        {
            "version": ["original", "depth-matched"],
            "arm": ["noise_pbs", "noise_pbs"],
            "specimen_auc": [0.00, 0.17],
            "permutation_p": [1.00, 0.71],
        }
    )
    depths = pd.DataFrame(
        {
            "condition": ["ER+ tumour", "Normal"],
            "depth": [6200, 2400],
            "included": [True, True],
        }
    )
    figure = depth_collapse_figure(summary, depths)
    figure.canvas.draw()
    panel = figure.axes[0].get_window_extent()
    bounds = figure.bbox
    for text in figure.axes[0].texts:
        box = text.get_window_extent()
        assert box.y0 >= panel.y0 - 1
        assert box.x0 >= bounds.x0
        assert box.x1 <= bounds.x1


def test_lp_figure_spans_the_cap_and_marks_alpha(tmp_path):
    curve = pd.DataFrame({"k": [0.0, 2.55, 42.7], "p_value": [0.90, 0.70, 0.22]})
    shares = pd.DataFrame(
        {
            "donor": ["B1", "N1"],
            "genotype": ["BRCA1", "WT"],
            "lp_share": [0.14, 0.30],
            "menopause_status": ["Post (oophorectomy)", "Pre"],
        }
    )
    path = tmp_path / "lp.png"
    plot_lp_expansion(curve, shares, path, alpha=0.01, signature_c=2.55, extreme_p=0.22)
    assert path.is_file()
    assert path.with_suffix(".pdf").is_file()
    figure = lp_expansion_figure(curve, shares, 0.01, 2.55, 0.22)
    axis = figure.axes[0]
    assert axis.get_xlim()[1] >= 42.7
    alpha_lines = [line.get_ydata()[0] for line in axis.lines if len(set(line.get_ydata())) == 1]
    assert 0.01 in alpha_lines


def test_palette_constants_are_the_iscb_notebook_values():
    text = Path("notebooks/GSE161529/iscb-scs/12-a-iscb-scs-volcano.ipynb").read_text(
        encoding="utf-8"
    )
    for color in ("#00274C", "#0072B2", "#009E73", "#D55E00", "#E69F00"):
        assert color in text
