"""Tests for signals_in_the_noise.analysis.cell_identity."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.cell_identity import (
    OTHER_TYPE,
    Signature,
    call_cell_types,
    score_signatures,
    signature_from_table,
)

SIGNATURES = {
    "lp": Signature(up=("KIT", "ELF5"), down=("KRT5",)),
    "basal": Signature(up=("KRT5", "ACTA2"), down=()),
}


def _adata(counts: list[list[float]], genes=("KIT", "ELF5", "KRT5", "ACTA2", "GAPDH")) -> AnnData:
    adata = AnnData(np.asarray(counts, dtype=float))
    adata.var_names = list(genes)
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    return adata


def test_signature_from_table_splits_by_fold_change_and_drops_excluded():
    table = pd.DataFrame(
        {
            "Symbol": ["KIT", "KRT14", "KRT5", None],
            "Average log fold-change": [1.2, 0.8, -0.5, 2.0],
        }
    )
    signature = signature_from_table(table, excluded_genes=["KRT14"])
    assert signature == Signature(up=("KIT",), down=("KRT5",))


def test_score_signatures_is_mean_log_difference():
    adata = _adata([[10, 0, 0, 0, 9990]])
    scores = score_signatures(adata, SIGNATURES)
    assert scores.loc["cell_0", "lp"] == pytest.approx(np.log1p(10) / 2)
    assert scores.loc["cell_0", "basal"] == pytest.approx(0.0)


def test_score_signatures_scores_each_cell_independently():
    base = [[5, 3, 0, 1, 100], [0, 0, 8, 6, 90]]
    alone = score_signatures(_adata(base), SIGNATURES)
    with_extra = score_signatures(_adata(base + [[50, 50, 50, 50, 1]]), SIGNATURES)
    pd.testing.assert_frame_equal(alone, with_extra.iloc[:2])


def test_score_signatures_requires_up_genes():
    with pytest.raises(ValueError, match="No up-regulated"):
        score_signatures(_adata([[1, 1, 1, 1, 1]]), {"x": Signature(up=("MISSING",), down=())})


def test_call_cell_types_uses_reference_scale_only():
    scores = pd.DataFrame(
        {"lp": [1.0, 0.0, 0.5, 9.0], "basal": [0.0, 1.0, 0.5, 9.5]},
        index=["a", "b", "c", "rescued"],
    )
    reference = pd.Series([True, True, True, False], index=scores.index)
    calls = call_cell_types(scores, reference)
    assert calls[["a", "b", "c"]].tolist() == ["lp", "basal", OTHER_TYPE]

    shifted = scores.copy()
    shifted.loc["rescued"] = [-50.0, 100.0]
    assert call_cell_types(shifted, reference)[["a", "b", "c"]].tolist() == ["lp", "basal", OTHER_TYPE]
