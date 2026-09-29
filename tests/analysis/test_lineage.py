"""Retained-cell lineage labels. The share has to ignore removed cells."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from signals_in_the_noise.analysis.lineage import retained_lineage_labels

_BARCODE = "AAACCTGAGCAAATCA-1"
_SIGNATURES = Path("resources/GSE161529/epithelial_cell_typing")
_FILES = {
    "basal": "41591_2009_BFnm2000_MOESM13_ESM.xls",
    "lp": "41591_2009_BFnm2000_MOESM14_ESM.xls",
    "ml": "41591_2009_BFnm2000_MOESM15_ESM.xls",
    "stromal": "41591_2009_BFnm2000_MOESM16_ESM.xls",
}


def test_lp_share_is_the_retained_fraction_and_removed_cells_are_not_scored():
    cells = pd.DataFrame(
        {
            "donor": ["A", "A", "A"],
            "barcode": [_BARCODE, "AAACCTGAGCAAATCA-2", "AAACCTGAGCAAATCA-3"],
            "is_noise": [0, 0, 1],
        }
    )
    counts = _matrix(["lib_" + _BARCODE, "AAACCTGAGCAAATCA-2", "AAACCTGAGCAAATCA-3"])
    donors = pd.DataFrame({"donor": ["A"], "specimen_id": ["s-a"]}).set_index("specimen_id")
    seen = {}

    def scorer(adata: AnnData) -> AnnData:
        seen["n_obs"] = adata.n_obs
        labeled = adata.copy()
        labeled.obs["predicted_type"] = np.where(labeled.obs_names == _BARCODE, "lp", "basal")
        return labeled

    before = counts.X.copy()
    result = retained_lineage_labels(cells, {"s-a": counts}, donors, scorer=scorer)
    assert seen["n_obs"] == 2
    assert np.array_equal(counts.X, before)
    assert result["shares"].loc["A", "lp_share"] == 0.5
    assert set(result["labels"]["barcode"]) == {_BARCODE, "AAACCTGAGCAAATCA-2"}


@pytest.mark.filterwarnings("ignore:Bitwise inversion:DeprecationWarning")
def test_lim_signatures_call_a_progenitor_cell_lp_and_a_flat_background_cell_other():
    lp_genes, other_genes = _signature_genes()
    background = [f"CTRL{i}" for i in range(50)]
    genes = lp_genes + other_genes + background
    matrix = np.zeros((2, len(genes)))
    matrix[0, : len(lp_genes)] = 100
    matrix[0, len(genes) - 50 :] = 1
    matrix[1, len(genes) - 50 :] = 10
    counts = AnnData(matrix, obs=pd.DataFrame(index=[_BARCODE, "AAACCTGAGCAAATCA-2"]))
    counts.var_names = genes
    cells = pd.DataFrame(
        {
            "donor": ["A", "A"],
            "barcode": [_BARCODE, "AAACCTGAGCAAATCA-2"],
            "is_noise": [0, 0],
        }
    )
    donors = pd.DataFrame({"donor": ["A"], "specimen_id": ["s-a"]}).set_index("specimen_id")
    result = retained_lineage_labels(cells, {"s-a": counts}, donors)
    types = result["labels"].set_index("barcode")["predicted_type"]
    assert types.loc[_BARCODE] == "lp"
    assert types.loc["AAACCTGAGCAAATCA-2"] == "other"
    assert result["shares"].loc["A", "lp_share"] == 0.5


def _matrix(barcodes: list[str]) -> AnnData:
    counts = AnnData(np.ones((len(barcodes), 2)))
    counts.obs_names = barcodes
    counts.var_names = ["G1", "G2"]
    return counts


def _signature_genes() -> tuple[list[str], list[str]]:
    upregulated = {}
    symbols = {}
    for name, filename in _FILES.items():
        table = pd.read_excel(_SIGNATURES / filename)
        table = table.loc[:, ["Symbol", "Average log fold-change"]].dropna()
        symbols[name] = set(table["Symbol"].astype(str))
        upregulated[name] = set(
            table.loc[table["Average log fold-change"] >= 0, "Symbol"].astype(str)
        )
    others = set().union(*(genes for name, genes in upregulated.items() if name != "lp"))
    lp_only = sorted(upregulated["lp"] - others)
    if len(lp_only) < 5:
        raise AssertionError("The LP signature has too few genes of its own to score against.")
    rest = sorted(set().union(*symbols.values()) - set(lp_only))
    return lp_only, rest
