"""Provisional per-cell epithelial identity from the Lim et al. (2009) signatures.

A stand-in for Step 2 (mapping onto Arm B clusters) that runs before Steps 1-3
exist. It keeps the property Step 2 requires: rescued cells cannot define what
an LP cell looks like.

- Each cell is scored on its own: counts are scaled to a fixed library size and
  the score is a mean-expression difference, so adding cells changes no other
  cell's score.
- Scores are standardised against a reference population (Arm B cells) only.

There is no ambient filter (Step 3), so this is not the identity used in the
main test.
"""

import logging
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from signals_in_the_noise.config import get_resources_path

logger = logging.getLogger(__name__)

LIM_SIGNATURE_FILES: dict[str, str] = {
    "basal": "GSE161529/epithelial_cell_typing/41591_2009_BFnm2000_MOESM13_ESM.xls",
    "lp": "GSE161529/epithelial_cell_typing/41591_2009_BFnm2000_MOESM14_ESM.xls",
    "ml": "GSE161529/epithelial_cell_typing/41591_2009_BFnm2000_MOESM15_ESM.xls",
    "stromal": "GSE161529/epithelial_cell_typing/41591_2009_BFnm2000_MOESM16_ESM.xls",
}
"""Lim et al. (2009) population signatures, relative to the resources directory."""

AMBIENT_CORE_MARKERS: tuple[str, ...] = ("ANKRD30A", "KRT14")
"""Genes in the shared ambient core; unreliable for mature luminal or basal identity."""

OTHER_TYPE = "other"


@dataclass(frozen=True)
class Signature:
    """Up- and down-regulated genes of one population signature."""

    up: tuple[str, ...]
    down: tuple[str, ...]


def signature_from_table(table: pd.DataFrame, *, excluded_genes: Iterable[str] = ()) -> Signature:
    """Split a Lim et al. table (``Symbol``, ``Average log fold-change``) into up and down genes."""
    table = table.loc[:, ["Symbol", "Average log fold-change"]].dropna()
    table = table.loc[~table["Symbol"].isin(set(excluded_genes))]
    up = table.loc[table["Average log fold-change"] >= 0, "Symbol"].unique()
    down = table.loc[table["Average log fold-change"] < 0, "Symbol"].unique()
    return Signature(up=tuple(up), down=tuple(down))


def load_lim_signatures(
    *, excluded_genes: Iterable[str] = AMBIENT_CORE_MARKERS
) -> dict[str, Signature]:
    """Read the :data:`LIM_SIGNATURE_FILES`, dropping ``excluded_genes``."""
    excluded_genes = tuple(excluded_genes)
    return {
        name: signature_from_table(
            pd.read_excel(get_resources_path(path)), excluded_genes=excluded_genes
        )
        for name, path in LIM_SIGNATURE_FILES.items()
    }


def score_signatures(
    adata: AnnData, signatures: Mapping[str, Signature], *, target_sum: float = 1e4
) -> pd.DataFrame:
    """Per-cell mean log expression of up genes minus mean of down genes.

    Args:
        adata: Raw counts with gene symbols as ``var_names``.
        signatures: Signature name to :class:`Signature`.
        target_sum: Library size every cell is scaled to before ``log1p``.

    Returns:
        One row per cell (indexed by ``obs_names``), one column per signature.

    Raises:
        ValueError: If a signature has no up-regulated genes in ``adata``.
    """
    counts = sparse.csr_matrix(adata.X, dtype=float)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals), where=totals > 0)
    logged = sparse.diags(scale) @ counts
    logged.data = np.log1p(logged.data)
    logged = logged.tocsc()

    genes = pd.Index(adata.var_names)

    def mean_expression(symbols: tuple[str, ...]) -> np.ndarray:
        columns = np.flatnonzero(genes.isin(symbols))
        if columns.size == 0:
            return np.zeros(adata.n_obs)
        return np.asarray(logged[:, columns].mean(axis=1)).ravel()

    scores = {}
    for name, signature in signatures.items():
        if not genes.isin(signature.up).any():
            raise ValueError(f"No up-regulated {name} genes found in adata.var_names.")
        scores[name] = mean_expression(signature.up) - mean_expression(signature.down)
    return pd.DataFrame(scores, index=adata.obs_names.astype(str))


def call_cell_types(scores: pd.DataFrame, reference: pd.Series) -> pd.Series:
    """Assign each cell the signature with the highest score, standardised on reference cells.

    Scores are z-scored with the mean and SD of ``reference`` cells only, so
    non-reference cells cannot shift the calls. Cells whose best z-score is not
    positive are :data:`OTHER_TYPE`.

    Args:
        scores: Output of :func:`score_signatures`.
        reference: Boolean mask aligned with ``scores`` selecting reference cells.
    """
    reference_scores = scores.loc[reference.to_numpy()]
    z = (scores - reference_scores.mean()) / reference_scores.std(ddof=0)
    calls = z.idxmax(axis=1)
    calls[z.max(axis=1) <= 0] = OTHER_TYPE
    return calls
