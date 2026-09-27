"""Quality-control rules that decide which cells are "noise".

The published GSE161529 rule uses per-sample thresholds (Pal et al. lowered the
gene threshold for low-coverage samples). The alternatives here let the same
flagging logic run under a study-wide uniform rule or a per-sample adaptive
(MAD) rule, and after equalizing sequencing depth by binomial thinning.
"""

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, fields

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from scipy.stats import median_abs_deviation

from signals_in_the_noise.utils.log import get_logger

logger = get_logger(__name__)

QC_METRICS: tuple[str, ...] = ("pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts")
"""Cell-level QC metrics used to define PBS subtypes."""

NOISE_FLAG_COLUMNS: tuple[str, ...] = (
    "is_low_num_genes",
    "is_high_num_genes",
    "is_high_mito",
    "is_high_total_count",
)
"""Individual QC flags; a cell is noise when any of them is set."""

DEFAULT_PERCENT_TOP: tuple[int, ...] = (50, 100, 200, 500)
"""Scanpy's default ``percent_top`` for :func:`scanpy.pp.calculate_qc_metrics`."""

POPULATIONS: tuple[str, ...] = ("all", "retained", "noise")
"""Cell populations relative to a noise flag: every cell, cells kept by QC, cells discarded."""


def select_population(
    obs: pd.DataFrame, population: str, noise_column: str = "is_noise"
) -> pd.DataFrame:
    """Return the ``"all"``, ``"retained"`` (flag 0) or ``"noise"`` (flag 1) cells of ``obs``."""
    if population == "all":
        return obs
    if population == "retained":
        return obs[obs[noise_column] == 0]
    if population == "noise":
        return obs[obs[noise_column] == 1]
    raise ValueError(f"Unknown population {population!r}; expected one of {POPULATIONS}.")


@dataclass(frozen=True)
class QcThresholds:
    """Absolute QC cutoffs for one sample, in the units used by Pal et al.

    Attributes:
        mito_upper: Maximum mitochondrial read fraction (0–1).
        genes_lower: Cells with ``n_genes_by_counts <= genes_lower`` are noise.
        genes_upper: Cells with ``n_genes_by_counts > genes_upper`` are noise.
        total_upper: Cells with ``total_counts >= total_upper`` are noise.
    """

    mito_upper: float
    genes_lower: float
    genes_upper: float
    total_upper: float

    UNS_KEYS = {
        "mito_upper": "qc_mito_upper",
        "genes_lower": "qc_genes_lower",
        "genes_upper": "qc_genes_upper",
        "total_upper": "qc_total_upper",
    }

    @classmethod
    def from_uns(cls, uns: Mapping) -> "QcThresholds":
        """Read the published per-sample thresholds from ``adata.uns``."""
        return cls(**{name: float(uns[key]) for name, key in cls.UNS_KEYS.items()})


def annotate_qc_metrics(
    adata: AnnData, *, percent_top: Iterable[int] | None = DEFAULT_PERCENT_TOP
) -> AnnData:
    """Flag mitochondrial genes and write scanpy QC metrics to ``adata`` in place."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None if percent_top is None else list(percent_top),
        inplace=True,
    )
    return adata


def flag_noise(obs: pd.DataFrame, thresholds: QcThresholds) -> pd.DataFrame:
    """Return the four QC flags plus ``is_noise`` (all ``int``) for each cell.

    Requires ``n_genes_by_counts``, ``pct_counts_mt`` and ``total_counts`` in ``obs``.
    """
    flags = pd.DataFrame(
        {
            "is_low_num_genes": obs["n_genes_by_counts"] <= thresholds.genes_lower,
            "is_high_num_genes": obs["n_genes_by_counts"] > thresholds.genes_upper,
            "is_high_mito": obs["pct_counts_mt"] / 100 > thresholds.mito_upper,
            "is_high_total_count": obs["total_counts"] >= thresholds.total_upper,
        },
        index=obs.index,
    )
    flags["is_noise"] = flags.any(axis=1)
    return flags.astype(int)


def modal_thresholds(thresholds: Iterable[QcThresholds]) -> QcThresholds:
    """Collapse per-sample thresholds into one study-wide rule using each field's mode.

    Ties resolve to the smallest modal value.
    """
    frame = pd.DataFrame([vars(t) for t in thresholds])
    if frame.empty:
        raise ValueError("At least one QcThresholds is required.")
    return QcThresholds(**{f.name: float(frame[f.name].mode().min()) for f in fields(QcThresholds)})


def _mad_bound(values: pd.Series, n_mads: float, sign: int) -> float:
    mad = median_abs_deviation(values, scale="normal", nan_policy="omit")
    if not np.isfinite(mad) or mad == 0:
        return float(sign * np.inf)
    return float(np.nanmedian(values) + sign * n_mads * mad)


def mad_thresholds(obs: pd.DataFrame, n_mads: float = 3.0) -> QcThresholds:
    """Derive per-sample adaptive thresholds, ``median ± n_mads * MAD``.

    Gene and library-size bounds are computed on the log1p scale and mapped
    back to counts; the mitochondrial bound is computed on the raw fraction.
    A metric with zero MAD gets an infinite bound, so it flags nothing.
    """
    log_genes = obs["log1p_n_genes_by_counts"]
    return QcThresholds(
        mito_upper=_mad_bound(obs["pct_counts_mt"] / 100, n_mads, +1),
        genes_lower=float(np.expm1(_mad_bound(log_genes, n_mads, -1))),
        genes_upper=float(np.expm1(_mad_bound(log_genes, n_mads, +1))),
        total_upper=float(np.expm1(_mad_bound(obs["log1p_total_counts"], n_mads, +1))),
    )


def relabel_noise(
    specimen_obs: Mapping[str, pd.DataFrame],
    rule: Callable[[pd.DataFrame], QcThresholds],
    column: str,
) -> dict[str, pd.DataFrame]:
    """Return copies of each specimen's ``obs`` with ``column`` set to the noise flag from ``rule``.

    Args:
        specimen_obs: Specimen identifier to cell-level QC frame.
        rule: Called with one specimen's ``obs``; returns that specimen's thresholds.
        column: Name of the new ``int`` noise column.
    """
    relabelled = {}
    for specimen_id, obs in specimen_obs.items():
        obs = obs.copy()
        obs[column] = flag_noise(obs, rule(obs))["is_noise"]
        relabelled[specimen_id] = obs
    return relabelled


def thin_counts(X, fraction: float, rng: np.random.Generator) -> sparse.csr_matrix:
    """Binomially thin a raw count matrix so each UMI survives with ``fraction``.

    Thinning mimics sequencing the same library less deeply, unlike per-cell
    downsampling which would erase between-cell depth differences.

    Raises:
        ValueError: If ``fraction`` is outside ``(0, 1]`` or ``X`` is not integer counts.
    """
    if not 0 < fraction <= 1:
        raise ValueError(f"fraction must be in (0, 1], got {fraction}.")
    thinned = sparse.csr_matrix(X, copy=True)
    counts = np.rint(thinned.data)
    if not np.allclose(counts, thinned.data):
        raise ValueError("thin_counts requires raw integer counts.")
    if fraction < 1:
        thinned.data = rng.binomial(counts.astype(np.int64), fraction).astype(thinned.dtype)
        thinned.eliminate_zeros()
    return thinned


def library_depth(adata: AnnData, mask: np.ndarray | pd.Series | None = None) -> float:
    """Median total counts per cell, optionally over a subset of cells."""
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    if mask is not None:
        totals = totals[np.asarray(mask, dtype=bool)]
    return float(np.median(totals))


def depth_matched_qc_obs(
    adata: AnnData,
    target_depth: float,
    rng: np.random.Generator,
    *,
    reference_mask: np.ndarray | pd.Series | None = None,
) -> tuple[pd.DataFrame, float]:
    """Thin ``adata`` to ``target_depth`` median counts and recompute QC metrics.

    Samples already at or below the target are left unthinned.

    Args:
        adata: Sample with raw counts in ``X``.
        target_depth: Desired median total counts over ``reference_mask`` cells.
        rng: Random generator used for thinning.
        reference_mask: Cells whose median depth is matched; defaults to all.

    Returns:
        ``(obs, fraction)`` — freshly computed QC metrics indexed like
        ``adata.obs`` and the thinning fraction that was applied.
    """
    fraction = min(1.0, target_depth / library_depth(adata, reference_mask))
    thinned = AnnData(
        X=thin_counts(adata.X, fraction, rng),
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    annotate_qc_metrics(thinned, percent_top=None)
    logger.debug("thinned sample by fraction %.3f", fraction)
    return thinned.obs, fraction
