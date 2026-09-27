"""Execute every qc-confound notebook end to end on a small synthetic cohort.

The real GSE161529 loader is replaced by synthetic specimens with count
matrices and Pal-style metadata (including a patient with two same-titled
specimens and a patient with both a Normal and an ER+ specimen). Permutation
and draw counts are shrunk so the whole series runs in seconds; the goal is to
catch errors in the notebook code, not to reproduce any result.
"""

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from anndata import AnnData  # noqa: E402

import signals_in_the_noise.config as config  # noqa: E402
import signals_in_the_noise.modeling.experiments as experiments  # noqa: E402
import signals_in_the_noise.preprocessing.gse161529 as gse161529  # noqa: E402
from signals_in_the_noise.preprocessing.qc import (  # noqa: E402
    QcThresholds,
    annotate_qc_metrics,
    flag_noise,
)

NOTEBOOK_DIR = Path(__file__).parents[2] / "notebooks" / "GSE161529" / "qc-confound"
NOTEBOOKS = sorted(NOTEBOOK_DIR.glob("*.ipynb"))

# (filename, cancer_type, patient)
SPECIMENS = [
    ("GSM1_N-0019-total", "Normal", "0019"),
    ("GSM2_N-0092-total", "Normal", "0092"),
    ("GSM3_N-0123-total", "Normal", "0123"),
    ("GSM4_N-MH0064-total", "Normal", "0064"),
    ("GSM5_ER-0001", "ER+ tumour", "0001"),
    ("GSM6_ER-0029-7C", "ER+ tumour", "0029"),
    ("GSM7_ER-0029-9C", "ER+ tumour", "0029"),
    ("GSM8_ER-MH0064-T", "ER+ tumour", "0064"),
    ("GSM9_ER-0125", "ER+ tumour", "0125"),
    ("GSM10_HER2-0308", "HER2+ tumour", "0308"),
]


def _synthetic_specimen(index: int, filename: str, cancer_type: str, patient: str) -> AnnData:
    rng = np.random.default_rng(index)
    n_cells, n_genes, n_mt = 700, 120, 4
    depth_scale = 1.0 if cancer_type == "Normal" else 1.6
    damaged = rng.random(n_cells) < 0.2
    library = rng.lognormal(mean=np.log(400 * depth_scale), sigma=0.6, size=n_cells)
    library[damaged] /= 6
    rates = rng.dirichlet(np.ones(n_genes))
    mito = np.clip(rng.beta(2, 20, size=n_cells), 0.001, 0.9)
    mito[damaged] = rng.uniform(0.3, 0.8, size=damaged.sum())
    rates_cells = np.tile(rates, (n_cells, 1))
    rates_cells[:, :n_mt] *= (mito / rates[:n_mt].sum())[:, None]
    rates_cells /= rates_cells.sum(axis=1, keepdims=True)
    X = sp.csr_matrix(rng.poisson(library[:, None] * rates_cells).astype(np.float32))

    adata = AnnData(X)
    genes = [f"MT-{i}" for i in range(n_mt)] + [f"GENE_{i}" for i in range(n_genes - n_mt)]
    adata.var_names = genes
    adata.obs_names = [f"{filename}_{i}" for i in range(n_cells)]
    adata.obs["adata-filename"] = f"{filename}.h5ad"
    annotate_qc_metrics(adata, percent_top=None)

    genes = adata.obs["n_genes_by_counts"]
    thresholds = QcThresholds(
        mito_upper=float(np.quantile(adata.obs["pct_counts_mt"] / 100, 0.9)),
        genes_lower=float(np.quantile(genes, 0.3 if index % 3 else 0.2)),
        genes_upper=float(np.quantile(genes, 0.99)),
        total_upper=float(np.quantile(adata.obs["total_counts"], 0.99)),
    )
    for column, values in flag_noise(adata.obs, thresholds).items():
        adata.obs[column] = values
    adata.uns.update(
        {
            "title": f"{cancer_type} Total cells from Patient {patient}",
            "cancer_type": cancer_type,
            "cell_population": "Total",
            "menopause_status": "Pre" if index % 2 else "Post",
            "num_genes_before": int((X.sum(axis=0) > 0).sum()),
            **{key: getattr(thresholds, name) for name, key in QcThresholds.UNS_KEYS.items()},
        }
    )
    return adata


class _FakeGSE161529:
    def __init__(self):
        self.objects = {
            f"{name}.h5ad": _synthetic_specimen(i, name, cancer_type, patient)
            for i, (name, cancer_type, patient) in enumerate(SPECIMENS)
        }


def _code_cells(path: Path) -> list[str]:
    cells = json.loads(path.read_text(encoding="utf-8"))["cells"]
    return ["".join(c["source"]) for c in cells if c["cell_type"] == "code"]


@pytest.fixture()
def synthetic_study(monkeypatch, tmp_path):
    monkeypatch.setattr(gse161529, "GSE161529", _FakeGSE161529)
    monkeypatch.setattr(config, "get_data_path", lambda subpath=None: tmp_path / (subpath or ""))
    for key, value in {"min_reference_cells": 5, "n_minority": 3, "n_draws": 2}.items():
        monkeypatch.setitem(experiments.PUBLISHED_PROTOCOL, key, value)
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: plt.close("all"))
    return tmp_path


def test_notebooks_are_present():
    assert [p.name[:2] for p in NOTEBOOKS] == ["00", "01", "02", "03", "04"]


@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=[p.stem for p in NOTEBOOKS])
def test_notebook_executes_on_synthetic_cohort(notebook, synthetic_study):
    namespace: dict = {"__name__": "__main__"}
    for index, source in enumerate(_code_cells(notebook)):
        source = re.sub(r"^N_PERMUTATIONS = \d+", "N_PERMUTATIONS = 3", source, flags=re.M)
        try:
            exec(compile(source, f"{notebook.name}[cell {index}]", "exec"), namespace)
        except Exception as error:
            raise AssertionError(f"{notebook.name} code cell {index} failed:\n{source}") from error
        finally:
            plt.close("all")

    written = sorted(p.name for p in (synthetic_study / "processed" / "qc_confound").glob("*.csv"))
    assert written, "notebook wrote no result files"
    assert all(p.startswith(notebook.name[:2]) for p in written)
    for csv in (synthetic_study / "processed" / "qc_confound").glob(f"{notebook.name[:2]}_*.csv"):
        assert not pd.read_csv(csv).empty
