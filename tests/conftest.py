"""Shared pytest fixtures."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


FIXTURE_DIR = Path(__file__).parent.parent / "data" / "fixtures"


def _qc_block(n: int, rng: np.random.Generator, *, coupled: bool, shift: float) -> pd.DataFrame:
    """QC metrics for ``n`` cells.

    When ``coupled``, one latent factor drives low depth, low gene count and
    high mitochondrial fraction together, which concentrates cells in the
    pbs-1 corner; otherwise the three metrics are independent.
    """
    latent = rng.normal(size=n)
    noise = rng.normal(size=(3, n))
    if coupled:
        depth = latent + 0.3 * noise[0]
        genes = latent + 0.3 * noise[1]
        mito = -latent + 0.3 * noise[2]
    else:
        depth, genes, mito = noise
    log_total = 7.0 + 0.5 * depth + shift
    log_genes = 6.0 + 0.4 * genes + shift
    return pd.DataFrame(
        {
            "pct_counts_mt": 100 / (1 + np.exp(-(mito - 2.0))),
            "log1p_total_counts": log_total,
            "log1p_n_genes_by_counts": log_genes,
            "total_counts": np.expm1(log_total),
            "n_genes_by_counts": np.expm1(log_genes),
        }
    )


def make_specimen_obs(
    n_noise: int = 400,
    n_retained: int = 400,
    *,
    noise_coupled: bool = False,
    retained_coupled: bool = False,
    shift: float = 0.0,
    seed: int = 0,
) -> pd.DataFrame:
    """Synthetic cell-level QC frame for one specimen with an ``is_noise`` column."""
    rng = np.random.default_rng(seed)
    noise = _qc_block(n_noise, rng, coupled=noise_coupled, shift=shift).assign(is_noise=1)
    retained = _qc_block(n_retained, rng, coupled=retained_coupled, shift=shift).assign(is_noise=0)
    obs = pd.concat([noise, retained], ignore_index=True)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    return obs


def make_world(
    kind: str, n_per_condition: int = 6, n_cells: int = 400, seed: int = 0
) -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    """Synthetic ER+/Normal cohort.

    ``"specific"``: only ER+ discarded cells have coupled QC metrics.
    ``"confounded"``: every ER+ cell is shifted and coupled, as a sample-level
    technical difference (depth, dissociation) would produce.
    ``"null"``: no difference between conditions.
    """
    specimen_obs, conditions = {}, {}
    for condition in ("ER+ tumour", "Normal"):
        er = condition == "ER+ tumour"
        for i in range(n_per_condition):
            specimen_id = f"{'ER' if er else 'N'}-{i:02d}"
            kwargs = {
                "n_noise": n_cells,
                "n_retained": n_cells,
                "seed": seed * 1000 + len(specimen_obs),
            }
            if kind == "specific":
                kwargs["noise_coupled"] = er
            elif kind == "confounded":
                kwargs.update(noise_coupled=er, retained_coupled=er, shift=1.0 if er else 0.0)
            elif kind != "null":
                raise ValueError(kind)
            specimen_obs[specimen_id] = make_specimen_obs(**kwargs)
            conditions[specimen_id] = condition
    return specimen_obs, conditions


@pytest.fixture()
def specimen_obs_factory():
    """Return :func:`make_specimen_obs`."""
    return make_specimen_obs


@pytest.fixture()
def world_factory():
    """Return :func:`make_world`."""
    return make_world


@pytest.fixture()
def fixture_dir() -> Path:
    """Return the path to the small test-fixture data directory."""
    return FIXTURE_DIR


@pytest.fixture()
def preprocessor_config_dict() -> dict:
    """Return a minimal valid PreprocessorConfig as a plain dict."""
    return {
        "data_loaded": False,
        "annotations_loaded": False,
        "annotations_applied": False,
        "custom": [],
    }


@pytest.fixture()
def config_json_file(tmp_path, preprocessor_config_dict) -> Path:
    """Write a PreprocessorConfig JSON file to a temp directory and return its path."""
    config_file = tmp_path / "test_study.json"
    config_file.write_text(json.dumps(preprocessor_config_dict), encoding="utf-8")
    return config_file
