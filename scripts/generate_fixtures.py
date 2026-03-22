"""Generate minimal 10x Genomics fixture files for functional tests.

Creates a small but structurally valid dataset under data/fixtures/tenx/:

    data/fixtures/tenx/
        FIXTURE_features.tsv.gz        <- shared gene list (30 genes)
        raw/
            SAMPLEONE-barcodes.tsv.gz  <- 5 cells
            SAMPLEONE-matrix.mtx.gz
            SAMPLETWO-barcodes.tsv.gz  <- 5 cells
            SAMPLETWO-matrix.mtx.gz

The files are intentionally tiny (< 5 KB each) so they can be committed
to git and make tests run in milliseconds.

Run this script any time the fixture schema needs to change:

    python scripts/generate_fixtures.py
"""

import gzip
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = REPO_ROOT / "data" / "fixtures" / "tenx"

N_GENES = 30
N_CELLS = 5

GENE_IDS = [f"ENSG{i:011d}" for i in range(N_GENES)]
GENE_NAMES = [f"GENE{i:04d}" for i in range(N_GENES)]

SAMPLES = {
    "SAMPLEONE": [
        "AAACATACAACCAC-1",
        "AAACATTGAGCTAC-1",
        "AAACATTGATCAGC-1",
        "AAACCGTGCTTCCG-1",
        "AAACCGTGTATGCG-1",
    ],
    "SAMPLETWO": [
        "AAACATACAAGGGT-1",
        "AAACATACAGAGCC-1",
        "AAACATACAGCATA-1",
        "AAACATACAGCCTG-1",
        "AAACATACAGTAAG-1",
    ],
}


def _write_gz(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)


def _make_count_matrix(n_genes: int, n_cells: int, seed: int) -> np.ndarray:
    """Return a dense integer count matrix (genes × cells) with ~30 % non-zeros."""
    rng = np.random.default_rng(seed)
    counts = rng.integers(1, 20, size=(n_genes, n_cells))
    mask = rng.random(size=(n_genes, n_cells)) < 0.70
    counts[mask] = 0
    return counts


def _write_mtx(path: Path, matrix: np.ndarray) -> None:
    """Write a Matrix Market file (integer, coordinate format, 1-indexed) gzipped."""
    n_genes, n_cells = matrix.shape
    nonzero_rows, nonzero_cols = np.nonzero(matrix)
    nnz = len(nonzero_rows)

    lines = [
        "%%MatrixMarket matrix coordinate integer general",
        "%",
        f"{n_genes} {n_cells} {nnz}",
    ]
    for r, c in zip(nonzero_rows, nonzero_cols):
        lines.append(f"{r + 1} {c + 1} {matrix[r, c]}")

    _write_gz(path, "\n".join(lines) + "\n")


def write_features(dest: Path) -> None:
    """Write the shared features (gene list) file."""
    rows = [
        f"{gid}\t{gname}\tGene Expression"
        for gid, gname in zip(GENE_IDS, GENE_NAMES)
    ]
    _write_gz(dest, "\n".join(rows) + "\n")


def write_sample(raw_dir: Path, sample_id: str, barcodes: list[str], seed: int) -> None:
    """Write barcodes and matrix files for one sample."""
    _write_gz(
        raw_dir / f"{sample_id}-barcodes.tsv.gz",
        "\n".join(barcodes) + "\n",
    )
    matrix = _make_count_matrix(N_GENES, len(barcodes), seed=seed)
    _write_mtx(raw_dir / f"{sample_id}-matrix.mtx.gz", matrix)


def main() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    raw_dir = FIXTURE_DIR / "raw"
    raw_dir.mkdir(exist_ok=True)

    write_features(FIXTURE_DIR / "FIXTURE_features.tsv.gz")

    for seed, (sample_id, barcodes) in enumerate(SAMPLES.items()):
        write_sample(raw_dir, sample_id, barcodes, seed=seed * 17)

    print(f"Fixtures written to: {FIXTURE_DIR}")
    print(f"  Genes  : {N_GENES}")
    print(f"  Samples: {list(SAMPLES)}")
    print(f"  Cells per sample: {N_CELLS}")

    total_bytes = sum(f.stat().st_size for f in FIXTURE_DIR.rglob("*.gz"))
    print(f"  Total size: {total_bytes:,} bytes")


if __name__ == "__main__":
    main()
