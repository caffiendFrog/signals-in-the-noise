<h1>
  <img src="images/project_logo_transparent.png" alt="Signals in the noise project logo" width="100" style="vertical-align: middle; margin-right: 0.5rem;">
  Signals in the Noise
</h1>

The graded final report for this project can be found here: https://github.com/caffiendFrog/portfolio/blob/main/academics/SIADS-699/.

## Contents

* [Introduction](#introduction)
  * [Motivation](#motivation)
  * [Objective](#objective)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [1 — Clone the repository](#1--clone-the-repository)
  * [2 — Create the conda environment](#2--create-the-conda-environment)
  * [3 — Install the package](#3--install-the-package)
  * [4 — Download the dataset](#4--download-the-dataset)
  * [5 — Verify the setup](#5--verify-the-setup)
  * [6 — Run the notebooks](#6--run-the-notebooks)
* [Repository Structure](#repository-structure)
* [Configuration](#configuration)
* [References](#references)

---

## Introduction

Single-cell RNA sequencing (scRNA-seq) has played a pivotal role in advancing the understanding of biology by enabling researchers to measure gene expression at the resolution of individual cells. Through scRNA-seq analyses, researchers have created comprehensive cell atlases and identified rare and/or previously unrecognized cellular subpopulations. Unlike bulk RNA sequencing (bulk RNA-seq), which uses whole tissue or bulk-sorted cells as inputs, scRNA-seq further breaks down the tissue samples into individual cells as inputs [1].

A necessary byproduct of this level of resolution is a dramatic increase in the number of observations, often by 3 to 4 orders of magnitude, resulting in significantly more data for downstream analysis. Another challenge is that the process of tagging mRNA may incorrectly tag mRNA from multiple cells with the same barcode or fail to tag anything at all. These two challenges highlight the importance of verifying the quality of the reads and filtering out noise. This is commonly done by calculating metrics such as total number of genes, percentage of genes that are for mitochondria, and total number of barcodes (cells) that contain a gene. Thresholds are then determined for the dataset and cells that fall outside the threshold are filtered out from further analysis [2].

_[Back to Top](#contents)_

### Motivation

This project is motivated by the causal ambiguity of identifying thresholds for quality control (QC) metrics in the pre-processing workflow. Specifically, thresholds for scRNA-seq are set using biological assumptions, while those same or related assumptions are being evaluated by scRNA-seq. One such biological assumption is that cells with higher total RNA are metabolically healthy. As a result, the QC process often prioritizes these cells, while treating cells with low total RNA counts as technical artifacts to be filtered out [3, 4]. This approach, while effective for minimizing noise from ambient RNA contamination, risks eliminating biologically meaningful signals.

| Feature to Threshold      | Filtered by QC Metric | Targeted by DDR | pbs-2                         |
|---------------------------|------------------------|-----------------|-------------------------------|
| Low total RNA content     | ✅ pbs-1               | ⚠️ Depends      | ✅ Viable but quiet cell       |
| High total RNA content    | ✅ Degraded cell       | ⚠️ Depends      | ✅ Limited active gene expression |
| Low number of genes       | ✅ Technical artifact  | ⚠️ Depends      | ✅ Limited active gene expression |
| Low mitochondrial RNA %   | ❌ Not filtered out    | ✅ pbs-1         | ✅ Limited energy needs        |
| High mitochondrial RNA %  | ✅ pbs-1               | ✅ pbs-1         | ❌ Not pbs-2                   |

*Table 1. Summary of QC metric thresholds and how they correspond to different kinds of cells.*

**Legend:** ✅ characteristic • ❌ not a characteristic • ⚠️ might be a characteristic (context-dependent)

_[Back to Top](#contents)_

### Objective

The goal of this study is to perform a comparative scRNA-seq analysis of cells classified as biological signals ("real") versus those labeled as technical artifacts ("noise"), with the aim of evaluating whether current QC processes systematically exclude potentially informative cellular states.

This repository contains the framework used to perform the comparative analysis.

_[Back to Top](#contents)_

---

## Getting Started

### Prerequisites

| Requirement | Version |
|---|---|
| Operating system | Windows 10 / 11 |
| [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda | Any recent version |
| Git | Any recent version |

> **Note:** The environment and all path handling are Windows-specific. The conda
> environment definition pins every dependency to an exact version for reproducibility.

_[Back to Top](#contents)_

### 1 — Clone the repository

```powershell
git clone https://github.com/caffiendFrog/signals-in-the-noise.git
cd signals-in-the-noise
```

_[Back to Top](#contents)_

### 2 — Create the conda environment

All dependencies, including Python itself, are declared in `environment.yml`.
Create the environment once and activate it for every subsequent session.

```powershell
conda env create -f environment.yml
conda activate signals-in-the-noise
```

To update the environment after a dependency change (e.g. pulling new commits):

```powershell
conda env update -f environment.yml --prune
```

_[Back to Top](#contents)_

### 3 — Install the package

Install the `signals_in_the_noise` package in editable mode so that the source
code in `src/` is importable from the notebooks without any path manipulation:

```powershell
pip install -e .
```

_[Back to Top](#contents)_

### 4 — Download the dataset

This project uses [GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529) [5, 6] — a human breast cancer scRNA-seq dataset. The raw files are too large to store in the repository and must be downloaded separately.

> "Processed scRNA-seq and bulk RNA-seq data generated for this study are available as GEO series GSE161529 and GSE161892, respectively. Raw data are available on request, subject to approval by our institutional Data Access Committee (dataaccess@wehi.edu.au) to ensure preservation of patient confidentiality."

The download script fetches the tar archive and the shared features file directly from GEO, then extracts the archive into `data/raw/`:

```powershell
python scripts\download_data.py
```

After the script completes, `data/raw/` will contain:

```
data/raw/
    GSE161529_RAW/               ← per-sample barcode and matrix files
    GSE161529_features.tsv.gz    ← shared gene features file
```

You will also need to place the supplementary Excel files provided by the authors into `data/reference/`:

```
data/reference/
    GSE161529/
        table_supplementary_1.xlsx
        table_supplementary_2.xlsx
        table_ev_4.xlsx
```

_[Back to Top](#contents)_

### 5 — Verify the setup

Run the test suite to confirm the environment and package are working correctly.
The tests use small committed fixture files and do not require the full dataset.

```powershell
pytest
```

All tests should pass. Expected output:

```
72 passed in ~10s
```

_[Back to Top](#contents)_

### 6 — Run the notebooks

Start Jupyter from within the activated conda environment:

```powershell
jupyter lab
```

The analysis is organised into numbered notebooks under `notebooks/GSE161529/`.
Run them in order:

| Notebook | Description |
|---|---|
| `00-a-epi-cell-typing-reverse-engineer.ipynb` | Reverse-engineer the authors' epithelial cell typing approach |
| `01-epi-cell-typing-figure-1c.ipynb` | Reproduce Figure 1C |
| `02-epi-cell-typing-figure-1d.ipynb` | Reproduce Figure 1D |
| `03-epi-cell-typing-figure-1e.ipynb` | Reproduce Figure 1E |
| `04-epi-cell-typing-figure-1h.ipynb` | Reproduce Figure 1H |
| `05-epi-cell-typing-empty-cells.ipynb` | Analyze cells labeled as empty / noise |
| `06-a-epi-cell-typing-associated-pathways-raw.ipynb` | Pathway analysis — raw expression |
| `06-b-epi-cell-typing-associated-pathways-avg.ipynb` | Pathway analysis — averaged expression |
| `06-c-epi-cell-typing-associated-pathways-gsea.ipynb` | Gene set enrichment analysis |
| `07-epi-cell-typing-total-populations.ipynb` | Compare total cell populations |
| `08-findings-eda.ipynb` | Exploratory data analysis of findings |

Each notebook imports from the `signals_in_the_noise` package installed in step 3.
To enable logging output in a notebook, add this to the first cell:

```python
from signals_in_the_noise.utils.logging_config import setup_logging
setup_logging()
```

_[Back to Top](#contents)_

---

## Repository Structure

```
signals-in-the-noise/
├── src/
│   └── signals_in_the_noise/
│       ├── config.py                   # project-wide path constants and helper functions
│       ├── analysis/
│       │   ├── noise_phenotypes.py     # noise-cell phenotype annotation logic
│       │   └── statistics.py          # statistical comparison helpers
│       ├── io/
│       │   ├── gmt.py                  # GMT gene-set file parser
│       │   └── tenx.py                 # 10x Genomics file reconstitution and AnnData loading
│       ├── preprocessing/
│       │   ├── base.py                 # Preprocessor base class and PreprocessorConfig dataclass
│       │   └── gse161529.py            # GSE161529-specific preprocessing subclass
│       └── utils/
│           ├── log.py                  # module-level logger factory (no handler configuration)
│           ├── logging_config.py       # setup_logging() for entry points and notebooks
│           └── visualization.py        # matplotlib figure/axes grid helper
├── tests/
│   ├── conftest.py                     # shared pytest fixtures
│   ├── analysis/
│   │   ├── test_noise_phenotypes.py
│   │   └── test_statistics.py
│   ├── functional/
│   │   └── test_tenx_functional.py    # end-to-end tests using committed fixture files
│   ├── io/
│   │   ├── test_gmt.py
│   │   └── test_tenx.py
│   ├── preprocessing/
│   │   ├── test_base.py
│   │   └── test_gse161529.py
│   ├── utils/
│   │   ├── test_log.py
│   │   ├── test_logging_config.py
│   │   └── test_visualization.py
│   └── test_config.py
├── notebooks/
│   └── GSE161529/                      # numbered analysis notebooks (run in order)
├── scripts/
│   ├── download_data.py                # downloads and extracts GSE161529 raw data
│   └── generate_fixtures.py            # regenerates committed test fixture files
├── data/
│   ├── raw/                            # git-ignored; downloaded dataset files
│   ├── processed/                      # git-ignored; intermediate pipeline outputs
│   ├── reference/                      # git-ignored; supplementary Excel files
│   └── fixtures/
│       └── tenx/                       # small committed files used by functional tests
├── resources/
│   └── GSE161529/
│       ├── epithelial_cell_typing/     # gene signature Excel files (Lim et al. 2009)
│       ├── gene_set_enrichment_analysis/  # MSigDB Hallmark GMT files
│       └── table_*.xlsx                # supplementary tables committed to git
├── environment.yml                     # conda environment — single source of truth for dependencies
├── pyproject.toml                      # build metadata and pytest/ruff configuration
└── README.md
```

_[Back to Top](#contents)_

---

## Configuration

All project-wide paths are resolved in `src/signals_in_the_noise/config.py`. It exposes two constants and two helper functions:

| Name | Type | Description |
|---|---|---|
| `PROJECT_ROOT` | `Path` | Repository root (resolved at import time via `__file__`) |
| `DATA_DIRECTORY` | `Path` | `<PROJECT_ROOT>/data/` |
| `RESOURCES_DIRECTORY` | `Path` | `<PROJECT_ROOT>/resources/` |
| `get_data_path(subpath)` | function | Returns `DATA_DIRECTORY / subpath`, or `DATA_DIRECTORY` if `subpath` is `None` |
| `get_resources_path(subpath)` | function | Returns `RESOURCES_DIRECTORY / subpath`, or `RESOURCES_DIRECTORY` if `subpath` is `None` |

There are no hardcoded absolute paths anywhere in the source code. Every file I/O operation resolves its path through one of these helpers.

_[Back to Top](#contents)_

---

## References

1. Lafzi A, Moutinho C, Picelli S, Heyn H. Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies. Nature protocols. London: Nature Publishing Group UK; 2018;13(12):2742–2757.
2. Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Molecular systems biology. London: Nature Publishing Group UK; 2019;15(6):e8746-n/a.
3. Young, Matthew D, Behjati, Sam. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. Gigascience. United States: Oxford University Press; 2020;9(12).
4. Cheng, Sophia K. Signals in the Noise: Uncovering the Biological Signatures of Ghost Cell Profiles in Human Breast Cancer. Dec 2025. Data Science for Social Good, University of Michigan, student paper.
5. Yeh, Albert C, Ramaswamy, Sridhar. Mechanisms of Cancer Cell Dormancy--Another Hallmark of Cancer? Cancer research (Chicago, Ill). United States; 2015;75(23):5014–5022.
6. Abad, Etna, Graifer, Dmitry, Lyakhovich, Alex. DNA damage response and resistance of cancer stem cells. Cancer letters. Ireland: Elsevier B.V; 2020;474:106–117.

_[Back to Top](#contents)_
