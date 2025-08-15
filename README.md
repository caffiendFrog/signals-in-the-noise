<h1>
  <img src="images/project_logo_transparent.png" alt="Signals in the noise project logo" width="100" style="vertical-align: middle; margin-right: 0.5rem;">
  Signals in the Noise
</h1>

## Contents
* [Overview](#overview)
  * [Motivation](#motivation) 
  * [Objective](#objective)
* [Getting Started](#getting-started)
  * [Repository Structure](#repository-structure) 
  * [Running Jupyter Notebooks](#running-jupyter-notebooks)
* [References](#references)

## Introduction
Single-cell RNA sequencing (scRNA-seq) has played a pivotal role in advancing the understanding of biology by enabling researchers to measure gene expression at the resolution of individual cells. Through scRNA-seq analyses, researchers have created comprehensive cell atlases and identified rare and/or previously unrecognized cellular subpopulations. Unlike bulk RNA sequencing (bulk RNA-seq), which uses whole tissue or or bulk-sorted cells as inputs,  scRNA-seq further breaks down the tissue samples into individual cells as inputs [1].

A necessary byproduct of this level of resolution is a dramatic increase in the number of observations, often by 3 to 4 orders of magnitude, resulting in significantly more data for downstream analysis. Another challenge is that the process of tagging mRNA may incorrectly tag mRNA from multiple cells with the same barcode or fail to tag anything at all. These two challenges highlight the importance of verifying the quality of the reads and filtering out noise. This is commonly done by calculating metrics such as total number of genes, percentage of genes that are for mitochondria, and total number of barcodes (cells) that contain a gene. Thresholds are then determined for the dataset and cells that fall outside the threshold are filtered out from further analysis [2].

_[Back to Top](#contents)_

### Motivation
This project is motivated by the causal ambiguity of identifying thresholds for quality control (QC) metrics in the pre-processing workflow. Specifically, thresholds for scRNA-seq are set using biological assumptions, while those same or related assumptions are being evaluated by scRNA-seq. One such biological assumption is that cells with higher total RNA are metabolically healthy. As a result, the QC process often prioritizes these cells, while treating cells with low total RNA counts as technical artifacts to be filtered out [3, 4]. This approach, while effective for minimizing noise from ambient RNA contamination, risks eliminating biologically meaningful signals. 

| Feature to Threshold      | Filtered by QC Metric | Targeted by DDR | Dormant Cells                 |
|---------------------------|------------------------|-----------------|-------------------------------|
| Low total RNA content     | ✅ Damaged cell        | ⚠️ Depends      | ✅ Viable but quiet cell       |
| High total RNA content    | ✅ Degraded cell       | ⚠️ Depends      | ✅ Limited active gene expression |
| Low number of genes       | ✅ Technical artifact  | ⚠️ Depends      | ✅ Limited active gene expression |
| Low mitochondrial RNA %   | ❌ Not filtered out    | ✅ Damaged cell  | ✅ Limited energy needs        |
| High mitochondrial RNA %  | ✅ Damaged cell        | ✅ Damaged cell  | ❌ Not dormant                 |

*Table 1. Summary of QC metric thresholds and how they correspond to different kinds of cells.*

**Legend:** ✅ characteristic • ❌ not a characteristic • ⚠️ might be a characteristic (context-dependent)

_[Back to Top](#contents)_

### Objective
The goal of this study is to perform a comparative scRNA-seq analysis of cells classified as biological signals (“real”) versus those labeled as technical artifacts (“noise”), with the aim of evaluating whether current QC processes systematically exclude potentially informative cellular states.

This repository contains the framework used to perform comparative analysis.

_[Back to Top](#contents)_
## Getting Started
1. Prerequisites
   * Python 3.12 or higher
   * `pip` installed
2. Clone the repository
    ```bash
       git clone https://github.com/caffiendFrog/signals-in-the-noise.git
       cd signals-in-the-noise
    ```
3. Activate virtual environment for isolation
   * Windows (CMD or Powershell)
    ```bash
        python -m venv .venv
      .venv\Scripts\activate
    ```
  * macOS/Linux
    ```bash
        python3 -m venv .venv
        source .venv/bin/activate
    ```
  * PyCharm (_verified on PyCharm 2025.1.2, Windows 11 Home_)
    * Project Settings -> Python Interpreters -> Add Python Interpreter -> Local Interpreter
      * Select `.venv` (matching above activation environment)<img src="images/pycharm_screenshot.png" width="500" alt="screenshot of adding python interpreter to pycharm 2025.1.2"/>

4. Install runtime dependencies
    ```bash
        python .\bin\install_dependencies.py
    ```
    * NOTE: this has not been verified on a Mac, but should work if  

5. Install package in editable mode
    ```bash
        pip install -e .
    ```
   * _This will allow using the source code in the jupyter notebooks._

6. Download the datasets
    ```bash
        python .\bin\download_datasets.py
    ```
   * We will be using [GSE161529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529) [5, 6]. The datasets are prohibitively large to store in GitHub. Datasets can be downloaded directly from the Gene Expression Omnibus (GEO) or by using the provided python script which will download the files to the `assets` directory and expand the `tar` file for the patient samples. The authors provide the following statement about the data:
     * > “Processed scRNA‐seq and bulk RNA‐seq data generated for this study are available as GEO series GSE161529 and GSE161892, respectively. Raw data are available on request, subject to approval by our institutional Data Access Committee (dataaccess@wehi.edu.au) to ensure preservation of patient confidentiality.”
   * Using the python script will ensure compatibility with the rest of the downstream workflow (e.g. file naming conventions and locations).

_[Back to Top](#contents)_

### Repository Structure

* `bin`
  * Scripts to install dependencies and download datasets
* `data`
  * Data in various stages of preprocessing
  * Marked as `.gitignore` due to the size of data
  * Recommended to mark this directory as excluded from indexing in IDE
* `images`
  * Images used in documentation 
* `notebook`
  * Jupyter notebooks (the analysis)
* `resources`
  * Additonal resources for the data
  * Recommended to mark this directory as excluded from indexing in IDE.
* `src\signals_in_the-noise`
  * `preprocessing`
    * Source code for preprocessing data.
  * `utilities`
    * Source code for utility functions

_[Back to Top](#contents)_

### Running Jupyter Notebooks
Jupyter notebooks must be started within the virtual environment.
* Window (CMD or Powershell), macOS/Linux
  * After activating the virtualenv, you should see `(.venv)` prefixed to your command line.
  * Start jupyter notebooks you normally would.
* PyCharm 2025.1.2
  * If you have the pro version or on a trial that allows interacting directly with jupyter notebooks, be sure to select your activated virtual environment as the interpreter.

_[Back to Top](#contents)_

## References
1. Lafzi A, Moutinho C, Picelli S, Heyn H. Tutorial: guidelines for the experimental design of single-cell RNA sequencing studies. Nature protocols. London: Nature Publishing Group UK; 2018;13(12):2742–2757.
2. Luecken MD, Theis FJ. Current best practices in single‐cell RNA‐seq analysis: a tutorial. Molecular systems biology. London: Nature Publishing Group UK; 2019;15(6):e8746-n/a.
3. Young, Matthew D, Behjati, Sam. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. Gigascience. United States: Oxford University Press; 2020;9(12).
4. Cheng, Sophia K. Signals in the Noise: Uncovering the Biological Signatures of Ghost Cell Profiles in Human Breast Cancer. Dec 2025. Data Science for Social Good, University of Michigan, student paper. 
5. Yeh, Albert C, Ramaswamy, Sridhar. Mechanisms of Cancer Cell Dormancy--Another Hallmark of Cancer? Cancer research (Chicago, Ill). United States; 2015;75(23):5014–5022. 
6. Abad, Etna, Graifer, Dmitry, Lyakhovich, Alex. DNA damage response and resistance of cancer stem cells. Cancer letters. Ireland: Elsevier B.V; 2020;474:106–117.

_[Back to Top](#contents)_
