#!/usr/bin/env python

import os

import subprocess
import sys
from pathlib import Path
import tarfile
from urllib import request


project_root = Path(__file__).resolve().parent.parent
assets_directory = project_root / "assets"
raw_data_file_path =  os.path.join(assets_directory, "GSE161529_RAW.tar")
features_file_path = os.path.join(assets_directory, "GSE161529_features.tsv.gz")

raw_data_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161529&format=file"
features_url = (
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE161529&format=file&file=GSE161529%5Ffeatures%2Etsv%2Egz"
)

for url, file_path, extract_to in (
    (raw_data_url, raw_data_file_path, "GSE161529_RAW"),
    (features_url, features_file_path, None),
):
    if not os.path.exists(file_path):
        print(f"Downloading raw data from {url} to {file_path}...")
        request.urlretrieve(url, file_path)
        if extract_to:
            extract_to_path = os.path.join(assets_directory, extract_to)
            if not os.path.exists(extract_to_path):
                mode = "r:gz" if file_path.endswith(".gz") else "r:"
                print(f"Extracting {file_path} to {extract_to_path} with {mode}...")
                with tarfile.open(file_path, mode) as tar:
                    tar.extractall(extract_to_path)
        print(f"Done!")
    else:
        print(f"File {file_path} already exists.")

