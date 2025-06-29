import logging
import os
import re
import shutil
from collections import defaultdict
from enum import StrEnum
from pathlib import Path
import scanpy as sc

from signals_in_the_noise.utilities.storage import get_data_path

L = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[logging.StreamHandler(),]
)

class DirectoryType(StrEnum):
    # SINGLE = "single"
    MULTIPLE = "multiple"


class TenX:
    """
    Utility class for reconstituting 10x Genomics raw data into a sparse anndata object.

    The directory for a single sample (one set of cells) is expected to be:
        - <sample identifiers>-barcodes.tsv.gz
        - <sample identifiers>-features.tsv.gz (or -genes.tsv.gz)
        - <sample identifiers>-matrix.tsv.gz

    The directory for multiple samples (many sets of cells) is expected to be:
        - <sample identifiers>-barcodes.tsv.gz
        - <sample identifiers>-matrix.tsv.gz

    With a single corresponding features/genes file for the samples:
        - <study identifier>_features.tsv.gz
    """
    def __init__(self, directory: str, directory_type: DirectoryType, *, features_filename: str=None):
        """
        :param directory: directory with files to reconstitute.
        :param directory_type: type of directory.
        :param features_filename: filename of genomics features file, required if directory type is MULTIPLE.
        """
        self.directory = Path(directory)
        self.directory_type = directory_type

        match self.directory_type:
            case DirectoryType.MULTIPLE:
                if not features_filename:
                    raise ValueError("features_filename required for directory type MULTIPLE")
                if not os.path.exists(features_filename):
                    raise FileNotFoundError(f"Required file not found: {features_filename}")
                self.features_path = Path(features_filename)

                if not self.features_path.name.endswith("_features.tsv.gz"):
                    raise ValueError("features_filename is not in the expected format, '_features.tsv.gz'")

                self.study_id = self.features_path.stem.replace("_features.tsv", "")
                self.multiple_adata = []

                # Create target directory for re-organized files
                self.study_directory = get_data_path(self.study_id)
                self.study_directory.mkdir(parents=True, exist_ok=True)
            case _:
                raise ValueError(f"Invalid directory_type {directory_type}")

    def load_data(self, *, cache=True):
        """

        :param cache: True will save the adata to disk after loading
        :return:
        """
        match self.directory_type:
            case DirectoryType.MULTIPLE:
                self.load_multiple_adata(cache=cache)

    def load_multiple_adata(self, *, cache: bool):
        cache_directory = None
        if cache:
            cache_directory = get_data_path(f"{self.study_id}_adata_cache")
            cache_directory.mkdir(parents=True, exist_ok=True)
        # Regex to match sample_identifier from filenames
        pattern = re.compile(r"^(?P<sample_id>.+?)-(barcodes\.tsv|matrix\.mtx)\.gz$")

        # Build a dictionary of sample_identifier to list of barcodes & matrix files.
        samples_to_files = defaultdict(list)
        for filename in os.listdir(self.directory):
            match = pattern.match(filename)
            if not match:
                continue
            samples_to_files[match.group("sample_id")].append(filename)

        # Copy the files
        missing_targets = defaultdict(list)
        for sample_identifier, filenames in samples_to_files.items():
            sample_dir = f"{self.study_directory}/{sample_identifier}"
            os.makedirs(sample_dir, exist_ok=True)
            for filename in filenames:
                target_path = None
                source_path = f"{self.directory}/{filename}"
                if 'barcodes' in filename:
                    target_path = f"{sample_dir}/barcodes.tsv.gz"
                elif 'matrix' in filename:
                    target_path = f"{sample_dir}/matrix.mtx.gz"
                if target_path:
                    shutil.copy2(source_path, target_path)
                else:
                    missing_targets[sample_identifier].append(filename)

            # Everyone gets the same features file
            target_path = f"{sample_dir}/features.tsv.gz"
            shutil.copy2(self.features_path, target_path)

            # Read the 10x file as AnnData
            if sample_identifier not in missing_targets:
                L.info(f"Reading {sample_identifier} as AnnData object.")
                adata = sc.read_10x_mtx(path=sample_dir)
                self.multiple_adata.append(adata)
                if cache:
                    L.info(f"...caching object.")
                    adata.write_h5ad(f"{cache_directory}/{sample_identifier}.h5ad")
            else:
                L.warning(f"Skipping {sample_identifier}, unable to determine target paths for {missing_targets[sample_identifier]}")
