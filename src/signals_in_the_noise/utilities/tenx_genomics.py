import os
import re
import shutil
from collections import defaultdict
from enum import StrEnum


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
        self.directory = directory
        self.directory_type = directory_type

        match self.directory_type:
            case DirectoryType.MULTIPLE:
                if not features_filename:
                    raise ValueError("features_filename required for directory type MULTIPLE")
                if not os.path.exists(features_filename):
                    raise FileNotFoundError(f"Required file not found: {features_filename}")
                self.features_filename = features_filename
                self.multiple_adata = []

                # Regex to match study_identifier from the features file
                pattern = re.compile(r"^\.\./data/(?P<study_id>.+?)_features\.tsv\.gz$")
                match = pattern.match(self.features_filename)
                self.study_id = match.group("study_id")

                # Create target directory for re-organized files
                parent_dir = os.path.dirname(self.directory)
                self.study_directory = f"{parent_dir}/{self.study_id}"
                os.makedirs(self.study_directory, exist_ok=True)
            case _:
                raise ValueError(f"Invalid directory_type {directory_type}")

    def load_data(self):
        match self.directory_type:
            case DirectoryType.MULTIPLE:
                self.load_multiple_adata()

    def load_multiple_adata(self):
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
        for sample_identifier, filenames in samples_to_files.items():
            sample_dir = f"{self.study_directory}/{sample_identifier}"
            os.makedirs(sample_dir, exist_ok=True)
            for filename in filenames:
                source_path = f"{self.directory}/{filename}"
                if 'barcodes' in filename:
                    target_path = f"{sample_dir}/barcodes.tsv.gz"
                elif 'matrix' in filename:
                    target_path = f"{sample_dir}/matrix.mtx.gz"
                shutil.copy2(source_path, target_path)

            # Everyone gets the same features file
            source_path = f"{self.study_directory}_features.tsv.gz"
            target_path = f"{sample_dir}/features.tsv.gz"
            shutil.copy2(source_path, target_path)
