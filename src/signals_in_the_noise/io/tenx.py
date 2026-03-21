import re
import shutil
from collections import defaultdict
from enum import StrEnum
from pathlib import Path

import scanpy as sc

from signals_in_the_noise.config import get_data_path
from signals_in_the_noise.utils.log import get_logger

logger = get_logger(__name__)


class DirectoryType(StrEnum):
    MULTIPLE = "multiple"


class TenX:
    """Utility class for reconstituting 10x Genomics raw data into a sparse AnnData object.

    The directory for multiple samples (many sets of cells) is expected to contain:
        - <sample identifier>-barcodes.tsv.gz
        - <sample identifier>-matrix.mtx.gz

    With a single shared features file:
        - <study identifier>_features.tsv.gz
    """

    def __init__(
        self,
        directory: str,
        directory_type: DirectoryType,
        *,
        features_filename: str = None,
    ):
        """Initialise a TenX loader.

        Args:
            directory: Path to the directory containing raw 10x files.
            directory_type: Layout type of the directory.
            features_filename: Path to the features TSV file; required for MULTIPLE layout.

        Raises:
            ValueError: If required arguments are missing or the features filename format is wrong.
            FileNotFoundError: If ``features_filename`` does not exist on disk.
        """
        self.directory = Path(directory)
        self.directory_type = directory_type

        match self.directory_type:
            case DirectoryType.MULTIPLE:
                if not features_filename:
                    raise ValueError("features_filename required for directory type MULTIPLE")
                features_path = Path(features_filename)
                if not features_path.exists():
                    raise FileNotFoundError(f"Required file not found: {features_filename}")
                self.features_path = features_path

                if not self.features_path.name.endswith("_features.tsv.gz"):
                    raise ValueError(
                        "features_filename is not in the expected format, '_features.tsv.gz'"
                    )

                self.study_id = self.features_path.stem.replace("_features.tsv", "")
                self.multiple_adata = []

                self.study_directory = get_data_path(self.study_id)
                self.study_directory.mkdir(parents=True, exist_ok=True)
            case _:
                raise ValueError(f"Invalid directory_type {directory_type}")

    @property
    def cache_directory_name(self) -> Path:
        """Return the path of the h5ad cache directory for this study."""
        return get_data_path(f"{self.study_id}_adata_cache")

    def load_adata(self) -> None:
        """Load AnnData objects from the h5ad cache directory."""
        cache_directory = self.cache_directory_name
        for file in cache_directory.iterdir():
            logger.info(f"Reading {file} as AnnData object.")
            adata = sc.read_h5ad(file)
            adata.obs["adata-filename"] = file.name
            self.multiple_adata.append(adata)

    def load_data(self, *, cache: bool = True) -> None:
        """Load raw 10x data and optionally cache each sample as h5ad.

        Args:
            cache: When True, each loaded AnnData object is written to disk.
        """
        match self.directory_type:
            case DirectoryType.MULTIPLE:
                self.load_multiple_adata(cache=cache)

    def load_multiple_adata(self, *, cache: bool) -> None:
        """Reconstitute and load all samples in a multi-sample directory.

        Args:
            cache: When True, AnnData objects are written to h5ad files in the cache directory.
        """
        cache_directory = None
        if cache:
            cache_directory = self.cache_directory_name
            cache_directory.mkdir(parents=True, exist_ok=True)

        samples_to_files = self._samples_to_file_dictionary()
        self._reconstitute_ten_x_file_structure(samples_to_files, cache_directory)

    def _samples_to_file_dictionary(self) -> dict:
        """Build a mapping of sample IDs to their barcode and matrix filenames.

        Returns:
            Dictionary mapping sample identifier strings to lists of matching filenames.
        """
        pattern = re.compile(r"^(?P<sample_id>.+?)-(barcodes\.tsv|matrix\.mtx)\.gz$")
        samples_to_files: dict = defaultdict(list)
        for path in self.directory.iterdir():
            match = pattern.match(path.name)
            if not match:
                continue
            samples_to_files[match.group("sample_id")].append(path.name)
        return samples_to_files

    def _reconstitute_ten_x_file_structure(
        self, samples_to_files: dict, cache_directory: Path | None
    ) -> None:
        """Reorganise raw files into per-sample 10x layout and load as AnnData.

        Args:
            samples_to_files: Mapping of sample IDs to their source filenames.
            cache_directory: Directory to write h5ad cache files; None disables caching.
        """
        missing_targets: dict = defaultdict(list)
        skipped_files = []

        for sample_identifier, filenames in samples_to_files.items():
            sample_dir = self.study_directory / sample_identifier
            sample_dir.mkdir(parents=True, exist_ok=True)

            cached_path = (
                cache_directory / f"{sample_identifier}.h5ad" if cache_directory else None
            )
            if cached_path and cached_path.exists():
                logger.info(
                    f"Skipping {sample_identifier} because it already exists as an `.h5ad` file."
                )
                skipped_files.append(cached_path)
                continue

            for filename in filenames:
                target_path = None
                source_path = self.directory / filename
                if "barcodes" in filename:
                    target_path = sample_dir / "barcodes.tsv.gz"
                elif "matrix" in filename:
                    target_path = sample_dir / "matrix.mtx.gz"
                if target_path:
                    shutil.copy2(source_path, target_path)
                else:
                    missing_targets[sample_identifier].append(filename)

            features_target = sample_dir / "features.tsv.gz"
            shutil.copy2(self.features_path, features_target)

            if sample_identifier not in missing_targets:
                logger.info(f"Reading {sample_identifier} as AnnData object.")
                adata_filename = f"{sample_identifier}.h5ad"
                adata = sc.read_10x_mtx(path=str(sample_dir))
                adata.obs["adata-filename"] = adata_filename
                self.multiple_adata.append(adata)
                if cache_directory:
                    adata_path = cache_directory / adata_filename
                    logger.info("...caching object.")
                    adata.write_h5ad(adata_path)
            else:
                logger.warning(
                    f"Skipping {sample_identifier}, unable to determine target paths for "
                    f"{missing_targets[sample_identifier]}"
                )

        for file in skipped_files:
            logger.info(f"Loading cached object from file {file}")
            adata = sc.read_h5ad(file)
            adata.obs["adata-filename"] = file.name
            self.multiple_adata.append(adata)
