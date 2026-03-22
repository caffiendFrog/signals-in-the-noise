import re
import shutil
from collections import defaultdict
from pathlib import Path

import scanpy as sc

from signals_in_the_noise.config import get_data_path
from signals_in_the_noise.utils.log import get_logger

logger = get_logger(__name__)


class TenX:
    """Utility class for reconstituting 10x Genomics raw data into a sparse AnnData object.

    The directory is expected to contain one pair of files per sample:
        - <sample identifier>-barcodes.tsv.gz
        - <sample identifier>-matrix.mtx.gz

    With a single shared features file (named ``<study identifier>_features.tsv.gz``)
    supplied separately.
    """

    def __init__(self, directory: str, *, features_filename: str):
        """Initialize a TenX loader.

        Args:
            directory: Path to the directory containing the raw per-sample files.
            features_filename: Path to the shared features TSV file.
                Must end with ``_features.tsv.gz``.

        Raises:
            FileNotFoundError: If ``features_filename`` does not exist on disk.
            ValueError: If the features filename is not in the expected format.
        """
        self.directory = Path(directory)

        features_path = Path(features_filename)
        if not features_path.exists():
            raise FileNotFoundError(f"Required file not found: {features_filename}")
        if not features_path.name.endswith("_features.tsv.gz"):
            raise ValueError(
                "features_filename is not in the expected format, '_features.tsv.gz'"
            )

        self.features_path = features_path
        self.study_id = features_path.stem.replace("_features.tsv", "")
        self.multiple_adata = []

        self.study_directory = get_data_path(self.study_id)
        self.study_directory.mkdir(parents=True, exist_ok=True)

    @property
    def cache_directory_name(self) -> Path:
        """Return the path of the h5ad cache directory for this study."""
        return get_data_path(f"{self.study_id}_adata_cache")

    def load_adata(self) -> None:
        """Load AnnData objects from the h5ad cache directory.

        Does nothing if the cache directory does not yet exist.
        """
        cache_directory = self.cache_directory_name
        if not cache_directory.exists():
            logger.info(f"Cache directory does not exist, nothing to load: {cache_directory}")
            return
        for file in cache_directory.iterdir():
            logger.info(f"Reading {file} as AnnData object.")
            adata = sc.read_h5ad(file)
            adata.obs["adata-filename"] = file.name
            self.multiple_adata.append(adata)

    def load_data(self, *, cache: bool = True) -> None:
        """Load raw 10x data, reorganize into per-sample directories, and read as AnnData.

        Args:
            cache: When True, each loaded AnnData object is written to disk as h5ad.
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
        """Reorganize raw files into per-sample 10x layout and load as AnnData.

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
                source_path = self.directory / filename
                if "barcodes" in filename:
                    shutil.copy2(source_path, sample_dir / "barcodes.tsv.gz")
                elif "matrix" in filename:
                    shutil.copy2(source_path, sample_dir / "matrix.mtx.gz")
                else:
                    missing_targets[sample_identifier].append(filename)

            shutil.copy2(self.features_path, sample_dir / "features.tsv.gz")

            if sample_identifier not in missing_targets:
                logger.info(f"Reading {sample_identifier} as AnnData object.")
                adata_filename = f"{sample_identifier}.h5ad"
                adata = sc.read_10x_mtx(path=str(sample_dir))
                adata.obs["adata-filename"] = adata_filename
                self.multiple_adata.append(adata)
                if cache_directory:
                    logger.info("...caching object.")
                    adata.write_h5ad(cache_directory / adata_filename)
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
