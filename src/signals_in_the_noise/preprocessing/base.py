import json
import random
from collections import defaultdict
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from signals_in_the_noise.config import get_data_path, get_resources_path
from signals_in_the_noise.utils.log import get_logger

logger = get_logger(__name__)


@dataclass
class PreprocessorConfig:
    """Track which preprocessing steps have been completed and persisted.

    Attributes:
        data_loaded: True once raw data has been loaded.
        annotations_loaded: True once annotation metadata has been merged in.
        annotations_applied: True once derived annotation columns have been written.
        custom: List of arbitrary step names recorded by subclasses.
    """

    data_loaded: bool
    annotations_loaded: bool
    annotations_applied: bool
    custom: List[str] = field(default_factory=list)

    def to_json(self, path: Path, indent: int = 2) -> None:
        """Serialise this config to a JSON file.

        Args:
            path: Destination file path (parent directories are created if needed).
            indent: JSON indentation level.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as f:
            logger.debug(f"Writing config to {path}.")
            json.dump(asdict(self), f, indent=indent, ensure_ascii=False)

    @classmethod
    def from_json(cls, path: Path) -> "PreprocessorConfig":
        """Deserialise a PreprocessorConfig from a JSON file.

        Args:
            path: Path to the JSON file.

        Returns:
            A new PreprocessorConfig instance populated from the file.
        """
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        return cls(**data)


class Preprocessor:
    """Base class for preprocessing single-cell RNA-seq data."""

    def __init__(self, study_id: str):
        """Initialise the preprocessor and load or create its config.

        Args:
            study_id: Identifier for the study, used to derive file paths.
        """
        self.cache_directory_path: Path | str = ""
        self.objects: dict = defaultdict(AnnData)
        self.random_kwargs: dict = defaultdict(str)

        self.STUDY_ID = study_id
        self.config_path = Path(get_data_path(f"{study_id}.json"))
        self.config_path.parent.mkdir(parents=True, exist_ok=True)

        if self.config_path.exists():
            self.config = PreprocessorConfig.from_json(self.config_path)
        else:
            self.config = PreprocessorConfig(False, False, False)

    @property
    def is_data_loaded(self) -> bool:
        """True if raw data has been loaded and recorded in the config."""
        return self.config.data_loaded

    def data_loaded(self) -> None:
        """Mark raw data as loaded and persist the config."""
        self.config.data_loaded = True
        self._save_config()

    @property
    def is_annotations_loaded(self) -> bool:
        """True if annotation metadata has been merged and recorded."""
        return self.config.annotations_loaded

    def annotations_loaded(self) -> None:
        """Mark annotations as loaded and persist the config."""
        self.config.annotations_loaded = True
        self._save_config()

    @property
    def is_annotations_applied(self) -> bool:
        """True if derived annotation columns have been applied."""
        return self.config.annotations_applied

    def annotations_applied(self) -> None:
        """Mark annotations as applied and persist the config."""
        self.config.annotations_applied = True
        self._save_config()

    def has_custom(self, value: str) -> bool:
        """Return True if ``value`` has been recorded as a custom step.

        Args:
            value: Step name to check.
        """
        return value in self.config.custom

    def add_custom(self, value: str) -> None:
        """Record a custom step name and persist the config.

        Args:
            value: Step name to record.
        """
        self.config.custom.append(value)
        self._save_config()

    def _save_config(self) -> None:
        self.config.to_json(self.config_path)

    def cache_adata_object(self, adata: AnnData, filename: str) -> None:
        """Write an AnnData object to the cache directory if one is set.

        Args:
            adata: AnnData object to write.
            filename: Filename (relative to cache directory) to write to.
        """
        if self.cache_directory_path:
            adata.write(self.cache_directory_path / filename)

    def get_dataset(self, filename: str) -> AnnData:
        """Return a copy of the stored dataset for a given filename.

        Args:
            filename: Key under which the AnnData was stored.

        Returns:
            A copy of the dataset if found, otherwise an empty AnnData.
        """
        actual = self.objects.get(filename, None)
        if actual is not None:
            return actual.copy()
        return AnnData()

    @property
    def random_seed(self) -> int:
        """Return the random seed, generating one on first access."""
        if "random_state" not in self.random_kwargs:
            self.random_kwargs["random_state"] = random.randint(1, 100)
        return self.random_kwargs.get("random_state", None)

    def cache_raw_gene_expression(
        self, adata: AnnData, genes_of_interest: dict, *, in_place: bool = False
    ) -> AnnData:
        """Store raw gene expression for selected genes as observation columns.

        Captures expression before downstream filtering removes genes from the dataset.

        Args:
            adata: Input AnnData object.
            genes_of_interest: Mapping of obs column name to gene name.
            in_place: When True, modifies ``adata`` directly; otherwise operates on a copy.

        Returns:
            AnnData with new obs columns for each gene of interest.
        """
        var_names_lower = {name.lower(): name for name in adata.var_names}
        if not in_place:
            adata = adata.copy()

        for obs_name, gene_name in genes_of_interest.items():
            if gene_name.lower() in var_names_lower:
                gene_expression = adata[:, gene_name].X
                if not isinstance(gene_expression, np.ndarray):
                    gene_expression = gene_expression.toarray()
                gene_expression = gene_expression.flatten()
            else:
                gene_expression = np.zeros(adata.n_obs)

            adata.obs[obs_name] = gene_expression

        return adata

    def check_adata_for_genes(self, adata: AnnData, genes_to_check: list) -> list:
        """Return the subset of gene names absent from the dataset (case-insensitive).

        Args:
            adata: AnnData object to check.
            genes_to_check: Gene names to look up.

        Returns:
            List of gene names that were not found in ``adata.var_names``.
        """
        var_names_lower = {name.lower(): name for name in adata.var_names}
        missing = [gene for gene in genes_to_check if gene.lower() not in var_names_lower]
        if missing:
            logger.info(f"{len(missing)} missing out of {len(genes_to_check)}")
        return missing

    def score_gene_signature_expression(
        self,
        adata: AnnData,
        gene_signature_filenames: Dict[str, str],
        *,
        log_normalize: bool,
        hvg_only: bool,
        hvg_flavor: str,
    ) -> AnnData:
        """Score the dataset for gene signature expression.

        Reads each signature from an Excel file, separates up- and down-regulated
        genes, scores each set independently, and stores the difference as a new
        obs column named ``score_<signature>``.

        Args:
            adata: AnnData object to score.
            gene_signature_filenames: Mapping of signature name to Excel filename
                (relative to the study resources directory).
            log_normalize: When True, apply total-count normalisation and log1p.
            hvg_only: When True, filter to highly variable genes before scoring.
            hvg_flavor: Flavor string for ``sc.pp.highly_variable_genes``.

        Returns:
            AnnData with a new obs column per gene signature.
        """
        if log_normalize:
            logger.info("Log normalizing dataset...")
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
        if hvg_only:
            logger.info("Filtering dataset to highly variable genes...")
            sc.pp.highly_variable_genes(adata, flavor=hvg_flavor)
            adata.raw = adata.copy()
            adata = adata[:, adata.var["highly_variable"]].copy()

        for gene_signature, filename in gene_signature_filenames.items():
            signature_df = pd.read_excel(get_resources_path(self.STUDY_ID + "/" + filename))
            signature_df = signature_df.loc[:, ["Symbol", "Average log fold-change"]].dropna()
            upregulated_genes = (
                signature_df.loc[signature_df["Average log fold-change"] >= 0, "Symbol"]
                .unique()
                .tolist()
            )
            downregulated_genes = (
                signature_df.loc[signature_df["Average log fold-change"] < 0, "Symbol"]
                .unique()
                .tolist()
            )
            actual_up_genes = [gene for gene in upregulated_genes if gene in adata.var_names]
            actual_down_genes = [gene for gene in downregulated_genes if gene in adata.var_names]
            sc.tl.score_genes(
                adata, gene_list=actual_up_genes, score_name="up_score",
                random_state=self.random_seed,
            )
            sc.tl.score_genes(
                adata, gene_list=actual_down_genes, score_name="down_score",
                random_state=self.random_seed,
            )
            adata.obs[f"score_{gene_signature}"] = (
                adata.obs["up_score"] - adata.obs["down_score"]
            )
            adata.obs.drop(columns=["up_score", "down_score"], inplace=True)

        return adata
