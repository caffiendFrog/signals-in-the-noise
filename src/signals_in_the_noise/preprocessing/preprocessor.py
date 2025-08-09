import json
import random

from collections import defaultdict
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List

import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData

from signals_in_the_noise.utilities.logging import get_logger
from signals_in_the_noise.utilities.storage import get_data_path
from signals_in_the_noise.utilities.storage import get_resources_path

L = get_logger(__name__)


@dataclass
class PreprocessorConfig:
    """
    Configuration used to track steps performed during preprocessing

    """
    data_loaded: bool
    annotations_loaded: bool
    annotations_applied: bool
    custom: List[str] = field(default_factory=list)

    def to_json(self, path: Path, indent: int = 2) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as f:
            L.debug(f"Writing config to {path}.")
            json.dump(asdict(self), f, indent=indent, ensure_ascii=False)

    @classmethod
    def from_json(cls, path: Path) -> "PreprocessorConfig":
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        return cls(**data)


class Preprocessor:
    """
    Base class for preprocessing data
    """
    def __init__(self, study_id: str):
        self.cache_directory_path = ''
        self.objects = defaultdict(AnnData)
        self.random_kwargs = defaultdict(str)

        self.STUDY_ID = study_id
        self.config_path = Path(get_data_path(f"{study_id}.json"))
        self.config_path.parent.mkdir(parents=True, exist_ok=True)

        if self.config_path.exists():
            self.config = PreprocessorConfig.from_json(self.config_path)
        else:
            self.config = PreprocessorConfig(False, False, False)

    @property
    def is_data_loaded(self) -> bool:
        return self.config.data_loaded

    def data_loaded(self) -> None:
        self.config.data_loaded = True
        self._save_config()

    @property
    def is_annotations_loaded(self) -> bool:
        return self.config.annotations_loaded


    def annotations_loaded(self) -> None:
        self.config.annotations_loaded = True
        self._save_config()

    @property
    def is_annotations_applied(self) -> bool:
        return self.config.annotations_applied

    def annotations_applied(self) -> None:
        self.config.annotations_applied = True
        self._save_config()

    def has_custom(self, value: str) -> bool:
        return value in self.config.custom

    def add_custom(self, value: str) -> None:
        self.config.custom.append(value)
        self._save_config()

    def _save_config(self):
        self.config.to_json(self.config_path)

    def cache_adata_object(self, adata: AnnData, filename: str):
        if self.cache_directory_path:
            adata.write(self.cache_directory_path / filename)

    def get_dataset(self, filename):
        """
        Returns a copy of the dataset
        :param filename:
        :return: a copy of the dataset if it exists, otherwise an empty AnnData object.
        """
        actual = self.objects.get(filename, None)
        if actual:
            return actual.copy()
        return AnnData()

    @property
    def random_seed(self):
        if 'random_state' not in self.random_kwargs:
            self.random_kwargs['random_state'] = random.randint(1, 100)
        return self.random_kwargs.get('random_state', None)

    def cache_raw_gene_expression(self, adata, genes_of_interest, *, in_place=False):
        """
        Store the raw gene expression for selected genes before downstream analysis and
        dimensionality reduction filters them out of the dataset.

        :param adata: AnnData
        :param genes_of_interest: list of genes
        :param in_place: True to modify in place (e.g. don't make a copy of the given adata)
        :return: AnnData
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

    def check_adata_for_genes(self, adata, genes_to_check):
        """
        Utility method to check if the given genes exist in the data set.

        Check is case insensitive.

        :param adata: AnnData
        :param genes_to_check: list of strings that are gene names to check.
        :return:
        """
        var_names_lower = {name.lower(): name for name in adata.var_names}
        missing = []
        for gene in genes_to_check:
            if gene.lower() not in var_names_lower:
                missing.append(gene)

        if missing:
            L.info(f"{len(missing)} missing out of {len(genes_to_check)}")
        return missing

    def score_gene_signature_expression(
            self,
            adata: AnnData,
            gene_signature_filenames: Dict[str, str],
            *,
            log_normalize: bool,
            hvg_only: bool,
            hvg_flavor: str,
    ):
        """
        Scores the given dataset for the gene signatures.

        Currently only supports signature files in excel and formatted for gse161529

        :param adata: AnnData object to score for gene expression
        :param gene_signature_filenames: gene name to its expression signature filename
        :param in_place: True to modify adata in place
        :param log_normalize: True if adata needs to be log-normalized
        :param hvg_only: True if adata needs to be filtered to highly variable genes
        :param hvg_flavor: One of the flavors for sc.pp.highly_variable_genes
        :return: Dataset that now has a new column for each gene signature containing the average expression of the gene for each cell.
        """
        if log_normalize:
            print("Log normalizing dataset...")
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
        if hvg_only:
            print("Filtering dataset to highly variable genes...")
            sc.pp.highly_variable_genes(adata, flavor=hvg_flavor)
            # save a static snapshot of the data before filtering out
            adata.raw = adata.copy()
            adata = adata[:, adata.var['highly_variable']].copy()

        for gene_signature, filename in gene_signature_filenames.items():
            # read gene signature from excel file
            signature_df = pd.read_excel(get_resources_path(self.STUDY_ID + '/' + filename))
            genes = signature_df.loc[:, 'Symbol'].dropna().unique().tolist()
            # only return the genes that exist in target dataset
            actual_genes = [gene for gene in genes if gene in adata.var_names]
            sc.tl.score_genes(adata, gene_list=actual_genes, score_name=f'score_{gene_signature}', random_state=43)

        return adata
