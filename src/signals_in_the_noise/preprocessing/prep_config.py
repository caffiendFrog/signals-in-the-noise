import json
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path

from anndata import AnnData

from signals_in_the_noise.utilities.logging import get_logger
from signals_in_the_noise.utilities.storage import get_data_path

L = get_logger(__name__)


@dataclass
class PreprocessorConfig:
    """
    Configuration used to track steps performed during preprocessing

    Might not need "cached" options, but including for now just in case.
    """
    data_loaded: bool
    data_loaded_cached: bool
    annotations_loaded: bool
    annotations_loaded_cached: bool
    annotations_applied: bool
    annotations_applied_cached: bool

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
    config: PreprocessorConfig
    cache_directory_path: Path
    objects: defaultdict

    def __init__(self, study_id: str):
        self.STUDY_ID = study_id
        self.config_path = Path(get_data_path(f"{study_id}.json"))
        self.config_path.parent.mkdir(parents=True, exist_ok=True)

        if self.config_path.exists():
            self.config = PreprocessorConfig.from_json(self.config_path)
        else:
            self.config = PreprocessorConfig(False, False, False, False, False, False)

    @property
    def is_data_loaded(self) -> bool:
        return self.config.data_loaded

    def data_loaded(self) -> None:
        self.config.data_loaded = True
        self.config.data_loaded_cached = True
        self._save_config()

    @property
    def is_data_loaded_cached(self) -> bool:
        return self.config.data_loaded_cached

    @property
    def is_annotations_loaded(self) -> bool:
        return self.config.annotations_loaded

    @property
    def is_annotations_loaded_cached(self) -> bool:
        return self.config.annotations_loaded_cached

    def annotations_loaded(self) -> None:
        self.config.annotations_loaded = True
        self.config.annotations_loaded_cached = True
        self._save_config()

    @property
    def is_annotations_applied(self) -> bool:
        return self.config.annotations_applied

    @property
    def is_annotations_applied_cached(self) -> bool:
        return self.config.annotations_applied_cached

    def annotations_applied(self) -> None:
        self.config.annotations_applied = True
        self.config.annotations_applied_cached = True
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
