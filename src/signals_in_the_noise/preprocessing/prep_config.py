import json
from dataclasses import dataclass, asdict
from pathlib import Path

from signals_in_the_noise.utilities.logging import get_logger

L = get_logger(__name__)


@dataclass
class PrepConfig:
    """
    Configuration used to track steps performed during preprocessing
    """
    data_loaded: bool
    data_loaded_cached: bool
    annotations_added: bool
    annotations_added_cached: bool

    def to_json(self, path: Path, indent: int = 2) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as f:
            L.debug(f"Writing config to {path}.")
            json.dump(asdict(self), f, indent=indent, ensure_ascii=False)

    @classmethod
    def from_json(cls, path: Path) -> "PrepConfig":
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        return cls(**data)


class Prep:
    """
    Base class for preprocessing data
    """
    config: PrepConfig

    def __init__(self, config_filename: str):
        self.config_path = Path(config_filename)
        self.config_path.parent.mkdir(parents=True, exist_ok=True)

        if self.config_path.exists():
            self.config = PrepConfig.from_json(self.config_path)
        else:
            self.config = PrepConfig(False, False, False, False)

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
    def is_annotations_added(self) -> bool:
        return self.config.annotations_added

    @property
    def is_annotations_added_cached(self) -> bool:
        return self.config.annotations_added_cached

    def _save_config(self):
        self.config.to_json(self.config_path)