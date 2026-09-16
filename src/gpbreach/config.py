"""Run configuration loading.

A config is a YAML file in ``configs/``. Configs are a protected path: changing
one changes what a run means (CLAUDE.md Section 2.4).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


def repo_root() -> Path:
    """Repository root, i.e. the directory containing ``src/`` and ``configs/``."""
    return Path(__file__).resolve().parents[2]


def _resolve(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else (repo_root() / p).resolve()


@dataclass(frozen=True)
class Config:
    run_id: str
    breach: dict[str, Any]
    dam: dict[str, Any]
    components: dict[str, Any]
    parameters: dict[str, Any]
    sampling: dict[str, Any]
    time: dict[str, Any]
    path: Path | None = None
    raw: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, path: str | Path) -> Config:
        path = _resolve(path)
        with path.open() as fh:
            data = yaml.safe_load(fh)
        missing = {"run_id", "breach", "dam", "components", "parameters", "sampling", "time"} - set(data)
        if missing:
            raise ValueError(f"{path.name} is missing required sections: {sorted(missing)}")
        return cls(
            run_id=data["run_id"],
            breach=data["breach"],
            dam=data["dam"],
            components=data["components"],
            parameters=data["parameters"],
            sampling=data["sampling"],
            time=data["time"],
            path=path,
            raw=data,
        )

    @property
    def dam_path(self) -> Path:
        return _resolve(self.dam["processed"])

    def build_path(self, key: str) -> Path:
        return _resolve(self.dam["build"][key])

    def nominal(self, name: str) -> float:
        try:
            return float(self.parameters[name]["nominal"])
        except KeyError as exc:
            raise KeyError(f"parameter {name!r} has no 'nominal' value in {self.path}") from exc

    @property
    def years(self) -> list[int]:
        t = self.time
        return list(range(int(t["start_year"]), int(t["end_year"]) + 1, int(t.get("step_years", 1))))
