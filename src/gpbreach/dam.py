"""The ``DamObject``: the small, engine-facing summary of one breach site.

Per CLAUDE.md Section 3, all spatial work happens in ``preprocess/``. The engine
never opens a raster; it only sees a DamObject.

For this site the "collapse space early" rule does *not* reduce to a 1D profile:
the owner's breach criterion is a connectivity test across the 2D dam (see
docs/baseline_method.md). The DamObject therefore carries two co-registered
grids plus the two seed cells whose connection defines the breach.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class DamObject:
    """Static inputs for one breach site.

    Attributes
    ----------
    name:
        Identifier, e.g. ``"breach1_lgp"``.
    surface:
        Glacier surface elevation at ``base_year`` (m), NaN outside the domain.
    bed:
        Bed/bedrock elevation (m), NaN outside the domain.
    seed_lake, seed_target:
        ``(row, col)`` cells, one in each of the two basins whose connection
        constitutes the breach. See D-005.
    base_year:
        Calendar year of ``surface``.
    transform, crs, vertical_datum:
        Geo-referencing and the confirmed vertical datum of both elevation grids.
    """

    name: str
    surface: np.ndarray
    bed: np.ndarray
    seed_lake: tuple[int, int]
    seed_target: tuple[int, int]
    base_year: int
    transform: tuple[float, float, float, float, float, float]
    crs: str | None = None
    vertical_datum: str | None = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.surface.shape != self.bed.shape:
            raise ValueError(
                f"surface {self.surface.shape} and bed {self.bed.shape} are not co-registered"
            )
        for label, seed in (("seed_lake", self.seed_lake), ("seed_target", self.seed_target)):
            r, c = seed
            if not (0 <= r < self.surface.shape[0] and 0 <= c < self.surface.shape[1]):
                raise ValueError(f"{label} {seed} is outside the grid {self.surface.shape}")

    @property
    def valid(self) -> np.ndarray:
        """Cells where both surface and bed are defined."""
        return np.isfinite(self.surface) & np.isfinite(self.bed)

    @property
    def thickness(self) -> np.ndarray:
        """Ice thickness at ``base_year`` (m)."""
        return self.surface - self.bed

    @property
    def cell_area_m2(self) -> float:
        return abs(self.transform[1] * self.transform[5])

    def xy(self, rowcol: tuple[int, int]) -> tuple[float, float]:
        x0, dx, _, y0, _, dy = self.transform
        r, c = rowcol
        return x0 + (c + 0.5) * dx, y0 + (r + 0.5) * dy

    def save(self, path: str | Path) -> Path:
        """Write to a compressed ``.npz``."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            name=self.name,
            surface=self.surface.astype(np.float32),
            bed=self.bed.astype(np.float32),
            seed_lake=np.asarray(self.seed_lake),
            seed_target=np.asarray(self.seed_target),
            base_year=self.base_year,
            transform=np.asarray(self.transform),
            crs=self.crs or "",
            vertical_datum=self.vertical_datum or "",
            meta=np.asarray(repr(self.meta)),
        )
        return path

    @classmethod
    def load(cls, path: str | Path) -> DamObject:
        import ast

        with np.load(path, allow_pickle=False) as z:
            return cls(
                name=str(z["name"]),
                surface=z["surface"].astype(np.float64),
                bed=z["bed"].astype(np.float64),
                seed_lake=tuple(int(v) for v in z["seed_lake"]),
                seed_target=tuple(int(v) for v in z["seed_target"]),
                base_year=int(z["base_year"]),
                transform=tuple(float(v) for v in z["transform"]),
                crs=str(z["crs"]) or None,
                vertical_datum=str(z["vertical_datum"]) or None,
                meta=ast.literal_eval(str(z["meta"])),
            )
