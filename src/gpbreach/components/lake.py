"""Lake level models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..dam import DamObject


@dataclass(frozen=True)
class FixedLevel:
    """Lake surface held at a constant elevation for all time. [D-007]

    The baseline uses 110 m for the upper (LGP) lake and 17 m for Alsek Lake,
    both recovered exactly from the legacy QGIS rasters. Whether each is a fixed
    sill, an ice spillway, or a river stage is an open question the owner must
    settle (docs/decisions.md D-007).
    """

    name: str = "fixed"

    def evaluate(self, t_years: np.ndarray, dam: DamObject, level: np.ndarray) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t_years, dtype=float))
        lv = np.atleast_1d(np.asarray(level, dtype=float))
        return np.broadcast_to(lv[:, None], (lv.size, t.size)).copy()


LAKE_MODELS = {FixedLevel.name: FixedLevel}
