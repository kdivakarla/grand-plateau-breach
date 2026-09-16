"""Bed elevation models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..dam import DamObject


@dataclass(frozen=True)
class StaticBed:
    """The owner's ``Bedrock_Aligned.tif`` used as-is, with no perturbation. [D-006]

    The baseline treats the bed as exactly known. Spatially correlated bed error
    is Phase 3 (CLAUDE.md Section 9) and will arrive as a sibling implementation
    here; it must not be added to this one.
    """

    name: str = "static"

    def evaluate(self, dam: DamObject, params: dict | None = None) -> np.ndarray:
        return dam.bed


BED_MODELS = {StaticBed.name: StaticBed}
