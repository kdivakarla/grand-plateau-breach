"""Swappable physics components, selected by name from YAML config.

Every implementation's docstring cites the equation it applies and the
``docs/decisions.md`` entry that authorised it (CLAUDE.md Section 3).
"""

from .bed import BED_MODELS, StaticBed
from .criterion import CRITERIA, FlotationConnectivity, flotation_field
from .lake import LAKE_MODELS, FixedLevel
from .surface import SURFACE_MODELS, LinearThinning

__all__ = [
    "BED_MODELS",
    "CRITERIA",
    "LAKE_MODELS",
    "SURFACE_MODELS",
    "FixedLevel",
    "FlotationConnectivity",
    "LinearThinning",
    "StaticBed",
    "flotation_field",
]
