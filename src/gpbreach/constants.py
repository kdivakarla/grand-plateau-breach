"""Physical constants.

Only this module may define densities or other physical constants. Any change
here is a physics change and requires owner approval (CLAUDE.md Section 2.2).
"""

from __future__ import annotations

#: Density of pure glacier ice (kg m-3). NOT used by the baseline criterion --
#: the baseline prescribes the flotation fraction directly. See D-002.
RHO_ICE: float = 917.0

#: Density of fresh water (kg m-3). See D-002.
RHO_WATER: float = 1000.0

#: Flotation fraction implied by the density ratio above. Provided for
#: comparison only; the baseline uses ``LEGACY_FLOTATION_FRACTION``. See D-002.
DENSITY_FLOTATION_FRACTION: float = RHO_ICE / RHO_WATER

#: Flotation fraction reproduced from the owner's QGIS rasters. The legacy
#: method compared ice above lake level against a fraction r = 0.1 of ice below
#: lake level, which is algebraically f = 1 / (1 + r). See D-001 and D-003.
LEGACY_THICKNESS_RATIO: float = 0.1
LEGACY_FLOTATION_FRACTION: float = 1.0 / (1.0 + LEGACY_THICKNESS_RATIO)


def ratio_to_fraction(r: float) -> float:
    """Convert the legacy above/below thickness ratio ``r`` to flotation fraction ``f``."""
    return 1.0 / (1.0 + r)


def fraction_to_ratio(f: float) -> float:
    """Convert flotation fraction ``f`` to the legacy above/below thickness ratio ``r``."""
    return (1.0 - f) / f
