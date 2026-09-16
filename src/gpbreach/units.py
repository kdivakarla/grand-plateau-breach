"""Unit and datum helpers.

All elevations and thicknesses are metres. Conversions to metres of water
equivalent live here and nowhere else (CLAUDE.md Section 5).
"""

from __future__ import annotations

from .constants import RHO_ICE, RHO_WATER


class DatumMismatch(ValueError):
    """Raised when elevations from two sources do not share a vertical datum."""


def assert_same_datum(*datums: str | None) -> str:
    """Return the common vertical datum, or raise.

    ``None`` means "unverified" and is always an error: the owner confirms every
    datum, Claude Code never assumes one (CLAUDE.md Section 5).
    """
    if any(d is None for d in datums):
        raise DatumMismatch(
            "At least one dataset has an unverified vertical datum. Record it in "
            "docs/data_sources.md before combining these elevations."
        )
    unique = set(datums)
    if len(unique) != 1:
        raise DatumMismatch(f"Vertical datums differ: {sorted(unique)}")
    return unique.pop()


def ice_to_water_equivalent(thickness_m: float) -> float:
    """Convert an ice thickness (m) to metres of water equivalent head."""
    return thickness_m * RHO_ICE / RHO_WATER


def water_equivalent_to_ice(head_m_we: float) -> float:
    """Convert metres of water equivalent head to an equivalent ice thickness (m)."""
    return head_m_we * RHO_WATER / RHO_ICE
