"""LGP→GPL cascade: pass search, drainage path, lake curves, inflow hydrograph.

Implements the owner's *Pass Search & Inflow Hydrograph Build Spec* (2026-10-03).
See docs/phase2_cascade_plan.md for the build order and docs/decisions.md for the
parameter choices.
"""

from . import grid

__all__ = ["grid"]
