"""Breach criteria.

``FlotationConnectivity`` is the baseline: the owner's QGIS method, ported
exactly. Its derivation is worth stating in full because it is what makes the
Monte Carlo framework cheap.

Legacy statement (owner's words)
--------------------------------
Split the ice dam at lake level. The part above lake level thins at ``dhdt``.
The dam is breachable where the remaining ice above lake level falls below a
fraction ``r`` of the ice below lake level. Breach occurs in the first year a
*connected pathway* of such cells spans the dam.

    H_above(t) = S + dhdt * t - L        H_below = L - B
    buoyant  <=>  H_above(t) <= r * H_below

Algebra
-------
    S + dhdt*t - L <= r * (L - B)
    S + dhdt*t + r*B - (1 + r)*L <= 0

Divide through by (1 + r) and write ``f = 1 / (1 + r)``:

    f*(S + dhdt*t) + (1 - f)*B <= L                                        (1)

which is exactly hydrostatic flotation, ``h_lake >= f * H_ice`` with
``h_lake = L - B`` and ``H_ice = S + dhdt*t - B``. So the owner's ratio ``r``
and the flotation fraction ``f`` are the same parameter in different clothes:

    f = 1 / (1 + r)          r = (1 - f) / f
    r = 0.1  <->  f = 0.909091

Why this shape matters
----------------------
Define the **static flotation field** on the base-year grids:

    Phi_f(x) = f * S(x) + (1 - f) * B(x)                                   (2)

Then (1) reads ``Phi_f(x) <= L - f*dhdt*t``: the buoyant region is a sub-level
set of a *time-independent* field under a threshold that rises linearly. The
breach is the first connection of the two basins, so it happens exactly when
the threshold reaches the saddle value ``Phi_crit(f)`` between them:

    t_breach = (Phi_crit(f) - L) / (-f * dhdt)                             (3)

``Phi_crit`` depends on ``f`` and on the bed/surface grids -- and on nothing
else. Lake level and thinning rate enter only through the closed form (3).
Sampling over ``L`` and ``dhdt`` is therefore free, and sampling over ``f``
needs only a one-dimensional precomputed curve (see ``analysis.fcurve``). Only
a perturbed *bed* forces the saddle to be recomputed.

See D-001, D-003, D-005.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..dam import DamObject
from ..saddle import locate_saddle, saddle_threshold


def flotation_field(surface: np.ndarray, bed: np.ndarray, f: float) -> np.ndarray:
    """Static flotation field ``Phi_f = f*S + (1-f)*B`` -- equation (2)."""
    if not 0.0 < f <= 1.0:
        raise ValueError(f"flotation fraction f must be in (0, 1], got {f!r}")
    return f * surface + (1.0 - f) * bed


@dataclass(frozen=True)
class FlotationConnectivity:
    """Breach when a connected buoyant pathway first links the two basins. [D-001]

    Parameters
    ----------
    connectivity:
        4 or 8 neighbour adjacency. The baseline dates are identical either way.
    """

    connectivity: int = 8
    name: str = "flotation_connectivity"

    def critical_threshold(
        self,
        surface: np.ndarray,
        bed: np.ndarray,
        f: float,
        seed_a: tuple[int, int],
        seed_b: tuple[int, int],
    ) -> float:
        """Saddle value ``Phi_crit(f)`` between the two seeds (m)."""
        return saddle_threshold(flotation_field(surface, bed, f), seed_a, seed_b, self.connectivity)

    def critical_threshold_for_dam(self, dam: DamObject, f: float) -> float:
        return self.critical_threshold(dam.surface, dam.bed, f, dam.seed_lake, dam.seed_target)

    @staticmethod
    def breach_time(
        phi_crit: np.ndarray, level: np.ndarray, dhdt: np.ndarray, f: np.ndarray
    ) -> np.ndarray:
        """Years after ``base_year`` until breach -- equation (3). Vectorised.

        Returns ``inf`` where the dam never breaches (no thinning, or the basins
        are already connected and the lake is below the saddle for all time).
        """
        phi_crit = np.asarray(phi_crit, dtype=float)
        level = np.asarray(level, dtype=float)
        dhdt = np.asarray(dhdt, dtype=float)
        f = np.asarray(f, dtype=float)
        rate = -f * dhdt  # m yr-1 at which the threshold rises; > 0 while thinning
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (phi_crit - level) / rate
        t = np.where(rate <= 0, np.inf, t)
        return np.where(np.isfinite(phi_crit), np.maximum(t, 0.0), np.inf)

    @staticmethod
    def threshold_at(
        level: np.ndarray, dhdt: np.ndarray, f: np.ndarray, t_years: np.ndarray
    ) -> np.ndarray:
        """Buoyancy threshold ``L - f*dhdt*t`` at time ``t`` (m)."""
        return np.asarray(level, dtype=float) - np.asarray(f, dtype=float) * np.asarray(
            dhdt, dtype=float
        ) * np.asarray(t_years, dtype=float)

    def buoyant_mask(
        self, dam: DamObject, f: float, level: float, dhdt: float, t_years: float
    ) -> np.ndarray:
        """Cells meeting the buoyancy condition at time ``t_years``."""
        phi = flotation_field(dam.surface, dam.bed, f)
        return np.isfinite(phi) & (phi <= self.threshold_at(level, dhdt, f, t_years))

    def breach_cells(self, dam: DamObject, f: float, phi_crit: float) -> np.ndarray:
        """``(n, 2)`` cells where the pathway completes -- the breach point."""
        return locate_saddle(
            flotation_field(dam.surface, dam.bed, f),
            phi_crit,
            dam.seed_lake,
            dam.seed_target,
            self.connectivity,
        )

    def legacy_margin(
        self, dam: DamObject, f: float, level: float, dhdt: float, t_years: float
    ) -> np.ndarray:
        """The owner's QGIS raster, reproduced cell-for-cell.

        ``Buoyancy_State_<year>.tif`` = ``S + r*B - (1 + r)*L - (-dhdt)*t`` with
        ``r = (1-f)/f``; negative means buoyant. This is ``(1/f)`` times the
        margin used internally, kept only so the port can be checked against the
        owner's files (docs/baseline_method.md).
        """
        r = (1.0 - f) / f
        return dam.surface + r * dam.bed - (1.0 + r) * level + dhdt * t_years


CRITERIA = {FlotationConnectivity.name: FlotationConnectivity}
