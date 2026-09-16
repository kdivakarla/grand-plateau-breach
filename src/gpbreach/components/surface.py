"""Surface elevation models: how the glacier surface evolves after ``base_year``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..dam import DamObject


@dataclass(frozen=True)
class LinearThinning:
    """Spatially uniform linear surface lowering. [D-004]

        S(x, t) = S(x, t0) + dhdt * (t - t0)

    ``dhdt`` is negative for thinning. The baseline uses -9.05 m yr-1, derived by
    the owner from OIB/DEM transect differences and recovered exactly from the
    legacy QGIS rasters (docs/baseline_method.md).

    No elevation dependence, no seasonality, no dynamic response: the baseline
    reproduces the owner's method, it does not improve it (CLAUDE.md Section 8.1).
    """

    name: str = "linear_thinning"

    def evaluate(self, t_years: np.ndarray, dam: DamObject, dhdt: np.ndarray) -> np.ndarray:
        """Surface lowering (m, negative) after ``t_years``, shaped ``(n_real, n_times)``."""
        t = np.atleast_1d(np.asarray(t_years, dtype=float))
        d = np.atleast_1d(np.asarray(dhdt, dtype=float))
        return d[:, None] * t[None, :]

    def elevation_offset_rate(self, dhdt: np.ndarray) -> np.ndarray:
        """Rate at which the surface field shifts (m yr-1), used by the engine."""
        return np.asarray(dhdt, dtype=float)


SURFACE_MODELS = {LinearThinning.name: LinearThinning}
