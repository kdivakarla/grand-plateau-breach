"""Precomputed ``Phi_crit(f)`` curve.

The saddle search is the only expensive operation in the model, and it depends
on exactly one sampled parameter: the flotation fraction ``f``. Tabulating
``Phi_crit`` on a grid of ``f`` once turns a Monte Carlo over ``f``, lake level
and thinning rate into pure vectorised arithmetic (components.criterion,
equation 3).

This is a performance device, not a physics change: on the ``f`` grid points the
tabulated values are exact, and ``spacing`` is chosen small enough that
interpolation error is far below the parameter uncertainty. Report the
interpolation residual with ``FCurve.max_interp_error`` before trusting it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..components.criterion import FlotationConnectivity
from ..dam import DamObject


@dataclass(frozen=True)
class FCurve:
    """Tabulated saddle value against flotation fraction."""

    f: np.ndarray
    phi_crit: np.ndarray
    dam_name: str
    connectivity: int

    def __call__(self, f: np.ndarray) -> np.ndarray:
        f = np.asarray(f, dtype=float)
        if np.any(f < self.f[0]) or np.any(f > self.f[-1]):
            raise ValueError(
                f"f outside tabulated range [{self.f[0]:.4f}, {self.f[-1]:.4f}]; "
                "widen the curve rather than extrapolating."
            )
        return np.interp(f, self.f, self.phi_crit)

    def max_interp_error(
        self, dam: DamObject, n_checks: int = 9, rng: np.random.Generator | None = None
    ) -> float:
        """Largest |interpolated - exact| at points between the tabulated nodes."""
        rng = rng or np.random.default_rng(0)
        probes = rng.uniform(self.f[0], self.f[-1], n_checks)
        crit = FlotationConnectivity(connectivity=self.connectivity)
        exact = np.array([crit.critical_threshold_for_dam(dam, float(p)) for p in probes])
        return float(np.max(np.abs(self(probes) - exact)))

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            path,
            f=self.f,
            phi_crit=self.phi_crit,
            dam_name=self.dam_name,
            connectivity=self.connectivity,
        )
        return path

    @classmethod
    def load(cls, path: str | Path) -> FCurve:
        with np.load(path) as z:
            return cls(z["f"], z["phi_crit"], str(z["dam_name"]), int(z["connectivity"]))


def build_fcurve(
    dam: DamObject,
    f_min: float = 0.80,
    f_max: float = 0.98,
    spacing: float = 0.005,
    connectivity: int = 8,
) -> FCurve:
    """Tabulate ``Phi_crit`` over ``[f_min, f_max]``."""
    n = round((f_max - f_min) / spacing) + 1
    grid = np.linspace(f_min, f_max, n)
    crit = FlotationConnectivity(connectivity=connectivity)
    values = np.array([crit.critical_threshold_for_dam(dam, float(v)) for v in grid])
    return FCurve(grid, values, dam.name, connectivity)
