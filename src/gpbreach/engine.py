"""The engine: margin over time, first crossing per realization.

Generic over realizations; the baseline is simply ``n_realizations = 1``.

The expensive part is the saddle value ``Phi_crit(f)``, which depends only on
``f`` and the (static) bed and surface grids -- not on lake level or thinning
rate. The engine therefore computes it once per *distinct* ``f`` and then solves
for every realization in closed form, vectorised (see components.criterion).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .components.criterion import FlotationConnectivity, flotation_field
from .config import Config
from .dam import DamObject
from .sampling import build_sample_table

log = logging.getLogger(__name__)


@dataclass
class RunResult:
    run_id: str
    dam_name: str
    samples: pd.DataFrame
    margin: pd.DataFrame  # index realization, columns year; m (positive = holds)
    timeseries: pd.DataFrame  # per-year diagnostics for realization 0
    phi_crit: dict[float, float]  # f -> saddle value (m)
    breach_cells: np.ndarray = field(default_factory=lambda: np.empty((0, 2), int))
    meta: dict = field(default_factory=dict)

    @property
    def breach_year_continuous(self) -> pd.Series:
        return self.samples["breach_year_continuous"]

    @property
    def breach_year(self) -> pd.Series:
        """Legacy reporting: first annual raster showing a connected pathway."""
        return self.samples["breach_year"]


def run(cfg: Config, dam: DamObject | None = None) -> RunResult:
    """Execute one run and return per-realization breach dates plus diagnostics."""
    dam = dam if dam is not None else DamObject.load(cfg.dam_path)
    criterion = FlotationConnectivity(connectivity=int(cfg.components.get("connectivity", 8)))

    samples = build_sample_table(cfg)
    f = samples["flotation_fraction"].to_numpy()
    level = samples["lake_level_m"].to_numpy()
    dhdt = samples["dhdt_m_per_yr"].to_numpy()

    # --- expensive step, once per distinct flotation fraction --------------
    # Grouping is on the exact float value, never a rounded one: the saddle is a
    # specific cell value of Phi_f, and downstream code locates that cell by
    # equality. A perturbed f must be used identically in both places.
    phi_crit_by_f: dict[float, float] = {}
    for f_value in np.unique(f):
        phi_crit_by_f[float(f_value)] = criterion.critical_threshold_for_dam(dam, float(f_value))
        log.info("Phi_crit(f=%.12g) = %.4f m", f_value, phi_crit_by_f[float(f_value)])
    if len(phi_crit_by_f) > 50:
        log.warning(
            "%d distinct flotation fractions required %d saddle solves; use "
            "analysis.fcurve for ensembles of this size.",
            len(phi_crit_by_f),
            len(phi_crit_by_f),
        )
    phi_crit = np.array([phi_crit_by_f[float(v)] for v in f])

    # --- closed-form solve, vectorised over realizations -------------------
    t_breach = criterion.breach_time(phi_crit, level, dhdt, f)
    year_continuous = dam.base_year + t_breach
    year_annual = np.where(np.isfinite(year_continuous), np.ceil(year_continuous), np.inf)

    samples = samples.assign(
        phi_crit_m=phi_crit,
        t_breach_yr=t_breach,
        breach_year_continuous=year_continuous,
        breach_year=year_annual,
    )

    # --- margin over time at the saddle ------------------------------------
    years = np.asarray(cfg.years, dtype=float)
    t = years - dam.base_year
    margin = phi_crit[:, None] + (f * dhdt)[:, None] * t[None, :] - level[:, None]
    margin_df = pd.DataFrame(margin, index=samples.index, columns=[int(y) for y in years])

    timeseries = _diagnostics(
        dam, criterion, years, float(f[0]), float(level[0]), float(dhdt[0]), float(phi_crit[0])
    )
    cells = criterion.breach_cells(dam, float(f[0]), float(phi_crit[0]))

    return RunResult(
        run_id=cfg.run_id,
        dam_name=dam.name,
        samples=samples,
        margin=margin_df,
        timeseries=timeseries,
        phi_crit=phi_crit_by_f,
        breach_cells=cells,
        meta={
            "connectivity": criterion.connectivity,
            "base_year": dam.base_year,
            "n_realizations": len(samples),
        },
    )


def _diagnostics(
    dam: DamObject,
    criterion: FlotationConnectivity,
    years: np.ndarray,
    f: float,
    level: float,
    dhdt: float,
    phi_crit: float,
) -> pd.DataFrame:
    """Year-by-year values at the breach point, for the baseline report."""
    phi = flotation_field(dam.surface, dam.bed, f)
    cells = criterion.breach_cells(dam, f, phi_crit)
    if cells.size:
        r, c = cells[0]
        at_saddle = True
    else:  # the basins never connect within the domain
        log.warning("no saddle cell found; reporting diagnostics at the lake seed instead")
        r, c = dam.seed_lake
        at_saddle = False
    s0, b0 = float(dam.surface[r, c]), float(dam.bed[r, c])
    t = years - dam.base_year
    surface = s0 + dhdt * t
    thickness = surface - b0
    lake_depth = level - b0
    rows = pd.DataFrame(
        {
            "year": years.astype(int),
            "surface_m": surface,
            "bed_m": b0,
            "ice_thickness_m": thickness,
            "lake_depth_m": lake_depth,
            "flotation_depth_m": f * thickness,
            "margin_m": f * thickness - lake_depth,
            "threshold_m": criterion.threshold_at(level, dhdt, f, t),
            "phi_crit_m": phi_crit,
            "connected": criterion.threshold_at(level, dhdt, f, t) >= phi_crit,
        }
    )
    rows.attrs["at_saddle"] = at_saddle
    rows.attrs["breach_cell"] = (int(r), int(c))
    rows.attrs["breach_xy"] = dam.xy((int(r), int(c)))
    rows.attrs["phi_at_cell"] = float(phi[r, c])
    return rows
