"""Project the ice surface forward in time.

Produces a surface raster for any future year, which then plugs into the existing
machinery as just another entry in the config's ``surfaces:`` list -- so the pass
search, path, volumes and hydrograph all run on it unchanged.

Three thinning models
---------------------
``uniform``
    One rate everywhere: ``z_s(t) = z_s(t0) + dhdt * (t - t0)``. This is D-004's
    -9.05 m/yr. Because the shift is constant, the whole head field shifts by
    ``kappa * dhdt * dt`` and the **pass does not move** -- the onset year is then
    available in closed form and projecting surfaces is unnecessary.

``by_elevation``
    A different rate per elevation band, from Loso's lidar (D-026). Thinning is
    strongest in the 200-300 m band, so the surface does *not* shift uniformly and
    **the pass can move**. This is the case where projected surfaces earn their
    keep.

``accelerating``
    ``dhdt(t) = dhdt0 + accel * (t - t0)``, integrating to a quadratic. D-026 found
    the rate rising from about -5 m/yr (2005-09) to -12 m/yr (2017-20); a linear
    extrapolation of the long-term mean therefore predicts a *later* breach than
    the recent trend implies.

Surfaces are clipped to the ice mask: lake and land cells are left untouched,
because thinning a lake surface is meaningless.
"""

from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)


def uniform(zs: np.ndarray, ice: np.ndarray, dhdt: float, dt_years: float) -> np.ndarray:
    """Constant rate everywhere on the ice (m/yr, negative for thinning)."""
    return np.where(ice, zs + dhdt * dt_years, zs)


def by_elevation(
    zs: np.ndarray,
    ice: np.ndarray,
    bands: list[tuple[float, float, float]],
    dt_years: float,
    default_dhdt: float = 0.0,
) -> np.ndarray:
    """Rate varying with elevation.

    ``bands`` is a list of ``(z_low, z_high, dhdt)``; the rate is chosen by the
    cell's **current** surface elevation. Cells outside every band take
    ``default_dhdt``.

    The lookup uses the surface at the start of the step, so for long projections
    call this in annual steps rather than one jump -- a cell that thins out of its
    band should pick up the next band's rate.
    """
    rate = np.full(zs.shape, default_dhdt, dtype=float)
    for lo, hi, r in bands:
        rate = np.where(ice & (zs >= lo) & (zs < hi), r, rate)
    return np.where(ice, zs + rate * dt_years, zs)


def accelerating(
    zs: np.ndarray, ice: np.ndarray, dhdt0: float, accel: float, dt_years: float
) -> np.ndarray:
    """Rate changing linearly in time: integral of ``dhdt0 + accel*t`` over ``dt``."""
    drop = dhdt0 * dt_years + 0.5 * accel * dt_years**2
    return np.where(ice, zs + drop, zs)


def project(
    zs: np.ndarray,
    ice: np.ndarray,
    years: list[int],
    base_year: int,
    model: str = "uniform",
    step_years: float = 1.0,
    **kw,
) -> dict[int, np.ndarray]:
    """Surfaces for each requested year, stepped annually so rates can update.

    Returns ``{year: surface}``. Thinning never raises the surface below the bed --
    that is the caller's concern, but a warning is logged if a surface goes
    negative in thickness against a supplied ``bed``.
    """
    bed = kw.pop("bed", None)
    out: dict[int, np.ndarray] = {}
    cur = zs.copy()
    year = base_year
    for target in sorted(years):
        while year < target:
            dt = min(step_years, target - year)
            if model == "uniform":
                cur = uniform(cur, ice, kw["dhdt"], dt)
            elif model == "by_elevation":
                cur = by_elevation(cur, ice, kw["bands"], dt, kw.get("default_dhdt", 0.0))
            elif model == "accelerating":
                elapsed = year - base_year
                rate = kw["dhdt0"] + kw["accel"] * elapsed
                cur = uniform(cur, ice, rate, dt)
            else:
                raise ValueError(f"unknown model {model!r}")
            year += dt
        out[target] = cur.copy()
        if bed is not None:
            thin = np.isfinite(bed) & ice & (cur - bed < 0)
            if thin.any():
                log.warning("%d cells have negative thickness by %d", int(thin.sum()), target)
    return out
