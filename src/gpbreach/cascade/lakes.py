"""Stage-area-volume curves and the released volume V_w.

V_w drives everything downstream: the empirical relations go roughly as
``V^0.66``, so a factor-two error in volume is a ~1.6x error in peak discharge.
It therefore deserves more care than a single number.

How far LGP actually drains
---------------------------
Not to GPL's level. The drainage path has a **bed sill** -- the highest bed
elevation along the route -- and LGP cannot drain below it, because once the lake
falls to the sill there is no longer a downhill path. So

    h_final = max(h_GPL, bed_sill)

and on the 2018 geometry the sill (58.71 m) controls, not GPL (27.6 m): 58.6 m of
drawdown rather than 89.7 m.

Bathymetry
----------
Volume below the waterline needs the lake floor. Where survey points exist,
``interpolate_bathymetry`` grids them. Where they do not, ``vertical_wall_bottom``
gives a strict **upper bound** by assuming the walls are vertical, so the lake
keeps its surface area all the way down. A real basin narrows with depth, so the
true volume is always smaller. Use the bound to bracket, never as an estimate.
"""

from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)


def vertical_wall_bottom(mask: np.ndarray, floor_elevation: float) -> np.ndarray:
    """Bottom elevation for a prism: every lake cell at ``floor_elevation`` (m).

    Yields the **largest** volume any basin with this outline could hold between
    two stages, because area never shrinks with depth. The honest use is as an
    upper bound while real bathymetry is missing.
    """
    return np.where(mask, float(floor_elevation), np.nan)


def conic_bottom(mask: np.ndarray, h_lake: float, max_depth: float) -> np.ndarray:
    """Bottom elevation for a cone-like basin, as a companion lower bound.

    Depth falls off linearly with normalised distance from the shore, so the
    deepest point sits at the centre of the lake and the margins are shallow.
    Crude, but it brackets the prism from the other side: a real basin usually
    lies between the two.
    """
    from scipy.ndimage import distance_transform_edt

    dist = distance_transform_edt(mask)
    if dist.max() <= 0:
        return np.where(mask, h_lake, np.nan)
    depth = max_depth * (dist / dist.max())
    return np.where(mask, h_lake - depth, np.nan)


def interpolate_bathymetry(
    xy: np.ndarray,
    depths: np.ndarray,
    mask: np.ndarray,
    h_lake: float,
    transform,
    method: str = "linear",
    fill_with_max: bool = True,
) -> np.ndarray:
    """Grid sounding depths into a bottom-elevation raster (m).

    Parameters
    ----------
    xy: ``(n, 2)`` survey positions in the raster CRS (m).
    depths: measured water depths (m, positive down).
    mask: lake cells to fill.
    h_lake: lake surface elevation on the survey date (m); bottom = ``h_lake - depth``.
    method: passed to ``scipy.interpolate.griddata``.
    fill_with_max: outside the soundings' convex hull ``griddata`` returns NaN.
        True fills those cells with the deepest sounding (conservative for volume);
        False leaves them NaN so the caller sees the gap.

    Censored soundings (recorded as "> 275 m") must already be resolved to a
    number by the caller -- this function cannot know what the cap means.
    """
    from scipy.interpolate import griddata

    rows, cols = np.nonzero(mask)
    x = transform.c + (cols + 0.5) * transform.a
    y = transform.f + (rows + 0.5) * transform.e
    bottom_pts = float(h_lake) - np.asarray(depths, dtype=float)
    grid = griddata(np.asarray(xy, dtype=float), bottom_pts, (x, y), method=method)

    n_nan = int(np.isnan(grid).sum())
    if n_nan and fill_with_max:
        log.info(
            "%d of %d lake cells lie outside the sounding hull; filled with the deepest sounding",
            n_nan,
            grid.size,
        )
        grid = np.where(np.isnan(grid), bottom_pts.min(), grid)
    out = np.full(mask.shape, np.nan)
    out[rows, cols] = grid
    return out


def stage_curve(bottom: np.ndarray, mask: np.ndarray, levels: np.ndarray, cell_area_m2: float):
    """Area and volume at each stage.

    ``A(h)`` is the area of lake cells whose floor lies below ``h``; ``V(h)`` is
    the water volume between the floor and ``h``, accumulated from the lake floor
    up. Returns a DataFrame with ``h`` (m), ``A_m2``, ``A_km2``, ``V_m3``, ``V_Mm3``.
    """
    import pandas as pd

    valid = mask & np.isfinite(bottom)
    b = bottom[valid]
    rows = []
    for h in np.asarray(levels, dtype=float):
        wet = b < h
        area = float(wet.sum()) * cell_area_m2
        vol = float(np.sum(h - b[wet])) * cell_area_m2 if wet.any() else 0.0
        rows.append(
            {"h": float(h), "A_m2": area, "A_km2": area / 1e6, "V_m3": vol, "V_Mm3": vol / 1e6}
        )
    return pd.DataFrame(rows)


def volume_at(curve, h: float) -> float:
    """Interpolate V (m^3) at stage ``h`` from a stage curve."""
    return float(np.interp(h, curve["h"].to_numpy(), curve["V_m3"].to_numpy()))


def drainage_floor(h_gpl: float, bed_sill: float) -> tuple[float, str]:
    """Lowest stage LGP can reach, and which constraint sets it."""
    if bed_sill >= h_gpl:
        return float(bed_sill), "bed_sill"
    return float(h_gpl), "receiving_lake"


def released_volume(curve, h_start: float, h_gpl: float, bed_sill: float) -> dict:
    """Water released as LGP falls from ``h_start`` to its drainage floor.

    Returns the volume, the floor, which constraint set it, and -- because the
    difference matters for interpretation -- what would have drained had the sill
    not been there.
    """
    h_final, control = drainage_floor(h_gpl, bed_sill)
    v_start, v_final = volume_at(curve, h_start), volume_at(curve, h_final)
    vw = max(v_start - v_final, 0.0)
    # Only meaningful if the curve actually extends below the receiving lake;
    # otherwise the prism floor sits at the sill and the answer is a tautology.
    floor = float(curve["h"].to_numpy().min())
    v_nosill = max(v_start - volume_at(curve, h_gpl), 0.0) if floor <= h_gpl else None
    return {
        "h_start_m": float(h_start),
        "h_final_m": h_final,
        "controlled_by": control,
        "drawdown_m": float(h_start - h_final),
        "V_start_m3": v_start,
        "V_final_m3": v_final,
        "Vw_m3": vw,
        "Vw_Mm3": vw / 1e6,
        "Vw_km3": vw / 1e9,
        "Vw_if_no_sill_Mm3": None if v_nosill is None else v_nosill / 1e6,
        "sill_withholds_Mm3": None if v_nosill is None else (v_nosill - vw) / 1e6,
        "sill_withholds_note": None
        if v_nosill is not None
        else "unknown: the stage curve does not extend below the receiving lake (needs bathymetry)",
    }


def connected_stage_area(
    dem: np.ndarray,
    seed_mask: np.ndarray,
    levels: np.ndarray,
    cell_area_m2: float,
    base_level: float,
    connectivity: int = 8,
):
    """Area, and volume **added**, as a lake rises above ``base_level``.

    At each stage, cells below that stage are flood-filled from the existing lake,
    so a basin that merely sits below the level but is cut off by higher ground is
    not counted.

    ``base_level`` matters more than it looks. The project's surface DEM has the
    lakes masked out, so lake cells are NaN. Treating them as "no data" makes the
    existing lake surface contribute nothing, and GPL's 45 km² -- by far the
    largest term -- silently vanishes from the volume. Seed cells are therefore
    assigned an effective elevation of ``base_level``: the water surface already
    there. Volume is then what the flood *adds*, which is what the routing stage
    needs.
    """
    import pandas as pd
    from scipy import ndimage

    struct = (
        np.ones((3, 3), bool)
        if connectivity == 8
        else np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool)
    )
    # Existing lake surface stands in for the missing bathymetry.
    effective = np.where(seed_mask, float(base_level), dem)
    rows = []
    for h in np.asarray(levels, dtype=float):
        below = (np.isfinite(effective) & (effective < h)) | seed_mask
        lab, _ = ndimage.label(below, structure=struct)
        keep = set(np.unique(lab[seed_mask])) - {0}
        conn = np.isin(lab, list(keep)) if keep else np.zeros_like(below)
        area = float(conn.sum()) * cell_area_m2
        depth = np.where(conn & np.isfinite(effective), h - effective, 0.0)
        vol = float(np.sum(np.maximum(depth, 0.0))) * cell_area_m2
        rows.append(
            {"h": float(h), "A_m2": area, "A_km2": area / 1e6, "V_m3": vol, "V_Mm3": vol / 1e6}
        )
    return pd.DataFrame(rows)


def detect_spill(curve, area_jump_frac: float = 0.10, outlier_factor: float = 8.0) -> dict | None:
    """First stage at which the lake's area jumps, i.e. it spills into a new basin.

    A stage curve for a closed basin grows smoothly. A step change means the water
    has topped a divide and flooded an adjoining depression. That matters twice
    over: the new basin can swallow a large volume for almost no further rise, so
    interpolating a stage from a volume **across** the step is meaningless; and
    beyond it the single-basin assumption behind the curve no longer holds.

    Detection is an **outlier** test, not a fixed threshold. Near the apex of a
    smoothly widening basin the relative area growth per step is naturally large
    -- a cone gains 20 % of its area in a step close to its floor -- so a bare
    cutoff flags ordinary basins. A step must exceed both ``area_jump_frac`` and
    ``outlier_factor`` times the median step to count.

    Returns ``None`` if the curve is smooth, otherwise the stage, the area jump,
    and the capacity of the original basin up to that point.
    """
    h = curve["h"].to_numpy()
    a = curve["A_m2"].to_numpy()
    v = curve["V_m3"].to_numpy()
    if h.size < 3:
        return None
    rel = np.diff(a) / np.maximum(a[:-1], 1.0)
    if rel.size == 0:
        return None
    i = int(np.argmax(rel))
    typical = float(np.median(rel))
    if rel[i] < area_jump_frac or rel[i] < outlier_factor * typical:
        return None
    return {
        "spill_stage_m": float(h[i + 1]),
        "stage_before_m": float(h[i]),
        "area_before_km2": float(a[i] / 1e6),
        "area_after_km2": float(a[i + 1] / 1e6),
        "area_jump_frac": float(rel[i]),
        "median_step_frac": typical,
        "capacity_before_spill_Mm3": float(v[i] / 1e6),
        "volume_absorbed_by_step_Mm3": float((v[i + 1] - v[i]) / 1e6),
    }
