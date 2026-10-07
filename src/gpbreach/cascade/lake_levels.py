"""Lake-surface elevation time series from ICESat-2 ATL06.

Pulls quality-screened ATL06 heights inside each lake outline and reduces them to
one elevation per lake per overpass.

Three screens, and all three earn their place
---------------------------------------------
1. **Product quality.** ``atl06_quality_summary == 0`` and agreement with the
   granule's own reference DEM. On a cloud-affected granule the unscreened heights
   can be hundreds of metres wrong while still looking like ordinary numbers --
   the 2026-07-19 granule has **zero** usable segments over this area.
2. **Inside the lake.** Points are selected by polygon, with an inward buffer so
   shoreline and calving-front returns do not contaminate the water surface.
3. **Flatness.** A lake surface is level, so scatter is a tell. Points more than
   ``mad_k`` robust deviations from the median are dropped, and an epoch whose
   residual scatter stays high is flagged rather than reported as a clean value.

The outlines are fixed at their 2018 extent (D-015) while the lakes themselves
move over the record, which is what the inward buffer is really protecting
against. Treat an epoch with few points or high scatter as suspect.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ..io import atl06

log = logging.getLogger(__name__)

#: class code -> label, matching configs/cascade_data.yaml
DEFAULT_LAKES = {2: "LGP", 3: "GPL", 4: "Alsek"}


def robust_level(h: np.ndarray, mad_k: float = 3.0) -> dict:
    """Median elevation with a robust outlier cut.

    Returns the level, its scatter, and how many points survived -- all three
    matter, because a confident-looking median over four points on a windy
    shoreline is not a lake level.
    """
    h = np.asarray(h, dtype=float)
    h = h[np.isfinite(h)]
    if h.size == 0:
        return {"n": 0, "level_m": np.nan, "mad_m": np.nan, "std_m": np.nan, "n_rejected": 0}
    med = float(np.median(h))
    mad = float(np.median(np.abs(h - med))) * 1.4826  # ~sigma for a normal
    # mad == 0 means every point is identical; keep them all rather than reject all.
    keep = np.abs(h - med) <= mad_k * mad if mad > 0 else np.ones(h.shape, bool)
    kept = h[keep]
    if kept.size == 0:
        kept, keep = h, np.ones(h.shape, bool)
    return {
        "n": int(kept.size),
        "n_rejected": int((~keep).sum()),
        "level_m": float(np.median(kept)),
        "mad_m": float(np.median(np.abs(kept - np.median(kept))) * 1.4826),
        "std_m": float(np.std(kept)),
    }


def adaptive_buffer(poly, buffer_m: float, min_area_frac: float = 0.65):
    """Shrink ``poly`` inward, but never past ``min_area_frac`` of its area.

    A fixed inward buffer is not scale-free: -150 m costs Alsek (75.6 km²) a
    tenth of its area and LGP (3.4 km²) a half, which is why LGP returned a usable
    level on 1 overpass of 43 while Alsek managed 18. The buffer is there to keep
    shorelines and calving fronts out of the sample, so it is backed off for small
    lakes rather than applied blindly.
    """
    if not buffer_m:
        return poly
    target = poly.area * min_area_frac

    full = poly.buffer(buffer_m)
    if not full.is_empty and full.area >= target:
        return full

    # Area shrinks monotonically as the buffer grows more negative, so bisect for
    # the most aggressive buffer that still leaves `target`. `ok` always satisfies
    # the constraint and `fail` never does.
    fail, ok = buffer_m, 0.0
    for _ in range(30):
        mid = (fail + ok) / 2.0
        g = poly.buffer(mid)
        if not g.is_empty and g.area >= target:
            ok = mid
        else:
            fail = mid
    out = poly.buffer(ok)
    log.info(
        "buffer eased from %.0f m to %.0f m to keep %.0f%% of the lake",
        buffer_m,
        ok,
        100 * min_area_frac,
    )
    return out if not out.is_empty else poly


def levels_from_granule(
    path: str | Path,
    lakes,
    lake_names: dict[int, str] | None = None,
    buffer_m: float = -150.0,
    mad_k: float = 3.0,
    max_scatter_m: float = 1.0,
) -> list[dict]:
    """One row per lake present in this granule.

    ``lakes`` is a GeoDataFrame with a ``class`` column. ``buffer_m`` is negative
    to shrink the polygons inward.
    """
    import geopandas as gpd

    names = lake_names or DEFAULT_LAKES
    path = Path(path)
    info = atl06.granule_info(path)

    minx, miny, maxx, maxy = lakes.to_crs(4326).total_bounds
    df = atl06.read_segments(path, bbox=(minx, miny, maxx, maxy), quality_only=True)
    if df.empty:
        raw = atl06.read_segments(path, bbox=(minx, miny, maxx, maxy), quality_only=False)
        log.warning("%s: %d segments in the box, none passed quality", path.name, len(raw))
        return [
            {
                "granule": path.name,
                "date": info.time_start[:10],
                "rgt": info.rgt,
                "cycle": info.cycle,
                "lake": name,
                "n": 0,
                "level_m": np.nan,
                "mad_m": np.nan,
                "std_m": np.nan,
                "n_rejected": 0,
                "flag": "no quality-passing segments",
            }
            for name in names.values()
        ]

    pts = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=4326).to_crs(
        lakes.crs
    )

    rows = []
    for code, name in names.items():
        sel = lakes[lakes["class"] == code]
        if sel.empty:
            continue
        poly = adaptive_buffer(sel.geometry.union_all(), buffer_m)
        inside = pts[pts.within(poly)]
        stats = robust_level(inside["h_li"].to_numpy(), mad_k=mad_k)
        flag = None
        if stats["n"] == 0:
            flag = "no points inside the outline"
        elif stats["n"] < 5:
            flag = f"only {stats['n']} points"
        elif stats["mad_m"] > max_scatter_m:
            flag = f"scatter {stats['mad_m']:.2f} m exceeds {max_scatter_m} m"
        geoid = (
            float(np.median(inside["geoid_h"]))
            if stats["n"] and np.isfinite(inside["geoid_h"]).any()
            else np.nan
        )
        rows.append(
            {
                "granule": path.name,
                "date": info.time_start[:10],
                "rgt": info.rgt,
                "cycle": info.cycle,
                "lake": name,
                **stats,
                "geoid_h_m": geoid,
                "level_egm2008_m": stats["level_m"] - geoid,
                "flag": flag,
            }
        )
    return rows


def build_timeseries(rows: list[dict]):
    """Tidy the per-granule rows into a sorted DataFrame."""
    import pandas as pd

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["date"] = pd.to_datetime(df["date"])
    df["year"] = df["date"].dt.year
    df["doy"] = df["date"].dt.dayofyear
    df["usable"] = df["flag"].isna()
    return df.sort_values(["lake", "date"]).reset_index(drop=True)
