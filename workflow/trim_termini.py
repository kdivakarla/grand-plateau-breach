#!/usr/bin/env python
"""Trim RGI outlines back to a current terminus position.

RGI 7.0 outlines for these four glaciers carry ``src_date = 2010-09-15``, so they
describe 2010 ice. Sixteen years of retreat later they overstate the glaciers at
the calving fronts. This cuts them back using a polygon **you** digitise from
current imagery.

Deciding where the current terminus is, is a physics/input judgement and belongs
to the owner (CLAUDE.md Section 2.2). This script does not detect a terminus and
does not guess one: it applies the cut you supply, reports how much area moved,
and leaves the original untouched.

Two ways to supply the cut, whichever is easier to draw:

  --mode difference  (default)  CUT = the deglaciated area: current lake and
                                bare ground where 2010 ice used to be.
                                Result = RGI outline minus CUT.

  --mode intersection           CUT = the current ice extent.
                                Result = RGI outline clipped to CUT.

Usage::

    python workflow/trim_termini.py --cut gis/current_lakes_2026.gpkg
    python workflow/trim_termini.py --cut gis/ice_extent_2026.gpkg --mode intersection

Output goes to ``data/interim/``; inputs are never modified.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IN = ROOT / "data/interim/rgi7_lake_terminating_4.gpkg"
OUT_DIR = ROOT / "data/interim"


def trim(
    glaciers: gpd.GeoDataFrame, cut: gpd.GeoDataFrame, mode: str = "difference"
) -> gpd.GeoDataFrame:
    """Return ``glaciers`` trimmed by ``cut``, reprojected to the glacier CRS."""
    if cut.crs != glaciers.crs:
        cut = cut.to_crs(glaciers.crs)
    merged = cut.geometry.union_all()

    out = glaciers.copy()
    out["area_before_km2"] = out.geometry.area / 1e6
    if mode == "difference":
        out["geometry"] = out.geometry.difference(merged)
    elif mode == "intersection":
        out["geometry"] = out.geometry.intersection(merged)
    else:
        raise ValueError(f"mode must be 'difference' or 'intersection', got {mode!r}")

    # A cut can shatter a tongue into slivers; keep only the main body but say so.
    out["n_parts"] = [
        len(g.geoms) if g.geom_type == "MultiPolygon" else (0 if g.is_empty else 1)
        for g in out.geometry
    ]
    out["area_after_km2"] = out.geometry.area / 1e6
    out["area_lost_km2"] = out.area_before_km2 - out.area_after_km2
    out["pct_lost"] = 100.0 * out.area_lost_km2 / out.area_before_km2
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--cut",
        required=True,
        help="polygon layer you digitised (lake/deglaciated, or current ice)",
    )
    ap.add_argument("--glaciers", default=str(DEFAULT_IN))
    ap.add_argument("--mode", choices=("difference", "intersection"), default="difference")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    g = gpd.read_file(args.glaciers)
    cut = gpd.read_file(args.cut)
    print(f"glaciers : {len(g)} features, CRS {g.crs.to_string()}")
    print(f"cut      : {len(cut)} features, CRS {cut.crs.to_string()}  mode={args.mode}\n")

    out = trim(g, cut, args.mode)

    print(f"{'rgi_id':26s} {'before':>9s} {'after':>9s} {'lost':>9s} {'%':>6s}  parts")
    for _, r in out.iterrows():
        print(
            f"{r.rgi_id:26s} {r.area_before_km2:9.2f} {r.area_after_km2:9.2f} "
            f"{r.area_lost_km2:9.2f} {r.pct_lost:6.2f}  {r.n_parts}"
        )
        if r.n_parts > 1:
            print("      ^ cut left multiple pieces -- inspect in QGIS before using")
        if r.pct_lost > 25:
            print("      ^ lost >25% of area; check the cut layer is what you meant")

    dest = Path(args.out) if args.out else OUT_DIR / "rgi7_lake_terminating_4_trimmed.gpkg"
    out.to_file(dest, layer="trimmed", driver="GPKG")
    print(f"\nwrote {dest}")
    print("Original outlines are unchanged; compare the two layers in QGIS.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
