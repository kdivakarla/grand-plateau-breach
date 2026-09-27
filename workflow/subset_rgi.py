#!/usr/bin/env python
"""Cut a local subset out of the RGI 7.0 Alaska inventory.

The full region-01 shapefile is 27,509 glaciers and 162 MB, which is awkward to
pan around in QGIS. This writes just the Grand Plateau neighbourhood, reprojected
to the project CRS (D-009) so it drops straight into the existing QGIS project.

    python workflow/subset_rgi.py [--buffer-km 20]

Source data in data/raw/ is read-only and untouched; output goes to data/interim/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import geopandas as gpd
import rasterio
from pyproj import CRS, Transformer

ROOT = Path(__file__).resolve().parents[1]
SHP = ROOT / "data/raw/RGI7_Alaska/RGI2000-v7.0-G-01_alaska.shp"
REFERENCE_RASTER = ROOT / "../2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif"
OUT_DIR = ROOT / "data/interim"

#: The four glaciers in contact with the three lakes of the system.
#: Identified by combining size, terminus elevation and proximity: all four
#: terminate between 19 and 24 m, consistent with calving into lakes whose
#: surfaces sit near 17-25 m. Owner to confirm against imagery (D-012).
LAKE_TERMINATING = {
    "RGI2000-v7.0-G-01-27357": "Grand Plateau Glacier -> Grand Plateau Lake",
    "RGI2000-v7.0-G-01-17002": "Grand Plateau North (unnamed) -> GP / upper lake",
    "RGI2000-v7.0-G-01-16980": "Alsek Glacier -> Alsek Lake (NE)",
    "RGI2000-v7.0-G-01-16987": "unnamed 115 km2 -> Alsek Lake (NE)",
}
#: Kept for the earlier two-glacier product.
GRAND_PLATEAU = {
    k: v
    for k, v in LAKE_TERMINATING.items()
    if k in ("RGI2000-v7.0-G-01-17002", "RGI2000-v7.0-G-01-27357")
}


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--buffer-km", type=float, default=20.0, help="halo around the analysis domain (default 20)"
    )
    args = ap.parse_args()

    with rasterio.open(REFERENCE_RASTER) as src:
        albers, bounds = src.crs, src.bounds
    buf = args.buffer_km * 1000

    inv = Transformer.from_crs(albers, CRS.from_epsg(4326), always_xy=True)
    corners = [
        inv.transform(x, y)
        for x in (bounds.left - buf, bounds.right + buf)
        for y in (bounds.bottom - buf, bounds.top + buf)
    ]
    lons, lats = [c[0] for c in corners], [c[1] for c in corners]
    bbox = (min(lons), min(lats), max(lons), max(lats))

    g = gpd.read_file(SHP, bbox=bbox).to_crs(albers)

    # Flag the ones that actually overlap the modelled domain, not just the halo.
    domain = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries.from_wkt(
            [
                f"POLYGON(({bounds.left} {bounds.bottom},{bounds.right} {bounds.bottom},"
                f"{bounds.right} {bounds.top},{bounds.left} {bounds.top},"
                f"{bounds.left} {bounds.bottom}))"
            ]
        ),
        crs=albers,
    )
    g["in_domain"] = g.intersects(domain.geometry.iloc[0])
    g["dist_to_domain_km"] = g.distance(domain.geometry.iloc[0]) / 1000.0
    g["is_grand_plateau"] = g.rgi_id.isin(GRAND_PLATEAU)
    g["gp_label"] = g.rgi_id.map(GRAND_PLATEAU)
    g["is_lake_terminating"] = g.rgi_id.isin(LAKE_TERMINATING)
    g["lake_label"] = g.rgi_id.map(LAKE_TERMINATING)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"rgi7_grand_plateau_{int(args.buffer_km)}km.gpkg"
    g.to_file(out, layer="rgi7_glaciers", driver="GPKG")

    gp = g[g.is_grand_plateau]
    gp_out = OUT_DIR / "rgi7_grand_plateau_dams.gpkg"
    gp.to_file(gp_out, layer="dam_glaciers", driver="GPKG")

    four = g[g.is_lake_terminating]
    four_out = OUT_DIR / "rgi7_lake_terminating_4.gpkg"
    four.to_file(four_out, layer="lake_terminating", driver="GPKG")
    missing = set(LAKE_TERMINATING) - set(four.rgi_id)
    if missing:
        raise SystemExit(f"expected glaciers not found in the subset: {sorted(missing)}")

    print(f"CRS out   : {albers.to_string()}")
    print(f"halo      : {args.buffer_km:.0f} km around the analysis domain")
    print(f"wrote     : {out}   ({len(g)} glaciers)")
    print(f"            {gp_out}   ({len(gp)} dam glaciers)")
    print(f"            {four_out}   ({len(four)} lake-terminating)")
    print(f"  overlapping the domain itself : {int(g.in_domain.sum())}")
    print(f"  total area in subset          : {g.area_km2.sum():,.1f} km2")
    print("\n  the four lake-terminating glaciers:")
    for _, r in four.sort_values("area_km2", ascending=False).iterrows():
        print(
            f"    {r.rgi_id}  {r.area_km2:7.2f} km2  z {r.zmin_m:5.1f}-{r.zmax_m:.0f} m  "
            f"src {str(r.src_date)[:10]}  {r.lake_label}"
        )
    print("\n  RGI src_date is 2010 for all four, so these outlines predate the")
    print("  current terminus by ~16 yr. Trim with: python workflow/trim_termini.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
