"""Co-register one hindcast event's rasters onto a common grid.

    python -m gpbreach.cascade.build_event_grid configs/events/knik.yaml
    python -m gpbreach.cascade.build_event_grid configs/events/*.yaml --res 20

The Grand Plateau ``build_grid`` is tied to that site's file layout. This is the
general version: it reads an event config, puts the bed and surface on one grid at
the event's own CRS and resolution, and reports what is still missing before the
pass search can run.

Beds and surfaces routinely disagree. Across the three events supplied, bed
rasters are IceBoost at 100 m while surfaces are 2 m DEMs, and Kaskawulsh's two
arrive in **different CRSs** (EPSG:3413 against EPSG:32607). Resolving that is
this module's whole job.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import rasterio
import yaml
from rasterio.warp import transform_bounds

from ..config import repo_root
from .grid import TargetGrid, regrid, write

log = logging.getLogger(__name__)


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def build(cfg: dict, res_override: float | None = None) -> dict:
    """Put every raster in one event on a shared grid. Returns a status report."""
    crs = rasterio.crs.CRS.from_string(cfg["grid"]["crs"])
    res = res_override or float(cfg["grid"]["res_m"])
    buf = float(cfg["grid"].get("buffer_m", 0))

    # Grid extent: where bed and surface actually overlap -- outside it the head
    # surface is undefined anyway, so a larger grid only adds empty cells.
    boxes = []
    for entry in cfg["beds"] + cfg["surfaces"]:
        p = _resolve(entry["source_path"])
        if not p.exists():
            continue
        with rasterio.open(p) as src:
            b = src.bounds if src.crs == crs else transform_bounds(src.crs, crs, *src.bounds)
            boxes.append(tuple(b))
    if len(boxes) < 2:
        raise SystemExit("need at least a bed and a surface")
    arr = np.array(boxes)
    left, bottom = arr[:, 0].max() - buf, arr[:, 1].max() - buf
    right, top = arr[:, 2].min() + buf, arr[:, 3].min() + buf
    if right <= left or top <= bottom:
        raise SystemExit("bed and surface do not overlap")

    width = int(np.ceil((right - left) / res))
    height = int(np.ceil((top - bottom) / res))
    target = TargetGrid(rasterio.Affine(res, 0, left, 0, -res, top), (height, width), crs)

    out_dir = _resolve(cfg["grid"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    report = {"event": cfg["event_id"], "grid": target.describe(), "layers": {}}

    for entry in cfg["beds"] + cfg["surfaces"]:
        src = _resolve(entry["source_path"])
        if not src.exists():
            report["layers"][entry["id"]] = "source missing"
            continue
        arr2 = regrid(src, target)
        dest = _resolve(entry["path"])
        dest.parent.mkdir(parents=True, exist_ok=True)
        write(dest, arr2, target)
        n = int(np.isfinite(arr2).sum())
        report["layers"][entry["id"]] = {
            "valid_cells": n,
            "valid_pct": round(100 * n / target.n_cells, 1),
            "min": round(float(np.nanmin(arr2)), 1) if n else None,
            "max": round(float(np.nanmax(arr2)), 1) if n else None,
            "output": str(dest),
        }

    # Thickness is the first thing that shows whether the two fit together.
    bed_p = _resolve(cfg["beds"][0]["path"])
    srf_p = _resolve(cfg["surfaces"][0]["path"])
    if bed_p.exists() and srf_p.exists():
        with rasterio.open(bed_p) as s1, rasterio.open(srf_p) as s2:
            zb, zs = s1.read(1).astype(float), s2.read(1).astype(float)
        zb[np.abs(zb) > 1e30] = np.nan
        zs[np.abs(zs) > 1e30] = np.nan
        H = zs - zb
        ok = np.isfinite(H)
        report["thickness"] = {
            "cells_with_both": int(ok.sum()),
            "median_m": round(float(np.nanmedian(H)), 1) if ok.any() else None,
            "max_m": round(float(np.nanmax(H)), 1) if ok.any() else None,
            "negative_cells": int(np.nansum(H < 0)),
            "negative_pct": round(100 * float(np.nansum(H < 0)) / max(int(ok.sum()), 1), 1),
        }
        write(out_dir / "thickness.tif", H.astype("float32"), target)

    missing = []
    if not cfg["observed"].get("date"):
        missing.append("observed breach date")
    if not (cfg["outlines"] or {}).get("lakes"):
        missing.append("lake outlines (polygons with a 'class' field)")
    if (cfg["lake_levels"] or {}).get("draining_lake") is None:
        missing.append("draining-lake level at the breach")
    if (cfg["lake_levels"] or {}).get("receiving_lake") is None:
        missing.append("receiving-lake level")
    report["blocking"] = missing
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("configs", nargs="+")
    ap.add_argument("--res", type=float, default=None)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")

    for cpath in args.configs:
        cfg = yaml.safe_load(Path(cpath).read_text())
        print("=" * 72)
        print(f"{cfg['event_id']}")
        try:
            rep = build(cfg, args.res)
        except SystemExit as exc:
            print(f"  FAILED: {exc}")
            continue
        print(f"  grid: {rep['grid']}")
        for name, info in rep["layers"].items():
            if isinstance(info, str):
                print(f"    {name:28s} {info}")
            else:
                print(
                    f"    {name:28s} {info['valid_pct']:5.1f}% valid  "
                    f"z {info['min']}..{info['max']}"
                )
        t = rep.get("thickness")
        if t:
            print(
                f"    thickness: median {t['median_m']} m, max {t['max_m']} m, "
                f"{t['negative_pct']}% negative ({t['negative_cells']:,} cells)"
            )
        if rep["blocking"]:
            print("  BLOCKING for calibration:")
            for m in rep["blocking"]:
                print(f"    - {m}")
        else:
            print("  ready for the pass search")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
