"""Write projected ice-surface rasters for future years.

    python -m gpbreach.cascade.run_project_surface --years 2025 2030 2035
    python -m gpbreach.cascade.run_project_surface --years 2030 --model by_elevation
    python -m gpbreach.cascade.run_project_surface --years 2030 --model accelerating \\
        --dhdt0 -9.05 --accel -0.3

Outputs go to ``data/standard/`` and the matching ``surfaces:`` entries are
printed, ready to paste into ``configs/cascade_data.yaml``. Once pasted, every
other step runs on them unchanged.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import rasterio
import yaml

from ..config import repo_root
from .grid import TargetGrid, write
from .project_surface import project

log = logging.getLogger(__name__)

#: Loso 2009-2020 means per elevation band, GP North (D-026).
LOSO_BANDS = [(100.0, 200.0, -7.06), (200.0, 300.0, -9.38), (300.0, 400.0, -8.00)]


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--config", default=str(repo_root() / "configs" / "cascade_data.yaml"))
    ap.add_argument("--years", type=int, nargs="+", required=True)
    ap.add_argument(
        "--model", choices=("uniform", "by_elevation", "accelerating"), default="uniform"
    )
    ap.add_argument("--dhdt", type=float, default=-9.05, help="uniform rate, m/yr")
    ap.add_argument("--dhdt0", type=float, default=-9.05, help="accelerating: initial rate")
    ap.add_argument(
        "--accel", type=float, default=-0.3, help="accelerating: change in rate, m/yr per yr"
    )
    ap.add_argument("--base-surface", default=None)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    cfg = yaml.safe_load(Path(args.config).read_text())
    entry = (
        next(s for s in cfg["surfaces"] if s["id"] == args.base_surface)
        if args.base_surface
        else cfg["surfaces"][0]
    )
    base_year = int(str(entry.get("date", "2018"))[:4])

    with rasterio.open(_resolve(entry["path"])) as src:
        zs = src.read(1).astype(np.float64)
        if src.nodata is not None:
            zs[zs == src.nodata] = np.nan
        target = TargetGrid(src.transform, src.shape, src.crs)
    zs[np.abs(zs) > 1e30] = np.nan
    with rasterio.open(_resolve(cfg["classes"])) as src:
        classes = src.read(1)
    bed = None
    bid = cfg.get("default_bed")
    if bid:
        bp = next((b["path"] for b in cfg["beds"] if b["id"] == bid), None)
        if bp and _resolve(bp).exists():
            with rasterio.open(_resolve(bp)) as src:
                bed = src.read(1).astype(np.float64)
            bed[np.abs(bed) > 1e30] = np.nan

    ice = classes == cfg["class_codes"]["ice"]
    print(f"base surface {entry['id']} ({base_year}), {int(ice.sum()):,} ice cells")
    print(f"model: {args.model}")

    kw: dict = {"bed": bed}
    if args.model == "uniform":
        kw["dhdt"] = args.dhdt
        print(f"  uniform {args.dhdt:+.2f} m/yr")
    elif args.model == "by_elevation":
        kw["bands"] = LOSO_BANDS
        kw["default_dhdt"] = args.dhdt
        print("  per-elevation rates (Loso 2009-2020, D-026):")
        for lo, hi, r in LOSO_BANDS:
            print(f"    {lo:.0f}-{hi:.0f} m : {r:+.2f} m/yr")
        print(f"    elsewhere   : {args.dhdt:+.2f} m/yr")
    else:
        kw["dhdt0"], kw["accel"] = args.dhdt0, args.accel
        print(f"  dh/dt starts {args.dhdt0:+.2f} m/yr, changing {args.accel:+.2f} per yr")

    out = project(zs, ice, args.years, base_year, model=args.model, **kw)
    out_dir = repo_root() / "data" / "standard"
    entries = []
    print(f"\n{'year':>6s} {'mean ice surface':>18s} {'change':>9s}  file")
    base_mean = float(np.nanmean(zs[ice]))
    for yr, arr in sorted(out.items()):
        sid = f"proj_{yr}_{args.model}"
        path = out_dir / f"surf_{sid}.tif"
        write(path, arr.astype("float32"), target)
        m = float(np.nanmean(arr[ice]))
        print(f"{yr:6d} {m:18.2f} {m - base_mean:+9.2f}  {path.name}")
        entries.append(
            {
                "id": sid,
                "path": f"data/standard/{path.name}",
                "date": f"{yr}-09-03",
                "vertical": entry.get("vertical"),
                "source": f"projected from {entry['id']} by {args.model}",
            }
        )

    print("\npaste into configs/cascade_data.yaml under surfaces::\n")
    print(yaml.dump(entries, sort_keys=False, default_flow_style=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
