"""Build the extended analysis grid and regrid every cascade input onto it.

    python -m gpbreach.cascade.build_grid --buffer-m 500
    python -m gpbreach.cascade.build_grid --dry-run     # report the grid only

Reads from ``data/processed/`` and writes to ``data/standard/``; inputs are never
modified. Every run writes ``grid_manifest.json`` recording the target grid, the
source of each layer, its resampling method, and its coverage.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import rasterio

from ..config import repo_root
from .grid import (
    TargetGrid,
    aligned_grid,
    coverage_report,
    rasterize_classes,
    regrid,
    union_bounds,
    verify_against,
    write,
)

log = logging.getLogger(__name__)

IN_DIR = repo_root() / "data" / "processed"
OUT_DIR = repo_root() / "data" / "standard"

#: Reference grid: defines the CRS, the cell size and the alignment phase.
REFERENCE = "bed_ellip.tif"

#: layer name -> (source file, categorical?)
LAYERS = {
    "bed_ellip": ("bed_ellip.tif", False),
    "surf_2018": ("surf_2m.tif", False),
}

#: Extra named beds from the config, regridded onto the canonical grid. A bed may
#: arrive in any CRS or resolution -- IceBoost is EPSG:32608 at 100 m against the
#: project's ESRI:102247 at 5 m -- so this is the single place that is handled.
CONFIG = repo_root() / "configs" / "cascade_data.yaml"

#: Classes are rebuilt from the outlines rather than resampled -- see
#: grid.rasterize_classes for why resampling cannot extend them.
ICE_VECTOR = "ice_2018.gpkg"
LAKES_VECTOR = "lakes_2018.gpkg"
ORIGINAL_CLASSES = "classes.tif"

#: Vector layers whose extent the grid must contain.
EXTENT_VECTORS = ("lakes_2018.gpkg",)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--buffer-m",
        type=float,
        default=500.0,
        help="halo beyond the union of inputs (default 500 m)",
    )
    ap.add_argument(
        "--res", type=float, default=None, help="target cell size in m (default: the reference's)"
    )
    ap.add_argument("--dry-run", action="store_true", help="report the grid, write nothing")
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    ref_path = IN_DIR / REFERENCE
    if not ref_path.exists():
        raise SystemExit(f"reference grid {ref_path} not found")

    with rasterio.open(ref_path) as ref:
        bounds = union_bounds(
            raster_paths=[ref_path],
            vector_paths=[IN_DIR / v for v in EXTENT_VECTORS if (IN_DIR / v).exists()],
            crs=ref.crs,
            buffer_m=args.buffer_m,
        )
        target = aligned_grid(bounds, ref, res=args.res)
        ref_grid = TargetGrid(ref.transform, ref.shape, ref.crs)

    print("reference :", ref_grid.describe())
    print("target    :", target.describe())
    growth = target.n_cells / ref_grid.n_cells
    print(f"            {growth:.1f}x the reference cell count, halo {args.buffer_m:g} m")
    # Alignment is the property that makes the extension non-destructive.
    off_x = (target.transform.c - ref_grid.transform.c) / target.res[0]
    off_y = (target.transform.f - ref_grid.transform.f) / target.res[1]
    print(
        f"            phase offset {off_x:+.6f}, {off_y:+.6f} cells "
        f"({'ALIGNED' if max(abs(off_x % 1), abs(off_y % 1)) < 1e-9 else 'NOT ALIGNED'})"
    )
    if args.dry_run:
        return 0

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "created_utc": datetime.now(UTC).isoformat(),
        "reference": str(ref_path),
        "buffer_m": args.buffer_m,
        "target": {
            "width": target.shape[1],
            "height": target.shape[0],
            "res_m": target.res[0],
            "bounds": [round(v, 3) for v in target.bounds],
            "crs": target.crs.to_string(),
            "transform": list(target.transform)[:6],
        },
        "layers": {},
    }

    classes = None
    for name, (src_name, categorical) in LAYERS.items():
        src = IN_DIR / src_name
        if not src.exists():
            log.warning("skipping %s: %s not found", name, src)
            manifest["layers"][name] = {"status": "missing", "source": str(src)}
            continue
        log.info(
            "regridding %s from %s (%s)", name, src_name, "nearest" if categorical else "bilinear"
        )
        arr = regrid(src, target, categorical=categorical)
        dest = write(out_dir / f"{name}.tif", arr, target)
        if categorical:
            classes = arr
            valid = int((arr != 0).sum())
        else:
            valid = int(np.isfinite(arr).sum())
        manifest["layers"][name] = {
            "status": "ok",
            "source": str(src),
            "resampling": "nearest" if categorical else "bilinear",
            "output": str(dest),
            "valid_cells": valid,
            "valid_pct": round(100.0 * valid / target.n_cells, 2),
        }
        print(
            f"  wrote {dest.name:16s} {valid:>12,} valid cells "
            f"({100.0 * valid / target.n_cells:5.1f}% of grid)"
        )

    # --- classes: rasterise from the polygons, then check against the original ---
    ice_p, lakes_p = IN_DIR / ICE_VECTOR, IN_DIR / LAKES_VECTOR
    if ice_p.exists() or lakes_p.exists():
        log.info("rasterising classes from outlines (ice then lakes)")
        classes = rasterize_classes(
            target,
            ice_path=ice_p if ice_p.exists() else None,
            lakes_path=lakes_p if lakes_p.exists() else None,
        )
        dest = write(out_dir / "classes.tif", classes, target, nodata=0)
        manifest["layers"]["classes"] = {
            "status": "ok",
            "method": "rasterised from outlines",
            "sources": [str(ice_p), str(lakes_p)],
            "output": str(dest),
        }
        print(f"  wrote {dest.name:16s} rasterised from outlines")
        orig = IN_DIR / ORIGINAL_CLASSES
        if orig.exists():
            v = verify_against(classes, target, orig)
            manifest["classes_verification"] = v
            if v.get("status") == "compared":
                print(
                    f"    agreement with original classes.tif over its footprint: "
                    f"{v['agree_pct']:.3f}% of {v['overlap_cells']:,} cells"
                )
                if v["agree_pct"] < 99.0:
                    print(
                        "    WARNING: below 99% -- the rebuild may have changed the "
                        "classification, not just extended it. Inspect before use."
                    )

    # --- named beds from the config -----------------------------------------
    import yaml

    if CONFIG.exists():
        cfg = yaml.safe_load(CONFIG.read_text())
        for entry in cfg.get("beds", []):
            src = entry.get("source_path")
            if not src:
                continue
            src_p = repo_root() / src
            dest_name = Path(entry["path"]).name
            if not src_p.exists():
                log.warning("bed %s: %s not found", entry["id"], src_p)
                continue
            if (
                (out_dir / dest_name).exists()
                and src_p.samefile(IN_DIR / Path(src).name)
                and dest_name == f"{Path(src).stem}.tif"
                and entry["id"] == "millan"
            ):
                pass
            log.info("regridding bed %s from %s", entry["id"], src_p.name)
            arr = regrid(src_p, target, categorical=False)
            write(out_dir / dest_name, arr, target)
            n = int(np.isfinite(arr).sum())
            manifest["layers"][f"bed_{entry['id']}"] = {
                "status": "ok",
                "source": str(src_p),
                "output": str(out_dir / dest_name),
                "resampling": "bilinear",
                "valid_cells": n,
                "valid_pct": round(100.0 * n / target.n_cells, 2),
            }
            print(
                f"  wrote {dest_name:22s} {n:>12,} valid cells "
                f"({100.0 * n / target.n_cells:5.1f}% of grid)  [bed: {entry['id']}]"
            )

    if classes is not None:
        lakes = IN_DIR / "lakes_2018.gpkg"
        cov = coverage_report(classes, target, lakes if lakes.exists() else None)
        manifest["lake_coverage"] = cov
        print("\n  lake capture after extension:")
        for name, e in cov.items():
            if "captured_pct" in e:
                print(
                    f"    {name:6s} {e['in_grid_km2']:8.2f} km2 of "
                    f"{e['digitised_km2']:8.2f} km2  = {e['captured_pct']:5.1f}%"
                )
            else:
                print(f"    {name:6s} {e['in_grid_km2']:8.2f} km2")

    mpath = out_dir / "grid_manifest.json"
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"\n  wrote {mpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
