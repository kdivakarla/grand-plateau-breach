"""Run the head surface and pass search for one surface epoch.

    python -m gpbreach.cascade.run_seal
    python -m gpbreach.cascade.run_seal --k 0.95 --connectivity 4

Writes ``head_<id>.tif``, ``margin_lgp_<id>.tif``, ``spill_<id>.tif`` and
``seal_summary_<id>.json`` to the configured output directory.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import rasterio
import yaml

from ..config import repo_root
from .grid import TargetGrid, write
from .seal import competing_outlets, head_surface, masked_smooth, overtopping_search, pass_search

log = logging.getLogger(__name__)

DATA_CFG = repo_root() / "configs" / "cascade_data.yaml"
PARAMS_CFG = repo_root() / "configs" / "cascade_params.yaml"


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def _git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "-C", str(repo_root()), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-config", default=str(DATA_CFG))
    ap.add_argument("--params-config", default=str(PARAMS_CFG))
    ap.add_argument("--surface-id", default=None, help="which surfaces entry to use")
    ap.add_argument("--k", type=float, default=None, help="override the flotation factor")
    ap.add_argument("--connectivity", type=int, default=None, choices=(4, 8))
    ap.add_argument("--smoothing-m", type=float, default=None)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    data = yaml.safe_load(Path(args.data_config).read_text())
    params = yaml.safe_load(Path(args.params_config).read_text())
    k = args.k if args.k is not None else params["k"]
    conn = args.connectivity if args.connectivity is not None else params["connectivity"]
    sigma = args.smoothing_m if args.smoothing_m is not None else params["smoothing_m"]

    surfaces = data["surfaces"]
    entry = (
        next(s for s in surfaces if s["id"] == args.surface_id) if args.surface_id else surfaces[0]
    )
    sid = entry["id"]

    codes = data["class_codes"]
    levels_by_name = data["lake_levels"]
    lake_levels = {codes[n]: v for n, v in levels_by_name.items() if n in codes}

    with rasterio.open(_resolve(data["bed"]["path"])) as src:
        zb = src.read(1).astype(np.float64)
        target = TargetGrid(src.transform, src.shape, src.crs)
        if src.nodata is not None:
            zb[zb == src.nodata] = np.nan
    with rasterio.open(_resolve(entry["path"])) as src:
        zs = src.read(1).astype(np.float64)
        if src.nodata is not None:
            zs[zs == src.nodata] = np.nan
        if src.shape != target.shape:
            raise SystemExit(
                f"surface {src.shape} is not on the bed grid {target.shape}; "
                "run gpbreach.cascade.build_grid first"
            )
    with rasterio.open(_resolve(data["classes"])) as src:
        classes = src.read(1)
        if src.shape != target.shape:
            raise SystemExit("classes raster is not on the bed grid")
    zb[np.abs(zb) > 1e30] = np.nan
    zs[np.abs(zs) > 1e30] = np.nan

    print(f"surface   : {sid}  ({entry.get('date', '?')})")
    print(f"grid      : {target.describe()}")
    print(
        f"k = {k}  (f = kappa = {k * 917 / 1000:.6f})   connectivity {conn}   smoothing {sigma} m"
    )

    if sigma and sigma > 0:
        zs = masked_smooth(zs, classes == codes["ice"], sigma, target.res[0])

    head, info = head_surface(zb, zs, classes, lake_levels, k=k)
    print(
        f"\nhead: {info['ice_cells_with_data']:,} ice cells with data of "
        f"{info['ice_cells']:,} ({100 * info['ice_cells_with_data'] / info['ice_cells']:.1f}%); "
        f"{info['negative_thickness_cells']:,} negative thickness"
    )

    res = pass_search(head, classes, params["source_class"], params["target_class"], conn)
    h_lgp = levels_by_name.get("lgp")

    print(f"\n{'=' * 64}")
    if not res.reached_target:
        print("  NO ROUTE from GPL to LGP -- the seal cannot open on this geometry")
    else:
        print(f"  h_pass        {res.h_pass:10.3f} m (ellipsoidal, metres of water)")
        if h_lgp is not None:
            delta = res.h_pass - h_lgp
            thin = delta / (k * 917 / 1000)
            print(f"  h_LGP         {h_lgp:10.3f} m")
            print(
                f"  Delta         {delta:+10.3f} m   "
                f"({'SEALED' if delta > 0 else 'OPEN -- would already drain'})"
            )
            print(f"  thinning needed at the pass {thin:8.2f} m")
        if res.pass_cell:
            x, y = target.transform * (res.pass_cell[1] + 0.5, res.pass_cell[0] + 0.5)
            print(f"  pass cell     {res.pass_cell}  map ({x:.1f}, {y:.1f})")
    print(f"{'=' * 64}")

    over = overtopping_search(
        zs, classes, lake_levels, params["source_class"], params["target_class"], conn
    )
    if over.reached_target:
        print(
            f"  overtopping route would need LGP at {over.h_pass:.1f} m "
            f"(vs {res.h_pass:.1f} m for flotation)"
        )

    comp = competing_outlets(
        head, classes, {"gpl": codes["gpl"], "alsek": codes["alsek"]}, params["target_class"], conn
    )
    print("  competing outlets:")
    for name, e in comp.items():
        if e["h_pass"] is None:
            print(f"    {name:6s} {e['status']}")
        else:
            print(f"    {name:6s} h_pass {e['h_pass']:.3f} m")

    out_dir = _resolve(data["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    write(
        out_dir / f"head_{sid}.tif",
        np.where(np.isinf(head), np.nan, head).astype("float32"),
        target,
    )
    spill_out = np.where(np.isinf(res.spill), np.nan, res.spill).astype("float32")
    write(out_dir / f"spill_{sid}.tif", spill_out, target)
    if h_lgp is not None:
        margin = np.where((classes == codes["ice"]) & np.isfinite(head), head - h_lgp, np.nan)
        write(out_dir / f"margin_lgp_{sid}.tif", margin.astype("float32"), target)

    summary = {
        "run": {
            "created_utc": datetime.now(UTC).isoformat(),
            "git_commit": _git_commit(),
            "data_config": str(args.data_config),
            "params_config": str(args.params_config),
        },
        "surface": entry,
        "params": {"k": k, "kappa": k * 917 / 1000, "connectivity": conn, "smoothing_m": sigma},
        "lake_levels": levels_by_name,
        "head_info": info,
        "pass": {
            "reached_target": bool(res.reached_target),
            "h_pass_m": None if not np.isfinite(res.h_pass) else res.h_pass,
            "delta_m": None
            if (h_lgp is None or not np.isfinite(res.h_pass))
            else res.h_pass - h_lgp,
            "thinning_needed_m": None
            if (h_lgp is None or not np.isfinite(res.h_pass))
            else (res.h_pass - h_lgp) / (k * 917 / 1000),
            "pass_cell": list(res.pass_cell) if res.pass_cell else None,
            "search": res.info,
        },
        "overtopping": {"h_pass_m": over.h_pass if over.reached_target else None},
        "competing_outlets": comp,
    }
    spath = out_dir / f"seal_summary_{sid}.json"
    spath.write_text(json.dumps(summary, indent=2, default=str))
    print(
        f"\n  wrote {out_dir}/head_{sid}.tif, spill_{sid}.tif, "
        f"margin_lgp_{sid}.tif\n  wrote {spath}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
