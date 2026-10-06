"""Stage curves for LGP and GPL, and the released volume V_w.

    python -m gpbreach.cascade.run_lakes
    python -m gpbreach.cascade.run_lakes --lgp-max-depth 120

LGP bathymetry is not available, so V_w is reported as a **bracket** between a
vertical-wall upper bound and a conic lower bound rather than as a single number.
Supply soundings to replace both.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import rasterio
import yaml

from ..config import repo_root
from .lakes import (
    conic_bottom,
    connected_stage_area,
    detect_spill,
    released_volume,
    stage_curve,
    vertical_wall_bottom,
)

log = logging.getLogger(__name__)


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-config", default=str(repo_root() / "configs" / "cascade_data.yaml"))
    ap.add_argument(
        "--bed-sill",
        type=float,
        default=None,
        help="m; default reads path_summary from the last run_path",
    )
    ap.add_argument(
        "--lgp-max-depth",
        type=float,
        default=None,
        help="m; deepest point for the conic lower bound "
        "(default: drawdown, i.e. just deep enough to drain)",
    )
    ap.add_argument("--gpl-rise-max", type=float, default=30.0)
    ap.add_argument("--step", type=float, default=0.5, help="stage step, m")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    data = yaml.safe_load(Path(args.data_config).read_text())
    codes, levels = data["class_codes"], data["lake_levels"]
    h_lgp, h_gpl = levels["lgp"], levels["gpl"]
    out_dir = _resolve(data["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    sill = args.bed_sill
    if sill is None:
        cands = sorted(out_dir.glob("path_summary_*.json"))
        if not cands:
            raise SystemExit("no path_summary_*.json; run run_path or pass --bed-sill")
        sill = json.loads(cands[-1].read_text())["summary"]["bed_sill_m"]
        print(f"bed sill {sill:.2f} m  (from {cands[-1].name})")

    with rasterio.open(_resolve(data["classes"])) as src:
        classes = src.read(1)
        cell = abs(src.transform.a * src.transform.e)
    surf_path = _resolve(data["surfaces"][0]["path"])
    with rasterio.open(surf_path) as src:
        dem = src.read(1).astype(np.float64)
        if src.nodata is not None:
            dem[dem == src.nodata] = np.nan
    dem[np.abs(dem) > 1e30] = np.nan

    lgp = classes == codes["lgp"]
    gpl = classes == codes["gpl"]
    a_lgp = lgp.sum() * cell / 1e6
    a_gpl = gpl.sum() * cell / 1e6
    print(f"LGP {a_lgp:.3f} km2   GPL {a_gpl:.3f} km2   cell {cell:.0f} m2")

    h_final = max(h_gpl, sill)
    drawdown = h_lgp - h_final
    print(
        f"\nLGP drains {h_lgp:.2f} -> {h_final:.2f} m "
        f"({'bed sill' if sill >= h_gpl else 'GPL'} controls), drawdown {drawdown:.2f} m"
    )

    # --- LGP: bracket, because there is no bathymetry ---------------------
    max_depth = args.lgp_max_depth if args.lgp_max_depth is not None else drawdown
    stages = np.arange(h_lgp - max_depth - 1.0, h_lgp + args.step, args.step)
    bounds = {}
    for name, bottom in (
        ("upper_vertical_walls", vertical_wall_bottom(lgp, h_lgp - max_depth)),
        ("lower_conic", conic_bottom(lgp, h_lgp, max_depth)),
    ):
        curve = stage_curve(bottom, lgp, stages, cell)
        curve.to_csv(out_dir / f"lgp_hypsometry_{name}.csv", index=False)
        bounds[name] = released_volume(curve, h_lgp, h_gpl, sill)

    print("\nreleased volume V_w -- a BRACKET, not an estimate (no LGP bathymetry):")
    for name, r in bounds.items():
        print(f"  {name:22s} {r['Vw_Mm3']:8.1f} Mm3 = {r['Vw_km3']:.3f} km3")
    lo = bounds["lower_conic"]["Vw_Mm3"]
    hi = bounds["upper_vertical_walls"]["Vw_Mm3"]
    print(f"  -> V_w between {lo:.0f} and {hi:.0f} Mm3 ({lo / 1e3:.2f}-{hi / 1e3:.2f} km3)")
    w = bounds["upper_vertical_walls"]["sill_withholds_Mm3"]
    print(
        f"  water the sill withholds below {sill:.1f} m: "
        + (f"{w:.0f} Mm3" if w is not None else "unknown without bathymetry")
    )

    spec = 500.0
    print(
        f"\n  spec expectation ~{spec:.0f} Mm3. Upper bound is "
        f"{hi / spec:.2f}x that, so the expectation looks high -- "
        f"{spec * 1e6 / (drawdown * 1e6):.2f} km2 of lake would be needed "
        f"against {a_lgp:.2f} km2 digitised."
    )

    # --- GPL: how far it rises taking the flood ---------------------------
    gpl_stages = np.arange(h_gpl, h_gpl + args.gpl_rise_max + args.step, args.step)
    gpl_curve = connected_stage_area(dem, gpl, gpl_stages, cell, base_level=h_gpl)
    gpl_curve.to_csv(out_dir / "gpl_hypsometry.csv", index=False)
    spill = detect_spill(gpl_curve)
    print("\nGPL response:")
    if spill is None:
        for v in (lo, hi):
            rise = float(
                np.interp(v * 1e6, gpl_curve["V_m3"].to_numpy(), gpl_curve["h"].to_numpy())
            )
            print(f"  rises to {rise:.2f} m (+{rise - h_gpl:.2f} m) taking {v:.0f} Mm3")
    else:
        cap = spill["capacity_before_spill_Mm3"]
        print(
            f"  own basin holds {cap:.0f} Mm3, rising to {spill['stage_before_m']:.2f} m "
            f"(+{spill['stage_before_m'] - h_gpl:.2f} m)"
        )
        print(
            f"  at {spill['spill_stage_m']:.2f} m it SPILLS into an adjoining basin: "
            f"area {spill['area_before_km2']:.1f} -> {spill['area_after_km2']:.1f} km2, "
            f"absorbing a further {spill['volume_absorbed_by_step_Mm3']:.0f} Mm3"
        )
        for label, v in (("lower", lo), ("upper", hi)):
            verdict = (
                "fits in GPL's own basin"
                if v <= cap
                else "EXCEEDS it -- water spills into the adjoining basin"
            )
            print(f"  V_w {label} bound {v:.0f} Mm3: {verdict}")
        print(
            "  -> do not read a stage off the curve across the step; it is a "
            "discontinuity, not a slope"
        )

    summary = {
        "created_utc": datetime.now(UTC).isoformat(),
        "lake_levels": {"lgp": h_lgp, "gpl": h_gpl},
        "areas_km2": {"lgp": a_lgp, "gpl": a_gpl},
        "bed_sill_m": sill,
        "h_final_m": h_final,
        "drawdown_m": drawdown,
        "bathymetry": "ABSENT -- V_w is bracketed by prism and cone bounds",
        "Vw_bounds": bounds,
        "Vw_Mm3_range": [lo, hi],
        "gpl_spill": spill,
        "spec_expectation_Mm3": spec,
    }
    (out_dir / "volume_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("\n  wrote lgp_hypsometry_*.csv, gpl_hypsometry.csv, volume_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
