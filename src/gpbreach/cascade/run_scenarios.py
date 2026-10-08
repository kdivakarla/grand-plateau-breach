"""Run a matrix of component swaps and compare what each one changes.

    # every bed against every thinning rate
    python -m gpbreach.cascade.run_scenarios --beds millan iceboost --thinning all

    # sensitivity to the flotation factor, one bed
    python -m gpbreach.cascade.run_scenarios --beds millan --k 0.9 0.95 1.0

    # list what is available to swap
    python -m gpbreach.cascade.run_scenarios --list

Writes ``scenario_comparison.csv`` and a markdown report. Geometry solves are
cached in ``data/interim/scenario_cache/``, so repeating a matrix with extra
thinning rates costs nothing.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import yaml

from ..config import repo_root
from .inflow_empirical import RELATIONS
from .scenarios import expand, run_scenario

log = logging.getLogger(__name__)
CFG = repo_root() / "configs" / "cascade_data.yaml"


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def _load_grids(cfg: dict, bed_id: str, surface_id: str) -> dict:
    bed = next(b for b in cfg["beds"] if b["id"] == bed_id)
    surf = next(s for s in cfg["surfaces"] if s["id"] == surface_id)

    def read(path):
        with rasterio.open(_resolve(path)) as src:
            a = src.read(1).astype(np.float64)
            if src.nodata is not None:
                a[a == src.nodata] = np.nan
            a[np.abs(a) > 1e30] = np.nan
            return a, src.transform, src.shape

    zb, transform, _ = read(bed["path"])
    zs, _, s2 = read(surf["path"])
    with rasterio.open(_resolve(cfg["classes"])) as src:
        classes = src.read(1)
        cell = abs(src.transform.a * src.transform.e)
    if not (zb.shape == zs.shape == classes.shape):
        raise SystemExit(
            f"grids disagree: bed {zb.shape}, surface {s2}, "
            f"classes {classes.shape} -- run build_grid"
        )
    return {
        "bed": zb,
        "surface": zs,
        "classes": classes,
        "transform": transform,
        "cell_area_m2": cell,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--config", default=str(CFG))
    ap.add_argument("--beds", nargs="*", default=None)
    ap.add_argument("--surfaces", nargs="*", default=None)
    ap.add_argument("--thinning", nargs="*", default=None)
    ap.add_argument("--k", nargs="*", type=float, default=[1.0])
    ap.add_argument("--connectivity", nargs="*", type=int, default=[8])
    ap.add_argument("--base-year", type=int, default=2018)
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    cfg = yaml.safe_load(Path(args.config).read_text())
    beds = {b["id"]: b for b in cfg["beds"]}
    surfaces = {s["id"]: s for s in cfg["surfaces"]}
    thinnings = {t["id"]: t for t in cfg["thinning"]}

    if args.list:
        print("swappable components in", args.config)
        for label, items in (("beds", beds), ("surfaces", surfaces), ("thinning", thinnings)):
            print(f"\n  {label}:")
            for k, v in items.items():
                extra = (
                    f"{v['dhdt_m_per_yr']:+g} m/yr"
                    if "dhdt_m_per_yr" in v
                    else Path(v.get("path", "")).name
                )
                print(f"    {k:24s} {extra}")
                print(f"      {v.get('source', '')}")
        return 0

    bed_ids = args.beds or [cfg.get("default_bed", next(iter(beds)))]
    surf_ids = args.surfaces or [next(iter(surfaces))]
    thin_ids = (
        list(thinnings) if args.thinning == ["all"] else args.thinning or [next(iter(thinnings))]
    )
    for name, chosen, avail in (
        ("bed", bed_ids, beds),
        ("surface", surf_ids, surfaces),
        ("thinning", thin_ids, thinnings),
    ):
        bad = set(chosen) - set(avail)
        if bad:
            raise SystemExit(f"unknown {name}(s): {sorted(bad)}. Available: {sorted(avail)}")

    codes = cfg["class_codes"]
    levels = cfg["lake_levels"]
    lake_levels = {codes[n]: v for n, v in levels.items() if n in codes and v is not None}

    scenarios = expand(bed_ids, surf_ids, thin_ids, tuple(args.k), tuple(args.connectivity))
    print(
        f"{len(scenarios)} scenario(s): {len(bed_ids)} bed x {len(surf_ids)} surface "
        f"x {len(thin_ids)} thinning x {len(args.k)} k\n"
    )

    results, grids_cache = [], {}
    for sc in scenarios:
        key = (sc.bed_id, sc.surface_id)
        if key not in grids_cache:
            grids_cache[key] = _load_grids(cfg, sc.bed_id, sc.surface_id)
        dhdt = thinnings[sc.thinning_id]["dhdt_m_per_yr"]
        r = run_scenario(
            sc,
            grids_cache[key],
            lake_levels,
            codes,
            dhdt,
            base_year=args.base_year,
            use_cache=not args.no_cache,
        )
        if r is None:
            print(f"  {sc.name}: NO ROUTE -- skipped")
            continue
        results.append(r)
        print(f"  {sc.name}")

    if not results:
        raise SystemExit("no scenario produced a result")

    df = pd.DataFrame([r.as_row() for r in results])
    out = _resolve(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "scenario_comparison.csv", index=False)

    print("\n" + "=" * 100)
    print("TIMING  — set by bed, surface and k through h_pass; then scaled by dh/dt")
    print(
        f"{'bed':10s} {'surface':22s} {'thin':13s} {'k':>5s} {'h_pass':>9s} "
        f"{'Delta':>9s} {'thin need':>10s} {'onset':>9s}"
    )
    for r in results:
        s = r.scenario
        print(
            f"{s.bed_id:10s} {s.surface_id:22s} {s.thinning_id:13s} {s.k:5.2f} "
            f"{r.h_pass_m:9.2f} {r.delta_m:+9.2f} {r.thinning_needed_m:10.1f} "
            f"{r.onset_year:9.1f}"
        )

    print("\n" + "=" * 100)
    print("FLOOD SIZE — the two readings of 'the sill' are not interchangeable")
    print(
        "  along-path : highest bed on the flotation route (conservative; assumes a fixed conduit)"
    )
    print("  bed-only   : lowest bed ridge on ANY route (what the bed raster shows)")
    print()
    print(
        f"{'bed':10s} {'sill def':11s} {'sill':>9s} {'h_final':>9s} {'drawdn':>8s} "
        f"{'controls':>16s} {'V_w Mm3':>15s} {'Qp tunnel m3/s':>17s}"
    )
    seen = set()
    for r in results:
        if r.scenario.bed_id in seen:
            continue
        seen.add(r.scenario.bed_id)
        rows = [
            (
                "along-path",
                r.bed_sill_m,
                r.h_final_m,
                r.drawdown_m,
                r.controlled_by,
                r.vw_lo_mm3,
                r.vw_hi_mm3,
            ),
            (
                "bed-only",
                r.bed_pass_m,
                r.h_final_bedpass_m,
                r.drawdown_bedpass_m,
                r.controlled_by_bedpass,
                r.vw_lo_bedpass_mm3,
                r.vw_hi_bedpass_mm3,
            ),
        ]
        for lbl, sill, hf, dd, ctrl, lo, hi in rows:
            tun = RELATIONS["walder_costa_tunnel"]
            qlo = tun.peak_discharge(max(lo, 0.0)) if np.isfinite(lo) else float("nan")
            qhi = tun.peak_discharge(max(hi, 0.0)) if np.isfinite(hi) else float("nan")
            print(
                f"{r.scenario.bed_id:10s} {lbl:11s} {sill:9.2f} {hf:9.2f} {dd:8.2f} "
                f"{ctrl:>16s} {lo:6.0f}-{hi:<8.0f} {qlo:7.0f}-{qhi:<9.0f}"
            )
        for n in r.notes:
            print(f"{'':10s}   note: {n}")
        print()
    if len(set(df.bed)) > 1:
        print("\n" + "=" * 100)
        print("WHAT THE BED SWAP MOVED")
        base = df[df.bed == bed_ids[0]].iloc[0]
        for bid in bed_ids[1:]:
            alt = df[df.bed == bid].iloc[0]
            print(f"\n  {bed_ids[0]} -> {bid}")
            for lbl, col, unit in (
                ("h_pass", "h_pass_m", "m"),
                ("thinning needed", "thinning_needed_m", "m"),
                ("onset year", "onset_year", "yr"),
                ("bed sill", "bed_sill_m", "m"),
                ("drawdown", "drawdown_m", "m"),
                ("V_w upper", "vw_hi_mm3", "Mm3"),
            ):
                d = alt[col] - base[col]
                print(f"    {lbl:18s} {base[col]:10.2f} -> {alt[col]:10.2f} ({d:+.2f} {unit})")

    print(f"\n  wrote {out / 'scenario_comparison.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
