"""Command line interface.

    gpbreach preprocess configs/baseline_breach1_lgp.yaml
    gpbreach run        configs/baseline_breach1_lgp.yaml
    gpbreach seeds      configs/baseline_breach1_lgp.yaml --threshold 266.32
    gpbreach fcurve     configs/baseline_breach1_lgp.yaml
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

from .analysis.fcurve import build_fcurve
from .analysis.report import write_manifest, write_report
from .analysis.maps import breach_map
from .analysis.sensitivity import sensitivity_table
from .config import Config, repo_root
from .dam import DamObject
from .engine import run as run_engine
from .preprocess.build_dam import build_dam, find_basin_seeds


def _load_dam(cfg: Config) -> DamObject:
    if not cfg.dam_path.exists():
        raise SystemExit(
            f"{cfg.dam_path} does not exist. Run:\n"
            f"    gpbreach preprocess {cfg.path}"
        )
    return DamObject.load(cfg.dam_path)


def cmd_preprocess(args) -> int:
    cfg = Config.load(args.config)
    dam = build_dam(cfg)
    print(f"wrote {cfg.dam_path}")
    print(f"  grid           {dam.surface.shape}, {dam.cell_area_m2:.0f} m2 cells")
    print(f"  valid cells    {int(dam.valid.sum()):,}")
    print(f"  base year      {dam.base_year}")
    print(f"  vertical datum {dam.vertical_datum}")
    print(f"  seed (lake)    {dam.seed_lake} -> map {dam.meta['seed_lake_xy']}")
    print(f"  seed (target)  {dam.seed_target} -> map {dam.meta['seed_target_xy']}")
    return 0


def cmd_seeds(args) -> int:
    cfg = Config.load(args.config)
    dam = _load_dam(cfg)
    found = find_basin_seeds(
        dam.surface, dam.bed, f=cfg.nominal("flotation_fraction"),
        reference_threshold=args.threshold, n_basins=args.n_basins,
        connectivity=int(cfg.components.get("connectivity", 8)),
    )
    print(f"basins of Phi_f <= {args.threshold} m  (f = {cfg.nominal('flotation_fraction'):.6f})")
    for d in found:
        x, y = dam.xy(d["seed"])
        area = d["n_cells"] * dam.cell_area_m2 / 1e6
        print(f"  rank {d['rank']}: {area:8.3f} km2  seed {d['seed']}  "
              f"map ({x:.1f}, {y:.1f})  Phi_min {d['phi_min_m']:.2f} m")
    if args.csv:
        out = Path(args.csv)
        out.write_text(
            "rank,row,col,x,y,area_km2,phi_min_m\n"
            + "".join(
                f"{d['rank']},{d['seed'][0]},{d['seed'][1]},"
                f"{dam.xy(d['seed'])[0]:.3f},{dam.xy(d['seed'])[1]:.3f},"
                f"{d['n_cells'] * dam.cell_area_m2 / 1e6:.6f},{d['phi_min_m']:.4f}\n"
                for d in found
            )
        )
        print(f"wrote {out}  (load in QGIS as delimited text, CRS {dam.crs and 'as per dam'})")
    return 0


def cmd_run(args) -> int:
    cfg = Config.load(args.config)
    dam = _load_dam(cfg)
    result = run_engine(cfg, dam)
    run_dir = repo_root() / "results" / "runs" / cfg.run_id
    write_manifest(cfg, result, run_dir)
    report = write_report(cfg, result, run_dir)
    result.samples.to_csv(run_dir / "samples.csv")
    try:
        result.samples.to_parquet(run_dir / "samples.parquet")
    except (ImportError, ValueError):
        pass  # pyarrow not installed; the CSV is the fallback record

    s = result.samples.iloc[0]
    print(f"{cfg.run_id}")
    print(f"  Phi_crit                {s['phi_crit_m']:.4f} m")
    print(f"  breach (continuous)     {s['breach_year_continuous']:.3f}")
    print(f"  breach (annual maps)    {int(s['breach_year'])}")
    if len(result.samples) > 1:
        q = result.samples["breach_year_continuous"].quantile([0.05, 0.5, 0.95])
        print(f"  ensemble 5/50/95        {q.iloc[0]:.2f} / {q.iloc[1]:.2f} / {q.iloc[2]:.2f}")
    print(f"  report                  {report}")
    return 0


def cmd_fcurve(args) -> int:
    cfg = Config.load(args.config)
    dam = _load_dam(cfg)
    curve = build_fcurve(dam, args.f_min, args.f_max, args.spacing,
                         int(cfg.components.get("connectivity", 8)))
    out = repo_root() / "data" / "processed" / f"fcurve_{dam.name}.npz"
    curve.save(out)
    nominal = cfg.nominal("flotation_fraction")
    level, dhdt = cfg.nominal("lake_level_m"), cfg.nominal("dhdt_m_per_yr")
    print(f"wrote {out}  ({curve.f.size} nodes, f in [{args.f_min}, {args.f_max}])")
    print(f"  max interpolation error {curve.max_interp_error(dam):.4f} m")
    print("\n  f        r        Phi_crit(m)   breach year")
    for f_value in np.linspace(args.f_min, args.f_max, 10):
        phi = float(curve(f_value))
        yr = dam.base_year + (phi - level) / (-f_value * dhdt)
        mark = "  <- nominal" if abs(f_value - nominal) < args.spacing / 2 else ""
        print(f"  {f_value:.4f}  {(1 - f_value) / f_value:.4f}  {phi:10.3f}  {yr:11.2f}{mark}")
    return 0


def cmd_map(args) -> int:
    cfg = Config.load(args.config)
    dam = _load_dam(cfg)
    out = breach_map(
        dam,
        f=cfg.nominal("flotation_fraction"),
        level=cfg.nominal("lake_level_m"),
        dhdt=cfg.nominal("dhdt_m_per_yr"),
        connectivity=int(cfg.components.get("connectivity", 8)),
        out_dir=repo_root() / "results" / "runs" / cfg.run_id,
        label_source=cfg.breach.get("source_label", "source lake"),
        label_target=cfg.breach.get("target_label", "target lake"),
    )
    print(f"{cfg.run_id}: breach point {out['crest']} -> "
          f"map {dam.xy(out['crest'])}, Phi_crit {out['phi_crit']:.4f} m")
    for key in ("map", "points", "pathway"):
        print(f"  {key:9s} {out[key]}")
    print("  load the two CSVs in QGIS (delimited text, x/y fields, CRS = the raster's)")
    return 0


def cmd_sensitivity(args) -> int:
    cfg = Config.load(args.config)
    dam = _load_dam(cfg)
    table = sensitivity_table(
        dam,
        f=cfg.nominal("flotation_fraction"),
        level=cfg.nominal("lake_level_m"),
        dhdt=cfg.nominal("dhdt_m_per_yr"),
        connectivity=int(cfg.components.get("connectivity", 8)),
    )
    print(f"{cfg.run_id}: breach year sensitivity (central differences)\n")
    print(f"{'parameter':<20} {'+/- step':>9} {'year(-)':>9} {'year(+)':>9} "
          f"{'span (yr)':>10} {'per unit':>12}")
    for _, r in table.iterrows():
        print(f"{r['parameter']:<20} {r['half_step']:>9.3g} {r['year_minus']:>9.2f} "
              f"{r['year_plus']:>9.2f} {r['year_span_over_step']:>10.2f} "
              f"{r['d_year_per_unit']:>12.4f}")
    if args.csv:
        table.to_csv(args.csv, index=False)
        print(f"\nwrote {args.csv}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="gpbreach", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("preprocess", help="build the dam object from rasters")
    sp.add_argument("config")
    sp.set_defaults(func=cmd_preprocess)

    sp = sub.add_parser("seeds", help="list the major basins and their seed cells")
    sp.add_argument("config")
    sp.add_argument("--threshold", type=float, required=True, help="reference Phi threshold (m)")
    sp.add_argument("--n-basins", type=int, default=4)
    sp.add_argument("--csv", help="write the basin seeds to this CSV for QGIS")
    sp.set_defaults(func=cmd_seeds)

    sp = sub.add_parser("run", help="run the engine and write the report")
    sp.add_argument("config")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("fcurve", help="tabulate Phi_crit against flotation fraction")
    sp.add_argument("config")
    sp.add_argument("--f-min", type=float, default=0.80)
    sp.add_argument("--f-max", type=float, default=0.98)
    sp.add_argument("--spacing", type=float, default=0.005)
    sp.set_defaults(func=cmd_fcurve)

    sp = sub.add_parser("map", help="render the breach pathway and export it for QGIS")
    sp.add_argument("config")
    sp.set_defaults(func=cmd_map)

    sp = sub.add_parser("sensitivity", help="breach-year sensitivity to each parameter")
    sp.add_argument("config")
    sp.add_argument("--csv")
    sp.set_defaults(func=cmd_sensitivity)

    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
