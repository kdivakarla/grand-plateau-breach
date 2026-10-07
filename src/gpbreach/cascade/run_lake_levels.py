"""Build an ICESat-2 lake-elevation time series for LGP, GPL and Alsek.

    # see what exists, download nothing
    python -m gpbreach.cascade.run_lake_levels --list

    # melt season only, a few granules to start
    python -m gpbreach.cascade.run_lake_levels --months 5 6 7 8 9 --limit 6

    # the whole melt-season record (~41 granules, ~0.7 GB)
    python -m gpbreach.cascade.run_lake_levels --months 5 6 7 8 9 --yes

Granules are cached in ``data/raw/ICESat-2/ATL06/`` and never re-downloaded.
Writes ``lake_levels.csv``, ``lake_levels_summary.csv`` and ``lake_levels.png``.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import yaml

from ..config import repo_root
from ..io.nsidc import download_granule, search_granules
from .lake_levels import build_timeseries, levels_from_granule

log = logging.getLogger(__name__)

CACHE = repo_root() / "data" / "raw" / "ICESat-2" / "ATL06"
COLORS = {"LGP": "#2a78d6", "GPL": "#eb6834", "Alsek": "#1baf7a"}
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def _plot(df, out_png: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    good = df[df.usable]
    fig, axes = plt.subplots(
        2, 1, figsize=(11, 7.4), facecolor=SURFACE, gridspec_kw={"hspace": 0.33}
    )
    fig.subplots_adjust(left=0.085, right=0.985, top=0.86, bottom=0.09)

    ax = axes[0]
    for lake, sub in good.groupby("lake"):
        ax.errorbar(
            sub["date"],
            sub["level_m"],
            yerr=sub["mad_m"],
            fmt="o-",
            ms=5,
            lw=1.3,
            capsize=2.5,
            color=COLORS.get(lake, INK2),
            label=lake,
        )
    bad = df[~df.usable & df.n.gt(0)]
    if not bad.empty:
        ax.scatter(
            bad["date"],
            bad["level_m"],
            marker="x",
            s=28,
            c=INK2,
            alpha=0.55,
            label="flagged",
            zorder=1,
        )
    ax.set_ylabel("lake surface (m, WGS84 ellipsoidal)", color=INK2)
    ax.set_title("Full record", color=INK, fontsize=10.5, loc="left", pad=8)
    ax.legend(
        frameon=True, facecolor=SURFACE, edgecolor="#d8d7d2", fontsize=9, labelcolor=INK, ncol=4
    )

    ax = axes[1]
    for lake, sub in good.groupby("lake"):
        ax.scatter(
            sub["doy"], sub["level_m"], s=34, color=COLORS.get(lake, INK2), label=lake, alpha=0.85
        )
    ax.set_xlabel("day of year", color=INK2)
    ax.set_ylabel("lake surface (m)", color=INK2)
    ax.set_title("Seasonal view — every year overlaid", color=INK, fontsize=10.5, loc="left", pad=8)

    for a in axes:
        a.set_facecolor(SURFACE)
        a.tick_params(colors=INK2, labelsize=9)
        for s in a.spines.values():
            s.set_color("#d8d7d2")
        a.grid(alpha=0.25)

    fig.text(
        0.012,
        0.965,
        "ICESat-2 ATL06 lake surface elevations — Grand Plateau",
        color=INK,
        fontsize=13.5,
        ha="left",
        va="top",
    )
    fig.text(
        0.012,
        0.906,
        f"{len(good)} usable of {len(df)} lake-overpass pairs  ·  "
        f"error bars are robust scatter (MAD)  ·  x marks flagged epochs",
        color=INK2,
        fontsize=9.3,
        ha="left",
        va="top",
    )
    fig.savefig(out_png, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return out_png


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-config", default=str(repo_root() / "configs" / "cascade_data.yaml"))
    ap.add_argument(
        "--months",
        type=int,
        nargs="*",
        default=None,
        help="keep only these months, e.g. 5 6 7 8 9 for the melt season",
    )
    ap.add_argument("--start", default=None, help="ISO date")
    ap.add_argument("--end", default=None, help="ISO date")
    ap.add_argument("--limit", type=int, default=None, help="cap how many granules to use")
    ap.add_argument(
        "--buffer-m",
        type=float,
        default=-150.0,
        help="inward buffer on the lake outlines (default -150 m)",
    )
    ap.add_argument("--max-scatter-m", type=float, default=1.0)
    ap.add_argument("--list", action="store_true", help="search only, download nothing")
    ap.add_argument("--yes", action="store_true", help="skip the download confirmation")
    ap.add_argument("--cache", default=str(CACHE))
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import geopandas as gpd

    data = yaml.safe_load(Path(args.data_config).read_text())
    lakes = gpd.read_file(_resolve(data["outlines"]["lakes"]))
    names = {
        v: k.upper() if k != "alsek" else "Alsek"
        for k, v in data["class_codes"].items()
        if k in ("lgp", "gpl", "alsek")
    }
    bbox = tuple(lakes.to_crs(4326).total_bounds)
    print(f"lakes: {', '.join(names.values())}   bbox {tuple(round(v, 3) for v in bbox)}")

    granules = search_granules(
        "ATL06",
        bbox,
        start=args.start,
        end=args.end,
        months=tuple(args.months) if args.months else None,
    )
    if args.limit:
        granules = granules[-args.limit :]
    if not granules:
        raise SystemExit("no granules matched")
    total = sum(g.size_mb or 16.0 for g in granules)
    print(
        f"{len(granules)} granules, {granules[0].date} to {granules[-1].date}, "
        f"~{total / 1024:.2f} GB"
    )

    if args.list:
        for g in granules:
            print(f"  {g.date}  RGT {g.rgt:4d} cyc {g.cycle:2d}  {g.name}")
        return 0

    cache = Path(args.cache)
    have = sum(1 for g in granules if (cache / g.name).exists())
    todo = len(granules) - have
    if todo and not args.yes:
        gb = total * todo / len(granules) / 1024
        print(f"\n{have} already cached, {todo} to download (~{gb:.2f} GB)")
        if input("proceed? [y/N] ").strip().lower() not in ("y", "yes"):
            print("aborted; nothing downloaded")
            return 1

    rows = []
    for i, g in enumerate(granules, 1):
        try:
            path = download_granule(g, cache)
        except Exception as exc:
            log.warning("skipping %s: %s", g.name, exc)
            continue
        print(f"  [{i}/{len(granules)}] {g.date}  {g.name}", end="  ", flush=True)
        try:
            got = levels_from_granule(
                path, lakes, names, buffer_m=args.buffer_m, max_scatter_m=args.max_scatter_m
            )
        except Exception as exc:
            # One unreadable granule must not cost the other forty-two.
            print(f"-> FAILED: {type(exc).__name__}", flush=True)
            log.warning("%s could not be read: %s", g.name, exc)
            rows.append(
                {
                    "granule": g.name,
                    "date": g.date,
                    "rgt": g.rgt,
                    "cycle": g.cycle,
                    "lake": None,
                    "n": 0,
                    "level_m": np.nan,
                    "flag": f"unreadable: {type(exc).__name__}",
                }
            )
            continue
        rows.extend(got)
        ok = [r for r in got if r.get("flag") is None]
        print(f"-> {len(ok)} usable lake(s)" if ok else "-> none usable", flush=True)

    df = build_timeseries(rows)
    if df.empty:
        raise SystemExit("no rows produced")

    out = _resolve(data["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "lake_levels.csv", index=False)

    good = df[df.usable]
    print(f"\n{len(good)} usable of {len(df)} lake-overpass pairs")
    if not good.empty:
        summary = (
            good.groupby("lake")
            .agg(
                n_epochs=("level_m", "size"),
                first=("date", "min"),
                last=("date", "max"),
                median_m=("level_m", "median"),
                min_m=("level_m", "min"),
                max_m=("level_m", "max"),
                typical_scatter_m=("mad_m", "median"),
            )
            .reset_index()
        )
        summary["range_m"] = summary["max_m"] - summary["min_m"]
        summary.to_csv(out / "lake_levels_summary.csv", index=False)
        print(summary.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
        png = _plot(df, out / "lake_levels.png")
        print(f"\n  wrote lake_levels.csv, lake_levels_summary.csv, {png.name}")
    for r in df[~df.usable].itertuples():
        log.info("flagged %s %s: %s", r.date.date(), r.lake, r.flag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
