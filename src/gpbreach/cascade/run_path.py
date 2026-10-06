"""Trace the LGP→GPL drainage path and write its profile.

    python -m gpbreach.cascade.run_path
    python -m gpbreach.cascade.run_path --n-mode mean_seal --no-refine

Consumes the head surface from ``run_seal`` (recomputed here so the parent array
is available) and writes ``path_<id>.gpkg``, ``profile_<id>.csv``,
``path_summary_<id>.json`` and ``profile_<id>.png``.
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
from .path import (
    monotonicity,
    profile,
    refine_downstream,
    summarize,
    trace_path,
    write_path_gpkg,
)
from .seal import head_surface, masked_smooth, pass_search

log = logging.getLogger(__name__)

# Categorical slots 1-3 of the reference palette: these three clear the
# all-pairs CVD and normal-vision floors on a light surface.
BED, SURF, HEAD = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def _plot(prof, summary, h_lgp, h_gpl, out_png: Path, title: str) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = prof["s"].to_numpy() / 1000.0
    fig, ax = plt.subplots(figsize=(10, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    ax.fill_between(s, prof["z_b"], prof["z_s"], color=SURF, alpha=0.13, lw=0, label="ice column")
    ax.plot(s, prof["z_s"], color=SURF, lw=2, label="ice surface")
    ax.plot(s, prof["z_b"], color=BED, lw=2, label="bed")
    ax.plot(s, prof["head"], color=HEAD, lw=2, label="hydraulic head")
    ax.axhline(h_lgp, color=INK, lw=1.3, ls="--")
    ax.axhline(h_gpl, color=INK2, lw=1.1, ls=":")
    ax.annotate(
        f"LGP {h_lgp:g} m",
        (s[-1], h_lgp),
        xytext=(-4, 5),
        textcoords="offset points",
        ha="right",
        color=INK,
        fontsize=9,
    )
    ax.annotate(
        f"GPL {h_gpl:g} m",
        (s[-1], h_gpl),
        xytext=(-4, -13),
        textcoords="offset points",
        ha="right",
        color=INK2,
        fontsize=9,
    )

    sp = summary.l_total_m - summary.l_pass_m
    ax.plot([sp / 1000.0], [summary.h_at_pass_m], "o", ms=11, mfc="none", mec=INK, mew=2, zorder=6)
    ax.annotate(
        f"pass  {summary.h_at_pass_m:.1f} m\nN = {summary.n_at_pass_m:+.1f} m",
        (sp / 1000.0, summary.h_at_pass_m),
        xytext=(-14, -46),
        textcoords="offset points",
        ha="right",
        color=INK,
        fontsize=9.5,
        fontweight="bold",
        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8),
        bbox=dict(boxstyle="round,pad=0.3", fc=SURFACE, ec=INK2, lw=0.6, alpha=0.92),
    )
    ax.plot([summary.bed_sill_s_m / 1000.0], [summary.bed_sill_m], "v", ms=9, color=BED, zorder=6)
    ax.annotate(
        f"bed sill {summary.bed_sill_m:.1f} m",
        (summary.bed_sill_s_m / 1000.0, summary.bed_sill_m),
        xytext=(6, -16),
        textcoords="offset points",
        color=BED,
        fontsize=9,
    )

    ax.set_xlabel("distance from LGP along the drainage path (km)", color=INK2)
    ax.set_ylabel("elevation / head (m, WGS84 ellipsoidal)", color=INK2)
    ax.set_title(title, color=INK, fontsize=12.5, loc="left", pad=20)
    ax.text(
        0.0,
        1.015,
        f"L_pass {summary.l_pass_m / 1000:.2f} km  ·  "
        f"bed sill {summary.bed_sill_m:.1f} m caps how far LGP can drain  ·  "
        f"flotation zone {summary.flotation_length_m / 1000:.2f} km",
        transform=ax.transAxes,
        color=INK2,
        fontsize=9.5,
        va="bottom",
    )
    ax.tick_params(colors=INK2, labelsize=9)
    for sp_ in ax.spines.values():
        sp_.set_color("#d8d7d2")
    ax.grid(alpha=0.25)
    ax.legend(
        frameon=True,
        facecolor=SURFACE,
        edgecolor="#d8d7d2",
        fontsize=9,
        labelcolor=INK,
        loc="upper right",
        ncol=1,
        framealpha=0.93,
    )
    ax.margins(x=0.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return out_png


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-config", default=str(repo_root() / "configs" / "cascade_data.yaml"))
    ap.add_argument("--params-config", default=str(repo_root() / "configs" / "cascade_params.yaml"))
    ap.add_argument("--surface-id", default=None)
    ap.add_argument("--n-mode", default="at_pass", choices=("at_pass", "mean_seal", "min"))
    ap.add_argument(
        "--no-refine", action="store_true", help="keep the raw minimax path below the pass"
    )
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    data = yaml.safe_load(Path(args.data_config).read_text())
    params = yaml.safe_load(Path(args.params_config).read_text())
    entry = (
        next(s for s in data["surfaces"] if s["id"] == args.surface_id)
        if args.surface_id
        else data["surfaces"][0]
    )
    sid = entry["id"]
    codes = data["class_codes"]
    levels = data["lake_levels"]
    h_lgp, h_gpl = levels["lgp"], levels["gpl"]
    lake_levels = {codes[n]: v for n, v in levels.items() if n in codes}

    with rasterio.open(_resolve(data["bed"]["path"])) as src:
        zb = src.read(1).astype(np.float64)
        transform, crs = src.transform, src.crs
        if src.nodata is not None:
            zb[zb == src.nodata] = np.nan
    with rasterio.open(_resolve(entry["path"])) as src:
        zs = src.read(1).astype(np.float64)
        if src.nodata is not None:
            zs[zs == src.nodata] = np.nan
    with rasterio.open(_resolve(data["classes"])) as src:
        classes = src.read(1)
    zb[np.abs(zb) > 1e30] = np.nan
    zs[np.abs(zs) > 1e30] = np.nan

    sigma = params["smoothing_m"]
    if sigma and sigma > 0:
        zs = masked_smooth(zs, classes == codes["ice"], sigma, abs(transform.a))

    head, _ = head_surface(zb, zs, classes, lake_levels, k=params["k"])
    res = pass_search(
        head, classes, params["source_class"], params["target_class"], params["connectivity"]
    )
    if not res.reached_target:
        raise SystemExit("no route from GPL to LGP; nothing to trace")

    p = trace_path(res.parent, res.info["best_target_cell"])
    n_raw = len(p)
    if not args.no_refine:
        p = refine_downstream(
            p, head, classes, res.pass_cell, target_class=codes["gpl"], parent=res.parent
        )
    prof = profile(p, zb, zs, head, transform, h_lgp, h_gpl)
    summ = summarize(
        prof, classes, h_lgp, res.pass_cell, lgp_class=codes["lgp"], n_mode=args.n_mode
    )

    refined = (not args.no_refine) and len(p) != n_raw
    print(
        f"path: {n_raw} cells traced"
        + (
            f", {len(p)} after downstream refinement"
            if refined
            else " (unrefined minimax route)"
            if args.no_refine
            else " -- refinement did not reach GPL, keeping the minimax route"
        )
    )

    print(f"\n  total length        {summ.l_total_m / 1000:8.3f} km")
    print(f"  LGP margin -> GPL   {summ.l_margin_m / 1000:8.3f} km   (L_margin)")
    print(f"  pass -> GPL         {summ.l_pass_m / 1000:8.3f} km   (L_pass)")
    print(f"\n  bed sill            {summ.bed_sill_m:8.2f} m  at {summ.bed_sill_s_m / 1000:.3f} km")
    print(
        f"    -> LGP cannot drain below max(h_GPL, bed sill) = {max(h_gpl, summ.bed_sill_m):.2f} m"
    )
    print(f"  head at pass        {summ.h_at_pass_m:8.2f} m")
    print(
        f"  thickness at pass   {summ.thickness_at_pass_m:8.2f} m   "
        f"(min {summ.thickness_min_m:.1f}, mean {summ.thickness_mean_m:.1f})"
    )
    print(
        f"  flotation-zone len  {summ.flotation_length_m / 1000:8.3f} km "
        f"(head < h_LGP along the path)"
    )
    print(f"\n  N at pass           {summ.n_at_pass_m:+8.2f} m")
    print(f"  N mean over seal    {summ.n_mean_seal_m:+8.2f} m")
    print(f"  N min along path    {summ.n_min_m:+8.2f} m")
    print(f"  -> representative N {summ.n_representative_m:+8.2f} m  ({args.n_mode})")

    below = prof.iloc[
        int(
            np.flatnonzero(
                (prof["row"].to_numpy() == res.pass_cell[0])
                & (prof["col"].to_numpy() == res.pass_cell[1])
            )[0]
        ) :
    ]
    mono = monotonicity(below, "head")
    print(
        f"\n  head below the pass: {'monotonic' if mono['monotonic'] else 'NOT monotonic'} "
        f"-- rises at {mono['rising_steps']}/{mono['steps']} steps, "
        f"largest {mono['largest_rise_m']:+.2f} m, total {mono['total_rise_m']:.1f} m"
    )
    if not mono["monotonic"]:
        print("    (surface roughness; this is what smoothing_m is for -- currently 0)")

    out = _resolve(data["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    prof.to_csv(out / f"profile_{sid}.csv", index=False)
    attrs = summ.as_dict() | {
        "surface_id": sid,
        "k": params["k"],
        "connectivity": params["connectivity"],
    }
    write_path_gpkg(prof, crs, out / f"path_{sid}.gpkg", attrs)
    png = _plot(
        prof, summ, h_lgp, h_gpl, out / f"profile_{sid}.png", f"LGP to GPL drainage path — {sid}"
    )
    (out / f"path_summary_{sid}.json").write_text(
        json.dumps(
            {
                "created_utc": datetime.now(UTC).isoformat(),
                "surface": entry,
                "n_mode": args.n_mode,
                "refined_downstream": not args.no_refine,
                "n_cells": len(p),
                "summary": summ.as_dict(),
                "head_monotonicity_below_pass": mono,
                "pass_cell": list(res.pass_cell),
                "h_pass_m": res.h_pass,
            },
            indent=2,
            default=str,
        )
    )
    print(f"\n  wrote profile_{sid}.csv, path_{sid}.gpkg, path_summary_{sid}.json, {png.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
