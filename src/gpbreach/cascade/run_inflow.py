# ruff: noqa: RUF001
"""Empirical inflow hydrographs for the LGP→GPL flood.

    python -m gpbreach.cascade.run_inflow
    python -m gpbreach.cascade.run_inflow --vw-mm3 200

Reads V_w from ``volume_summary.json`` unless given. Because V_w is a bracket
(no LGP bathymetry, D-015), each relation is run at both ends and the result is
reported as a range.
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from ..config import repo_root
from .inflow_empirical import RELATIONS, hydrograph, summarise

log = logging.getLogger(__name__)

# Categorical slots 1-3; the reference relation is deliberately NOT a fourth hue.
COLORS = {
    "walder_costa_tunnel": "#2a78d6",
    "clague_mathews": "#1baf7a",
    "walder_costa_nontunnel": "#eb6834",
}
REF_COLOR = "#74848f"
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"


def _resolve(p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else repo_root() / q


def _plot(curves, out_png: Path, vw_lo: float, vw_hi: float, hw: float) -> Path:
    """Two panels, shared y-axis.

    Tunnel and non-tunnel drainage differ by ~14x in peak and ~20x in duration, so
    one linear time axis squashes the fast case and a log axis distorts the shapes
    into something that no longer looks like a hydrograph. Separate panels keep
    each curve honest; the shared y-axis is what lets the peaks be compared, which
    is the comparison that matters. Both panels enclose the same area, because
    every curve integrates to the same V_w.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fast = ["walder_costa_nontunnel", "froehlich"]
    slow = ["walder_costa_tunnel", "clague_mathews"]
    fig, axes = plt.subplots(
        1, 2, figsize=(11, 5.9), facecolor=SURFACE, sharey=True, gridspec_kw={"wspace": 0.06}
    )
    # Explicit margins rather than tight_layout: figure-level text is not part of
    # the tight layout calculation, so the two fight and the labels collide.
    fig.subplots_adjust(left=0.085, right=0.985, top=0.775, bottom=0.165)

    for ax, keys, title in (
        (axes[0], fast, "Non-tunnel — breach or mechanical failure"),
        (axes[1], slow, "Tunnel — subglacial conduit"),
    ):
        ax.set_facecolor(SURFACE)
        xmax = 0.0
        for key in keys:
            rel = RELATIONS[key]
            band = [c for c in curves if c["relation"] == key]
            if not band:
                continue
            colour = REF_COLOR if rel.reference_only else COLORS[key]
            style = {"ls": (0, (5, 3)), "lw": 1.6} if rel.reference_only else {"lw": 2.1}
            for c in band:
                hours = c["t"] / 3600.0
                xmax = max(xmax, hours[-1])
                ax.plot(
                    hours, c["q"], color=colour, alpha=0.95 if c["edge"] == "hi" else 0.40, **style
                )
            hi = next(c for c in band if c["edge"] == "hi")
            px = hi["t"][int(np.argmax(hi["q"]))] / 3600.0
            # Flip the label inward near the right margin, else it clips.
            right = px > 0.58 * xmax
            ax.annotate(
                rel.label.split(" —")[0],
                (px, hi["q"].max()),
                xytext=(-6 if right else 6, 6),
                textcoords="offset points",
                ha="right" if right else "left",
                color=colour,
                fontsize=9.5,
            )
        ax.set_xlim(0, xmax * 1.12)
        ax.set_title(title, color=INK, fontsize=10.5, loc="left", pad=10)
        ax.set_xlabel("hours since the seal opens", color=INK2)
        ax.tick_params(colors=INK2, labelsize=9)
        for sp in ax.spines.values():
            sp.set_color("#d8d7d2")
        ax.grid(alpha=0.25)

    axes[0].set_ylabel("inflow to Grand Plateau Lake  Q (m³/s)", color=INK2)
    fig.text(
        0.012,
        0.955,
        "Empirical inflow bracket — LGP → GPL",
        color=INK,
        fontsize=13.5,
        ha="left",
        va="top",
    )
    fig.text(
        0.012,
        0.892,
        f"V_w = {vw_lo:.0f}–{vw_hi:.0f} Mm³ · faint = lower bound, solid = upper · "
        f"shared y-axis · equal area under every curve",
        color=INK2,
        fontsize=9.3,
        ha="left",
        va="top",
    )
    fig.text(
        0.012,
        0.035,
        f"Dashed grey is Froehlich (2025), a moraine-breach regression shown for "
        f"reference at H_w = {hw:.0f} m. It does not apply to ice-dam drainage.",
        color=INK2,
        fontsize=8.8,
        ha="left",
        va="bottom",
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
        "--vw-mm3",
        type=float,
        nargs="*",
        default=None,
        help="override V_w; give one or two values (Mm3)",
    )
    ap.add_argument(
        "--hw-m",
        type=float,
        default=None,
        help="head for the Froehlich reference (default: the drawdown)",
    )
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    data = yaml.safe_load(Path(args.data_config).read_text())
    out_dir = _resolve(data["output_dir"])
    vs = out_dir / "volume_summary.json"

    hw = args.hw_m
    if args.vw_mm3:
        bounds = sorted(args.vw_mm3) if len(args.vw_mm3) > 1 else [args.vw_mm3[0]] * 2
    else:
        if not vs.exists():
            raise SystemExit(f"{vs} not found; run run_lakes or pass --vw-mm3")
        summary = json.loads(vs.read_text())
        bounds = summary["Vw_Mm3_range"]
        if hw is None:
            hw = summary["drawdown_m"]
        print(f"V_w {bounds[0]:.1f}-{bounds[1]:.1f} Mm3 (from {vs.name})")
    if hw is None:
        hw = 100.0
    print(f"head for the Froehlich reference: {hw:.1f} m\n")

    rows, curves = [], []
    for edge, vw in (("lo", bounds[0]), ("hi", bounds[1])):
        df = summarise(vw, hw_m=hw)
        df.insert(0, "edge", edge)
        df.insert(1, "Vw_Mm3", vw)
        rows.append(df)
        for _, r in df.iterrows():
            t, q = hydrograph(r.Qp_m3s, vw * 1e6, r["shape"])
            curves.append({"relation": r.relation, "edge": edge, "t": t, "q": q})
    table = pd.concat(rows, ignore_index=True)

    print(
        f"{'relation':26s} {'mode':14s} {'Qp lo':>9s} {'Qp hi':>9s} {'dur lo':>8s} {'dur hi':>8s}"
    )
    for key, rel in RELATIONS.items():
        sub = table[table.relation == key]
        if sub.empty:
            continue
        lo = sub[sub.edge == "lo"].iloc[0]
        hi = sub[sub.edge == "hi"].iloc[0]
        tag = "  [reference only]" if rel.reference_only else ""
        print(
            f"{key:26s} {rel.mode:14s} {lo.Qp_m3s:9.0f} {hi.Qp_m3s:9.0f} "
            f"{lo.duration_h:7.1f}h {hi.duration_h:7.1f}h{tag}"
        )

    applicable = table[~table.reference_only]
    print(
        f"\n  applicable peak-discharge range: "
        f"{applicable.Qp_m3s.min():.0f} – {applicable.Qp_m3s.max():.0f} m3/s"
    )
    print(f"  worst volume error across all curves: {table.volume_error_frac.max():.2e}")
    tun = applicable[applicable["mode"] == "tunnel"].Qp_m3s
    non = applicable[applicable["mode"] == "non_tunnel"].Qp_m3s
    print(
        f"  tunnel {tun.min():.0f}–{tun.max():.0f} vs non-tunnel "
        f"{non.min():.0f}–{non.max():.0f} m3/s "
        f"— a factor of {non.max() / tun.min():.1f} between drainage modes,"
    )
    print("    which dwarfs the V_w uncertainty. The mode is the decision that matters.")

    long = pd.concat(
        [
            pd.DataFrame(
                {"relation": c["relation"], "edge": c["edge"], "t_s": c["t"], "Q_m3s": c["q"]}
            )
            for c in curves
        ],
        ignore_index=True,
    )
    long.to_csv(out_dir / "inflow_empirical.csv", index=False)
    table.to_csv(out_dir / "inflow_empirical_summary.csv", index=False)
    png = _plot(curves, out_dir / "inflow_empirical.png", bounds[0], bounds[1], hw)
    (out_dir / "inflow_empirical_summary.json").write_text(
        json.dumps(
            {
                "created_utc": datetime.now(UTC).isoformat(),
                "Vw_Mm3_bounds": bounds,
                "hw_m": hw,
                "relations": table.to_dict("records"),
            },
            indent=2,
            default=str,
        )
    )
    print(f"\n  wrote inflow_empirical.csv, inflow_empirical_summary.csv/json, {png.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
