"""Run provenance and the baseline report (CLAUDE.md Sections 5 and 8.3 Step 7)."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ..config import Config, repo_root
from ..engine import RunResult


def _checksum(path: Path, limit: int = 64 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    read = 0
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
            read += len(chunk)
            if read >= limit:
                h.update(b"TRUNCATED")
                break
    return h.hexdigest()


def _git_state() -> dict:
    try:
        root = repo_root()
        commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                                    capture_output=True, text=True, check=True).stdout.strip())
        return {"commit": commit, "dirty": dirty}
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {"commit": None, "dirty": None, "note": "not a git repository"}


def write_manifest(cfg: Config, result: RunResult, run_dir: Path) -> Path:
    """Record everything needed to reproduce the run."""
    run_dir.mkdir(parents=True, exist_ok=True)
    inputs = {}
    for key in ("surface", "bed"):
        try:
            p = cfg.build_path(key)
            inputs[p.name] = {"path": str(p), "sha256": _checksum(p), "bytes": p.stat().st_size}
        except (KeyError, FileNotFoundError):
            continue
    manifest = {
        "run_id": cfg.run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(cfg.path),
        "config": cfg.raw,
        "git": _git_state(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "seed": int(cfg.sampling["seed"]),
        "inputs": inputs,
        "phi_crit_m": {str(k): v for k, v in result.phi_crit.items()},
        "n_realizations": len(result.samples),
    }
    path = run_dir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, default=str))
    return path


def write_report(cfg: Config, result: RunResult, run_dir: Path) -> Path:
    """Write ``baseline_report.md`` with the year-by-year table and a margin plot."""
    run_dir.mkdir(parents=True, exist_ok=True)
    ts = result.timeseries
    s = result.samples.iloc[0]
    cell = ts.attrs.get("breach_cell")
    x, y = ts.attrs.get("breach_xy", (float("nan"), float("nan")))

    lines = [
        f"# Baseline report -- {cfg.run_id}",
        "",
        f"- Breach site: **{cfg.breach.get('description', result.dam_name)}**",
        f"- Realizations: {len(result.samples)}",
        f"- Flotation fraction f = {s['flotation_fraction']:.6f} "
        f"(legacy ratio r = {(1 - s['flotation_fraction']) / s['flotation_fraction']:.4f})",
        f"- Lake level L = {s['lake_level_m']:.2f} m",
        f"- Surface change dh/dt = {s['dhdt_m_per_yr']:.3f} m yr-1",
        f"- Saddle value Phi_crit = **{s['phi_crit_m']:.4f} m**",
        "",
        f"## Breach date",
        "",
        f"- Continuous crossing: **{s['breach_year_continuous']:.3f}**",
        f"- First annual map showing a connected pathway: **{int(s['breach_year'])}**",
        f"- Breach point: cell {cell}, map coordinates ({x:.1f}, {y:.1f})",
        f"- Candidate breach cells at the saddle: {len(result.breach_cells)}",
        "",
        "## Year by year at the breach point",
        "",
        "`margin = f * H_ice - h_lake` in metres; positive means the dam holds.",
        "",
        "| year | surface (m) | ice thickness (m) | lake depth (m) | f*H (m) | margin (m) | threshold (m) | connected |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for _, r in ts.iterrows():
        lines.append(
            f"| {int(r['year'])} | {r['surface_m']:.2f} | {r['ice_thickness_m']:.2f} | "
            f"{r['lake_depth_m']:.2f} | {r['flotation_depth_m']:.2f} | {r['margin_m']:.2f} | "
            f"{r['threshold_m']:.2f} | {'yes' if r['connected'] else 'no'} |"
        )

    if len(result.samples) > 1:
        yrs = result.samples["breach_year_continuous"]
        lines += [
            "",
            "## Ensemble",
            "",
            f"- median {yrs.median():.2f}, 5th {yrs.quantile(0.05):.2f}, "
            f"95th {yrs.quantile(0.95):.2f}",
        ]

    plot = _plot_margin(result, run_dir)
    if plot:
        lines += ["", "## Margin vs time", "", f"![margin]({plot.name})"]

    path = run_dir / "baseline_report.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def _plot_margin(result: RunResult, run_dir: Path) -> Path | None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        return None

    years = np.asarray(result.margin.columns, dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    values = result.margin.to_numpy()
    if values.shape[0] == 1:
        ax.plot(years, values[0], lw=2, color="#1f4e79", label="margin")
    else:
        lo, mid, hi = np.percentile(values, [5, 50, 95], axis=0)
        ax.fill_between(years, lo, hi, color="#1f4e79", alpha=0.25, label="5-95%")
        ax.plot(years, mid, lw=2, color="#1f4e79", label="median")
    ax.axhline(0.0, color="#c00000", lw=1.2, ls="--", label="flotation")
    yr = float(result.samples["breach_year_continuous"].iloc[0])
    if np.isfinite(yr):
        ax.axvline(yr, color="#c00000", lw=1, alpha=0.6)
        ax.annotate(f"{yr:.2f}", (yr, 0), textcoords="offset points", xytext=(6, 8),
                    color="#c00000")
    ax.set_xlabel("year")
    ax.set_ylabel("margin  f·H − h_lake  (m)")
    ax.set_title(f"{result.run_id}: dam margin at the breach point")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = run_dir / "margin_vs_time.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
