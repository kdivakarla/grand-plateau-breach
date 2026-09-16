"""Breach-pathway map, for confirming D-005 against the owner's own reading.

A static figure plus QGIS-loadable CSVs. The point of both is the same question:
*is the route the model opens the same stretch of dam you identified by eye?*

Colour: the flotation field is sequential (one hue, light to dark). The three
overlays -- target basin, source basin, pathway -- take categorical slots 1-3 of
the reference palette, which are the slots documented as clearing the all-pairs
CVD and normal-vision floors on a light surface. Every overlay is also directly
labelled, so identity is never carried by colour alone.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import ndimage

from ..components.criterion import FlotationConnectivity, flotation_field
from ..dam import DamObject
from ..saddle import _structure, minimax_path

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
TARGET = "#2a78d6"  # categorical slot 1 -- basin the water drains into
SOURCE = "#eb6834"  # slot 2 -- basin that drains
PATH = "#1baf7a"  # slot 3 -- breach pathway


def breach_map(
    dam: DamObject,
    f: float,
    level: float,
    dhdt: float,
    connectivity: int = 8,
    out_dir: str | Path = ".",
    label_source: str = "source lake",
    label_target: str = "target lake",
) -> dict:
    """Render the pathway map and write the companion CSVs. Returns the paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    crit = FlotationConnectivity(connectivity=connectivity)
    phi = np.where(dam.valid, flotation_field(dam.surface, dam.bed, f), np.nan)
    T = crit.critical_threshold_for_dam(dam, f)
    path = minimax_path(
        np.where(np.isfinite(phi), phi, np.inf), T, dam.seed_lake, dam.seed_target, connectivity
    )
    crest = tuple(path[np.argmax(phi[path[:, 0], path[:, 1]])]) if path.size else dam.seed_lake

    # Basins just *below* the saddle: the two pools about to merge.
    below = np.isfinite(phi) & (phi < T)
    labels, _ = ndimage.label(below, structure=_structure(connectivity))
    src = labels == labels[dam.seed_lake]
    tgt = labels == labels[dam.seed_target]

    x0, dx, _, y0, _, dy = dam.transform
    nrows, ncols = phi.shape
    extent = [
        x0 / 1000,
        (x0 + ncols * dx) / 1000,
        (y0 + nrows * dy) / 1000,
        y0 / 1000,
    ]  # km, north up

    fig, ax = plt.subplots(figsize=(8.4, 8.0), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.imshow(
        phi,
        extent=extent,
        origin="upper",
        cmap="Greys",
        vmin=np.nanpercentile(phi, 1),
        vmax=np.nanpercentile(phi, 99),
        interpolation="nearest",
    )
    for mask, colour in ((tgt, TARGET), (src, SOURCE)):
        ax.imshow(
            np.where(mask, 1.0, np.nan),
            extent=extent,
            origin="upper",
            cmap=matplotlib.colors.ListedColormap([colour]),
            alpha=0.55,
            interpolation="nearest",
        )
    if path.size:
        xs = np.array([dam.xy(tuple(rc))[0] for rc in path]) / 1000
        ys = np.array([dam.xy(tuple(rc))[1] for rc in path]) / 1000
        ax.plot(xs, ys, color=PATH, lw=1.8, ls=(0, (6, 3)), solid_capstyle="round", zorder=5)

    def mark(cell, text, colour, dxy=(0.25, 0.25)):
        mx, my = dam.xy(cell)
        ax.plot(mx / 1000, my / 1000, "o", ms=9, mfc=colour, mec=SURFACE, mew=2, zorder=6)
        ax.annotate(
            text,
            (mx / 1000 + dxy[0], my / 1000 + dxy[1]),
            color=INK,
            fontsize=10,
            fontweight="medium",
            zorder=7,
            bbox=dict(boxstyle="round,pad=0.25", fc=SURFACE, ec="none", alpha=0.85),
        )

    mark(dam.seed_lake, label_source, SOURCE)
    mark(dam.seed_target, label_target, TARGET)
    mx, my = dam.xy(crest)
    ax.plot(mx / 1000, my / 1000, "o", ms=15, mfc="none", mec=INK, mew=2.2, zorder=8)
    ax.annotate(
        f"breach point\nΦ = {T:.1f} m",
        (mx / 1000 + 0.3, my / 1000 - 0.55),
        color=INK,
        fontsize=10,
        fontweight="semibold",
        zorder=9,
        bbox=dict(boxstyle="round,pad=0.3", fc=SURFACE, ec=INK_2, lw=0.6, alpha=0.9),
    )

    year = (
        dam.base_year
        + crit.breach_time(np.array([T]), np.array([level]), np.array([dhdt]), np.array([f]))[0]
    )
    ax.set_title(
        f"{dam.name}: breach pathway   —   {np.ceil(year):.0f}  (crossing {year:.2f})",
        color=INK,
        fontsize=13,
        pad=22,
        loc="left",
    )
    # Be explicit about what the shaded regions are: sub-level sets of the
    # flotation field just below the saddle, not digitised lake outlines. Any
    # basin with a lower saddle has already merged into them by this threshold.
    ax.text(
        0.0,
        1.012,
        f"buoyant regions at Φ = {T:.1f} m (just below breach); the route is one "
        f"valid path — its crest sets the date",
        transform=ax.transAxes,
        color=INK_2,
        fontsize=9.5,
        va="bottom",
    )
    ax.set_xlabel("easting (km)", color=INK_2, fontsize=10)
    ax.set_ylabel("northing (km)", color=INK_2, fontsize=10)
    ax.tick_params(colors=INK_2, labelsize=9)
    for s in ax.spines.values():
        s.set_color("#d8d7d2")
    ax.set_aspect("equal")
    ax.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="s",
                ls="",
                ms=10,
                mfc=TARGET,
                mec="none",
                label=f"{label_target} basin",
            ),
            Line2D(
                [],
                [],
                marker="s",
                ls="",
                ms=10,
                mfc=SOURCE,
                mec="none",
                label=f"{label_source} basin",
            ),
            Line2D([], [], color=PATH, lw=2.4, label="breach pathway"),
            Line2D(
                [],
                [],
                marker="o",
                ls="",
                ms=10,
                mfc="none",
                mec=INK,
                mew=2,
                label="breach point (saddle)",
            ),
        ],
        loc="lower left",
        frameon=True,
        facecolor=SURFACE,
        edgecolor="#d8d7d2",
        fontsize=9,
        labelcolor=INK,
    )
    fig.tight_layout()
    png = out_dir / f"breach_map_{dam.name}.png"
    fig.savefig(png, dpi=160, facecolor=SURFACE)
    plt.close(fig)

    pts = out_dir / f"breach_points_{dam.name}.csv"
    rows = [
        ("seed_source", *dam.seed_lake),
        ("seed_target", *dam.seed_target),
        ("breach_point", *crest),
    ]
    pts.write_text(
        "name,row,col,x,y,phi_m\n"
        + "".join(
            f"{n},{r},{c},{dam.xy((r, c))[0]:.3f},{dam.xy((r, c))[1]:.3f},{phi[r, c]:.4f}\n"
            for n, r, c in rows
        )
    )
    line = out_dir / f"breach_pathway_{dam.name}.csv"
    line.write_text(
        "order,row,col,x,y,phi_m\n"
        + "".join(
            f"{i},{r},{c},{dam.xy((r, c))[0]:.3f},{dam.xy((r, c))[1]:.3f},{phi[r, c]:.4f}\n"
            for i, (r, c) in enumerate(map(tuple, path))
        )
    )
    return {
        "map": png,
        "points": pts,
        "pathway": line,
        "phi_crit": T,
        "breach_year": year,
        "crest": crest,
    }
