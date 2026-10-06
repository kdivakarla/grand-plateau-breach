"""Drainage path from LGP to GPL, and the along-path profile.

Turns the pass search into one route plus the handful of numbers the flood model
needs (spec §path.py).

Two segments, found differently
-------------------------------
*Above the pass* the route is the minimax path: the one the priority flood found,
which is the cheapest way to reach the pass at all. *Below the pass* the minimax
path is **not unique** -- once water is over the pass, any descending route works
and the flood's parent links pick an arbitrary one. ``refine_downstream``
replaces that segment with steepest descent on the head surface, which is where
water actually goes once the seal opens.

Direction convention
--------------------
Paths run **LGP -> GPL**, the direction water flows. The pass search seeds at GPL
and stores parent links pointing back toward it, so following parents from an LGP
cell already yields that order.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

_D8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))


def trace_path(parent: np.ndarray, start: tuple[int, int]) -> np.ndarray:
    """Follow parent links from ``start`` to the seed, as an ``(n, 2)`` array.

    ``start`` should be the target-lake cell the search reached (``PassResult``'s
    ``info["best_target_cell"]``). The returned order is start-first, i.e.
    LGP -> GPL for the standard seeding.
    """
    cells = [(int(start[0]), int(start[1]))]
    limit = parent.shape[0] * parent.shape[1] + 1
    while len(cells) < limit:
        r, c = cells[-1]
        pr, pc = parent[r, c]
        if pr < 0 or pc < 0:
            break
        cells.append((int(pr), int(pc)))
    else:  # pragma: no cover - a cycle would mean a corrupt parent array
        raise RuntimeError("parent chain did not terminate; possible cycle")
    return np.array(cells, dtype=int)


def steepest_descent(
    head: np.ndarray,
    start: tuple[int, int],
    classes: np.ndarray,
    stop_class: int,
    max_steps: int | None = None,
    forbidden: set[tuple[int, int]] | None = None,
    parent: np.ndarray | None = None,
) -> np.ndarray:
    """D8 steepest descent on ``head`` from ``start`` until ``stop_class`` is hit.

    Pure descent does not work on a real potential field. Two things stop it:
    local minima, and **plateaus** in the flood-filled surface where many cells
    share one bottleneck value (on the 2018 geometry, nine consecutive cells sit
    at 232.714 m). A strictly-downhill rule stalls on both.

    ``parent`` supplies the fallback. When no neighbour is strictly lower, the
    walk steps to the pass search's parent cell, which by construction lies one
    step closer to the flood's source. Descent is therefore physical wherever the
    head genuinely falls, and flood-guided only where it does not, and it always
    terminates at the source.

    ``forbidden`` cells are never entered; ``refine_downstream`` uses it to
    exclude the upstream path, since at the pass the neighbour back toward LGP is
    often lower than the one toward GPL.
    """
    nrows, ncols = head.shape
    max_steps = max_steps or (nrows * ncols)
    cells = [(int(start[0]), int(start[1]))]
    seen = {cells[0]}
    if forbidden:
        seen |= set(forbidden)
    fallbacks = 0

    for _ in range(max_steps):
        r, c = cells[-1]
        if classes[r, c] == stop_class:
            if fallbacks:
                log.info("descent used the parent fallback %d times (plateaus/pits)", fallbacks)
            return np.array(cells, dtype=int)

        best, best_h = None, head[r, c]
        for dr, dc in _D8:
            rr, cc = r + dr, c + dc
            if not (0 <= rr < nrows and 0 <= cc < ncols) or (rr, cc) in seen:
                continue
            h = head[rr, cc]
            if np.isfinite(h) and h < best_h:
                best, best_h = (rr, cc), h

        if best is None and parent is not None:
            pr, pc = parent[r, c]
            if pr >= 0 and (int(pr), int(pc)) not in seen:
                best = (int(pr), int(pc))
                fallbacks += 1

        if best is None:
            log.warning(
                "descent stalled at %s before reaching class %d "
                "(no lower neighbour, no usable parent)",
                (r, c),
                stop_class,
            )
            return np.array(cells, dtype=int)
        cells.append(best)
        seen.add(best)

    log.warning("descent hit the step limit")
    return np.array(cells, dtype=int)


def refine_downstream(
    path: np.ndarray,
    head: np.ndarray,
    classes: np.ndarray,
    pass_cell: tuple[int, int],
    target_class: int = 3,
    routing_surface: np.ndarray | None = None,
    parent: np.ndarray | None = None,
) -> np.ndarray:
    """Replace the below-pass segment with steepest descent to ``target_class``.

    ``path`` runs LGP -> GPL. Everything up to ``pass_cell`` is kept; the
    remainder is re-derived by descent. If the pass is not on the path, the path
    is returned unchanged.

    **Refinement is attempted, never forced.** If the descent fails to reach
    ``target_class`` the original path is returned unchanged and a warning is
    logged. Silently returning a truncated route would be far worse than keeping
    the minimax one: on the 2018 geometry the descent stops 63 m below the pass,
    which would have reported ``L_pass`` as 0.06 km instead of 5.6 km and fed
    that straight into the conduit model.

    Pass ``parent=PassResult.parent`` to let the walk fall back to the flood tree
    on plateaus; see :func:`steepest_descent`.
    """
    idx = np.flatnonzero((path[:, 0] == pass_cell[0]) & (path[:, 1] == pass_cell[1]))
    if idx.size == 0:
        log.warning("pass cell %s is not on the traced path; leaving it unrefined", pass_cell)
        return path
    cut = int(idx[0])
    upstream = {(int(r), int(c)) for r, c in path[:cut]}
    surface = head if routing_surface is None else routing_surface
    below = steepest_descent(
        surface, pass_cell, classes, target_class, forbidden=upstream, parent=parent
    )
    if classes[tuple(below[-1])] != target_class:
        log.warning(
            "downstream refinement did not reach class %d (stopped after %d cells); "
            "keeping the unrefined minimax path",
            target_class,
            len(below),
        )
        return path
    return np.vstack([path[:cut], below])


def monotonicity(prof, column: str = "head") -> dict:
    """How far ``column`` departs from monotonically decreasing along the path.

    A perfectly smooth potential field would fall monotonically below the pass.
    Real ones do not, and the size of the departure is a useful signal: it is the
    roughness that defeats strict steepest descent, and it is what the spec's
    ``smoothing_m`` knob exists to control.
    """
    v = prof[column].to_numpy()
    d = np.diff(v)
    rises = d[d > 1e-9]
    return {
        "steps": int(d.size),
        "rising_steps": int(rises.size),
        "rising_fraction": float(rises.size / d.size) if d.size else 0.0,
        "largest_rise_m": float(rises.max()) if rises.size else 0.0,
        "total_rise_m": float(rises.sum()) if rises.size else 0.0,
        "monotonic": bool(rises.size == 0),
    }


def profile(
    path: np.ndarray,
    zb: np.ndarray,
    zs: np.ndarray,
    head: np.ndarray,
    transform,
    h_lgp: float,
    h_gpl: float,
    res: float | None = None,
):
    """Sample bed, surface, head and effective pressures along ``path``.

    Returns a DataFrame with, per cell: ``s`` distance along the path (m),
    ``x``/``y`` map coordinates (m), ``z_b``, ``z_s``, ``H`` thickness (m),
    ``head`` (m of water), and ``N_lake``/``N_gpl`` = head minus each lake level
    (m of water). ``N_lake`` is the same quantity as the margin raster.
    """
    import pandas as pd

    res = res or abs(transform.a)
    r, c = path[:, 0], path[:, 1]
    step = np.hypot(np.diff(r.astype(float)), np.diff(c.astype(float))) * res
    s = np.concatenate([[0.0], np.cumsum(step)])
    # Explicit arithmetic rather than `transform * (col, row)`, which is deprecated.
    x = transform.c + (c + 0.5) * transform.a
    y = transform.f + (r + 0.5) * transform.e
    zbv, zsv, hv = zb[r, c], zs[r, c], head[r, c]
    return pd.DataFrame(
        {
            "row": r,
            "col": c,
            "s": s,
            "x": np.asarray(x),
            "y": np.asarray(y),
            "z_b": zbv,
            "z_s": zsv,
            "H": zsv - zbv,
            "head": hv,
            "N_lake": hv - h_lgp,
            "N_gpl": hv - h_gpl,
        }
    )


@dataclass
class PathSummary:
    """Scalars the flood model consumes. Lengths in m, heads in m of water."""

    l_total_m: float
    l_margin_m: float
    l_pass_m: float
    bed_sill_m: float
    bed_sill_s_m: float
    flotation_length_m: float
    h_at_pass_m: float
    thickness_at_pass_m: float
    thickness_min_m: float
    thickness_mean_m: float
    n_representative_m: float
    n_representative_mode: str
    n_at_pass_m: float
    n_mean_seal_m: float
    n_min_m: float

    def as_dict(self) -> dict:
        return {
            k: (None if v is None or (isinstance(v, float) and not np.isfinite(v)) else v)
            for k, v in self.__dict__.items()
        }


def summarize(
    prof,
    classes: np.ndarray,
    h_lgp: float,
    pass_cell: tuple[int, int] | None = None,
    lgp_class: int = 2,
    n_mode: str = "at_pass",
) -> PathSummary:
    """Reduce a profile to the scalars M4-M6 need.

    ``bed_sill`` is the highest bed elevation along the route: LGP cannot drain
    below it, so it caps the released volume. ``n_mode`` selects the
    representative effective pressure -- ``at_pass`` (default, per D-014),
    ``mean_seal`` (mean over the sealed segment) or ``min`` (minimum along path).
    """
    cls = classes[prof["row"].to_numpy(), prof["col"].to_numpy()]
    s = prof["s"].to_numpy()

    # Where the route leaves LGP: the last cell still in the lake.
    in_lgp = np.flatnonzero(cls == lgp_class)
    s_margin = s[in_lgp[-1]] if in_lgp.size else 0.0

    i_pass = None
    if pass_cell is not None:
        hit = np.flatnonzero(
            (prof["row"].to_numpy() == pass_cell[0]) & (prof["col"].to_numpy() == pass_cell[1])
        )
        if hit.size:
            i_pass = int(hit[0])
    if i_pass is None:
        i_pass = int(np.nanargmax(prof["head"].to_numpy()))

    zb = prof["z_b"].to_numpy()
    sill_i = int(np.nanargmax(zb))
    head = prof["head"].to_numpy()
    n_lake = prof["N_lake"].to_numpy()
    sealed = head > h_lgp
    thick = prof["H"].to_numpy()

    n_at_pass = float(n_lake[i_pass])
    n_mean_seal = float(np.nanmean(n_lake[sealed])) if sealed.any() else float("nan")
    n_min = float(np.nanmin(n_lake))
    n_rep = {"at_pass": n_at_pass, "mean_seal": n_mean_seal, "min": n_min}.get(n_mode)
    if n_rep is None:
        raise ValueError(f"n_mode must be at_pass|mean_seal|min, got {n_mode!r}")

    return PathSummary(
        l_total_m=float(s[-1]),
        l_margin_m=float(s[-1] - s_margin),
        l_pass_m=float(s[-1] - s[i_pass]),
        bed_sill_m=float(zb[sill_i]),
        bed_sill_s_m=float(s[sill_i]),
        flotation_length_m=float(np.sum(np.diff(s)[~sealed[:-1]])) if s.size > 1 else 0.0,
        h_at_pass_m=float(head[i_pass]),
        thickness_at_pass_m=float(thick[i_pass]),
        thickness_min_m=float(np.nanmin(thick)),
        thickness_mean_m=float(np.nanmean(thick)),
        n_representative_m=float(n_rep),
        n_representative_mode=n_mode,
        n_at_pass_m=n_at_pass,
        n_mean_seal_m=n_mean_seal,
        n_min_m=n_min,
    )


def write_path_gpkg(prof, crs, out_path: str | Path, attributes: dict | None = None) -> Path:
    """Write the route as a single LineString carrying the summary attributes."""
    import geopandas as gpd
    from shapely.geometry import LineString

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    line = LineString(list(zip(prof["x"], prof["y"], strict=True)))
    attrs = {
        k: [v]
        for k, v in (attributes or {}).items()
        if isinstance(v, (int, float, str, bool)) or v is None
    }
    gpd.GeoDataFrame(attrs or {"id": [1]}, geometry=[line], crs=crs).to_file(
        out_path, layer="drainage_path", driver="GPKG"
    )
    return out_path
