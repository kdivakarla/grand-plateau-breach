"""Hydraulic head surface and the LGP→GPL pass search.

The seal between LGP and GPL opens when LGP's level reaches the **pass head**
``h_pass``: the lowest achievable "highest head" on any subglacial route between
the two lakes. Below that, every route crosses somewhere the ice still holds.

Head under ice, in metres of water (spec §seal.py)::

    h_phi = z_b + k * (rho_i / rho_w) * max(z_s - z_b, 0)

Relation to the Phase 1 baseline
--------------------------------
Writing ``kappa = k * rho_i / rho_w`` and expanding, this is

    h_phi = kappa * z_s + (1 - kappa) * z_b

which is exactly ``components.criterion.flotation_field`` with ``f = kappa``. The
spec's flotation factor ``k`` and decision D-003's flotation fraction ``f`` are
one parameter in different units:

    f = k * rho_i / rho_w          k = f / (rho_i / rho_w)

The baseline's ``f = 0.9091`` is ``k = 0.9914``. ``test_cascade_seal.py`` asserts
this correspondence against ``saddle.saddle_threshold`` so the two halves of the
project cannot drift apart.

Why a priority flood rather than the existing saddle search
-----------------------------------------------------------
``saddle.saddle_threshold`` returns the same ``h_pass`` but finds it by bisecting
over sorted levels, which gives no parent links. ``path.py`` needs the route, so
this module uses a minimax priority flood (Dijkstra with ``max`` in place of
``+``) that yields the spill field and a parent array in one pass.
"""

from __future__ import annotations

import heapq
import logging
from dataclasses import dataclass, field

import numpy as np

from ..constants import RHO_ICE, RHO_WATER

log = logging.getLogger(__name__)

#: Neighbour offsets.
_OFFSETS_8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
_OFFSETS_4 = ((-1, 0), (0, -1), (0, 1), (1, 0))


def kappa_from_k(k: float) -> float:
    """Flotation factor ``k`` -> flotation fraction ``f`` (= kappa). Dimensionless."""
    return k * RHO_ICE / RHO_WATER


def k_from_kappa(kappa: float) -> float:
    """Inverse of :func:`kappa_from_k`. Dimensionless."""
    return kappa * RHO_WATER / RHO_ICE


def masked_smooth(z: np.ndarray, mask: np.ndarray, sigma_m: float, res: float) -> np.ndarray:
    """Gaussian-smooth ``z`` over ``mask`` only, by normalised convolution.

    Parameters
    ----------
    z: field to smooth (m), NaN allowed off the mask.
    mask: boolean, True where ``z`` is meaningful.
    sigma_m: smoothing length in **metres**; 0 returns ``z`` unchanged.
    res: cell size in metres.

    Returns NaN off the mask. Smoothing ``z*m`` and dividing by smoothed ``m``
    stops off-mask zeros from bleeding into the result.
    """
    if sigma_m <= 0:
        return np.where(mask, z, np.nan)
    from scipy.ndimage import gaussian_filter

    sigma_px = sigma_m / res
    zz = np.where(mask, z, 0.0).astype(float)
    mm = mask.astype(float)
    num = gaussian_filter(zz, sigma_px, mode="constant", cval=0.0)
    den = gaussian_filter(mm, sigma_px, mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    return np.where(mask & (den > 1e-12), out, np.nan)


def head_surface(
    zb: np.ndarray,
    zs: np.ndarray,
    classes: np.ndarray,
    lake_levels: dict[int, float],
    k: float = 1.0,
) -> tuple[np.ndarray, dict]:
    """Hydraulic head in metres of water, with land and undefined cells as barriers.

    Ice cells (class 1) take the Shreve head; lake cells take their lake level;
    everything else — land, and any ice cell missing bed or surface — becomes
    ``+inf`` so the search cannot route through it.

    Returns ``(head, info)``; ``info`` counts the barrier cells, so an ice area
    silently lost to missing data shows up rather than being absorbed.
    """
    kappa = kappa_from_k(k)
    head = np.full(zb.shape, np.inf, dtype=np.float64)

    ice = classes == 1
    defined = ice & np.isfinite(zb) & np.isfinite(zs)
    thickness = np.where(defined, zs - zb, 0.0)
    negative = int((defined & (thickness < 0)).sum())
    head[defined] = zb[defined] + kappa * np.maximum(thickness[defined], 0.0)

    for value, level in lake_levels.items():
        if level is None:
            continue
        head[classes == value] = float(level)

    info = {
        "k": k,
        "kappa": kappa,
        "ice_cells": int(ice.sum()),
        "ice_cells_with_data": int(defined.sum()),
        "ice_cells_missing_data": int((ice & ~defined).sum()),
        "negative_thickness_cells": negative,
        "barrier_cells": int(np.isinf(head).sum()),
        "lake_levels": {int(k_): v for k_, v in lake_levels.items()},
    }
    if info["ice_cells_missing_data"]:
        log.warning(
            "%d ice cells lack bed or surface and are barriers", info["ice_cells_missing_data"]
        )
    if negative:
        log.warning("%d ice cells have negative thickness; clipped at 0", negative)
    return head, info


@dataclass
class PassResult:
    """Outcome of one pass search. Heads and spills are metres of water."""

    h_pass: float
    pass_cell: tuple[int, int] | None
    spill: np.ndarray
    parent: np.ndarray
    reached_target: bool
    info: dict = field(default_factory=dict)


def pass_search(
    head: np.ndarray,
    classes: np.ndarray,
    source_class: int = 3,
    target_class: int = 2,
    connectivity: int = 8,
) -> PassResult:
    """Minimax priority flood from ``source_class`` to ``target_class``.

    ``spill[c]`` is the lowest achievable maximum head on any route from the
    source to cell ``c``. Target cells are set to ``-inf`` head before the search
    so that the value stored at them is the barrier between the lakes and does
    not include the target lake's own level. ``h_pass`` is the minimum spill over
    all target cells.

    Returns ``inf`` for ``h_pass`` when no route exists.
    """
    if connectivity == 8:
        offsets = _OFFSETS_8
    elif connectivity == 4:
        offsets = _OFFSETS_4
    else:
        raise ValueError(f"connectivity must be 4 or 8, got {connectivity!r}")

    work = head.copy()
    target = classes == target_class
    source = classes == source_class
    if not source.any():
        raise ValueError(f"no cells of source class {source_class}")
    if not target.any():
        raise ValueError(f"no cells of target class {target_class}")
    work[target] = -np.inf

    nrows, ncols = work.shape
    spill = np.full(work.shape, np.inf, dtype=np.float64)
    parent = np.full((*work.shape, 2), -1, dtype=np.int32)

    heap: list[tuple[float, int, int]] = []
    for r, c in zip(*np.nonzero(source), strict=True):
        v = work[r, c]
        if np.isinf(v) and v > 0:
            continue
        spill[r, c] = v
        heap.append((float(v), int(r), int(c)))
    heapq.heapify(heap)

    popped = 0
    while heap:
        v, r, c = heapq.heappop(heap)
        if v > spill[r, c]:
            continue  # stale entry
        popped += 1
        for dr, dc in offsets:
            rr, cc = r + dr, c + dc
            if not (0 <= rr < nrows and 0 <= cc < ncols):
                continue
            hn = work[rr, cc]
            if np.isinf(hn) and hn > 0:
                continue  # barrier
            cand = v if v > hn else hn
            if cand < spill[rr, cc]:
                spill[rr, cc] = cand
                parent[rr, cc] = (r, c)
                heapq.heappush(heap, (float(cand), rr, cc))

    target_spill = spill[target]
    finite = target_spill[np.isfinite(target_spill)]
    if finite.size == 0:
        return PassResult(
            float("inf"),
            None,
            spill,
            parent,
            False,
            {"cells_visited": popped, "reason": "no route to target"},
        )

    h_pass = float(finite.min())
    # The pass cell is the highest-head cell on the route to the best target cell.
    tr, tc = _best_target_cell(spill, target)
    pass_cell = _highest_on_path(work, parent, (tr, tc))
    return PassResult(
        h_pass,
        pass_cell,
        spill,
        parent,
        True,
        {
            "cells_visited": popped,
            "target_cells": int(target.sum()),
            "best_target_cell": [int(tr), int(tc)],
        },
    )


def _best_target_cell(spill: np.ndarray, target: np.ndarray) -> tuple[int, int]:
    masked = np.where(target, spill, np.inf)
    idx = int(np.argmin(masked))
    return divmod(idx, spill.shape[1])


def _highest_on_path(
    head: np.ndarray, parent: np.ndarray, start: tuple[int, int]
) -> tuple[int, int] | None:
    """Cell with the greatest head on the parent chain from ``start``."""
    best, best_h = None, -np.inf
    r, c = start
    seen = 0
    limit = head.size + 1
    while r >= 0 and c >= 0 and seen < limit:
        h = head[r, c]
        if np.isfinite(h) and h > best_h:
            best_h, best = h, (int(r), int(c))
        pr, pc = parent[r, c]
        if pr < 0:
            break
        r, c = int(pr), int(pc)
        seen += 1
    return best


def overtopping_search(
    zs: np.ndarray,
    classes: np.ndarray,
    lake_levels: dict[int, float],
    source_class: int = 3,
    target_class: int = 2,
    connectivity: int = 8,
) -> PassResult:
    """Same search on **ice surface elevation**, for comparison with flotation onset.

    Answers "how high would LGP have to rise to spill over the ice?" — an upper
    bound that the flotation route should beat.
    """
    surf = np.where(classes == 1, zs, np.inf).astype(np.float64)
    surf[~np.isfinite(surf)] = np.inf
    for value, level in lake_levels.items():
        if level is not None:
            surf[classes == value] = float(level)
    return pass_search(surf, classes, source_class, target_class, connectivity)


def competing_outlets(
    head: np.ndarray,
    classes: np.ndarray,
    candidates: dict[str, int],
    target_class: int = 2,
    connectivity: int = 8,
) -> dict[str, dict]:
    """Pass head from each candidate receiving lake to the target lake.

    The lake giving LGP the lowest pass is the one it drains to first.
    """
    out: dict[str, dict] = {}
    for name, source_class in candidates.items():
        if not (classes == source_class).any():
            out[name] = {"status": "absent from grid", "h_pass": None}
            continue
        try:
            res = pass_search(head, classes, source_class, target_class, connectivity)
        except ValueError as exc:
            out[name] = {"status": str(exc), "h_pass": None}
            continue
        out[name] = {
            "status": "ok",
            "h_pass": res.h_pass,
            "pass_cell": list(res.pass_cell) if res.pass_cell else None,
            "reached": res.reached_target,
        }
    return out
