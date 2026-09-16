"""Minimum-bottleneck ("saddle") thresholds on a raster.

The breach criterion is a connectivity test on a sub-level set of a static
field (docs/baseline_method.md). The quantity that matters is therefore the
*saddle value*

    Phi_crit = min over paths p from seed A to seed B of  max over cells in p of Phi

i.e. the lowest threshold at which the two seeds first belong to one connected
component of ``{Phi <= T}``. Everything downstream is algebra on that scalar.

The value is found exactly -- not to a tolerance -- by bisecting over the sorted
array of distinct cell values, so the answer is always an actual cell value.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

#: 8-connectivity matches how the owner read connectivity off the QGIS rasters.
STRUCT_8 = np.ones((3, 3), dtype=bool)
STRUCT_4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)


def _structure(connectivity: int) -> np.ndarray:
    if connectivity == 8:
        return STRUCT_8
    if connectivity == 4:
        return STRUCT_4
    raise ValueError(f"connectivity must be 4 or 8, got {connectivity!r}")


def connected_at(
    field: np.ndarray,
    threshold: float,
    seed_a: tuple[int, int],
    seed_b: tuple[int, int],
    connectivity: int = 8,
) -> bool:
    """Are both seeds in one connected component of ``{field <= threshold}``?"""
    labels, _ = ndimage.label(field <= threshold, structure=_structure(connectivity))
    la, lb = labels[seed_a], labels[seed_b]
    return bool(la != 0 and la == lb)


def saddle_threshold(
    field: np.ndarray,
    seed_a: tuple[int, int],
    seed_b: tuple[int, int],
    connectivity: int = 8,
) -> float:
    """Exact minimum threshold at which ``seed_a`` and ``seed_b`` connect.

    ``field`` may contain NaN for cells outside the domain; those never join a
    component. Returns ``inf`` if the seeds never connect within the domain.
    """
    work = np.where(np.isfinite(field), field, np.inf)
    finite = work[np.isfinite(work)]
    if finite.size == 0:
        return float("inf")
    levels = np.unique(finite)

    # The seeds themselves must be inside the domain.
    if not (np.isfinite(work[seed_a]) and np.isfinite(work[seed_b])):
        return float("inf")
    if not connected_at(work, levels[-1], seed_a, seed_b, connectivity):
        return float("inf")

    lo, hi = -1, levels.size - 1  # lo: known-disconnected index, hi: known-connected
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if connected_at(work, levels[mid], seed_a, seed_b, connectivity):
            hi = mid
        else:
            lo = mid
    return float(levels[hi])


def locate_saddle(
    field: np.ndarray,
    threshold: float,
    seed_a: tuple[int, int],
    seed_b: tuple[int, int],
    connectivity: int = 8,
) -> np.ndarray:
    """Cells that complete the connection at ``threshold`` -- the breach point.

    Returns an ``(n, 2)`` array of ``(row, col)``. These are the cells whose
    value equals the saddle threshold and which touch both basins immediately
    below it, so they can be exported and checked in QGIS.
    """
    work = np.where(np.isfinite(field), field, np.inf)
    struct = _structure(connectivity)
    below = work < threshold
    labels, _ = ndimage.label(below, structure=struct)
    la, lb = labels[seed_a], labels[seed_b]
    if la == 0 or lb == 0:
        return np.empty((0, 2), dtype=int)
    touch_a = ndimage.binary_dilation(labels == la, structure=struct)
    touch_b = ndimage.binary_dilation(labels == lb, structure=struct)
    at_level = work == threshold
    strict = at_level & touch_a & touch_b
    if strict.any():
        return np.argwhere(strict)
    # The bridge can need more than one cell at the saddle level (exact ties in
    # the field). Fall back to the cells at that level touching either basin.
    return np.argwhere(at_level & (touch_a | touch_b))


def minimax_path(
    field: np.ndarray,
    threshold: float,
    seed_a: tuple[int, int],
    seed_b: tuple[int, int],
    connectivity: int = 8,
) -> np.ndarray:
    """A path from ``seed_a`` to ``seed_b`` staying at or below ``threshold``.

    With ``threshold = saddle_threshold(...)`` this is the breach pathway: the
    route the buoyant zone opens along, and the thing to check against the
    stretch of dam the owner identified by eye. Breadth-first, so the path is
    also the shortest such route in cell counts.

    Returns an ``(n, 2)`` array of ``(row, col)``, empty if no path exists.
    """
    from collections import deque

    work = np.where(np.isfinite(field), field, np.inf)
    passable = work <= threshold
    if not (passable[seed_a] and passable[seed_b]):
        return np.empty((0, 2), dtype=int)

    nrows, ncols = work.shape
    if connectivity == 8:
        steps = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    else:
        steps = [(-1, 0), (0, -1), (0, 1), (1, 0)]

    came_from = np.full((nrows, ncols, 2), -1, dtype=np.int32)
    seen = np.zeros((nrows, ncols), dtype=bool)
    seen[seed_a] = True
    queue = deque([seed_a])
    while queue:
        r, c = queue.popleft()
        if (r, c) == seed_b:
            break
        for dr, dc in steps:
            rr, cc = r + dr, c + dc
            if 0 <= rr < nrows and 0 <= cc < ncols and passable[rr, cc] and not seen[rr, cc]:
                seen[rr, cc] = True
                came_from[rr, cc] = (r, c)
                queue.append((rr, cc))
    if not seen[seed_b]:
        return np.empty((0, 2), dtype=int)

    path = [seed_b]
    while path[-1] != seed_a:
        path.append(tuple(came_from[path[-1]]))
    return np.array(path[::-1], dtype=int)
