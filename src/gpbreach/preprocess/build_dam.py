"""Build a ``DamObject`` from the owner's rasters.

This is the only place that opens a GeoTIFF. It runs once per breach site; the
engine consumes the resulting ``.npz`` (CLAUDE.md Section 3).
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy import ndimage

from ..components.criterion import flotation_field
from ..config import Config
from ..dam import DamObject
from ..io import read_raster
from ..saddle import _structure

log = logging.getLogger(__name__)


def find_basin_seeds(
    surface: np.ndarray,
    bed: np.ndarray,
    f: float,
    reference_threshold: float,
    n_basins: int = 3,
    connectivity: int = 8,
) -> list[dict]:
    """Locate the deepest cell of each major basin of the flotation field.

    The basins are the connected components of ``{Phi_f <= reference_threshold}``.
    Each is represented by its minimum-``Phi`` cell, which is a stable interior
    point: it does not move when ``f``, the lake level, or the thinning rate
    change, so the same seeds remain valid across a Monte Carlo ensemble.

    Returned largest-first. Which basins correspond to which lakes is for the
    owner to confirm (D-005).
    """
    phi = np.where(
        np.isfinite(surface) & np.isfinite(bed), flotation_field(surface, bed, f), np.inf
    )
    labels, _ = ndimage.label(phi <= reference_threshold, structure=_structure(connectivity))
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    out = []
    for rank, label in enumerate(np.argsort(sizes)[::-1][:n_basins], start=1):
        if sizes[label] == 0:
            break
        cells = np.argwhere(labels == label)
        values = phi[cells[:, 0], cells[:, 1]]
        r, c = cells[np.argmin(values)]
        out.append(
            {
                "rank": rank,
                "seed": (int(r), int(c)),
                "n_cells": int(sizes[label]),
                "phi_min_m": float(values.min()),
            }
        )
    return out


def build_dam(cfg: Config, out_path: str | Path | None = None) -> DamObject:
    """Read the configured rasters, resolve the seeds, and write the dam object."""
    build = cfg.dam["build"]
    surface_r = read_raster(cfg.build_path("surface"))
    bed_r = read_raster(cfg.build_path("bed"))
    if surface_r.grid_key() != bed_r.grid_key():
        raise ValueError(
            "surface and bed rasters are not co-registered:\n"
            f"  {surface_r.path.name}: {surface_r.grid_key()}\n"
            f"  {bed_r.path.name}: {bed_r.grid_key()}\n"
            "Resampling is a physics/input decision -- ask the owner (CLAUDE.md 2.2)."
        )

    surface = surface_r.values
    bed = bed_r.values
    valid = np.isfinite(surface) & np.isfinite(bed)
    surface = np.where(valid, surface, np.nan)
    bed = np.where(valid, bed, np.nan)

    # `2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif` carries no nodata flag and uses 0.0
    # as fill over 56% of its extent. Those cells sit outside the bedrock domain,
    # so they never reach the analysis -- but assert it rather than trust it, since
    # a wider bed raster would silently admit 0 m "ground". See D-008.
    fill = build.get("surface_fill_sentinel", 0.0)
    if fill is not None:
        leaked = int(np.count_nonzero(valid & (surface_r.values == fill)))
        if leaked:
            raise ValueError(
                f"{surface_r.path.name}: {leaked:,} cells inside the analysis domain "
                f"equal the fill sentinel {fill}. These are not real elevations. "
                "Mask them or set dam.build.surface_fill_sentinel to null after "
                "confirming with the owner that the value is genuine."
            )

    seeds_cfg = build["seeds"]
    method = seeds_cfg.get("method", "explicit")
    discovered: list[dict] = []
    if method == "explicit":
        seed_lake = tuple(int(v) for v in seeds_cfg["seed_lake"])
        seed_target = tuple(int(v) for v in seeds_cfg["seed_target"])
    elif method == "auto_largest_basins":
        discovered = find_basin_seeds(
            surface,
            bed,
            f=cfg.nominal("flotation_fraction"),
            reference_threshold=float(seeds_cfg["reference_threshold_m"]),
            n_basins=int(seeds_cfg.get("n_basins", 3)),
            connectivity=int(cfg.components.get("connectivity", 8)),
        )
        by_rank = {d["rank"]: d["seed"] for d in discovered}
        seed_lake = by_rank[int(seeds_cfg.get("lake_basin_rank", 1))]
        seed_target = by_rank[int(seeds_cfg.get("target_basin_rank", 2))]
    else:
        raise ValueError(f"unknown seed method {method!r}")

    dam = DamObject(
        name=cfg.breach["name"],
        surface=surface,
        bed=bed,
        seed_lake=seed_lake,
        seed_target=seed_target,
        base_year=int(build["base_year"]),
        transform=surface_r.transform,
        crs=surface_r.crs,
        vertical_datum=build.get("vertical_datum"),
        meta={
            "surface_raster": surface_r.path.name,
            "bed_raster": bed_r.path.name,
            "n_valid_cells": int(valid.sum()),
            "seed_method": method,
            "discovered_basins": discovered,
            "seed_lake_xy": None,
            "seed_target_xy": None,
        },
    )
    dam.meta["seed_lake_xy"] = dam.xy(dam.seed_lake)
    dam.meta["seed_target_xy"] = dam.xy(dam.seed_target)

    out_path = Path(out_path) if out_path else cfg.dam_path
    dam.save(out_path)
    log.info("wrote %s (%d valid cells)", out_path, valid.sum())
    return dam
