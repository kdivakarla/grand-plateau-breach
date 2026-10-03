"""Build one canonical analysis grid and put every input on it.

Why this exists
---------------
The cascade analysis needs all rasters cell-for-cell co-registered, and it needs a
footprint that *contains the lakes*. Neither holds for the inputs as prepared:

* the 12.5 x 12.7 km footprint truncates all three lakes -- Alsek to 1.5 % of its
  area, LGP to 25 %, GPL to 52 % -- so stage-area-volume curves, and therefore the
  released volume and the whole hydrograph, cannot be computed on it;
* ``head_2018`` and ``margin_2018_lgp`` are 2497 columns at 4.999597 x 5.000177 m
  with a 2.1 m origin offset, so they are not aligned to the bed grid;
* ``surf_2m`` and ``surf_filled`` have **anisotropic** pixels, 2.382442 x
  1.288162 m, from a reprojection that did not fix the output resolution.

This module defines a target grid and resamples everything onto it.

Alignment rule
--------------
The target grid keeps the **phase** of the reference grid: every target cell edge
falls on a reference cell edge, so the extended product is pixel-compatible with
the existing 5 m products and with the Phase 1 baseline rasters. The origin moves
only by whole multiples of the cell size.

What this does not do
---------------------
Extending the footprint does not create data. The bed is only known inside the
original clip and stays NaN outside it; that is correct and is reported rather
than filled. Deciding how to extend the bed is a physics/input decision for the
owner (CLAUDE.md Section 2.2).
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.enums import Resampling
from rasterio.warp import reproject

log = logging.getLogger(__name__)

#: Resampling per kind of layer. Categorical data must never be interpolated.
CONTINUOUS = Resampling.bilinear
CATEGORICAL = Resampling.nearest


@dataclass(frozen=True)
class TargetGrid:
    """A regular grid: transform, shape and CRS.

    Attributes
    ----------
    transform: affine mapping (col, row) -> (x, y) in CRS units (metres).
    shape: ``(height, width)`` in cells.
    crs: coordinate reference system.
    """

    transform: Affine
    shape: tuple[int, int]
    crs: rasterio.crs.CRS

    @property
    def res(self) -> tuple[float, float]:
        """Cell size (x, y) in metres."""
        return (abs(self.transform.a), abs(self.transform.e))

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        """``(left, bottom, right, top)`` in metres."""
        h, w = self.shape
        left, top = self.transform.c, self.transform.f
        return (left, top - h * abs(self.transform.e), left + w * abs(self.transform.a), top)

    @property
    def n_cells(self) -> int:
        return self.shape[0] * self.shape[1]

    def describe(self) -> str:
        left, bottom, right, top = self.bounds
        return (
            f"{self.shape[1]} x {self.shape[0]} cells @ {self.res[0]:g} m  "
            f"({(right - left) / 1000:.2f} x {(top - bottom) / 1000:.2f} km)  "
            f"origin ({left:.1f}, {top:.1f})  {self.n_cells / 1e6:.1f} M cells"
        )


def aligned_grid(
    bounds: tuple[float, float, float, float],
    reference: rasterio.DatasetReader,
    res: float | None = None,
) -> TargetGrid:
    """Smallest grid containing ``bounds`` that keeps ``reference``'s phase.

    Parameters
    ----------
    bounds: ``(left, bottom, right, top)`` in metres, in the reference CRS.
    reference: open raster whose origin phase and CRS are preserved.
    res: target cell size in metres; defaults to the reference's (square) size.

    The returned grid's edges coincide with reference cell edges, so resampling
    between the two is a pure crop/pad with no interpolation where they overlap.
    """
    rx, ry = abs(reference.transform.a), abs(reference.transform.e)
    if res is None:
        if not np.isclose(rx, ry):
            raise ValueError(
                f"reference has anisotropic pixels ({rx:g} x {ry:g} m); pass res= explicitly"
            )
        res = rx

    ox, oy = reference.transform.c, reference.transform.f
    left, bottom, right, top = bounds
    # Snap outward onto the reference phase.
    new_left = ox + np.floor((left - ox) / res) * res
    new_top = oy + np.ceil((top - oy) / res) * res
    width = int(np.ceil((right - new_left) / res))
    height = int(np.ceil((new_top - bottom) / res))
    return TargetGrid(
        transform=Affine(res, 0.0, new_left, 0.0, -res, new_top),
        shape=(height, width),
        crs=reference.crs,
    )


def union_bounds(
    raster_paths: Iterable[str | Path] = (),
    vector_paths: Iterable[str | Path] = (),
    crs: rasterio.crs.CRS | None = None,
    buffer_m: float = 0.0,
) -> tuple[float, float, float, float]:
    """Union of every input's extent, reprojected to ``crs``, plus a halo."""
    import geopandas as gpd

    boxes = []
    for p in raster_paths:
        with rasterio.open(p) as src:
            b = src.bounds
            if crs is not None and src.crs != crs:
                from rasterio.warp import transform_bounds

                b = transform_bounds(src.crs, crs, *b)
            boxes.append(tuple(b))
    for p in vector_paths:
        g = gpd.read_file(p)
        if crs is not None and g.crs != crs:
            g = g.to_crs(crs)
        boxes.append(tuple(g.total_bounds))
    if not boxes:
        raise ValueError("no inputs given")
    arr = np.array(boxes, dtype=float)
    return (
        arr[:, 0].min() - buffer_m,
        arr[:, 1].min() - buffer_m,
        arr[:, 2].max() + buffer_m,
        arr[:, 3].max() + buffer_m,
    )


def regrid(
    src_path: str | Path, target: TargetGrid, categorical: bool = False, dtype: str = "float32"
) -> np.ndarray:
    """Resample one raster onto ``target``.

    Continuous layers are bilinear and come back as float with NaN nodata;
    categorical layers are nearest-neighbour and keep their integer coding.
    """
    with rasterio.open(src_path) as src:
        if categorical:
            out = np.zeros(target.shape, dtype=src.dtypes[0])
            dst_nodata = src.nodata if src.nodata is not None else 0
            out[:] = dst_nodata
            reproject(
                source=rasterio.band(src, 1),
                destination=out,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=target.transform,
                dst_crs=target.crs,
                src_nodata=src.nodata,
                dst_nodata=dst_nodata,
                resampling=CATEGORICAL,
            )
            return out

        out = np.full(target.shape, np.nan, dtype=dtype)
        reproject(
            source=rasterio.band(src, 1),
            destination=out,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=target.transform,
            dst_crs=target.crs,
            src_nodata=src.nodata,
            dst_nodata=np.nan,
            resampling=CONTINUOUS,
        )
    # GDAL float32 nodata sentinels survive resampling as huge magnitudes.
    out[np.abs(out) > 1e30] = np.nan
    return out


def write(
    path: str | Path, data: np.ndarray, target: TargetGrid, nodata: float | int | None = None
) -> Path:
    """Write ``data`` on ``target`` as a tiled, compressed GeoTIFF."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if nodata is None:
        nodata = np.nan if np.issubdtype(data.dtype, np.floating) else 0
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=target.shape[0],
        width=target.shape[1],
        count=1,
        dtype=data.dtype,
        crs=target.crs,
        transform=target.transform,
        nodata=nodata,
        compress="deflate",
        tiled=True,
        blockxsize=512,
        blockysize=512,
        predictor=2 if np.issubdtype(data.dtype, np.floating) else 1,
    ) as dst:
        dst.write(data, 1)
    return path


def coverage_report(
    classes: np.ndarray,
    target: TargetGrid,
    lake_paths: str | Path | None = None,
    class_names: dict[int, str] | None = None,
) -> dict:
    """Fraction of each digitised lake now captured by the grid."""
    import geopandas as gpd

    names = class_names or {2: "LGP", 3: "GPL", 4: "Alsek"}
    cell_km2 = target.res[0] * target.res[1] / 1e6
    out = {}
    polys = {}
    if lake_paths is not None:
        g = gpd.read_file(lake_paths)
        if g.crs != target.crs:
            g = g.to_crs(target.crs)
        for _, row in g.iterrows():
            polys[int(row["class"])] = row.geometry.area / 1e6
    for value, name in names.items():
        in_grid = float((classes == value).sum() * cell_km2)
        entry = {"class": value, "in_grid_km2": round(in_grid, 3)}
        if value in polys:
            entry["digitised_km2"] = round(polys[value], 3)
            entry["captured_pct"] = round(100.0 * in_grid / polys[value], 1)
        out[name] = entry
    return out


def rasterize_classes(
    target: TargetGrid,
    ice_path: str | Path | None = None,
    lakes_path: str | Path | None = None,
    ice_value: int = 1,
    land_value: int = 0,
) -> np.ndarray:
    """Build the class raster from the digitised outlines, on ``target``.

    Resampling an existing class raster onto a larger grid cannot extend it --
    the new area is simply nodata. The classes must be rasterised from the
    polygons instead, which is what this does.

    Coding follows the spec: 0 land, 1 ice, 2 LGP, 3 GPL, 4 Alsek. **Lakes are
    burned after ice and therefore win** where outlines overlap, since a cell
    covered by water is water regardless of the glacier outline around it.

    This re-expresses the owner's own polygons on a bigger canvas; it does not
    reinterpret them. ``verify_against`` checks that claim.
    """
    import geopandas as gpd
    from rasterio.features import rasterize

    out = np.full(target.shape, land_value, dtype="uint8")

    if ice_path is not None:
        ice = gpd.read_file(ice_path)
        if ice.crs != target.crs:
            ice = ice.to_crs(target.crs)
        rasterize(
            ((geom, ice_value) for geom in ice.geometry),
            out=out,
            transform=target.transform,
            all_touched=False,
        )

    if lakes_path is not None:
        lakes = gpd.read_file(lakes_path)
        if lakes.crs != target.crs:
            lakes = lakes.to_crs(target.crs)
        if "class" not in lakes.columns:
            raise ValueError(f"{lakes_path} has no 'class' field")
        # Burn in ascending class order so the result is deterministic.
        pairs = sorted(
            ((int(r["class"]), r.geometry) for _, r in lakes.iterrows()), key=lambda t: t[0]
        )
        rasterize(
            ((geom, value) for value, geom in pairs),
            out=out,
            transform=target.transform,
            all_touched=False,
        )
    return out


def verify_against(
    new: np.ndarray, target: TargetGrid, original_path: str | Path, sample_limit: int = 2_000_000
) -> dict:
    """Compare a rebuilt class raster with the original over their overlap.

    Returns agreement statistics and a per-class confusion summary. Disagreement
    is expected only at polygon edges (half-cell rasterisation differences); a
    large or structured disagreement means the rebuild changed the owner's
    classification and must be reported, not accepted.
    """
    with rasterio.open(original_path) as src:
        orig = src.read(1)
        orig_t, orig_shape = src.transform, src.shape
        res = abs(orig_t.a)
        # Index of the original's origin within the target grid.
        col0 = round((orig_t.c - target.transform.c) / res)
        row0 = round((target.transform.f - orig_t.f) / res)

    h, w = orig_shape
    sub = new[row0 : row0 + h, col0 : col0 + w]
    if sub.shape != orig.shape:
        return {
            "status": "no overlap or shape mismatch",
            "new_sub_shape": list(sub.shape),
            "orig_shape": list(orig.shape),
        }

    same = sub == orig
    stats = {
        "status": "compared",
        "overlap_cells": int(same.size),
        "agree_cells": int(same.sum()),
        "agree_pct": round(100.0 * float(same.mean()), 3),
        "per_class": {},
    }
    for v in np.unique(np.concatenate([np.unique(orig), np.unique(sub)])):
        in_orig = int((orig == v).sum())
        in_new = int((sub == v).sum())
        stats["per_class"][int(v)] = {
            "original_cells": in_orig,
            "rebuilt_cells": in_new,
            "delta_cells": in_new - in_orig,
        }
    if same.size <= sample_limit:
        stats["disagree_cells"] = int((~same).sum())
    return stats
