"""Grid extension: alignment, rasterisation and the faithfulness check.

Synthetic cases with known answers, per the spec's testing convention, plus a
check against the real extended grid when it has been built.
"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio
from affine import Affine

from gpbreach.cascade.grid import TargetGrid, aligned_grid, rasterize_classes, verify_against
from gpbreach.config import repo_root

CRS = rasterio.crs.CRS.from_string("ESRI:102247")
STANDARD = repo_root() / "data" / "standard"


class _FakeRef:
    """Minimal stand-in for an open raster, for the pure alignment maths."""

    def __init__(self, origin=(1000.0, 5000.0), res=5.0, shape=(100, 100)):
        self.transform = Affine(res, 0.0, origin[0], 0.0, -res, origin[1])
        self.shape = shape
        self.crs = CRS


def test_target_keeps_reference_phase() -> None:
    """Every target edge must fall on a reference edge, else the extension resamples."""
    ref = _FakeRef(origin=(1000.0, 5000.0), res=5.0)
    t = aligned_grid((317.0, 1234.0, 2611.0, 6789.0), ref)
    off_x = (t.transform.c - ref.transform.c) / 5.0
    off_y = (t.transform.f - ref.transform.f) / 5.0
    assert abs(off_x - round(off_x)) < 1e-9
    assert abs(off_y - round(off_y)) < 1e-9


def test_target_contains_requested_bounds() -> None:
    ref = _FakeRef()
    want = (317.0, 1234.0, 2611.0, 6789.0)
    left, bottom, right, top = aligned_grid(want, ref).bounds
    assert left <= want[0] and bottom <= want[1]
    assert right >= want[2] and top >= want[3]


def test_grid_is_snapped_outward_not_inward() -> None:
    """Snapping must never clip the requested area, even by one cell."""
    ref = _FakeRef(origin=(0.0, 100.0), res=10.0)
    t = aligned_grid((1.0, 1.0, 99.0, 99.0), ref)
    left, bottom, right, top = t.bounds
    assert left <= 1.0 and bottom <= 1.0 and right >= 99.0 and top >= 99.0


def test_anisotropic_reference_is_refused_without_explicit_res() -> None:
    """surf_2m has 2.382 x 1.288 m pixels; guessing a cell size from that is wrong."""
    ref = _FakeRef()
    ref.transform = Affine(2.382, 0.0, 0.0, 0.0, -1.288, 0.0)
    with pytest.raises(ValueError, match="anisotropic"):
        aligned_grid((0.0, 0.0, 10.0, 10.0), ref)
    assert aligned_grid((0.0, 0.0, 10.0, 10.0), ref, res=5.0).res == (5.0, 5.0)


def _write_vectors(tmp_path):
    import geopandas as gpd
    from shapely.geometry import box

    ice = gpd.GeoDataFrame(geometry=[box(0, 0, 100, 100)], crs=CRS)
    # a lake wholly inside the ice, to prove lakes win the overlap
    lakes = gpd.GeoDataFrame({"class": [2]}, geometry=[box(20, 20, 60, 60)], crs=CRS)
    ip, lp = tmp_path / "ice.gpkg", tmp_path / "lakes.gpkg"
    ice.to_file(ip, driver="GPKG")
    lakes.to_file(lp, driver="GPKG")
    return ip, lp


def test_lakes_overwrite_ice(tmp_path) -> None:
    target = TargetGrid(Affine(10.0, 0.0, 0.0, 0.0, -10.0, 100.0), (10, 10), CRS)
    ip, lp = _write_vectors(tmp_path)
    cls = rasterize_classes(target, ice_path=ip, lakes_path=lp)
    assert set(np.unique(cls)) <= {0, 1, 2}
    assert (cls == 2).sum() > 0, "lake was not burned"
    assert (cls == 1).sum() > 0, "ice was fully overwritten"
    # the lake centre must be lake, not ice
    assert cls[target.shape[0] // 2, target.shape[1] // 2] == 2


def test_rasterize_without_lakes_is_ice_and_land(tmp_path) -> None:
    target = TargetGrid(Affine(10.0, 0.0, 0.0, 0.0, -10.0, 100.0), (10, 10), CRS)
    ip, _ = _write_vectors(tmp_path)
    cls = rasterize_classes(target, ice_path=ip, lakes_path=None)
    assert set(np.unique(cls)) <= {0, 1}


@pytest.mark.skipif(
    not (STANDARD / "classes.tif").exists(),
    reason="extended grid not built; run gpbreach.cascade.build_grid",
)
def test_extended_classes_reproduce_the_original() -> None:
    """The rebuild must extend the owner's classification, never reinterpret it."""
    original = repo_root() / "data" / "processed" / "classes.tif"
    if not original.exists():
        pytest.skip("original classes.tif not present")
    with rasterio.open(STANDARD / "classes.tif") as src:
        target = TargetGrid(src.transform, src.shape, src.crs)
        new = src.read(1)
    stats = verify_against(new, target, original)
    assert stats["status"] == "compared"
    assert stats["agree_pct"] > 99.9, f"only {stats['agree_pct']}% agreement with the original"


@pytest.mark.skipif(
    not (STANDARD / "grid_manifest.json").exists(), reason="extended grid not built"
)
def test_extended_grid_captures_all_three_lakes() -> None:
    import json

    m = json.loads((STANDARD / "grid_manifest.json").read_text())
    cov = m.get("lake_coverage", {})
    assert cov, "manifest has no lake coverage"
    for name, entry in cov.items():
        if "captured_pct" in entry:
            assert entry["captured_pct"] > 99.0, f"{name} only {entry['captured_pct']}% captured"
