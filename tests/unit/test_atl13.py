"""ATL13 reader.

The granule lives in ``data/raw/`` and is gitignored, so every test that needs
it skips when it is absent. The datum assertions are the point: they fail if a
future product version changes what the file declares, which is exactly the
change that must not pass unnoticed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from gpbreach.config import repo_root
from gpbreach.io.atl13 import (
    granule_info,
    navd88_from_ellipsoid,
    read_segments,
    summarise_water_bodies,
)

GRANULE_DIR = repo_root() / "data" / "raw" / "ICESat-Inland-lake-elevations"
# Grand Plateau / Alsek region
BBOX = (-139.0, 58.6, -137.2, 59.6)


def _granule() -> Path:
    found = sorted(GRANULE_DIR.glob("ATL13_*.h5")) if GRANULE_DIR.is_dir() else []
    if not found:
        pytest.skip(f"no ATL13 granule in {GRANULE_DIR}")
    return found[0]


def test_granule_identifies_as_atl13() -> None:
    info = granule_info(_granule())
    assert info.product == "ATL13"
    assert "ATL13" in info.doi


def test_water_surface_is_declared_on_the_ellipsoid() -> None:
    """ht_water_surf must say WGS84 ellipsoid -- the basis of every conversion."""
    info = granule_info(_granule())
    desc = info.datums["ht_water_surf"].lower()
    assert "wgs84" in desc.replace("-", "") or "wgs 84" in desc
    assert "ellipsoid" in desc
    assert info.units["ht_water_surf"] == "meters"


def test_orthometric_height_is_egm2008_not_navd88() -> None:
    """The trap this module exists to prevent: ht_ortho is EGM2008."""
    info = granule_info(_granule())
    desc = info.datums["ht_ortho"].upper()
    assert "EGM2008" in desc
    assert "NAVD" not in desc


def test_geoid_is_mean_tide() -> None:
    info = granule_info(_granule())
    assert "mean-tide" in info.datums["segment_geoid"].lower()


def test_ortho_equals_ellipsoid_minus_geoid() -> None:
    """Internal consistency: ht_ortho == ht_water_surf - segment_geoid."""
    df = read_segments(_granule(), bbox=BBOX)
    if df.empty:
        pytest.skip("granule does not cover the Grand Plateau region")
    residual = df["ht_water_surf"] - df["segment_geoid"] - df["ht_ortho"]
    assert np.abs(residual).max() < 0.01


def test_fill_values_become_nan() -> None:
    """ATL13 uses 127 as the int8 fill; it must not survive as a real flag."""
    df = read_segments(_granule(), bbox=BBOX)
    if df.empty:
        pytest.skip("granule does not cover the Grand Plateau region")
    for col in ("ice_flag", "qf_cloud", "qf_ice"):
        if col in df:
            assert not (df[col] == 127).any(), f"{col} still carries the fill value"


def test_bbox_actually_restricts() -> None:
    df = read_segments(_granule(), bbox=BBOX)
    if df.empty:
        pytest.skip("granule does not cover the Grand Plateau region")
    assert df["segment_lat"].between(BBOX[1], BBOX[3]).all()
    assert df["segment_lon"].between(BBOX[0], BBOX[2]).all()


def test_navd88_conversion_is_a_subtraction() -> None:
    assert navd88_from_ellipsoid(32.745, 7.483) == pytest.approx(25.262, abs=1e-3)


def test_summary_reports_the_egm2008_navd88_gap() -> None:
    """The ~1 m gap between the two orthometric systems must be visible."""
    df = read_segments(_granule(), bbox=BBOX)
    if df.empty:
        pytest.skip("granule does not cover the Grand Plateau region")
    table = summarise_water_bodies(df, geoid12_height_m=7.483)
    assert not table.empty
    gap = table["egm2008_minus_navd88_m"].abs()
    # EGM2008 and GEOID12A disagree here; if this ever collapses to zero the
    # conversion has been short-circuited somewhere.
    assert (gap > 0.5).all()
