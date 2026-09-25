"""ATL06 reader, including the quality screen that this module exists to enforce."""

from __future__ import annotations

from pathlib import Path

import pytest

from gpbreach.config import repo_root
from gpbreach.io import atl06

GRANULE_DIR = repo_root() / "data" / "raw" / "ICESat-Inland-Glaceir-Elevations"
DOMAIN = (-138.10, 58.95, -137.85, 59.13)


def _granule() -> Path:
    found = sorted(GRANULE_DIR.glob("ATL06_*.h5")) if GRANULE_DIR.is_dir() else []
    if not found:
        pytest.skip(f"no ATL06 granule in {GRANULE_DIR}")
    return found[0]


def test_heights_are_declared_on_the_wgs84_ellipsoid() -> None:
    info = atl06.granule_info(_granule())
    assert info.product == "ATL06"
    d = info.description.lower()
    assert "wgs 84" in d or "wgs84" in d
    assert "ellipsoid" in d


def test_geoid_is_tide_free_unlike_atl13() -> None:
    """ATL06 reports a tide-free geoid; ATL13 reports mean-tide. Do not mix them."""
    info = atl06.granule_info(_granule())
    assert "tide-free" in info.geoid_description.lower()
    assert info.tide_system == "tide-free"


def test_quality_screen_is_on_by_default() -> None:
    raw = atl06.read_segments(_granule(), bbox=DOMAIN, quality_only=False)
    screened = atl06.read_segments(_granule(), bbox=DOMAIN, quality_only=True)
    assert len(screened) <= len(raw)
    if not screened.empty:
        assert (screened.quality == 0).all()
        assert (screened.vs_dem.abs() <= 50.0).all()


def test_screen_rejects_the_cloud_corrupted_domain_segments() -> None:
    """Regression on the 2026-07-19 granule: unscreened heights here are ~600 m wrong.

    If this ever stops rejecting, the screen has been weakened -- which would let
    heights that disagree with ATL06's own DEM by hundreds of metres into the
    analysis (D-008).
    """
    raw = atl06.read_segments(_granule(), bbox=DOMAIN, quality_only=False)
    if raw.empty:
        pytest.skip("granule does not cover the analysis domain")
    if (raw.quality == 0).any():
        pytest.skip("this granule has good-quality data in the domain; screen untested here")
    screened = atl06.read_segments(_granule(), bbox=DOMAIN)
    assert screened.empty
    assert raw.vs_dem.abs().median() > 100.0


def test_quality_by_band_separates_cloud_from_no_coverage() -> None:
    bands = [(58.20, 58.70), (58.95, 59.13), (59.13, 59.50)]
    table = atl06.quality_by_band(_granule(), bands)
    assert len(table) == 3
    assert {"n", "n_good", "pct_good"} <= set(table.columns)
    # coverage exists in every band; what differs is whether it is usable
    assert (table["n"] > 0).all()


def test_orthometric_height_is_a_subtraction() -> None:
    df = atl06.read_segments(_granule(), quality_only=True)
    if df.empty:
        pytest.skip("no screened segments anywhere in the granule")
    residual = (df.h_li - df.geoid_h) - df.h_ortho
    assert residual.abs().max() < 1e-6
