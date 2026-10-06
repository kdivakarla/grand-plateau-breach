"""lakes.py: the analytic cone from the spec's test table, plus the sill logic."""

from __future__ import annotations

import numpy as np
import pytest
from affine import Affine

from gpbreach.cascade.lakes import (
    conic_bottom,
    connected_stage_area,
    detect_spill,
    drainage_floor,
    interpolate_bathymetry,
    released_volume,
    stage_curve,
    vertical_wall_bottom,
    volume_at,
)

RES = 2.0
CELL = RES * RES
T = Affine(RES, 0.0, 0.0, 0.0, -RES, 0.0)


def _cone(n=420, slope=0.5, z0=0.0):
    """Inverted cone: bottom(r) = z0 + slope*r. At stage h, A = pi*r^2, V = A*(h-z0)/3."""
    c = n / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.hypot((xx - c) * RES, (yy - c) * RES)
    bottom = z0 + slope * r
    mask = r <= (n / 2.0 - 2) * RES
    return bottom, mask, slope, z0


def test_cone_area_and_volume_match_analytic() -> None:
    """Spec pass condition: A(h) and V(h) within 1% of the analytic cone."""
    bottom, mask, slope, z0 = _cone()
    for h in (40.0, 70.0, 100.0):
        curve = stage_curve(bottom, mask, [h], CELL)
        r_h = (h - z0) / slope
        a_exact = np.pi * r_h**2
        v_exact = a_exact * (h - z0) / 3.0
        a = curve["A_m2"].iloc[0]
        v = curve["V_m3"].iloc[0]
        assert abs(a - a_exact) / a_exact < 0.01, f"area off at h={h}: {a:.0f} vs {a_exact:.0f}"
        assert abs(v - v_exact) / v_exact < 0.01, f"volume off at h={h}: {v:.0f} vs {v_exact:.0f}"


def test_volume_is_monotonic_in_stage() -> None:
    bottom, mask, _, _ = _cone()
    curve = stage_curve(bottom, mask, np.linspace(10, 100, 10), CELL)
    assert (np.diff(curve["V_m3"].to_numpy()) > 0).all()
    assert (np.diff(curve["A_m2"].to_numpy()) >= 0).all()


def test_empty_at_or_below_the_floor() -> None:
    bottom, mask, _, z0 = _cone()
    curve = stage_curve(bottom, mask, [z0 - 1.0], CELL)
    assert curve["V_m3"].iloc[0] == 0.0
    assert curve["A_m2"].iloc[0] == 0.0


# --- the drainage floor ---------------------------------------------------


def test_bed_sill_controls_when_it_sits_above_the_receiving_lake() -> None:
    h, who = drainage_floor(h_gpl=27.6, bed_sill=58.71)
    assert h == pytest.approx(58.71)
    assert who == "bed_sill"


def test_receiving_lake_controls_when_the_sill_is_lower() -> None:
    h, who = drainage_floor(h_gpl=27.6, bed_sill=10.0)
    assert h == pytest.approx(27.6)
    assert who == "receiving_lake"


def test_sill_withholds_water() -> None:
    """A sill above the receiving lake must reduce V_w, and the result says by how much."""
    mask = np.ones((60, 60), bool)
    bottom = vertical_wall_bottom(mask, 0.0)
    curve = stage_curve(bottom, mask, np.linspace(0, 120, 121), CELL)
    with_sill = released_volume(curve, h_start=117.3, h_gpl=27.6, bed_sill=58.71)
    assert with_sill["controlled_by"] == "bed_sill"
    assert with_sill["Vw_m3"] < with_sill["Vw_if_no_sill_Mm3"] * 1e6
    assert with_sill["sill_withholds_Mm3"] > 0
    # prism: V_w should be area x drawdown
    area = mask.sum() * CELL
    assert with_sill["Vw_m3"] == pytest.approx(area * (117.3 - 58.71), rel=0.01)


def test_released_volume_never_negative() -> None:
    mask = np.ones((20, 20), bool)
    curve = stage_curve(vertical_wall_bottom(mask, 0.0), mask, np.linspace(0, 50, 51), CELL)
    out = released_volume(curve, h_start=10.0, h_gpl=27.6, bed_sill=40.0)
    assert out["Vw_m3"] == 0.0


# --- bounds ---------------------------------------------------------------


def test_prism_bounds_the_cone_from_above() -> None:
    """Vertical walls must give more volume than any narrowing basin."""
    bottom, mask, _, _ = _cone()
    h = 80.0
    cone_v = stage_curve(bottom, mask, [h], CELL)["V_m3"].iloc[0]
    floor = np.nanmin(bottom[mask])
    prism_v = stage_curve(vertical_wall_bottom(mask, floor), mask, [h], CELL)["V_m3"].iloc[0]
    assert prism_v > cone_v


def test_conic_bottom_is_shallower_at_the_margin() -> None:
    mask = np.zeros((80, 80), bool)
    mask[10:70, 10:70] = True
    b = conic_bottom(mask, h_lake=100.0, max_depth=50.0)
    # the shallowest cell is one cell in from the shore, so slightly below h_lake
    assert np.nanmax(b[mask]) == pytest.approx(100.0, abs=2.0)  # shore
    assert np.nanmin(b[mask]) == pytest.approx(50.0, abs=2.0)  # centre


# --- bathymetry gridding --------------------------------------------------


def test_interpolated_bathymetry_honours_the_soundings() -> None:
    mask = np.ones((40, 40), bool)
    xy = np.array([[10.0, -10.0], [70.0, -10.0], [10.0, -70.0], [70.0, -70.0], [40.0, -40.0]])
    depths = np.array([10.0, 10.0, 10.0, 10.0, 50.0])
    b = interpolate_bathymetry(xy, depths, mask, h_lake=100.0, transform=T)
    assert np.isfinite(b[mask]).all()
    assert np.nanmin(b) >= 100.0 - depths.max() - 1e-6
    assert np.nanmax(b) <= 100.0 - depths.min() + 1e-6


# --- connected rise -------------------------------------------------------


def test_connected_rise_excludes_a_separate_basin() -> None:
    dem = np.full((40, 40), 100.0)
    dem[5:15, 5:15] = 0.0  # the lake
    dem[25:35, 25:35] = 0.0  # a separate pit, below the stage but not connected
    seed = np.zeros((40, 40), bool)
    seed[5:15, 5:15] = True
    out = connected_stage_area(dem, seed, [10.0], CELL, base_level=0.0)
    assert out["A_m2"].iloc[0] == pytest.approx(100 * CELL, rel=0.01)


def test_connected_rise_grows_with_stage() -> None:
    dem = np.full((40, 40), 50.0)
    dem[18:22, 18:22] = 0.0
    dem[15:25, 15:25] = np.minimum(dem[15:25, 15:25], 20.0)
    seed = np.zeros((40, 40), bool)
    seed[18:22, 18:22] = True
    out = connected_stage_area(dem, seed, [10.0, 30.0], CELL, base_level=0.0)
    assert out["A_m2"].iloc[1] > out["A_m2"].iloc[0]


def test_volume_at_interpolates() -> None:
    import pandas as pd

    c = pd.DataFrame({"h": [0.0, 10.0], "V_m3": [0.0, 100.0]})
    assert volume_at(c, 5.0) == pytest.approx(50.0)


def test_masked_lake_still_contributes_volume() -> None:
    """The surface DEM masks out lakes; the lake must not vanish from the volume."""
    dem = np.full((40, 40), 100.0)
    seed = np.zeros((40, 40), bool)
    seed[10:30, 10:30] = True
    dem[seed] = np.nan  # lakes masked, as in surf_2018
    out = connected_stage_area(dem, seed, [5.0], CELL, base_level=0.0)
    expected = seed.sum() * CELL * 5.0  # area x rise
    assert out["V_m3"].iloc[0] == pytest.approx(expected, rel=0.01)


def test_sill_withholding_is_unknown_without_deep_bathymetry() -> None:
    mask = np.ones((20, 20), bool)
    curve = stage_curve(
        vertical_wall_bottom(mask, 58.71), mask, np.linspace(58.71, 117.3, 60), CELL
    )
    out = released_volume(curve, 117.3, h_gpl=27.6, bed_sill=58.71)
    assert out["sill_withholds_Mm3"] is None
    assert "bathymetry" in out["sill_withholds_note"]


def test_detect_spill_finds_the_step() -> None:
    """A lake topping a divide into a second basin must be flagged, not smoothed over."""
    dem = np.full((60, 60), 100.0)
    dem[10:20, 10:20] = 0.0  # main basin
    dem[40:55, 40:55] = 0.0  # adjoining basin, deep
    dem[20:40, 20:40] = 12.0  # the divide between them
    seed = np.zeros((60, 60), bool)
    seed[10:20, 10:20] = True
    curve = connected_stage_area(dem, seed, np.arange(0.0, 25.0, 1.0), CELL, base_level=0.0)
    spill = detect_spill(curve)
    assert spill is not None
    assert spill["area_after_km2"] > spill["area_before_km2"]
    assert 10.0 <= spill["spill_stage_m"] <= 16.0


def test_detect_spill_returns_none_for_a_smooth_basin() -> None:
    bottom, mask, _, _ = _cone()
    curve = stage_curve(bottom, mask, np.linspace(20, 100, 40), CELL)
    assert detect_spill(curve) is None
