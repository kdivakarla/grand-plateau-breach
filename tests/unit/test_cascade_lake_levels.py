"""lake_levels.py — robust aggregation and scale-aware buffering."""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Point

from gpbreach.cascade.lake_levels import adaptive_buffer, build_timeseries, robust_level


def test_median_ignores_a_few_wild_points() -> None:
    """An iceberg or a shoreline return must not move the level."""
    h = np.concatenate(
        [np.full(200, 30.0) + np.random.default_rng(0).normal(0, 0.05, 200), [55.0, 61.0, -12.0]]
    )
    out = robust_level(h)
    assert out["level_m"] == pytest.approx(30.0, abs=0.05)
    assert out["n_rejected"] == 3
    assert out["n"] == 200


def test_scatter_is_reported_not_hidden() -> None:
    rng = np.random.default_rng(1)
    tight = robust_level(np.full(100, 30.0) + rng.normal(0, 0.05, 100))
    loose = robust_level(np.full(100, 30.0) + rng.normal(0, 3.0, 100))
    assert tight["mad_m"] < 0.2
    assert loose["mad_m"] > 1.0


def test_empty_input_is_not_a_level() -> None:
    out = robust_level(np.array([]))
    assert out["n"] == 0
    assert np.isnan(out["level_m"])


def test_all_identical_points_survive() -> None:
    """MAD is zero here; the filter must not reject everything."""
    out = robust_level(np.full(20, 42.0))
    assert out["n"] == 20
    assert out["level_m"] == pytest.approx(42.0)


def test_buffer_is_eased_for_a_small_lake() -> None:
    """A fixed buffer is not scale-free; small lakes must keep usable area."""
    big = Point(0, 0).buffer(5000.0)  # ~78 km2
    small = Point(0, 0).buffer(1000.0)  # ~3 km2
    big_b = adaptive_buffer(big, -150.0)
    small_b = adaptive_buffer(small, -150.0)
    assert big_b.area / big.area == pytest.approx(0.94, abs=0.02)  # full buffer applied
    assert small_b.area / small.area >= 0.64  # eased to the floor


def test_buffer_never_drops_below_the_area_floor() -> None:
    small = Point(0, 0).buffer(400.0)
    out = adaptive_buffer(small, -350.0, min_area_frac=0.65)
    assert out.area / small.area >= 0.64


def test_zero_buffer_is_a_no_op() -> None:
    poly = Point(0, 0).buffer(1000.0)
    assert adaptive_buffer(poly, 0.0).area == pytest.approx(poly.area)


def test_timeseries_flags_drive_usability() -> None:
    rows = [
        {
            "granule": "a",
            "date": "2022-07-01",
            "lake": "GPL",
            "n": 100,
            "level_m": 30.0,
            "flag": None,
        },
        {
            "granule": "b",
            "date": "2022-08-01",
            "lake": "GPL",
            "n": 2,
            "level_m": 31.0,
            "flag": "only 2 points",
        },
    ]
    df = build_timeseries(rows)
    assert df["usable"].tolist() == [True, False]
    assert df["doy"].iloc[0] == 182
