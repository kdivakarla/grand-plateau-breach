"""project_surface.py — the three thinning models."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.cascade.project_surface import accelerating, by_elevation, project, uniform


def test_uniform_shifts_only_the_ice() -> None:
    zs = np.array([[100.0, 200.0, 300.0]])
    ice = np.array([[True, False, True]])
    out = uniform(zs, ice, -9.05, 10.0)
    assert out[0, 0] == pytest.approx(100.0 - 90.5)
    assert out[0, 1] == pytest.approx(200.0), "non-ice must not be thinned"
    assert out[0, 2] == pytest.approx(300.0 - 90.5)


def test_uniform_shift_is_constant_so_the_pass_cannot_move() -> None:
    """The reason projected surfaces are unnecessary for a uniform rate."""
    rng = np.random.default_rng(0)
    zs = rng.uniform(100, 500, (40, 40))
    ice = np.ones(zs.shape, bool)
    out = uniform(zs, ice, -9.05, 7.0)
    diff = out - zs
    assert diff.std() == pytest.approx(0.0, abs=1e-12)


def test_by_elevation_applies_the_right_band() -> None:
    zs = np.array([[150.0, 250.0, 350.0]])
    ice = np.ones(zs.shape, bool)
    bands = [(100.0, 200.0, -5.0), (200.0, 300.0, -10.0), (300.0, 400.0, -2.0)]
    out = by_elevation(zs, ice, bands, 1.0)
    assert out[0, 0] == pytest.approx(145.0)
    assert out[0, 1] == pytest.approx(240.0)
    assert out[0, 2] == pytest.approx(348.0)


def test_by_elevation_falls_back_outside_every_band() -> None:
    zs = np.array([[50.0]])
    ice = np.ones(zs.shape, bool)
    out = by_elevation(zs, ice, [(100.0, 200.0, -5.0)], 1.0, default_dhdt=-1.0)
    assert out[0, 0] == pytest.approx(49.0)


def test_annual_stepping_lets_a_cell_change_band() -> None:
    """A cell thinning out of its band must pick up the next band's rate."""
    zs = np.array([[205.0]])
    ice = np.ones(zs.shape, bool)
    bands = [(100.0, 200.0, -1.0), (200.0, 300.0, -10.0)]
    stepped = project(zs, ice, [2020], 2018, "by_elevation", bands=bands)[2020]
    one_jump = by_elevation(zs, ice, bands, 2.0)
    assert stepped[0, 0] > one_jump[0, 0], "stepping must slow once out of the fast band"


def test_accelerating_matches_the_integral() -> None:
    zs = np.array([[100.0]])
    ice = np.ones(zs.shape, bool)
    out = accelerating(zs, ice, dhdt0=-9.0, accel=-0.5, dt_years=4.0)
    assert out[0, 0] == pytest.approx(100.0 + (-9.0 * 4 + 0.5 * -0.5 * 16))


def test_accelerating_thins_faster_than_constant() -> None:
    zs = np.array([[500.0]])
    ice = np.ones(zs.shape, bool)
    acc = project(zs, ice, [2030], 2018, "accelerating", dhdt0=-9.05, accel=-0.3)[2030]
    flat = project(zs, ice, [2030], 2018, "uniform", dhdt=-9.05)[2030]
    assert acc[0, 0] < flat[0, 0]


def test_project_returns_every_requested_year() -> None:
    zs = np.full((5, 5), 300.0)
    ice = np.ones(zs.shape, bool)
    out = project(zs, ice, [2025, 2030, 2020], 2018, "uniform", dhdt=-9.05)
    assert sorted(out) == [2020, 2025, 2030]
    assert out[2030].mean() < out[2025].mean() < out[2020].mean()


def test_unknown_model_rejected() -> None:
    with pytest.raises(ValueError, match="unknown model"):
        project(np.zeros((2, 2)), np.ones((2, 2), bool), [2020], 2018, "magic")
