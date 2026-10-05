"""seal.py: toy grids with known answers, plus the tie back to Phase 1."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.cascade.seal import (
    competing_outlets,
    head_surface,
    k_from_kappa,
    kappa_from_k,
    masked_smooth,
    overtopping_search,
    pass_search,
)
from gpbreach.components.criterion import flotation_field
from gpbreach.constants import LEGACY_FLOTATION_FRACTION, RHO_ICE, RHO_WATER
from gpbreach.saddle import saddle_threshold

LAND, ICE, LGP, GPL, ALSEK = 0, 1, 2, 3, 4


def _corridor(barrier_head: float = 50.0):
    """GPL on the left, LGP on the right, one ice corridor with a known high point."""
    classes = np.full((5, 9), LAND, dtype=np.uint8)
    classes[2, 0] = GPL
    classes[2, 8] = LGP
    classes[2, 1:8] = ICE
    head = np.full((5, 9), np.inf)
    head[2, 1:8] = 10.0
    head[2, 4] = barrier_head  # the pass
    head[2, 0] = 5.0  # GPL level
    head[2, 8] = -np.inf  # set by the search anyway
    return head, classes


def test_pass_head_is_the_barrier() -> None:
    head, classes = _corridor(50.0)
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    assert r.reached_target
    assert r.h_pass == pytest.approx(50.0)
    assert r.pass_cell == (2, 4)


def test_lower_of_two_routes_wins() -> None:
    classes = np.full((5, 9), LAND, dtype=np.uint8)
    classes[1, 1:8] = ICE
    classes[3, 1:8] = ICE
    classes[1:4, 0] = GPL
    classes[1:4, 8] = LGP
    head = np.full((5, 9), np.inf)
    head[1, 1:8] = 10.0
    head[3, 1:8] = 10.0
    head[1, 4] = 70.0  # high route
    head[3, 4] = 30.0  # low route -- the answer
    head[1:4, 0] = 5.0
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    assert r.h_pass == pytest.approx(30.0)
    assert r.pass_cell == (3, 4)


def test_land_wall_blocks_completely() -> None:
    head, classes = _corridor()
    classes[2, 4] = LAND
    head[2, 4] = np.inf
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    assert not r.reached_target
    assert np.isinf(r.h_pass)
    assert r.pass_cell is None


def test_missing_data_is_a_barrier_not_a_hole() -> None:
    """An ice cell without bed or surface must not be routed through."""
    zb = np.zeros((3, 5))
    zs = np.full((3, 5), 100.0)
    classes = np.full((3, 5), LAND, dtype=np.uint8)
    classes[1, :] = ICE
    classes[1, 0] = GPL
    classes[1, 4] = LGP
    zs[1, 2] = np.nan
    head, info = head_surface(zb, zs, classes, {GPL: 5.0, LGP: 50.0}, k=1.0)
    assert np.isinf(head[1, 2])
    assert info["ice_cells_missing_data"] == 1
    assert not pass_search(head, classes, GPL, LGP, connectivity=4).reached_target


def test_eight_connectivity_crosses_a_diagonal_four_cannot() -> None:
    classes = np.full((3, 3), LAND, dtype=np.uint8)
    classes[0, 0] = GPL
    classes[2, 2] = LGP
    classes[1, 1] = ICE
    head = np.full((3, 3), np.inf)
    head[0, 0] = 1.0
    head[1, 1] = 10.0
    assert pass_search(head, classes, GPL, LGP, 8).h_pass == pytest.approx(10.0)
    assert not pass_search(head, classes, GPL, LGP, 4).reached_target


# --- head surface ---------------------------------------------------------


def test_head_formula_matches_the_spec() -> None:
    zb = np.array([[0.0, -100.0]])
    zs = np.array([[500.0, 400.0]])
    classes = np.array([[ICE, ICE]], dtype=np.uint8)
    head, _ = head_surface(zb, zs, classes, {}, k=1.0)
    expected = zb + 1.0 * (RHO_ICE / RHO_WATER) * (zs - zb)
    assert head == pytest.approx(expected)


def test_negative_thickness_is_clipped_and_counted() -> None:
    zb = np.array([[100.0]])
    zs = np.array([[50.0]])  # surface below bed
    classes = np.array([[ICE]], dtype=np.uint8)
    head, info = head_surface(zb, zs, classes, {}, k=1.0)
    assert head[0, 0] == pytest.approx(100.0)  # max(zs-zb,0) = 0
    assert info["negative_thickness_cells"] == 1


def test_lake_cells_take_their_level() -> None:
    zb = np.zeros((1, 3))
    zs = np.full((1, 3), 100.0)
    classes = np.array([[LGP, ICE, GPL]], dtype=np.uint8)
    head, _ = head_surface(zb, zs, classes, {LGP: 117.3, GPL: 27.6}, k=1.0)
    assert head[0, 0] == pytest.approx(117.3)
    assert head[0, 2] == pytest.approx(27.6)


# --- the tie to Phase 1 ---------------------------------------------------


def test_head_surface_equals_the_flotation_field() -> None:
    """h_phi == Phi_f with f = k*rho_i/rho_w. The two halves of the project agree."""
    rng = np.random.default_rng(0)
    zb = rng.uniform(-400, 200, (60, 60))
    zs = zb + rng.uniform(0, 600, (60, 60))
    classes = np.full(zb.shape, ICE, dtype=np.uint8)
    for k in (1.0, 0.95, 0.90):
        head, _ = head_surface(zb, zs, classes, {}, k=k)
        assert head == pytest.approx(flotation_field(zs, zb, kappa_from_k(k)), abs=1e-9)


def test_k_and_f_are_inverses() -> None:
    for k in (0.8, 0.9, 1.0):
        assert k_from_kappa(kappa_from_k(k)) == pytest.approx(k)
    assert k_from_kappa(LEGACY_FLOTATION_FRACTION) == pytest.approx(0.991375, abs=1e-6)


def test_pass_search_agrees_with_saddle_threshold() -> None:
    """The priority flood must return exactly what Phase 1's saddle search does.

    Different algorithms (flood vs. bisection over sorted levels), same answer --
    this is what stops the cascade and baseline halves from drifting apart.
    """
    rng = np.random.default_rng(1)
    zb = rng.uniform(-200, 100, (40, 40))
    zs = zb + rng.uniform(0, 400, (40, 40))
    classes = np.full(zb.shape, ICE, dtype=np.uint8)
    classes[0, 0] = GPL
    classes[-1, -1] = LGP
    k = 0.95
    head, _ = head_surface(zb, zs, classes, {GPL: -1e6, LGP: 1e6}, k=k)

    r = pass_search(head, classes, GPL, LGP, connectivity=8)
    field = flotation_field(zs, zb, kappa_from_k(k))
    expected = saddle_threshold(field, (0, 0), (zb.shape[0] - 1, zb.shape[1] - 1), 8)
    assert r.h_pass == pytest.approx(expected, abs=1e-9)


# --- smoothing ------------------------------------------------------------


def test_zero_sigma_is_a_no_op_on_the_mask() -> None:
    z = np.arange(25.0).reshape(5, 5)
    mask = np.ones((5, 5), dtype=bool)
    mask[0, 0] = False
    out = masked_smooth(z, mask, 0.0, 5.0)
    assert np.isnan(out[0, 0])
    assert out[mask] == pytest.approx(z[mask])


def test_smoothing_does_not_pull_in_offmask_zeros() -> None:
    """Normalised convolution: a constant field stays constant near the mask edge."""
    z = np.full((9, 9), 100.0)
    mask = np.zeros((9, 9), dtype=bool)
    mask[3:6, 3:6] = True
    out = masked_smooth(z, mask, 10.0, 5.0)
    assert out[mask] == pytest.approx(100.0, abs=1e-6)
    assert np.isnan(out[0, 0])


# --- other searches -------------------------------------------------------


def test_overtopping_uses_surface_not_head() -> None:
    zs = np.full((3, 5), 200.0)
    classes = np.full((3, 5), LAND, dtype=np.uint8)
    classes[1, :] = ICE
    classes[1, 0] = GPL
    classes[1, 4] = LGP
    zs[1, 2] = 300.0
    r = overtopping_search(zs, classes, {GPL: 5.0, LGP: 50.0}, GPL, LGP, 4)
    assert r.h_pass == pytest.approx(300.0)


def test_competing_outlets_ranks_receivers() -> None:
    classes = np.full((5, 9), LAND, dtype=np.uint8)
    classes[2, 4] = LGP
    classes[2, 0] = GPL
    classes[2, 8] = ALSEK
    classes[2, 1:4] = ICE
    classes[2, 5:8] = ICE
    head = np.full((5, 9), np.inf)
    head[2, 1:4] = 60.0  # toward GPL
    head[2, 5:8] = 20.0  # toward Alsek -- cheaper
    head[2, 0] = 5.0
    head[2, 8] = 5.0
    out = competing_outlets(head, classes, {"gpl": GPL, "alsek": ALSEK}, LGP, 4)
    assert out["gpl"]["h_pass"] == pytest.approx(60.0)
    assert out["alsek"]["h_pass"] == pytest.approx(20.0)


def test_competing_outlets_reports_absent_lake() -> None:
    classes = np.full((3, 3), LAND, dtype=np.uint8)
    classes[1, 0] = GPL
    classes[1, 2] = LGP
    head = np.full((3, 3), np.inf)
    head[1, :] = 10.0
    out = competing_outlets(head, classes, {"alsek": ALSEK}, LGP, 4)
    assert out["alsek"]["h_pass"] is None
    assert "absent" in out["alsek"]["status"]
