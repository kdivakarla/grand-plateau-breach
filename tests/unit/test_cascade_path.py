"""path.py: toy grids with known routes, per the spec's test table."""

from __future__ import annotations

import numpy as np
import pytest
from affine import Affine

from gpbreach.cascade.path import (
    monotonicity,
    profile,
    refine_downstream,
    steepest_descent,
    summarize,
    trace_path,
)
from gpbreach.cascade.seal import pass_search

LAND, ICE, LGP, GPL = 0, 1, 2, 3
T = Affine(5.0, 0.0, 0.0, 0.0, -5.0, 0.0)


def _corridor():
    """GPL left, LGP right, one ice corridor with the pass at column 4."""
    classes = np.full((5, 9), LAND, dtype=np.uint8)
    classes[2, 0] = GPL
    classes[2, 8] = LGP
    classes[2, 1:8] = ICE
    head = np.full((5, 9), np.inf)
    #        col:   1     2     3     4(pass) 5     6     7
    head[2, 1:8] = np.array([10.0, 15.0, 20.0, 50.0, 45.0, 42.0, 40.0])
    head[2, 0] = 5.0
    return head, classes


def test_path_runs_from_lgp_to_gpl() -> None:
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    assert tuple(p[0]) == (2, 8), "path must start in LGP"
    assert classes[tuple(p[-1])] == GPL, "path must end in GPL"


def test_path_is_contiguous() -> None:
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    steps = np.abs(np.diff(p.astype(int), axis=0)).max(axis=1)
    assert (steps == 1).all(), "consecutive path cells must be neighbours"


def test_head_is_non_increasing_below_the_pass() -> None:
    """Spec pass condition: once over the pass, water only goes downhill."""
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    p = refine_downstream(p, head, classes, r.pass_cell, target_class=GPL)
    idx = np.flatnonzero((p[:, 0] == r.pass_cell[0]) & (p[:, 1] == r.pass_cell[1]))
    below = head[p[int(idx[0]) :, 0], p[int(idx[0]) :, 1]]
    below = below[np.isfinite(below)]
    assert (np.diff(below) <= 1e-9).all(), f"head rises below the pass: {below}"


def test_steepest_descent_reaches_the_target() -> None:
    head, classes = _corridor()
    p = steepest_descent(head, (2, 4), classes, GPL)
    assert classes[tuple(p[-1])] == GPL


def test_steepest_descent_reports_a_pit() -> None:
    classes = np.full((3, 5), LAND, dtype=np.uint8)
    classes[1, :] = ICE
    classes[1, 0] = GPL
    head = np.full((3, 5), np.inf)
    head[1, :] = [5.0, 99.0, 1.0, 50.0, 60.0]  # start at col 3, descend to the pit at col 2
    p = steepest_descent(head, (1, 3), classes, GPL)
    assert classes[tuple(p[-1])] != GPL, "should stall in the pit, not reach GPL"


def test_refine_escapes_pits_via_the_parent_fallback() -> None:
    """Raw head descent stalls in a pit; the parent fallback still reaches GPL."""
    head, classes = _corridor()
    head[2, 2] = 2.0  # carve a pit between the pass and GPL
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])

    raw = refine_downstream(p, head, classes, r.pass_cell, target_class=GPL)
    # refinement failed, so the ORIGINAL path comes back -- never a truncated one
    assert np.array_equal(raw, p)
    assert classes[tuple(raw[-1])] == GPL, "fallback must still end in GPL"

    guided = refine_downstream(p, head, classes, r.pass_cell, target_class=GPL, parent=r.parent)
    assert classes[tuple(guided[-1])] == GPL


def test_refine_leaves_path_alone_if_pass_absent() -> None:
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    same = refine_downstream(p, head, classes, (0, 0), target_class=GPL)
    assert np.array_equal(p, same)


def test_profile_distance_accounts_for_diagonals() -> None:
    p = np.array([[0, 0], [1, 1], [1, 2]])
    zb = np.zeros((3, 3))
    zs = np.full((3, 3), 100.0)
    head = np.full((3, 3), 50.0)
    prof = profile(p, zb, zs, head, T, h_lgp=117.3, h_gpl=27.6, res=5.0)
    assert prof["s"].iloc[1] == pytest.approx(5.0 * np.sqrt(2))
    assert prof["s"].iloc[2] == pytest.approx(5.0 * np.sqrt(2) + 5.0)


def test_profile_effective_pressures() -> None:
    p = np.array([[1, 1]])
    zb = np.zeros((3, 3))
    zs = np.full((3, 3), 100.0)
    head = np.full((3, 3), 150.0)
    prof = profile(p, zb, zs, head, T, h_lgp=117.3, h_gpl=27.6)
    assert prof["N_lake"].iloc[0] == pytest.approx(150.0 - 117.3)
    assert prof["N_gpl"].iloc[0] == pytest.approx(150.0 - 27.6)
    assert prof["H"].iloc[0] == pytest.approx(100.0)


def test_summary_bed_sill_is_the_high_point() -> None:
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    zb = np.zeros((5, 9))
    zb[2, 3] = 42.0  # a bed high point on the route
    zs = np.full((5, 9), 300.0)
    prof = profile(p, zb, zs, head, T, h_lgp=117.3, h_gpl=27.6)
    s = summarize(prof, classes, h_lgp=117.3, pass_cell=r.pass_cell)
    assert s.bed_sill_m == pytest.approx(42.0)


def test_summary_n_modes_differ_and_are_selectable() -> None:
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    zb = np.zeros((5, 9))
    zs = np.full((5, 9), 300.0)
    prof = profile(p, zb, zs, head, T, h_lgp=20.0, h_gpl=5.0)
    at_pass = summarize(prof, classes, 20.0, r.pass_cell, n_mode="at_pass")
    minimum = summarize(prof, classes, 20.0, r.pass_cell, n_mode="min")
    assert at_pass.n_representative_m == at_pass.n_at_pass_m
    assert minimum.n_representative_m == pytest.approx(minimum.n_min_m)
    assert minimum.n_min_m <= at_pass.n_at_pass_m
    with pytest.raises(ValueError, match="n_mode"):
        summarize(prof, classes, 20.0, r.pass_cell, n_mode="nonsense")


def test_summary_lengths_are_ordered() -> None:
    head, classes = _corridor()
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    zb = np.zeros((5, 9))
    zs = np.full((5, 9), 300.0)
    prof = profile(p, zb, zs, head, T, h_lgp=117.3, h_gpl=27.6)
    s = summarize(prof, classes, h_lgp=117.3, pass_cell=r.pass_cell)
    assert 0 < s.l_pass_m <= s.l_margin_m <= s.l_total_m


def test_failed_refinement_returns_the_original_path() -> None:
    """A truncated route would silently corrupt L_pass; refusing is the only safe move."""
    head, classes = _corridor()
    head[2, 2] = 2.0  # pit between the pass and GPL
    r = pass_search(head, classes, GPL, LGP, connectivity=4)
    p = trace_path(r.parent, r.info["best_target_cell"])
    out = refine_downstream(p, head, classes, r.pass_cell, target_class=GPL)
    assert np.array_equal(out, p)


def test_monotonicity_detects_rises() -> None:
    import pandas as pd

    falling = pd.DataFrame({"head": [10.0, 8.0, 5.0, 1.0]})
    assert monotonicity(falling)["monotonic"]
    bumpy = pd.DataFrame({"head": [10.0, 8.0, 9.0, 1.0]})
    m = monotonicity(bumpy)
    assert not m["monotonic"]
    assert m["rising_steps"] == 1
    assert m["largest_rise_m"] == pytest.approx(1.0)
