"""calibrate.py — inverting an observed breach for the flotation factor."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.cascade.calibrate import Calibration, h_pass_at, implied_k, summarise

CODES = {"land": 0, "ice": 1, "lgp": 2, "gpl": 3, "alsek": 4}


def _dam(n=50, seed=3):
    cl = np.full((n, n), CODES["ice"], dtype=np.uint8)
    cl[2:9, 2:9] = CODES["gpl"]
    cl[n - 9 : n - 2, n - 9 : n - 2] = CODES["lgp"]
    rng = np.random.default_rng(seed)
    zb = np.full((n, n), -50.0) + rng.normal(0, 8, (n, n))
    zs = zb + 200.0 + rng.normal(0, 15, (n, n))
    zs[n // 2 - 2 : n // 2 + 2, :] += 90.0
    return zb, zs, cl, {CODES["gpl"]: 10.0, CODES["lgp"]: 120.0}


@pytest.mark.parametrize("k_true", [0.70, 0.85, 1.00])
def test_round_trip_recovers_the_flotation_factor(k_true: float) -> None:
    """The inversion must return the k that generated the pass head."""
    zb, zs, cl, lv = _dam()
    hp = h_pass_at(k_true, zb, zs, cl, CODES, lv)
    c = implied_k(zb, zs, cl, CODES, lv, h_lake=hp, refine=8)
    assert c.status == "ok"
    assert c.implied_k == pytest.approx(k_true, abs=0.01)


def test_pass_head_increases_with_k() -> None:
    """Monotonicity is what makes the one-dimensional root find valid."""
    zb, zs, cl, lv = _dam()
    hs = [h_pass_at(k, zb, zs, cl, CODES, lv) for k in (0.6, 0.8, 1.0)]
    assert hs[0] < hs[1] < hs[2]


def test_lake_above_every_barrier_is_flagged() -> None:
    zb, zs, cl, lv = _dam()
    c = implied_k(zb, zs, cl, CODES, lv, h_lake=1e4, n_table=5, refine=2)
    assert c.status == "above_range"
    assert np.isnan(c.implied_k)


def test_lake_below_every_barrier_is_flagged() -> None:
    zb, zs, cl, lv = _dam()
    c = implied_k(zb, zs, cl, CODES, lv, h_lake=-1e4, n_table=5, refine=2)
    assert c.status == "below_range"
    assert np.isnan(c.implied_k)


def test_k_above_one_carries_a_warning() -> None:
    """A lake past flotation without draining points at the inputs, not the ice."""
    zb, zs, cl, lv = _dam()
    hp = h_pass_at(1.15, zb, zs, cl, CODES, lv)
    c = implied_k(zb, zs, cl, CODES, lv, h_lake=hp, k_min=0.8, k_max=1.4, refine=8)
    if c.status == "ok" and c.implied_k > 1.0:
        assert "flotation" in c.note


def test_kappa_is_reported_alongside_k() -> None:
    zb, zs, cl, lv = _dam()
    hp = h_pass_at(0.9, zb, zs, cl, CODES, lv)
    c = implied_k(zb, zs, cl, CODES, lv, h_lake=hp, refine=6)
    assert c.implied_kappa == pytest.approx(c.implied_k * 0.917, rel=1e-6)


def test_summarise_pools_events_into_a_prior() -> None:
    res = {
        "a": Calibration(0.90, 0.825, 100.0, 110.0, (0.89, 0.91)),
        "b": Calibration(0.80, 0.734, 90.0, 105.0, (0.79, 0.81)),
        "c": Calibration(float("nan"), float("nan"), 0.0, 0.0, (0, 0), status="no_route"),
    }
    df = summarise(res)
    assert len(df) == 3
    assert df.attrs["n_events"] == 2
    assert df.attrs["k_mean"] == pytest.approx(0.85)
    assert df.attrs["k_min"] == pytest.approx(0.80)
