"""scenarios.py — the caching contract and the propagation chain."""

from __future__ import annotations

import numpy as np
import pytest
from affine import Affine

from gpbreach.cascade.scenarios import Scenario, expand, volume_bracket

CODES = {"land": 0, "ice": 1, "lgp": 2, "gpl": 3, "alsek": 4}
T = Affine(5.0, 0.0, 0.0, 0.0, -5.0, 0.0)


def test_thinning_swap_reuses_the_geometry_solve() -> None:
    """The expensive solve must not depend on the thinning rate."""
    a = Scenario("millan", "s2018", "fast")
    b = Scenario("millan", "s2018", "slow")
    assert a.geometry_key == b.geometry_key
    assert a.name != b.name


@pytest.mark.parametrize(
    "field,value",
    [
        ("bed_id", "iceboost"),
        ("surface_id", "other"),
        ("k", 0.9),
        ("connectivity", 4),
        ("smoothing_m", 50.0),
    ],
)
def test_every_geometry_input_invalidates_the_cache(field: str, value) -> None:
    """Anything that changes h_pass must change the key, or a stale answer is served."""
    base = Scenario("millan", "s2018", "fast")
    alt = Scenario(**{**base.__dict__, field: value})
    assert base.geometry_key != alt.geometry_key, f"{field} did not invalidate the cache"


def test_expand_is_the_full_product() -> None:
    s = expand(["b1", "b2"], ["s1"], ["t1", "t2", "t3"], ks=(0.9, 1.0))
    assert len(s) == 2 * 1 * 3 * 2
    assert len({x.name for x in s}) == len(s)


def _classes(n=40):
    c = np.full((n, n), CODES["ice"], dtype=np.uint8)
    c[2:8, 2:8] = CODES["lgp"]
    c[30:38, 30:38] = CODES["gpl"]
    return c


def test_sill_above_the_lake_means_no_drainage() -> None:
    """A sill above the lake surface must yield V_w = 0, not a negative volume."""
    out = volume_bracket(_classes(), CODES, 25.0, h_lgp=117.3, h_gpl=27.6, bed_sill=185.4)
    assert out["drawdown"] == 0.0
    assert out["lo"] == 0.0 and out["hi"] == 0.0
    assert out["controlled_by"] == "bed_sill"


def test_lower_sill_releases_more_water() -> None:
    high = volume_bracket(_classes(), CODES, 25.0, 117.3, 27.6, bed_sill=90.0)
    low = volume_bracket(_classes(), CODES, 25.0, 117.3, 27.6, bed_sill=40.0)
    assert low["drawdown"] > high["drawdown"]
    assert low["hi"] > high["hi"]


def test_prism_bound_is_above_the_conic_bound() -> None:
    out = volume_bracket(_classes(), CODES, 25.0, 117.3, 27.6, bed_sill=58.71)
    assert out["hi"] > out["lo"] > 0


def test_receiving_lake_controls_when_the_sill_is_low() -> None:
    out = volume_bracket(_classes(), CODES, 25.0, 117.3, 27.6, bed_sill=10.0)
    assert out["controlled_by"] == "receiving_lake"
    assert out["h_final"] == pytest.approx(27.6)
