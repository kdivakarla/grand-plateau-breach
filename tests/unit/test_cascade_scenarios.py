"""scenarios.py — the caching contract and the propagation chain."""

from __future__ import annotations

import numpy as np
import pytest
from affine import Affine

from gpbreach.cascade.scenarios import Scenario, bed_pass, expand, volume_bracket

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


def test_bed_pass_finds_the_low_route_not_the_high_one() -> None:
    """The bed-only pass must ignore where the flotation route happens to go.

    Two corridors join the lakes: one over a 180 m ridge, one over a 20 m ridge.
    The bed-only minimax must return 20 whichever the head would pick.
    """
    n = 20
    cl = np.full((n, n), CODES["land"], dtype=np.uint8)
    cl[0:3, 0:3] = CODES["gpl"]
    cl[17:20, 17:20] = CODES["lgp"]
    bed = np.full((n, n), np.nan)
    bed[0:3, 0:3] = -50.0
    bed[17:20, 17:20] = -50.0
    bed[1, 3:18] = 180.0  # high corridor, across then down
    bed[1:18, 17] = 180.0
    bed[3:18, 1] = 20.0  # low corridor, down then across
    bed[17, 1:18] = 20.0
    assert bed_pass(bed, cl, CODES, connectivity=8) == pytest.approx(20.0)


def test_bed_pass_is_inf_when_no_route_exists() -> None:
    n = 20
    cl = np.full((n, n), CODES["land"], dtype=np.uint8)
    cl[1:4, 1:4] = CODES["gpl"]
    cl[16:19, 16:19] = CODES["lgp"]
    bed = np.full((n, n), np.nan)
    bed[1:4, 1:4] = -10.0
    bed[16:19, 16:19] = -10.0  # no bed anywhere between them
    assert np.isinf(bed_pass(bed, cl, CODES))


def test_the_two_sills_give_different_volumes() -> None:
    """The point of reporting both: they are not interchangeable."""
    cls = _classes()
    along = volume_bracket(cls, CODES, 25.0, 117.3, 27.6, bed_sill=58.71)
    bedonly = volume_bracket(cls, CODES, 25.0, 117.3, 27.6, bed_sill=11.08)
    assert bedonly["drawdown"] > along["drawdown"]
    assert bedonly["hi"] > along["hi"]
    assert along["controlled_by"] == "bed_sill"
    assert bedonly["controlled_by"] == "receiving_lake"


def test_cache_key_is_versioned() -> None:
    """A cached payload that predates bed_pass_m must not be served.

    The key carries a schema version; without it an old cache answers a new
    question with a field missing.
    """
    import hashlib

    unversioned = hashlib.sha256(b"millan|s|1|8|0").hexdigest()[:16]
    assert Scenario("millan", "s", "t").geometry_key != unversioned
