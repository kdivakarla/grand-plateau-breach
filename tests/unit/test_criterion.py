"""The flotation algebra: ratio form, flotation form, and the closed-form solve."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.constants import (
    LEGACY_FLOTATION_FRACTION,
    LEGACY_THICKNESS_RATIO,
    fraction_to_ratio,
    ratio_to_fraction,
)
from gpbreach.components.criterion import FlotationConnectivity, flotation_field


def test_ratio_and_fraction_are_inverses() -> None:
    for r in (0.05, 0.1, 1 / 9, 0.2):
        assert fraction_to_ratio(ratio_to_fraction(r)) == pytest.approx(r)


def test_legacy_ratio_is_not_exactly_f_090() -> None:
    """r = 0.1 implies f = 0.9091, not 0.90. Recorded as D-003, not silently fixed."""
    assert LEGACY_FLOTATION_FRACTION == pytest.approx(1 / 1.1)
    assert LEGACY_FLOTATION_FRACTION != pytest.approx(0.90, abs=1e-4)
    assert ratio_to_fraction(1 / 9) == pytest.approx(0.90)
    assert LEGACY_THICKNESS_RATIO == 0.1


def test_ratio_and_flotation_criteria_agree() -> None:
    """H_above <= r*H_below  <=>  f*S_t + (1-f)*B <= L, for the same parameters."""
    rng = np.random.default_rng(0)
    S = rng.uniform(0, 600, 5000)
    B = rng.uniform(-400, 200, 5000)
    L, dhdt, t = 110.0, -9.05, 19.0
    r = LEGACY_THICKNESS_RATIO
    f = ratio_to_fraction(r)

    S_t = S + dhdt * t
    ratio_form = (S_t - L) <= r * (L - B)
    flotation_form = flotation_field(S_t, B, f) <= L
    assert np.array_equal(ratio_form, flotation_form)

    # and the flotation statement h_lake >= f*H
    h_lake = L - B
    H = S_t - B
    assert np.array_equal(ratio_form, h_lake >= f * H)


def test_static_field_shifts_linearly_with_time() -> None:
    """Phi_f(S + dhdt*t) == Phi_f(S) + f*dhdt*t -- the basis of the closed form."""
    rng = np.random.default_rng(1)
    S, B = rng.uniform(0, 600, 100), rng.uniform(-400, 200, 100)
    f, dhdt, t = 0.9090909090909091, -9.05, 7.0
    assert flotation_field(S + dhdt * t, B, f) == pytest.approx(
        flotation_field(S, B, f) + f * dhdt * t
    )


def test_breach_time_closed_form() -> None:
    crit = FlotationConnectivity()
    phi, L, dhdt, f = 271.5517938787287, 110.0, -9.05, 1 / 1.1
    t = crit.breach_time(np.array([phi]), np.array([L]), np.array([dhdt]), np.array([f]))
    assert 2010 + t[0] == pytest.approx(2029.63613, abs=1e-4)
    # at that instant the threshold has exactly reached the saddle
    assert crit.threshold_at(L, dhdt, f, t[0]) == pytest.approx(phi)


def test_breach_time_is_vectorised_over_realizations() -> None:
    crit = FlotationConnectivity()
    n = 1000
    rng = np.random.default_rng(2)
    phi = np.full(n, 271.55)
    L = rng.normal(110.0, 5.0, n)
    dhdt = rng.normal(-9.05, 1.0, n)
    f = rng.uniform(0.85, 0.95, n)
    t = crit.breach_time(phi, L, dhdt, f)
    assert t.shape == (n,)
    # spot-check one element against the scalar formula
    i = 17
    assert t[i] == pytest.approx((phi[i] - L[i]) / (-f[i] * dhdt[i]))


def test_no_thinning_never_breaches() -> None:
    crit = FlotationConnectivity()
    t = crit.breach_time(np.array([271.55]), np.array([110.0]), np.array([0.0]),
                         np.array([0.9]))
    assert np.isinf(t[0])


def test_higher_lake_level_breaches_sooner() -> None:
    crit = FlotationConnectivity()
    args = (np.array([271.55, 271.55]), np.array([100.0, 120.0]),
            np.array([-9.05, -9.05]), np.array([1 / 1.1, 1 / 1.1]))
    t = crit.breach_time(*args)
    assert t[1] < t[0]


def test_invalid_flotation_fraction_rejected() -> None:
    with pytest.raises(ValueError):
        flotation_field(np.zeros(3), np.zeros(3), 0.0)
    with pytest.raises(ValueError):
        flotation_field(np.zeros(3), np.zeros(3), 1.5)
