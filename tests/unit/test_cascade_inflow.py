"""inflow_empirical.py — the spec's expected values and the volume constraint."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.cascade.inflow_empirical import (
    RELATIONS,
    hydrograph,
    summarise,
    unit_shape,
)

V500 = 500.0


def test_spec_expected_peak_discharges() -> None:
    """Spec test table: ~2,780 / ~4,820 / ~16,900 m3/s at V = 500 Mm3."""
    assert RELATIONS["walder_costa_tunnel"].peak_discharge(V500) == pytest.approx(2780, rel=0.01)
    assert RELATIONS["clague_mathews"].peak_discharge(V500) == pytest.approx(4820, rel=0.01)
    assert RELATIONS["walder_costa_nontunnel"].peak_discharge(V500) == pytest.approx(
        16900, rel=0.01
    )


def test_froehlich_needs_a_head_and_reproduces_the_reference_point() -> None:
    f = RELATIONS["froehlich"]
    with pytest.raises(ValueError, match="hw_m"):
        f.peak_discharge(V500)
    assert f.peak_discharge(V500, hw_m=100.0) == pytest.approx(51000, rel=0.02)


def test_froehlich_is_flagged_reference_only() -> None:
    """It is a moraine-breach regression; nothing may silently treat it as applicable."""
    assert RELATIONS["froehlich"].reference_only is True
    assert all(not r.reference_only for k, r in RELATIONS.items() if k != "froehlich")


@pytest.mark.parametrize("key", list(RELATIONS))
def test_every_curve_integrates_to_the_released_volume(key: str) -> None:
    """Spec pass condition: each curve integrates to V_w."""
    rel = RELATIONS[key]
    qp = rel.peak_discharge(V500, hw_m=100.0)
    t, q = hydrograph(qp, V500 * 1e6, rel.shape)
    assert np.trapezoid(q, t) == pytest.approx(V500 * 1e6, rel=1e-9)


@pytest.mark.parametrize("key", list(RELATIONS))
def test_peak_is_preserved(key: str) -> None:
    rel = RELATIONS[key]
    qp = rel.peak_discharge(V500, hw_m=100.0)
    _, q = hydrograph(qp, V500 * 1e6, rel.shape)
    assert q.max() == pytest.approx(qp, rel=1e-9)


def test_triangle_base_time_matches_the_spec_formula() -> None:
    """Spec: non-tunnel base time t_b = 2 V_w / Q_p."""
    qp = RELATIONS["walder_costa_nontunnel"].peak_discharge(V500)
    t, _ = hydrograph(qp, V500 * 1e6, "triangle")
    assert t[-1] == pytest.approx(2 * V500 * 1e6 / qp, rel=1e-6)


def test_tunnel_shape_rises_slowly_and_falls_abruptly() -> None:
    """The defining asymmetry of tunnel drainage."""
    u, g = unit_shape("lognormal_reversed")
    peak = int(np.argmax(g))
    rise, fall = u[peak], 1.0 - u[peak]
    assert rise > fall, f"peak at u={u[peak]:.2f}; rise should be the longer limb"
    assert rise > 0.7


def test_tunnel_drainage_lasts_longer_than_a_breach() -> None:
    """Same volume, far smaller peak, so the flood is drawn out."""
    df = summarise(V500, hw_m=100.0)
    tunnel = df[df.relation == "walder_costa_tunnel"].iloc[0]
    breach = df[df.relation == "walder_costa_nontunnel"].iloc[0]
    assert tunnel.Qp_m3s < breach.Qp_m3s
    assert tunnel.duration_h > breach.duration_h


def test_shapes_are_normalised_and_nonnegative() -> None:
    for shape in ("triangle", "lognormal_reversed"):
        u, g = unit_shape(shape)
        assert g.min() >= 0.0
        assert g.max() == pytest.approx(1.0)
        assert u[0] == pytest.approx(0.0) and u[-1] == pytest.approx(1.0)


def test_triangle_volume_is_independent_of_peak_position() -> None:
    """Area is half base times height wherever the apex sits."""
    qp = 1000.0
    vols = []
    for pf in (0.2, 0.5, 0.8):
        tt, q = hydrograph(qp, 1e8, "triangle", peak_fraction=pf)
        vols.append(np.trapezoid(q, tt))
    assert np.allclose(vols, 1e8, rtol=1e-9)


def test_zero_volume_gives_no_flood() -> None:
    _, q = hydrograph(0.0, 0.0, "triangle")
    assert q.max() == 0.0
    assert RELATIONS["clague_mathews"].peak_discharge(0.0) == 0.0


def test_unknown_shape_rejected() -> None:
    with pytest.raises(ValueError, match="unknown shape"):
        unit_shape("sawtooth")


def test_scaling_is_roughly_two_thirds_power() -> None:
    """Qp ~ V^2/3 for tunnel drainage (Ng & Bjornsson 2003), as the spec notes."""
    r = RELATIONS["walder_costa_tunnel"]
    ratio = r.peak_discharge(1000.0) / r.peak_discharge(100.0)
    assert ratio == pytest.approx(10**0.66, rel=1e-6)
    assert 0.6 < np.log10(ratio) < 0.7
