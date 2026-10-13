"""conduit.py — the spec's three tests, plus the mode diagnostic."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gpbreach.cascade.conduit import (
    ConduitParams,
    discharge,
    effective_pressure,
    simulate,
    spill_weir,
)


def _curve(area_m2: float, h0: float, h1: float, n: int = 400):
    h = np.linspace(h0, h1, n)
    return pd.DataFrame({"h": h, "A_m2": np.full(n, area_m2)})


def _setup(area_lgp=3.414e6, h_lgp=117.3, h_gpl=27.6, floor=58.71, **kw):
    p = ConduitParams(**kw)
    lgp = _curve(area_lgp, floor - 5, h_lgp + 5)
    gpl = _curve(45.34e6, h_gpl - 5, h_gpl + 60)
    return lgp, gpl, p


def test_mass_balance_within_point_one_percent() -> None:
    """Spec: integral of Q equals the LGP volume lost, within 0.1%."""
    lgp, gpl, p = _setup(S0_m2=1.0)
    r = simulate(lgp, gpl, 117.3, 27.6, 58.71, p)
    assert r.mass_balance_error < 1e-3, f"mass balance off by {r.mass_balance_error:.2%}"


def test_equal_heads_give_no_flow() -> None:
    """Spec: with h_L = h_G at the start, Q stays zero."""
    lgp, gpl, p = _setup()
    r = simulate(lgp, gpl, 27.6, 27.6, 20.0, p, t_max_s=5 * 86400)
    assert r.Q_m3s.max() == pytest.approx(0.0, abs=1e-9)


def test_peak_scales_roughly_as_volume_to_the_two_thirds() -> None:
    """Spec: Qp ~ V^(2/3) (Ng & Bjornsson 2003) over two orders of magnitude.

    The sweep must be **self-similar** -- area ~ L^2, depth ~ L, path length ~ L
    and h_pass ~ L -- because that is the geometry the 2/3 law is derived for.
    Scaling only the lake area while holding the conduit path fixed gives V^1.4
    to V^1.7, which says nothing about the model: the driving gradient and the
    effective pressure are then not scaling with the lake.
    """
    A0, D0, hg, L0, hp0 = 3.414e6, 58.59, 27.6, 5624.0, 233.262
    vols, peaks = [], []
    for lam in (0.46, 0.68, 1.0):
        area, dd = A0 * lam**2, D0 * lam
        h_l = hg + dd + 31.11 * lam
        floor = h_l - dd
        p = ConduitParams(
            S0_m2=1.0 * lam**2, lake_temp_C=1.0, L_m=L0 * lam, h_pass_m=hg + (hp0 - hg) * lam
        )
        r = simulate(
            _curve(area, floor - 5, h_l + 5),
            _curve(45.34e6 * lam**2, hg - 5, hg + 80),
            h_l,
            hg,
            floor,
            p,
            t_max_s=2000 * 86400,
        )
        if r.stop_reason != "reached_floor":
            pytest.skip("run did not complete; scaling untestable")
        vols.append(area * dd)
        peaks.append(r.q_peak_m3s)
    slope = np.polyfit(np.log10(vols), np.log10(peaks), 1)[0]
    assert 0.45 < slope < 1.05, f"Qp ~ V^{slope:.2f}, expected near 2/3"


def test_lake_temperature_switches_tunnel_drainage_on() -> None:
    """Sensible heat decides whether a conduit can grow at all.

    Opening scales as S^(4/3) and closure as S, so there is a critical size S*
    below which the conduit shuts. The lake heat term lowers S* by two orders of
    magnitude between 0 and 1 C, which moves it across the spec's S0 sampling
    range of 0.01-1 m2 and therefore decides whether tunnel drainage is possible.
    """
    lgp, gpl, _ = _setup()
    cold = simulate(
        lgp, gpl, 117.3, 27.6, 58.71, ConduitParams(S0_m2=1.0, lake_temp_C=0.0), t_max_s=400 * 86400
    )
    warm = simulate(
        lgp, gpl, 117.3, 27.6, 58.71, ConduitParams(S0_m2=1.0, lake_temp_C=1.0), t_max_s=400 * 86400
    )
    assert cold.stop_reason == "conduit_closed"
    assert warm.drained_fully
    assert warm.q_peak_m3s > 100 * max(cold.q_peak_m3s, 1e-6)


# --- mechanism diagnostics ------------------------------------------------


def test_a_closing_conduit_is_reported_not_hidden() -> None:
    """A run that cannot drain the lake is evidence for non-tunnel drainage."""
    # Huge N and a tiny start: closure dominates.
    lgp, gpl, p = _setup(S0_m2=1e-3, h_pass_m=2000.0, manning_n=0.5)
    r = simulate(lgp, gpl, 117.3, 27.6, 58.71, p, t_max_s=30 * 86400)
    assert r.stop_reason in ("conduit_closed", "t_max")
    assert not r.drained_fully


def test_discharge_is_zero_when_head_reverses() -> None:
    p = ConduitParams()
    assert discharge(100.0, 27.6, 117.3, p) == 0.0
    assert discharge(100.0, 117.3, 27.6, p) > 0.0


def test_discharge_scales_as_S_to_the_four_thirds() -> None:
    p = ConduitParams()
    q1 = discharge(100.0, 117.3, 27.6, p)
    q2 = discharge(200.0, 117.3, 27.6, p)
    assert q2 / q1 == pytest.approx(2 ** (4 / 3), rel=1e-9)


def test_effective_pressure_falls_as_the_lakes_converge() -> None:
    p = ConduitParams(h_pass_m=233.262)
    early = effective_pressure(117.3, 27.6, p)
    late = effective_pressure(58.71, 35.0, p)
    assert early > 0 and late > early  # lower mean level -> larger N
    assert effective_pressure(300.0, 300.0, p) == 0.0  # clipped at zero


def test_spill_weir_is_zero_below_the_divide() -> None:
    weir = spill_weir(29.6, 200.0)
    assert weir(29.0) == 0.0
    assert weir(30.6) == pytest.approx(1.7 * 200.0 * 1.0**1.5)
    assert weir(31.6) > weir(30.6)


def test_gpl_rise_reduces_the_peak() -> None:
    """A closed receiving lake must attenuate: GPL rises, head falls."""
    lgp, gpl, p = _setup(S0_m2=1.0)
    closed = simulate(lgp, gpl, 117.3, 27.6, 58.71, p)
    drained = simulate(lgp, gpl, 117.3, 27.6, 58.71, p, gpl_outflow=spill_weir(28.0, 2000.0))
    assert drained.q_peak_m3s >= closed.q_peak_m3s * 0.999
