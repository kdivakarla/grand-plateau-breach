"""Physical inflow hydrograph: a Nye/Clarke-type conduit coupled to both lakes.

Why this module tightens the answer
-----------------------------------
The empirical relations (M5) take only ``V_w`` and return a peak, and which
relation applies depends on the drainage mode -- a choice worth a **factor of 17**
in ``Q_p`` against the 1.9x that volume uncertainty contributes. This model does
not take the mode as input. It integrates the conduit and the mode falls out:

* **opening** scales as ``Q ~ S^(4/3)`` -- superlinear in conduit area;
* **closure** scales as ``S * N^3`` -- linear in area.

Above a threshold size opening outruns closure and the conduit runs away, which
*is* tunnel drainage. Below it the conduit pinches shut and the lake cannot drain
through a tunnel at all -- which is itself evidence that the real system drains
some other way (sheet flow, marginal, or mechanical failure of the dam).

So a run that **fails** to empty the lake is as informative as one that succeeds.
``stop_reason`` reports which happened.

What it cannot settle
---------------------
Non-tunnel drainage includes mechanical failure of the ice dam, which is not
conduit physics. A successful tunnel run shows a tunnel *can* carry the flood; it
does not prove the dam will not fail some other way first.

Equations (spec section conduit.py)::

    A_L(h_L) dh_L/dt = Q_in,L - Q
    A_G(h_G) dh_G/dt = Q + Q_in,G - O(h_G)
    Q        = S^(4/3) / (n' (4 pi)^(1/3)) * ((h_L - h_G)/L)^(1/2)
    dS/dt    = Q rho_w g (h_L-h_G)/L / (rho_i L_f)      <- viscous dissipation
             + rho_w c_w Q T_L / (rho_i L_f L)          <- lake sensible heat
             - 2 A S (N/n)^n                            <- creep closure
    N        = rho_w g max(h_pass - (h_L + h_G)/2, 0)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

log = logging.getLogger(__name__)

RHO_W, RHO_I, G = 1000.0, 917.0, 9.81
L_F, C_W = 3.34e5, 4186.0
FOUR_PI_CBRT = (4.0 * np.pi) ** (1.0 / 3.0)


@dataclass
class ConduitParams:
    """Conduit and flow-law parameters. See docs/decisions.md for provenance."""

    manning_n: float = 0.08
    S0_m2: float = 0.1
    lake_temp_C: float = 0.5
    A_glen: float = 2.4e-24
    n_glen: int = 3
    L_m: float = 5624.0  # path length, L_pass by default (D-014)
    h_pass_m: float = 233.262
    q_in_lgp: float = 0.0  # background inflow to LGP, m3/s
    q_in_gpl: float = 0.0


@dataclass
class ConduitResult:
    t_s: np.ndarray
    Q_m3s: np.ndarray
    h_L_m: np.ndarray
    h_G_m: np.ndarray
    S_m2: np.ndarray
    stop_reason: str
    q_peak_m3s: float
    t_peak_s: float
    volume_released_m3: float
    duration_s: float
    mass_balance_error: float = float("nan")
    meta: dict = field(default_factory=dict)

    @property
    def drained_fully(self) -> bool:
        return self.stop_reason == "reached_floor"


def discharge(S: float, h_L: float, h_G: float, p: ConduitParams) -> float:
    """Manning flow in a full circular conduit (m^3/s). Zero if the head reverses."""
    dh = h_L - h_G
    if dh <= 0 or S <= 0:
        return 0.0
    return S ** (4.0 / 3.0) / (p.manning_n * FOUR_PI_CBRT) * np.sqrt(dh / p.L_m)


def effective_pressure(h_L: float, h_G: float, p: ConduitParams) -> float:
    """N in Pa, referenced to the mean of the two lake levels (spec conduit.py)."""
    return RHO_W * G * max(p.h_pass_m - 0.5 * (h_L + h_G), 0.0)


def _area_fn(curve):
    """Interpolator for lake area against stage (m^2), flat outside the range."""
    h = np.asarray(curve["h"], dtype=float)
    a = np.asarray(curve["A_m2"], dtype=float)
    order = np.argsort(h)
    h, a = h[order], a[order]

    def area(stage: float) -> float:
        return float(np.interp(stage, h, a, left=a[0], right=a[-1]))

    return area


def simulate(
    lgp_curve,
    gpl_curve,
    h_lgp0: float,
    h_gpl0: float,
    h_floor: float,
    p: ConduitParams,
    gpl_outflow=None,
    t_max_s: float = 60 * 86400,
    rtol: float = 1e-6,
    atol: float = 1e-8,
    max_step_s: float | None = None,
) -> ConduitResult:
    """Integrate the coupled system until the lake drains or the conduit closes.

    ``gpl_outflow`` is a callable ``O(h_G) -> m^3/s``; ``None`` means a closed
    receiving lake, which is conservative because GPL then rises further and
    reduces the head driving the flood.
    """
    from scipy.integrate import solve_ivp

    A_L, A_G = _area_fn(lgp_curve), _area_fn(gpl_curve)
    out = gpl_outflow or (lambda h: 0.0)

    def rhs(_t, y):
        h_L, h_G, S = y
        S = max(S, 1e-9)
        Q = discharge(S, h_L, h_G, p)
        N = effective_pressure(h_L, h_G, p)
        dh = max(h_L - h_G, 0.0)
        melt_visc = Q * RHO_W * G * (dh / p.L_m) / (RHO_I * L_F)
        melt_heat = RHO_W * C_W * Q * p.lake_temp_C / (RHO_I * L_F * p.L_m)
        close = 2.0 * p.A_glen * S * (N / p.n_glen) ** p.n_glen
        return [
            (p.q_in_lgp - Q) / max(A_L(h_L), 1.0),
            (Q + p.q_in_gpl - out(h_G)) / max(A_G(h_G), 1.0),
            melt_visc + melt_heat - close,
        ]

    def ev_floor(_t, y):
        return y[0] - h_floor

    ev_floor.terminal, ev_floor.direction = True, -1

    def ev_closed(_t, y):
        return y[2] - 1e-4

    ev_closed.terminal, ev_closed.direction = True, -1

    sol = solve_ivp(
        rhs,
        (0.0, t_max_s),
        [h_lgp0, h_gpl0, p.S0_m2],
        method="LSODA",
        events=[ev_floor, ev_closed],
        rtol=rtol,
        atol=atol,
        dense_output=True,
        max_step=max_step_s or np.inf,
    )

    t, (h_L, h_G, S) = sol.t, sol.y
    Q = np.array([discharge(s, a, b, p) for s, a, b in zip(S, h_L, h_G, strict=True)])
    if sol.t_events[0].size:
        reason = "reached_floor"
    elif sol.t_events[1].size:
        reason = "conduit_closed"
    else:
        reason = "t_max"

    released = float(np.trapezoid(Q, t)) if t.size > 1 else 0.0
    # Independent check: water leaving LGP by stage change must match the integral.
    v_lost = 0.0
    if t.size > 1:
        hh = np.linspace(h_L[-1], h_L[0], 200)
        v_lost = float(np.trapezoid([A_L(x) for x in hh], hh))
    err = abs(released - v_lost) / v_lost if v_lost > 0 else float("nan")
    ipk = int(np.argmax(Q)) if Q.size else 0

    return ConduitResult(
        t_s=t,
        Q_m3s=Q,
        h_L_m=h_L,
        h_G_m=h_G,
        S_m2=S,
        stop_reason=reason,
        q_peak_m3s=float(Q.max()) if Q.size else 0.0,
        t_peak_s=float(t[ipk]) if t.size else 0.0,
        volume_released_m3=released,
        duration_s=float(t[-1]) if t.size else 0.0,
        mass_balance_error=err,
        meta={
            "S_final_m2": float(S[-1]),
            "S_max_m2": float(S.max()),
            "h_L_final_m": float(h_L[-1]),
            "h_G_final_m": float(h_G[-1]),
            "n_steps": int(t.size),
            "solver_message": sol.message,
        },
    )


def spill_weir(spill_stage_m: float, width_m: float, C: float = 1.7):
    """GPL outflow once it tops its divide: a broad-crested weir over the front.

    Placeholder until the routing stage exists (spec). Below ``spill_stage_m`` the
    outflow is zero, which is conservative: GPL rises, the head falls, and the
    flood is attenuated rather than exaggerated.
    """

    def outflow(h_G: float) -> float:
        head = max(h_G - spill_stage_m, 0.0)
        return C * width_m * head**1.5

    return outflow
