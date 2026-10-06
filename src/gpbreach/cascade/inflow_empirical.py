"""Empirical peak discharge and hydrograph shapes from the released volume alone.

These relations are the sanity check on the physical conduit model (M6): they use
nothing but V_w, so when the two disagree it is the conduit model that has to
explain itself.

Two drainage modes, and they differ by more than an order of magnitude
---------------------------------------------------------------------
A **subglacial tunnel** enlarges by melting, so discharge climbs slowly and stops
abruptly when the lake empties. A **non-tunnel** failure -- a breach, a marginal
spillway, mechanical collapse of the dam -- releases water far faster. Which mode
applies is a physical judgement, not something the volume decides, so both are
reported and neither is preferred here.

``Froehlich (2025)`` is a **moraine-breach** regression and does not apply to
ice-dam drainage. The spec includes it as a labelled reference point only, and
``REFERENCE_ONLY`` marks it so plots and tables cannot quietly promote it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Relation:
    """An empirical peak-discharge relation. ``V`` in Mm^3, ``Q_p`` in m^3 s^-1."""

    key: str
    label: str
    coefficient: float
    v_exponent: float
    mode: str  # "tunnel" | "non_tunnel" | "moraine_breach"
    citation: str
    hw_exponent: float = 0.0  # 0 unless the relation also uses a head
    reference_only: bool = False

    def peak_discharge(self, v_mm3: float, hw_m: float | None = None) -> float:
        """Q_p (m^3 s^-1) for a released volume ``v_mm3`` (Mm^3)."""
        if v_mm3 <= 0:
            return 0.0
        q = self.coefficient * v_mm3**self.v_exponent
        if self.hw_exponent:
            if hw_m is None:
                raise ValueError(f"{self.key} needs hw_m (head, m)")
            q *= hw_m**self.hw_exponent
        return float(q)

    @property
    def shape(self) -> str:
        """Hydrograph shape implied by the drainage mode."""
        return "lognormal_reversed" if self.mode == "tunnel" else "triangle"


RELATIONS: dict[str, Relation] = {
    "walder_costa_tunnel": Relation(
        "walder_costa_tunnel",
        "Walder & Costa (1996), tunnel",
        46.0,
        0.66,
        "tunnel",
        "Walder & Costa (1996)",
    ),
    "clague_mathews": Relation(
        "clague_mathews", "Clague & Mathews (1973)", 75.0, 0.67, "tunnel", "Clague & Mathews (1973)"
    ),
    "walder_costa_nontunnel": Relation(
        "walder_costa_nontunnel",
        "Walder & Costa (1996), non-tunnel",
        1100.0,
        0.44,
        "non_tunnel",
        "Walder & Costa (1996)",
    ),
    "froehlich": Relation(
        "froehlich",
        "Froehlich (2025) — moraine breach, reference only",
        142.0,
        0.378,
        "moraine_breach",
        "Froehlich (2025)",
        hw_exponent=0.766,
        reference_only=True,
    ),
}


def unit_shape(
    shape: str, n: int = 2001, sigma: float = 0.5, peak_fraction: float = 0.5
) -> tuple[np.ndarray, np.ndarray]:
    """Dimensionless hydrograph ``g(u)`` on ``u`` in [0, 1], scaled to a peak of 1.

    ``lognormal_reversed`` — a lognormal rises sharply and decays slowly; reversed
    in time it rises slowly and falls abruptly, which is the tunnel signature.

    ``triangle`` — ``peak_fraction`` sets where the apex sits. The spec does not
    specify it; 0.5 is the neutral choice and the released volume does not depend
    on it, since a triangle's area is half base times height wherever the apex is.
    """
    if shape == "triangle":
        u = np.linspace(0.0, 1.0, n)
        p = float(np.clip(peak_fraction, 1e-6, 1 - 1e-6))
        g = np.where(u <= p, u / p, (1.0 - u) / (1.0 - p))
        return u, np.clip(g, 0.0, 1.0)

    if shape == "lognormal_reversed":
        # Capture the tail out to where the pdf has fallen well below the peak.
        x = np.linspace(1e-3, 4.0, n)
        pdf = np.exp(-((np.log(x)) ** 2) / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi))
        g = pdf / pdf.max()
        g = g[::-1]  # slow rise, abrupt fall
        u = (x - x[0]) / (x[-1] - x[0])
        return u, g

    raise ValueError(f"unknown shape {shape!r}")


def hydrograph(
    q_peak: float,
    v_w_m3: float,
    shape: str,
    n: int = 2001,
    sigma: float = 0.5,
    peak_fraction: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Discharge against time, honouring both the peak **and** the volume.

    The dimensionless shape fixes the form; the duration is then solved so the
    integral equals ``v_w_m3`` exactly::

        D = V_w / (Q_p * integral of g over u)

    Returns ``(t_seconds, Q_m3s)``. A final numerical rescale removes the residual
    left by discretising the shape, so the curve integrates to V_w to machine
    precision rather than merely close to it.
    """
    if q_peak <= 0 or v_w_m3 <= 0:
        return np.zeros(2), np.zeros(2)
    u, g = unit_shape(shape, n=n, sigma=sigma, peak_fraction=peak_fraction)
    integral = float(np.trapezoid(g, u))
    duration = v_w_m3 / (q_peak * integral)
    t = u * duration
    q = g * q_peak
    # Numerical rescale: the spec asks for the area to equal V_w exactly.
    actual = float(np.trapezoid(q, t))
    if actual > 0:
        t = t * (v_w_m3 / actual)
    return t, q


def summarise(
    v_w_mm3: float,
    hw_m: float | None = None,
    keys=None,
    sigma: float = 0.5,
    peak_fraction: float = 0.5,
):
    """One row per relation: Q_p, duration, time to peak, and the volume check."""
    import pandas as pd

    rows = []
    for key in keys or RELATIONS:
        rel = RELATIONS[key]
        try:
            qp = rel.peak_discharge(v_w_mm3, hw_m)
        except ValueError as exc:
            log.warning("skipping %s: %s", key, exc)
            continue
        t, q = hydrograph(qp, v_w_mm3 * 1e6, rel.shape, sigma=sigma, peak_fraction=peak_fraction)
        vol = float(np.trapezoid(q, t))
        rows.append(
            {
                "relation": key,
                "label": rel.label,
                "mode": rel.mode,
                "reference_only": rel.reference_only,
                "Qp_m3s": qp,
                "shape": rel.shape,
                "duration_h": float(t[-1]) / 3600.0,
                "time_to_peak_h": float(t[int(np.argmax(q))]) / 3600.0,
                "volume_check_Mm3": vol / 1e6,
                "volume_error_frac": abs(vol - v_w_mm3 * 1e6) / (v_w_mm3 * 1e6),
            }
        )
    return pd.DataFrame(rows)
