"""Invert an observed breach for the flotation factor k.

A forecast runs the model forward and asks *when* the seal opens. A hindcast
already knows the answer, so it runs the model **backwards** and asks *what value
of k would have opened the seal when it actually did.*

The inversion
-------------
At the moment of breach the lake level equals the pass head::

    h_pass(k) = h_lake(observed)

``h_pass`` rises monotonically with ``k`` -- more ice weight means a higher
barrier -- so this is a one-dimensional root find. Each evaluation is a full pass
search (~35 s on a 30 M-cell grid), so the curve is tabulated once and the root
taken by interpolation, then refined by a few exact evaluations near the answer.

Reading the result
------------------
``k = 1`` is pure flotation: the lake floated its dam exactly. Hindcasts of past
river-piracy events are expected to give ``k < 1`` -- lakes tend to drain *before*
reaching full flotation, and the shortfall tends to be larger for thicker dams.
A k above 1 means the lake exceeded flotation without draining, which usually
indicates a problem in the inputs rather than in the ice.

What this does **not** calibrate: the drainage mode, the sill definition, or the
thinning rate. Those need the observed hydrograph, volume, and surface history
respectively.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from .seal import head_surface, kappa_from_k, pass_search

log = logging.getLogger(__name__)


@dataclass
class Calibration:
    """Result of inverting one observed breach."""

    implied_k: float
    implied_kappa: float
    h_lake_m: float
    h_pass_at_k1_m: float
    bracket: tuple[float, float]
    table: list[tuple[float, float]] = field(default_factory=list)
    status: str = "ok"
    note: str = ""

    def as_row(self) -> dict:
        return {
            "implied_k": self.implied_k,
            "implied_kappa": self.implied_kappa,
            "h_lake_m": self.h_lake_m,
            "h_pass_at_k1_m": self.h_pass_at_k1_m,
            "k_lo": self.bracket[0],
            "k_hi": self.bracket[1],
            "status": self.status,
            "note": self.note,
        }


def h_pass_at(
    k: float, zb, zs, classes, codes: dict, lake_levels: dict, connectivity: int = 8
) -> float:
    """Pass head for one flotation factor (m). One full search."""
    head, _ = head_surface(zb, zs, classes, lake_levels, k=k)
    res = pass_search(head, classes, codes["gpl"], codes["lgp"], connectivity)
    return float(res.h_pass) if res.reached_target else float("inf")


def implied_k(
    zb,
    zs,
    classes,
    codes: dict,
    lake_levels: dict,
    h_lake: float,
    k_min: float = 0.5,
    k_max: float = 1.2,
    n_table: int = 9,
    refine: int = 6,
    connectivity: int = 8,
) -> Calibration:
    """Solve ``h_pass(k) = h_lake`` for the observed breach.

    ``lake_levels`` sets the lake cells in the head surface; the *target* level is
    ``h_lake``, which for a hindcast is the level at the moment of drainage.
    """
    ks = np.linspace(k_min, k_max, n_table)
    table: list[tuple[float, float]] = []
    for k in ks:
        hp = h_pass_at(float(k), zb, zs, classes, codes, lake_levels, connectivity)
        table.append((float(k), hp))
        log.info("k=%.3f -> h_pass=%.3f", k, hp)

    finite = [(k, h) for k, h in table if np.isfinite(h)]
    if len(finite) < 2:
        return Calibration(
            float("nan"),
            float("nan"),
            h_lake,
            float("nan"),
            (k_min, k_max),
            table,
            "no_route",
            "the pass search found no route at any k",
        )

    kk = np.array([k for k, _ in finite])
    hh = np.array([h for _, h in finite])
    h_at_1 = float(np.interp(1.0, kk, hh))

    if hh.min() > h_lake:
        return Calibration(
            float("nan"),
            float("nan"),
            h_lake,
            h_at_1,
            (k_min, k_max),
            table,
            "below_range",
            f"h_pass exceeds the lake level even at k={kk[0]:.2f}; "
            "the implied k is below the search range",
        )
    if hh.max() < h_lake:
        return Calibration(
            float("nan"),
            float("nan"),
            h_lake,
            h_at_1,
            (k_min, k_max),
            table,
            "above_range",
            f"h_pass is below the lake level even at k={kk[-1]:.2f}; "
            "the seal should already have been open",
        )

    # bracket, then refine with exact evaluations
    i = int(np.searchsorted(hh, h_lake))
    lo, hi = kk[max(i - 1, 0)], kk[min(i, kk.size - 1)]
    for _ in range(refine):
        mid = 0.5 * (lo + hi)
        hp = h_pass_at(float(mid), zb, zs, classes, codes, lake_levels, connectivity)
        table.append((float(mid), hp))
        if not np.isfinite(hp):
            break
        if hp < h_lake:
            lo = mid
        else:
            hi = mid
    k_hat = 0.5 * (lo + hi)
    return Calibration(
        float(k_hat),
        float(kappa_from_k(k_hat)),
        h_lake,
        h_at_1,
        (float(lo), float(hi)),
        table,
        "ok",
        ""
        if k_hat <= 1.0
        else "k > 1 means the lake passed flotation without draining; "
        "suspect the inputs rather than the ice",
    )


def summarise(results: dict[str, Calibration]):
    """Pool calibrations across events into a prior for k."""
    import pandas as pd

    rows = []
    for name, c in results.items():
        r = c.as_row()
        r["event"] = name
        rows.append(r)
    df = pd.DataFrame(rows)
    good = df[df.status == "ok"]
    if not good.empty:
        df.attrs["k_mean"] = float(good.implied_k.mean())
        df.attrs["k_sd"] = float(good.implied_k.std())
        df.attrs["k_min"] = float(good.implied_k.min())
        df.attrs["k_max"] = float(good.implied_k.max())
        df.attrs["n_events"] = len(good)
    return df
