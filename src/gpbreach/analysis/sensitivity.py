"""Local sensitivity of the breach date to each parameter.

This is the bridge to the Phase 3 uncertainty work (CLAUDE.md Section 9): before
choosing distributions, it is worth knowing how far the answer moves per unit of
each input. Nothing here assumes a distribution -- it reports derivatives, which
are properties of the model, not of anyone's belief about the inputs.

From components.criterion equation (3), with ``m = -dhdt > 0``:

    year = t0 + (Phi_crit(f) - L) / (f * m)

so, holding the saddle location fixed,

    d(year)/dL  = -1 / (f*m)
    d(year)/dm  = -(year - t0) / m
    d(year)/dS  = +1 / m           (uniform surface bias; dPhi/dS = f)
    d(year)/dB  = +(1-f) / (f*m)   (uniform bed bias;     dPhi/dB = 1-f)

``d(year)/df`` is evaluated numerically because the saddle *location* can switch
routes as ``f`` changes, which the closed form does not capture.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..components.criterion import FlotationConnectivity
from ..dam import DamObject


def breach_year(dam: DamObject, criterion: FlotationConnectivity, f: float, level: float,
                dhdt: float, bed_offset: float = 0.0, surface_offset: float = 0.0) -> float:
    """Continuous breach year, optionally with uniform offsets applied to the grids."""
    phi = criterion.critical_threshold(
        dam.surface + surface_offset, dam.bed + bed_offset, f, dam.seed_lake, dam.seed_target
    )
    t = criterion.breach_time(np.array([phi]), np.array([level]), np.array([dhdt]),
                              np.array([f]))[0]
    return float(dam.base_year + t)


def sensitivity_table(dam: DamObject, f: float, level: float, dhdt: float,
                      connectivity: int = 8, steps: dict[str, float] | None = None) -> pd.DataFrame:
    """Central-difference sensitivity of the breach year to each input.

    ``steps`` gives the half-step used for each parameter; the defaults are
    deliberately generous so that a route switch in the saddle shows up rather
    than hiding inside a tiny perturbation.
    """
    crit = FlotationConnectivity(connectivity=connectivity)
    steps = steps or {
        "flotation_fraction": 0.02,
        "lake_level_m": 10.0,
        "dhdt_m_per_yr": 1.0,
        "bed_offset_m": 20.0,
        "surface_offset_m": 10.0,
    }
    base = breach_year(dam, crit, f, level, dhdt)

    def probe(**kw) -> float:
        args = {"f": f, "level": level, "dhdt": dhdt, "bed_offset": 0.0, "surface_offset": 0.0}
        args.update(kw)
        return breach_year(dam, crit, **args)

    rows = []
    for name, h in steps.items():
        if name == "flotation_fraction":
            hi, lo = probe(f=f + h), probe(f=f - h)
        elif name == "lake_level_m":
            hi, lo = probe(level=level + h), probe(level=level - h)
        elif name == "dhdt_m_per_yr":
            hi, lo = probe(dhdt=dhdt + h), probe(dhdt=dhdt - h)
        elif name == "bed_offset_m":
            hi, lo = probe(bed_offset=h), probe(bed_offset=-h)
        else:
            hi, lo = probe(surface_offset=h), probe(surface_offset=-h)
        rows.append(
            {
                "parameter": name,
                "half_step": h,
                "year_minus": lo,
                "year_base": base,
                "year_plus": hi,
                "d_year_per_unit": (hi - lo) / (2 * h),
                "year_span_over_step": hi - lo,
            }
        )
    table = pd.DataFrame(rows)
    table["abs_span"] = table["year_span_over_step"].abs()
    return table.sort_values("abs_span", ascending=False).drop(columns="abs_span").reset_index(
        drop=True
    )
