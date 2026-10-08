"""Swap components and see what each change does to the cascade.

A scenario is one choice of **bed x surface x thinning x flotation factor**. This
module runs the whole chain for a scenario and reports every intermediate, so a
swap can be traced rather than just compared at the end.

The cost of a swap is not uniform -- and that shapes the whole design
----------------------------------------------------------------------
Only the *saddle search* is expensive (~100 s on the 29 M-cell grid). What it
depends on decides what a swap costs:

| swap | changes h_pass? | changes V_w? | cost |
|---|---|---|---|
| **bed** | yes | yes, via the bed sill | full re-solve |
| **surface** | yes | no -- the sill is bed-only | full re-solve |
| **flotation factor k** | yes | no | full re-solve, but 1-D in k |
| **thinning dh/dt** | **no** | no | **free** |

Thinning is free because it only sets how fast the threshold rises. The geometry
question "how much thinning is needed" is answered once per bed x surface x k;
dividing by a rate is arithmetic. So the engine caches on
``(bed, surface, k, connectivity)`` and a thinning sweep costs one solve, not one
per rate. This is the same structure as Phase 1, where ``Phi_crit`` depended only
on ``f`` and the grids.

The propagation chain
---------------------
Two branches run from the same swap, and a bed change moves **both**::

    bed, surface, k  ->  h_pass  ->  Delta = h_pass - h_LGP
                                 ->  thinning needed = Delta / kappa
                     ->  (with dh/dt)  ->  years to breach  ->  onset year

    bed              ->  sill      ->  h_final = max(h_GPL, sill)
                                   ->  drawdown  ->  V_w  ->  Q_p

Two sills, and they are not interchangeable
-------------------------------------------
"The sill" is ambiguous, and the two readings differ by 48 m on the Millan bed and
152 m on IceBoost:

``bed_sill_m`` -- the highest bed along the **flotation path**. That route is
chosen to minimise hydraulic potential, not bed elevation, so it will climb over a
bed ridge where the ice is thick. This is the spec's definition and it is the
conservative one: it assumes the flood drains through the conduit that opened and
that the route cannot migrate.

``bed_pass_m`` -- the lowest bed ridge on **any** route between the lakes,
regardless of where the seal opens. This is what you see reading the bed raster in
QGIS. It assumes that once drainage is under way water finds the lowest route.

The truth is somewhere between, so both are reported and ``V_w`` is computed under
each. Which to believe is a physics decision for the owner, not something to bury
in a definition.

A surface or thinning change moves only the timing. A bed change moves the timing
*and* the flood size, which is why bed is the swap to be most careful about.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..config import repo_root
from .lakes import conic_bottom, released_volume, stage_curve, vertical_wall_bottom
from .path import profile, summarize, trace_path
from .seal import head_surface, kappa_from_k, masked_smooth, pass_search

log = logging.getLogger(__name__)

CACHE_DIR = repo_root() / "data" / "interim" / "scenario_cache"


@dataclass(frozen=True)
class Scenario:
    """One combination of swappable components."""

    bed_id: str
    surface_id: str
    thinning_id: str
    k: float = 1.0
    connectivity: int = 8
    smoothing_m: float = 0.0

    @property
    def name(self) -> str:
        return (
            f"{self.bed_id}__{self.surface_id}__{self.thinning_id}"
            f"__k{self.k:g}__c{self.connectivity}"
        )

    @property
    def geometry_key(self) -> str:
        """Identity of the expensive part: everything except the thinning rate.

        The ``v2`` prefix is a schema version, not decoration: v2 added
        ``bed_pass_m`` to the cached payload, and without the bump a v1 cache
        would still be served and silently answer with that field missing.
        Bump it whenever the cached dict gains or changes a field.
        """
        raw = (
            f"v2|{self.bed_id}|{self.surface_id}|{self.k:.10g}|"
            f"{self.connectivity}|{self.smoothing_m:g}"
        )
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


@dataclass
class ScenarioResult:
    """Every intermediate in the chain, so a swap can be traced step by step."""

    scenario: Scenario
    h_pass_m: float
    delta_m: float
    thinning_needed_m: float
    pass_cell: tuple[int, int] | None
    bed_sill_m: float
    bed_pass_m: float
    l_pass_m: float
    l_total_m: float
    n_at_pass_m: float
    h_final_m: float
    drawdown_m: float
    controlled_by: str
    vw_lo_mm3: float
    vw_hi_mm3: float
    dhdt_m_per_yr: float
    years_to_breach: float
    onset_year: float
    # The bed-only reading of the sill, carried alongside the along-path one.
    h_final_bedpass_m: float = float("nan")
    drawdown_bedpass_m: float = float("nan")
    controlled_by_bedpass: str = ""
    vw_lo_bedpass_mm3: float = float("nan")
    vw_hi_bedpass_mm3: float = float("nan")
    qp_tunnel_lo: float = float("nan")
    qp_tunnel_hi: float = float("nan")
    qp_nontunnel_lo: float = float("nan")
    qp_nontunnel_hi: float = float("nan")
    notes: list[str] = field(default_factory=list)

    def as_row(self) -> dict:
        d = {
            "scenario": self.scenario.name,
            "bed": self.scenario.bed_id,
            "surface": self.scenario.surface_id,
            "thinning": self.scenario.thinning_id,
            "k": self.scenario.k,
        }
        d.update({k: v for k, v in self.__dict__.items() if k not in ("scenario", "notes")})
        d["notes"] = "; ".join(self.notes)
        return d


def _cache_path(scenario: Scenario) -> Path:
    return CACHE_DIR / f"geom_{scenario.geometry_key}.json"


def solve_geometry(
    scenario: Scenario,
    zb,
    zs,
    classes,
    transform,
    lake_levels: dict,
    codes: dict,
    use_cache: bool = True,
) -> dict:
    """The expensive half: head surface, pass search, path and bed sill.

    Cached on ``geometry_key``, which deliberately excludes the thinning rate --
    nothing here depends on it.
    """
    cache = _cache_path(scenario)
    if use_cache and cache.exists():
        log.info("cache hit for %s/%s/k=%g", scenario.bed_id, scenario.surface_id, scenario.k)
        return json.loads(cache.read_text())

    surf = zs
    if scenario.smoothing_m:
        surf = masked_smooth(zs, classes == codes["ice"], scenario.smoothing_m, abs(transform.a))
    # Bed-only pass: independent of the head, so it is computed once alongside it.
    bp = bed_pass(zb, classes, codes, scenario.connectivity)
    head, info = head_surface(zb, surf, classes, lake_levels, k=scenario.k)
    res = pass_search(head, classes, codes["gpl"], codes["lgp"], scenario.connectivity)
    out: dict = {"reached": bool(res.reached_target), "head_info": info, "bed_pass_m": bp}
    if not res.reached_target:
        out.update({"h_pass_m": None, "bed_sill_m": None})
    else:
        p = trace_path(res.parent, res.info["best_target_cell"])
        prof = profile(
            p, zb, surf, head, transform, lake_levels[codes["lgp"]], lake_levels[codes["gpl"]]
        )
        summ = summarize(
            prof, classes, lake_levels[codes["lgp"]], res.pass_cell, lgp_class=codes["lgp"]
        )
        out.update(
            {
                "h_pass_m": float(res.h_pass),
                "pass_cell": [int(v) for v in res.pass_cell] if res.pass_cell else None,
                "bed_sill_m": summ.bed_sill_m,
                "l_pass_m": summ.l_pass_m,
                "l_total_m": summ.l_total_m,
                "n_at_pass_m": summ.n_at_pass_m,
                "thickness_at_pass_m": summ.thickness_at_pass_m,
            }
        )
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(out, indent=2, default=str))
    return out


def bed_pass(bed: np.ndarray, classes: np.ndarray, codes: dict, connectivity: int = 8) -> float:
    """Lowest bed ridge on any route between the two lakes (m).

    Both lake footprints are set to ``-inf`` -- water is already there, so they are
    not barriers -- and the minimax then runs over the intervening bed. Seeding is
    region-to-region via ``pass_search``, not cell-to-cell, so the answer is the
    true lowest ridge rather than the ridge between two arbitrary points.
    """
    field = np.where(np.isfinite(bed), bed, np.inf)
    field = np.where((classes == codes["lgp"]) | (classes == codes["gpl"]), -np.inf, field)
    res = pass_search(field, classes, codes["gpl"], codes["lgp"], connectivity)
    return float(res.h_pass) if res.reached_target else float("inf")


def volume_bracket(
    classes,
    codes: dict,
    cell_area_m2: float,
    h_lgp: float,
    h_gpl: float,
    bed_sill: float,
    step: float = 0.5,
) -> dict:
    """V_w bounds for a given sill. Prism above, cone below (D-015)."""
    lgp = classes == codes["lgp"]
    h_final = max(h_gpl, bed_sill)
    drawdown = max(h_lgp - h_final, 0.0)
    if drawdown <= 0:
        return {
            "lo": 0.0,
            "hi": 0.0,
            "h_final": h_final,
            "drawdown": 0.0,
            "controlled_by": "bed_sill" if bed_sill >= h_gpl else "receiving_lake",
        }
    stages = np.arange(h_lgp - drawdown - 1.0, h_lgp + step, step)
    out = {}
    for key, bottom in (
        ("hi", vertical_wall_bottom(lgp, h_lgp - drawdown)),
        ("lo", conic_bottom(lgp, h_lgp, drawdown)),
    ):
        curve = stage_curve(bottom, lgp, stages, cell_area_m2)
        out[key] = released_volume(curve, h_lgp, h_gpl, bed_sill)["Vw_Mm3"]
    out.update(
        {
            "h_final": h_final,
            "drawdown": drawdown,
            "controlled_by": "bed_sill" if bed_sill >= h_gpl else "receiving_lake",
        }
    )
    return out


def run_scenario(
    scenario: Scenario,
    grids: dict,
    lake_levels: dict,
    codes: dict,
    dhdt: float,
    base_year: int = 2018,
    use_cache: bool = True,
) -> ScenarioResult | None:
    """Run the full chain for one scenario."""
    from .inflow_empirical import RELATIONS

    geom = solve_geometry(
        scenario,
        grids["bed"],
        grids["surface"],
        grids["classes"],
        grids["transform"],
        lake_levels,
        codes,
        use_cache=use_cache,
    )
    notes = []
    if not geom.get("reached"):
        log.warning("%s: no route from GPL to LGP", scenario.name)
        return None

    h_lgp, h_gpl = lake_levels[codes["lgp"]], lake_levels[codes["gpl"]]
    h_pass = geom["h_pass_m"]
    kappa = kappa_from_k(scenario.k)
    delta = h_pass - h_lgp
    thinning = delta / kappa
    years = thinning / abs(dhdt) if dhdt else float("inf")
    if delta <= 0:
        notes.append("already open at this geometry")
        years = 0.0

    vol = volume_bracket(
        grids["classes"], codes, grids["cell_area_m2"], h_lgp, h_gpl, geom["bed_sill_m"]
    )
    # The same volume calculation under the bed-only reading of the sill.
    bp = geom.get("bed_pass_m", float("nan"))
    if np.isfinite(bp):
        vol_bp = volume_bracket(grids["classes"], codes, grids["cell_area_m2"], h_lgp, h_gpl, bp)
    else:
        vol_bp = {
            "lo": float("nan"),
            "hi": float("nan"),
            "h_final": float("nan"),
            "drawdown": float("nan"),
            "controlled_by": "no bed route",
        }

    if vol["drawdown"] <= 0:
        notes.append("along-path sill at or above the lake surface: no drainage on that route")
    if np.isfinite(bp) and vol_bp["drawdown"] > vol["drawdown"] + 1.0:
        notes.append(
            f"a lower bed route exists at {bp:.1f} m; drawdown would be "
            f"{vol_bp['drawdown']:.1f} m instead of {vol['drawdown']:.1f} m"
        )

    qp = {}
    for key in ("walder_costa_tunnel", "walder_costa_nontunnel"):
        rel = RELATIONS[key]
        qp[key] = (rel.peak_discharge(vol["lo"]), rel.peak_discharge(vol["hi"]))

    return ScenarioResult(
        scenario=scenario,
        h_pass_m=h_pass,
        delta_m=delta,
        thinning_needed_m=thinning,
        pass_cell=tuple(geom["pass_cell"]) if geom.get("pass_cell") else None,
        bed_sill_m=geom["bed_sill_m"],
        l_pass_m=geom.get("l_pass_m", float("nan")),
        l_total_m=geom.get("l_total_m", float("nan")),
        n_at_pass_m=geom.get("n_at_pass_m", float("nan")),
        h_final_m=vol["h_final"],
        drawdown_m=vol["drawdown"],
        controlled_by=vol["controlled_by"],
        vw_lo_mm3=vol["lo"],
        vw_hi_mm3=vol["hi"],
        bed_pass_m=bp,
        h_final_bedpass_m=vol_bp["h_final"],
        drawdown_bedpass_m=vol_bp["drawdown"],
        controlled_by_bedpass=vol_bp["controlled_by"],
        vw_lo_bedpass_mm3=vol_bp["lo"],
        vw_hi_bedpass_mm3=vol_bp["hi"],
        dhdt_m_per_yr=dhdt,
        years_to_breach=years,
        onset_year=base_year + years,
        qp_tunnel_lo=qp["walder_costa_tunnel"][0],
        qp_tunnel_hi=qp["walder_costa_tunnel"][1],
        qp_nontunnel_lo=qp["walder_costa_nontunnel"][0],
        qp_nontunnel_hi=qp["walder_costa_nontunnel"][1],
        notes=notes,
    )


def expand(beds, surfaces, thinnings, ks=(1.0,), connectivity=(8,)) -> list[Scenario]:
    """Cartesian product of the component choices."""
    return [
        Scenario(b, s, t, k, c)
        for b in beds
        for s in surfaces
        for t in thinnings
        for k in ks
        for c in connectivity
    ]
