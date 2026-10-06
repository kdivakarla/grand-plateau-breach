# Phase 2 plan — LGP→GPL pass search and inflow hydrograph

Claude Code's implementation plan in response to
*Grand Plateau LGP→GPL Breach: Pass Search & Inflow Hydrograph Build Spec*
(Krishna, 2026-10-03). **Draft for owner review — no code written yet.**

The spec is detailed and I am not proposing changes to its physics. What follows
is (1) what I verified on disk, (2) three findings that change the build order,
(3) an architecture recommendation, (4) milestones, (5) what blocks me.

---

## 1. Input audit

| Spec input | On disk | Status |
|---|---|---|
| `bed_ellip.tif` | — | **MISSING** |
| `bed_filled.tif` | `data/processed/bed_filled.tif` | present, 2496×2543 @ 5 m |
| `surf_50m.tif` | — | **MISSING** |
| `surf_filled.tif` | — | **MISSING** |
| `classes.tif` | `data/processed/classes.tif` | present, 2496×2543 @ 5 m, classes 0–4 as specified |
| `ice_2018.gpkg` | `data/processed/ice_2018.gpkg` | present |
| `lakes_2018.gpkg` | `data/processed/lakes_2018.gpkg` | present, `class` field correct (2/3/4) |
| LGP bathymetry | — | not located; format still open |
| lake levels | — | config values not yet supplied |
| *(extra)* | `data/processed/head_2018.tif` | present, **5240×5337 @ 2.38 m** |

Also new since the last session: `arcticDEM/` (2 m DEM, hillshade, Alsek outline,
2018-09-02/05) and `Geoid_height_conversion/` (`us_noaa_g2012ba0.tif` = GEOID12B,
`reprojected_geoid_undulation.tif`, `bedrock_height_corrected.tif`). These look
like the provenance behind the prepared rasters, which closes most of D-008.

**The ice surface is the blocking gap.** Everything in `seal.py` needs `z_s`, and
no surface raster is on disk at any grid. The spec says to read originals rather
than `_filled`, but for the bed only `_filled` exists.

---

## 2. Three findings that change the build order

### 2.1 `seal.py` is largely already built — and it is the same mathematics as Phase 1

The spec's head surface is **algebraically identical** to the flotation field the
baseline already uses. Writing `κ = k·ρ_i/ρ_w`:

```
spec      h_φ  = z_b + k·(ρ_i/ρ_w)·(z_s − z_b)
gpbreach  Φ_f  = f·z_s + (1 − f)·z_b
```

Expand the first and they are the same expression with **`f = κ = k·ρ_i/ρ_w`**.
Verified numerically over 10⁵ random points: max difference 2×10⁻¹³ m at
k = 1.0, 0.95, 0.90.

And the spec's `pass_search` — "lowest possible highest head on any route" — is a
minimax/bottleneck search, which is exactly what `src/gpbreach/saddle.py`
already does and what `tests/unit/test_saddle.py` already covers:

| spec function | existing equivalent |
|---|---|
| `pass_search` → `h_pass` | `saddle.saddle_threshold()` |
| `trace_path` | `saddle.minimax_path()` |
| pass cell | `saddle.locate_saddle()` |
| `head_surface` | `components.criterion.flotation_field()` |

**Consequence worth noting:** the baseline's `f = 0.9091` corresponds to
**k = 0.9914** — the 2030/2035 dates are already very nearly a pure-flotation
result. A hindcast-calibrated `k < 1` would push both dates *later*, not earlier.
The spec's `k` and decision D-003's `f` are one parameter, so calibrating `k`
feeds the existing baseline directly rather than being a separate exercise.

What is genuinely new in `seal.py`: class-based barriers (land = +inf, lake cells
take their lake level), and a priority-flood implementation returning the **spill
and parent arrays** that `path.py` needs. My current `saddle_threshold` bisects
over sorted levels, which returns the threshold but not the parent links — so the
priority flood should replace it internally while keeping the tested behaviour.

### 2.2 All three lakes are badly truncated by the current grid

The grid is 12.5 × 12.7 km. Clipping the digitised 2018 outlines against it:

| lake | digitised | inside grid | captured |
|---|---|---|---|
| Alsek | 75.61 km² | 1.11 km² | **1.5 %** |
| LGP | 3.41 km² | 0.87 km² | **25.4 %** |
| GPL | 45.34 km² | 23.59 km² | **52.0 %** |

All three extend past the grid edge. **`lakes.py` cannot produce correct
stage–area–volume curves on this grid**, and since the empirical relations go as
roughly `V^0.66`, a factor-2 error in `V_w` is a ~1.6× error in `Q_p` — which
propagates into every downstream result.

So **expanding the grid is a prerequisite for milestone 4, not a tidy-up.** It
likely also matters for `seal.py`: the LGP→GPL route may leave and re-enter the
current footprint. This is the same truncation already flagged under D-005 and
evidenced by the RGI outlines in D-011, now with a second, sharper consequence.

### 2.3 `head_2018.tif` is on a different grid from `bed_filled` and `classes`

`head_2018.tif` is 5240×5337 at **2.38 m**; `bed_filled.tif` and `classes.tif`
are 2496×2543 at **5 m**. The spec requires everything on one grid. Before any
comparison of my `h_pass` against your QGIS margin map is meaningful, we need to
agree which grid is canonical — the 5 m bed grid, or a finer ArcticDEM-native grid.

---

## 3. Architecture: build inside `gpbreach`, not a parallel repo

The spec proposes a new tree (`grand_plateau/`, `src/gp/`). I recommend building
these as modules inside the existing package instead —
`src/gpbreach/cascade/{seal,path,lakes,inflow_empirical,conduit}.py` — because
the existing repo already supplies, tested and working:

- `io/raster.py` — three-backend reader, NaN nodata, grid-identity checks
- `config.py` — YAML loading with the `nominal`/distribution structure the Monte
  Carlo stage will need
- `constants.py`, `units.py` — ρ_i, ρ_w, and the datum-assertion helper
- `saddle.py` — the pass search itself (§2.1)
- `io/atl06.py`, `io/atl13.py` — lake-level observations, quality-screened
- the conda environment, ruff/pytest config, run manifests with git commit and
  input checksums (the spec's "JSON sidecar" requirement, already implemented in
  `analysis/report.py`)
- 45 passing tests including the 2030/2035 regression

A parallel repo duplicates all of that and splits the project in two, with the
`k`↔`f` link spanning the gap. The spec's module names and `run(cfg)` convention
can be kept exactly; only the parent directory changes. The spec's
`data/standard/` maps onto the existing `data/processed/`.

**If you prefer the separate repo, say so and I will follow the spec as written** —
this is a structural preference, not a correctness issue.

---

## 4. Milestones

Following the spec's build order, with prerequisites inserted. One milestone per
session, each confirmed against its tests before the next.

**M0 — Grid and inputs — DONE 2026-10-03.**
`gpbreach.cascade.grid` + `build_grid`. Canonical grid is `bed_ellip.tif`'s:
5 m, ESRI:102247, extended to 5018 x 5854 cells (25.09 x 29.27 km, 29.4 M) with a
500 m halo, phase-aligned to the original (integer cell offset, so existing 5 m
products stay pixel-compatible). All three lakes now 100 % captured. Classes are
**rasterised from the outlines**, not resampled — resampling cannot extend a
raster — and agree with the original `classes.tif` at 100.000 % over its
footprint. Outputs in `data/standard/` with `grid_manifest.json`.

Outstanding from M0, for the owner:
- The **bed covers only 30.4 %** of ice cells on the extended grid (it is still
  the original clip). A pass-search route between LGP and GPL does exist within
  that coverage (verified: one connected traversable component, 110.9 km²), so
  M2 is not blocked — but `competing_outlets` toward Alsek probably is.
- **2,237 cells have negative ice thickness** (surface below bed, to −33.9 m).
  `max(z_s − z_b, 0)` absorbs them, but they indicate bed/surface inconsistency
  worth a look.
- `head_2018.tif` / `margin_2018_lgp.tif` are 2497 cols at 4.999597 x 5.000177 m
  with a 2.1 m origin offset, and the surface rasters have anisotropic
  2.382 x 1.288 m pixels. All are regularised onto the canonical grid now, but
  regenerate them from `data/standard/` before using them as cross-checks.

**M1 — Scaffold.** `cascade/` package, `config/data.yaml` + `params.yaml` in the
spec's schema, `io` extensions (class raster reader, grid-consistency assertion
across bed/surface/classes), output sidecars, pytest wiring. Deliverable: a
config that loads and validates, with a clear error naming any missing input.

**M2 — `seal.py` — DONE 2026-10-05.** See D-013. 16 toy tests plus the
equivalence test; 2018 run gives h_pass 233.262 m, Δ +115.96 m, pass 277 m from
the Phase 1 breach point. Found a 6.7 m datum mismatch in the QGIS margin map.

*(original scope)* Priority-flood pass search returning spill + parent arrays;
`head_surface` reusing `flotation_field`; class barriers; `masked_smooth`;
`overtopping_search`; `competing_outlets`. Toy-grid tests first (spec's table),
then: **an equivalence test asserting the new head/pass code reproduces
`saddle_threshold` on the Phase 1 rasters at `f = k·ρ_i/ρ_w`** — this ties the two
halves of the project together and protects both. Then the 2018 run, with Δ and
pass location compared against your QGIS margin map.

**M3 — `path.py` — DONE 2026-10-06.** See D-014. L_pass 5.624 km, bed sill
58.71 m controls the drawdown. Downstream refinement implemented but disabled by
its own safety check; needs smoothing to work.

*(original scope)* `trace_path`, `refine_downstream` (D8 on head), `profile`,
`summarize`. Outputs: path GeoPackage, profile CSV/PNG, summary JSON.

**M4 — `lakes.py`.** *Requires M0's extended grid.* Bathymetry interpolation,
stage curves, `released_volume`, `gpl_stage_area`. Check `V_w` against your
~0.5 km³ expectation and explain any gap.

**M5 — `inflow_empirical.py`.** The four relations, both hydrograph shapes,
volume-preserving rescale. Spec gives exact expected values at V = 500 Mm³, so
this is strongly testable. Froehlich labelled reference-only.

**M6 — `conduit.py`.** Coupled two-lake ODE, `solve_ivp` with Radau/LSODA,
events, the three spec tests (mass balance, zero-head, V^2/3 scaling). Provisional
2018-geometry run overlaid on the bracket.

**M7 — Hooks.** Entry points for projected surfaces by year, the photogrammetry
DEM as a new `surfaces` entry, Monte Carlo wrapper, GPL routing/outlet stubs.

---

## 5. What blocks me

**Hard blockers (cannot start M2 without):**

1. **The ice surface raster.** Nothing on disk. Which file, which grid?
2. **Canonical grid** — 5 m bed grid or ~2.4 m ArcticDEM grid (§2.3)?
3. **Lake levels** for LGP, GPL, Alsek in ellipsoidal metres, with source and date.
   Note this is also open decision D-007/D-008 for the baseline; whatever you
   supply here should settle both.

**Blocks M4:**

4. **Extended grid** covering the full lakes (§2.2).
5. **LGP bathymetry** — file location, format, interpolation method.

**Can default, flag if you disagree:** `k = 1.0` with `f = k·ρ_i/ρ_w` wired to the
existing constants; connectivity 8 with 4 as a check; `smoothing_m = 0`;
representative N at the pass; `L_pass` for conduit length. These match the spec's
defaults — I will use them and record each in `decisions.md`, but per the spec's
guardrail I will not silently change any of them later.

**Deferred, not needed before M6:** n′ range, S0 range, lake temperature,
background inflows, placeholder `O(h_G)`, photogrammetry metadata.

**Minor:** `lakes_2018.gpkg` row `gp_outline` has a stale `path` attribute
pointing at the LGP source layer. Geometries and `class` values are correct;
cosmetic only, but worth fixing at source so provenance stays readable.

---

## 6. Suggested first move

M0 is yours and blocks everything. If it would help, the fastest unblock is for
me to write the grid-extension and surface-preparation script so M0 becomes a
scripted, re-runnable step rather than manual QGIS work — that also makes it
reproducible when the photogrammetry DEM lands. Say the word and that becomes M0.
