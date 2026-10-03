# Grand Plateau LGP→GPL Breach: Pass Search & Inflow Hydrograph Build Spec

Oct 3, 2026 · @Krishna

## Purpose and context

Build a Python pipeline that finds when and where the subglacial seal between Little Grand Plateau Lake (LGP) and Grand Plateau Lake (GPL) opens, then produces the LGP→GPL inflow hydrograph.

- **The cascade.** LGP sits \~100 m above GPL behind Grand Plateau Glacier. LGP's level is treated as constant; the dam thins until it reaches flotation and LGP drains under the ice into GPL. That flood loads GPL's boulder-armored outlet, which may fail and turn GPL into an ocean bay.
- **Scope of this spec.** Seal geometry and onset, drainage path and profile, lake curves and released volume, and the inflow hydrograph (empirical bracket plus a physical conduit model).
- **Later stages, not built yet.** Routing through GPL, outlet hydraulics and armor stability, breach development, and the Monte Carlo likelihood. Outputs here must plug into those stages.
- **Key physics choice.** Flotation is expressed through the Shreve hydraulic potential with a factor k. k = 1 is pure flotation; hindcasts suggest pure flotation is an upper bound, so k < 1 is allowed.
- **What not to use for this phase.** Froehlich (2025) is a moraine-breach regression and does not apply to ice-dam drainage. It is reserved for the later GPL outlet breach.

## Inputs and conventions

All inputs are already prepared in QGIS on one grid (the bed grid) and in WGS84 ellipsoidal heights.

| Input | File | Notes |
| --- | --- | --- |
| Bed elevation | `bed_ellip.tif` (also `bed_filled.tif`) | IfSAR minus Millan thickness; NAVD88 → ellipsoid via GEOID12B plus stable-ground offset; values below −500 m removed and interpolated; nodata outside glacier |
| Ice surface 2018 | `surf_50m.tif` (also `surf_filled.tif`) | ArcticDEM WV02 + WV03 strips, 2018-09-03; lakes masked, clipped to ice, Average-resampled to the bed grid |
| Classes | `classes.tif` | 0 land, 1 ice, 2 LGP, 3 GPL, 4 Alsek; land may be nodata |
| Outlines | `ice_2018.gpkg`, `lakes_2018.gpkg` | Digitized from Sentinel-2, 2018-09-05; lakes carry a `class` field |
| Lake levels | config values | LGP, GPL, Alsek in ellipsoidal metres; IfSAR + N + same offset as the bed |
| LGP bathymetry | plumb-line points | Depths reported as ">275" were set to 275 m; format to confirm |

Conventions:

- **Read originals, not the `_filled` rasters.** In Python, read `bed_ellip.tif` and `surf_50m.tif` with nodata as NaN. The `_filled` versions use 0 as a placeholder for the QGIS Raster Calculator and must not be treated as elevations.
- **Units.** SI throughout; elevations in metres; time in seconds internally, reported in hours.
- **Constants.** ρ\_i = 917 kg m⁻³, ρ\_w = 1000 kg m⁻³, g = 9.81 m s⁻², ρ\_i/ρ\_w = 0.917.
- **Known caveats.** The two ArcticDEM strips disagree by 0.3–2 m in their overlap; the 2018 surface is unsmoothed; a \~1 m residual datum offset may remain.

## Repository layout and coding conventions

Every data source and parameter lives in config files, so a new DEM (such as the upcoming photogrammetry surface) is added by editing config, not code.

```
grand_plateau/
├── config/          data.yaml, params.yaml
├── data/raw/        original inputs (read-only)
├── data/standard/   QGIS-prepared rasters and outlines (inputs above)
├── data/derived/    outputs of this pipeline, one folder per run
├── src/gp/          io.py, seal.py, path.py, lakes.py,
│                    inflow_empirical.py, conduit.py
├── tests/           synthetic test cases
└── notebooks/       exploration and figures
```

Example `data.yaml` entries:

```yaml
surfaces:
  - id: arcticdem_20180903
    path: data/standard/surf_50m.tif
    date: 2018-09-03
    vertical: ellipsoid_wgs84
bed: {path: data/standard/bed_ellip.tif, vertical: ellipsoid_wgs84}
classes: data/standard/classes.tif
lake_levels: {lgp: null, gpl: null, alsek: null}   # metres, ellipsoidal
```

Example `params.yaml` entries:

```yaml
k: 1.0                # flotation factor; hindcast-calibrated range later
smoothing_m: 0        # masked Gaussian sigma in metres; 0 = off
connectivity: 8
conduit: {manning_n: 0.08, S0_m2: 0.1, lake_temp_C: 0.5, A: 2.4e-24, n_glen: 3}
```

Conventions for Claude Code:

- **Stack.** Python 3.11+, numpy, scipy, pandas, rasterio, geopandas, shapely, pyyaml, matplotlib; numba optional for the pass search.
- **Structure.** Each module exposes pure functions plus a `run(cfg)` that reads inputs and writes outputs; run with `python -m gp.<module> --config config/data.yaml`.
- **Outputs.** GeoTIFF, GeoPackage, CSV, and JSON summaries. Each output gets a JSON sidecar with the config used, input file paths, and the git commit.
- **Docstrings state units** for every argument and return value.
- **Guardrails.** Never write to `data/raw/` or `data/standard/`. Ask before changing a physical assumption (k, connectivity, the N definition, flow law) rather than silently picking one.

## Module: seal.py (head surface and pass search)

`seal.py` computes the hydraulic head surface and finds the pass head h\_pass, the lowest possible "highest head" on any subglacial route from LGP to GPL. The seal opens when h\_LGP ≥ h\_pass.

Head under ice, in metres of water:

```latex
h_\phi = z_b + k\,\frac{\rho_i}{\rho_w}\,\max(z_s - z_b,\ 0)
```

Functions:

1. `load_inputs(cfg)` returns bed, surface, classes (NaN for nodata), the affine transform, and cell size.
2. `masked_smooth(z, mask, sigma_m, res)` smooths ice cells only, by normalized convolution: `gaussian(z·m) / gaussian(m)`. Returns NaN off the mask. `sigma_m = 0` returns the input unchanged.
3. `head_surface(zb, zs, classes, lake_levels, k)` returns the head array. Ice cells use the formula above; lake cells take their lake level; land and ice cells with NaN inputs become barriers (`+inf`), and NaN ice cells are counted and logged.
4. `pass_search(head, classes, source_class=3, target_class=2, connectivity=8)` runs a minimax priority flood (heapq):
   - Seed every source-class cell with its head.
   - Pop the lowest; for each non-barrier neighbour, candidate = max(current, neighbour head); keep it if lower than the stored value and record the parent index.
   - Set target-class cells to −inf before the search so their stored value equals h\_pass.
   - Return the spill array, parent array, h\_pass, and the pass cell (highest-head cell on the traced path).
5. `overtopping_search(...)` runs the same search on ice surface elevation instead of head, for comparison with the flotation onset.
6. `competing_outlets(...)` reruns the search seeded from Alsek (class 4) and reports which receiving lake gives LGP the lowest pass.

Outputs per surface: `head_<id>.tif`, `margin_lgp_<id>.tif` (head − h\_LGP on ice), `spill_<id>.tif`, and `seal_summary_<id>.json` with h\_pass, Δ = h\_pass − h\_LGP, pass coordinates, and thinning needed at the pass (Δ divided by k·ρ\_i/ρ\_w).

Notes:

- **Connectivity.** 8-connectivity lets water pass diagonally between two land cells that touch only at a corner. Default to 8 and expose 4 as an option for a conservative check.
- **Speed.** A few hundred thousand cells runs in seconds with heapq. Add numba only if per-year loops get slow.
- **Check against QGIS.** For 2018, Δ must be positive and close to the value read off the QGIS margin map.

## Module: path.py (drainage path and along-path profile)

`path.py` turns the pass search into one drainage path and the handful of numbers the flood model needs.

Functions:

1. `trace_path(parent, classes)` follows parent links from the LGP cell reached by the search back to GPL.
2. `refine_downstream(head, pass_cell, classes)` replaces the segment below the pass with steepest descent (D8) on the head surface to GPL. Water flows downhill in potential once the seal opens, while the minimax path below the pass is not unique.
3. `profile(path, zb, zs, head, h_lgp, h_gpl)` samples along the path and returns a table:
   - distance s (m), x, y
   - z\_b, z\_s, H, head (m)
   - N\_lake = head − h\_LGP and N\_gpl = head − h\_GPL (effective pressure in metres of water, bracketing the conduit's water pressure)
4. `summarize(profile)` reports:
   - L\_margin (LGP margin to GPL) and L\_pass (pass to GPL), in m
   - bed sill = max z\_b along the path; this sets the lowest level LGP can drain to
   - flotation-zone length (where head < h\_LGP)
   - H at the pass, and min / mean H along the path
   - representative N (default: N\_lake at the pass; alternatives in Open decisions)

Outputs: `path_<id>.gpkg` (LineString with summary attributes), `profile_<id>.csv`, `path_summary_<id>.json`, and `profile_<id>.png` (bed, surface, head and both lake levels against distance).

## Module: lakes.py (stage–area–volume curves and released volume)

`lakes.py` builds both lakes' stage curves and computes the released volume Vw, which the notes put near 0.5 km³.

Functions:

1. `interpolate_bathymetry(points, lake_mask, h_lake)` converts plumb-line depths to bottom elevation (h\_LGP − depth) and grids them on the analysis grid (linear TIN by default). Capped ">275 m" points likely don't matter, because drainage removes at most \~100 m from the top.
2. `stage_curve(bottom, mask, levels)` returns A(h) and V(h) by counting cells below each stage, with V accumulated from the lake floor.
3. `released_volume(curve, h_lgp, bed_sill, h_gpl)` sets h\_final = max(h\_GPL, bed sill) and returns Vw = V(h\_LGP) − V(h\_final). It also flags whether the bed sill or GPL controls.
4. `gpl_stage_area(dem, gpl_mask, h_gpl, rise_max=30)` gives GPL's area for stages above today's level. At each stage it counts cells below that stage that are connected to the lake (flood fill), so separate basins are excluded.

Outputs: `lgp_hypsometry.csv` (h, A, V), `gpl_hypsometry.csv`, and `volume_summary.json` (h\_final, the controlling constraint, Vw in m³ and Mm³).

Check: plot A(h) and V(h) for both lakes; Vw should be the same order as 0.5 km³ or the difference should be explainable.

## Module: inflow\_empirical.py (empirical bracket)

`inflow_empirical.py` gives a first inflow range from Vw alone and serves as the sanity check on the conduit model.

| Relation | Peak discharge Qp (m³/s), V in Mm³ | Drainage mode | Qp at V = 500 Mm³ |
| --- | --- | --- | --- |
| Walder & Costa (1996), tunnel | 46 V^0.66 | Subglacial tunnel | \~2,800 |
| Clague & Mathews (1973) | 75 V^0.67 | Subglacial tunnel | \~4,800 |
| Walder & Costa (1996), non-tunnel | 1100 V^0.44 | Breach, marginal, or mechanical failure | \~17,000 |
| Froehlich (2025), reference only | 142 Hw^0.766 V^0.378 | Moraine breach | \~51,000 |

Hydrograph shapes:

- **Tunnel drainage.** Slow rise, abrupt fall. Use a reversed lognormal (σ = 0.5), as in published hazard studies.
- **Non-tunnel drainage.** Triangle with base time t\_b = 2 Vw / Qp.
- **Volume check.** Rescale each curve numerically so its area equals Vw exactly.

Outputs: `inflow_empirical.csv` (t and Q for each relation) and `inflow_empirical.png`.

The Froehlich row is a reference point only and should be labelled that way in plots.

## Module: conduit.py (physical inflow hydrograph)

`conduit.py` solves a lumped Nye/Clarke-type conduit model coupled to both lakes and returns the LGP→GPL hydrograph. Both lakes are solved together, because GPL's rise lowers the head driving the flood.

State: LGP level h\_L, GPL level h\_G, conduit cross-section S.

```latex
A_{L}(h_L)\,\frac{dh_L}{dt} = Q_{in,L} - Q
```

```latex
A_{G}(h_G)\,\frac{dh_G}{dt} = Q + Q_{in,G} - O(h_G)
```

```latex
Q = \frac{S^{4/3}}{n'\,(4\pi)^{1/3}}\left(\frac{h_L - h_G}{L}\right)^{1/2}
```

```latex
\frac{dS}{dt} = \frac{Q\,\rho_w g\,(h_L - h_G)/L}{\rho_i L_f} + \frac{\rho_w c_w\,Q\,T_L}{\rho_i L_f\,L} - 2A\,S\left(\frac{N}{n}\right)^{n}
```

Terms and defaults:

- **Flow law.** Manning flow in a full circular conduit; n′ sampled around 0.03–0.12 (to confirm against the literature).
- **Melt.** The first term is wall melting by viscous heat. The second is lake sensible heat T\_L, applied as a lumped approximation over length L; check it against the Clarke (2003) form.
- **Creep closure.** A = 2.4 × 10⁻²⁴ Pa⁻³ s⁻¹, n = 3. N = ρ\_w g × max(head\_pass − (h\_L + h\_G)/2, 0), updated as the lakes change.
- **GPL outflow O(h\_G).** A placeholder rating or constant base outflow for now; the routing stage replaces it later.
- **Initial state.** h\_L = h\_LGP, h\_G = h\_GPL, S0 sampled log-uniform over 0.01–1 m².
- **Constants.** L\_f = 3.34 × 10⁵ J kg⁻¹, c\_w = 4186 J kg⁻¹ K⁻¹, L from `path_summary` (L\_pass by default).

Solver:

- `solve_ivp` with Radau or LSODA, rtol 1e-6, and `max_step` small enough to resolve the peak.
- Events stop the run when h\_L reaches h\_final, when S closes below a small threshold after the peak, or at t\_max.
- Report whether the run ended at the sill or by closure, since closure means partial drainage.

Outputs: `inflow_conduit_<scenario>.csv` (t, Q, h\_L, h\_G, S), and a summary with Qp, time to peak, duration, released volume, and the stop reason. Overlay each hydrograph on the empirical bracket in one plot.

Later extension (stub only): the Walder & Costa (1996) breach-widening model for the mechanical-failure scenario.

## Testing and validation

Each module ships with synthetic tests that have known answers, plus QA figures for the real 2018 run.

| Module | Test | Pass condition |
| --- | --- | --- |
| seal | Toy grid with one known saddle between two basins | h\_pass and pass cell match exactly |
| seal | Land wall fully separating the lakes | h\_pass = +inf, no path |
| path | Path on the toy grid | Ends in GPL; head non-increasing below the pass |
| lakes | Synthetic cone-shaped basin | A(h) and V(h) within 1% of analytic values |
| inflow\_empirical | Vw = 500 Mm³ | \~2,780 / \~4,820 / \~16,900 m³/s; each curve integrates to Vw |
| conduit | Mass balance | ∫Q dt equals the LGP volume lost within 0.1% |
| conduit | h\_L = h\_G at start | Q stays zero |
| conduit | Sweep Vw over two orders of magnitude | Qp scales roughly as V^(2/3) (Ng & Björnsson 2003) |

QA outputs for each run:

- Margin map PNG with the path and pass marked; compare with the QGIS margin map.
- Profile plot of bed, surface, head and lake levels along the path.
- Hydrograph plot overlaying conduit runs on the empirical bracket.
- A short `run_report.md` listing inputs, parameters, h\_pass, Δ, Vw, Qp and any warnings (NaN ice cells, closure-limited runs).

## Build order

Build in this order, one milestone per Claude Code session, and confirm each against its tests before starting the next.

1. **Scaffold.** Repo layout, `data.yaml` and `params.yaml`, `io.py` (raster read with NaN nodata, grid checks, output sidecars), pytest setup.
2. **seal.py.** Toy-grid tests first, then the 2018 run. Compare Δ and the pass location with the QGIS margin map.
3. **path.py.** Path, profile table and profile plot for 2018.
4. **lakes.py.** Needs the bathymetry format confirmed. Produce both hypsometries and Vw.
5. **inflow\_empirical.py.** Bracket hydrographs from Vw.
6. **conduit.py.** Tests first, then a provisional run using 2018 geometry as if the seal opened today. Overlay on the bracket.
7. **Hooks for later.** Leave clean entry points for projected surfaces by year (the thinning module), the photogrammetry DEM as a new `surfaces` entry, the Monte Carlo wrapper, and the GPL routing and outlet modules.

The 2018 hydrograph is provisional. The real onset geometry comes from projected surfaces at the onset year.

## Open decisions

These need an answer from Krishna before or during the build. Claude Code should ask rather than guess.

- [ ] Lake levels for LGP, GPL and Alsek in ellipsoidal metres, with sources and dates
- [ ] LGP bathymetry format (raw points or an interpolated surface) and the interpolation method
- [ ] Smoothing length for the ice surface (start at 0, then test a few ice thicknesses)
- [ ] k range, calibrated from the three hindcast events
- [ ] Connectivity: 8 (default) or 4
- [ ] Representative N: at the pass, mean over the seal segment, or minimum along the path
- [ ] Path length used in the conduit model: L\_pass or L\_margin
- [ ] LGP water temperature and seasonal background inflows
- [ ] Conduit roughness n′ range and initial size S0 range
- [ ] Placeholder GPL outflow O(h\_G) for the conduit phase, until the routing module exists
- [ ] Timeline and metadata for the photogrammetry DEM (date, CRS, vertical datum, ground control)
- [ ] Workflow manager: plain module runs for now, or Snakemake from the start
