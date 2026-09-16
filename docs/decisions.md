# Decision log

Every physics, input, or methodological choice. Claude Code appends entries
**only after owner approval** (CLAUDE.md Section 6).

All entries below are seeded as `STATUS: TO CONFIRM` — they record what the
baseline currently does, recovered from the owner's own outputs. None has been
approved yet, and nothing here was chosen by Claude Code: each value was
measured from the legacy rasters or stated by the owner.

---

## D-001 — Breach criterion: flotation plus a connected pathway
- Date: 2026-09-16
- Decided by: owner (reconstructed; **STATUS: TO CONFIRM**)
- Context: what counts as a breach.
- Decision: a cell is buoyant when ice above lake level falls to a fraction `r`
  of ice below lake level; a breach occurs the first year a connected pathway of
  buoyant cells links the two lake basins. Implemented as
  `components/criterion.py: FlotationConnectivity`.
- Affects: `components/criterion.py`, `saddle.py`, both baseline configs.
- Expected effect: reproduces 2030 and 2035 exactly.

## D-002 — Densities are not used by the baseline
- Date: 2026-09-16
- Decided by: owner (reconstructed; **STATUS: TO CONFIRM**)
- Context: CLAUDE.md Section 5 asks for margins in m w.e., which presumes the
  flotation fraction comes from `rho_ice / rho_water`.
- Decision: the baseline prescribes `f` directly, so `RHO_ICE = 917` and
  `RHO_WATER = 1000` are defined in `constants.py` but unused. Margin is
  reported in metres as `f*H - h_lake`, algebraically the same comparison.
- Note for the owner: the density ratio implies `f = 0.917`, which is *higher*
  than both 0.90 and the 0.9091 actually used. Worth deciding whether `f` is
  meant to be a density ratio, or a calibrated dam-strength parameter that
  absorbs ice fracture and bed friction (the hindcast phase would inform this).

## D-003 — Flotation fraction: `r = 0.1`, i.e. `f = 0.909091`
- Date: 2026-09-16
- Decided by: owner (reconstructed; **STATUS: TO CONFIRM**)
- Context: the owner described the target as "lake level 90% of ice dam
  thickness", then implemented "above-lake ice < 10% of below-lake ice".
- Options: `r = 0.1` (what the rasters do, `f = 0.9091`); `r = 1/9 = 0.1111`
  (what `f = 0.90` exactly requires).
- Decision: keep `r = 0.1` as built. Not "fixed" (CLAUDE.md Section 8.1).
- Expected effect: `f = 0.90` would move both dates ≈0.4 yr earlier and change
  neither integer year. `gpbreach fcurve` shows the full `f` dependence.

## D-004 — Surface change: uniform -9.05 m/yr, linear
- Date: 2026-09-16
- Decided by: owner (reconstructed; **STATUS: TO CONFIRM**)
- Context: the owner recalled "about -9 m/year" from OIB/DEM transects; the
  rasters step by exactly 9.05 m/yr.
- Decision: `-9.05 m/yr`, spatially uniform, linear, no seasonality.
- Affects: `components/surface.py`, both configs.
- **This is the parameter the answer is most sensitive to by a wide margin**
  (`gpbreach sensitivity`: ~2.2–2.8 yr per 1 m/yr, versus 0.2 yr per 0.01 of `f`).
  Replacing it with a calibrated model is the seasonal-thinning work already
  under way with PyGEM/OGGM.
- **2026-09-16, owner:** an updated rate is coming; stay on the current value for
  now. Kept at **-9.05**, not -9.00. The owner refers to the rate as "-9 m/yr",
  but the published rasters encode 9.05 exactly, and the cell-for-cell
  regression test depends on it. For the record, -9.00 would give 2029.745 and
  2034.739 — the same integer years, ~0.11 yr later. When the new rate lands,
  change it in the configs, expect `test_matches_legacy_rasters_every_year` to
  fail by design, and re-baseline deliberately rather than loosening the test.

## D-005 — Basin seeds and which pair defines each breach
- Date: 2026-09-16
- Decided by: owner — **CONFIRMED by description 2026-09-16**, pending a QGIS check
- Context: the legacy connectivity judgement was made visually in QGIS and is not
  recorded in any file. The port needs two explicit seed cells per breach.
- Decision (provisional): seeds are the minimum-`Phi` cell of each of the three
  dominant basins — A (1903, 649), B (122, 22), C (2486, 2450). Breach 1 is A–C
  with L = 110 m; breach 2 is A–B with L = 17 m.
- **Owner's independent description of the two pathways** (2026-09-16):
  - breach 1, LGP→GP: "the slope that is west of the lake which the glacier wraps
    up and around";
  - breach 2, Alsek→GP: "the slope to the south of the lake where the glacier
    wraps around".
  Both match the saddles the model finds, which were located before the owner
  described them:
  - breach 1 crest at (908752.5, 1112597.5) — on the **west** flank of the upper
    (LGP) basin;
  - breach 2 crest at (903622.5, 1118227.5) — at the **south** tip of the Alsek
    basin, where the two lobes pinch together.
  `gpbreach map` renders both (`results/runs/<run_id>/breach_map_*.png`) and
  exports `breach_points_*.csv` / `breach_pathway_*.csv` for overlay in QGIS.
- **Residual concern:** the seeds for the Alsek and upper-LGP basins sit within
  ~25 cells of the raster edge, so those basins may be clipped by the domain.
  This does not affect the saddles found (both crests are well inside the
  domain), but a wider clip could in principle expose a lower pass. Worth one
  QGIS check against the full-extent DEM.
- Affects: both configs, `preprocess/build_dam.py`.
- Expected effect: this is what selects the saddle, hence the date. If a seed is
  in the wrong basin, the date is wrong.

## D-006 — Bed treated as exactly known
- Date: 2026-09-16
- Decided by: owner (reconstructed; **STATUS: TO CONFIRM**)
- Context: `Bedrock_Aligned.tif` is used as-is.
- Decision: no perturbation, no error model in the baseline. Spatially
  correlated bed error is Phase 3.
- Note: a *uniform* bed bias barely moves the answer (0.011 yr per m, because
  only the `(1-f)` weight touches the bed). A *spatially correlated* error could
  matter much more, since it can move the saddle to a different pass. Do not
  infer from the uniform-bias number that bed uncertainty is negligible.
- Provenance of `Bedrock_Aligned.tif` is not documented anywhere in the repo —
  see `docs/data_sources.md`.

## D-007 — Lake levels: fixed at 110 m and 17 m
- Date: 2026-09-16
- Decided by: owner (reconstructed; **STATUS: TO CONFIRM**)
- Context: both values recover exactly from the rasters, to 4 dp.
- Decision: constant in time for each breach.
- Open: is each a fixed sill, an ice spillway that falls as the dam thins, or a
  river stage? A spillway that tracks the thinning surface would *delay* breach;
  a fixed sill is the more conservative choice. Sensitivity is
  -0.12 yr per metre of lake level.

## D-008 — Vertical datum: the two grids are on different references
- Date: 2026-09-16
- Decided by: **owner decision still required**, but the facts are now established
- Context: `2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif` and `Bedrock_Aligned.tif` are
  differenced cell-by-cell, which is only valid on a common vertical datum.
  Neither file declares one.

### Evidence gathered 2026-09-16
1. **Neither raster carries a vertical CRS.** `gdalinfo` shows a 2D projected CRS
   only. The Loso lidar tiles do declare `VERTCRS`, but as `"unknown"`.
2. **The IFSAR surface is NAVD88 orthometric.** USGS states that Alaska IFSAR is
   delivered in NAD83 horizontal / NAVD88 vertical, metres, at 5 m post spacing —
   which matches this grid exactly.
3. **`Bedrock_Aligned.tif` is tied to the OIB radar bed.** Sampled at all 932 OIB
   points inside the domain, `bed_h_m` minus `Bedrock_Aligned` has a median of
   **+0.04 m** (IQR -0.32 to +0.54, sd 1.30). The raster was evidently built by
   aligning to those picks — hence the name — so it inherits the OIB datum.
   IceBridge L2 products are normally WGS84 ellipsoid heights; **this is the one
   fact still to confirm, from the NSIDC user guide for `IRUAFHF1B`.**
4. **The offset is a constant.** GEOID12A geoid height over the domain is
   +6.66 m (NW corner) to +6.86 m (SE corner), i.e. **+6.80 m ± 0.10 m** — a
   0.21 m gradient across 12 km. A single scalar correction is sufficient; no
   geoid grid is needed.

### Effect on the answer (measured, not estimated)
With a 6.80 m mismatch, depending on which dataset is the odd one out:

| case | shift | integer year |
|---|---|---|
| bed ellipsoidal, surface and lake level orthometric | **-0.075 yr** | 2030 / 2035 — unchanged |
| surface ellipsoidal, bed and lake level orthometric | -0.751 yr | 2029 / 2034 |
| lake level on the other datum from both grids | -0.826 yr | 2029 / 2034 |

The **first row is the likely configuration** given the evidence above: the bed
came from OIB, the surface is IFSAR, and the lake levels were presumably read off
the IFSAR-based QGIS project. In that case the datum mismatch is worth 0.08 yr
and changes neither headline date.

### What is still open
- Confirm the OIB vertical datum from the NSIDC user guide (item 3).
- Confirm which datum the lake levels of 110 m and 17 m were read in. This is the
  single input with the largest datum leverage (-0.12 yr per metre).
- Decide whether to correct the bed by -6.80 m for internal consistency. This is
  a physics/input change and needs owner approval; the baseline currently applies
  **no** correction, matching the owner's original rasters.

### Aside — a fill-value trap in the IFSAR raster
`2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif` has **no nodata flag set** and uses
`0.0` as fill for 3,570,204 of its 6,347,328 cells (56%). None of them fall
inside the current analysis domain, which is masked by the bedrock raster's
proper nodata, so the baseline is unaffected — verified, and `preprocess` now
asserts it. But any future widening of the domain, or any use of this raster on
its own, must mask zeros explicitly or it will read 0 m fill as real ground.

## D-009 — Horizontal CRS stays Alaska Albers (ESRI:102247)
- Date: 2026-09-16
- Decided by: owner — **CONFIRMED 2026-09-16** ("the Alaska Albers datum is the
  one to use")
- Context: CLAUDE.md Section 5 specifies EPSG:3338 (NAD83 Alaska Albers). All
  project rasters are ESRI:102247 (NAD83(CORS96) Alaska Albers) — the same
  projection on a different realisation of the datum.
- Decision: no reprojection. Every input is already co-registered on one 5 m
  grid, and reprojecting would resample the owner's data (CLAUDE.md Section 2.2).
- Note: the vector layers are mixed — `outletStrengthCalculations/lakePolygon.gpkg`
  and `OIB-and-Hig-depths-merged.gpkg` are EPSG:32607 (UTM 7N), while
  `outletOutline.gpkg` is in Albers but lies entirely *outside* the raster extent.
- **Consequence for later phases:** the owner plans to lean on OIB lidar, which
  is UTM 7N. Albers stays the analysis CRS, so `ingest/` must reproject UTM →
  Albers on the way in, and must never reproject the analysis grid itself. The
  two Albers variants in play differ only in datum realisation — ESRI:102247 is
  NAD83(CORS96), EPSG:3338 is NAD83 — a sub-metre horizontal difference, but it
  should be stated explicitly per dataset rather than assumed equivalent.
- **Horizontal CRS is settled; vertical datum is not.** See D-008, which the
  OIB plan makes more urgent, not less.

## D-010 — Interpolating `Phi_crit(f)` for large ensembles
- Date: 2026-09-16
- Decided by: Claude Code (performance only; **STATUS: TO CONFIRM**)
- Context: a saddle solve takes ~0.5 s, so a 10,000-member ensemble over `f`
  would take over an hour.
- Decision: `analysis/fcurve.py` tabulates `Phi_crit` on an `f` grid and
  interpolates. Exact at the nodes; measured interpolation error is 0.09 m
  (≈0.01 yr) at 0.005 spacing.
- This is a numerical device, not a physics change, but it is recorded here
  because it sits between the sampler and the answer. The engine does **not**
  use it by default — it solves exactly, once per distinct `f`.
