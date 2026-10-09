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

### ICESat-2 ATL13 evidence, 2026-09-25
Granule `ATL13_20260719041126_05373201_007_01.h5` (v007, acquired 2026-07-19),
read with `gpbreach atl13`. Its own dataset attributes declare:

| field | declared reference |
|---|---|
| `ht_water_surf` | **WGS84 ellipsoid** |
| `ht_ortho` | **EGM2008** orthometric |
| `segment_geoid` | EGM2008 geoid, **mean-tide** system |

**EGM2008 is not NAVD88.** At 59.2318 N, 138.0970 W the granule's EGM2008 geoid
is 8.679 m while NOAA's GEOID12A is 7.483 m — a **1.196 m** disagreement. Using
`ht_ortho` as if it were NAVD88 would inject ~1.2 m of lake-level error, i.e.
~0.14 yr of breach date. Two further caveats on any conversion: ICESat-2 is
ITRF2014 while GEOID12 expects NAD83 (order a metre in Alaska, needs HTDP/NCAT
to do properly), and ATL13's geoid is mean-tide while NAVD88 is tide-free
(`segment_geoid_free2mean` = -0.155 m here).

**The granule is the right one; ATL13 is the wrong product.** The file is named
for RGT 0537 but spans ~6.3 h of orbit and actually contains RGTs 537-540. **RGT
540 crosses the site**: fitted from its own nearby segments, beam pair gt1
passes **0.41 km** west of Grand Plateau Lake and **0.47 km** west of the modelled
Alsek basin; gt3 passes 1.4 km from the upper LGP lake. CMR independently
confirms ATL06/ATL03 granules on RGT 0540 intersect a tight box on the lakes.

Yet ATL13 reports **zero water segments on all six beams** between 58.95 N and
59.20 N — precisely the band holding Grand Plateau Lake (59.0346 N) and the
Alsek basin (59.1187 N). The nearest reported water is 20 km away.

The cause is ATL13's reference water mask. Every water body it does report here
carries `inland_water_body_source = 1` (**HydroLAKES**), a static global lake
database that does not contain these ice-dammed proglacial lakes. No ATL13
granule on any RGT will help until that mask includes them.

**Use ATL06 (land ice height) instead**, which is masked to glacier ice rather
than to a lake inventory and therefore does cover them. The matching granule for
the same day and track is `ATL06_20260719090121_05403202_007_01.h5` (RGT 0540,
cycle 32). ATL03 is the fallback if ATL06 turns out to be masked out over open
water. ATL06 heights are expected to be WGS84 ellipsoid with a geoid field
alongside, but **verify with the same read-it-from-the-file approach** rather
than assuming — that is what `gpbreach atl13` does and what an ATL06 reader
should do too.

**ATL06 outcome, 2026-09-25.** `ATL06_20260719090121_05403202_007_01.h5`
(RGT 0540, cycle 32) was obtained. Its declared reference, read from the file:
`h_li` is **WGS 84 ellipsoid, ITRF2014**, and `dem/geoid_h` is EGM2008 **in the
tide-free system** — the opposite tide convention to ATL13's mean-tide
`segment_geoid`. NAVD88/GEOID12 is tide-free, so ATL06 matches it and ATL13 does
not; each product carries its own `geoid_free2mean` and they apply it in
opposite directions.

**ATL06 works over water here.** Over ATL13's water body 5001, on the same track
and the same day, ATL06 returns 159/167 and 134/145 good-quality segments and
agrees with ATL13 to **0.10–0.25 m** (32.50 / 32.64 m vs 32.745 m ellipsoid).
That cross-validates both products and removes any doubt that the land-ice
product can measure a lake surface.

**But this granule is unusable over the study area.** Of the 1382 segments
between 58.95 N and 59.13 N — the analysis domain — **zero** pass
`atl06_quality_summary`, and the median disagreement with ATL06's own reference
DEM is **597 m**. Band statistics show why: 52% good north of 59.13 N and 92%
good south of 58.70 N, but 0% good across 58.70–59.13 N, with `cloud_flg_atm`
rising from 0 to 1 over exactly that interval. It is cloud, not terrain, not
water, and not the track. A different date is required; no processing choice
recovers these heights.

**Consequence for the code.** `gpbreach.io.atl06` screens on
`atl06_quality_summary == 0` **and** `|h_li - dem_h| <= 50 m` by default. Without
that screen these segments read as plausible numbers and would have implied a
~700 m elevation error. A regression test asserts they are rejected.

**Available record.** CMR lists 77 ATL06 granules intersecting the lakes on RGTs
0060, 0502, 0540 and 0982, spanning 2018-10-31 to 2026-07-19 at roughly nine per
year. That is enough for a lake-level time series, which would settle the datum
question and also feed the Breach 1 to Breach 2 cascade work directly.

**Note on the earlier comparison.** Sampling the IFSAR DSM at the ICESat points
of water body 5001 returns exactly 21.000 m everywhere (sd 0.000) — an IFSAR
water-flattening constant, not a measurement. No datum offset can be read from
it, and body 5001 is 9.2 km north of the raster anyway, identity unconfirmed.

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

## D-011 — RGI 7.0 Alaska glacier outlines added as an input
- Date: 2026-09-27
- Decided by: owner — **approved by direct request 2026-09-27**
- Context: the owner asked for RGI v7 outlines near Grand Plateau for manual
  analysis. Downloading a dataset and writing to `data/raw/` are both
  stop-and-ask items (CLAUDE.md Section 2.2); the request is the approval.
- Decision: downloaded `RGI2000-v7.0-G-01_alaska.zip` (83 MB, 27,509 glaciers)
  from the NSIDC DAAC to `data/raw/RGI7_Alaska/`, using the Earthdata
  credentials already in the owner's `~/.netrc`. Nothing else in `data/raw/`
  was touched. `workflow/subset_rgi.py` writes a reprojected local subset to
  `data/interim/`; it never writes to `data/raw/`.
- Citation: RGI 7.0 Consortium, 2023. Randolph Glacier Inventory Version 7.0.
  NSIDC. doi:10.5067/f6jmovy5navz
- Source CRS EPSG:4326; the subset is reprojected to ESRI:102247 per D-009.

### The two dam glaciers
| RGI id | name | area | z range |
|---|---|---|---|
| `RGI2000-v7.0-G-01-17002` | unnamed in RGI; Loso's "Grand Plateau North" | 182.55 km² | 19–3219 m |
| `RGI2000-v7.0-G-01-27357` | **Grand Plateau Glacier** | 237.25 km² | 22–4643 m |

Both are `term_type = 9`. Neighbours within 20 km include Alsek (`-16980`),
Melbern (`-17004`), Fairweather (`-27011`) and Ferris (`-27355`); 234 glaciers
totalling 2,367 km² are in the subset.

### What this says about D-005
Clipping the two RGI polygons against the project rasters:

| glacier | inside the analysis domain | fraction of the glacier |
|---|---|---|
| `-17002` | 14.54 km² | **8 %** |
| `-27357` | 54.08 km² | **23 %** |

Their combined 68.6 km² is essentially the whole modelled domain (69.3 km² of
valid bedrock), which confirms the domain *is* these two glaciers and nothing
else. But **77–92 % of each glacier lies outside the clip.** That is direct
evidence for the residual concern recorded under D-005: the flotation basins are
truncated by the raster extent, and a wider clip could in principle expose a
lower saddle than the ones driving the 2030 and 2035 dates. The RGI outlines now
give an objective boundary against which to test that, should the owner want the
rasters re-clipped.

## D-012 — The four lake-terminating glaciers
- Date: 2026-09-27
- Decided by: owner identified them by eye; **selection below is Claude Code's
  match to RGI ids and needs owner confirmation against imagery**
- Context: the owner wants the four glaciers in contact with the three lakes —
  Grand Plateau, the glacier north of it, and two calving into the NE arm of
  Alsek Lake.
- Method: filtered the 234-glacier subset to area > 5 km² and zmin < 200 m, then
  ranked by polygon distance to each lake reference (upper LGP lake polygon,
  Grand Plateau basin seed, and the ICESat-2 track across Alsek Lake).

| rgi_id | name | area | zmin | lake |
|---|---|---|---|---|
| `-01-27357` | Grand Plateau Glacier | 237.25 km² | 21.9 m | Grand Plateau Lake |
| `-01-17002` | unnamed (Loso: GP North) | 182.55 km² | 19.0 m | GP / upper lake |
| `-01-16987` | unnamed | 115.28 km² | 21.6 m | Alsek Lake (NE) |
| `-01-16980` | Alsek Glacier | 107.86 km² | 23.8 m | Alsek Lake (NE) |

All four terminate between **19.0 and 23.8 m**, consistent with calving into
lakes whose surfaces sit near 17–25 m — the strongest evidence for the
selection. `-01-16966` (40 km²) is closer to Alsek Lake than two of the four but
has zmin 63 m, well above the water surface, so it was rejected as not calving.
Confirm the NE-arm pair against Sentinel-2 before relying on it.

- **All four carry `src_date = 2010-09-15`.** Despite the "RGI2000" name these
  outlines are 2010 imagery — the same epoch as the IFSAR DEM, which is
  convenient, but ~16 years stale at the calving fronts.
- `workflow/trim_termini.py` applies an owner-digitised cut; it neither detects
  nor guesses a terminus, since that is an input decision (Section 2.2).

## D-013 — Cascade seal module, and a datum mismatch in the QGIS margin map
- Date: 2026-10-05
- Decided by: owner supplied the lake levels; the finding below is Claude Code's
  and needs the owner's call on which product to correct.
- Lake levels supplied 2026-10-05, WGS84 ellipsoidal: **LGP 117.3 m, GPL 27.6 m**
  (Alsek still outstanding, which blocks `competing_outlets`). Their difference,
  89.7 m, is consistent with the spec's "LGP sits ~100 m above GPL".

### 2018 result (k = 1.0, connectivity 8, no smoothing)
| quantity | value |
|---|---|
| h_pass | **233.262 m** |
| h_LGP | 117.300 m |
| Δ = h_pass − h_LGP | **+115.962 m — sealed** |
| thinning needed at the pass | **126.46 m** |
| pass cell | (4663, 3710) → map (908597.5, 1112827.5) |
| overtopping alternative | 253.1 m, so flotation opens first by ~20 m |

The pass sits **277 m** from the Phase 1 breach point for breach 1
(908752.5, 1112597.5). Those come from independent data — 2010 IFSAR with
`Bedrock_Aligned` in NAVD88 versus 2018 ArcticDEM with `bed_ellip` in ellipsoidal
heights — so the agreement is a genuine cross-validation of both.

### The finding: `margin_2018_lgp.tif` mixes vertical datums
Comparing the pipeline's output with the QGIS rasters over 2.48 M cells:

- **head surfaces agree**: median difference −0.003 m (IQR −0.070 to +0.061). The
  pipeline reproduces the QGIS head surface.
- **margins differ by a constant −6.703 m** (IQR −6.77 to −6.64).

Since the heads agree, the whole offset is in the lake level: the QGIS margin was
built with **h_LGP ≈ 110.60 m**, an orthometric value, against an **ellipsoidal**
head surface. The gap is the geoid height, ~6.8 m here.

**Consequence.** The QGIS margin map overstates how sealed the dam is by 6.7 m,
i.e. about 0.74 yr at −9.05 m/yr. The pipeline run is the internally consistent
one. Recommend regenerating `margin_2018_lgp.tif` from `data/derived/`.

### This also closes an open question under D-008
D-008 asked which datum the baseline's lake levels of 110.0 m and 17.0 m were
read in. The new ellipsoidal LGP level is 117.3 m; minus the ~6.8 m geoid that is
110.5 m, against the baseline's 110.0 m. **The baseline lake levels were NAVD88**,
consistent with the IFSAR surface they were read from — so the baseline is
internally consistent after all, and the "lake levels on the other datum" row of
the D-008 table (−0.83 yr) does not apply. The remaining D-008 item is the bed.

## D-014 — Drainage path, and why the downstream refinement is not used
- Date: 2026-10-06
- Decided by: Claude Code (implementation); **owner to confirm the two choices below**

### 2018 path (k = 1.0, connectivity 8, smoothing 0)
| quantity | value |
|---|---|
| total length LGP → GPL | 10.835 km |
| L_pass (pass → GPL) | **5.624 km** |
| bed sill | **58.71 m** at 1.471 km |
| ⇒ LGP can drain only to max(h_GPL, sill) | **58.71 m** — the sill controls, not GPL |
| thickness at pass | 336.26 m (min 98.8, mean 268.5) |
| flotation-zone length | 2.848 km |
| N at pass / mean over seal / min | +115.96 / +89.07 / −96.44 m |

**[SUPERSEDED by D-019 — this holds only under the along-path reading of the\nsill; under the bed-only reading GPL controls.]** LGP drains from 117.3 m
to 58.71 m, not to 27.6 m — about 58.6 m of drawdown rather than 89.7 m. That
feeds M4 directly and makes `V_w` materially smaller than a GPL-controlled
estimate would suggest.

### Choice 1 — representative N defaults to `at_pass`
The spec left this open. Default is N at the pass (+115.96 m); `--n-mode`
offers `mean_seal` (+89.07) and `min` (−96.44). The three differ enough to
matter for the conduit closure term, so this needs an owner decision before M6.

### Choice 2 — the spec's `refine_downstream` is implemented but does not run
The spec asks for steepest descent on head below the pass, because the minimax
path there is not unique. On the 2018 geometry that cannot work:

- a raw D8 descent from the pass **stalls after 63 m** in a local minimum;
- descending the flood-filled surface instead stalls too, on **plateaus** — nine
  consecutive cells share the value 232.714 m, and a strictly-downhill rule has
  nowhere to go;
- a fallback that follows the flood's parent links escapes plateaus but then
  stalls where the walk strays to the LGP side of the pass, whose parent chain
  runs back through forbidden upstream cells.

The underlying reason is that the head field is rough: below the pass it rises
at **288 of 899 steps**, by up to 5.15 m and 103.4 m in total. Strict descent is
not well posed on it.

`refine_downstream` therefore **verifies it reached GPL and returns the original
path unchanged if not**, logging a warning. That is not cosmetic: the truncated
route reported `L_pass` as 0.063 km instead of 5.624 km, an 89-fold error that
would have gone straight into the conduit model as its length.

So the reported path is the minimax route. It reaches GPL, descends monotonically
in spill by construction, and falls from 233.26 m to 27.60 m in head.

**For the owner:** the cleanest fix is smoothing. `smoothing_m` is still 0, which
the spec flagged as an open decision ("start at 0, then test a few ice
thicknesses"). Smoothing the surface over a few ice thicknesses should make the
head field well posed for descent and let the refinement run as specified.

### Note on the pass
Head varies by only 8 m (225.2–233.3) over the 3 km upstream of the pass, so the
pass is the crest of a broad near-plateau rather than a sharp saddle. ``h_pass``
is robust — it is a maximum — but the pass *cell* could move within that plateau
under small changes to the surface, so do not over-interpret its exact location.

## D-015 — Lake curves, and two things the volume estimate runs into
- Date: 2026-10-06
- Decided by: Claude Code (implementation); **two owner decisions needed**

### There is no LGP bathymetry in the repository
`OIB-and-Hig-depths-merged.gpkg` holds **23** points with a depth value. Sampled
against the 2018 classes they fall on **GPL (5)** and **ice (18)** — **none in
LGP**. Depths run 0–275 m with 275 the censored cap. So the plumb-line survey
covers Grand Plateau Lake, not the upper lake, and `interpolate_bathymetry` has
nothing to grid for LGP.

V_w is therefore reported as a **bracket**, not a number:

| bound | basin shape | V_w |
|---|---|---|
| upper | vertical walls (area constant with depth) | **200.0 Mm³** |
| lower | conic (deepest at centre) | **77.9 Mm³** |

A real basin lies between them. Drawdown is 117.3 → 58.71 m = 58.59 m, the bed
sill controlling (D-014).

**[D-019 adds a second V_w answer: under the bed-only sill it is 119–306 Mm³.]**\n\n**The spec's ~0.5 km³ expectation is superseded — RESOLVED 2026-10-06.** The
figure is reproducible, and the decomposition shows exactly what it assumed. Two
LGP outlines exist in the repo for the same lake (centroids 913 m apart; 98.8 % of
the 2018 outline lies inside the other):

| outline | area | drain to sill 58.7 m | drain to GPL 27.6 m |
|---|---|---|---|
| `lakes_2018.gpkg` class 2 | 3.414 km² | **200 Mm³** | 306 Mm³ |
| `outletStrengthCalculations/lakePolygon.gpkg` | 5.083 km² | 298 Mm³ | **456 Mm³** |

456 Mm³ ≈ 0.46 km³ is the ~0.5 km³ figure: it used the **larger outline**,
drainage **all the way to GPL** rather than to the bed sill, and **vertical
walls**. The bed sill is new information from D-014 that an earlier estimate
could not have had.

**Owner decision, 2026-10-06: use `lakes_2018`.** Reason given — it is the only
outline epoch-consistent with the rest of the inputs, which are all 2018: the
ArcticDEM ice surface (2018-09-03) and the lake outlines digitised from a 2018
Sentinel-2 scene (2018-09-05). `lakePolygon.gpkg` carries no date or provenance.
Consistency of epoch across layers outranks any single layer's extent, because
mixing a lake outline from one year with an ice surface from another silently
mis-states the dam geometry.

So **V_w = 78–200 Mm³ stands**, and the ~0.5 km³ figure should be retired from the
spec rather than reconciled.

`lakePolygon.gpkg` is **superseded**; it is referenced nowhere in the code or
configs, and `data_sources.md` now says so.

### GPL spills into an adjoining basin before it can absorb the flood
Its stage curve is not smooth. GPL's own basin holds **68 Mm³**, rising 1.50 m to
29.10 m. Between 29.1 and 29.6 m the connected area jumps **45.4 → 61.8 km²** and
absorbs a further **309 Mm³** for almost no extra rise: the lake tops a divide
into an adjoining depression. It does **not** connect to Alsek or LGP at these
stages.

Both V_w bounds (78 and 200 Mm³) exceed GPL's own 68 Mm³ capacity, so on this
geometry the flood spills laterally. **Owner decision 2:** whether the routing
stage should treat that adjoining basin as part of GPL, or as a separate
receiving body.

Reading a stage off the curve across that step is meaningless — it is a
discontinuity, not a slope. `detect_spill` finds it and the runner refuses to
interpolate through it. Detection is an outlier test rather than a fixed
threshold, because a smoothly widening basin gains a large *fraction* of its area
per step near its floor and a bare cutoff flags ordinary cones.

### A bug worth recording
The project's surface DEM has lakes masked to NaN. The first version of
`connected_stage_area` read depth from that DEM, so GPL's 45 km² — much the
largest term — contributed **zero** volume. Seed cells now take an effective
elevation of the current lake surface, and a test asserts a masked lake still
contributes its area times the rise.

## D-016 — Empirical inflow bracket
- Date: 2026-10-06
- Decided by: Claude Code (implementation); **one owner decision needed**
- All four of the spec's relations reproduce its expected values at V = 500 Mm³:
  2,780 / 4,824 / 16,941 m³/s, and Froehlich 50,641 — the last confirming that the
  spec's ~51,000 assumed **H_w = 100 m**, the nominal LGP–GPL head. Every curve
  integrates to V_w to machine precision (worst error 2×10⁻¹⁶).

### Result over the actual V_w bracket (78–200 Mm³, H_w = 58.6 m)
| relation | mode | Q_p (m³/s) | duration |
|---|---|---|---|
| Walder & Costa, tunnel | tunnel | 815 – 1,519 | 96 – 133 h |
| Clague & Mathews | tunnel | 1,388 – 2,611 | 57 – 77 h |
| Walder & Costa, non-tunnel | non-tunnel | 7,475 – 11,320 | 5.8 – 9.8 h |
| *Froehlich (reference only)* | *moraine breach* | *16,650 – 23,783* | *2.6 – 4.7 h* |

**The drainage mode matters far more than V_w.** Tunnel gives 815–2,611 m³/s;
non-tunnel gives 7,475–11,320 — a factor of **13.9**, against the ~2.6× spread
that the entire volume uncertainty contributes. Resolving whether this is a
tunnel or a non-tunnel failure is worth more than any improvement to the
bathymetry. **Owner decision.**

Which mode applies is physical, not statistical, so the code reports both and
prefers neither. Froehlich is flagged `reference_only` in the data structure, and
a test asserts the flag, so no table or plot can quietly promote a moraine-breach
regression into the applicable set.

### Plotting note
Tunnel and non-tunnel differ by ~14× in peak and ~20× in duration. A single
linear time axis squashes the fast case; a log axis distorts the curves into
shapes that no longer read as hydrographs. The figure therefore uses two panels
with a **shared y-axis** — peaks comparable, timescales honest, and equal area
under every curve.

## D-017 — ICESat-2 lake-level time series
- Date: 2026-10-07
- Decided by: Claude Code (implementation); **two owner decisions below**
- `gpbreach.io.nsidc` (CMR search + Earthdata download via `~/.netrc`) and
  `gpbreach.cascade.lake_levels`. 43 melt-season ATL06 granules, 2019–2026,
  cached in `data/raw/ICESat-2/ATL06/` (~0.7 GB, gitignored).

**ATL06, not ATL13.** ATL13 is the inland-water product and the obvious choice,
but its HydroLAKES reference mask does not contain these ice-dammed lakes, so it
returns nothing over all three regardless of track (D-008). ATL06 covers them and
agrees with ATL13 to 0.10–0.25 m where both exist.

### Result: 33 usable of 129 lake-overpass pairs
| lake | epochs | median | range | typical scatter |
|---|---|---|---|---|
| Alsek | 18 | 32.21 m | 27.88–32.93 (**5.05 m**) | 0.21 m |
| GPL | 14 | 30.35 m | 29.77–31.00 (1.23 m) | 0.11 m |
| LGP | **1** | 119.75 m | — | 0.07 m |

All WGS84 ellipsoidal (ITRF2014); `level_egm2008_m` carries the geoid-corrected
form. Rejected rows are retained with a `flag` saying why, so what was discarded
is visible rather than merely absent.

### 1. Alsek has a strong melt-season signal — and fills the `alsek: null` gap
27.88 m (1 May 2019) → 32.56 m (30 Jul 2019); 29.72 m (26 May 2022) → 32.16 m
(25 Jun 2022). It rises ~4–5 m through the season, far more than GPL's 1.2 m
total spread. **Owner decision: what value to set for `lake_levels.alsek`** — a
melt-season median of ~32.2 m, or a date-matched value, given the variation is
larger than the quantity itself in some comparisons.

### 2. Both lakes read ~2.5–2.75 m above the configured levels
| lake | config | ICESat-2 median | difference |
|---|---|---|---|
| GPL | 27.6 m | 30.35 m | **+2.75 m** |
| LGP | 117.3 m | 119.75 m | **+2.45 m** |

Two independent lakes agreeing to within 0.3 m points at a **systematic offset**,
not coincidence. The leading candidate is the NAD83 ↔ ITRF2014 ellipsoid
difference already flagged under D-008, which is order a metre or two in Alaska;
ICESat-2 is ITRF2014 and the provenance of the configured values is not recorded.
**Owner decision: resolve before these levels feed Δ**, since lake level enters
the margin one-for-one (D-013) and 2.6 m is ~0.3 yr of thinning.

### 3. LGP is genuinely under-sampled, and buffering is not the cause
Only **6 of 43** granules place any segment in LGP's bounding box, and just 2
survive quality screening. At 3.414 km² it is a small target that the tracks
mostly miss. An adaptive buffer was added — a fixed −150 m inward buffer costs
Alsek a tenth of its area but LGP a half, which is not scale-free — and it did not
change the count. **LGP's level cannot be monitored this way**; ATL03 photons or
a different mission would be needed.

### Bug fixed
`atl06.granule_info` assumed beam `gt1l` exists. `ATL06_20220825050258` carries
only the gt2 and gt3 pairs, and the assumption aborted a 43-granule run at
granule 23. `granule_info` now searches the beams for one that has the field, the
runner catches per-granule failures so one bad file cannot cost the other
forty-two, and a regression test covers it.

## D-018 — Component swapping framework, and the Millan/IceBoost bed disagreement
- Date: 2026-10-08
- Decided by: Claude Code (implementation); **bed choice is an open owner question**
- `configs/cascade_data.yaml` now declares **named** beds, surfaces and thinning
  rates; `cascade.scenarios` runs any combination and `run_scenarios` compares
  them. See `docs/component_swapping.md`.

Solves are cached on `(bed, surface, k, connectivity)` — deliberately excluding
the thinning rate, which cannot move `h_pass`. A thinning sweep therefore costs
one solve rather than one per rate.

### Millan vs IceBoost on the 2018 surface
| | h_pass | Δ | thinning needed | onset | bed sill | drawdown | V_w |
|---|---|---|---|---|---|---|---|
| Millan | 233.26 | +115.96 | 126.5 m | 2032.0 | 58.71 | 58.59 | 78–200 Mm³ |
| IceBoost | 245.00 | +127.70 | 139.3 m | 2033.4 | **185.39** | **0.00** | **0** |

**Timing is insensitive (+1.4 yr); flood size is not.** IceBoost's sill sits above
LGP's own surface, which taken literally means no drainage at all. Along the 2018
path IceBoost holds the bed above 58 m for 5.8 km against Millan's 21 m — mean
difference +105 m, sd 62 m. This confirms the owner's observation that the sill is
a far broader feature in IceBoost.

A V_w of zero is not credible for a lake that drains, so this is a reason to
**doubt the IceBoost ingest or its datum**, not yet a reason to prefer either bed.
`vertical: ellipsoid_wgs84` for IceBoost is marked TO CONFIRM.

### The arbiter exists but is not yet usable
`glathida_thickness/` holds 1,947 IceBridge UAF radar soundings over both
glaciers. GlaThiDa's `elevation` field does **not** match our ice surface (median
+111 m, IQR +95 to +236), so `bed = elevation − thickness` is not defensible as
ground truth yet. An earlier comparison built on that assumption put both beds
~200 m below radar; that number is **withdrawn** — it reflects the unresolved
field definition. Settling what `elevation` means in IRUAFHF2 (Tober et al. 2025)
is the highest-value open item here.

## D-019 — Two sills, and a correction to D-014, D-015 and D-018
- Date: 2026-10-08
- Decided by: owner spotted the discrepancy; **which sill governs is still an
  owner decision**

The owner read the IceBoost bed in QGIS, saw a pass near 20 m, and asked why the
analysis reported 185 m. The analysis was wrong — or rather, it was answering a
different question than the one it appeared to answer.

### The two quantities
| | millan | iceboost |
|---|---|---|
| **along-path**: highest bed on the flotation route | 58.71 m | 185.39 m |
| **bed-only**: lowest bed ridge on *any* route | **11.08 m** | **33.85 m** |

The flotation route is chosen to minimise hydraulic potential, not bed elevation,
so it will climb over a high bed ridge wherever the ice is thick. Reporting its
maximum bed as "the sill" implied it was the lowest bed route. It is not, and the
two differ by 48 m on Millan and 152 m on IceBoost.

### What that corrects
- **D-018's "IceBoost gives V_w = 0, no flood" is withdrawn.** Under the bed-only
  reading IceBoost's sill is 33.85 m, an 83.45 m drawdown and V_w of 111–285 Mm³.
  The suspicion recorded there that V_w = 0 was implausible was right; the cause
  was this definition, not the data or its datum.
- **D-015's V_w of 78–200 Mm³ is one of two answers, not the answer.** Millan's
  bed-only pass of 11.08 m lies *below* GPL at 27.6 m, so the receiving lake would
  control rather than any sill: drawdown 89.7 m and V_w 119–306 Mm³.
- **D-014's "the bed sill controls the drawdown, not GPL" holds only under the
  along-path reading.** Under the bed-only reading it is false for Millan.

### Both are now reported
`scenarios.bed_pass` computes the bed-only minimax region-to-region, and
`run_scenarios` prints both readings with V_w and Q_p under each:

```
bed        sill def      sill  h_final  drawdn        controls    V_w Mm3   Qp tunnel
millan     along-path   58.71    58.71   58.59        bed_sill     78-200    815-1519
millan     bed-only     11.08    27.60   89.70  receiving_lake    119-306   1079-2012
iceboost   along-path  185.39   185.39    0.00        bed_sill      0-0         0-0
iceboost   bed-only     33.85    33.85   83.45        bed_sill    111-285   1029-1918
```

### The physics question, unresolved
*Along-path* assumes the flood drains through the conduit that opened and that the
route cannot migrate — conservative, smaller flood. *Bed-only* assumes water finds
the lowest available route once drainage is under way — larger flood. The spec
specifies the first. The truth is likely between, and the gap is large enough
(Millan 78–200 against 119–306 Mm³) that it should be chosen deliberately rather
than inherited from a definition.

### Also fixed
The scenario cache key now carries a schema version. Without it a cache written
before `bed_pass_m` existed would still have been served, answering the new
question with the field missing. A test asserts the versioned key differs from the
unversioned one.
\n

## D-020 — IceBoost datum, the spill basin, and the armor threshold
- Date: 2026-10-09

### IceBoost is ellipsoidal — CONFIRMED by owner
`configs/cascade_data.yaml` updated from TO CONFIRM. This closes gap A4 in
`docs/path_to_failure_probability.md`: the +105 m mean offset from Millan along
the path is a genuine model disagreement, not a datum artifact.

### The spill basin drains — it does not pond
The D-015 "adjoining basin that absorbs 309 Mm³" was an **artifact of the DEM**.

`surf_2m` is clipped to ice with lakes masked, so the coastal lowland around GPL
is full of holes. A flood fill on it finds an apparently enclosed depression
reaching −46 m that never spills, even at a 60 m stage over a 30 × 28 km window.
A closed depression 53 m below sea level beside the ocean is not credible.

Repeating the test on the **IFSAR DSM**, which is not ice-clipped, reverses it:
at GPL's *current* level (27.6 m ellipsoidal = 20.8 m NAVD88) the lake is already
connected to a 281 km² region that reaches sea level and every tile edge. That is
simply GPL's existing drainage to the Alsek River and Dry Bay.

**Method, for repeating this anywhere:**
1. Use a DEM that is **not** clipped to ice. Ice-clipped products cannot answer
   drainage questions off the glacier.
2. **Mask the void sentinel.** The IFSAR DSM carries −10000 with no nodata flag
   set. Read naively it is "valid 100%", and those voids connect everything —
   the first run reported 555 km² draining to the sea for the wrong reason.
3. Coarsen by **block minimum**, not averaging. Averaging closes the narrow
   channels that decide connectivity.
4. Flood fill from the lake and ask whether the region reaches a known outlet or
   the sea.

**Consequence.** GPL is not a two-reservoir system with a large hidden store. The
"spill" at ~29.6 m is GPL overtopping onto the coastal plain — and the outlet
polygon lies **inside** that footprint, so the lateral spill is plausibly the
outflow route itself rather than a separate sink. Routing should treat it as
conveyance, not storage. Whether the flood then **bypasses** the boulder outlet,
which would leave the armor unloaded, is now the open question.

### Armor threshold from D50 = 0.961 m (owner, 2026-10-09)
`tau_c = theta_c (rho_s − rho_w) g D50`, with rho_s = 2650 kg m⁻³:

| theta_c | tau_c | depth-slope product needed |
|---|---|---|
| 0.030 | 467 Pa | 0.0476 m |
| 0.045 | 700 Pa | 0.0714 m |
| 0.060 | 933 Pa | 0.0951 m |

Flow depth required to mobilise, by slope:

| slope | theta_c 0.03 | 0.045 | 0.06 | Lamb et al. (2008) |
|---|---|---|---|---|
| 0.005 | 9.51 m | 14.27 m | 19.03 m | — |
| 0.010 | 4.76 m | 7.14 m | 9.51 m | — |
| 0.020 | 2.38 m | 3.57 m | 4.76 m | 4.47 m |
| 0.050 | 0.95 m | 1.43 m | 1.90 m | 2.25 m |
| 0.100 | 0.48 m | 0.71 m | 0.95 m | 1.34 m |

**The outlet slope now matters more than anything else in this calculation** — it
swings the required depth by a factor of 20 across the plausible range, far more
than the choice of `theta_c`. Measuring it is cheap and is the next thing to do.
Lamb et al.'s slope-dependent `theta_c ~ 0.15 S^0.25` is also tabulated, because
constant-Shields is optimistic on steep boulder beds; both should be reported.

