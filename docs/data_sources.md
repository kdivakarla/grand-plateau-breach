# Data sources

Provenance, datum, and units for every dataset. Cells marked `?` are for the
owner to fill in or confirm; Claude Code does not assume a datum
(CLAUDE.md Section 5).

## Currently used by the baseline

| Dataset | Use | Location | Format | Horizontal CRS | Vertical datum | Units | Notes |
|---|---|---|---|---|---|---|---|
| `2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif` | base-year glacier surface | repo parent | GeoTIFF | ESRI:102247 | **?** (D-008) | m | 2496x2543 @ 5 m; 6,347,328 valid cells; range 0–621 m |
| `Bedrock_Aligned.tif` | bed | repo parent | GeoTIFF | ESRI:102247 | **?** (D-008) | m | same grid; 2,773,050 valid cells; range -998–612 m. **Provenance undocumented — how was it derived, and "aligned" to what?** |
| `Buoyant_Maps_LGP_GP/Buoyancy_State_<2027..2032>.tif` | legacy output, regression target | repo parent | GeoTIFF | ESRI:102247 | n/a (margin) | m | negative = buoyant |
| `Buoyant_Maps_Alsek_GP/Buoyancy_State_<2032..2036>.tif` | legacy output, regression target | repo parent | GeoTIFF | ESRI:102247 | n/a (margin) | m | negative = buoyant |

The surface and bed grids are byte-for-byte co-registered: identical size,
origin (900605, 1122935), and 5 m pixels. `preprocess/build_dam.py` asserts this
and refuses to resample.

## Present in the repository, not yet used

| Dataset | Use | Location | Format | Horizontal CRS | Vertical datum | Units | Notes |
|---|---|---|---|---|---|---|---|
| `Lidar_Elevation_Loso/GrandPlateau_2p5m/*.tif` | surface history 2005–2020 | repo parent | GeoTIFF | ? | **?** | m | 9 epochs (2005, 2009, 2011x2, 2014x2, 2017, 2020x2). **Would constrain D-004 directly** |
| `Lidar_Elevation_Loso/DiffStats_GP/*.json` | dh/dt from cop30 2005–2020 | repo parent | JSON | — | — | m | north/south Grand Plateau |
| `OIB-and-Hig-depths-merged.gpkg` | surface/bed comparison | repo parent | GPKG | EPSG:32607 | **?** | m | 934 points; `ice_thickness_m`, `bed_h_m`, `surface_height_m` |
| `FULL-transects-with-avg-OIB-values.csv` | dh/dt source | repo parent | CSV | — | **?** | m | the transects behind the -9.05 m/yr figure |
| `outletStrengthCalculations/lakePolygon.gpkg` | lake geometry | repo parent | GPKG | EPSG:32607 | n/a | — | 1 polygon; would let D-005 seeds be defined from real lake outlines |
| `outletStrengthCalculations/outletOutline.gpkg` | outlet geometry | repo parent | GPKG | ESRI:102247 | n/a | — | 1 polygon, extent lies **outside** the raster domain |
| `first-attempt-clipped.tif` | ? | repo parent | GeoTIFF | ? | ? | m | unclear role |

## Later phases (not needed for Phase 1)

| Dataset | Use | Notes |
|---|---|---|
| Millan et al. (2022) ice thickness | bed prior | not in repo |
| NASA IceBridge UAF radar/lidar | bed correction, surface history | not in repo as raw HDF5 |
| ArcticDEM strips/mosaic | surface history | not in repo |
| Plumb-line bathymetry (Higman) | lake depth | `>275 m` values censored |
| ERA5 | daily forcing | phase 4 |
| PyGEM outputs | mass-balance ensemble | phase 5; owner emailed Loso 2026-08-25 |
