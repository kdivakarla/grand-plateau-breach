# legacy/

Read-only reference. Excluded from the package.

**There is no legacy script.** The 2030/2035 results were produced with the QGIS
raster calculator / Python console and the code was not saved. What survives is
the output: eleven `Buoyancy_State_<year>.tif` rasters in the parent directory.

The method was therefore recovered *from those rasters* rather than from source.
The recovery, its numerical evidence, and the resulting parameters are written up
in `docs/baseline_method.md` §2. In short: each raster is an exact linear
function of the two input grids with zero residual, which pins every parameter
(`r = 0.1`, `dh/dt = -9.05 m/yr`, base year 2010, `L = 110 m` and `L = 17 m`)
with no freedom remaining.

`tests/regression/test_baseline.py::test_matches_legacy_rasters_every_year`
compares the port against all eleven rasters cell-by-cell and serves as the
characterization test that CLAUDE.md Section 8.3 Step 2 asks for.

If the original QGIS expressions or console history do turn up, put them here —
they would settle D-005 (which basins the owner actually connected), which is the
one part of the method the rasters do not record.
