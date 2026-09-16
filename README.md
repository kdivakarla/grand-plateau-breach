# gpbreach — Grand Plateau breach forecasting

> **Preliminary research code — not peer-reviewed. Do not use for operational
> hazard decisions.**
>
> The dates below reproduce an earlier deterministic calculation; they are not a
> validated forecast. Open items, largest first:
>
> - the **thinning rate** is a single uniform linear trend, and it is the input
>   the answer is most sensitive to by a wide margin (≈2.2 yr per m/yr); an
>   updated estimate is pending;
> - **no uncertainty is attached** to the headline dates — they are one
>   realization, not a distribution;
> - the **vertical datum** of the surface and bed grids is not declared in either
>   file. The surface is NAVD88 per the USGS Alaska IFSAR spec; the bed is
>   aligned to OIB radar picks and so is probably WGS84 ellipsoid. The geoid
>   height here is a near-constant +6.80 m, which is worth **0.08 yr** in the
>   likely configuration and at most 0.83 yr in the worst one (D-008);
> - the **bed** is treated as exactly known; a *uniform* bias barely matters
>   (0.011 yr per m), but spatially correlated error could move the saddle.
>
> Every such item is tracked in [`docs/decisions.md`](docs/decisions.md). Read it
> before citing or reusing anything here.

Forecasts the timing of two glacier-dammed lake breaches at Grand Plateau
Glacier, southeast Alaska.

| | breach | baseline date |
|---|---|---|
| **Breach 1** | upper (LGP) lake → Grand Plateau Lake | **2030** (crossing 2029.64) |
| **Breach 2** | Alsek Lake → Grand Plateau Lake | **2035** (crossing 2034.60) |

Both reproduce the owner's original QGIS results exactly, and the regression
suite checks the port against all eleven `Buoyancy_State_<year>.tif` rasters
cell-by-cell.

## The idea in one equation

The breach criterion is hydrostatic flotation plus a connected pathway across the
dam. Because the criterion is linear in the surface and bed grids, it collapses
to a **static field** and a **rising threshold**:

```
Phi_f(x) = f*S(x) + (1-f)*B(x)          buoyant where  Phi_f <= L - f*(dh/dt)*t
```

Two basins connect the moment the threshold reaches the **saddle** (lowest pass)
between them, so the breach date is closed-form:

```
year = t0 + (Phi_crit(f) - L) / (-f * dh/dt)
```

`Phi_crit` depends only on `f` and the grids — not on lake level or thinning
rate. That is what makes the Monte Carlo cheap: sampling `L` and `dh/dt` is pure
arithmetic, sampling `f` needs only a precomputed 1-D curve, and only a perturbed
*bed* forces a new saddle solve. Full derivation in `docs/baseline_method.md`.

## Quick start

```bash
export PYTHONPATH=src            # until the env below exists and `pip install -e .` is run

python -m gpbreach.cli preprocess  configs/baseline_breach1_lgp.yaml
python -m gpbreach.cli run         configs/baseline_breach1_lgp.yaml
python -m gpbreach.cli sensitivity configs/baseline_breach1_lgp.yaml
python -m gpbreach.cli fcurve      configs/baseline_breach1_lgp.yaml
```

or `make -f workflow/Makefile all`.

## Environment

`environment.yml` is **proposed, not created** — making or changing an
environment needs owner approval (CLAUDE.md §2.2). Everything above runs today in
the existing `gpgn-318` env, with two fallbacks:

- no `rasterio` → preprocessing shells out to `gdal_translate` (found inside
  QGIS.app). One warning per raster; results are identical.
- no `pytest` → run `python workflow/run_tests_without_pytest.py`, which executes
  the real test functions with a small shim.

## How much does the answer move?

`gpbreach sensitivity`, breach 1 (breach 2 is very similar):

| parameter | ± step | span (yr) | per unit |
|---|---|---|---|
| `dh/dt` | 1 m/yr | **4.4** | 2.20 yr per m/yr |
| lake level | 10 m | 2.4 | -0.12 yr per m |
| uniform surface bias | 10 m | 2.2 | 0.11 yr per m |
| flotation fraction `f` | 0.02 | 0.8 | 0.21 yr per 0.01 |
| uniform bed bias | 20 m | 0.4 | 0.011 yr per m |

The thinning rate dominates. `f` — the parameter the method is named after — is
comparatively weak: the whole range `f` ∈ [0.82, 0.96] spans under three years.
A *uniform* bed bias is nearly irrelevant, but that says nothing about
spatially correlated bed error, which can move the saddle to a different pass.

## Layout

```
configs/      run definitions (protected)
docs/         decisions.md, data_sources.md, baseline_method.md
src/gpbreach/
  saddle.py     minimum-bottleneck threshold between two cells
  dam.py        DamObject: two grids + two seeds
  components/   surface, bed, lake, criterion (swappable physics, protected)
  engine.py     vectorised over realizations; deterministic = N of 1
  analysis/     fcurve, sensitivity, report
tests/        unit, regression (expected/ is protected)
legacy/       why there is no legacy script, and how the method was recovered
```

## Status

Phase 1 complete and passing. Ten decisions are logged in `docs/decisions.md`;
D-005 (which basins each breach connects) and D-009 (Alaska Albers as the
analysis CRS) are confirmed, the rest remain `TO CONFIRM`. **D-008 — the vertical
datum — is the open item most likely to change the answer**, and it becomes
urgent as soon as OIB lidar (WGS84 ellipsoid heights, UTM 7N) is mixed with the
IFSAR surface.

Phases 2–6 (hindcast harness, uncertainty one source at a time, ERA5 daily
forcing, PyGEM ensembles, Sobol analysis) are not started.

No license file is present, so the default is all-rights-reserved. Add one if you
want others to be able to reuse this.
