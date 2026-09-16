# Baseline method

The owner's deterministic method, as ported. Written by Claude Code from (a) the
owner's verbal description and (b) the legacy QGIS rasters themselves.
**Awaiting owner confirmation** (CLAUDE.md Section 8.3, Step 0).

## 1. Owner's description

> I used `Bedrock_Aligned` as the bedrock, `2010_IFSAR` as the base DEM, then
> calculated the thinning rates from transects with OIB to get about -9 m/year.
> I tried to find the time when the lake level is 90% of the ice dam thickness
> and set this as the threshold. Once there was a full pathway below this
> threshold, I could say the breach will occur. The way I did this was separating
> ice dam thickness into above lake level and below lake level. Then the above
> bit thins by 9 m per year, and once that is less than 10 percent of the ice
> below lake level, that meets our threshold.

## 2. What the legacy rasters actually contain

There was no legacy script to read, so the method was recovered from the output
rasters. Every `Buoyancy_State_<year>.tif` in `Buoyant_Maps_LGP_GP/` and
`Buoyant_Maps_Alsek_GP/` is an exact linear function of the two input grids.
Least squares over 300,000 sampled cells returns **zero residual** (max
|residual| < 1e-4 m, i.e. float32 round-off):

```
Buoyancy_State_Y = S + 0.1 * B + c_Y
```

with `S` = `2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif`, `B` = `Bedrock_Aligned.tif`,
and `c_Y` decreasing by exactly **9.05** per year. Writing
`c_Y = -(1.1*L + 9.05*(Y - 2010))` and solving for `L` gives the same value from
every year of both series, to 4 decimal places:

| series | folder | implied L | years present |
|---|---|---|---|
| Breach 1 (upper/LGP lake) | `Buoyant_Maps_LGP_GP` | **110.0 m** | 2027–2032 |
| Breach 2 (Alsek Lake) | `Buoyant_Maps_Alsek_GP` | **17.0 m** | 2032–2036 |

Both recover as exact round numbers, and their difference (93 m) matches the
"upper lake sits about 100 m above Grand Plateau Lake" description. So the
legacy parameters are fixed with no freedom left:

| quantity | value | note |
|---|---|---|
| base DEM year | 2010 | |
| surface change `dh/dt` | **-9.05 m/yr** | not -9.0; see D-004 |
| above/below thickness ratio `r` | **0.1** | equivalently `f = 0.909091`; see D-003 |
| lake level, breach 1 | 110.0 m | D-007 |
| lake level, breach 2 | 17.0 m | D-007 |
| grid | 2496 x 2543 at 5 m, ESRI:102247 | all inputs co-registered exactly |

The raster is **negative where the ice is buoyant**.

## 3. The criterion, and why it is flotation

Split the dam at lake level: `H_above(t) = S + dhdt*t - L`, `H_below = L - B`.
The cell is buoyant when `H_above <= r * H_below`. Expanding and dividing by
`(1 + r)`, with `f = 1/(1+r)`:

```
f * (S + dhdt*t) + (1 - f) * B  <=  L                                    (1)
```

which is exactly `h_lake >= f * H_ice` with `h_lake = L - B` and
`H_ice = S + dhdt*t - B`. The owner's ratio `r` and a flotation fraction `f` are
the same parameter:

```
f = 1 / (1 + r)        r = (1 - f) / f
r = 0.1  <->  f = 0.909091
```

**Note for D-003:** the owner's stated intent was `f = 0.90`, which corresponds
to `r = 1/9 = 0.1111`, not `r = 0.1`. The rasters use `r = 0.1`, i.e.
`f = 0.9091`. The discrepancy is small — using `f = 0.90` exactly moves both
breach dates about 0.4 yr earlier and does not change either integer year — but
it has been left as the owner built it, not "fixed" (CLAUDE.md Section 8.1).

## 4. The breach event is a connectivity test

A single buoyant cell is not a breach. The owner's rule is that a breach occurs
in the first year a **connected pathway** of buoyant cells spans the dam between
two lake basins. This is the one part of the method that is not in the rasters:
the owner read it off the maps by eye.

Restating (1) on the base-year grids, define the **static flotation field**

```
Phi_f(x) = f * S(x) + (1 - f) * B(x)                                     (2)
```

Then (1) is `Phi_f(x) <= L - f*dhdt*t`: the buoyant region is a sub-level set of
a *time-independent* field under a threshold rising linearly at `f*|dhdt|` per
year. Two basins first connect exactly when the threshold reaches the **saddle
value** `Phi_crit` between them — the lowest pass on any path from one to the
other. So the breach date is closed-form:

```
t_breach = (Phi_crit(f) - L) / (-f * dhdt)                               (3)
year     = 2010 + t_breach
```

and the annual map the owner would have read is the first integer year at or
after that crossing, `ceil(year)`.

This replaces "render a raster per year and look" with one scalar per breach
site, and it is what makes uncertainty analysis affordable: **`Phi_crit` depends
only on `f` and the two grids** — not on lake level, not on thinning rate. Those
enter only through (3), so sampling them is free. `Phi_crit(f)` is a
one-dimensional curve that can be tabulated once (`gpbreach fcurve`). Only a
perturbed *bed* forces the saddle to be recomputed.

## 5. Basins and seeds

Three basins of `Phi_f` dominate the domain. Each is represented by its
minimum-`Phi` cell, a stable interior point that does not move when `f`, lake
level, or thinning rate change:

| basin | seed (row, col) | map x, y (ESRI:102247) | interpretation — **TO CONFIRM (D-005)** |
|---|---|---|---|
| A | (1903, 649) | 903852.5, 1113417.5 | Grand Plateau Lake (central) |
| B | (122, 22) | 900717.5, 1122322.5 | Alsek Lake (NW) |
| C | (2486, 2450) | 912857.5, 1110502.5 | upper / LGP lake (SE) |

Saddle values are properties of the field alone and are independent of lake level:

| pair | `Phi_crit` | with L | breach year |
|---|---|---|---|
| A–C (breach 1) | 271.5518 m | 110.0 m | 2029.636 → **2030** |
| A–B (breach 2) | 219.4077 m | 17.0 m | 2034.602 → **2035** |

Both reproduce the owner's published dates.

## 6. Verification

`tests/regression/test_baseline.py` asserts, for both breaches:

1. the integer breach year (2030 / 2035) exactly;
2. `Phi_crit` and the continuous crossing to 1e-6 m / 1e-3 yr;
3. the located breach cell;
4. that 4- and 8-connectivity give the identical answer;
5. **the ported `legacy_margin` against every one of the owner's eleven
   `Buoyancy_State_<year>.tif` rasters, cell for cell** (max |difference| below
   1e-3 m over ~2.77 million cells per year);
6. that applying the connectivity test to the owner's *own* rasters returns the
   same first-connected year.

## 7. Known departures from the brief

- CLAUDE.md Section 3 says "collapse space early" to 1D dam objects. That is not
  possible here: the criterion is a 2D connectivity test, so the `DamObject`
  carries two grids plus seeds. Space is instead collapsed one step later, at
  `Phi_crit`, and the engine never opens a raster.
- CLAUDE.md Section 5 asks for margins in m w.e. The baseline `f` is prescribed,
  not derived from densities, so the margin is reported in metres as
  `f*H - h_lake`, which is algebraically the same comparison. See D-002.
