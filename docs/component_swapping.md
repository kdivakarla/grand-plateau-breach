# Swapping components

How to run the cascade with a different bed, surface or thinning rate, and what
each swap actually changes.

## The one-liner

```bash
python -m gpbreach.cascade.run_scenarios --list                      # what's available
python -m gpbreach.cascade.run_scenarios --beds millan iceboost      # compare two beds
python -m gpbreach.cascade.run_scenarios --beds millan --k 0.9 0.95 1.0
python -m gpbreach.cascade.run_scenarios --beds millan iceboost --thinning all
```

Output: a comparison table on stdout and `data/derived/scenario_comparison.csv`.

## Adding a component

Everything is declared in `configs/cascade_data.yaml`; no code changes.

```yaml
beds:
  - id: my_new_bed
    path: data/standard/bed_my_new_bed.tif      # where build_grid will put it
    source_path: data/processed/whatever.tif     # where it arrives, any CRS/resolution
    vertical: ellipsoid_wgs84
    source: one line on provenance
```

Then `python -m gpbreach.cascade.build_grid` regrids it onto the canonical 5 m
ESRI:102247 grid. Surfaces and thinning entries work the same way. A thinning
entry is just an id and a rate.

## What each swap costs, and why

Only the saddle search is expensive (~35 s per distinct geometry here). What it
depends on determines the cost:

| swap | moves `h_pass`? | moves `V_w`? | cost |
|---|---|---|---|
| **bed** | yes | **yes**, via the bed sill | full re-solve |
| **surface** | yes | no — the sill is bed-only | full re-solve |
| **flotation `k`** | yes | no | full re-solve, 1-D in `k` |
| **thinning `dh/dt`** | **no** | no | **free** |

Thinning is free because it only sets how fast the threshold rises. "How much
thinning is needed" is a geometry question answered once; dividing by a rate is
arithmetic. Solves are cached on `(bed, surface, k, connectivity)` in
`data/interim/scenario_cache/`, so adding thinning rates to a matrix costs nothing.

## How a swap propagates

Two branches, and **a bed change moves both**:

```
bed, surface, k ──► h_pass ──► Δ = h_pass − h_LGP
                           ──► thinning needed = Δ / κ
                   + dh/dt ──► years to breach ──► onset year

bed             ──► bed sill ──► h_final = max(h_GPL, sill)
                             ──► drawdown ──► V_w ──► Q_p
```

A surface or thinning swap moves **only the timing**. A bed swap moves the timing
*and* the flood size — which is why the bed deserves the most scrutiny.

## Worked example: Millan vs IceBoost

```
TIMING
bed        h_pass    Δ        thinning needed   onset
millan     233.26   +115.96        126.5        2032.0
iceboost   245.00   +127.70        139.3        2033.4

FLOOD SIZE
bed        sill     h_final  drawdown  V_w (Mm³)   Q_p tunnel
millan     58.71     58.71     58.59    78–200      815–1519
iceboost  185.39    185.39      0.00       0–0           0–0
```

**Timing barely cares: +1.4 yr. Flood size depends on which sill you mean.**

Those `V_w` figures use the **along-path** sill. Under the **bed-only** sill —
the lowest bed ridge on any route, which is what you see reading the raster in
QGIS — Millan gives 119–306 Mm³ and IceBoost 111–285 Mm³, and neither blocks
drainage. `run_scenarios` prints both; see D-019.

IceBoost puts the *along-path* sill at 185.4 m, above LGP's surface, which taken
literally says no drainage on that route. Its *bed-only* sill is 33.85 m, which
permits an 83 m drawdown. The two models disagree about the drainage corridor far
more than they disagree about the dam — but so do the two definitions of the sill,
and that gap is the larger one.

The owner's observation is confirmed and then some: along the 2018 path, IceBoost
holds the bed above 58 m for **5.8 km**, where Millan does so for **21 m**. Mean
difference +105 m, sd 62 m.

**This is not yet a reason to prefer either bed.** IceBoost's vertical reference
is still marked TO CONFIRM, and the +105 m mean offset from Millan along the path
is large enough to want independent checking before either is trusted.

## The validation that is still open

`glathida_thickness/` holds 1,947 IceBridge UAF radar soundings covering both
Grand Plateau glaciers — in principle the arbiter. It is not usable yet: GlaThiDa's
`elevation` field does not match our ice surface (median +111 m, IQR +95 to +236),
so `bed = elevation − thickness` is not currently a defensible ground truth. An
earlier comparison on that assumption suggested both beds sit ~200 m below radar;
**disregard it** — it reflects the unresolved field definition, not the beds.

Resolving what `elevation` means in the IRUAFHF2 product (Tober et al. 2025,
doi:10.5067/225AMW1MQB4O) would turn 1,947 ground-truth points into a direct test
of which bed to trust — the single highest-value open item for this question.
