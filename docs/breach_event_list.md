# Breach event list — input for gpbreach Phase 2 (hindcast + calibration)

One block per lake-breach event you want in the analysis. This becomes the
`events/` registry described in the Phase 2 design, so the fields here map
one-to-one onto that schema.

## How to fill this in

- **List every event you want, even with no data.** A row with nothing but a date
  and a citation is still useful — it shapes the schema, and `data_status: none`
  is a valid state. Don't pre-filter to the events that happen to be ready.
- **Leave blanks blank.** An empty field is information ("not known yet"); a
  guessed number is not. Anything you're unsure of, mark `?`.
- **Rough is fine.** "mid-May, some year around 2016" beats omitting the event.
- Add or drop fields freely — if a field doesn't fit an event, say so in `notes`.

### Field reference

| field | meaning |
|---|---|
| `id` | short slug, lowercase, no spaces — becomes the directory name |
| `status` | `observed` (historic, we know the answer) or `forecast` |
| `date` / `precision` | when it breached; precision = `day` \| `season` \| `year` |
| `source` | paper, report, or DOI. `?` if you know the event but not the reference |
| `data_status` | `complete` \| `partial` \| `none` — your call, gates whether it produces output |
| `surface` | is there a DEM, which one, what epoch, **how far from the breach date** |
| `bed` | `none` \| `millan` \| `millan+radar` \| `local` — and which flat-bed variants are plausible |
| `lake_level` | elevation if known, **and its vertical datum**; or "needs estimating" |
| `dhdt` | thinning rate if known, with source — often unnecessary, see below |
| `notes` | spillway outlets, progressive vs. sudden drainage, repeat events, anything odd |

**On `surface`:** this is the field that most affects effort. If the DEM is close
in time to the breach, dh/dt nearly drops out of the calibration — so "2016 DEM,
event in 2016" is worth far more than a dense surface history from a decade away.

**On datums:** every elevation needs its vertical datum recorded. Getting this
wrong silently shifted our own Grand Plateau answer (see `docs/decisions.md`
D-008 in the repo) — it's the single easiest way to corrupt a multi-event
comparison.

---

## Summary table

Fill this in as you go; it's the at-a-glance view.

| id | name | status | date | data_status |
|---|---|---|---|---|
| `alsek-grand-plateau` | Alsek Lake → Grand Plateau Lake | forecast | 2035 | complete |
| `lgp-grand-plateau` | Upper (LGP) lake → Grand Plateau Lake | forecast | 2030 | complete |
|  |  |  |  |  |
|  |  |  |  |  |
|  |  |  |  |  |

---

## Events

### alsek-grand-plateau  *(worked example — already complete, nothing to do here)*

```
id:          alsek-grand-plateau
name:        Alsek Lake drains into Grand Plateau Lake
status:      forecast
location:    Grand Plateau Glacier, Glacier Bay NP&P, SE Alaska
date:        2035            precision: year   (continuous crossing 2034.60)
source:      this project; AGU abstract in supplemental_info/
data_status: complete

surface:     2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif, epoch 2010, 5 m
             datum NAVD88 (per USGS Alaska IFSAR spec; not declared in the file)
             offset from breach date: ~25 yr  -> dh/dt matters a lot here
bed:         local — Bedrock_Aligned.tif, aligned to OIB radar picks
             (median +0.04 m vs. OIB bed_h_m at 932 points)
             datum: inherits OIB, probably WGS84 ellipsoid — D-008, unconfirmed
lake_level:  17.0 m   datum: ? (presumed same as IFSAR surface — needs confirming)
dhdt:        -9.05 m/yr, uniform linear, from OIB/DEM transects. Revision pending.
notes:       2-D connectivity criterion, not a flat bed. Saddle Phi_crit =
             219.4077 m; breach point (903622.5, 1118227.5). This is the one
             event where the full saddle machinery runs.
```

### lgp-grand-plateau  *(worked example — already complete)*

```
id:          lgp-grand-plateau
name:        Upper (LGP) lake drains into Grand Plateau Lake
status:      forecast
location:    Grand Plateau Glacier, Glacier Bay NP&P, SE Alaska
date:        2030            precision: year   (continuous crossing 2029.64)
source:      this project
data_status: complete

surface:     as above (shared grid)
bed:         as above (shared grid)
lake_level:  110.0 m  datum: ? — same question as above
dhdt:        -9.05 m/yr (shared)
notes:       Saddle Phi_crit = 271.5518 m; breach point (908752.5, 1112597.5).
             Lake sits ~93 m above Grand Plateau Lake.
```

---

### <copy this block for each new event>

```
id:
name:
status:          observed | forecast
location:
date:                            precision:   day | season | year
source:
data_status:     complete | partial | none

surface:         which DEM? what epoch? how far from the breach date?
                 vertical datum:
bed:             none | millan | millan+radar | local
                 if flat-bed: what defines the elevation?
lake_level:      value + datum, or "needs estimating"
dhdt:            value + source, or "unknown"
notes:
```

---

### historic event 1

```
id:
name:
status:          observed
date:                            precision:
source:
data_status:
surface:
bed:
lake_level:
dhdt:
notes:
```

### historic event 2

```
id:
name:
status:          observed
date:                            precision:
source:
data_status:
surface:
bed:
lake_level:
dhdt:
notes:
```

### historic event 3

```
id:
name:
status:          observed
date:                            precision:
source:
data_status:
surface:
bed:
lake_level:
dhdt:
notes:
```

---

## Open questions for you (answer inline, or ignore if not yet known)

1. **Flat-bed variants** — how should these be defined? Options I can see:
   at the lake-floor elevation; at a fixed depth below the lake surface;
   at the minimum of a Millan-derived bed. Which, and how many per event?

2. **Holdout validation** — this needs enough events to be meaningful. With 3
   historic events, calibrating on 2 and testing on 1 is weak. Worth saying now
   whether the list is likely to reach ~6+, since it changes the plan.

3. Any events you want **deliberately excluded**, and why? (Worth recording —
   an exclusion rule stated up front is defensible; one applied later isn't.)

---

*Once this is filled in I'll turn it into `docs/phase2_plan.md` in the repo,
with the event count settled and the holdout question decided, for your approval
before any code gets written.*
