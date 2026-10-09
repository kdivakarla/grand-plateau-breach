# What is needed to say how likely the outlet is to fail

Gap analysis, 2026-10-09. Written against the state at commit `d2c86e5`.

---

## 1. Where the analysis actually stands

The chain runs from ice geometry to a peak inflow. It does **not** yet reach the
outlet, and the inflow it produces spans a factor of **17**:

| bed | sill | mode | Q_p lo | Q_p hi |
|---|---|---|---|---|
| millan | along-path | tunnel | 816 | 1,519 |
| millan | along-path | non-tunnel | 7,480 | 11,320 |
| millan | bed-only | tunnel | 1,078 | 2,011 |
| millan | bed-only | non-tunnel | 9,008 | 13,649 |
| iceboost | bed-only | tunnel | 1,030 | 1,918 |
| iceboost | bed-only | non-tunnel | 8,736 | 13,229 |

**816 – 13,649 m³/s.** A failure probability computed on that would be precise
and meaningless: it would mostly report which branch was assumed.

### The spread is not statistical
Most of it comes from **choices that are still open**, not from measurement error.
Collapsing each one:

| decision | best-case factor remaining |
|---|---|
| **drainage mode** (tunnel vs non-tunnel) | **1.8×** |
| sill definition (along-path vs bed-only) | 13.3× |
| which bed (Millan vs IceBoost) | 12.8× |
| bathymetry, i.e. V_w within a branch | 1.9× |

**One decision dominates everything.** Settling tunnel versus non-tunnel takes
the peak discharge from a factor of 17 to a factor of 1.8. Nothing else comes
close, and no amount of Monte Carlo over parameters substitutes for it — sampling
`k` and `dh/dt` while the drainage mode is unresolved produces a confident
distribution around the wrong branch.

---

## 2. What is missing, in priority order

### A. Resolve the structural choices — these are not distributions

| # | question | why it is not a sampling problem | what resolves it |
|---|---|---|---|
| A1 | **tunnel or non-tunnel drainage?** | a binary physical mechanism, 17× → 1.8× | **conduit model (M6)**, and the hindcasts |
| A2 | **which sill governs?** | a modelling definition (D-019) | owner decision; or make the conduit model route the flood and find out |
| A3 | ~~does the spill basin pond or drain?~~ | **RESOLVED D-020: it drains.** New question: does the flood *bypass* the outlet? | — |
| A4 | **which bed?** | IceBoost datum **confirmed ellipsoidal** (D-020), so the +105 m offset is real | validate against radar |

A1 is the one to do first, and it is already specified: **`conduit.py` (M6)** solves
melt-opening against creep-closure and *predicts* whether a tunnel can enlarge fast
enough. It turns A1 from an assumption into a result. It is the single highest-value
piece of unbuilt work.

**The hindcasts matter here too.** The spec's Phase 2 — three past river-piracy
events — would calibrate `k` *and* show empirically which mode these systems
exhibit. It has never been started, and it is the only external check the project
has.

### B. Data that does not exist yet — without it there is no failure criterion

| # | needed | for | status |
|---|---|---|---|
| B1 | outlet **crest elevation** and **width** | the rating curve `Q = C·L·H^1.5` | extractable — `outletOutline.gpkg` touches GPL and is inside the extended grid |
| B2 | ~~boulder size~~ **D50 = 0.961 m** supplied 2026-10-09 | critical shear — thresholds now computed, D-020 | **done**; D84 would still help |
| B3 | ~~outlet slope~~ **S = 0.16** supplied 2026-10-09 (D-021) | applied shear | **done**; cross-section/width still needed |
| B4 | LGP **bathymetry** | tightens V_w | missing; worth only 1.9× — **do this last** |

B2 is the hard blocker. Without a grain size there is no threshold to compare
against, and no failure criterion at all.

### C. Physics modules not yet built

1. **`conduit.py` (M6)** — coupled two-lake ODE. Also answers A1.
2. **GPL routing** — and it must be a **two-reservoir** model. GPL tops its divide
   1.7–23 h into the flood in *every* scenario, always inside the first quarter,
   after which stage is pinned near 29.6 m while a 16.6 km² basin absorbs ~309 Mm³.
   A single level-pool reservoir will badly overstate the stage the outlet sees.
3. **Outlet rating** — and it cannot stay fixed. Passing even the *smallest*
   flood needs ~1.1 m of head over a 400 m crest, and GPL can only rise ~1.5 m
   before spilling. Every larger case requires the outlet to enlarge, which is
   the failure being predicted. The rating must couple to erosion.
4. **Armor stability** — applied versus critical shear, with two formulations
   (Shields depth-slope, and a stream-power or critical-unit-discharge form that
   travels better on steep boulder beds) reported side by side.

### D. Only then, the probabilistic layer

Once the branches are closed and the modules exist:

- distributions on what remains genuinely uncertain — `k`, `dh/dt`, lake levels,
  bed error, `C`, `θ_c`, grain size, conduit roughness and initial size;
- the Monte Carlo wrapper, seeded, one row per realization (the engine is already
  built for this: deterministic is N = 1);
- **a failure criterion that is not binary.** Armor failure is progressive —
  incipient motion, then winnowing, then incision, then runaway. Report
  `P(mobilisation)` and, conditional on that, `P(runaway)`. A single
  fails/doesn't-fail number hides the mechanism.

---

## 3. What can honestly be said, and when

**Now (updated 2026-10-09):** with D50 and slope in hand, a first branch-wise
verdict exists — under non-tunnel drainage the armor mobilises with a large
margin; under tunnel drainage it sits near threshold (D-021). Peak inflow is
still known only to a factor of ~17, and that same unresolved branch is what
separates the two verdicts.

**After A1 + A3 and modules C1–C2:** a conditional statement — *"given tunnel
drainage, GPL's stage reaches X and the outlet sees Q_out of Y"*. Useful, and
honest about what it is conditional on.

**After B1–B3 and C3–C4:** a deterministic verdict on mobilisation under each
branch. This is probably the most defensible product the project can reach
without new field data.

**After D:** a probability — but one whose credibility rests entirely on the
hindcasts. Without them there is no validation of the chain against any observed
event, and the distribution would be an expression of assumed inputs rather than
a tested forecast.

**A recommendation.** Aim for the conditional, branch-wise answer first and treat
the unconditional probability as a later product. It is more useful to a hazard
assessment to say *"the outlet mobilises under non-tunnel drainage but holds under
tunnel drainage, and here is which one the physics supports"* than to report a
single number that averages over a question nobody has answered.

---

## 4. Suggested order

1. **A3** — spill basin: ponds or drains. Cheap, and gates everything downstream.
2. **C1 / A1** — build `conduit.py`; it resolves the dominant unknown.
3. **B1, B3** — outlet crest, width, slope from existing data.
4. **C2** — two-reservoir routing.
5. **B2** — boulder size. Field data or photogrammetry; the long pole.
6. **C3, C4** — coupled rating and armor stability.
7. **Phase 2 hindcasts** — in parallel from the start; they calibrate `k`, test the
   drainage mode, and are the only external validation available.
8. **D** — Monte Carlo, last.
