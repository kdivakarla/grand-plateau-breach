# data/

Gitignored. Nothing here is edited by Claude Code (CLAUDE.md Section 2.2).

- `raw/` — read-only originals. **Currently empty.** The baseline reads the
  owner's rasters in place from the parent `Grand_Plateau_Master/` directory,
  via the `../`-relative paths in `configs/*.yaml`. Those paths work on the
  owner's machine but **not in a fresh clone** — the rasters are ~250 MB and are
  deliberately not in git. To run from a clone, put
  `2010_IFSAR_GLACIER_SURFACE_CLIPPED.tif` and `Bedrock_Aligned.tif` here (copy
  or symlink) and point `dam.build.surface` / `dam.build.bed` at `data/raw/`.
- `interim/` — reprojected/standardized intermediates. The only directory Claude
  Code may delete from, along with `results/scratch/`.
- `processed/` — dam objects (`dam_*.npz`) and flotation-fraction curves
  (`fcurve_*.npz`), both rebuildable with `gpbreach preprocess` / `gpbreach fcurve`.
