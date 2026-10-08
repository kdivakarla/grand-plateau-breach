# QGIS integration

Three ways to work with QGIS, in increasing order of effort.

## 1. Just open the outputs (nothing to install)

Everything the pipeline writes is already QGIS-ready: GeoTIFF and GeoPackage in
**ESRI:102247**, the project CRS, on the canonical 5 m grid. Drag them in from
`data/standard/`, `data/derived/` and `data/interim/`.

## 2. Load them all, styled, from the Python Console

`load_cascade_layers.py` adds every cascade output as a styled layer in one go —
the margin raster gets a diverging ramp pinned at 0, which is the flotation
threshold and the only hard break on that layer.

Edit `REPO` at the top, then: **Plugins → Python Console → Show Editor → Open
Script → Run**. It reports anything missing so you know which pipeline step to run.

## 3. Run the analysis from QGIS (Processing script)

The pipeline is a normal Python package, so a Processing script can shell out to
it. The catch is that QGIS ships its own Python, which does **not** have the
`gpbreach` environment. Call the conda environment's interpreter explicitly:

```python
import subprocess
PY = "/opt/miniconda3/envs/gpbreach/bin/python"
subprocess.run([PY, "-m", "gpbreach.cascade.run_scenarios",
                "--beds", "millan", "iceboost"],
               cwd="/path/to/grand-plateau-breach", check=True)
```

Then reload the layers with the script above.

**Worth knowing:** QGIS's bundled GDAL and the conda environment's are different
builds. Mixing them in one process causes hard-to-diagnose crashes, which is why
this shells out rather than importing `gpbreach` into QGIS directly. Keep the two
Pythons separate and pass data between them as files.
