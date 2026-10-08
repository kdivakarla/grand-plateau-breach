"""Load the cascade outputs into QGIS, styled.

Paste into the QGIS Python Console (Plugins -> Python Console), or run with
Plugins -> Python Console -> Show Editor -> Open Script.

Everything the pipeline writes is already QGIS-ready -- GeoTIFF and GeoPackage in
ESRI:102247, the project CRS -- so this only saves the clicking and applies
sensible styling. Nothing here modifies the data.
"""

from pathlib import Path

from qgis.core import (
    QgsColorRampShader,
    QgsProject,
    QgsRasterLayer,
    QgsRasterShader,
    QgsSingleBandPseudoColorRenderer,
    QgsVectorLayer,
)
from qgis.PyQt.QtGui import QColor

# --- EDIT THIS to your checkout -------------------------------------------
REPO = Path("/Users/krishnadivakarla/work/Research/Glaciers/Grand_Plateau_Master"
            "/grand-plateau-breach")
# ---------------------------------------------------------------------------

STANDARD = REPO / "data" / "standard"
DERIVED = REPO / "data" / "derived"
INTERIM = REPO / "data" / "interim"


def _diverging(layer, low, mid, high, colors=("#ca0020", "#f7f7f7", "#2a78d6")):
    """Diverging ramp pinned so the midpoint lands on a meaningful value."""
    shader = QgsRasterShader()
    ramp = QgsColorRampShader()
    ramp.setColorRampType(QgsColorRampShader.Interpolated)
    ramp.setColorRampItemList([
        QgsColorRampShader.ColorRampItem(low, QColor(colors[0]), f"{low:.0f}"),
        QgsColorRampShader.ColorRampItem(mid, QColor(colors[1]), f"{mid:.0f}"),
        QgsColorRampShader.ColorRampItem(high, QColor(colors[2]), f"{high:.0f}"),
    ])
    shader.setRasterShaderFunction(ramp)
    layer.setRenderer(QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader))
    layer.triggerRepaint()


def load() -> None:
    project = QgsProject.instance()
    loaded, missing = [], []

    # Rasters: (path, display name, styling)
    rasters = [
        (DERIVED / "margin_lgp_arcticdem_20180903.tif", "margin vs LGP (2018)", "margin"),
        (DERIVED / "head_arcticdem_20180903.tif", "hydraulic head (2018)", None),
        (STANDARD / "bed_ellip.tif", "bed — Millan", None),
        (STANDARD / "bed_iceboost.tif", "bed — IceBoost", None),
        (STANDARD / "surf_2018.tif", "ice surface 2018", None),
        (STANDARD / "classes.tif", "classes", None),
    ]
    for path, name, style in rasters:
        if not path.exists():
            missing.append(path.name)
            continue
        lyr = QgsRasterLayer(str(path), name)
        if not lyr.isValid():
            missing.append(f"{path.name} (invalid)")
            continue
        if style == "margin":
            # Zero is flotation: the only value on this layer worth a hard break.
            _diverging(lyr, -100.0, 0.0, 300.0)
        project.addMapLayer(lyr)
        loaded.append(name)

    vectors = [
        (DERIVED / "path_arcticdem_20180903.gpkg", "drainage path (2018)"),
        (REPO / "data" / "processed" / "lakes_2018.gpkg", "lakes 2018"),
        (REPO / "data" / "processed" / "ice_2018.gpkg", "ice 2018"),
        (INTERIM / "rgi7_lake_terminating_4.gpkg", "RGI 7 — the four dam glaciers"),
    ]
    for path, name in vectors:
        if not path.exists():
            missing.append(path.name)
            continue
        lyr = QgsVectorLayer(str(path), name, "ogr")
        if lyr.isValid():
            project.addMapLayer(lyr)
            loaded.append(name)
        else:
            missing.append(f"{path.name} (invalid)")

    print(f"loaded {len(loaded)} layer(s):")
    for n in loaded:
        print("   ", n)
    if missing:
        print(f"\nnot found ({len(missing)}) — run the relevant pipeline step:")
        for n in missing:
            print("   ", n)


load()
