"""Read-only GeoTIFF access.

Three backends are tried in order so the pipeline runs both in a fully
provisioned environment and in the current ``gpgn-318`` environment, which has
no geospatial stack installed:

1. ``rasterio``      -- preferred; declared in environment.yml
2. ``osgeo.gdal``    -- equivalent, if GDAL python bindings are present
3. ``gdal_translate``-- subprocess fallback that converts the raster to a flat
   ENVI binary read with ``numpy.fromfile``. Used only when neither library is
   importable; it finds ``gdal_translate`` on PATH or inside QGIS.app.

Rasters are never written by this module. ``data/raw`` is read-only
(CLAUDE.md Section 2.2).
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

_QGIS_BIN = Path("/Applications/QGIS.app/Contents/MacOS")


@dataclass(frozen=True)
class RasterData:
    """A single-band raster loaded as float64 with nodata set to NaN."""

    values: np.ndarray
    transform: tuple[float, float, float, float, float, float]
    crs: str | None
    path: Path

    @property
    def shape(self) -> tuple[int, int]:
        return self.values.shape

    def xy(self, row: np.ndarray | int, col: np.ndarray | int):
        """Map (row, col) to projected (x, y) cell centres."""
        x0, dx, _, y0, _, dy = self.transform
        return x0 + (np.asarray(col) + 0.5) * dx, y0 + (np.asarray(row) + 0.5) * dy

    def grid_key(self) -> tuple:
        """Identity of the grid, for asserting that two rasters are co-registered."""
        return (self.shape, tuple(round(v, 6) for v in self.transform))


def read_raster(path: str | Path) -> RasterData:
    """Load band 1 of ``path`` as float64 with nodata represented as NaN."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    for backend in (_read_rasterio, _read_gdal, _read_gdal_translate):
        try:
            return backend(path)
        except ImportError:
            continue
    raise RuntimeError(
        f"No raster backend available to read {path}. Install rasterio "
        "(see environment.yml) or make gdal_translate available on PATH."
    )


def _finalise(values: np.ndarray, nodata: float | None) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if nodata is not None and np.isfinite(nodata):
        values = np.where(values == nodata, np.nan, values)
    # GDAL's float32 nodata sentinels (+/-3.4e38) survive as huge magnitudes.
    return np.where(np.abs(values) > 1e30, np.nan, values)


def _read_rasterio(path: Path) -> RasterData:
    import rasterio

    with rasterio.open(path) as src:
        t = src.transform
        return RasterData(
            values=_finalise(src.read(1), src.nodata),
            transform=(t.c, t.a, t.b, t.f, t.d, t.e),
            crs=str(src.crs) if src.crs else None,
            path=path,
        )


def _read_gdal(path: Path) -> RasterData:
    from osgeo import gdal

    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    band = ds.GetRasterBand(1)
    gt = ds.GetGeoTransform()
    return RasterData(
        values=_finalise(band.ReadAsArray(), band.GetNoDataValue()),
        transform=(gt[0], gt[1], gt[2], gt[3], gt[4], gt[5]),
        crs=ds.GetProjection() or None,
        path=path,
    )


def _find_gdal_translate() -> str:
    exe = shutil.which("gdal_translate")
    if exe:
        return exe
    candidate = _QGIS_BIN / "gdal_translate"
    if candidate.exists():
        return str(candidate)
    raise ImportError("gdal_translate not found")


def _read_gdal_translate(path: Path) -> RasterData:
    exe = _find_gdal_translate()
    log.warning(
        "Reading %s via the gdal_translate subprocess fallback; install rasterio "
        "for a supported path.",
        path.name,
    )
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "band1.bin"
        subprocess.run(
            [exe, "-q", "-of", "ENVI", "-ot", "Float64", "-b", "1", str(path), str(out)],
            check=True,
            capture_output=True,
        )
        hdr = (out.with_suffix(".hdr")).read_text()
        nx = int(re.search(r"^samples\s*=\s*(\d+)", hdr, re.M).group(1))
        ny = int(re.search(r"^lines\s*=\s*(\d+)", hdr, re.M).group(1))
        values = np.fromfile(out, dtype="<f8").reshape(ny, nx)
        mc = re.search(r"map info = \{([^}]*)\}", hdr)
        parts = [p.strip() for p in mc.group(1).split(",")] if mc else []
        # ENVI map info: name, ref_x(1-based), ref_y(1-based), easting, northing, dx, dy, ...
        x0 = float(parts[3]) - (float(parts[1]) - 1) * float(parts[5])
        y0 = float(parts[4]) + (float(parts[2]) - 1) * float(parts[6])
        transform = (x0, float(parts[5]), 0.0, y0, 0.0, -float(parts[6]))
        cs = re.search(r"coordinate system string = \{(.*)\}", hdr, re.S)
    return RasterData(
        values=_finalise(values, None),
        transform=transform,
        crs=cs.group(1).strip() if cs else None,
        path=path,
    )
