"""Read ICESat-2 ATL13 inland water surface height granules.

ATL13 gives lake surface elevations, which is the one input D-008 leaves open:
the vertical datum of the lake levels used in the baseline.

**The datum answer is in the file, and this module reads it rather than
asserting it.** For the v007 granule in `data/raw/`, the dataset attributes say:

    ht_water_surf   "Water surface height ... with reference to WGS84 ellipsoid"
    ht_ortho        "Orthometric height EGM2008 converted from ellipsoidal height."
    segment_geoid   "Applicable mean-tide system geoid value ..."

So ATL13 is natively **WGS84 ellipsoid**, with an **EGM2008** orthometric height
provided alongside.

The trap
--------
``ht_ortho`` is EGM2008. The project's IFSAR surface is **NAVD88** (D-008), and
NAVD88 in Alaska is realised through GEOID12A/12B, *not* EGM2008. The two geoid
models disagree by ~1.2 m at Grand Plateau. Treating ``ht_ortho`` as NAVD88
therefore injects a systematic ~1 m lake-level error, which is ~0.15 yr of breach
date at the project's sensitivity of -0.12 yr per metre.

To get NAVD88, take the *ellipsoid* height and subtract a GEOID12 value from
NOAA (see ``navd88_from_ellipsoid``). Two caveats on that conversion, both
recorded in D-008:

* GEOID12 expects **NAD83** ellipsoid heights; ICESat-2 is on ITRF2014. The two
  differ by of order a metre in Alaska. An exact conversion needs NOAA's HTDP or
  NCAT; this module does not pretend to do it.
* ATL13's geoid is **mean-tide**; NAVD88/GEOID12 is **tide-free**. The granule
  carries ``segment_geoid_free2mean`` for that step.

Nothing here writes to ``data/raw/`` -- it is read-only (CLAUDE.md Section 2.2).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

BEAMS = ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")

#: Fields pulled per segment. Keep this narrow: an ATL13 granule is ~60 MB and
#: most of it is quality flags nobody has asked for yet.
SEGMENT_FIELDS = (
    "segment_lat",
    "segment_lon",
    "ht_water_surf",
    "ht_ortho",
    "segment_geoid",
    "segment_geoid_free2mean",
    "stdev_water_surf",
    "err_ht_water_surf",
    "inland_water_body_id",
    "inland_water_body_type",
    "ice_flag",
    "qf_cloud",
    "qf_ice",
    "delta_time",
    "cycle_number",
    "rgt",
)

#: Height fields whose declared datum we report verbatim from the file.
DATUM_FIELDS = (
    "ht_water_surf",
    "ht_ortho",
    "segment_geoid",
    "segment_geoid_free2mean",
    "segment_dem_ht",
)


@dataclass
class GranuleInfo:
    """Identity and declared datums of one ATL13 granule."""

    path: Path
    product: str
    version: str
    doi: str
    time_start: str
    time_end: str
    date_created: str
    datums: dict[str, str] = field(default_factory=dict)
    units: dict[str, str] = field(default_factory=dict)

    def describe(self) -> str:
        lines = [
            f"{self.product} {self.version}   {self.path.name}",
            f"  acquired  {self.time_start}  ..  {self.time_end}",
            f"  created   {self.date_created}",
            f"  doi       {self.doi}",
            "",
            "  Declared vertical references (read from the file, not assumed):",
        ]
        for name, desc in self.datums.items():
            unit = self.units.get(name, "?")
            lines.append(f"    {name:26s} [{unit}]")
            lines.append(f"        {desc}")
        return "\n".join(lines)


def _attr(obj, key: str, default: str = "") -> str:
    v = obj.attrs.get(key, default)
    return v.decode() if isinstance(v, bytes) else str(v)


def granule_info(path: str | Path) -> GranuleInfo:
    """Read the granule's identity and the declared datum of every height field."""
    import h5py

    path = Path(path)
    with h5py.File(path, "r") as f:
        info = GranuleInfo(
            path=path,
            product=_attr(f, "identifier_product_type", "?"),
            version=_attr(f, "identifier_product_format_version", "?"),
            doi=_attr(f, "identifier_product_doi", "?"),
            time_start=_attr(f, "time_coverage_start", "?"),
            time_end=_attr(f, "time_coverage_end", "?"),
            date_created=_attr(f, "date_created", "?"),
        )
        for name in DATUM_FIELDS:
            key = f"{BEAMS[0]}/{name}"
            if key in f:
                info.datums[name] = _attr(f[key], "description", "(no description)")
                info.units[name] = _attr(f[key], "units", "?")
    return info


def read_segments(
    path: str | Path,
    bbox: tuple[float, float, float, float] | None = None,
    beams: tuple[str, ...] = BEAMS,
    drop_fill: bool = True,
) -> pd.DataFrame:
    """Read water-surface segments into a DataFrame.

    Parameters
    ----------
    bbox:
        ``(lon_min, lat_min, lon_max, lat_max)`` in degrees. ``None`` reads the
        whole granule, which is global and large.
    drop_fill:
        Replace each field's declared ``_FillValue`` with NaN. ATL13 uses 127 as
        the int8 fill, which is easily mistaken for a real flag value -- an
        ``ice_flag`` of 127 means "not set", not "ice".
    """
    import h5py

    frames = []
    with h5py.File(path, "r") as f:
        for beam in beams:
            if beam not in f:
                continue
            g = f[beam]
            lat, lon = g["segment_lat"][:], g["segment_lon"][:]
            if bbox is not None:
                lon0, lat0, lon1, lat1 = bbox
                keep = (lat >= lat0) & (lat <= lat1) & (lon >= lon0) & (lon <= lon1)
            else:
                keep = np.ones(lat.shape, dtype=bool)
            if not keep.any():
                continue
            idx = np.flatnonzero(keep)
            data = {"beam": beam}
            for name in SEGMENT_FIELDS:
                if name not in g:
                    continue
                values = g[name][:][idx]
                if drop_fill:
                    fill = g[name].attrs.get("_FillValue")
                    if fill is not None:
                        values = np.where(values == fill, np.nan, values.astype(np.float64))
                data[name] = values
            frames.append(pd.DataFrame(data))
    if not frames:
        return pd.DataFrame(columns=["beam", *SEGMENT_FIELDS])
    return pd.concat(frames, ignore_index=True)


def navd88_from_ellipsoid(ellipsoid_height_m, geoid12_height_m: float):
    """Convert an ellipsoid height to approximate NAVD88.

    ``geoid12_height_m`` comes from NOAA for the site, e.g.::

        https://geodesy.noaa.gov/api/geoid/ght?lat=59.02&lon=-137.94&model=12

    Approximate: ICESat-2 is ITRF2014 while GEOID12 expects NAD83, and the tide
    systems differ. Use for reconnaissance, not for a published number -- see
    this module's docstring and D-008.
    """
    return np.asarray(ellipsoid_height_m, dtype=float) - float(geoid12_height_m)


def summarise_water_bodies(df: pd.DataFrame, geoid12_height_m: float | None = None) -> pd.DataFrame:
    """One row per water body: location, elevation in each available datum."""
    if df.empty:
        return pd.DataFrame()
    out = []
    for wbid, s in df.groupby("inland_water_body_id"):
        row = {
            "water_body_id": int(wbid),
            "n_segments": len(s),
            "beams": ",".join(sorted(s["beam"].unique())),
            "lat": s["segment_lat"].mean(),
            "lon": s["segment_lon"].mean(),
            "ellipsoid_m": s["ht_water_surf"].median(),
            "egm2008_m": s["ht_ortho"].median(),
            "geoid_egm2008_m": s["segment_geoid"].median(),
            "spread_m": s["ht_water_surf"].max() - s["ht_water_surf"].min(),
            "sd_m": s["ht_water_surf"].std(),
        }
        if geoid12_height_m is not None:
            row["navd88_est_m"] = row["ellipsoid_m"] - geoid12_height_m
            row["egm2008_minus_navd88_m"] = row["egm2008_m"] - row["navd88_est_m"]
        out.append(row)
    return pd.DataFrame(out).sort_values("n_segments", ascending=False).reset_index(drop=True)
