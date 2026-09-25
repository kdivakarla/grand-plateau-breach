"""Read ICESat-2 ATL06 land-ice height granules.

ATL06 is the product to use over the Grand Plateau lakes. ATL13, the purpose-built
inland-water product, masks to **HydroLAKES**, which does not contain these
ice-dammed proglacial lakes, so it returns nothing here however good the track
(D-008). ATL06 is masked to land ice instead and does cover them.

Declared references, read from the granule rather than assumed:

    root/description    "land-ice surface heights (above the WGS 84 ellipsoid,
                         ITRF2014 reference frame)"
    dem/geoid_h         "Geoid height above WGS-84 reference ellipsoid ...
                         in the tide-free system"

Note the tide-system difference from ATL13, which reports its geoid in the
**mean-tide** system. Both carry a ``geoid_free2mean`` field for the conversion,
and they apply it in opposite directions -- read each field's description.
NAVD88/GEOID12 is tide-free, so ATL06's convention matches it and ATL13's does not.

Quality filtering is not optional
---------------------------------
ATL06 returns a height for many segments that are badly wrong. In the 2026-07-19
granule, **every one of the 733 segments inside the analysis domain** failed
``atl06_quality_summary`` and disagreed with ATL06's own reference DEM by 500 m
or more. Taken at face value they imply a 700 m elevation error; they are simply
cloud-corrupted. ``read_segments`` therefore applies two independent screens by
default:

1. ``atl06_quality_summary == 0`` -- the product's own flag.
2. ``|h_li - dem_h| <= dem_tolerance_m`` -- a sanity check against the reference
   DEM the granule carries. This catches gross outliers that slip the flag.

Both are defaults, not obligations: pass ``quality_only=False`` to inspect what
is being rejected. Do not disable them to obtain data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

BEAMS = ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")

_FIELDS = {
    "lat": "latitude",
    "lon": "longitude",
    "h_li": "h_li",
    "h_li_sigma": "h_li_sigma",
    "sigma_geo_h": "sigma_geo_h",
    "quality": "atl06_quality_summary",
    "segment_id": "segment_id",
    "delta_time": "delta_time",
    "geoid_h": "dem/geoid_h",
    "geoid_free2mean": "dem/geoid_free2mean",
    "dem_h": "dem/dem_h",
    "cloud_flg_atm": "geophysical/cloud_flg_atm",
    "bsnow_conf": "geophysical/bsnow_conf",
}


@dataclass
class ATL06Info:
    """Identity and declared vertical reference of one ATL06 granule."""

    path: Path
    product: str
    version: str
    doi: str
    rgt: int | None
    cycle: int | None
    time_start: str
    description: str
    geoid_description: str
    tide_system: str
    extras: dict[str, str] = field(default_factory=dict)

    def describe(self) -> str:
        import textwrap

        out = [
            f"{self.product} {self.version}   {self.path.name}",
            f"  RGT {self.rgt}  cycle {self.cycle}   acquired {self.time_start}",
            f"  doi {self.doi}",
            "",
            "  Height reference (from the granule):",
        ]
        out += ["      " + ln for ln in textwrap.wrap(self.description, 88)]
        out += ["", "  Geoid field (dem/geoid_h):"]
        out += ["      " + ln for ln in textwrap.wrap(self.geoid_description, 88)]
        out += ["", f"  Tide system: {self.tide_system}"]
        return "\n".join(out)


def _attr(obj, key: str, default: str = "") -> str:
    v = obj.attrs.get(key, default)
    return v.decode() if isinstance(v, bytes) else str(v)


def granule_info(path: str | Path) -> ATL06Info:
    """Read identity and the declared vertical reference."""
    import h5py

    path = Path(path)
    with h5py.File(path, "r") as f:
        desc = _attr(f, "description", "")
        geoid_desc = _attr(f[f"{BEAMS[0]}/land_ice_segments/dem/geoid_h"], "description", "")
        rgt = cycle = None
        if "orbit_info/rgt" in f:
            rgt = int(f["orbit_info/rgt"][0])
        if "orbit_info/cycle_number" in f:
            cycle = int(f["orbit_info/cycle_number"][0])
        tide = "tide-free" if "tide-free" in geoid_desc.lower() else "UNSTATED - check"
        return ATL06Info(
            path=path,
            product=_attr(f, "identifier_product_type", "?"),
            version=_attr(f, "identifier_product_format_version", "?"),
            doi=_attr(f, "identifier_product_doi", "?"),
            rgt=rgt,
            cycle=cycle,
            time_start=_attr(f, "time_coverage_start", "?"),
            description=desc,
            geoid_description=geoid_desc,
            tide_system=tide,
        )


def read_segments(
    path: str | Path,
    bbox: tuple[float, float, float, float] | None = None,
    beams: tuple[str, ...] = BEAMS,
    quality_only: bool = True,
    dem_tolerance_m: float = 50.0,
) -> pd.DataFrame:
    """Read land-ice segments, quality-screened by default.

    See the module docstring: on a cloud-affected granule the unscreened heights
    can be wrong by hundreds of metres while still looking like plausible numbers.

    Adds ``h_ortho = h_li - geoid_h`` (EGM2008, tide-free) and ``vs_dem``, the
    disagreement with the granule's own reference DEM.
    """
    import h5py

    frames = []
    with h5py.File(path, "r") as f:
        for beam in beams:
            base = f"{beam}/land_ice_segments"
            if base not in f:
                continue
            g = f[base]
            data = {"beam": beam}
            for out_name, key in _FIELDS.items():
                if key not in g:
                    continue
                v = g[key][:].astype(np.float64)
                v[np.abs(v) > 1e30] = np.nan
                data[out_name] = v
            df = pd.DataFrame(data)
            if bbox is not None:
                lon0, lat0, lon1, lat1 = bbox
                df = df[df.lat.between(lat0, lat1) & df.lon.between(lon0, lon1)]
            frames.append(df)

    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out = out[np.isfinite(out.h_li)]
    out["h_ortho"] = out.h_li - out.geoid_h
    out["vs_dem"] = out.h_li - out.dem_h

    if quality_only:
        n0 = len(out)
        keep = out.quality == 0
        if dem_tolerance_m is not None:
            keep &= out.vs_dem.abs() <= dem_tolerance_m
        out = out[keep]
        if n0 and not len(out):
            log.warning(
                "all %d segments were rejected by quality screening -- this granule "
                "is likely cloud-affected over the requested area",
                n0,
            )
        else:
            log.info("quality screening kept %d of %d segments", len(out), n0)
    return out.reset_index(drop=True)


def quality_by_band(
    path: str | Path, bands: list[tuple[float, float]], beams: tuple[str, ...] = BEAMS
) -> pd.DataFrame:
    """Fraction of segments passing quality in each latitude band.

    Use this before trusting -- or discarding -- a granule: it separates "the
    track missed" from "the data is there but cloud-wiped", which look identical
    if you only ever read the screened output.
    """
    raw = read_segments(path, beams=beams, quality_only=False)
    rows = []
    for lat0, lat1 in bands:
        s = raw[raw.lat.between(lat0, lat1)]
        good = s[(s.quality == 0) & (s.vs_dem.abs() <= 50.0)]
        rows.append(
            {
                "lat_min": lat0,
                "lat_max": lat1,
                "n": len(s),
                "n_good": len(good),
                "pct_good": 100.0 * len(good) / len(s) if len(s) else np.nan,
                "cloud_flg_median": s.cloud_flg_atm.median() if len(s) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def navd88_from_ellipsoid(h_li, geoid12_height_m: float):
    """Approximate NAVD88 from ATL06 ellipsoid height. See io.atl13 for caveats."""
    return np.asarray(h_li, dtype=float) - float(geoid12_height_m)
