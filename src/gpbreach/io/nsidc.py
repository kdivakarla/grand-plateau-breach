"""Search and download ICESat-2 granules from NSIDC.

Search uses NASA's CMR API and needs no credentials. Download goes through
Earthdata Login, which does need them: the script reads ``~/.netrc`` rather than
ever handling a password itself.

Why ATL06 and not ATL13
-----------------------
ATL13 is the inland-water product and is the obvious choice for lake levels, but
it does not work at Grand Plateau. Its reference water mask is **HydroLAKES**,
which does not contain these ice-dammed proglacial lakes, so ATL13 returns zero
water segments over all three of them regardless of how good the ground track is
(D-008). ATL06 is masked to land ice instead and does cover them; over a lake
ATL13 *does* carry, the two agree to 0.10-0.25 m, so ATL06 measures water surfaces
perfectly well here.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)

CMR_GRANULES = "https://cmr.earthdata.nasa.gov/search/granules.json"


@dataclass(frozen=True)
class Granule:
    """One granule's identity and where to fetch it."""

    name: str
    time_start: str
    url: str
    size_mb: float | None = None

    @property
    def date(self) -> str:
        return self.time_start[:10]

    @property
    def rgt(self) -> int | None:
        """Reference ground track, parsed from the ICESat-2 filename convention."""
        parts = self.name.split("_")
        return int(parts[2][:4]) if len(parts) > 2 and len(parts[2]) >= 8 else None

    @property
    def cycle(self) -> int | None:
        parts = self.name.split("_")
        return int(parts[2][4:6]) if len(parts) > 2 and len(parts[2]) >= 8 else None


def search_granules(
    short_name: str,
    bbox: tuple[float, float, float, float],
    start: str | None = None,
    end: str | None = None,
    months: tuple[int, ...] | None = None,
    page_size: int = 500,
) -> list[Granule]:
    """Find granules intersecting ``bbox``.

    ``bbox`` is ``(lon_min, lat_min, lon_max, lat_max)``. ``months`` keeps only
    granules whose acquisition falls in those months, which is how a melt-season
    series is selected -- CMR itself has no such filter.
    """
    import requests

    params = {
        "short_name": short_name,
        "bounding_box": ",".join(str(v) for v in bbox),
        "page_size": page_size,
        "sort_key": "start_date",
    }
    if start or end:
        params["temporal"] = f"{start or ''},{end or ''}"

    r = requests.get(CMR_GRANULES, params=params, timeout=60)
    r.raise_for_status()
    out = []
    for e in r.json()["feed"]["entry"]:
        href = next(
            (
                link["href"]
                for link in e.get("links", [])
                if link.get("href", "").endswith(".h5") and "s3://" not in link.get("href", "")
            ),
            None,
        )
        if href is None:
            continue
        if months and int(e["time_start"][5:7]) not in months:
            continue
        size = e.get("granule_size")
        out.append(
            Granule(
                name=e["producer_granule_id"],
                time_start=e["time_start"],
                url=href,
                size_mb=float(size) if size else None,
            )
        )
    out.sort(key=lambda g: g.time_start)
    return out


def _check_credentials() -> None:
    netrc = Path(os.environ.get("NETRC", Path.home() / ".netrc"))
    if not netrc.exists():
        raise SystemExit(
            f"{netrc} not found. Earthdata Login is required to download from NSIDC.\n"
            "Create one with:\n"
            "  machine urs.earthdata.nasa.gov login <user> password <pass>\n"
            "then: chmod 600 ~/.netrc"
        )
    if "urs.earthdata.nasa.gov" not in netrc.read_text():
        raise SystemExit(f"{netrc} has no urs.earthdata.nasa.gov entry.")


def download_granule(granule: Granule, dest_dir: str | Path, overwrite: bool = False) -> Path:
    """Fetch one granule, skipping it if already present.

    Earthdata answers with a redirect chain to its login service; ``requests``
    applies the ``~/.netrc`` credentials on the way through. Downloads land in a
    ``.part`` file and are renamed on success, so an interrupted run never leaves
    a truncated granule that later looks complete.
    """
    import requests

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / granule.name
    if dest.exists() and not overwrite:
        log.info("have %s", granule.name)
        return dest

    _check_credentials()
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.Session() as s:
        resp = s.get(granule.url, stream=True, timeout=300, allow_redirects=True)
        if resp.status_code in (401, 403):
            raise SystemExit(
                f"Earthdata rejected the credentials for {granule.name} "
                f"(HTTP {resp.status_code}). Check ~/.netrc, and that the account has "
                "accepted the NSIDC data-use agreement at urs.earthdata.nasa.gov."
            )
        resp.raise_for_status()
        with tmp.open("wb") as fh:
            for chunk in resp.iter_content(1 << 20):
                fh.write(chunk)
    tmp.rename(dest)
    log.info("downloaded %s (%.1f MB)", granule.name, dest.stat().st_size / 1e6)
    return dest
