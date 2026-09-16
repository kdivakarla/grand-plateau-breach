"""Baseline regression: the pipeline must reproduce the owner's QGIS results.

Three levels of assertion, weakest to strongest:

1. the headline breach year (2030 / 2035);
2. the continuous crossing and the saddle value;
3. the full legacy raster, cell for cell, for every year the owner produced.

(3) is the real test. It compares against
``Buoyant_Maps_*/Buoyancy_State_<year>.tif`` directly, so it fails if any part
of the port -- geometry, constants, sign, or time step -- drifts.

If these fail, do not tune anything. Follow CLAUDE.md Section 8.4.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from gpbreach.components.criterion import FlotationConnectivity
from gpbreach.config import Config, repo_root
from gpbreach.dam import DamObject
from gpbreach.engine import run

EXPECTED_PATH = Path(__file__).parent / "expected" / "baseline.yaml"
CASES = yaml.safe_load(EXPECTED_PATH.read_text())


def _result(case: dict):
    cfg = Config.load(case["config"])
    dam_path = cfg.dam_path
    if not dam_path.exists():
        pytest.skip(f"{dam_path} missing; run `gpbreach preprocess {cfg.path}` first")
    return cfg, DamObject.load(dam_path), run(cfg, DamObject.load(dam_path))


@pytest.mark.parametrize("name", sorted(CASES))
def test_breach_year(name: str) -> None:
    case = CASES[name]
    _, _, result = _result(case)
    assert int(result.breach_year.iloc[0]) == case["breach_year"]


@pytest.mark.parametrize("name", sorted(CASES))
def test_continuous_crossing_and_saddle(name: str) -> None:
    case = CASES[name]
    tol = case["tolerance"]
    _, _, result = _result(case)
    s = result.samples.iloc[0]
    assert s["phi_crit_m"] == pytest.approx(case["phi_crit_m"], abs=tol["phi_crit_m"])
    assert s["breach_year_continuous"] == pytest.approx(
        case["breach_year_continuous"], abs=tol["breach_year_continuous"]
    )


@pytest.mark.parametrize("name", sorted(CASES))
def test_breach_point_location(name: str) -> None:
    case = CASES[name]
    _, _, result = _result(case)
    assert result.timeseries.attrs["at_saddle"] is True
    assert list(result.timeseries.attrs["breach_cell"]) == case["breach_cell"]


@pytest.mark.parametrize("name", sorted(CASES))
def test_connectivity_choice_does_not_change_the_answer(name: str) -> None:
    """4- and 8-connectivity give the same baseline; the choice is not load-bearing."""
    case = CASES[name]
    cfg, dam, _ = _result(case)
    f = cfg.nominal("flotation_fraction")
    values = {
        c: FlotationConnectivity(connectivity=c).critical_threshold_for_dam(dam, f) for c in (4, 8)
    }
    assert values[4] == pytest.approx(values[8], abs=1e-9)


@pytest.mark.parametrize("name", sorted(CASES))
def test_matches_legacy_rasters_every_year(name: str) -> None:
    """Cell-for-cell equality with the owner's QGIS output, for every year."""
    from gpbreach.io import read_raster

    case = CASES[name]
    folder = (repo_root() / case["legacy_rasters"]).resolve()
    if not folder.is_dir():
        pytest.skip(f"legacy rasters not available at {folder}")

    cfg, dam, _ = _result(case)
    criterion = FlotationConnectivity(connectivity=int(cfg.components.get("connectivity", 8)))
    f = cfg.nominal("flotation_fraction")
    level = cfg.nominal("lake_level_m")
    dhdt = cfg.nominal("dhdt_m_per_yr")
    tol = case["tolerance"]["legacy_raster_max_abs_diff_m"]

    worst = 0.0
    for year in case["legacy_years"]:
        path = folder / f"Buoyancy_State_{year}.tif"
        assert path.exists(), f"missing {path}"
        legacy = read_raster(path).values
        ours = criterion.legacy_margin(dam, f, level, dhdt, year - dam.base_year)
        both = np.isfinite(legacy) & np.isfinite(ours)
        assert both.sum() > 1_000_000, f"{path.name}: too few comparable cells"
        diff = float(np.max(np.abs(legacy[both] - ours[both])))
        worst = max(worst, diff)
        assert diff <= tol, f"{path.name}: max |legacy - ported| = {diff:.6g} m > {tol} m"

    # And the sign convention agrees: the owner's raster is negative where buoyant.
    assert worst <= tol


@pytest.mark.parametrize("name", sorted(CASES))
def test_legacy_raster_connectivity_gives_the_same_year(name: str) -> None:
    """Reading the breach year off the owner's own rasters reproduces the answer."""
    from gpbreach.io import read_raster

    case = CASES[name]
    folder = (repo_root() / case["legacy_rasters"]).resolve()
    if not folder.is_dir():
        pytest.skip(f"legacy rasters not available at {folder}")

    cfg, dam, _ = _result(case)
    criterion = FlotationConnectivity(connectivity=int(cfg.components.get("connectivity", 8)))
    connected_years = []
    for year in case["legacy_years"]:
        legacy = read_raster(folder / f"Buoyancy_State_{year}.tif").values
        # The owner's raster is negative where buoyant, so the sub-level set is
        # {legacy <= 0}; a breach is a connected path between the two seeds.
        field = np.where(np.isfinite(legacy), legacy, np.inf)
        from gpbreach.saddle import connected_at

        if connected_at(field, 0.0, dam.seed_lake, dam.seed_target, criterion.connectivity):
            connected_years.append(year)

    assert connected_years, "no legacy raster shows a connected pathway"
    assert min(connected_years) == case["breach_year"]
