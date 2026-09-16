"""Realization sampling.

Deterministic runs are Monte Carlo with ``n_realizations: 1`` -- there is no
separate deterministic code path (CLAUDE.md Section 3). All randomness flows
from ``sampling.seed``; nothing here calls an unseeded RNG.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .config import Config


def _draw(rng: np.random.Generator, spec: dict[str, Any], n: int) -> np.ndarray:
    """Draw ``n`` values for one parameter.

    With no ``distribution`` key the nominal value is repeated, which is what
    makes an ``n_realizations: 1`` run reproduce the deterministic baseline.
    """
    nominal = float(spec["nominal"])
    dist = spec.get("distribution")
    if dist is None:
        return np.full(n, nominal)
    kind = dist["type"]
    if kind == "normal":
        values = rng.normal(dist.get("mean", nominal), dist["sd"], n)
    elif kind == "uniform":
        values = rng.uniform(dist["low"], dist["high"], n)
    elif kind == "triangular":
        values = rng.triangular(dist["low"], dist.get("mode", nominal), dist["high"], n)
    else:
        raise ValueError(f"unknown distribution type {kind!r}")
    lo, hi = dist.get("clip_low"), dist.get("clip_high")
    if lo is not None or hi is not None:
        values = np.clip(
            values, lo if lo is not None else -np.inf, hi if hi is not None else np.inf
        )
    return values


def build_sample_table(cfg: Config) -> pd.DataFrame:
    """One row per realization, one column per sampled parameter."""
    n = int(cfg.sampling["n_realizations"])
    rng = np.random.default_rng(int(cfg.sampling["seed"]))
    table = pd.DataFrame(
        {name: _draw(rng, spec, n) for name, spec in cfg.parameters.items()},
        index=pd.RangeIndex(n, name="realization"),
    )
    table.attrs["seed"] = int(cfg.sampling["seed"])
    return table
