"""Sampling determinism and the N=1 identity."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.config import Config
from gpbreach.sampling import build_sample_table

CONFIG = "configs/baseline_breach1_lgp.yaml"


def test_n1_reproduces_nominal_values() -> None:
    cfg = Config.load(CONFIG)
    table = build_sample_table(cfg)
    assert len(table) == 1
    for name in cfg.parameters:
        assert table[name].iloc[0] == pytest.approx(cfg.nominal(name))


def test_same_seed_same_draws() -> None:
    cfg = Config.load(CONFIG)
    spec = dict(cfg.raw)
    spec["sampling"] = {"n_realizations": 500, "seed": 42}
    spec["parameters"] = dict(spec["parameters"])
    spec["parameters"]["lake_level_m"] = {
        "nominal": 110.0,
        "distribution": {"type": "normal", "sd": 5.0},
    }
    c = Config(spec["run_id"], spec["breach"], spec["dam"], spec["components"],
               spec["parameters"], spec["sampling"], spec["time"], cfg.path, spec)
    a, b = build_sample_table(c), build_sample_table(c)
    assert np.array_equal(a.to_numpy(), b.to_numpy())
    assert a["lake_level_m"].std() > 0
