"""Saddle solver correctness on fields whose answer is known by hand."""

from __future__ import annotations

import numpy as np
import pytest

from gpbreach.saddle import connected_at, locate_saddle, saddle_threshold


def test_single_gap_is_the_saddle() -> None:
    """Two basins separated by a wall with one low notch."""
    field = np.full((5, 7), 100.0)
    field[2, 0:3] = 0.0     # basin A
    field[2, 4:7] = 0.0     # basin B
    field[2, 3] = 42.0      # the notch
    assert saddle_threshold(field, (2, 0), (2, 6), connectivity=4) == pytest.approx(42.0)


def test_lowest_of_several_passes_wins() -> None:
    """Two routes between the basins; the minimax path takes the lower notch."""
    field = np.full((5, 7), 100.0)
    field[1:4, 0:3] = 0.0   # basin A, spanning rows 1-3 so both notches are reachable
    field[1:4, 4:7] = 0.0   # basin B
    field[1, 3] = 70.0      # high notch
    field[3, 3] = 30.0      # low notch -- this is the saddle
    assert saddle_threshold(field, (1, 0), (1, 6), connectivity=4) == pytest.approx(30.0)


def test_result_is_an_actual_cell_value() -> None:
    rng = np.random.default_rng(0)
    field = rng.normal(size=(40, 40))
    t = saddle_threshold(field, (0, 0), (39, 39))
    assert np.any(field == t)


def test_threshold_is_tight() -> None:
    """Connected exactly at the saddle, disconnected just below it."""
    rng = np.random.default_rng(1)
    field = rng.normal(size=(30, 30))
    t = saddle_threshold(field, (0, 0), (29, 29))
    assert connected_at(field, t, (0, 0), (29, 29))
    below = field[field < t].max()
    assert not connected_at(field, below, (0, 0), (29, 29))


def test_nan_cells_never_connect() -> None:
    field = np.full((5, 7), 100.0)
    field[2, 0:3] = 0.0
    field[2, 4:7] = 0.0
    field[:, 3] = np.nan    # an impassable gap
    assert saddle_threshold(field, (2, 0), (2, 6), connectivity=4) == np.inf


def test_eight_connectivity_can_cross_a_diagonal() -> None:
    field = np.full((3, 3), 100.0)
    field[0, 0] = field[1, 1] = field[2, 2] = 0.0
    assert saddle_threshold(field, (0, 0), (2, 2), connectivity=8) == pytest.approx(0.0)
    assert saddle_threshold(field, (0, 0), (2, 2), connectivity=4) == pytest.approx(100.0)


def test_locate_saddle_finds_the_notch() -> None:
    field = np.full((5, 7), 100.0)
    field[2, 0:3] = 0.0
    field[2, 4:7] = 0.0
    field[2, 3] = 42.0
    cells = locate_saddle(field, 42.0, (2, 0), (2, 6), connectivity=4)
    assert cells.tolist() == [[2, 3]]
