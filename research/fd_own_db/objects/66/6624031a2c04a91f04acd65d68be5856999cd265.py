"""Tests for the fixed-margins G1 null and explicitly projected hypervolume."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pytest

from fd.g1_fixed_margins import (
    _bh, _curveball_once, _draw_cohort, _margins, _p_value, run_g1,
)
from fd.hypervolume import _area90, _grid_kernels, run_hypervolume
from fd.pipeline import REF_INPUT, SCOPES


def _write(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _fixture(root: Path, g0: Path) -> None:
    # Three distinct colour/form categories in a well-spread common PCoA.
    _write(root / REF_INPUT,
           ["profile_id","flower_primary_color","floral_form","PCoA1","PCoA2"],
           [{"profile_id": i, "flower_primary_color": c, "floral_form": f,
             "PCoA1": x, "PCoA2": y}
            for i, c, f, x, y in (
                (0, "red", "bell", -1, -1),
                (1, "blue", "bell", 1, -1),
                (2, "red", "tubular", -1, 1),
                (3, "blue", "tubular", 1, 1),
            )])
    for scope in SCOPES:
        base = root / "floral_space_20261106" / scope
        for island, ids in (("A", (0, 1)), ("B", (0, 2)),
                            ("C", (1, 3)), ("D", (2, 3))):
            pass
        _write(base / "island_profile_carriers.csv.gz",
               ["island_id","profile_id","n_species"],
               [{"island_id": island, "profile_id": p, "n_species": 50}
                for island, ids in (("A", (0,1)), ("B", (0,2)),
                                    ("C", (1,3)), ("D", (2,3))) for p in ids])
        _write(g0 / f"{scope}_islands.csv",
               ["island_id","analysis_regime","n_names","n_profiles"],
               [{"island_id": island, "analysis_regime": "northern_midlatitude",
                 "n_names": 100, "n_profiles": 2} for island in ("A","B","C","D")])
        _write(g0 / f"{scope}_pairs.csv",
               ["island_a","island_b","analysis_regime"],
               [{"island_a": a, "island_b": b,
                 "analysis_regime": "northern_midlatitude"} for a,b in
                (("A","B"),("A","C"),("A","D"),("B","C"),("B","D"),("C","D"))])


def test_curveball_preserves_row_column_margins() -> None:
    rng = np.random.default_rng(31)
    rows = [{0,1}, {1,2}, {0,2}, {1,3}, {2,3}]
    original = _margins(rows)
    for _ in range(1800):
        _curveball_once(rows, rng)
    assert _margins(rows) == original


def test_pair_selection_is_deterministic_and_observation_blind() -> None:
    pairs = [("A","B"),("A","C"),("B","C"),("A","D"),("B","D")]
    x = _draw_cohort(pairs, seed=141, local_n=3, pair_n=2)
    y = _draw_cohort(pairs, seed=141, local_n=3, pair_n=2)
    assert x == y
    assert len(x[1]) == 2 and len(x[2]) == 3


def test_adjustment_and_inclusive_ties() -> None:
    assert np.allclose(_bh([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])
    assert _p_value(0.1, np.array([0.1]*15)) == 1.0


def test_gaussian_single_center_hdr_area_near_analytic() -> None:
    coords = np.array([[0.0, 0.0], [2.0, 2.0]])
    kernel, area = _grid_kernels(coords, 0.35, 88)
    numerical, mass = _area90(kernel[:1], area)
    # Two-dimensional Gaussian 90% HDR area = pi*h^2 * 2ln(10).
    analytical = np.pi * 0.35**2 * (2.0 * np.log(10.0))
    assert numerical[0] == pytest.approx(analytical, rel=0.10)
    assert mass[0] > 0.994


def test_frozen_g0_to_g1_and_projected_hv_contract(tmp_path: Path) -> None:
    inp, g0 = tmp_path / "input", tmp_path / "g0"
    _fixture(inp, g0)
    g1 = run_g1(inp, g0, tmp_path / "g1", local_n=3, pair_n=3,
                draws_per_chain=12, burn_trades_per_island=4,
                thin_trades_per_island=1)
    assert g1["observed_test_slots"] == 32
    assert (tmp_path / "g1" / "g1_draws.csv").is_file()
    assert g1["conditioned_on"][0].startswith("island profile")
    hv = run_hypervolume(inp, g0, tmp_path / "hv", grid_n=60)
    assert hv["reference_profiles"] == 4
    assert hv["n_islands"] == 16
    assert hv["two_axis_projection_loss"] > 0.8
    with (tmp_path / "hv" / "all_hypervolumes.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert rows[0]["observed_reference_occupancy_fraction"] == "0.5"
    assert float(rows[0]["kde90_area_pc12_equal_h035"]) > 0
