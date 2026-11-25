from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from fd.core import ProfileSpace, distance, load_carriers, load_profiles
from fd.pipeline import PAIR_INPUT, REF_INPUT, run


PROFILES = {
    0: (frozenset({"red"}), frozenset({"bell"})),
    1: (frozenset({"blue"}), frozenset({"bell"})),
    2: (frozenset({"red"}), frozenset({"tubular"})),
    3: (frozenset({"blue"}), frozenset({"tubular"})),
}


def write_csv(path: Path, headers: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def test_original_distance_and_multistate_non_cartesian() -> None:
    space = ProfileSpace(PROFILES)
    assert distance(PROFILES[0], PROFILES[1]) == pytest.approx(0.5)
    assert distance(PROFILES[0], PROFILES[2]) == pytest.approx(0.5)
    assert distance(PROFILES[0], PROFILES[3]) == pytest.approx(1)
    multistate = (frozenset({"blue", "red"}), frozenset({"bell"}))
    assert distance(multistate, PROFILES[0]) == pytest.approx(0.25)
    assert len(ProfileSpace({0: PROFILES[0], 1: multistate}).ids) == 2


def test_rao_q_and_profile_equal_weight() -> None:
    space = ProfileSpace(PROFILES)
    row = space.island({0: 1, 1: 1})
    assert row["rao_q"] == pytest.approx(0.25)
    assert row["rao_q_without_replacement"] == pytest.approx(0.5)
    assert row["rao_q_profile_equal"] == pytest.approx(0.25)
    assert space.island({0: 7})["rao_q"] == pytest.approx(0)
    assert space.island({0: 7})["rao_q_without_replacement"] == pytest.approx(0)
    with pytest.raises(ValueError, match="Unknown"):
        space.island({999: 1})
    with pytest.raises(ValueError, match="Empty"):
        space.island({})


def test_pair_functional_dissimilarity_distinct_from_cross_distance() -> None:
    space = ProfileSpace(PROFILES)
    a = {0: 1, 1: 1}
    b = {2: 1, 3: 1}
    score = space.pair(a, b)
    assert score["profile_jaccard"] == 0
    assert score["functional_nn_dissimilarity"] == pytest.approx(0.5)
    assert score["mean_cross_profile_distance"] == pytest.approx(0.75)
    assert space.pair(a, a)["functional_nn_dissimilarity"] == pytest.approx(0)
    assert space.pair(a, a)["mean_cross_profile_distance"] == pytest.approx(0.25)


def fixture_inputs(root: Path) -> None:
    write_csv(root / REF_INPUT,
              ["profile_id", "flower_primary_color", "floral_form", "PCoA1"],
              [{"profile_id": i, "flower_primary_color": next(iter(p[0])),
                "floral_form": next(iter(p[1])), "PCoA1": 9000 + i}
               for i, p in PROFILES.items()])
    write_csv(root / PAIR_INPUT,
              ["encoding", "evidence_scope", "analysis_regime", "island_a", "island_b"],
              [{"encoding": "single_state" if "single" in s else "full_state_sets",
                "evidence_scope": s, "analysis_regime": "northern_midlatitude",
                "island_a": "A", "island_b": "B"}
               for s in ("all", "direct", "all_single_state", "direct_single_state")])
    for scope in ("all", "direct", "all_single_state", "direct_single_state"):
        base = root / "floral_space_20261106" / scope
        write_csv(base / "island_metrics.csv",
                  ["island_id", "analysis_regime", "n_species", "n_profiles",
                   "joint_coverage", "area_km2", "distance_to_continent_km"],
                  [{"island_id": x, "analysis_regime": "northern_midlatitude",
                    "n_species": 2, "n_profiles": 2, "joint_coverage": 0.75,
                    "area_km2": 50, "distance_to_continent_km": 100}
                   for x in ("A", "B")])
        write_csv(base / "island_profile_carriers.csv.gz",
                  ["island_id", "profile_id", "n_species"],
                  [{"island_id": x, "profile_id": i, "n_species": 1}
                   for x, profiles in (("A", (0, 1)), ("B", (2, 3)))
                   for i in profiles])


def test_pipeline_public_aggregate_end_to_end(tmp_path: Path) -> None:
    root = tmp_path / "inputs"
    fixture_inputs(root)
    manifest = run(root, tmp_path / "outputs", min_names=2)
    assert manifest["reference_profiles"] == 4
    assert manifest["summaries"]["all"]["n_islands_eligible"] == 2
    assert manifest["summaries"]["direct"]["n_pairs"] == 1
    with (tmp_path / "outputs" / "all_pairs.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert float(rows[0]["functional_nn_dissimilarity"]) == pytest.approx(0.5)
    assert float(rows[0]["mean_cross_profile_distance"]) == pytest.approx(0.75)
    assert manifest["method"]["genus_conditioning"] is False
    assert manifest["method"]["shared_taxa_removed"] is False


def test_pipeline_support_gate_excludes_observation_not_zeros(tmp_path: Path) -> None:
    root = tmp_path / "inputs"
    fixture_inputs(root)
    result = run(root, tmp_path / "outputs", min_names=3)
    assert result["summaries"]["all"]["n_islands_eligible"] == 0
    assert result["summaries"]["all"]["n_pairs"] == 0


def test_input_validation(tmp_path: Path) -> None:
    profile_path = tmp_path / "profiles.csv"
    write_csv(profile_path, ["profile_id", "flower_primary_color", "floral_form"],
              [{"profile_id": 0, "flower_primary_color": "red|red", "floral_form": "bell"}])
    with pytest.raises(ValueError, match="Duplicated state"):
        load_profiles(profile_path)
    carriers_path = tmp_path / "carriers.csv.gz"
    write_csv(carriers_path, ["island_id", "profile_id", "n_species"],
              [{"island_id": "A", "profile_id": 0, "n_species": 1},
               {"island_id": "A", "profile_id": 0, "n_species": 1}])
    with pytest.raises(ValueError, match="Duplicate island/profile"):
        load_carriers(carriers_path)
