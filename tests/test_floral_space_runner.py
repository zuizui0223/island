"""Reproduce a tiny four-scope analysis from frozen ledger-shaped inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest


REPOSITORY = Path(__file__).resolve().parents[1]
RUNNER = REPOSITORY / "scripts/run_floral_space_exploration.py"
TRAITS = ["flower_primary_color", "floral_form", "floral_symmetry", "tube_depth_class"]
NAMES = ["Alpha a", "Alpha b", "Beta c", "Gamma d", "Delta e", "Alpha", "Alpha sp."]


def _row(species, trait, state, quality="medium"):
    return {
        "accepted_species": species,
        "trait_name": trait,
        "normalized_value": state,
        "quality": quality,
        "source_groups": "fixture",
        "source_lineages": "fixture:one",
        "resolution_status": "resolved",
    }


def _fixture(tmp_path):
    root = tmp_path / "inputs"
    snapshot = root / "progressive/snapshot"
    source = root / "trait_source"
    flora_directory = root / "progressive/fixed/canonical/input"
    for directory in (snapshot, source, flora_directory):
        directory.mkdir(parents=True)
    rows = []
    for name, colour, form in [
        ("Alpha a", "white", "open_radial"),
        ("Alpha b", "red_pink", "tubular"),
        ("Beta c", "white|red_pink", "bell_campanulate"),
        ("Alpha", "white", "spurred"),
        ("Alpha sp.", "white", "spurred"),
    ]:
        rows.extend([_row(name, TRAITS[0], colour), _row(name, TRAITS[1], form)])
    rows.extend([
        _row("Alpha a", TRAITS[2], "actinomorphic"),
        _row("Alpha a", TRAITS[3], "absent_or_open"),
    ])
    direct = [row for row in rows
              if (row["accepted_species"], row["trait_name"]) != ("Alpha b", TRAITS[1])]
    low = [
        _row("Gamma d", TRAITS[0], "yellow_orange", "low"),
        _row("Gamma d", TRAITS[1], "open_radial", "low"),
    ]
    pd.DataFrame(rows).to_csv(source / "direct_species_trait_ledger.csv.gz", index=False)
    pd.DataFrame(direct).to_csv(snapshot / "chapter1_trait_ledger_direct_only.csv.gz", index=False)
    pd.DataFrame(direct + low).to_csv(snapshot / "chapter1_trait_ledger_all_analysis.csv.gz", index=False)
    # This is an axis-grain table: fourteen rows must yield seven name units.
    pd.DataFrame([
        {"accepted_species": name, "axis": axis}
        for name in NAMES for axis in ("floral_signal", "floral_architecture")
    ]).to_csv(source / "species_axis_coverage.csv.gz", index=False)
    pd.DataFrame([
        {"island_id": island, "accepted_species": name}
        for island, names in [
            ("one", ["Alpha a", "Alpha b", "Beta c", "Delta e", "Alpha", "Alpha sp."]),
            ("two", ["Alpha a", "Gamma d", "Delta e"]),
        ]
        for name in names
    ]).to_csv(flora_directory / "chapter1_status_flora.csv.gz", index=False)
    geography = tmp_path / "geography.csv"
    pd.DataFrame({"island_id": ["one", "two"], "analysis_regime": ["r", "r"]}).to_csv(
        geography, index=False
    )
    return root, geography


def _run(root, geography, output, *, expected=7):
    assert RUNNER.is_file(), "reproducible floral-space runner has not been implemented"
    return subprocess.run([
        sys.executable, str(RUNNER), "--input-root", str(root),
        "--output", str(output), "--geography", str(geography),
        "--expected-species", str(expected), "--min-species", "1",
        "--rarefaction-size", "1", "--pair-count", "1", "--seed", "17",
    ], cwd=REPOSITORY, capture_output=True, text=True)


def test_runner_preserves_population_and_uses_one_reference_space(tmp_path):
    root, geography = _fixture(tmp_path)
    output = tmp_path / "results"
    result = _run(root, geography, output)
    assert result.returncode == 0, result.stderr

    prepared = {}
    for scope in ("all", "direct", "all_single_state", "direct_single_state"):
        frame = pd.read_csv(root / f"{scope}_profiles.csv.gz")
        assert frame.columns.tolist() == [
            "accepted_species", *TRAITS, *[f"{trait}__quality" for trait in TRAITS]
        ]
        assert len(frame) == 7
        assert set(frame.accepted_species) == set(NAMES)
        assert not frame.accepted_species.duplicated().any()
        prepared[scope] = frame.set_index("accepted_species")
        assert prepared[scope].loc["Delta e", TRAITS].isna().all()
    assert prepared["all"].loc["Beta c", TRAITS[0]] == "red_pink|white"
    assert prepared["all"].loc["Alpha b", TRAITS[1]] == "tubular"
    assert prepared["direct"].loc["Gamma d", TRAITS].isna().all()
    assert prepared["all"].loc["Gamma d", f"{TRAITS[0]}__quality"] == "low"
    for scope in ("all_single_state", "direct_single_state"):
        assert prepared[scope].loc["Beta c", TRAITS[:2]].isna().all()
        assert prepared[scope].loc["Beta c", f"{TRAITS[0]}__quality"] == "medium"

    atlas = None
    for scope, complete in {"all": 4, "direct": 3, "all_single_state": 3,
                            "direct_single_state": 2}.items():
        summary = json.loads((output / scope / "summary.json").read_text())
        assert summary["upstream_population_name_units"] == 7
        assert summary["n_profile_population_species"] == 5
        assert summary["n_excluded_uninomial_names"] == 1
        assert summary["n_excluded_unresolved_names"] == 1
        assert summary["n_complete_species"] == complete
        assert summary["n_profiles"] == complete
        assert summary["n_reference_profiles"] == 4
        coordinates = pd.read_csv(output / scope / "profile_coordinates.csv")
        coordinates = coordinates.drop(columns="n_species_scope")
        if atlas is None:
            atlas = coordinates
        else:
            pd.testing.assert_frame_equal(coordinates, atlas)
        inputs = json.loads((output / scope / "input_manifest.json").read_text())
        assert inputs["space_profiles"]["sha256"] == hashlib.sha256(
            (root / "all_profiles.csv.gz").read_bytes()
        ).hexdigest()

    audit = json.loads((output / "data_quality_audit.json").read_text())
    assert audit["population_species"] == 7
    assert audit["selection"]["raw_direct_cells_restored"] == 1
    coverage = pd.read_csv(output / "joint_coverage.csv")
    assert len(coverage) == 12
    assert set(coverage.population_species) == {7}
    primary = coverage.loc[coverage.trait_set.eq("colour_form")].set_index("scope")
    assert primary.n_species.to_dict() == {
        "snapshot_all": 5, "snapshot_direct": 4, "restored_all": 6, "restored_direct": 5
    }
    assert primary.loc["restored_all", "fraction"] == pytest.approx(6 / 7)
    provenance = json.loads((output / "provenance_manifest.json").read_text())
    assert provenance["population_species"] == 7
    raw = provenance["sources"]["raw_direct"]
    assert raw["sha256"] == hashlib.sha256(
        (root / "trait_source/direct_species_trait_ledger.csv.gz").read_bytes()
    ).hexdigest()
    assert raw["matches_reference_sha256"] is False
    assert provenance["upstream_reference_artifacts"]["trait_source"]["run_id"] == 34191508045
    assert provenance["upstream_reference_artifacts"]["progressive"]["run_id"] == 34232450884
    assert provenance["software"]["python"]
    assert {"numpy", "pandas", "scipy", "PyYAML"} <= set(provenance["software"]["libraries"])
    assert provenance["status"] == "complete"


def test_runner_rejects_wrong_population_before_writing_outputs(tmp_path):
    root, geography = _fixture(tmp_path)
    output = tmp_path / "results"
    result = _run(root, geography, output, expected=8)
    assert result.returncode != 0
    assert "Expected 8" in result.stderr and "found 7" in result.stderr
    assert not output.exists()
    assert not (root / "all_profiles.csv.gz").exists()
