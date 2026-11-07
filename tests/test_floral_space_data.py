"""Protect the scientific input contract, independently of ordination code."""

from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path

import pandas as pd
import pytest
import yaml


TRAITS = (
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "tube_depth_class",
)
COLUMNS = [
    "accepted_species",
    "trait_name",
    "normalized_value",
    "quality",
    "source_groups",
    "source_lineages",
    "resolution_status",
    "evidence_scope",
]


def _row(species, trait, value, *, quality="medium", lineage="source:one", **extra):
    return {
        "accepted_species": species,
        "trait_name": trait,
        "normalized_value": value,
        "quality": quality,
        "source_groups": "reviewed_source",
        "source_lineages": lineage,
        "resolution_status": "resolved",
        "evidence_scope": "validated_low" if quality == "low" else "direct",
        **extra,
    }


def _inputs(tmp_path, raw, snapshot_direct, snapshot_low=()):
    paths = []
    for name, rows in [
        ("raw.csv.gz", raw),
        ("snapshot_all.csv.gz", [*snapshot_direct, *snapshot_low]),
        ("snapshot_direct.csv.gz", snapshot_direct),
    ]:
        path = tmp_path / name
        pd.DataFrame(rows, columns=COLUMNS).to_csv(path, index=False)
        paths.append(path)
    ontology = {
        "traits": {
            "flower_primary_color": {
                "allowed_values": [
                    "white", "red_pink", "blue_purple", "yellow_orange",
                    "other_described", "multicolored_variable", "unresolved",
                ]
            },
            "floral_form": {
                "allowed_values": [
                    "open_radial", "tubular", "bell_campanulate",
                    "other_described", "unresolved",
                ]
            },
            "floral_symmetry": {
                "allowed_values": ["actinomorphic", "zygomorphic", "unresolved"]
            },
            "tube_depth_class": {
                "allowed_values": [
                    "absent_or_open", "shallow", "intermediate", "deep", "unresolved"
                ]
            },
        }
    }
    ontology_path = tmp_path / "ontology.yml"
    ontology_path.write_text(yaml.safe_dump(ontology), encoding="utf-8")
    return [*paths, ontology_path]


def _load(*args):
    # An assertion, rather than an import-collection error, records the initial RED.
    spec = importlib.util.find_spec("island_v2.floral_space_data")
    assert spec is not None, "floral-space evidence adapter has not been implemented"
    return importlib.import_module("island_v2.floral_space_data").load_profiles(*args)


def test_keeps_multistate_sets_in_one_species_profile(tmp_path: Path):
    # Cartesian expansion would invent four flowers and multiply the species count.
    rows = [
        _row("Alpha one", TRAITS[0], "white|red_pink|white"),
        _row("Alpha one", TRAITS[1], "tubular|open_radial"),
    ]
    all_profiles, direct, _ = _load(*_inputs(tmp_path, rows, rows))
    assert len(all_profiles) == len(direct) == 1
    assert direct.iloc[0][TRAITS[0]] == "red_pink|white"
    assert direct.iloc[0][TRAITS[1]] == "open_radial|tubular"


def test_retains_snapshot_direct_and_restores_missing_raw_traits(tmp_path: Path):
    # Axis-level quality selection must not discard an independent Medium form.
    snapshot = [_row("Alpha one", TRAITS[0], "white", quality="high")]
    raw = [*snapshot, _row("Alpha one", TRAITS[1], "tubular")]
    low = [
        _row("Alpha one", TRAITS[1], "open_radial", quality="low"),
        _row("Alpha one", TRAITS[2], "actinomorphic", quality="low"),
    ]
    all_profiles, direct, audit = _load(*_inputs(tmp_path, raw, snapshot, low))
    assert all_profiles.iloc[0][TRAITS[1]] == "tubular"
    assert direct.iloc[0][TRAITS[1]] == "tubular"
    assert all_profiles.iloc[0][TRAITS[2]] == "actinomorphic"
    assert pd.isna(direct.iloc[0][TRAITS[2]])
    assert direct.iloc[0][f"{TRAITS[1]}__quality"] == "medium"
    assert direct.iloc[0][f"{TRAITS[1]}__record_source"] == "raw_direct_recovery"
    assert all_profiles.iloc[0][f"{TRAITS[2]}__evidence_scope"] == "validated_low"
    assert audit["selection"]["raw_direct_cells_restored"] == 1
    assert audit["selection"]["low_cells_blocked_by_direct_key"] == 1


def test_preserves_existing_broader_colour_evidence_without_unioning(tmp_path: Path):
    # Taking the raw row would erase a documented colour and its independent source.
    raw = [_row("Alpha one", TRAITS[0], "white", lineage="source:one")]
    snapshot = [
        _row("Alpha one", TRAITS[0], "white|yellow_orange", lineage="source:one|source:two")
    ]
    _, direct, audit = _load(*_inputs(tmp_path, raw, snapshot))
    assert direct.iloc[0][TRAITS[0]] == "white|yellow_orange"
    assert direct.iloc[0][f"{TRAITS[0]}__source_lineages"] == "source:one|source:two"
    assert audit["comparison"]["state_disagreements"] == 1
    assert audit["comparison"]["snapshot_state_and_lineage_supersets"] == 1


def test_unexplained_direct_disagreement_fails_closed(tmp_path: Path):
    # Unioning contradictory signatures would manufacture a variable phenotype.
    raw = [_row("Alpha one", TRAITS[0], "white", lineage="source:one")]
    snapshot = [_row("Alpha one", TRAITS[0], "red_pink", lineage="source:two")]
    with pytest.raises(ValueError, match="unexplained direct disagreement"):
        _load(*_inputs(tmp_path, raw, snapshot))


@pytest.mark.parametrize("value", ["", "unresolved", "other_described", "white|multicolored_variable"])
def test_unspecified_cells_are_missing_and_cannot_be_replaced_by_low(tmp_path: Path, value):
    # Unknown/generic colour descriptions must not create artificial floral types.
    raw = [_row("Alpha one", TRAITS[0], value)]
    low = [_row("Alpha one", TRAITS[0], "blue_purple", quality="low")]
    all_profiles, direct, audit = _load(*_inputs(tmp_path, raw, [], low))
    assert pd.isna(all_profiles.iloc[0][TRAITS[0]])
    assert pd.isna(direct.iloc[0][TRAITS[0]])
    assert audit["selection"]["low_cells_blocked_by_direct_key"] == 1


def test_does_not_treat_an_out_of_ontology_label_as_a_flower_form(tmp_path: Path):
    # Known frozen Low label permutations must not enter the ordination as states.
    low = [_row("Alpha one", TRAITS[1], "actinomorphic", quality="low")]
    all_profiles, _, audit = _load(*_inputs(tmp_path, [], [], low))
    assert pd.isna(all_profiles.iloc[0][TRAITS[1]])
    assert audit["input_validation"]["snapshot_all"]["invalid_ontology_cells"] == 1


def test_rejects_duplicate_species_trait_keys(tmp_path: Path):
    # A many-to-many pivot/join can silently multiply carriers or choose a record.
    row = _row("Alpha one", TRAITS[0], "white")
    with pytest.raises(ValueError, match="duplicate.*accepted_species.*trait_name"):
        _load(*_inputs(tmp_path, [row, row], [row]))


def test_excluded_direct_record_cannot_be_promoted(tmp_path: Path):
    # Excluded source conflicts remain unresolved, even if a normalized value exists.
    raw = [_row("Alpha one", TRAITS[0], "white", resolution_status="excluded")]
    all_profiles, direct, audit = _load(*_inputs(tmp_path, raw, []))
    assert pd.isna(all_profiles.iloc[0][TRAITS[0]])
    assert pd.isna(direct.iloc[0][TRAITS[0]])
    assert audit["input_validation"]["raw_direct"]["nonresolved_cells"] == 1


def test_joint_coverage_requires_every_requested_trait(tmp_path: Path):
    # Marginal completeness must never be reported as simultaneous completeness.
    snapshot = [
        _row("Alpha one", TRAITS[0], "white"),
        _row("Alpha one", TRAITS[1], "tubular"),
        _row("Beta two", TRAITS[0], "red_pink"),
    ]
    _, _, audit = _load(*_inputs(tmp_path, snapshot, snapshot))
    assert audit["coverage"]["restored_direct"]["colour_form"] == 1
    assert audit["coverage"]["restored_direct"]["colour_form_symmetry"] == 0
    assert audit["coverage"]["restored_direct"]["colour_form_symmetry_depth"] == 0
