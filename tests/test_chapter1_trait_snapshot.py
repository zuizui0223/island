import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_trait_snapshot import (
    AXES,
    build_snapshot,
    materialize_long_ledger,
    transition_audit,
    validate_species_axis,
)


def _frame() -> pd.DataFrame:
    rows = []
    for species in ("Alpha one", "Beta two"):
        for axis in AXES:
            rows.append(
                {
                    "accepted_species": species,
                    "axis": axis,
                    "trait_composition": "",
                    "trait_names": "",
                    "source_groups": "",
                    "source_lineages": "",
                    "quality": "",
                }
            )
    frame = pd.DataFrame(rows)

    def set_cell(species, axis, composition, quality, groups="source_a", lineages="line_a"):
        mask = frame["accepted_species"].eq(species) & frame["axis"].eq(axis)
        frame.loc[mask, "trait_composition"] = composition
        frame.loc[mask, "trait_names"] = "|".join(
            item.split("=", 1)[0] for item in composition.split("|")
        )
        frame.loc[mask, "source_groups"] = groups
        frame.loc[mask, "source_lineages"] = lineages
        frame.loc[mask, "quality"] = quality

    set_cell(
        "Alpha one",
        "flower_colour",
        'flower_primary_color=["red_pink"]',
        "high",
    )
    set_cell(
        "Alpha one",
        "floral_structural_complexity",
        'floral_form=["tubular"]|tube_depth_class=["deep"]',
        "medium",
    )
    set_cell(
        "Alpha one",
        "reproductive_assurance",
        'self_incompatibility=["SI"]',
        "medium",
    )
    set_cell(
        "Beta two",
        "reproductive_assurance",
        'self_incompatibility=["SC"]',
        "low",
    )
    return frame


def test_snapshot_materializes_direct_and_all_ledgers(tmp_path: Path):
    source = tmp_path / "coverage.csv.gz"
    _frame().to_csv(source, index=False, compression="gzip")
    out = tmp_path / "snapshot"
    manifest = build_snapshot(
        species_axis_csv=source,
        output_dir=out,
        expected_species=2,
        source_run_id="123",
        source_artifact_name="wave-test",
    )

    all_ledger = pd.read_csv(out / "chapter1_trait_ledger_all_analysis.csv.gz")
    direct = pd.read_csv(out / "chapter1_trait_ledger_direct_only.csv.gz")
    assert len(all_ledger) == 5
    assert len(direct) == 4
    assert set(direct["quality"]) == {"high", "medium"}
    assert "low" in set(all_ledger["quality"])
    si = direct.loc[
        direct["trait_name"].eq("self_incompatibility")
        & direct["accepted_species"].eq("Alpha one"),
        "normalized_value",
    ].item()
    assert si == "SI"
    assert manifest["coverage"]["global_species_axis_fill_fraction"] == pytest.approx(4 / 6)
    assert manifest["analysis_eligibility_rule"].startswith("response-specific")
    saved = json.loads((out / "chapter1_trait_snapshot_manifest.json").read_text())
    assert saved["source_workflow_run_id"] == "123"


def test_validated_low_is_not_direct():
    frame = validate_species_axis(_frame(), expected_species=2)
    all_ledger = materialize_long_ledger(frame, direct_only=False)
    direct = materialize_long_ledger(frame, direct_only=True)
    beta = all_ledger.loc[all_ledger["accepted_species"].eq("Beta two")]
    assert len(beta) == 1
    assert beta.iloc[0]["evidence_scope"] == "validated_low"
    assert "Beta two" not in set(direct["accepted_species"])


def test_reported_low_without_composition_is_retained_but_not_materialized(tmp_path: Path):
    frame = _frame()
    mask = frame["accepted_species"].eq("Beta two") & frame["axis"].eq("flower_colour")
    frame.loc[mask, "quality"] = "low"
    source = tmp_path / "coverage.csv.gz"
    frame.to_csv(source, index=False, compression="gzip")
    out = tmp_path / "snapshot"
    manifest = build_snapshot(
        species_axis_csv=source,
        output_dir=out,
        expected_species=2,
        source_run_id="wave",
        source_artifact_name="wave",
    )
    ledger = pd.read_csv(out / "chapter1_trait_ledger_all_analysis.csv.gz")
    assert not (
        ledger["accepted_species"].eq("Beta two")
        & ledger["axis"].eq("flower_colour")
    ).any()
    assert manifest["coverage"]["source_reported_low_without_trait_composition"] == 1
    colour = manifest["coverage"]["by_axis"]["flower_colour"]
    assert colour["source_reported_filled_species"] == 2
    assert colour["resolved_species"] == 1


def test_direct_quality_without_composition_is_rejected():
    frame = _frame()
    mask = frame["accepted_species"].eq("Beta two") & frame["axis"].eq("flower_colour")
    frame.loc[mask, "quality"] = "medium"
    with pytest.raises(ValueError, match="direct-quality"):
        validate_species_axis(frame, expected_species=2)


def test_cross_axis_trait_is_rejected():
    frame = _frame()
    mask = frame["accepted_species"].eq("Alpha one") & frame["axis"].eq("flower_colour")
    frame.loc[mask, "trait_composition"] = 'floral_form=["tubular"]'
    frame.loc[mask, "trait_names"] = "floral_form"
    with pytest.raises(ValueError, match="outside the declared axis"):
        validate_species_axis(frame, expected_species=2)


def test_transition_audit_marks_new_resolution_and_revision():
    previous = validate_species_axis(_frame(), expected_species=2)
    current_raw = _frame()
    new_mask = current_raw["accepted_species"].eq("Beta two") & current_raw["axis"].eq(
        "flower_colour"
    )
    current_raw.loc[new_mask, "trait_composition"] = 'flower_primary_color=["white_plain"]'
    current_raw.loc[new_mask, "trait_names"] = "flower_primary_color"
    current_raw.loc[new_mask, "source_groups"] = "source_b"
    current_raw.loc[new_mask, "source_lineages"] = "line_b"
    current_raw.loc[new_mask, "quality"] = "medium"
    revision_mask = current_raw["accepted_species"].eq("Alpha one") & current_raw["axis"].eq(
        "flower_colour"
    )
    current_raw.loc[revision_mask, "trait_composition"] = 'flower_primary_color=["blue_purple"]'
    current = validate_species_axis(current_raw, expected_species=2)
    audit, counts = transition_audit(previous, current)
    assert counts["newly_resolved"] == 1
    assert counts["value_revision"] == 1
    assert len(audit) == 6


def test_revision_fails_closed_without_explicit_override(tmp_path: Path):
    previous = tmp_path / "previous.csv.gz"
    current = tmp_path / "current.csv.gz"
    _frame().to_csv(previous, index=False, compression="gzip")
    revised = _frame()
    mask = revised["accepted_species"].eq("Alpha one") & revised["axis"].eq("flower_colour")
    revised.loc[mask, "trait_composition"] = 'flower_primary_color=["blue_purple"]'
    revised.to_csv(current, index=False, compression="gzip")

    with pytest.raises(ValueError, match="allow-evidence-revision"):
        build_snapshot(
            species_axis_csv=current,
            previous_species_axis_csv=previous,
            output_dir=tmp_path / "out",
            expected_species=2,
            source_run_id="new",
            source_artifact_name="new-wave",
        )
