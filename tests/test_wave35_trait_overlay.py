from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

import island_v2.wave35_trait_overlay as overlay_module
from island_v2.wave35_trait_overlay import AXES, build_wave35_overlay


@pytest.fixture(autouse=True)
def small_universe(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(overlay_module, "EXPECTED_SPECIES", 4)
    monkeypatch.setattr(overlay_module, "EXPECTED_CELLS", 12)


def _write_inputs(root: Path, *, wrong_low_axis: bool = False) -> tuple[Path, Path]:
    species = [f"Genus species{i}" for i in range(4)]
    baseline = pd.DataFrame(
        [(name, axis, "", "", "", "", "") for name in species for axis in AXES],
        columns=[
            "accepted_species",
            "axis",
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ],
    )
    baseline.loc[
        (baseline.accepted_species == species[0]) & (baseline.axis == "flower_colour"),
        ["trait_composition", "trait_names", "source_groups", "source_lineages", "quality"],
    ] = [
        'flower_primary_color=["red_pink"]',
        "flower_primary_color",
        "prior",
        "doi:prior",
        "high",
    ]
    baseline.loc[
        (baseline.accepted_species == species[1])
        & (baseline.axis == "reproductive_assurance"),
        ["trait_composition", "trait_names", "source_groups", "source_lineages", "quality"],
    ] = [
        'self_incompatibility=["SI"]',
        "self_incompatibility",
        "prior_low",
        "validated-low:prior",
        "low",
    ]
    baseline_path = root / "baseline.csv.gz"
    baseline.to_csv(baseline_path, index=False)

    input_dir = root / "inputs"
    input_dir.mkdir()
    direct = pd.DataFrame(
        [
            {
                "accepted_species": species[0],
                "axis": "flower_colour",
                "trait_name": "flower_primary_color",
                "resolution_status": "resolved",
                "quality": "high",
                "state_set": '["red_pink", "white"]',
                "source_groups": "floraweb_biolflor_bulk",
                "source_lineages": "floraweb:one",
            },
            {
                "accepted_species": species[1],
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "quality": "high",
                "state_set": '["SC"]',
                "source_groups": "glopl_2025_trait_compilation",
                "source_lineages": "glopl:one",
            },
            {
                "accepted_species": species[2],
                "axis": "floral_structural_complexity",
                "trait_name": "floral_symmetry",
                "resolution_status": "resolved",
                "quality": "medium",
                "state_set": '["actinomorphic"]',
                "source_groups": "floraweb_biolflor_bulk",
                "source_lineages": "floraweb:two",
            },
        ]
    )
    direct.to_csv(input_dir / "wave35_resolved_direct_species_trait.csv.gz", index=False)

    low_axis = "flower_colour" if wrong_low_axis else "reproductive_assurance"
    low = pd.DataFrame(
        [
            {
                "accepted_species": species[3],
                "genus": "Genus",
                "axis": low_axis,
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "quality": "low",
                "family_inference_used": "False",
                "global_fallback_used": "False",
                "source_lineage": 'validated-low:Genus:self_incompatibility:["SC"]',
            }
        ]
    )
    low.to_csv(input_dir / "wave35_candidate_validated_low_species_trait.csv.gz", index=False)
    rules = pd.DataFrame(
        [
            {
                "genus": "Genus",
                "axis": low_axis,
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "eligible": "True",
                "n_direct_species": "3",
                "species_loo_accuracy": "1.0",
                "lineage_loo_accuracy": "1.0",
            }
        ]
    )
    rules.to_csv(input_dir / "wave35_provider_touched_new_rule_audit.csv.gz", index=False)
    pd.DataFrame([{"comparison_status": "invalidated_new_counterexample_or_validation"}]).to_csv(
        input_dir / "wave35_old_low_comparison.csv.gz", index=False
    )
    hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in input_dir.iterdir()
        if path.is_file()
    }
    (input_dir / "wave35_source_manifest.json").write_text(
        json.dumps(
            {
                "contract": "wave35_reproduction_morphology_source_manifest_v1",
                "file_sha256": hashes,
            }
        ),
        encoding="utf-8",
    )
    return baseline_path, input_dir


def test_adds_direct_and_trait_specific_low_without_loss(tmp_path: Path) -> None:
    baseline, inputs = _write_inputs(tmp_path)
    output = tmp_path / "output"
    summary = build_wave35_overlay(
        baseline_csv=baseline,
        input_dir=inputs,
        output_dir=output,
    )

    assert summary["delta"]["gross_gain_species_axis"] == 2
    assert summary["delta"]["loss_species_axis"] == 0
    assert summary["delta"]["action_counts"] == {
        "direct_enrichment": 1,
        "direct_upgrade": 1,
        "direct_fill": 1,
        "validated_low_fill": 1,
    }
    assert summary["validated_low_audit"]["old_low_invalidated_in_strict_rebuild"] == 1
    result = pd.read_csv(output / "wave35_species_axis_coverage.csv.gz").fillna("")
    upgraded = result.loc[
        (result.accepted_species == "Genus species1")
        & (result.axis == "reproductive_assurance")
    ].iloc[0]
    inferred = result.loc[
        (result.accepted_species == "Genus species3")
        & (result.axis == "reproductive_assurance")
    ].iloc[0]
    assert upgraded.quality == "high"
    assert upgraded.trait_composition == 'self_incompatibility=["SC"]'
    assert inferred.quality == "low"
    assert inferred.trait_names == "self_incompatibility"


def test_rejects_genus_axis_join_that_drops_trait_identity(tmp_path: Path) -> None:
    baseline, inputs = _write_inputs(tmp_path, wrong_low_axis=True)
    with pytest.raises(ValueError, match="cross-trait or cross-axis"):
        build_wave35_overlay(
            baseline_csv=baseline,
            input_dir=inputs,
            output_dir=tmp_path / "output",
        )
