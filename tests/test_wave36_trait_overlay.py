from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

import island_v2.wave36_trait_overlay as overlay
from island_v2.wave35_trait_overlay import AXES


@pytest.fixture(autouse=True)
def small_universe(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(overlay, "EXPECTED_SPECIES", 3)
    monkeypatch.setattr(overlay, "EXPECTED_CELLS", 9)


def test_lossless_direct_and_trait_specific_low_overlay(tmp_path: Path) -> None:
    species = ["Alpha one", "Alpha two", "Alpha three"]
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
    baseline_path = tmp_path / "baseline.csv.gz"
    baseline.to_csv(baseline_path, index=False)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    pd.DataFrame(
        [
            {
                "accepted_species": species[0],
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "resolution_status": "resolved",
                "quality": "medium",
                "state_set": '["SI"]',
                "source_groups": "useful_plants_citation_remine",
                "source_lineages": "citation:one",
            }
        ]
    ).to_csv(inputs / "wave36_resolved_direct_species_trait.csv.gz", index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": species[1],
                "genus": "Alpha",
                "axis": "floral_structural_complexity",
                "trait_name": "floral_form",
                "inferred_state_set": '["bell_campanulate"]',
                "quality": "low",
                "family_inference_used": "False",
                "global_fallback_used": "False",
                "source_lineage": "validated-low:Alpha:floral_form",
            }
        ]
    ).to_csv(inputs / "wave36_candidate_validated_low_species_trait.csv.gz", index=False)
    pd.DataFrame(
        [
            {
                "genus": "Alpha",
                "axis": "floral_structural_complexity",
                "trait_name": "floral_form",
                "inferred_state_set": '["bell_campanulate"]',
                "eligible": "True",
                "n_direct_species": "3",
                "species_loo_accuracy": "1.0",
                "lineage_loo_accuracy": "1.0",
            }
        ]
    ).to_csv(inputs / "wave36_provider_touched_new_rule_audit.csv.gz", index=False)
    hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs.iterdir()}
    (inputs / "wave36_source_manifest.json").write_text(
        json.dumps(
            {
                "contract": "wave36_nhm_reproductive_remine_source_manifest_v1",
                "file_sha256": hashes,
            }
        ),
        encoding="utf-8",
    )
    summary = overlay.build_wave36_overlay(
        baseline_csv=baseline_path,
        input_dir=inputs,
        output_dir=tmp_path / "output",
    )
    assert summary["delta"]["gross_gain_species_axis"] == 2
    assert summary["delta"]["loss_species_axis"] == 0
    assert summary["delta"]["action_counts"] == {
        "direct_fill": 1,
        "validated_low_fill": 1,
    }
    assert summary["wave36_after"]["quality_counts"] == {
        "high": 0,
        "medium": 1,
        "low": 1,
    }
