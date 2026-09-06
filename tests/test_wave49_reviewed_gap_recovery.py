from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from island_v2.all_evidence_trait_audit import EVIDENCE_COLUMNS
from island_v2.wave49_reviewed_gap_recovery import (
    _load_manifest,
    _select_gap_evidence,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = (
    REPO_ROOT
    / "data/v2/staging/traits/wave49_reviewed_gap_recovery/source_manifest.json"
)


def test_wave49_manifest_pins_reviewed_source_receipts() -> None:
    manifest = _load_manifest(REPO_ROOT, MANIFEST)
    assert len(manifest["sources"]) == 11
    by_id = {source["source_id"]: source for source in manifest["sources"]}
    assert by_id["wfo_global_six"]["quality_overrides"] == {
        "flower_primary_color": "medium"
    }
    assert by_id["wfo_global_six"]["excluded_traits"] == {
        "tube_depth_class": "manual audit precision 0.65"
    }
    assert by_id["wfo_combined_high_yield"]["excluded_traits"] == {
        "flower_size_class": "manual audit precision 0.80"
    }


def _evidence_row(species: str, trait: str, quality: str = "high") -> dict[str, str]:
    row = {column: "" for column in EVIDENCE_COLUMNS}
    row.update(
        {
            "accepted_species": species,
            "axis": "flower_colour",
            "trait_name": trait,
            "normalized_value": "white",
            "quality": quality,
            "source_group": "official_flora",
            "source_provider": "fixture flora",
            "source_url": f"https://example.org/{species.replace(' ', '-')}",
            "source_record_id": species,
            "source_citation": "Fixture flora species treatment",
            "source_excerpt": "Flowers white.",
            "evidence_scope": "species_direct",
            "name_match_method": "accepted_name_exact",
            "source_lineage": f"fixture:{species}",
            "lineage_method": "original_treatment",
            "source_run_id": "fixture-run",
            "source_artifact": "fixture-artifact",
            "source_file": "fixture.csv",
            "acceptance_contract": "fixture_exact_quote_v1",
        }
    )
    return row


def test_select_gap_evidence_skips_completed_keys_and_downgrades(
    tmp_path: Path,
) -> None:
    evidence_path = tmp_path / "evidence.csv"
    pd.DataFrame(
        [
            _evidence_row("Alpha one", "flower_primary_color"),
            _evidence_row("Beta two", "flower_primary_color"),
        ]
    ).to_csv(evidence_path, index=False)
    manifest = {
        "sources": [
            {
                "source_id": "fixture",
                "evidence_path": "evidence.csv",
                "approved_traits": ["flower_primary_color"],
                "quality_overrides": {"flower_primary_color": "medium"},
                "review_precision": 0.91,
                "cultivar_contamination_rate": 0.0,
            }
        ]
    }
    selected, audit = _select_gap_evidence(
        repo_root=tmp_path,
        manifest=manifest,
        target_species={"Alpha one", "Beta two"},
        previous_direct_keys={("Alpha one", "flower_primary_color")},
    )
    assert selected[["accepted_species", "trait_name"]].to_records(
        index=False
    ).tolist() == [("Beta two", "flower_primary_color")]
    assert selected.iloc[0]["quality"] == "medium"
    assert audit.iloc[0]["eligible_target_rows"] == 2
    assert audit.iloc[0]["missing_species_trait_rows"] == 1


def test_select_gap_evidence_rejects_fuzzy_identity(tmp_path: Path) -> None:
    row = _evidence_row("Alpha one", "flower_primary_color")
    row["name_match_method"] = "fuzzy_match"
    pd.DataFrame([row]).to_csv(tmp_path / "evidence.csv", index=False)
    manifest = {
        "sources": [
            {
                "source_id": "fixture",
                "evidence_path": "evidence.csv",
                "approved_traits": None,
                "quality_overrides": {},
                "review_precision": 1.0,
                "cultivar_contamination_rate": 0.0,
            }
        ]
    }
    with pytest.raises(ValueError, match="fuzzy name match"):
        _select_gap_evidence(
            repo_root=tmp_path,
            manifest=manifest,
            target_species={"Alpha one"},
            previous_direct_keys=set(),
        )
