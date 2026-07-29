from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.open_web_pilot import (
    _language_from_evidence,
    _production_approved_traits,
    _source_type,
    _stratified_unique_page_audit,
)


def test_baseline_metadata_audit_exposes_missing_language_instead_of_guessing() -> None:
    row = pd.Series(
        {
            "source_url": "https://unknown.example/species",
            "source_group": "three_pass_medium",
            "source_provider": "specialist_site",
        }
    )
    assert _language_from_evidence(row) == "unknown_not_recorded"
    assert _source_type(row) == "mixed_public_web"


def test_workflows_and_entry_points_exist() -> None:
    root = Path(__file__).resolve().parents[1]
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    assert "island-v2-open-web-search" in pyproject
    assert "island-v2-open-web-evidence" in pyproject
    assert "island-v2-open-web-pilot" in pyproject


def test_manual_audit_sample_is_trait_stratified_and_page_unique() -> None:
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": f"{trait}-{index}",
                "trait_name": trait,
                "source_tier": "B",
                "source_url": f"https://example.test/{trait}/{index}",
            }
            for trait in ("floral_form", "flower_primary_color")
            for index in range(15)
        ]
    )
    audit = _stratified_unique_page_audit(candidates, 20)
    assert len(audit) == 20
    assert audit["source_url"].nunique() == 20
    assert audit.groupby("trait_name").size().to_dict() == {
        "floral_form": 10,
        "flower_primary_color": 10,
    }


def test_production_approves_only_sufficiently_audited_precise_traits() -> None:
    assert _production_approved_traits(
        {
            "floral_form": {"audited": 10, "precision": 0.9},
            "floral_symmetry": {"audited": 9, "precision": 1.0},
            "flower_size_class": {"audited": 20, "precision": 0.89},
        }
    ) == ["floral_form"]
