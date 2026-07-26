from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import rebuild_validated_low
from island_v2.open_web_multidomain import _hash_sample


def test_multidomain_audit_sample_is_hash_based_trait_stratified_and_unique() -> None:
    rows = []
    for trait in ("flower_primary_color", "floral_form"):
        for index in range(20):
            rows.append(
                {
                    "candidate_id": f"{trait}-{index}",
                    "trait_name": trait,
                    "source_url": f"https://example.test/{trait}/{index}",
                }
            )
    frame = pd.DataFrame(rows)
    first = _hash_sample(frame, sample_pages=20, seed="fixed")
    second = _hash_sample(frame, sample_pages=20, seed="fixed")
    assert first["candidate_id"].tolist() == second["candidate_id"].tolist()
    assert first["source_url"].nunique() == 20
    assert first.groupby("trait_name").size().to_dict() == {
        "floral_form": 10,
        "flower_primary_color": 10,
    }


def test_low_rebuild_applies_rule_by_trait_not_axis() -> None:
    coverage = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "axis": "floral_structural_complexity",
                "after_quality": quality,
            }
            for species, quality in (
                ("Alpha one", "high"),
                ("Alpha two", "medium"),
                ("Alpha three", "high"),
                ("Alpha target", ""),
            )
        ]
    )
    evidence = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "trait_name": "floral_symmetry",
                "normalized_value": "actinomorphic",
                "quality": quality,
                "evidence_scope": "species_direct",
                "source_lineage": f"url:{species}",
            }
            for species, quality in (
                ("Alpha one", "high"),
                ("Alpha two", "medium"),
                ("Alpha three", "high"),
            )
        ]
    )
    rules, low, conflicts, report = rebuild_validated_low(
        coverage=coverage,
        evidence=evidence,
        promoted=pd.DataFrame(),
    )
    assert conflicts.empty
    assert bool(rules.iloc[0]["eligible"])
    target = low.loc[low["accepted_species"].eq("Alpha target")]
    assert set(target["trait_name"]) == {"floral_symmetry"}
    assert report["rebuilt_validated_low_species_axis"] == 1


def test_pollen_vector_does_not_fill_structure_axis() -> None:
    coverage = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "axis": "floral_structural_complexity",
                "after_quality": quality,
            }
            for species, quality in (
                ("Beta one", "high"),
                ("Beta two", "high"),
                ("Beta three", "high"),
                ("Beta target", ""),
            )
        ]
    )
    evidence = pd.DataFrame(
        [
            {
                "accepted_species": species,
                "trait_name": "pollen_vector_mode",
                "normalized_value": "biotic",
                "quality": "high",
                "evidence_scope": "species_direct",
                "source_lineage": f"url:{species}",
            }
            for species in ("Beta one", "Beta two", "Beta three")
        ]
    )
    rules, low, conflicts, report = rebuild_validated_low(
        coverage=coverage,
        evidence=evidence,
        promoted=pd.DataFrame(),
    )
    assert rules.empty
    assert low.empty
    assert conflicts.empty
    assert report["rebuilt_validated_low_species_axis"] == 0


def test_precision_denominator_keeps_bad_accepts() -> None:
    source = Path("src/island_v2/open_web_finalize.py").read_text(encoding="utf-8")
    assert "len(accepted_audit) / len(reviewed)" in source
    assert "len(accepted_audit) + len(rejected)" not in source
