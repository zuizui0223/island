import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.targeted_support2_wave24_checkpoint import build_checkpoint

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_support2_wave24_checkpoint_20260824"
)
SOURCE_GROUP = "targeted_support2_wave24_checkpoint_20260824"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_wave24_has_unique_reviewed_trait_specific_rows() -> None:
    evidence = _read("targeted_support2_wave24_evidence_20260824.csv")
    audit = _read("targeted_support2_wave24_manual_audit_20260824.csv")

    assert len(evidence) == len(audit) == 22
    assert evidence["accepted_species"].nunique() == 21
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert not evidence["candidate_id"].duplicated().any()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["name_match_method"].eq("accepted_name_exact").all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_lineage"].nunique() == 21
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 15,
        "medium": 7,
    }
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["provenance_complete"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_wave24_keeps_rejected_candidates_visible_and_unpromoted() -> None:
    evidence = _read("targeted_support2_wave24_evidence_20260824.csv")
    rejected = _read("targeted_support2_wave24_rejected_candidates_20260824.csv")

    assert len(rejected) == 13
    assert rejected["decision"].eq("reject").all()
    assert rejected["decision_reason"].ne("").all()
    rejected_pairs = set(zip(rejected["accepted_species"], rejected["queried_trait"]))
    accepted_pairs = set(zip(evidence["accepted_species"], evidence["trait_name"]))
    assert not rejected_pairs & accepted_pairs
    assert "wrong floral organ" in set(rejected["decision_reason"])
    assert "strict species identity unresolved" in set(rejected["decision_reason"])


def test_wave24_combined_checkpoint_passes_common_review_gate(tmp_path: Path) -> None:
    build_checkpoint(tmp_path)
    evidence = pd.read_csv(
        tmp_path / "combined_curated_evidence_20260824.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        tmp_path / "combined_curated_manual_audit_20260824.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_980
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 22


def test_wave24_manifest_records_rule_ceiling_cost_and_guardrails() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave24.json").read_text(
            encoding="utf-8"
        )
    )

    assert manifest["fixed_integrated_baseline_run_id"] == 31488685866
    assert manifest["baseline_formal_run_id"] == 32702160934
    assert manifest["accepted_evidence_rows"] == 22
    assert manifest["accepted_species_trait"] == 22
    assert manifest["accepted_species"] == 21
    assert manifest["accepted_source_lineages"] == 21
    assert manifest["quality_counts"] == {"high": 15, "medium": 7}
    assert manifest["rejected_reviewed_candidates"] == 13
    assert len(manifest["targeted_support2_rules"]) == 22
    assert manifest["theoretical_rule_cells_touched"] == 203
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["guardrails"]["genus_axis_only_join"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
    assert manifest["guardrails"]["min_species_two_production"] is False


def test_wave24_workflow_verifies_fixed_goals_and_new_rules() -> None:
    workflow = Path(".github/workflows/open-web-review-promote.yml").read_text(
        encoding="utf-8"
    )

    assert "Verify the Wave 24 support-two checkpoint goals" in workflow
    assert 'assert len(targeted_new) >= 20' in workflow
    assert 'assert shared["coverage_increment_species_axis"] == 194' in workflow
    assert 'assert coverage["filled_species_axis"] >= 159443' in workflow
    assert 'assert reproduction["filled_species"] == 31637' in workflow
    assert 'assert coverage["species_by_filled_axis_count"]["0"] <= 14500' in workflow
    assert "wave24_support2_goal_verification.json" in workflow
    assert "targeted_support2_wave24_evidence_20260824.csv" in workflow
    assert "targeted_support2_wave24_manual_audit_20260824.csv" in workflow
    assert "combined_curated_evidence_20260824.csv" not in workflow
