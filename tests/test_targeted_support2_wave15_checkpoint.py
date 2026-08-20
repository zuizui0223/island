import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave15_checkpoint_20260820"
)
SOURCE_GROUP = "targeted_support2_wave15_checkpoint_20260820"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_wave15_has_100_unique_reviewed_direct_species_traits() -> None:
    evidence = _read("targeted_support2_wave15_evidence_20260820.csv")
    audit = _read("targeted_support2_wave15_manual_audit_20260820.csv")

    assert len(evidence) == len(audit) == 100
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert not evidence["candidate_id"].duplicated().any()
    assert set(evidence["evidence_scope"]) == {"species_direct", "synonym_direct"}
    assert set(evidence["evidence_quality"]) <= {"high", "medium"}
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_url"].str.startswith("http").all()
    assert evidence["source_lineage"].ne("").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["provenance_complete"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_wave15_preserves_trait_specific_contract() -> None:
    evidence = _read("targeted_support2_wave15_evidence_20260820.csv")
    assert not evidence["trait_name"].isin({"pollen_vector_mode", "reward_type"}).any()
    reproduction = evidence.loc[evidence["axis"].eq("reproductive_assurance")]
    assert set(reproduction["trait_name"]) <= {
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "mating_system",
        "cleistogamy",
    }
    assert not reproduction["trait_name"].eq("").any()
    stemonoporus = evidence.loc[evidence["accepted_species"].eq("Stemonoporus ceylanicus")].iloc[0]
    assert stemonoporus["evidence_scope"] == "synonym_direct"
    assert stemonoporus["name_match_method"] == "synonym_exact"
    assert stemonoporus["matched_page_name"] == "Stemonoporus Wightii"


def test_wave15_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260820.csv")
    audit = _read("combined_curated_manual_audit_20260820.csv")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_858
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 100


def test_wave15_manifest_records_artifacts_cost_and_fail_closed_rules() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave15.json").read_text(encoding="utf-8")
    )
    assert manifest["accepted_species_trait_rows"] == 100
    assert manifest["targeted_support2_rule_count"] == 20
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert sum(item["accepted_rows"] for item in manifest["source_artifacts"]) == 83
    assert len(manifest["source_artifacts"]) == 5
    assert manifest["guardrails"]["search_snippet_as_evidence"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
    assert manifest["guardrails"]["genus_axis_only_join"] is False
    assert manifest["guardrails"]["min_species_two_production"] is False
