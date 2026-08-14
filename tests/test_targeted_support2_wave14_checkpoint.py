import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_support2_wave14_checkpoint_20260814"
)
SOURCE_GROUP = "targeted_support2_wave14_checkpoint_20260814"


def test_wave14_rows_are_reviewed_species_direct_individual_traits() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "targeted_support2_wave14_evidence_20260814.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "targeted_support2_wave14_manual_audit_20260814.csv", dtype=str
    ).fillna("")

    assert len(evidence) == len(audit) == 8
    assert evidence["accepted_species"].nunique() == 5
    assert set(evidence["evidence_scope"]) == {"species_direct", "synonym_direct"}
    assert evidence["evidence_quality"].eq("high").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["inference_rule"].eq("").all()
    assert set(evidence["trait_name"]) == {
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "floral_symmetry",
    }
    assert not evidence["trait_name"].isin({"pollen_vector_mode", "reward_type"}).any()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_lineage"].nunique() == 3
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert audit["decision"].eq("accept").all()


def test_wave14_preserves_two_backbone_strict_synonym_identity() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "targeted_support2_wave14_evidence_20260814.csv", dtype=str
    ).fillna("")
    villosa = evidence.loc[evidence["accepted_species"].eq("Goeppertia villosa")]

    assert len(villosa) == 2
    assert villosa["matched_page_name"].eq("Calathea villosa").all()
    assert villosa["evidence_scope"].eq("synonym_direct").all()
    assert villosa["name_match_method"].eq("synonym_exact").all()
    assert villosa["name_resolution_lineage"].str.contains("gbif:").all()
    assert villosa["name_resolution_lineage"].str.contains("wfo:").all()


def test_wave14_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "combined_curated_evidence_20260814.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260814.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_759
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 8


def test_wave14_manifest_is_fail_closed_and_costed() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave14.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["manual_open_web_queries"] == 43
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["accepted_species_trait_rows"] == 8
    assert len(manifest["withheld_or_rejected"]) == 6
    assert manifest["guardrails"]["search_snippet_as_evidence"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
