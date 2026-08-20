import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave16_checkpoint_20260820"
)
SOURCE_GROUP = "targeted_support2_wave16_checkpoint_20260820"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_wave16_has_27_unique_reviewed_direct_species_traits() -> None:
    evidence = _read("targeted_support2_wave16_evidence_20260820.csv")
    audit = _read("targeted_support2_wave16_manual_audit_20260820.csv")

    assert len(evidence) == len(audit) == 27
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert not evidence["candidate_id"].duplicated().any()
    assert set(evidence["evidence_scope"]) == {"species_direct"}
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


def test_wave16_preserves_trait_contract_and_primary_evidence() -> None:
    evidence = _read("targeted_support2_wave16_evidence_20260820.csv").set_index(
        ["accepted_species", "trait_name"]
    )

    assert (
        not evidence.index.get_level_values("trait_name")
        .isin({"pollen_vector_mode", "reward_type"})
        .any()
    )
    rhamnus = evidence.loc[("Rhamnus alaternus", "mating_system")]
    assert rhamnus["normalized_value"] == "predominantly_outcrossing"
    assert rhamnus["source_lineage"] == ("institutional-species-page:up-mhnc:Rhamnus_alaternus")
    assert rhamnus["content_sha256"] == (
        "4763748cc8f47740c773c4a8e91e2f3a262fb5b7166a7ffd59215f396b2fbb19"
    )
    viscum = evidence.loc[("Viscum coloratum", "mating_system")]
    assert viscum["normalized_value"] == "predominantly_outcrossing"
    assert viscum["source_lineage"] == "doi:10.3390/biom15070974"
    assert viscum["content_sha256"] == (
        "95543993c10ec2a67dcc9fa0f54e27019064e642166dd520310be9ea6f58b37b"
    )


def test_wave16_snapshot_uses_three_completed_artifacts() -> None:
    source = _read("source_artifact_rows_wave16.csv")

    assert len(source) == 25
    assert set(source["source_artifact_id"]) == {
        "9014569071",
        "9018965291",
        "9020100238",
    }
    assert set(source["source_run_id"]) == {
        "31233181709",
        "31247162370",
        "31249406427",
    }
    assert not source["supporting_excerpt"].str.contains(r"(?i)solitary spikelet").any()
    assert not (
        source["domain"].eq("plants.ces.ncsu.edu")
        & source["trait_name"].eq("flower_size_class")
        & source["raw_value"].eq("less than 1 inch")
    ).any()


def test_wave16_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260820.csv")
    audit = _read("combined_curated_manual_audit_20260820.csv")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_885
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 27


def test_wave16_manifest_records_artifacts_cost_and_fail_closed_rules() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave16.json").read_text(encoding="utf-8")
    )

    assert manifest["accepted_species_trait_rows"] == 27
    assert manifest["targeted_support2_rule_count"] == 4
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert sum(item["accepted_rows"] for item in manifest["source_artifacts"]) == 25
    assert len(manifest["source_artifacts"]) == 3
    assert manifest["guardrails"]["search_snippet_as_evidence"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
    assert manifest["guardrails"]["genus_axis_only_join"] is False
    assert manifest["guardrails"]["min_species_two_production"] is False
    assert manifest["guardrails"]["solitary_spikelet_mapped_to_solitary_flower_display"] is False
