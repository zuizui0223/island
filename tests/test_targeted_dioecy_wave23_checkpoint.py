import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_dioecy_wave23_checkpoint_20260821"
)
SOURCE_GROUP = "targeted_dioecy_wave23_checkpoint_20260821"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_wave23_has_unique_reviewed_exact_species_dioecy_rows() -> None:
    evidence = _read("targeted_dioecy_wave23_evidence_20260821.csv")
    audit = _read("targeted_dioecy_wave23_manual_audit_20260821.csv")

    assert len(evidence) == len(audit) == 11
    assert evidence["accepted_species"].nunique() == 11
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert not evidence["candidate_id"].duplicated().any()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["name_match_method"].eq("accepted_name_exact").all()
    assert evidence["trait_name"].eq("mating_system").all()
    assert evidence["axis"].eq("reproductive_assurance").all()
    assert evidence["normalized_value"].eq("predominantly_outcrossing").all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_lineage"].ne("").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 9,
        "medium": 2,
    }
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["provenance_complete"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_wave23_preserves_shared_original_lineages() -> None:
    evidence = _read("targeted_dioecy_wave23_evidence_20260821.csv").set_index(
        "accepted_species"
    )

    asparagus = {
        evidence.loc["Asparagus cochinchinensis", "source_lineage"],
        evidence.loc["Asparagus schoberioides", "source_lineage"],
    }
    assert asparagus == {"doi:10.1371/journal.pone.0266376"}
    assert (
        evidence.loc["Myrsine australis", "source_lineage"]
        == "provider_treatment:nzpcn:myrsine-australis"
    )
    assert evidence["source_lineage"].nunique() == 10
    assert not evidence["trait_name"].isin(
        ["self_incompatibility", "autonomous_selfing_capacity"]
    ).any()


def test_wave23_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260821.csv")
    audit = _read("combined_curated_manual_audit_20260821.csv")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_958
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 11


def test_wave23_manifest_records_ceiling_cost_and_fail_closed_rules() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave23.json").read_text(
            encoding="utf-8"
        )
    )

    assert manifest["baseline_formal_run_id"] == 32449801624
    assert manifest["accepted_evidence_rows"] == 11
    assert manifest["accepted_species_trait"] == 11
    assert manifest["accepted_species"] == 11
    assert manifest["accepted_source_lineages"] == 10
    assert manifest["quality_counts"] == {"high": 9, "medium": 2}
    assert len(manifest["targeted_support2_rules"]) == 5
    assert manifest["theoretical_rule_cells_touched"] == 831
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["guardrails"]["genus_axis_only_join"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
    assert manifest["guardrails"]["dioecy_mapped_to_self_incompatibility"] is False
    assert manifest["guardrails"]["dioecy_mapped_to_autonomous_selfing"] is False
