import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave20_checkpoint_20260821"
)
SOURCE_GROUP = "targeted_support2_wave20_checkpoint_20260821"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_wave20_has_unique_reviewed_direct_species_traits() -> None:
    evidence = _read("targeted_support2_wave20_evidence_20260821.csv")
    audit = _read("targeted_support2_wave20_manual_audit_20260821.csv")

    assert len(evidence) == len(audit) == 24
    assert evidence["accepted_species"].nunique() == 16
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert not evidence["candidate_id"].duplicated().any()
    assert set(evidence["evidence_scope"]) == {"species_direct", "synonym_direct"}
    assert evidence["evidence_scope"].eq("synonym_direct").sum() == 3
    assert set(evidence["evidence_quality"]) == {"high", "medium"}
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_lineage"].ne("").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["provenance_complete"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_wave20_preserves_trait_specific_values_and_tier_limits() -> None:
    evidence = _read("targeted_support2_wave20_evidence_20260821.csv").set_index(
        ["accepted_species", "trait_name"]
    )

    assert evidence.loc[("Helicteres guazumifolia", "self_incompatibility")][
        "normalized_value"
    ] == "SC"
    assert evidence.loc[("Patersonia sericea", "tube_depth_class")][
        "normalized_value"
    ] == "deep"
    assert evidence.loc[("Amischotolype tenuis", "flower_primary_color")][
        "normalized_value"
    ] == "blue_purple|red_pink"
    assert evidence.loc[("Catesbaea parviflora", "floral_form")][
        "normalized_value"
    ] == "tubular"
    assert evidence.loc[("Critonia macropoda", "inflorescence_display")][
        "normalized_value"
    ] == "composite_display|umbel_corymb"
    assert evidence.loc[("Aeonium appendiculatum", "flower_size_class")][
        "normalized_value"
    ] == "medium"
    assert evidence.loc[("Arytera microphylla", "flower_primary_color")][
        "normalized_value"
    ] == "white"
    assert evidence.loc[("Paraboea barnettiae", "inflorescence_display")][
        "normalized_value"
    ] == "umbel_corymb"
    assert evidence.loc[("Xanthosoma brasiliense", "flower_primary_color")][
        "normalized_value"
    ] == "white|yellow_orange"
    nursery = evidence.loc[("Schizomeria ovata", "flower_primary_color")]
    assert nursery["source_tier"] == "C"
    assert nursery["evidence_quality"] == "medium"


def test_wave20_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260821.csv")
    audit = _read("combined_curated_manual_audit_20260821.csv")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_935
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 24


def test_wave20_manifest_records_ceiling_cost_and_fail_closed_rules() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave20.json").read_text(
            encoding="utf-8"
        )
    )

    assert manifest["baseline_formal_run_id"] == 32347770192
    assert manifest["accepted_evidence_rows"] == 24
    assert manifest["accepted_species_trait"] == 24
    assert manifest["accepted_species"] == 16
    assert len(manifest["targeted_support2_rules"]) == 14
    assert manifest["theoretical_rule_cells_touched"] == 179
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["guardrails"]["genus_axis_only_join"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
    assert (
        manifest["guardrails"]["cultivar_or_hybrid_transferred_to_wild_species"]
        is False
    )
