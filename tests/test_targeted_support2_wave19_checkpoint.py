import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/targeted_support2_wave19_checkpoint_20260820"
)
SOURCE_GROUP = "targeted_support2_wave19_checkpoint_20260820"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_wave19_has_five_unique_reviewed_direct_species_traits() -> None:
    evidence = _read("targeted_support2_wave19_evidence_20260820.csv")
    audit = _read("targeted_support2_wave19_manual_audit_20260820.csv")

    assert len(evidence) == len(audit) == 5
    assert evidence["accepted_species"].nunique() == 2
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert not evidence["candidate_id"].duplicated().any()
    assert set(evidence["evidence_scope"]) == {"species_direct"}
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


def test_wave19_preserves_trait_specific_and_multistate_values() -> None:
    evidence = _read("targeted_support2_wave19_evidence_20260820.csv").set_index(
        ["accepted_species", "trait_name"]
    )

    assert evidence.loc[("Lilium philadelphicum", "mating_system")][
        "normalized_value"
    ] == "predominantly_outcrossing"
    assert evidence.loc[("Pulicaria arabica", "floral_symmetry")][
        "normalized_value"
    ] == "actinomorphic|zygomorphic"
    assert evidence.loc[("Pulicaria arabica", "inflorescence_display")][
        "normalized_value"
    ] == "composite_display|umbel_corymb"
    assert evidence.loc[("Pulicaria arabica", "flower_primary_color")][
        "normalized_value"
    ] == "yellow_orange"
    assert evidence.loc[("Pulicaria arabica", "flower_size_class")][
        "normalized_value"
    ] == "small"


def test_wave19_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260820.csv")
    audit = _read("combined_curated_manual_audit_20260820.csv")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_911
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 5


def test_wave19_manifest_records_ceiling_cost_and_fail_closed_rules() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave19.json").read_text(
            encoding="utf-8"
        )
    )

    assert manifest["baseline_formal_run_id"] == 32344088168
    assert manifest["accepted_evidence_rows"] == 5
    assert manifest["accepted_species_trait"] == 5
    assert manifest["accepted_species"] == 2
    assert len(manifest["targeted_support2_rules"]) == 3
    assert manifest["theoretical_rule_cells_touched"] == 67
    assert manifest["formal_search_api_queries"] == 0
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["guardrails"]["genus_axis_only_join"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
    assert manifest["guardrails"]["paywall_or_access_control_bypassed"] is False
    assert (
        manifest["guardrails"][
            "cited_characterization_counted_as_independent_replication"
        ]
        is False
    )
