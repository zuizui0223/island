import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.ferrer_2024_powo_synonym_checkpoint import (
    MAPPING_AUDIT_SHA256,
    RESPONSES_SHA256,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

SOURCE_DIR = Path("data/v2/external/traits/ferrer_2024_two_backbone")
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "ferrer_2024_powo_synonym_checkpoint_20260821"
)
SOURCE_GROUP = "ferrer_2024_powo_synonym_checkpoint_20260821"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_ferrer_two_backbone_snapshots_are_immutable() -> None:
    assert (
        _file_sha256(SOURCE_DIR / "ferrer_two_backbone_responses_20260821.jsonl.gz")
        == RESPONSES_SHA256
    )
    assert (
        _file_sha256(SOURCE_DIR / "ferrer_two_backbone_mapping_audit_20260821.csv.gz")
        == MAPPING_AUDIT_SHA256
    )


def test_ferrer_wfo_powo_mappings_are_exact_and_family_consistent() -> None:
    mappings = pd.read_csv(
        SOURCE_DIR / "powo_wfo_selected_mappings_20260821.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")
    joined = mappings.merge(master[["accepted_species", "family"]], on="accepted_species")
    assert len(joined) == 3
    assert joined["source_family"].eq(joined["family"]).all()
    assert joined["wfo_classification_version"].eq("2026-06").all()
    assert joined["wfo_id"].str.startswith("wfo-").all()
    assert joined.apply(
        lambda row: "Synonym of: " + row.accepted_species in row.powo_visible_result,
        axis=1,
    ).all()


def test_ferrer_checkpoint_accepts_only_three_medium_synonym_direct_sc_rows() -> None:
    evidence = _read("ferrer_2024_powo_synonym_evidence_20260821.csv")
    audit = _read("ferrer_2024_powo_synonym_manual_audit_20260821.csv")
    assert len(evidence) == len(audit) == 3
    assert set(evidence["accepted_species"]) == {
        "Agalinis decemloba",
        "Cuphea flava",
        "Chaetogastra herbacea",
    }
    assert evidence["trait_name"].eq("self_incompatibility").all()
    assert evidence["normalized_value"].eq("SC").all()
    assert evidence["evidence_quality"].eq("medium").all()
    assert evidence["evidence_scope"].eq("synonym_direct").all()
    assert evidence["name_match_method"].eq("synonym_exact_two_backbone").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_lineage"].nunique() == 3
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()


def test_ferrer_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260821.csv")
    audit = _read("combined_curated_manual_audit_20260821.csv")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")
    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_947
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 3


def test_ferrer_manifest_does_not_claim_theoretical_rule_gain() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_ferrer_2024_powo.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["unique_unmatched_names_queried"] == 1_171
    assert manifest["wfo_requests"] == 1_171
    assert manifest["powo_manual_result_pages_reviewed"] == 3
    assert manifest["accepted_evidence_rows"] == 3
    assert manifest["theoretical_rule_candidate"]["potential_unresolved_cells"] == 12
    assert manifest["theoretical_rule_candidate"]["formal_gain_claimed"] is False
    assert manifest["guardrails"]["single_backbone_match_accepted"] is False
    assert manifest["guardrails"]["genus_axis_only_join"] is False
