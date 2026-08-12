from pathlib import Path

import pandas as pd

from island_v2.morphology_rule_unlock_checkpoint import (
    SOURCE_GROUP,
    reviewed_rows,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "morphology_rule_unlock_checkpoint_20260813"
)


def test_rows_are_direct_trait_specific_and_do_not_duplicate_pngtrees() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 9
    assert evidence[["accepted_species", "trait_name"]].duplicated().sum() == 0
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 6,
        "medium": 3,
    }
    assert evidence["trait_name"].value_counts().to_dict() == {
        "inflorescence_display": 4,
        "floral_form": 2,
        "flower_size_class": 2,
        "floral_symmetry": 1,
    }
    assert "Gonystylus macrophyllus" not in set(evidence["accepted_species"])
    assert "Sloanea dasycarpa" not in set(evidence["accepted_species"])
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_committed_combined_checkpoint_passes_individual_review_gate() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "combined_curated_evidence_20260813.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260813.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(audit) == len(accepted) == 1_746
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    morphology = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(morphology) == 9
    assert morphology["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert morphology["source_lineage"].nunique() == 9
