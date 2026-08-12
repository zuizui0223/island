from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import (
    validate_individually_reviewed_evidence,
)
from island_v2.rule_unlock_wave2_checkpoint import reviewed_rows

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812"
)


def test_wave2_rows_are_trait_specific_and_source_backed() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 14
    assert evidence["accepted_species"].nunique() == 8
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 14
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["source_excerpt"].str.len().gt(30).all()
    assert evidence["source_lineage"].nunique() == 9
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 12,
        "medium": 2,
    }

    reproductive = evidence.loc[evidence["axis"].eq("reproductive_assurance")]
    assert len(reproductive) == 12
    assert set(reproductive["trait_name"]) == {
        "autonomous_selfing_capacity",
        "self_incompatibility",
    }
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_zero_bagged_set_does_not_overwrite_hand_self_compatibility() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    for species in ["Hosta ventricosa", "Melastoma malabathricum"]:
        rows = evidence.loc[evidence["accepted_species"].eq(species)]
        observed = set(zip(rows["trait_name"], rows["normalized_value"]))
        assert ("autonomous_selfing_capacity", "absent") in observed
        assert ("self_incompatibility", "SC") in observed


def test_morphology_rows_are_species_direct_medium_only() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    morphology = evidence.loc[
        evidence["trait_name"].eq("floral_symmetry")
    ].set_index("accepted_species")

    assert set(morphology.index) == {
        "Drypetes assamica",
        "Melastoma malabathricum",
    }
    assert morphology["normalized_value"].eq("actinomorphic").all()
    assert morphology["evidence_quality"].eq("medium").all()


def test_committed_combined_checkpoint_passes_review_gate() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "combined_curated_evidence_20260812.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260812.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(audit) == len(accepted) == 799
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    wave = accepted.loc[accepted["source_group"].eq(
        "rule_unlock_wave2_checkpoint_20260812"
    )]
    assert len(wave) == 14
