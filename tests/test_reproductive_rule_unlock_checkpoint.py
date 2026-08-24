from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import (
    validate_individually_reviewed_evidence,
)
from island_v2.reproductive_rule_unlock_checkpoint import reviewed_rows

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "reproductive_rule_unlock_checkpoint_20260812"
)


def test_rows_preserve_reproductive_trait_boundaries() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 9
    assert evidence["accepted_species"].nunique() == 6
    assert evidence["trait_name"].value_counts().to_dict() == {
        "mating_system": 4,
        "autonomous_selfing_capacity": 3,
        "self_incompatibility": 2,
    }
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )
    assert evidence["source_lineage"].nunique() == 6
    assert evidence["evidence_quality"].eq("high").all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["source_excerpt"].str.len().gt(30).all()

    sonchus = evidence.loc[evidence["accepted_species"].eq("Sonchus microcephalus")]
    assert set(zip(sonchus["trait_name"], sonchus["normalized_value"])) == {
        ("autonomous_selfing_capacity", "autonomous"),
        ("self_incompatibility", "SC"),
    }
    fritillaria = evidence.loc[
        evidence["accepted_species"].eq("Fritillaria koidzumiana")
    ]
    assert set(zip(fritillaria["trait_name"], fritillaria["normalized_value"])) == {
        ("autonomous_selfing_capacity", "absent"),
        ("self_incompatibility", "mixed_or_variable"),
        ("mating_system", "predominantly_outcrossing"),
    }


def test_dioecy_does_not_cross_map_to_other_reproductive_traits() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    vitis = evidence.loc[evidence["accepted_species"].str.startswith("Vitis ")]

    assert len(vitis) == 2
    assert vitis["trait_name"].eq("mating_system").all()
    assert vitis["normalized_value"].eq("predominantly_outcrossing").all()
    assert not vitis["source_excerpt"].str.contains("autonomous selfing", case=False).any()


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

    assert len(evidence) == len(audit) == len(accepted) == 785
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    unlock = accepted.loc[accepted["source_group"].eq(
        "reproductive_rule_unlock_checkpoint_20260812"
    )]
    assert len(unlock) == 9
    assert unlock["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
