from pathlib import Path

import pandas as pd

from island_v2.kudo_2022_alpine_reproduction_checkpoint import (
    TABLE_ROWS,
    reviewed_rows,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "kudo_2022_alpine_reproduction_checkpoint_20260813"
)
MASTER = Path("data/v2/staging/gbif/collected/island_taxa.csv")
SYNONYMS = Path(
    "data/v2/external/traits/rodger_2021_pollinator_contribution/"
    "rodger_2021_two_backbone_synonyms_20260811.csv"
)


def test_all_source_rows_are_audited_and_only_species_rank_is_accepted() -> None:
    source_audit = pd.read_csv(
        CHECKPOINT / "kudo_2022_full_table_audit_20260813.csv", dtype=str
    ).fillna("")

    assert len(TABLE_ROWS) == len(source_audit) == 46
    assert source_audit["decision"].value_counts().to_dict() == {
        "accept": 23,
        "reject": 23,
    }
    rejected = source_audit.loc[source_audit["decision"].eq("reject")]
    assert rejected["decision_reason"].value_counts().to_dict() == {
        "rank_below_species_rejected": 8,
        "no_strict_target_master_match": 7,
        "no_trait_statement": 5,
        "family_conflict_rejected": 2,
        "source_name_typo_not_fuzzily_matched": 1,
    }


def test_trait_categories_are_not_substituted_across_reproductive_traits() -> None:
    evidence = pd.DataFrame(reviewed_rows(MASTER, SYNONYMS)).fillna("")

    assert len(evidence) == 23
    assert evidence["candidate_id"].is_unique
    assert evidence["evidence_quality"].eq("high").all()
    assert evidence["source_lineage"].nunique() == 1
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["source_excerpt"].str.contains("\t", regex=False).all()
    assert evidence["trait_name"].value_counts().to_dict() == {
        "mating_system": 15,
        "self_incompatibility": 7,
        "autonomous_selfing_capacity": 1,
    }

    autonomous = evidence.loc[
        evidence["accepted_species"].eq("Saussurea yanagisawae")
    ].iloc[0]
    assert autonomous["trait_name"] == "autonomous_selfing_capacity"
    assert autonomous["normalized_value"] == "autonomous"

    incompatible = evidence.loc[
        evidence["accepted_species"].eq("Arnica unalaschcensis")
    ].iloc[0]
    assert incompatible["trait_name"] == "self_incompatibility"
    assert incompatible["normalized_value"] == "SI"


def test_only_frozen_two_backbone_synonym_is_used() -> None:
    evidence = pd.DataFrame(reviewed_rows(MASTER, SYNONYMS)).fillna("")
    synonym = evidence.loc[
        evidence["name_match_method"].eq("exact_synonym")
    ]

    assert len(synonym) == 1
    row = synonym.iloc[0]
    assert row["accepted_species"] == "Kalmia procumbens"
    assert row["matched_page_name"] == "Loiseleuria procumbens"
    assert row["evidence_scope"] == "synonym_direct"
    assert row["trait_name"] == "self_incompatibility"
    assert row["normalized_value"] == "SI"


def test_committed_kudo_combined_checkpoint_passes_review_gate() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "combined_curated_evidence_20260813.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260813.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(MASTER, dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(audit) == len(accepted) == 1737
    assert evidence["candidate_id"].is_unique
    wave = accepted.loc[
        accepted["source_group"].eq(
            "kudo_2022_alpine_reproduction_checkpoint_20260813"
        )
    ]
    assert len(wave) == 23
    assert wave["evidence_scope"].value_counts().to_dict() == {
        "species_direct": 22,
        "synonym_direct": 1,
    }
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
