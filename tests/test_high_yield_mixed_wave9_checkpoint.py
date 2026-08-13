from pathlib import Path

import pandas as pd

from island_v2.high_yield_mixed_wave9_checkpoint import (
    SOURCE_GROUP,
    SOURCES,
    reviewed_rows,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "high_yield_mixed_wave9_checkpoint_20260814"
)


def test_rows_are_exact_trait_specific_direct_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 9
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 9
    assert set(evidence["trait_name"]) == {
        "autonomous_selfing_capacity",
        "floral_symmetry",
        "mating_system",
        "self_incompatibility",
    }
    assert set(evidence["evidence_scope"]) == {"species_direct"}
    assert set(evidence["language"]) == {"en"}
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_lineage"].nunique() == 6
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["inference_rule"].eq("").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_queue_potential_is_diagnostic_only() -> None:
    assert sum(int(source["potential"]) for source in SOURCES) == 84


def test_category_errors_are_avoided() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    values = {
        (row.accepted_species, row.trait_name): row.normalized_value
        for row in evidence.itertuples()
    }
    assert values[("Callistemon phoeniceus", "autonomous_selfing_capacity")] == (
        "mixed_or_variable"
    )
    assert values[("Triumfetta rhomboidea", "autonomous_selfing_capacity")] == (
        "delayed"
    )
    assert values[("Melaleuca quinquenervia", "mating_system")] == "mixed_mating"
    assert values[("Melaleuca quinquenervia", "self_incompatibility")] == "SC"


def test_strict_synonym_identity_is_explicit() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    row = evidence.loc[evidence["accepted_species"].eq("Sporobolus alterniflorus")].iloc[0]
    assert row["matched_page_name"] == "Spartina alterniflora"
    assert row["name_match_method"] == "synonym_exact"
    assert row["name_resolution_lineage"] == "powo_strict_synonym_to_master_accepted"


def test_committed_combined_checkpoint_passes_individual_review_gate() -> None:
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

    assert len(evidence) == len(audit) == len(accepted) == 1_952
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(wave) == 9
    assert wave["source_lineage"].nunique() == 6
    wave_audit = audit.loc[audit["candidate_id"].isin(wave["candidate_id"])]
    assert wave_audit["decision"].str.casefold().eq("accept").all()
    assert wave_audit["cultivar_contamination"].str.casefold().eq("false").all()
