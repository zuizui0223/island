from pathlib import Path

import pandas as pd

from island_v2.independent_lineage_recovery_checkpoint import (
    SOURCE_GROUP,
    reviewed_rows,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "independent_lineage_recovery_checkpoint_20260813"
)


def test_rows_are_direct_trait_specific_independent_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 20
    assert evidence[["accepted_species", "trait_name"]].duplicated().sum() == 0
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 16,
        "medium": 4,
    }
    assert evidence["trait_name"].value_counts().to_dict() == {
        "autonomous_selfing_capacity": 6,
        "mating_system": 7,
        "self_incompatibility": 2,
        "floral_symmetry": 2,
        "flower_size_class": 2,
        "inflorescence_display": 1,
    }
    assert evidence["source_lineage"].nunique() == 14
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_synonym_identity_and_underlying_lineage_are_explicit() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    erigeron = evidence.loc[evidence["accepted_species"].eq("Erigeron canadensis")]
    assert erigeron.iloc[0]["matched_page_name"] == "Conyza (=Erigeron) canadensis"
    assert erigeron.iloc[0]["name_match_method"] == "exact_synonym"

    lactuca = evidence.loc[evidence["accepted_species"].eq("Lactuca canadensis")]
    assert lactuca.iloc[0]["source_lineage"].startswith("citation:parrish_bazzaz_1979")

    new_erigeron = evidence.loc[
        evidence["accepted_species"].isin(
            {"Erigeron bonariensis", "Erigeron sumatrensis"}
        )
        & evidence["trait_name"].eq("autonomous_selfing_capacity")
    ]
    assert len(new_erigeron) == 2
    assert new_erigeron["name_match_method"].eq("exact_synonym").all()
    assert new_erigeron["source_lineage"].nunique() == 1


def test_multi_species_lineage_and_begonia_trait_separation_are_preserved() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    ficus = evidence.loc[
        evidence["accepted_species"].str.startswith("Ficus ")
        & evidence["trait_name"].eq("mating_system")
    ]
    assert len(ficus) == 5
    assert ficus["source_lineage"].eq("doi:10.1038/35607").all()
    assert ficus["normalized_value"].eq("predominantly_outcrossing").all()

    begonia = evidence.loc[evidence["accepted_species"].eq("Begonia integerrima")]
    assert set(zip(begonia["trait_name"], begonia["normalized_value"], strict=True)) == {
        ("self_incompatibility", "SI"),
    }
    assert "autonomous_selfing_capacity" not in set(begonia["trait_name"])
    assert begonia["source_lineage"].nunique() == 1

    casearia = evidence.loc[
        evidence["accepted_species"].eq("Casearia grandiflora")
    ]
    assert set(casearia["trait_name"]) == {
        "floral_symmetry",
        "flower_size_class",
    }
    assert casearia["source_lineage"].nunique() == 1
    assert casearia["content_sha256"].nunique() == 1


def test_bad_lolium_rigidum_row_is_explicitly_excluded() -> None:
    exclusions = pd.read_csv(
        "data/v2/staging/traits/open_web_pilot/"
        "direct_evidence_exclusions_20260811.csv",
        dtype=str,
    ).fillna("")
    match = exclusions.loc[
        exclusions["accepted_species"].eq("Lolium rigidum")
        & exclusions["trait_name"].eq("autonomous_selfing_capacity")
        & exclusions["normalized_value"].eq("autonomous")
    ]
    assert len(match) == 1
    assert "Avena sterilis" in match.iloc[0]["reason"]


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

    assert len(evidence) == len(audit) == len(accepted) == 1_766
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    recovered = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(recovered) == 20
    assert recovered["source_lineage"].nunique() == 14
    recovered_audit = audit.loc[audit["candidate_id"].isin(recovered["candidate_id"])]
    assert recovered_audit["reviewed_at_utc"].str.fullmatch(
        r"2026-08-13T(?:03:59:37|09:42:00)Z"
    ).all()
