import hashlib
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.reproductive_counterevidence_checkpoint import (
    SOURCE_GROUP,
    reviewed_rows,
)

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/reproductive_counterevidence_checkpoint_20260813"
)


def test_rows_keep_exact_reproductive_traits_and_counterevidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 9
    assert evidence[["accepted_species", "trait_name"]].duplicated().sum() == 0
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 5,
        "medium": 4,
    }
    assert evidence["trait_name"].value_counts().to_dict() == {
        "autonomous_selfing_capacity": 7,
        "mating_system": 2,
    }
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type", "self_incompatibility"}
    )

    allium = evidence.loc[evidence["accepted_species"].str.startswith("Allium ")]
    assert set(zip(allium["accepted_species"], allium["normalized_value"], strict=True)) == {
        ("Allium cernuum", "autonomous"),
        ("Allium oleraceum", "absent"),
    }
    liparis = evidence.loc[evidence["accepted_species"].str.startswith("Liparis ")]
    assert set(zip(liparis["accepted_species"], liparis["normalized_value"], strict=True)) == {
        ("Liparis liliifolia", "absent"),
        ("Liparis loeselii", "autonomous"),
    }


def test_synonym_and_lineage_contracts_are_explicit() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    calanthe = evidence.loc[evidence["accepted_species"].eq("Calanthe striata")]
    assert calanthe.iloc[0]["matched_page_name"] == "Calanthe sieboldii"
    assert calanthe.iloc[0]["name_match_method"] == "exact_synonym"
    assert calanthe.iloc[0]["normalized_value"] == "absent"

    hakea = evidence.loc[evidence["accepted_species"].eq("Hakea carinata")]
    assert hakea.iloc[0]["trait_name"] == "mating_system"
    assert hakea.iloc[0]["source_lineage"] == "doi:10.1071/BT97123"

    abutilon = evidence.loc[evidence["accepted_species"].eq("Abutilon theophrasti")]
    assert abutilon.iloc[0]["evidence_quality"] == "medium"
    assert abutilon.iloc[0]["source_lineage"].startswith("citation:andersen_1988")


def test_excerpt_backed_fingerprints_match_the_stored_quote() -> None:
    for row in reviewed_rows():
        if row["content_sha256_basis"] != "verified_exact_supporting_excerpt_utf8_bytes":
            continue
        assert hashlib.sha256(row["source_excerpt"].encode()).hexdigest() == row["content_sha256"]


def test_committed_combined_checkpoint_passes_individual_review_gate() -> None:
    evidence = pd.read_csv(CHECKPOINT / "combined_curated_evidence_20260813.csv", dtype=str).fillna(
        ""
    )
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260813.csv", dtype=str
    ).fillna("")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(audit) == len(accepted) == 1_875
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    recovered = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(recovered) == 9
    assert recovered["trait_name"].isin({"autonomous_selfing_capacity", "mating_system"}).all()
    recovered_audit = audit.loc[audit["candidate_id"].isin(recovered["candidate_id"])]
    assert recovered_audit["decision"].str.casefold().eq("accept").all()
    assert recovered_audit["cultivar_contamination"].str.casefold().eq("false").all()
