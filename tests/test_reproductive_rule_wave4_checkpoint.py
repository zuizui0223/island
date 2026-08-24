import hashlib
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.reproductive_rule_wave4_checkpoint import (
    SOURCE_EXCERPT,
    SOURCE_GROUP,
    SOURCE_LINEAGE,
    reviewed_rows,
)

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "reproductive_rule_wave4_checkpoint_20260814"
)


def test_rows_are_exact_species_direct_compatibility_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 6
    assert set(evidence["accepted_species"]) == {
        "Anacamptis papilionacea",
        "Himantoglossum adriaticum",
        "Ophrys apifera",
        "Ophrys bertolonii",
        "Orchis provincialis",
        "Serapias vomeracea",
    }
    assert "Orchis patens" not in set(evidence["accepted_species"])
    assert evidence["trait_name"].eq("self_incompatibility").all()
    assert evidence["normalized_value"].eq("SC").all()
    assert evidence["evidence_quality"].eq("high").all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_lineage"].eq(SOURCE_LINEAGE).all()
    assert evidence["source_lineage"].nunique() == 1
    assert not set(evidence["trait_name"]).intersection(
        {"autonomous_selfing_capacity", "mating_system", "cleistogamy"}
    )


def test_shared_excerpt_fingerprint_is_reproducible() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    expected = hashlib.sha256(SOURCE_EXCERPT.encode()).hexdigest()

    assert evidence["content_sha256"].eq(expected).all()
    assert evidence["content_sha256_basis"].eq(
        "verified_publisher_fulltext_excerpt_utf8_bytes"
    ).all()


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

    assert len(evidence) == len(audit) == len(accepted) == 1_911
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(wave) == 6
    assert wave["source_lineage"].nunique() == 1
    wave_audit = audit.loc[audit["candidate_id"].isin(wave["candidate_id"])]
    assert wave_audit["decision"].str.casefold().eq("accept").all()
    assert wave_audit["cultivar_contamination"].str.casefold().eq("false").all()
