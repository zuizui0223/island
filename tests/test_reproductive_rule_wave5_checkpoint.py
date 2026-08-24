from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.reproductive_rule_wave5_checkpoint import (
    SOURCE_GROUP,
    SOURCES,
    reviewed_rows,
)

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "reproductive_rule_wave5_checkpoint_20260814"
)


def test_rows_are_exact_trait_species_direct_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 7
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 7
    assert set(evidence["trait_name"]) == {
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "mating_system",
    }
    assert set(evidence["evidence_quality"]) == {"high", "medium"}
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_lineage"].nunique() == 7
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_queue_potential_is_diagnostic_only() -> None:
    assert sum(int(source["potential"]) for source in SOURCES) == 463
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    assert evidence["inference_rule"].eq("").all()


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

    assert len(evidence) == len(audit) == len(accepted) == 1_918
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(wave) == 7
    assert wave["source_lineage"].nunique() == 7
    wave_audit = audit.loc[audit["candidate_id"].isin(wave["candidate_id"])]
    assert wave_audit["decision"].str.casefold().eq("accept").all()
    assert wave_audit["cultivar_contamination"].str.casefold().eq("false").all()
