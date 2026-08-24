from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.visible_morphology_wave6_checkpoint import (
    SOURCE_GROUP,
    SOURCES,
    reviewed_rows,
)

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "visible_morphology_wave6_checkpoint_20260814"
)


def test_rows_are_exact_trait_species_or_synonym_direct_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 11
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 11
    assert set(evidence["trait_name"]) == {
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
        "tube_depth_class",
    }
    assert set(evidence["evidence_quality"]) == {"high", "medium"}
    assert set(evidence["evidence_scope"]) == {"species_direct", "synonym_direct"}
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_lineage"].nunique() == 10
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_multistate_and_queue_potential_are_conservative() -> None:
    assert sum(int(source["potential"]) for source in SOURCES) == 146
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    cassipourea = evidence.loc[evidence["accepted_species"].eq("Cassipourea elliptica")]
    assert cassipourea.iloc[0]["normalized_value"] == "solitary|umbel_corymb"
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

    assert len(evidence) == len(audit) == len(accepted) == 1_929
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(wave) == 11
    assert wave["source_lineage"].nunique() == 10
    wave_audit = audit.loc[audit["candidate_id"].isin(wave["candidate_id"])]
    assert wave_audit["decision"].str.casefold().eq("accept").all()
    assert wave_audit["cultivar_contamination"].str.casefold().eq("false").all()
