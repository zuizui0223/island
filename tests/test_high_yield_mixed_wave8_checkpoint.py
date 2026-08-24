from pathlib import Path

import pandas as pd

from island_v2.high_yield_mixed_wave8_checkpoint import (
    SOURCE_GROUP,
    SOURCES,
    reviewed_rows,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "high_yield_mixed_wave8_checkpoint_20260814"
)


def test_rows_are_exact_trait_specific_direct_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 7
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 7
    assert set(evidence["trait_name"]) == {
        "flower_primary_color",
        "flower_size_class",
        "mating_system",
        "self_incompatibility",
    }
    assert set(evidence["evidence_scope"]) == {"species_direct"}
    assert set(evidence["language"]) == {"en", "pt"}
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["source_lineage"].nunique() == 7
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["inference_rule"].eq("").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_queue_potential_is_diagnostic_only() -> None:
    assert sum(int(source["potential"]) for source in SOURCES) == 60


def test_pimenta_pdf_row_is_exact_species_and_flower_colour() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    pimenta = evidence.loc[
        evidence["accepted_species"].eq("Pimenta pseudocaryophyllus")
    ].iloc[0]
    assert pimenta["trait_name"] == "flower_primary_color"
    assert pimenta["normalized_value"] == "white"
    assert "flores brancas" in pimenta["source_excerpt"]
    assert pimenta["content_sha256_basis"] == "downloaded_official_pdf_bytes"


def test_counterevidence_is_preserved_not_cherry_picked() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    states = {
        (row.accepted_species, row.trait_name): row.normalized_value
        for row in evidence.itertuples()
    }
    assert states[("Lyonia lucida", "self_incompatibility")] == "SI"
    assert states[("Anthoxanthum odoratum", "self_incompatibility")] == "SI"
    assert states[("Oxytropis campestris", "self_incompatibility")] == "SC"
    oxytropis = evidence.loc[
        evidence["accepted_species"].eq("Oxytropis campestris")
    ].iloc[0]
    assert oxytropis["source_lineage"] == "isbn:978-1-63482-200-8:chapter3"


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

    assert len(evidence) == len(audit) == len(accepted) == 1_943
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(wave) == 7
    assert wave["source_lineage"].nunique() == 7
    wave_audit = audit.loc[audit["candidate_id"].isin(wave["candidate_id"])]
    assert wave_audit["decision"].str.casefold().eq("accept").all()
    assert wave_audit["cultivar_contamination"].str.casefold().eq("false").all()
