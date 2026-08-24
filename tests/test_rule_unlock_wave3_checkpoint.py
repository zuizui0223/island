from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.rule_unlock_wave3_checkpoint import reviewed_rows

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave3_checkpoint_20260813"
)


def test_wave3_rows_are_small_trait_specific_and_source_backed() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 7
    assert evidence["accepted_species"].nunique() == 6
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 7
    assert evidence["candidate_id"].is_unique
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["source_excerpt"].str.len().gt(40).all()
    assert evidence["source_lineage"].nunique() == 5
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 4,
        "medium": 3,
    }
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_bidens_uses_only_species_named_in_bagging_result() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[evidence["source_lineage"].eq("doi:10.1111/wbm.12142")]

    assert set(rows["accepted_species"]) == {
        "Bidens bipinnata",
        "Bidens biternata",
        "Bidens frondosa",
    }
    assert rows["trait_name"].eq("autonomous_selfing_capacity").all()
    assert rows["normalized_value"].eq("autonomous").all()
    assert rows["source_excerpt"].str.contains("high seed-set in bagged capitula").all()
    assert "Bidens pilosa" not in set(rows["accepted_species"])


def test_morphology_and_reproduction_stay_separate() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    coccothrinax = evidence.loc[
        evidence["accepted_species"].eq("Coccothrinax barbadensis")
    ]
    assert set(zip(coccothrinax["trait_name"], coccothrinax["normalized_value"])) == {
        ("flower_size_class", "small"),
        ("inflorescence_display", "raceme_spike_panicle"),
    }
    assert coccothrinax["evidence_quality"].eq("medium").all()

    bambusa = evidence.loc[evidence["accepted_species"].eq("Bambusa tulda")].iloc[0]
    assert bambusa["trait_name"] == "inflorescence_display"
    assert bambusa["normalized_value"] == "raceme_spike_panicle"
    assert bambusa["evidence_quality"] == "high"

    bulbophyllum = evidence.loc[
        evidence["accepted_species"].eq("Bulbophyllum tseanum")
    ].iloc[0]
    assert bulbophyllum["trait_name"] == "self_incompatibility"
    assert bulbophyllum["normalized_value"] == "SI"
    assert bulbophyllum["evidence_quality"] == "medium"


def test_committed_wave3_combined_checkpoint_passes_review_gate() -> None:
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

    assert len(evidence) == len(audit) == len(accepted) == 1714
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    wave = accepted.loc[
        accepted["source_group"].eq("rule_unlock_wave3_checkpoint_20260813")
    ]
    assert len(wave) == 7
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
