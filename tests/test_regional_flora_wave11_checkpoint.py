from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.regional_flora_wave11_checkpoint import SOURCE_GROUP, reviewed_rows

CHECKPOINT = Path("data/v2/staging/traits/open_web_pilot/regional_flora_wave11_checkpoint_20260814")
SOURCE_ROWS = CHECKPOINT / "regional_flora_wave11_source_rows_20260814.csv"


def test_wave11_is_reviewed_species_direct_evidence() -> None:
    evidence = pd.DataFrame(reviewed_rows(SOURCE_ROWS)).fillna("")

    assert len(evidence) == 759
    assert evidence["accepted_species"].nunique() == 729
    assert (
        evidence[["accepted_species", "trait_name", "source_lineage"]].drop_duplicates().shape[0]
        == 759
    )
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["name_match_method"].eq("accepted_name_exact").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert set(evidence["trait_name"]) == {
        "flower_primary_color",
        "floral_form",
        "flower_size_class",
        "inflorescence_display",
        "tube_depth_class",
    }
    assert not evidence["trait_name"].isin({"pollen_vector_mode", "reward_type"}).any()


def test_provider_quality_and_lineage_contract() -> None:
    evidence = pd.DataFrame(reviewed_rows(SOURCE_ROWS)).fillna("")

    assert evidence["source_provider"].value_counts().to_dict() == {
        "Endemia New Caledonia": 364,
        "New Zealand Plant Conservation Network": 283,
        "Flora of Zimbabwe": 112,
    }
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "medium": 476,
        "high": 283,
    }
    assert evidence["source_lineage"].nunique() == 729

    endemia = evidence.loc[evidence["source_provider"].eq("Endemia New Caledonia")]
    assert endemia["source_tier"].eq("B").all()
    assert endemia["evidence_quality"].eq("medium").all()
    assert endemia["source_url"].str.contains("/en/search?").all()
    assert endemia["source_citation"].str.contains("species fiche").all()

    nzpcn = evidence.loc[evidence["source_provider"].eq("New Zealand Plant Conservation Network")]
    assert nzpcn["source_tier"].eq("A").all()
    assert nzpcn["evidence_quality"].eq("high").all()
    assert (
        nzpcn["source_excerpt"]
        .str.replace("\r\n", "\n", regex=False)
        .str.startswith("Flower colours\n")
        .all()
    )
    assert nzpcn["wild_cultivated_cultivar_status"].str.startswith("native_species").all()


def test_fixed_master_species_and_family_are_exact() -> None:
    source = pd.read_csv(SOURCE_ROWS, dtype=str).fillna("")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")
    identity = master.set_index("accepted_species")["family"].to_dict()

    assert set(source["accepted_species"]).issubset(identity)
    assert all(
        identity[species] == family
        for species, family in zip(source["accepted_species"], source["family"], strict=True)
    )


def test_committed_combined_checkpoint_passes_review_gate() -> None:
    evidence = pd.read_csv(CHECKPOINT / "combined_curated_evidence_20260814.csv", dtype=str).fillna(
        ""
    )
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260814.csv", dtype=str
    ).fillna("")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_718
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique

    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    wave_audit = audit.loc[audit["candidate_id"].isin(wave["candidate_id"])]
    assert len(wave) == len(wave_audit) == 759
    assert wave_audit["decision"].str.casefold().eq("accept").all()
    assert wave_audit["species_identity_correct"].str.casefold().eq("true").all()
    assert wave_audit["value_correct"].str.casefold().eq("true").all()
    assert wave_audit["provenance_complete"].str.casefold().eq("true").all()
    assert wave_audit["cultivar_contamination"].str.casefold().eq("false").all()
