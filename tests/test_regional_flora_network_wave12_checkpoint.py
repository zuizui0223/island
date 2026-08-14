from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import validate_individually_reviewed_evidence
from island_v2.regional_flora_wave11_checkpoint import reviewed_rows

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "regional_flora_network_wave12_checkpoint_20260814"
)
SOURCE = CHECKPOINT / "regional_flora_network_wave12_source_rows_20260814.csv"
SOURCE_GROUP = "regional_flora_network_wave12_checkpoint_20260814"


def test_full_extraction_audit_is_fail_closed() -> None:
    audit = pd.read_csv(
        CHECKPOINT / "regional_flora_network_extraction_audit_20260814.csv",
        dtype=str,
    ).fillna("")
    accepted = audit["review_decision"].eq("accept")

    assert len(audit) == 37
    assert int(accepted.sum()) == 31
    assert int(audit["review_decision"].eq("reject").sum()) == 6
    assert accepted.mean() == 31 / 37
    assert not audit["cultivar_contamination"].str.casefold().eq("true").any()
    assert audit["source_run_id"].eq("31767694551").all()
    assert audit["source_artifact_id"].eq("9207019490").all()
    assert audit["source_url"].str.startswith("https://").all()
    assert audit["source_excerpt"].ne("").all()
    assert audit["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()

    rejected = set(
        audit.loc[~accepted, ["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    )
    assert rejected == {
        ("Alternanthera nodiflora", "flower_size_class"),
        ("Euphorbia zambesiana", "flower_primary_color"),
        ("Nervilia bicarinata", "flower_size_class"),
        ("Dovyalis lucida", "flower_size_class"),
        ("Thalassodendron ciliatum", "flower_size_class"),
        ("Zanthoxylum delagoense", "flower_size_class"),
    }


def test_accepted_source_rows_are_exact_species_direct_medium() -> None:
    source = pd.read_csv(SOURCE, dtype=str).fillna("")
    evidence = pd.DataFrame(
        reviewed_rows(SOURCE, source_group=SOURCE_GROUP)
    ).fillna("")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")
    family = master.set_index("accepted_species")["family"].to_dict()

    assert len(source) == len(evidence) == 31
    assert source["accepted_species"].nunique() == 26
    assert all(
        family[species] == source_family
        for species, source_family in zip(
            source["accepted_species"], source["family"], strict=True
        )
    )
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["evidence_quality"].eq("medium").all()
    assert evidence["inference_rule"].eq("").all()
    assert not evidence["trait_name"].isin({"pollen_vector_mode", "reward_type"}).any()
    assert set(evidence["trait_name"]) == {
        "flower_primary_color",
        "floral_form",
        "floral_symmetry",
        "flower_size_class",
        "inflorescence_display",
    }


def test_combined_checkpoint_passes_the_common_review_gate() -> None:
    evidence = pd.read_csv(CHECKPOINT / "combined_curated_evidence_20260814.csv", dtype=str).fillna(
        ""
    )
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260814.csv", dtype=str
    ).fillna("")
    master = pd.read_csv("data/v2/staging/gbif/collected/island_taxa.csv", dtype=str).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence,
        audit,
        master,
        Path("config/trait_ontology.yml"),
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_749
    wave = accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]
    assert len(wave) == 31
    assert wave["candidate_id"].is_unique
