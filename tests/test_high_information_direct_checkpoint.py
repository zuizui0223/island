from pathlib import Path

import pandas as pd

from island_v2.high_information_direct_checkpoint import reviewed_rows
from island_v2.open_web_finalize import validate_individually_reviewed_evidence


def test_reproductive_records_do_not_substitute_traits() -> None:
    rows = reviewed_rows()
    aloe = next(row for row in rows if row["accepted_species"] == "Aloe maculata")
    begonia = next(
        row for row in rows if row["accepted_species"] == "Begonia nelumbiifolia"
    )
    abutilon = [row for row in rows if row["accepted_species"] == "Abutilon indicum"]
    pohnpei = [
        row
        for row in rows
        if row["source_provider"] == "yomai_williams_2021_primary_article"
    ]
    assert aloe["trait_name"] == "autonomous_selfing_capacity"
    assert aloe["normalized_value"] == "autonomous"
    assert begonia["trait_name"] == "mating_system"
    assert begonia["normalized_value"] == "mixed_mating"
    assert {(row["trait_name"], row["normalized_value"]) for row in abutilon} == {
        ("self_incompatibility", "SC"),
        ("autonomous_selfing_capacity", "autonomous"),
    }
    assert len({row["source_lineage"] for row in abutilon}) == 1
    assert len(pohnpei) == 18
    assert {row["trait_name"] for row in pohnpei} == {
        "self_incompatibility",
        "autonomous_selfing_capacity",
    }
    mimosa_autonomy = next(
        row
        for row in pohnpei
        if row["accepted_species"] == "Mimosa diplotricha"
        and row["trait_name"] == "autonomous_selfing_capacity"
    )
    assert mimosa_autonomy["normalized_value"] == "absent"
    jasmine = next(
        row for row in rows if row["accepted_species"] == "Jasminum auriculatum"
    )
    assert jasmine["trait_name"] == "self_incompatibility"
    assert jasmine["normalized_value"] == "SC"
    assert "not autonomous selfing" in jasmine["raw_value"]
    assert not {row["trait_name"] for row in rows}.intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_committed_high_information_checkpoint_passes_review_gate() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "next_acquisition_checkpoint_20260811"
    )
    evidence = pd.read_csv(
        root / "combined_curated_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "combined_curated_manual_audit_20260811.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(accepted) == 662
    assert evidence["candidate_id"].is_unique
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    added = evidence.loc[
        evidence["source_group"].eq("high_leverage_direct_checkpoint")
        & evidence["accepted_species"].isin(
            {
                "Ficus religiosa",
                "Erythroxylum citrifolium",
                "Ternstroemia penangiana",
                "Aloe maculata",
                "Begonia nelumbiifolia",
                "Abutilon indicum",
                "Crotalaria pallida",
                "Dendrolobium umbellatum",
                "Leucaena leucocephala",
                "Mimosa diplotricha",
                "Senna alata",
                "Vigna marina",
                "Hibiscus tiliaceus",
                "Sida acuta",
                "Thespesia populnea",
                "Claoxylon taitense",
                "Osmoxylon novoguineense",
                "Globba campsophylla",
                "Berberis darwinii",
                "Jasminum auriculatum",
                "Allophylus cobbe",
                "Dicliptera sexangularis",
            }
        )
    ]
    assert len(added) == 39
    assert added["source_excerpt"].ne("").all()
    assert added["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
