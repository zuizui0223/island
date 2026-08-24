from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import (
    validate_individually_reviewed_evidence,
)
from island_v2.pohnpei_reproductive_checkpoint import reviewed_rows

SOURCE_TABLE = Path(
    "data/v2/external/traits/yomai_2021_pohnpei/hand_pollination_table_2_1.csv"
)
CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "pohnpei_reproductive_checkpoint_20260812"
)


def test_source_lineage_and_trait_boundaries() -> None:
    source = pd.read_csv(SOURCE_TABLE, dtype=str, keep_default_na=False)
    evidence = pd.DataFrame(reviewed_rows(source)).fillna("")

    assert len(source) == 36
    assert source["include_independent_chapter2"].eq("true").sum() == 25
    assert source["prior_article_lineage"].eq("true").sum() == 11
    assert len(evidence) == 49
    assert evidence["accepted_species"].nunique() == 25
    assert evidence["trait_name"].value_counts().to_dict() == {
        "self_incompatibility": 25,
        "autonomous_selfing_capacity": 24,
    }
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type", "mating_system"}
    )
    assert evidence["source_lineage"].eq(
        "repository:trace-tennessee:utk_graddiss:6964:chapter2"
    ).all()
    assert evidence["source_excerpt"].str.contains("Table 2.", regex=False).all()

    bidens = evidence.loc[evidence["accepted_species"].eq("Bidens pilosa")]
    assert set(bidens["trait_name"]) == {"self_incompatibility"}
    assert set(
        evidence.loc[
            evidence["accepted_species"].isin(
                {"Crotalaria pallida", "Dendrolobium umbellatum", "Vigna marina"}
            ),
            "accepted_species",
        ]
    ) == set()


def test_identity_and_autofertility_states_are_fail_closed() -> None:
    source = pd.read_csv(SOURCE_TABLE, dtype=str, keep_default_na=False)
    evidence = pd.DataFrame(reviewed_rows(source)).fillna("")

    synonyms = evidence.loc[
        evidence["accepted_species"].isin({"Grona heterocarpos", "Ludwigia octovalvis"})
    ]
    assert synonyms["evidence_scope"].eq("synonym_direct").all()
    assert synonyms["name_match_method"].eq("synonym_exact").all()
    same_document = evidence.loc[
        evidence["accepted_species"].isin(
            {"Centrosema molle", "Dalbergia candenatensis"}
        )
    ]
    assert same_document["name_match_method"].eq("accepted_name_exact").all()
    assert same_document["name_resolution_lineage"].str.startswith(
        "dissertation_table_1_species_list_"
    ).all()

    autonomy = evidence.loc[
        evidence["trait_name"].eq("autonomous_selfing_capacity")
    ].set_index("accepted_species")
    assert autonomy.at["Acacia confusa", "normalized_value"] == "autonomous"
    assert autonomy.at["Ipomoea indica", "normalized_value"] == "mixed_or_variable"
    assert autonomy.at["Senna occidentalis", "normalized_value"] == (
        "mixed_or_variable"
    )
    assert autonomy["normalized_value"].value_counts().to_dict() == {
        "autonomous": 20,
        "mixed_or_variable": 4,
    }


def test_committed_combined_checkpoint_passes_review_gate() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "combined_curated_evidence_20260812.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260812.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(audit) == len(accepted) == 776
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    pohnpei = accepted.loc[
        accepted["source_group"].eq("pohnpei_reproductive_checkpoint_20260812")
    ]
    assert len(pohnpei) == 49
    assert pohnpei["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
