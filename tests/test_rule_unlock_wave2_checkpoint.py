from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import (
    validate_individually_reviewed_evidence,
)
from island_v2.rule_unlock_wave2_checkpoint import reviewed_rows

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "rule_unlock_wave2_checkpoint_20260812"
)


def test_wave2_rows_are_trait_specific_and_source_backed() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert len(evidence) == 36
    assert evidence["accepted_species"].nunique() == 24
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 36
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["source_excerpt"].str.len().gt(30).all()
    assert evidence["source_lineage"].nunique() == 22
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 28,
        "medium": 8,
    }

    reproductive = evidence.loc[evidence["axis"].eq("reproductive_assurance")]
    assert len(reproductive) == 23
    assert set(reproductive["trait_name"]) == {
        "autonomous_selfing_capacity",
        "mating_system",
        "self_incompatibility",
    }
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_zero_bagged_set_does_not_overwrite_hand_self_compatibility() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    for species in ["Hosta ventricosa", "Melastoma malabathricum"]:
        rows = evidence.loc[evidence["accepted_species"].eq(species)]
        observed = set(zip(rows["trait_name"], rows["normalized_value"]))
        assert ("autonomous_selfing_capacity", "absent") in observed
        assert ("self_incompatibility", "SC") in observed


def test_morphology_rows_keep_explicit_species_descriptions_trait_specific() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    morphology = evidence.loc[
        evidence["trait_name"].eq("floral_symmetry")
    ].set_index("accepted_species")

    assert set(morphology.index) == {
        "Adenia cissampeloides",
        "Alangium salviifolium",
        "Corchorus olitorius",
        "Drypetes assamica",
        "Melastoma malabathricum",
        "Sideritis canariensis",
        "Tristaniopsis laurina",
    }
    assert morphology["normalized_value"].to_dict() == {
        "Adenia cissampeloides": "actinomorphic",
        "Alangium salviifolium": "actinomorphic",
        "Corchorus olitorius": "actinomorphic",
        "Drypetes assamica": "actinomorphic",
        "Melastoma malabathricum": "actinomorphic",
        "Sideritis canariensis": "zygomorphic",
        "Tristaniopsis laurina": "actinomorphic",
    }
    assert morphology.loc[
        ["Adenia cissampeloides", "Tristaniopsis laurina"],
        "evidence_quality",
    ].eq("high").all()
    assert morphology.loc[
        [
            "Alangium salviifolium",
            "Corchorus olitorius",
            "Drypetes assamica",
            "Melastoma malabathricum",
            "Sideritis canariensis",
        ],
        "evidence_quality",
    ].eq("medium").all()

    tube_depth = evidence.loc[
        evidence["accepted_species"].eq("Palaquium obovatum")
        & evidence["trait_name"].eq("tube_depth_class")
    ]
    assert len(tube_depth) == 1
    assert tube_depth.iloc[0]["normalized_value"] == "shallow"
    assert tube_depth.iloc[0]["evidence_quality"] == "high"

    rule_unlocks = evidence.loc[
        evidence["accepted_species"].isin(
            {"Polycarpaea corymbosa", "Boronia muelleri", "Benstonea foetida"}
        )
    ].set_index("accepted_species")
    assert rule_unlocks[["trait_name", "normalized_value"]].to_dict("index") == {
        "Polycarpaea corymbosa": {
            "trait_name": "inflorescence_display",
            "normalized_value": "umbel_corymb",
        },
        "Boronia muelleri": {
            "trait_name": "floral_form",
            "normalized_value": "open_radial",
        },
        "Benstonea foetida": {
            "trait_name": "inflorescence_display",
            "normalized_value": "raceme_spike_panicle",
        },
    }


def test_tristaniopsis_bagging_is_not_converted_to_self_incompatibility() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[evidence["accepted_species"].eq("Tristaniopsis laurina")]
    observed = set(zip(rows["trait_name"], rows["normalized_value"]))

    assert observed == {
        ("autonomous_selfing_capacity", "absent"),
        ("floral_symmetry", "actinomorphic"),
    }
    assert "self_incompatibility" not in set(rows["trait_name"])


def test_commelineae_thesis_keeps_si_and_mating_system_as_separate_traits() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    species = {
        "Commelina diffusa",
        "Dictyospermum montanum",
        "Floscopa scandens",
        "Murdannia nudiflora",
        "Rhopalephora scaberrima",
    }
    rows = evidence.loc[evidence["accepted_species"].isin(species)]

    assert len(rows.loc[rows["source_lineage"].eq(
        "study:veena-2020:commelineae-controlled-pollination"
    )]) == 10
    for accepted_species in species:
        observed = set(
            zip(
                rows.loc[rows["accepted_species"].eq(accepted_species), "trait_name"],
                rows.loc[rows["accepted_species"].eq(accepted_species), "normalized_value"],
            )
        )
        assert ("self_incompatibility", "SC") in observed
        assert ("mating_system", "mixed_mating") in observed


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

    assert len(evidence) == len(audit) == len(accepted) == 821
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    wave = accepted.loc[accepted["source_group"].eq(
        "rule_unlock_wave2_checkpoint_20260812"
    )]
    assert len(wave) == 36


def test_next_rule_unlocks_are_explicit_and_do_not_cross_traits() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["accepted_species"].isin(
            {"Pleioluma balansana", "Kunzea ericoides", "Corchorus olitorius"}
        )
    ].set_index("accepted_species")

    assert rows[["trait_name", "normalized_value", "evidence_quality"]].to_dict(
        "index"
    ) == {
        "Pleioluma balansana": {
            "trait_name": "flower_primary_color",
            "normalized_value": "white",
            "evidence_quality": "high",
        },
        "Kunzea ericoides": {
            "trait_name": "floral_form",
            "normalized_value": "open_radial",
            "evidence_quality": "medium",
        },
        "Corchorus olitorius": {
            "trait_name": "floral_symmetry",
            "normalized_value": "actinomorphic",
            "evidence_quality": "medium",
        },
    }
    assert not set(rows["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )
