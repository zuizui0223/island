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

    assert len(evidence) == 404
    assert evidence["accepted_species"].nunique() == 364
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 404
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["source_excerpt"].str.len().gt(30).all()
    assert evidence["source_lineage"].nunique() == 326
    assert evidence["evidence_quality"].value_counts().to_dict() == {
        "high": 46,
        "medium": 358,
    }

    reproductive = evidence.loc[evidence["axis"].eq("reproductive_assurance")]
    assert len(reproductive) == 75
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

    expected = {
        "Adenia cissampeloides",
        "Alangium salviifolium",
        "Celtis sinensis",
        "Corchorus olitorius",
        "Drypetes assamica",
        "Melastoma malabathricum",
        "Sideritis canariensis",
        "Suregada lanceolata",
        "Tristaniopsis laurina",
    }
    assert expected.issubset(set(morphology.index))
    assert morphology.loc[sorted(expected), "normalized_value"].to_dict() == {
        "Adenia cissampeloides": "actinomorphic",
        "Alangium salviifolium": "actinomorphic",
        "Celtis sinensis": "actinomorphic",
        "Corchorus olitorius": "actinomorphic",
        "Drypetes assamica": "actinomorphic",
        "Melastoma malabathricum": "actinomorphic",
        "Sideritis canariensis": "zygomorphic",
        "Suregada lanceolata": "actinomorphic",
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
            "Suregada lanceolata",
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

    assert len(evidence) == len(audit) == len(accepted) == 1189
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert audit["decision"].str.casefold().eq("accept").all()
    assert audit["cultivar_contamination"].str.casefold().eq("false").all()
    wave = accepted.loc[accepted["source_group"].eq(
        "rule_unlock_wave2_checkpoint_20260812"
    )]
    assert len(wave) == 404


def test_iisc_full_candidate_audit_is_precise_and_fail_closed() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[evidence["domain"].eq("indiaflora-ces.iisc.ac.in")]
    full_audit = pd.read_csv(
        CHECKPOINT / "india_flora_online_full_candidate_audit_20260812.csv",
        dtype=str,
    ).fillna("")

    assert len(full_audit) == 262
    assert full_audit["decision"].value_counts().to_dict() == {
        "accept": 250,
        "reject": 12,
    }
    assert len(rows) == 250
    assert rows["accepted_species"].nunique() == 232
    assert rows["evidence_quality"].eq("medium").all()
    assert rows["source_tier"].eq("A").all()
    assert rows["source_lineage"].nunique() == 232
    assert set(rows["trait_name"]) == {
        "flower_primary_color",
        "floral_form",
        "flower_size_class",
        "inflorescence_display",
    }
    assert full_audit["decision"].eq("accept").mean() >= 0.95
    assert full_audit["cultivar_contamination"].eq("true").mean() <= 0.02

    rejected_rows = full_audit.loc[full_audit["decision"].eq("reject")]
    rejected = set(rejected_rows["accepted_species"])
    rejected_keys = set(
        zip(rejected_rows["accepted_species"], rejected_rows["trait_name"])
    )
    accepted_keys = set(zip(rows["accepted_species"], rows["trait_name"]))
    assert {
        "Cordia grandis",
        "Croton tiglium",
        "Dinochloa andamanica",
        "Licuala peltata",
        "Nicotiana sanderae",
        "Spiraea bumalda",
    }.issubset(rejected)
    assert rejected_keys.isdisjoint(accepted_keys)


def test_bfis_bulk_snapshot_adds_only_explicit_novel_symmetry_records() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["content_sha256"].eq(
            "80f4d97075ebb66d41b5707b1ae4bc330dd41ea4df567b4e"
            "a19f63531aa88b4f"
        )
    ]

    assert len(rows) == rows["accepted_species"].nunique() == 49
    assert set(rows["trait_name"]) == {"floral_symmetry"}
    assert set(rows["normalized_value"]) == {"actinomorphic", "zygomorphic"}
    assert rows["evidence_quality"].eq("medium").all()
    assert rows["domain"].eq("bfis.bforest.gov.bd").all()
    assert rows["source_lineage"].nunique() == 49
    assert rows["source_excerpt"].str.contains(
        r"Species: .+; floral symmetry: (?:actinomorphic|zygomorphic)\.",
        regex=True,
    ).all()
    assert rows.loc[
        rows["accepted_species"].eq("Sophora wightii"), "normalized_value"
    ].tolist() == ["zygomorphic"]


def test_seventh_increment_excludes_conflicted_galapagos_row() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    assert "Chiococca alba" not in set(evidence["accepted_species"])
    galapagos = evidence.loc[
        evidence["source_url"].eq(
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC3489146/"
        )
    ]
    assert len(galapagos) == 13
    assert set(galapagos["evidence_quality"]) == {"medium"}
    assert set(galapagos["trait_name"]) == {
        "autonomous_selfing_capacity",
        "self_incompatibility",
    }


def test_encyclia_selfing_and_exclusion_keep_distinct_traits() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[evidence["accepted_species"].str.startswith("Encyclia ")]
    observed = set(
        zip(
            rows["accepted_species"],
            rows["trait_name"],
            rows["normalized_value"],
            rows["evidence_quality"],
        )
    )

    assert observed == {
        ("Encyclia phoenicea", "self_incompatibility", "SC", "medium"),
        ("Encyclia plicata", "self_incompatibility", "SC", "medium"),
        (
            "Encyclia tampensis",
            "autonomous_selfing_capacity",
            "absent",
            "high",
        ),
    }
    assert not rows.loc[rows["accepted_species"].eq("Encyclia tampensis")][
        "trait_name"
    ].eq("self_incompatibility").any()


def test_southern_ocean_compilation_uses_one_conservative_lineage() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["source_lineage"].eq(
            "compilation:doi:10.1093/aobpla/plv095:table-s1"
        )
    ]

    assert len(rows) == 30
    assert set(rows["trait_name"]) == {"self_incompatibility"}
    assert set(rows["evidence_quality"]) == {"medium"}
    partial = rows.loc[rows["accepted_species"].eq("Colobanthus affinis")]
    assert partial["normalized_value"].tolist() == ["mixed_or_variable"]


def test_juan_fernandez_rows_do_not_substitute_reproductive_traits() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[evidence["source_lineage"].eq("doi:10.2307/2657013")]

    observed = set(
        zip(
            rows["accepted_species"],
            rows["trait_name"],
            rows["normalized_value"],
        )
    )
    assert observed == {
        ("Berberis corymbosa", "self_incompatibility", "SI"),
        ("Wahlenbergia berteroi", "autonomous_selfing_capacity", "autonomous"),
        ("Wahlenbergia fernandeziana", "autonomous_selfing_capacity", "delayed"),
        ("Escallonia callcottiae", "autonomous_selfing_capacity", "autonomous"),
        ("Escallonia callcottiae", "mating_system", "mixed_mating"),
    }
    assert set(rows["evidence_quality"]) == {"high"}


def test_ninth_increment_keeps_visible_traits_species_direct() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["accepted_species"].isin(
            {
                "Pyrostria commersonii",
                "Quintinia acutifolia",
                "Suregada lanceolata",
            }
        )
        & evidence["source_record_id"].str.contains(
            "academie-reunion|vibrant-earth|suregada-lanceolata",
            regex=True,
        )
    ].set_index("accepted_species")

    assert rows[["trait_name", "normalized_value", "evidence_quality"]].to_dict(
        "index"
    ) == {
        "Pyrostria commersonii": {
            "trait_name": "flower_primary_color",
            "normalized_value": "white|yellow_orange",
            "evidence_quality": "medium",
        },
        "Quintinia acutifolia": {
            "trait_name": "flower_primary_color",
            "normalized_value": "white",
            "evidence_quality": "medium",
        },
        "Suregada lanceolata": {
            "trait_name": "floral_symmetry",
            "normalized_value": "actinomorphic",
            "evidence_quality": "medium",
        },
    }
    assert rows["evidence_scope"].eq("species_direct").all()
    assert not set(rows["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )
    assert rows.loc[
        "Quintinia acutifolia", "wild_cultivated_cultivar_status"
    ] == "species_level_horticultural_record_not_cultivar_limited"


def test_sixth_increment_unlocks_only_exact_traits_and_preserves_colour_sets() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")

    morphology = evidence.loc[
        evidence["accepted_species"].isin(
            {"Turraea cadetii", "Mosiera yamaniguensis"}
        )
    ]
    observed = set(
        zip(
            morphology["accepted_species"],
            morphology["trait_name"],
            morphology["normalized_value"],
        )
    )
    assert observed == {
        ("Turraea cadetii", "inflorescence_display", "umbel_corymb"),
        ("Turraea cadetii", "flower_primary_color", "white"),
        ("Mosiera yamaniguensis", "inflorescence_display", "solitary"),
        ("Mosiera yamaniguensis", "flower_primary_color", "white"),
    }

    colour = evidence.loc[
        evidence["accepted_species"].isin(
            {
                "Dichaetanthera rutenbergiana",
                "Acropogon bullatus",
                "Acropogon mesophilus",
                "Acropogon veillonii",
            }
        )
    ].set_index("accepted_species")
    assert colour["trait_name"].eq("flower_primary_color").all()
    assert colour["normalized_value"].to_dict() == {
        "Dichaetanthera rutenbergiana": (
            "green_brown_inconspicuous|yellow_orange"
        ),
        "Acropogon bullatus": "red_pink|yellow_orange",
        "Acropogon mesophilus": "blue_purple|yellow_orange",
        "Acropogon veillonii": "red_pink|yellow_orange",
    }
    assert colour["evidence_quality"].eq("high").all()
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


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


def test_third_rule_unlocks_are_species_direct_and_trait_specific() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["accepted_species"].isin(
            {"Jacobaea aquatica", "Carpinus laxiflora", "Citharexylum myrianthum"}
        )
    ].set_index("accepted_species")

    assert rows[["trait_name", "normalized_value", "evidence_quality"]].to_dict(
        "index"
    ) == {
        "Jacobaea aquatica": {
            "trait_name": "floral_form",
            "normalized_value": "composite_head",
            "evidence_quality": "medium",
        },
        "Carpinus laxiflora": {
            "trait_name": "inflorescence_display",
            "normalized_value": "raceme_spike_panicle",
            "evidence_quality": "high",
        },
        "Citharexylum myrianthum": {
            "trait_name": "inflorescence_display",
            "normalized_value": "raceme_spike_panicle",
            "evidence_quality": "high",
        },
    }
    assert rows["evidence_scope"].eq("species_direct").all()
    assert not set(rows["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_fourth_rule_unlocks_keep_autogamy_and_symmetry_trait_specific() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["accepted_species"].isin({"Hakea carinata", "Celtis sinensis"})
    ].set_index("accepted_species")

    assert rows[["trait_name", "normalized_value", "evidence_quality"]].to_dict(
        "index"
    ) == {
        "Hakea carinata": {
            "trait_name": "autonomous_selfing_capacity",
            "normalized_value": "autonomous",
            "evidence_quality": "high",
        },
        "Celtis sinensis": {
            "trait_name": "floral_symmetry",
            "normalized_value": "actinomorphic",
            "evidence_quality": "medium",
        },
    }
    assert rows["evidence_scope"].eq("species_direct").all()
    assert rows["retrieved_at_utc"].eq("2026-08-12T09:53:30Z").all()
    assert not set(rows["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )


def test_fifth_rule_unlocks_keep_horticulture_medium_and_government_db_high() -> None:
    evidence = pd.DataFrame(reviewed_rows()).fillna("")
    rows = evidence.loc[
        evidence["accepted_species"].isin(
            {"Phoenix roebelenii", "Gonystylus confusus"}
        )
    ].set_index("accepted_species")

    assert rows[["trait_name", "normalized_value", "evidence_quality"]].to_dict(
        "index"
    ) == {
        "Phoenix roebelenii": {
            "trait_name": "inflorescence_display",
            "normalized_value": "raceme_spike_panicle",
            "evidence_quality": "medium",
        },
        "Gonystylus confusus": {
            "trait_name": "flower_primary_color",
            "normalized_value": "yellow_orange",
            "evidence_quality": "high",
        },
    }
    assert rows.loc[
        "Phoenix roebelenii", "wild_cultivated_cultivar_status"
    ] == "species_level_horticultural_record_not_cultivar_limited"
    assert rows["evidence_scope"].eq("species_direct").all()
    assert not set(rows["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )
