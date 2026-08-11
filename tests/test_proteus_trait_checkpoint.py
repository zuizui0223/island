from __future__ import annotations

import pandas as pd

from island_v2.open_web_common import reviewed_source_package_evidence
from island_v2.proteus_trait_checkpoint import (
    CHARACTER_TRAIT,
    build_from_frames,
    deterministic_review_sample,
    exclude_resolved_direct_pairs,
    normalise_proteus_state,
)


def test_trait_specific_mapping_rejects_cross_trait_states() -> None:
    assert normalise_proteus_state(
        "82. Autonomous selfing (D1)",
        "82. Autonomous selfing (D1): (0) yes",
    ) == ("autonomous_selfing_capacity", "autonomous")
    assert normalise_proteus_state(
        "83. Self-incompatibility system (genetic) (D1)",
        "83. Self-incompatibility system (genetic) (D1): (0) self-compatible",
    ) == ("self_incompatibility", "SC")
    assert normalise_proteus_state(
        "80. Phenotypic mating system (D1)",
        "80. Phenotypic mating system (D1): (3) apomixis",
    ) == ("mating_system", "")
    assert normalise_proteus_state(
        "207. Symmetry of perianth (D1)",
        "207. Symmetry of perianth (D1): (3) disymmetric",
    ) == ("floral_symmetry", "")


def _row(index: int, species: str, character: str, state: str) -> dict[str, str]:
    return {
        "NDat": str(index),
        "Taxon": species,
        "Family": "Testaceae",
        "Character": character,
        "Reference": "Author 2025 (article)",
        "CharacterState": state,
        "Value": "",
        "Min": "",
        "Max": "",
        "Notes": "",
    }


def _references() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Citation": "Author 2025 (article)",
                "Type": "article",
                "Reference": "Author A. 2025. Test traits. doi:10.1234/example.1",
                "URL": "https://doi.org/10.1234/example.1",
            }
        ]
    )


def test_exact_name_family_gate_and_original_reference_lineage() -> None:
    data = pd.DataFrame(
        [
            _row(
                1,
                "Testus alpha",
                "82. Autonomous selfing (D1)",
                "82. Autonomous selfing (D1): (2) no",
            ),
            _row(
                2,
                "Testus beta",
                "83. Self-incompatibility system (genetic) (D1)",
                "83. Self-incompatibility system (genetic) (D1): (1) gametophytic self-incompatible",
            ),
            _row(
                3,
                "Testus absent",
                "82. Autonomous selfing (D1)",
                "82. Autonomous selfing (D1): (0) yes",
            ),
        ]
    )
    master = pd.DataFrame(
        [
            {"accepted_species": "Testus alpha", "family": "Testaceae"},
            {"accepted_species": "Testus beta", "family": "Testaceae"},
        ]
    )
    evidence, exclusions = build_from_frames(data, _references(), master)
    assert len(evidence) == 2
    assert set(evidence["normalized_value"]) == {"absent", "SI"}
    assert evidence["quality"].eq("high").all()
    assert evidence["source_lineage"].nunique() == 1
    assert evidence["source_lineage"].iloc[0] == "doi:10.1234/example.1"
    assert evidence["source_excerpt"].str.contains("Original reference key=").all()
    assert exclusions["reason"].tolist() == ["not_exact_accepted_species"]


def test_review_sample_and_production_gate_are_trait_specific() -> None:
    state_by_character = {
        "207. Symmetry of perianth (D1)": "207. Symmetry of perianth (D1): (0) actinomorphic (strictly polysymmetric)",
        "239. Main color of perianth at anthesis (D1)": "239. Main color of perianth at anthesis (D1): (3) white",
        "80. Phenotypic mating system (D1)": "80. Phenotypic mating system (D1): (2) mixed",
        "82. Autonomous selfing (D1)": "82. Autonomous selfing (D1): (1) variable",
        "83. Self-incompatibility system (genetic) (D1)": "83. Self-incompatibility system (genetic) (D1): (0) self-compatible",
    }
    rows: list[dict[str, str]] = []
    master_rows: list[dict[str, str]] = []
    index = 0
    for character in CHARACTER_TRAIT:
        for within_trait in range(42):
            index += 1
            species = f"Testus s{index}"
            rows.append(_row(index, species, character, state_by_character[character]))
            master_rows.append({"accepted_species": species, "family": "Testaceae"})
    evidence, exclusions = build_from_frames(
        pd.DataFrame(rows), _references(), pd.DataFrame(master_rows)
    )
    assert exclusions.empty
    audit = deterministic_review_sample(evidence)
    assert len(audit) == 200
    assert audit["trait_name"].value_counts().eq(40).all()
    assert audit["candidate_id"].nunique() == 200

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)
    assert summary["package_gate_passed"] is True
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert len(selected) == len(evidence)
    assert scopes["production_approved"].all()


def test_only_resolved_direct_pairs_are_excluded() -> None:
    evidence = pd.DataFrame(
        [
            {"accepted_species": "Testus alpha", "trait_name": "mating_system"},
            {"accepted_species": "Testus beta", "trait_name": "mating_system"},
            {"accepted_species": "Testus gamma", "trait_name": "mating_system"},
        ]
    )
    direct = pd.DataFrame(
        [
            {
                "accepted_species": "Testus alpha",
                "trait_name": "mating_system",
                "resolution_status": "resolved",
            },
            {
                "accepted_species": "Testus beta",
                "trait_name": "mating_system",
                "resolution_status": "unresolved",
            },
        ]
    )
    kept, excluded, pairs = exclude_resolved_direct_pairs(evidence, direct)
    assert kept["accepted_species"].tolist() == ["Testus beta", "Testus gamma"]
    assert excluded["accepted_species"].tolist() == ["Testus alpha"]
    assert pairs == {("Testus alpha", "mating_system")}
