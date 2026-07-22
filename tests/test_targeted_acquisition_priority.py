import pandas as pd

from island_v2.targeted_acquisition_priority import build_priority


def row(species: str, trait: str, value: str) -> dict[str, str]:
    return {
        "accepted_species": species,
        "trait_name": trait,
        "reported_value": value,
        "evidence_quality": "high",
        "evidence_scope": "species_direct",
    }


def test_uses_target_aligned_candidate_and_supportive_value():
    evidence = pd.DataFrame([
        row("Alpha one", "floral_form", "tubular"),
        row("Alpha two", "floral_form", "tubular"),
    ])
    unresolved = pd.DataFrame([
        {"accepted_species": "Alpha three", "axis": "floral_structural_complexity"},
        {"accepted_species": "Alpha four", "axis": "floral_structural_complexity"},
    ])
    hints = pd.DataFrame([
        {
            "accepted_species": "Alpha three",
            "trait_name": "flower_primary_color",
            "provisional_candidate_value": "red",
        },
        {
            "accepted_species": "Alpha four",
            "trait_name": "floral_form",
            "provisional_candidate_value": "tubular",
        },
    ])

    priority, targets, summary = build_priority(
        evidence, unresolved, candidate_hints=hints
    )

    assert list(targets["accepted_species"]) == ["Alpha four"]
    assert targets.iloc[0]["target_traits"] == "floral_form"
    assert targets.iloc[0]["target_value"] == "tubular"
    assert summary["candidate_hints_used"] is True
    assert len(priority) == 1


def test_rejects_groups_without_target_trait_candidates():
    evidence = pd.DataFrame([
        row("Beta one", "flower_primary_color", "white"),
        row("Beta two", "flower_primary_color", "white"),
    ])
    unresolved = pd.DataFrame([
        {"accepted_species": "Beta three", "axis": "flower_colour"}
    ])
    hints = pd.DataFrame([
        {
            "accepted_species": "Beta three",
            "trait_name": "floral_form",
            "candidate_value": "tubular",
        }
    ])

    priority, targets, summary = build_priority(
        evidence, unresolved, candidate_hints=hints
    )

    assert priority.empty
    assert targets.empty
    assert summary["rejected_groups_without_target_aligned_candidates"] == 1


def test_rejects_optimistically_insufficient_dominance():
    evidence = pd.DataFrame([
        row("Gamma one", "self_incompatibility", "SI"),
        row("Gamma two", "self_incompatibility", "SC"),
    ])
    unresolved = pd.DataFrame([
        {"accepted_species": "Gamma three", "axis": "reproductive_assurance"}
    ])
    hints = pd.DataFrame([
        {
            "accepted_species": "Gamma three",
            "trait_name": "self_incompatibility",
            "candidate_value": "SI",
        }
    ])

    priority, targets, summary = build_priority(
        evidence, unresolved, candidate_hints=hints
    )

    assert priority.empty
    assert targets.empty
    assert summary["rejected_groups_unable_to_meet_dominance"] == 1
