import pandas as pd

from island_v2.targeted_acquisition_priority import build_priority


def test_prioritizes_near_threshold_genera_by_leverage():
    evidence = pd.DataFrame([
        {"accepted_species": "Alpha one", "trait_name": "floral_form", "evidence_quality": "high", "evidence_scope": "species_direct"},
        {"accepted_species": "Alpha two", "trait_name": "floral_form", "evidence_quality": "medium", "evidence_scope": "species_direct"},
        {"accepted_species": "Beta one", "trait_name": "floral_form", "evidence_quality": "high", "evidence_scope": "species_direct"},
    ])
    unresolved = pd.DataFrame([
        {"accepted_species": "Alpha three", "axis": "floral_structural_complexity"},
        {"accepted_species": "Alpha four", "axis": "floral_structural_complexity"},
        {"accepted_species": "Alpha five", "axis": "floral_structural_complexity"},
        {"accepted_species": "Beta two", "axis": "floral_structural_complexity"},
        {"accepted_species": "Beta three", "axis": "floral_structural_complexity"},
    ])

    priority, targets, summary = build_priority(evidence, unresolved)

    assert list(priority["genus"]) == ["Alpha", "Beta"]
    assert priority.loc[0, "additional_species_needed"] == 1
    assert priority.loc[0, "leverage_score"] == 3
    assert set(targets.loc[targets["genus"].eq("Alpha"), "accepted_species"]) == {"Alpha five"}
    assert len(targets.loc[targets["genus"].eq("Beta")]) == 2
    assert summary["target_species_axis_tasks"] == 3
    assert summary["potential_unresolved_species_axis_coverage"] == 5


def test_excludes_no_evidence_and_already_eligible_genera():
    evidence = pd.DataFrame([
        {"accepted_species": f"Gamma {name}", "trait_name": "flower_primary_color", "evidence_quality": "high", "evidence_scope": "species_direct"}
        for name in ("one", "two", "three")
    ])
    unresolved = pd.DataFrame([
        {"accepted_species": "Gamma four", "axis": "flower_colour"},
        {"accepted_species": "Delta one", "axis": "flower_colour"},
    ])

    priority, targets, summary = build_priority(evidence, unresolved)

    assert priority.empty
    assert targets.empty
    assert summary["near_threshold_genus_axis_groups"] == 0
