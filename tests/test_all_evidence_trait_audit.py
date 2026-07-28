from __future__ import annotations

import json

import pandas as pd

from island_v2 import all_evidence_trait_audit as audit
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

ONTOLOGY = {
    "flower_primary_color": {
        "white",
        "red_pink",
        "yellow_orange",
        "blue_purple",
        "green_brown_inconspicuous",
        "multicolored_variable",
        "other_described",
    },
    "floral_form": {"tubular", "open_radial"},
    "floral_symmetry": {"actinomorphic", "zygomorphic"},
    "self_incompatibility": {"SI", "SC", "mixed_or_variable"},
    "autonomous_selfing_capacity": {"absent", "autonomous"},
    "mating_system": {
        "predominantly_outcrossing",
        "mixed_mating",
        "predominantly_selfing",
        "obligate_selfing",
    },
}


def row(
    species: str,
    trait: str,
    value: str,
    lineage: str,
    *,
    quality: str = "medium",
    source_group: str = "example",
    excerpt: str = "",
) -> dict[str, str]:
    axis = audit.trait_axis(trait)
    return {
        "accepted_species": species,
        "axis": axis,
        "trait_name": trait,
        "normalized_value": value,
        "quality": quality,
        "source_group": source_group,
        "source_provider": source_group,
        "source_url": f"https://example.org/{lineage}",
        "source_record_id": lineage,
        "source_citation": "Example source",
        "source_excerpt": excerpt or f"{trait}={value}",
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact",
        "source_lineage": lineage,
        "lineage_method": "test",
        "source_run_id": "1",
        "source_artifact": "test",
        "source_file": "test.csv",
        "acceptance_contract": "test",
    }


def test_lineage_dedup_is_trait_specific_not_axis_specific() -> None:
    evidence = pd.DataFrame(
        [
            row("Alpha one", "floral_form", "tubular", "page:1"),
            row("Alpha one", "floral_symmetry", "actinomorphic", "page:1"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, duplicates = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    assert len(lineages) == 2
    assert duplicates.empty
    assert set(lineages["trait_name"]) == {"floral_form", "floral_symmetry"}


def test_direct_conflicts_are_classified_without_row_order_selection() -> None:
    evidence = pd.DataFrame(
        [
            row("Colour one", "flower_primary_color", "red_pink", "colour:a"),
            row("Colour one", "flower_primary_color", "white", "colour:b"),
            row("Repro one", "self_incompatibility", "SC", "repro:a"),
            row("Repro one", "self_incompatibility", "SI", "repro:b"),
            row("Garden one", "floral_form", "tubular", "garden:wild"),
            row(
                "Garden one",
                "floral_form",
                "open_radial",
                "garden:cultivar",
                excerpt="cultivar selection has open radial flowers",
            ),
            row("Ontology one", "floral_form", "not_an_ontology_value", "bad:a"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    resolved, cell_audit = audit.resolve_direct_cells(lineages)

    colour = resolved.loc[
        resolved["accepted_species"].eq("Colour one")
    ].iloc[0]
    assert colour["classification"] == "true_multistate_variable"
    assert json.loads(colour["state_set"]) == ["red_pink", "white"]

    repro = cell_audit.loc[
        cell_audit["accepted_species"].eq("Repro one")
    ].iloc[0]
    assert repro["classification"] == "unresolved_direct_conflict"
    assert repro["resolution_status"] == "excluded"

    garden = resolved.loc[
        resolved["accepted_species"].eq("Garden one")
    ].iloc[0]
    assert garden["classification"] == "cultivar_contamination"
    assert garden["normalized_value"] == "tubular"

    ontology = cell_audit.loc[
        cell_audit["accepted_species"].eq("Ontology one")
    ].iloc[0]
    assert ontology["classification"] == "source_ontology_mismatch"


def test_legacy_colour_aliases_remain_multistate_and_internal_conflicts_are_explicit() -> None:
    states, invalid = audit.normalise_state_set(
        "flower_primary_color",
        "cream|violet|black|multicolored/variable",
        ONTOLOGY,
    )
    assert states == (
        "blue_purple",
        "multicolored_variable",
        "other_described",
        "white",
    )
    assert invalid == ()

    evidence = pd.DataFrame(
        [
            row("Delta one", "self_incompatibility", "SC", "same:lineage"),
            row("Delta one", "self_incompatibility", "SI", "same:lineage"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    resolved, cell_audit = audit.resolve_direct_cells(lineages)
    assert resolved.empty
    assert cell_audit.iloc[0]["classification"] == "unresolved_direct_conflict"
    assert cell_audit.iloc[0]["lineage_internal_conflicts"] == 1


def test_thresholds_and_application_remain_trait_specific() -> None:
    evidence_rows = [
        row(f"Beta {index}", "floral_form", value, f"beta:{index}")
        for index, value in enumerate(
            ["tubular", "tubular", "tubular", "open_radial"],
            start=1,
        )
    ]
    evidence_rows += [
        row(f"Gamma {index}", "floral_symmetry", "actinomorphic", f"gamma:{index}")
        for index in range(1, 3)
    ]
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    cells, _ = audit.resolve_direct_cells(lineages)
    cells["genus"] = cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    old_low = pd.DataFrame(
        columns=["accepted_species", "genus", "axis", "trait_name", "state_set"]
    )
    rules = audit.build_rule_audit(cells, lineages, old_low)

    beta_current = rules.loc[
        rules["setting"].eq("current_min3")
        & rules["genus"].eq("Beta")
        & rules["trait_name"].eq("floral_form")
    ].iloc[0]
    beta_relaxed = rules.loc[
        rules["setting"].eq("relaxed_min3")
        & rules["genus"].eq("Beta")
        & rules["trait_name"].eq("floral_form")
    ].iloc[0]
    assert not bool(beta_current["eligible"])
    assert bool(beta_relaxed["eligible"])

    gamma_min3 = rules.loc[
        rules["setting"].eq("current_min3")
        & rules["genus"].eq("Gamma")
    ].iloc[0]
    gamma_min2 = rules.loc[
        rules["setting"].eq("current_min2_diagnostic")
        & rules["genus"].eq("Gamma")
    ].iloc[0]
    assert not bool(gamma_min3["eligible"])
    assert bool(gamma_min2["eligible"])
    assert bool(gamma_min2["diagnostic_only"])

    master = pd.DataFrame(
        [
            {
                "accepted_species": "Gamma target",
                "genus": "Gamma",
                "n_islands": "3",
                "n_records": "10",
            }
        ]
    )
    inferred = audit.apply_genus_rules(
        master,
        cells,
        rules,
        "current_min2_diagnostic",
    )
    assert len(inferred) == 1
    assert inferred.iloc[0]["trait_name"] == "floral_symmetry"


def test_species_axis_coverage_uses_quality_precedence_and_exact_denominator() -> None:
    master = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "genus": "Alpha",
                "n_islands": "1",
                "n_records": "1",
            },
            {
                "accepted_species": "Beta two",
                "genus": "Beta",
                "n_islands": "1",
                "n_records": "1",
            },
        ]
    )
    direct = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "axis": "flower_colour",
                "trait_name": "flower_primary_color",
                "state_set": '["white"]',
                "normalized_value": "white",
                "quality": "high",
                "source_groups": "direct",
                "source_lineages": "lineage:a",
                "classification": "single_independent_lineage",
            }
        ]
    )
    low = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "axis": "flower_colour",
                "trait_name": "flower_primary_color",
                "state_set": '["red_pink"]',
                "normalized_value": "red_pink",
                "source_lineage": "rule:a",
            },
            {
                "accepted_species": "Beta two",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "state_set": '["SC"]',
                "normalized_value": "SC",
                "source_lineage": "rule:b",
            },
        ]
    )
    coverage = audit.species_axis_coverage(master, direct, low)
    snapshot = audit.coverage_snapshot(coverage)
    assert len(coverage) == 6
    assert snapshot["quality_counts"] == {"high": 1, "medium": 0, "low": 1}
    assert snapshot["unresolved_species_axis"] == 4


def test_acquisition_queue_ranks_cells_unlocked_per_required_acquisition() -> None:
    master = pd.DataFrame(
        [
            {
                "accepted_species": f"Alpha {index}",
                "genus": "Alpha",
                "n_islands": "1",
                "n_records": "1",
            }
            for index in range(100)
        ]
        + [
            {
                "accepted_species": f"Beta {index}",
                "genus": "Beta",
                "n_islands": "1",
                "n_records": "1",
            }
            for index in range(80)
        ]
    )
    coverage = master[["accepted_species"]].copy()
    coverage["axis"] = "reproductive_assurance"
    coverage["quality"] = ""
    rules = pd.DataFrame(
        [
            {
                "setting": "current_min3",
                "genus": "Alpha",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "n_direct_species": 1,
                "dominant_species": 1,
                "required_dominance": 0.95,
                "dominance": 1.0,
                "value_distribution": '{"[\\"SC\\"]":1}',
            },
            {
                "setting": "current_min3",
                "genus": "Beta",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "n_direct_species": 2,
                "dominant_species": 2,
                "required_dominance": 0.95,
                "dominance": 1.0,
                "value_distribution": '{"[\\"SC\\"]":2}',
            },
        ]
    )
    queue = audit.acquisition_queue(master, coverage, rules)
    assert queue.iloc[0]["genus"] == "Beta"
    assert queue.iloc[0]["acquisitions_needed_for_min3_dominance"] == 1
    assert queue.iloc[0]["potential_cells_per_acquisition"] == 80
