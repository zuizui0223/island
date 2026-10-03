import pandas as pd

from island_v2.chapter1_h1_three_axis_raw_state import (
    build_axis_state_counts,
    build_valid_state_ledger,
)


def _config():
    return {
        "axes": {
            "flower_colour": {"traits": ["flower_primary_color"]},
            "floral_structural_complexity": {
                "traits": ["floral_symmetry", "floral_form"]
            },
            "reproductive_assurance": {
                "traits": ["self_incompatibility", "mating_system"]
            },
        },
        "evidence_scopes": {
            "all_analysis_eligible": ["high", "medium", "low"],
            "direct_only": ["high", "medium"],
        },
    }


def _ontology():
    return {
        "traits": {
            "flower_primary_color": {
                "allowed_values": ["white", "red_pink", "unresolved"]
            },
            "floral_symmetry": {
                "allowed_values": ["actinomorphic", "zygomorphic", "unresolved"]
            },
            "floral_form": {
                "allowed_values": ["open_radial", "tubular", "unresolved"]
            },
            "self_incompatibility": {
                "allowed_values": ["SC", "SI", "unresolved"]
            },
            "mating_system": {
                "allowed_values": [
                    "predominantly_outcrossing",
                    "predominantly_selfing",
                    "unresolved",
                ]
            },
        }
    }


def test_single_component_axis_cell_is_retained():
    cells = pd.DataFrame(
        [
            {
                "accepted_species": "A one",
                "axis": "reproductive_assurance",
                "trait_composition": 'self_incompatibility=["SC"]',
                "quality": "high",
            },
            {
                "accepted_species": "B two",
                "axis": "reproductive_assurance",
                "trait_composition": 'mating_system=["predominantly_selfing"]',
                "quality": "medium",
            },
        ]
    )
    ledger, audit = build_valid_state_ledger(
        cells, _ontology(), _config(), evidence_scope="all_analysis_eligible"
    )
    assert set(ledger["accepted_species"]) == {"A one", "B two"}
    row = audit.loc[audit["axis"].eq("reproductive_assurance")].iloc[0]
    assert int(row["resolved_axis_cells"]) == 2
    assert int(row["ontology_valid_axis_cells"]) == 2


def test_multistate_is_retained_and_cross_field_state_is_filtered():
    cells = pd.DataFrame(
        [
            {
                "accepted_species": "A one",
                "axis": "flower_colour",
                "trait_composition": 'flower_primary_color=["white","red_pink"]',
                "quality": "high",
            },
            {
                "accepted_species": "B two",
                "axis": "floral_structural_complexity",
                "trait_composition": 'floral_symmetry=["actinomorphic","tubular"]',
                "quality": "high",
            },
        ]
    )
    ledger, audit = build_valid_state_ledger(
        cells, _ontology(), _config(), evidence_scope="all_analysis_eligible"
    )
    colour = ledger.loc[ledger["accepted_species"].eq("A one"), "state"].tolist()
    assert set(colour) == {"white", "red_pink"}
    structure = ledger.loc[ledger["accepted_species"].eq("B two"), "state"].tolist()
    assert structure == ["actinomorphic"]
    assert int(audit["ontology_invalid_only_axis_cells"].sum()) == 0


def test_axis_state_counts_use_trait_specific_denominators():
    cells = pd.DataFrame(
        [
            {
                "accepted_species": "A one",
                "axis": "reproductive_assurance",
                "trait_composition": 'self_incompatibility=["SC"]',
                "quality": "high",
            },
            {
                "accepted_species": "B two",
                "axis": "reproductive_assurance",
                "trait_composition": 'self_incompatibility=["SI"]|mating_system=["predominantly_selfing"]',
                "quality": "high",
            },
        ]
    )
    ledger, _ = build_valid_state_ledger(
        cells, _ontology(), _config(), evidence_scope="all_analysis_eligible"
    )
    flora = pd.DataFrame(
        {
            "island_id": ["I1", "I1"],
            "accepted_species": ["A one", "B two"],
        }
    )
    counts, _ = build_axis_state_counts(flora, ledger, stratum="all_observed")
    sc = counts.loc[counts["outcome"].eq("self_incompatibility::SC")].iloc[0]
    si = counts.loc[counts["outcome"].eq("self_incompatibility::SI")].iloc[0]
    mating = counts.loc[
        counts["outcome"].eq("mating_system::predominantly_selfing")
    ].iloc[0]
    assert int(sc["trials"]) == 2 and int(sc["successes"]) == 1
    assert int(si["trials"]) == 2 and int(si["successes"]) == 1
    assert int(mating["trials"]) == 1 and int(mating["successes"]) == 1
