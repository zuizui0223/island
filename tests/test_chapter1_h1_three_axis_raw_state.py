import pandas as pd

from island_v2.chapter1_h1_three_axis_raw_state import (
    build_axis_state_counts,
    build_valid_state_ledger,
    formal_state_support,
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


def test_unique_within_axis_trait_permutation_is_recovered():
    cells = pd.DataFrame(
        [
            {
                "accepted_species": "A one",
                "axis": "floral_structural_complexity",
                "trait_composition": (
                    'floral_symmetry=["raceme_spike_panicle"]|'
                    'inflorescence_display=["zygomorphic"]'
                ),
                "quality": "low",
            }
        ]
    )
    config = _config()
    config["axes"]["floral_structural_complexity"]["traits"].append(
        "inflorescence_display"
    )
    ontology = _ontology()
    ontology["traits"]["inflorescence_display"] = {
        "allowed_values": ["raceme_spike_panicle", "unresolved"]
    }
    ledger, audit = build_valid_state_ledger(
        cells,
        ontology,
        config,
        evidence_scope="all_analysis_eligible",
    )
    got = set(zip(ledger["trait_name"], ledger["state"], strict=False))
    assert got == {
        ("inflorescence_display", "raceme_spike_panicle"),
        ("floral_symmetry", "zygomorphic"),
    }
    row = audit.loc[
        audit["axis"].eq("floral_structural_complexity")
    ].iloc[0]
    assert int(row["ontology_repaired_axis_cells"]) == 1
    assert int(row["n_reassigned_state_memberships"]) == 2
    assert int(row["ontology_invalid_only_axis_cells"]) == 0


def test_ambiguous_invalid_state_is_not_silently_reassigned():
    cells = pd.DataFrame(
        [
            {
                "accepted_species": "A one",
                "axis": "reproductive_assurance",
                "trait_composition": (
                    'self_incompatibility=["absent"]|'
                    'autonomous_selfing_capacity=["mixed_or_variable"]'
                ),
                "quality": "high",
            }
        ]
    )
    config = _config()
    config["axes"]["reproductive_assurance"]["traits"].append(
        "autonomous_selfing_capacity"
    )
    ontology = _ontology()
    ontology["traits"]["self_incompatibility"]["allowed_values"].append(
        "mixed_or_variable"
    )
    ontology["traits"]["autonomous_selfing_capacity"] = {
        "allowed_values": ["absent", "mixed_or_variable", "unresolved"]
    }
    ledger, audit = build_valid_state_ledger(
        cells,
        ontology,
        config,
        evidence_scope="all_analysis_eligible",
    )
    # The invalid cross-field state is dropped; the already-valid current-trait
    # state is retained. No cross-field guess is made.
    assert set(zip(ledger["trait_name"], ledger["state"], strict=False)) == {
        ("autonomous_selfing_capacity", "mixed_or_variable")
    }
    row = audit.loc[
        audit["axis"].eq("reproductive_assurance")
    ].iloc[0]
    assert int(row["ontology_repaired_axis_cells"]) == 0
    assert int(row["ontology_invalid_only_axis_cells"]) == 0


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


def test_formal_state_support_excludes_single_species_state_without_dropping_other_states():
    ledger = pd.DataFrame(
        [
            {"accepted_species": "A one", "axis": "reproductive_assurance", "trait_name": "cleistogamy", "state": "obligate"},
            {"accepted_species": "A one", "axis": "reproductive_assurance", "trait_name": "self_incompatibility", "state": "SC"},
            {"accepted_species": "B two", "axis": "reproductive_assurance", "trait_name": "self_incompatibility", "state": "SC"},
            {"accepted_species": "C three", "axis": "reproductive_assurance", "trait_name": "self_incompatibility", "state": "SI"},
        ]
    )
    filtered, support = formal_state_support(
        ledger,
        minimum_unique_species=2,
        evidence_scope="all_analysis_eligible",
    )
    assert "cleistogamy" not in set(filtered["trait_name"])
    assert set(filtered["trait_name"]) == {"self_incompatibility"}
    obligate = support.loc[
        support["trait_name"].eq("cleistogamy")
        & support["state"].eq("obligate")
    ].iloc[0]
    assert int(obligate["n_unique_species"]) == 1
    assert bool(obligate["eligible_for_formal_axis_test"]) is False
