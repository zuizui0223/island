import pandas as pd

from island_v2.chapter1_context_input import build_island_trait_composition


def _status(flora):
    return pd.DataFrame(
        {
            "island_id": flora["island_id"],
            "accepted_species": flora["species"],
            "origin_status": ["native"] * len(flora),
            "endemic_status": ["nonendemic"] * len(flora),
        }
    )


def test_colour_and_form_multistates_expand_to_atomic_presence():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "species": ["A a", "A b"],
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "A a", "A b"],
            "trait_name": [
                "flower_primary_color",
                "flower_primary_color",
                "floral_form",
                "floral_form",
            ],
            "normalized_value": [
                "red_pink|white",
                "white",
                "tubular|salverform",
                "open_radial",
            ],
            "evidence_scope": ["species_direct"] * 4,
            "resolution_status": ["resolved"] * 4,
        }
    )
    composition, _, audit = build_island_trait_composition(flora, _status(flora), evidence)

    colour = composition.loc[composition["trait_name"].eq("flower_primary_color")]
    assert set(colour["trait_state"]) == {"red_pink", "white"}
    assert int(colour.loc[colour["trait_state"].eq("red_pink"), "successes"].iloc[0]) == 1
    assert int(colour.loc[colour["trait_state"].eq("white"), "successes"].iloc[0]) == 2
    assert set(colour["trials"]) == {2}

    form = composition.loc[composition["trait_name"].eq("floral_form")]
    assert set(form["trait_state"]) == {"tubular", "salverform", "open_radial"}
    assert set(form["trials"]) == {2}

    multistate = audit.loc[
        audit["accepted_species"].eq("A a")
        & audit["trait_name"].eq("flower_primary_color")
    ].iloc[0]
    assert multistate["multistate_policy"] == "atomic_presence"
    assert multistate["n_output_states"] == 2


def test_sc_si_multistate_is_one_mixed_state_not_double_counted():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1"],
            "species": ["A a", "A b", "A c"],
        }
    )
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "A c"],
            "trait_name": ["self_incompatibility"] * 3,
            "normalized_value": ["SC|SI", "SC", "SI"],
            "evidence_scope": ["species_direct"] * 3,
            "resolution_status": ["resolved"] * 3,
        }
    )
    composition, _, _ = build_island_trait_composition(flora, _status(flora), evidence)
    si = composition.loc[composition["trait_name"].eq("self_incompatibility")]
    assert set(si["trait_state"]) == {"SC", "SI", "mixed_or_variable"}
    assert set(si["successes"]) == {1}
    assert set(si["trials"]) == {3}


def test_conflicting_sc_and_si_sources_still_fail_closed():
    flora = pd.DataFrame({"island_id": ["i1"], "species": ["A a"]})
    evidence = pd.DataFrame(
        {
            "accepted_species": ["A a", "A a"],
            "trait_name": ["self_incompatibility", "self_incompatibility"],
            "normalized_value": ["SC", "SI"],
            "evidence_scope": ["species_direct", "species_direct"],
            "resolution_status": ["resolved", "resolved"],
        }
    )
    composition, _, audit = build_island_trait_composition(flora, _status(flora), evidence)
    assert composition.empty
    row = audit.iloc[0]
    assert not bool(row["resolved_for_primary"])
    assert row["n_distinct_signatures"] == 2
