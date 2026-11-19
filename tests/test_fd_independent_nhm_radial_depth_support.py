"""Synthetic original-source NHM symmetry/form plus Direct tube-depth gate tests.

No real GIFT plant names, source IDs, island IDs, providers or coordinates.
"""
import pytest

from run_fd_independent_nhm_radial_depth_support import (
    external_hybrid_donors,external_donors_for_focal,
)


def test_independent_A_and_original_Direct_depth_B_same_species_only():
    NHM={
        "Alpha mainland1":"actinomorphic",
        "Alpha mainland2":"zygomorphic",
        "Beta mainland1":"actinomorphic",
        "Beta mainland2":"actinomorphic",
        "Beta island1":"actinomorphic",
        "Gamma mainland1":"actinomorphic",
    }
    depth={
        "Alpha mainland1":"deep",
        "Alpha mainland2":"deep",
        "Beta mainland1":"shallow",
        "Beta mainland2":"deep",
        "Beta island1":"deep",
        # Gamma lacks independent original true tube depth.
    }
    gift=set(NHM)
    island={"Beta island1"}
    scored=external_hybrid_donors(
        NHM,depth,gift,island,frozenset({"actinomorphic"})
    )
    assert scored=={
        "Alpha mainland1":1,
        "Alpha mainland2":0,
        "Beta mainland1":0,
        "Beta mainland2":1,
    }
    assert not scored.keys() & island
    assert "Gamma mainland1" not in scored


def test_unknown_direct_depth_is_excluded_not_turned_into_false():
    out=external_hybrid_donors(
        {"Alpha mainland1":"actinomorphic"},
        {},
        {"Alpha mainland1"},set(),frozenset({"actinomorphic"})
    )
    assert out=={}


def test_more_than_one_external_positive_shape_state_allowed():
    rows=external_hybrid_donors(
        {"Beta mainland1":"bell_campanulate",
         "Beta mainland2":"tubular","Beta mainland3":"open_radial"},
        {name:"deep" for name in
         ("Beta mainland1","Beta mainland2","Beta mainland3")},
        {"Beta mainland1","Beta mainland2","Beta mainland3"},
        set(),frozenset({"tubular","bell_campanulate"})
    )
    assert list(sorted(rows.values()))==[0,1,1]


def test_requires_explicit_positive_external_trait_category():
    with pytest.raises(ValueError,match="explicit"):
        external_hybrid_donors(
            {"Alpha x":"actinomorphic"},
            {"Alpha x":"deep"},
            {"Alpha x"},set(),{"actinomorphic"}
        )


def test_all_original_island_species_removed_regardless_of_status():
    a={"Alpha introduced":"actinomorphic",
       "Alpha native":"actinomorphic",
       "Alpha unknown":"actinomorphic"}
    out=external_hybrid_donors(
        a,{name:"deep" for name in a},set(a),
        set(a),frozenset({"actinomorphic"})
    )
    assert out=={}


def test_focal_only_name_exclusion_keeps_external_species_from_other_islands():
    # A plant on island #1 can be a valid source reference for island #2
    # if and only if it is absent from the focal #2 accepted flora.
    eligible={"Alpha island1","Alpha mainland1","Beta mainland1"}
    selected_source={"Alpha island1","Alpha mainland1"}
    focal={"Beta island2"}
    q=external_donors_for_focal(selected_source,eligible,focal)
    assert q=={"Alpha island1","Alpha mainland1"}
    assert external_donors_for_focal(
        selected_source,eligible,{"Alpha island1"}
    )=={"Alpha mainland1"}
    assert not external_donors_for_focal(
        selected_source,eligible,{"Alpha island1"}
    )&{"Alpha island1"}
