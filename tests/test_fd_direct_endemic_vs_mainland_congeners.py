"""Synthetic actual Direct endemic species vs different source congeners.

These functions compare observed category codes; do not infer ancestors,
adaptive island populations, effective guilds or unobserved arrival.
"""
import pytest
from fd_direct_endemic_vs_mainland_congeners import (
    direct_plant_vs_actual_source_congeners
)


def test_exact_source_congener_signed_direction_and_absolute_disagreement():
    state={
        "A endemic":1,"A mainland1":0,"A mainland2":1,
        "B endemic":0,"B mainland1":1,
        "C endemic":1,"D mainland1":1}
    focal={"A endemic","B endemic","C endemic"}
    source={"A mainland1","A mainland2","B mainland1","D mainland1"}
    r=direct_plant_vs_actual_source_congeners(focal,source,state)
    assert len(r)==2
    by={x["_private_accepted_plant_name"]:x for x in r}
    assert by["A endemic"]["signed_plant_minus_actual_congener_category"]==pytest.approx(.5)
    assert by["A endemic"]["absolute_plant_congener_discordance"]==pytest.approx(.5)
    assert by["A endemic"]["real_source_congener_n"]==2
    assert by["A endemic"]["real_source_congener_has_variable_categories"] is True
    assert by["B endemic"]["signed_plant_minus_actual_congener_category"]==-1
    assert by["B endemic"]["absolute_plant_congener_discordance"]==1
    assert by["B endemic"]["real_source_congener_has_variable_categories"] is False
    assert "C endemic" not in by


def test_no_source_congener_does_not_become_zero_ancestral_trait():
    states={"A endemic":1,"B mainland":0}
    records=direct_plant_vs_actual_source_congeners(
        {"A endemic"},{"B mainland"},states)
    assert records==[]


def test_focal_accepted_name_and_missing_source_plant_label_leakage_fails_closed():
    states={"A endemic":1,"A mainland":0}
    with pytest.raises(ValueError):
        direct_plant_vs_actual_source_congeners(
            {"A endemic"},{"A endemic"},states)
    with pytest.raises(ValueError):
        direct_plant_vs_actual_source_congeners(
            {"A endemic"},{"A unknown"},states)
    with pytest.raises(ValueError):
        direct_plant_vs_actual_source_congeners(
            {"A unknown"},{"A mainland"},states)


def test_differing_native_endemic_source_states_do_not_require_other_traits():
    states={"A endemic1":1,"A mainland1":0,"A native2":0}
    a=direct_plant_vs_actual_source_congeners(
        {"A endemic1"},{"A mainland1"},states)
    b=direct_plant_vs_actual_source_congeners(
        {"A native2"},{"A mainland1"},states)
    assert len(a)==1 and len(b)==1
    assert a[0]["signed_plant_minus_actual_congener_category"]==1
    assert b[0]["signed_plant_minus_actual_congener_category"]==0
