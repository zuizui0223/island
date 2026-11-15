"""NHM CC BY 4.0 independent botanical-description flower source strict parsing.

All tests use artificial taxa and publicly generic trait words; NO protected
GIFT original plant/island data or source-region identities are included.
"""
from collections import defaultdict
import numpy as np
import pandas as pd
import pytest

from fd_nhm_external_flora_preflight import (
    accepted_exact_binomial,single_description_color,
    single_description_form,single_description_symmetry,
    nhm_external_categories,_coverage,_external_gift_source_coverage,
    EXTERNAL_COLUMNS,
)


def _row(taxon,**values):
    obj={x:"" for x in EXTERNAL_COLUMNS}
    obj["taxon"]=taxon
    obj.update(values)
    return obj


def test_species_name_acceptance_fails_closed_at_binomial():
    assert accepted_exact_binomial("  Example   valid  ")=="Example valid"
    assert accepted_exact_binomial("Example valid L.") is None
    assert accepted_exact_binomial("Example var. valid") is None
    assert accepted_exact_binomial("Example ×hybrid") is None
    assert accepted_exact_binomial("example valid") is None
    assert accepted_exact_binomial("") is None


def test_independent_flower_colour_refuses_mixed_or_ambiguous():
    assert single_description_color("whitish")=="white"
    assert single_description_color("blue-purple")=="blue_purple"
    assert single_description_color("rose-pink")=="red_pink"
    assert single_description_color("deep purple")=="blue_purple"
    assert single_description_color("yellowish")=="yellow_orange"
    assert single_description_color("greenish")=="green_brown_inconspicuous"
    assert single_description_color("red, white") is None
    assert single_description_color("green-yellow") is None
    assert single_description_color("unrecorded") is None
    assert single_description_color("no flower description") is None
    assert single_description_color("9 white") is None
    assert single_description_color("") is None


def test_form_and_symmetry_not_inferred_from_colour_or_pollen():
    assert single_description_form("campanulate")=="bell_campanulate"
    assert single_description_form("tubular")=="tubular"
    assert single_description_form("rotate")=="open_radial"
    assert single_description_form("ligulate, tubular") is None
    assert single_description_symmetry("actinomorphic")=="actinomorphic"
    assert single_description_symmetry("zygomorphic")=="zygomorphic"
    assert single_description_symmetry("radial") is None


def test_taxon_multisource_and_multiple_organs_require_consensus():
    d=pd.DataFrame([
        _row("Aster test",**{"flower colour":"pink","source":"ecoflora",
                             "flower shape":"campanulate"}),
        _row("Aster test",**{"petal colour":"rose","source":"eflora"}),
        _row("Beta valid",**{"flower colour":"white","source":"bhl"}),
        _row("Beta valid",**{"flower colour":"red","source":"ecoflora"}),
        _row("Gamma valid",**{"corolla colour":"yellowish","source":"bhl"}),
        _row("Delta valid",**{"corolla colour":"green",
                              "petal colour":"yellow","source":"eflora"}),
        _row("Epsilon valid subsp.",**{"flower colour":"white","source":"bhl"}),
    ])
    state,meta=nhm_external_categories(d)
    assert state["primary"]=={"Aster test":"red_pink"}
    assert state["expanded"]=={
        "Aster test":"red_pink","Gamma valid":"yellow_orange"}
    assert state["form"]=={"Aster test":"bell_campanulate"}
    assert state["symmetry"]=={}
    assert meta["n_external_rows"]==7
    assert meta["n_exact_taxon_binomial_names"]==4
    assert meta["n_species_with_documentation_conflict_by_axis"]["primary"]==1
    assert meta["n_species_with_documentation_conflict_by_axis"]["expanded"]==1


def test_exact_original_trait_agreement_never_fabricated_for_undocumented():
    q=_coverage({"A a":"white","B b":"red_pink"},
                {"A a":"white","B b":"yellow_orange","C c":"blue_purple"})
    assert q["n_exact_accepted_name_comparable"]==2
    assert q["n_equal_categories"]==1
    assert q["n_different_categories"]==1
    assert q["fraction_exact_name_agreement"]==pytest.approx(.5)
    assert _coverage({"A a":"white"},{"B b":"white"})[
        "fraction_exact_name_agreement"] is None


def test_source_donor_coverage_protects_all_original_target_names():
    original={"Island1":{"A first","B second"},"Island2":{"B third"}}
    # The target original-Direct codebook is NEVER a training donor source.
    measured={"__target_original":{
        "A first":"white","B second":"yellow_orange","B third":"white"},
        "A donor":"white","C donor":"blue_purple"}
    source={"S1":{"A donor"},"S2":{"C donor"}}
    near={"Island1":("S1","S2"),"Island2":("S2",)}
    r=_external_gift_source_coverage(original,["Island1","Island2"],
                                      measured,source,near)
    assert r["n_direct_original_target_carrier_occurrences"]==3
    assert r["n_original_target_carriers_with_near_source_external_same_genus"]==1
    assert r["n_original_floras_with_any_external_colour_source_species"]==2
    assert r["n_islands_with_any_exact_congeneric_external_source_donor"]==1
