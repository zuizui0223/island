"""Synthetic original taxon-disjoint NHM joint colour×form source support checks."""
import pandas as pd
import pytest
from fd_nhm_source_only_joint_colour_form_opportunity import (
    cross_corpus_joint,source_provenance_opportunity,
)


def _provider():
    return {
        "BHL_text":dict(
            primary={"A a":"white","B b":"yellow_orange","C c":"white","D d":"red_pink"},
            form={"A a":"tubular","C c":"open_radial","D d":"tubular"}),
        "FloraChina_text":dict(
            primary={"A a":"white","B b":"yellow_orange","C c":"blue_purple"},
            form={"A a":"tubular","B b":"open_radial"}),
    }


def test_exact_same_source_pair_restricted_vs_complementary_corpora():
    hierarchical,qa=cross_corpus_joint(_provider())
    same=hierarchical["same_published_corpus_both_traits"]
    joint=hierarchical["any_separate_corpora_consistent_both_traits"]
    twice=hierarchical["two_distinct_named_corpora_both_traits_agree"]
    assert same=={"A a":("white","tubular"),"D d":("red_pink","tubular")}
    assert joint=={
        "A a":("white","tubular"),
        "B b":("yellow_orange","open_radial"),
        "D d":("red_pink","tubular")}
    assert twice=={"A a":("white","tubular")}
    assert qa["n_conflicted_primary_colour_species"]==1
    assert qa["n_dual_corpus_form_species"]==1
    with pytest.raises(ValueError):
        cross_corpus_joint({})


def _flora():
    return {
        "I1":{"F sp1","F sp2","F sp3","G sp4","G sp5"},
        "I2":{"A sp1","A sp2","B sp3","B sp4","C sp5"},
    }


def test_true_global_external_name_exclusion_and_source_congeners():
    all_flora=_flora()
    donor={"F source1","F source2","B source1","I1 sp1"}
    gift=pd.DataFrame([
        {"entity_ID":"S1","work_species":"F source1"},
        {"entity_ID":"S1","work_species":"F source2"},
        {"entity_ID":"S2","work_species":"B source1"},
        {"entity_ID":"S2","work_species":"I1 sp1"},
    ])
    source=pd.DataFrame([
        dict(island_id=i,source_mode="geo_k5",source_rank=r,entity_ID=s)
        for i,s in (("I1","S1"),("I2","S2")) for r in (1,)
    ])
    # The donor species named I1 sp1 is deliberately excluded as an island
    # botanical member even if it were misreported as GIFT native source.
    all_flora["I1"].add("I1 sp1")
    target={sp for v in all_flora.values() for sp in v}
    # Same original accepted species label may be present on a GIFT list
    # but cannot count as a truly outside-the-island donor.
    q=source_provenance_opportunity(
        donor,all_flora,["I1","I2"],target,gift,source)
    assert q["n_source_only_names_with_both_traits"]==3
    assert q["n_target_islands_with_at_least_one_external_joint_source_donor"]==2
    assert q["n_target_islands_with_at_least_five_external_joint_source_donors"]==0
    assert q["n_direct_colour_form_target_plant_carrier_occurrences"]==11
    assert q["n_target_plant_carriers_with_nearest_external_joint_congener"]==5
    assert q["n_target_islands_with_any_source_congeneric_external_joint_donor"]==2


def test_source_shape_not_imputed_as_monotypic_or_exhaustive():
    all_flora=_flora()
    source=pd.DataFrame([dict(entity_ID="S",work_species="Y novel")])
    rank=pd.DataFrame([
        dict(island_id=i,source_mode="geo_k5",source_rank=1,entity_ID="S")
        for i in ("I1","I2")])
    q=source_provenance_opportunity(
        set(),all_flora,["I1","I2"],
        set().union(*all_flora.values()),source,rank)
    assert q["n_source_only_names_with_both_traits"]==0
    assert q["n_target_islands_with_at_least_one_external_joint_source_donor"]==0
    assert q["fraction_direct_joint_plant_carriers_with_source_genus_support"]==0
    with pytest.raises(ValueError):
        source_provenance_opportunity(
            set(),all_flora,["I1"],set(),
            source,rank)
