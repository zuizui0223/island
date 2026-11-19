"""Original GIFT geo/source TRAINING exclusion, not merely bootstrap blocks."""
import numpy as np
import pandas as pd
import pytest

from fd_global_geo_source_training_holdout import (
    original_full_flora_group_blacklists,source_training_ban,
    make_donor_models,score_paired_identical_original_species,
    original_group_equal_stress,MODES,
)


def _synthetic():
    flora={
        "focal":{"Alpha test0","Alpha test1","Alpha test2","Alpha test3","Alpha test4"},
        "nearby":{"Alpha nearby","Alpha local1","Alpha source_only"},
        "remote":{"Alpha local2","Alpha local3","Beta local1"},
        "same_source":{"Alpha local4","Alpha remote_source","Beta external"},
        "external":{"Alpha out1","Alpha out2","Alpha out3","Beta out1"},
    }
    blocks=dict(focal="B1",nearby="B1",remote="B2",
                same_source="B4",external="B3")
    sources=dict(focal="A1",nearby="A2",remote="A3",
                 same_source="A1",external="A4")
    states={x:("blue" if "test" in x or "local" in x
               else "red") for names in flora.values() for x in names}
    states.update({"Alpha out1":"red","Alpha out2":"red",
                   "Alpha out3":"red"})
    indices={"blue":0,"red":1}
    return flora,blocks,sources,states,indices


def test_blacklist_removes_ENTIRE_original_neighbor_flora_from_donor_training():
    flora,blocks,sources,_,_=_synthetic()
    g=original_full_flora_group_blacklists(flora,blocks,sources)
    sameblock=source_training_ban(
        "focal",flora["focal"],"B1","A1",g,"block_out")
    samesource=source_training_ban(
        "focal",flora["focal"],"B1","A1",g,"first_source_out")
    both=source_training_ban(
        "focal",flora["focal"],"B1","A1",g,"block_and_first_source_out")
    assert "Alpha local1" in sameblock
    assert "Alpha local1" not in samesource
    assert "Alpha local4" in samesource
    assert "Alpha local4" not in sameblock
    assert both==sameblock|samesource
    assert flora["focal"]<=sameblock&samesource
    assert g["n_islands_per_block"]["B1"]==2
    assert g["n_islands_per_source"]["A1"]==2


def test_missing_physical_block_or_source_group_never_falls_back_to_baseline():
    flora,blocks,sources,_,_=_synthetic()
    g=original_full_flora_group_blacklists(flora,blocks,sources)
    assert source_training_ban("focal",flora["focal"],None,"A1",g,"block_out") is None
    assert source_training_ban("focal",flora["focal"],"B1",None,g,
                               "first_source_out") is None
    assert source_training_ban("focal",flora["focal"],"B1",None,g,
                               "block_and_first_source_out") is None
    with pytest.raises(ValueError):
        source_training_ban("focal",flora["focal"],"B1","A1",g,"not_mode")


def test_recalculate_external_blind_prior_AFTER_geographic_and_source_donor_exclusion():
    _,_,_,states,index=_synthetic()
    local={"Alpha local1","Alpha local2","Alpha local3","Alpha local4"}
    external={"Alpha out1","Alpha out2","Alpha out3"}
    before=make_donor_models(local,external,set(),states,index)
    after=make_donor_models(local,external,{"Alpha local1","Alpha out1","Alpha out2"},
                            states,index)
    assert before["local_names"]==4 and before["external_names"]==3
    assert after["local_names"]==3 and after["external_names"]==1
    assert after["n_local_excluded"]==1 and after["n_external_excluded"]==2
    assert np.isclose(before["prior"].sum(),1)
    assert np.isclose(after["prior"].sum(),1)
    assert after["prior"][1]<before["prior"][1]
    with pytest.raises(ValueError,match="not disjoint"):
        make_donor_models(local,local,set(),states,index)


def test_score_baseline_and_strict_ON_IDENTICAL_focal_accepted_species():
    flora,_,_,states,index=_synthetic()
    local={"Alpha local1","Alpha local2","Alpha local3","Alpha local4"}
    external={"Alpha out1","Alpha out2","Alpha out3"}
    old=make_donor_models(local,external,set(),states,index)
    strict=make_donor_models(local,external,{"Alpha local1","Alpha local4",
                                              "Alpha out1"},states,index)
    score=score_paired_identical_original_species(
        flora["focal"],flora["focal"],states,index,old,strict)
    assert score["common_target_direct_carrier_count"]==5
    assert score["n_different_target_genera"]==1
    assert score["strict_local_source_names"]==2
    assert score["strict_external_source_names"]==2
    assert score["strict_minus_baseline_local_vs_external_Brier"]==pytest.approx(
        score["strict_local_vs_external_Brier"] -
        score["baseline_on_identical_carriers_local_vs_external_Brier"])
    flipped=dict(states)
    for sp in flora["focal"]:flipped[sp]="red"
    y=score_paired_identical_original_species(
        flora["focal"],flora["focal"],flipped,index,old,strict)
    assert np.isfinite(y["strict_local_vs_external_logloss"])
    assert score["strict_local_vs_external_Brier"]>0
    assert y["strict_local_vs_external_Brier"]<0
    too_strict=make_donor_models(local,external,set(local),states,index)
    assert too_strict is None
    assert score_paired_identical_original_species(
        flora["focal"],flora["focal"],states,index,old,too_strict) is None


def test_anonymous_source_physical_and_first_source_group_equal_scores():
    rows=[]
    for i in range(60):
        rows.append(dict(spatial_block=f"B{i//2}",
           source_cluster=f"E{i//3}",
           strict_local_vs_external_Brier=-.04,
           baseline_on_identical_carriers_local_vs_external_Brier=.03,
           strict_minus_baseline_local_vs_external_Brier=-.07,
           strict_local_vs_external_logloss=-.06,
           baseline_on_identical_carriers_local_vs_external_logloss=.02,
           strict_minus_baseline_local_vs_external_logloss=-.08,
           strict_other_region_genus_vs_blind_Brier=.015))
    df=pd.DataFrame(rows)
    for group in ("spatial_block","source_cluster"):
        a=original_group_equal_stress(df,group,seed=20261119,nboot=99)
        assert a["status"]=="ORIGINAL_STRICT_TRAINING_EXCLUSION_PAIRED_FOCAL_NAME_HOLDOUT"
        assert a["n_scored_island_floras"]==60
        assert a["strict_local_vs_external_Brier_group_equal_mean"]==pytest.approx(-.04)
        assert a["strict_minus_baseline_local_vs_external_Brier_group_resample_95"]==pytest.approx([-.07,-.07])
        assert a["strict_local_vs_external_logloss_group_resample_95"]==pytest.approx([-.06,-.06])
    insufficient=original_group_equal_stress(df.iloc[:9],
                "spatial_block",seed=3,nboot=99)
    assert insufficient["status"]=="NOT_IDENTIFIABLE_ORIGINAL_GROUP_SUPPORT"
    assert insufficient["group_resampling_not_performed"]


def test_full_original_source_flora_group_integrity():
    flora,blocks,sources,_,_=_synthetic()
    with pytest.raises(ValueError,match="not aligned"):
        original_full_flora_group_blacklists(flora,blocks,dict(sources,focal=None,
                                                               external=None,unknown="A1"))
    with pytest.raises(ValueError,match="empty"):
        original_full_flora_group_blacklists(dict(flora,extra=set()),dict(blocks,extra="B5"),
                                             dict(sources,extra="A5"))
