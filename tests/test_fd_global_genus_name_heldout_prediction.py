"""Synthetic source-name exclusion and genus-blind predictive control tests."""
import numpy as np
import pandas as pd
import pytest

from fd_global_genus_name_heldout_prediction import (
    compile_donors,predict_original_flora,block_bootstrap_paired_prediction)


def make_state():
    return {
        "Alpha focal":"white",
        "Alpha hidden":"white",
        "Alpha donor":"red",
        "Beta donor":"white",
        "Beta other":"red",
    }


def test_full_flora_taxon_name_out_including_NONcarrier_trait_names():
    states=make_state()
    category={"red":0,"white":1}
    donors=compile_donors(states,set(states),category)
    full={"Alpha focal","Alpha hidden","Nontrait species"}
    target={"Alpha focal"}
    z=predict_original_flora(full,target,states,donors)
    assert z["n_all_original_flora_names_excluded"]==2
    assert z["n_unique_external_same_ledger_species_donors_after_name_holdout"]==3
    assert z["share_carriers_with_independent_NAME_genus_donor"]==1.
    assert z["share_carriers_with_at_least_two_distinct_genus_profiles"]==0.
    # Changing the trait of the non-target but present island species can
    # NEVER change prediction because it must be excluded by full flora name.
    changed=dict(states,**{"Alpha hidden":"red"})
    p2=predict_original_flora(full,target,changed,
                             compile_donors(changed,set(changed),category))
    for key in ("blind_category_frequency_D2",
                "name_heldout_genus_category_frequency_D2",
                "name_heldout_genus_prediction_gain_D2"):
        assert z[key]==pytest.approx(p2[key])


def test_no_congener_donor_has_exact_prior_and_zero_improvement():
    state={"Unknown focal":"red","Alpha donor":"white","Beta donor":"red"}
    categories={"red":0,"white":1}
    donor=compile_donors(state,{"Alpha donor","Beta donor"},categories)
    z=predict_original_flora({"Unknown focal"},{"Unknown focal"},state,donor)
    assert z["share_carriers_with_independent_NAME_genus_donor"]==0.
    assert z["name_heldout_genus_prediction_gain_D2"]==pytest.approx(0.)
    assert z["blind_category_frequency_D2"]==pytest.approx(
        z["name_heldout_genus_category_frequency_D2"])


def test_category_pair_is_one_REAL_accepted_species_state():
    state={"A first":("white","tube"),"A second":("red","open"),
           "B first":("red","open")}
    categories={c:i for i,c in enumerate(sorted(set(state.values()),key=repr))}
    donors=compile_donors(state,set(state),categories)
    z=predict_original_flora({"A first"},{"A first"},state,donors)
    assert z["n_direct_carriers"]==1
    assert z["n_all_original_flora_names_excluded"]==1
    assert z["n_unique_external_same_ledger_species_donors_after_name_holdout"]==2


def test_negative_prediction_gains_remain_possible_and_not_filtered():
    state={"A test":"white","A wrong":"red","B common":"white"}
    donor=compile_donors(state,set(state),{"red":0,"white":1})
    z=predict_original_flora({"A test"},{"A test"},state,donor)
    assert z["name_heldout_genus_prediction_gain_D2"]<0


def test_cluster_equal_predictive_gain_and_source_identity_support():
    rows=[]
    for i in range(60):
        rows.append(dict(physical_block=f"Block{i//2}",
             name_heldout_genus_prediction_gain_D2=.025,
             blind_category_frequency_D2=.1,
             name_heldout_genus_category_frequency_D2=.075))
    q=block_bootstrap_paired_prediction(
        pd.DataFrame(rows),group="physical_block",n_boot=99)
    assert q["status"]=="ORIGINAL_GROUP_EQUAL_NAME_HELDOUT_GENUS_PREDICTION"
    assert q["n_scored_islands"]==60
    assert q["n_groups"]==30
    assert q["genus_vs_prior_group_equal_gain_95"]==pytest.approx([.025,.025])
    assert q["source_conditional_group_equal_relative_frequency_loss_gain"]==pytest.approx(.25)
    assert block_bootstrap_paired_prediction(pd.DataFrame(rows[:10]),
        group="physical_block",n_boot=99)["status"]==(
        "INSUFFICIENT_ORIGINAL_INDEPENDENT_SOURCE_CLUSTERS")


def test_prohibited_source_missing_states_and_impossible_heldout_roster():
    states=make_state()
    donors=compile_donors(states,set(states),{"red":0,"white":1})
    with pytest.raises(ValueError):
        predict_original_flora({"Alpha hidden"},{"Alpha focal"},states,donors)
    with pytest.raises(ValueError):
        predict_original_flora({"Alpha focal"},{"Alpha focal"},states,donors,
                               alpha=0)
