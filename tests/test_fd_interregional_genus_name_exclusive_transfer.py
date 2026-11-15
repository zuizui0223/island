"""Synthetic exact original regional source taxon exclusions and matched donors."""
import numpy as np
import pandas as pd
import pytest

from fd_interregional_genus_name_exclusive_transfer import (
    donor_source_region_ownership,hashed_donors_by_genus,
    reference_profile_prior,balanced_probabilities,
    island_heldout_comparison,group_equal_pair_stats,
)


def _toy():
    states={
      "A focal":"blue","A local1":"blue","A local2":"blue",
      "A outside1":"red","A outside2":"red",
      "B focal":"red","B local1":"red","B outside1":"blue",
      "C focal":"yellow"}
    idx={"blue":0,"red":1,"yellow":2}
    prior=reference_profile_prior({"A outside1","A outside2","B outside1"},
                                   states,idx)
    local=hashed_donors_by_genus(
        {"A focal","A local1","A local2","B focal","B local1"},seed=20261115)
    out=hashed_donors_by_genus(
        {"A outside1","A outside2","B outside1"},seed=20261115)
    return states,idx,prior,local,out


def test_source_region_exclusivity_strictly_removes_names_in_both():
    source={
        "northern_midlatitude":{"A local","D shared","B other"},
        "tropical":{"A outside","D shared"},
        "southern_extratropical":{"C outside","D shared"},
        "northern_high_latitude":{"E outside"},
        "geography_unresolved":{"A unclear"}}
    one,two=donor_source_region_ownership(source,"northern_midlatitude")
    assert one=={"A local","B other"}
    assert two=={"A outside","C outside","E outside","A unclear"}-{"A unclear"}
    assert not one&two
    assert "D shared" not in one|two
    assert "A unclear" not in two
    with pytest.raises(ValueError):
        donor_source_region_ownership(source,"not_real")


def test_equal_genus_training_name_counts_after_focal_full_flora_exclusion():
    states,idx,prior,local,other=_toy()
    l,e,q=balanced_probabilities(
        "A",{"A focal","A local1"},local,other,states,idx,prior)
    assert q["match_n"]==1
    assert q["local_names"]==1
    assert q["external_names"]==2
    assert l[0]>e[0]
    assert l[1]<e[1]
    assert np.isclose(l.sum(),1) and np.isclose(e.sum(),1)
    assert local["A"]==hashed_donors_by_genus(
        set(local["A"]+local["B"]),seed=20261115)["A"]


def test_focal_carriers_observed_and_unmatched_genus_fallback_not_zero():
    states,idx,prior,local,other=_toy()
    whole={"A focal","B focal","C focal"}
    target=whole
    z=island_heldout_comparison(whole,target,states,idx,prior,local,other)
    assert z["n_original_real_direct_carriers"]==3
    assert z["n_matched_genus_direct_carriers"]==2
    assert z["matched_carrier_fraction"]==pytest.approx(2/3)
    assert z["matched_genus_donor_names_sampled_per_side"]>=2
    assert np.isfinite(z["local_minus_external_proper_logloss_gain_per_comparable_plant"])
    assert np.isfinite(z["local_minus_external_proper_Brier_gain_per_comparable_plant"])
    assert np.isfinite(z["external_genus_vs_blind_proper_logloss_gain_per_comparable_plant"])
    assert np.isfinite(z["external_genus_vs_blind_proper_Brier_gain_per_comparable_plant"])
    assert z["local_exclusive_frequency_D2"]>=0
    assert z["external_exclusive_frequency_D2"]>=0


def test_no_matched_genus_keeps_local_external_predictions_identical():
    states,idx,prior,local,other=_toy()
    z=island_heldout_comparison({"C focal"},{"C focal"},states,idx,
                               prior,local,other)
    assert z["n_matched_genus_direct_carriers"]==0
    assert z["locally_exclusive_minus_externally_exclusive_frequency_D2_gain"]==pytest.approx(0)
    assert z["local_minus_external_proper_logloss_gain_per_comparable_plant"] is None
    assert z["external_genus_vs_blind_proper_logloss_gain_per_comparable_plant"] is None
    assert z["external_genus_vs_external_blind_frequency_D2_gain"]==pytest.approx(0)
    with pytest.raises(ValueError):
        island_heldout_comparison({"A local1"},{"A focal"},states,idx,prior,local,other)


def test_negative_and_positive_local_vs_external_evidence_not_filtered():
    states,idx,prior,local,other=_toy()
    z=island_heldout_comparison({"A focal"},{"A focal"},states,idx,prior,local,other)
    assert z["local_minus_external_proper_logloss_gain_per_comparable_plant"]>0
    flipped=dict(states)
    flipped["A focal"]="red"
    x=island_heldout_comparison({"A focal"},{"A focal"},flipped,idx,prior,local,other)
    assert x["local_minus_external_proper_logloss_gain_per_comparable_plant"]<0


def test_group_equal_original_block_contrast_and_source_support_gates():
    rows=[]
    for i in range(60):
        rows.append(dict(spatial_block=f"block{i//2}",
            locally_exclusive_minus_externally_exclusive_frequency_D2_gain=.03,
            local_minus_external_proper_logloss_gain_per_comparable_plant=.02,
            local_minus_external_proper_Brier_gain_per_comparable_plant=.01,
            external_genus_vs_blind_proper_logloss_gain_per_comparable_plant=.015,
            external_genus_vs_blind_proper_Brier_gain_per_comparable_plant=.008,
            external_genus_vs_external_blind_frequency_D2_gain=.019))
    q=group_equal_pair_stats(pd.DataFrame(rows),"spatial_block",
                              seed=20261115,nboot=99)
    assert q["status"]=="ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON"
    assert q["n_scored_islands"]==60
    assert q["n_clusters"]==30
    assert q["locally_exclusive_minus_externally_exclusive_frequency_D2_gain_group_equal_mean"]==pytest.approx(.03)
    assert q["local_minus_external_proper_logloss_gain_per_comparable_plant_group_resample_95"]==pytest.approx([.02,.02])
    assert q["external_genus_vs_blind_proper_logloss_gain_per_comparable_plant_group_equal_mean"]==pytest.approx(.015)
    no=group_equal_pair_stats(pd.DataFrame(rows[:12]),"spatial_block",
                              seed=20261115,nboot=99)
    assert no["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_ORIGINAL_GROUP_SUPPORT"
    assert no["n_scored_islands"]==12
