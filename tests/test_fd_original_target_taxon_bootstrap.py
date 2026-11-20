"""Synthetic source-private GIFT TARGET species/genus resampling and proper loss checks.

NO original plant, GIFT island, donor group or geographic ID is released.
"""
import numpy as np
import pytest

from fd_original_stage2_name_disjoint_flower_forecast import (
    TARGET_GAINS, target_species_gain_matrix,
    distinct_target_cluster_summary, taxon_dual_loss_positive,
)
from fd.two_stage_floral_forecast import proper_losses


def test_each_plant_name_proper_brier_and_log_gain_averages_back_to_source_island():
    codes=np.array([0,1,0,2,2,1,0],dtype=int)
    model=np.array([.45,.30,.25])
    refs={
        "regional_blind":np.array([.25,.50,.25]),
        "external_blind":np.array([.33,.33,.34]),
        "regional_genus_prior":np.array([.4,.4,.2]),
    }
    v=target_species_gain_matrix(codes,model,refs)
    assert v.shape==(len(codes),len(TARGET_GAINS))==(7,6)
    for j,base in enumerate(("regional_blind","external_blind","regional_genus_prior")):
        b=proper_losses(refs[base],codes)
        m=proper_losses(model,codes)
        np.testing.assert_allclose(v.mean(axis=0)[2*j:2*j+2],
                                   np.asarray(b)-np.asarray(m),atol=1e-12)


def test_a_single_plant_recorded_on_100_islands_counts_once_in_species_equal():
    accumulated={}
    for g in range(12):
        for i in range(3):
            name=f"Genus{g} species{i}"
            accumulated[name]=np.r_[np.full(len(TARGET_GAINS),.3),1.]
    frequent="Genus0 extremelyfrequent"
    accumulated[frequent]=np.r_[np.full(len(TARGET_GAINS),-30.0),100.0]
    # 36 species + 1 highly replicated = 37 true species, 12 genera.
    # The repeated species has 100 evaluations each showing -0.3, not
    # "100 independent negative plant phenotypes".
    result=distinct_target_cluster_summary(accumulated,nboot=299,seed=20261120)
    assert result["n_distinct_target_species"]==37
    assert result["n_distinct_target_genera"]==12
    assert result["n_species_island_occurrences"]==136
    assert result["target_accepted_names_absent_from_export"]
    assert result["species_equal"]["supported_for_inference"]
    assert result["genus_equal"]["supported_for_inference"]
    for baseline in ("regional_blind","external_blind","regional_genus_prior"):
        assert taxon_dual_loss_positive(result["species_equal"],baseline)
        assert taxon_dual_loss_positive(result["genus_equal"],baseline)
    first=f"environment_source_nearest_gain_vs_regional_blind_Brier"
    assert result["species_equal"][first]["group_equal_mean"]==pytest.approx(
        (36*.3 + (-30/100))/37)
    assert result["genus_equal"][first]["group_equal_mean"]>0
    assert "Genus0" not in str(result) and frequent not in str(result)


def test_no_original_target_species_or_genus_donor_support_is_not_zero_effect():
    empty=distinct_target_cluster_summary({},nboot=99)
    assert empty["n_distinct_target_species"]==0
    assert not empty["species_equal"]["supported_for_inference"]
    assert not empty["genus_equal"]["supported_for_inference"]
    assert not taxon_dual_loss_positive(empty["species_equal"],"regional_blind")
    with pytest.raises(ValueError,match="bootstrap"):
        distinct_target_cluster_summary({},nboot=10)


def test_target_category_ontology_guard_and_missing_source_probabilities():
    code=np.array([0,1,2])
    model=np.array([.3,.3,.4])
    a={"regional_blind":np.array([.5,.5]),
       "external_blind":np.array([.2,.3,.5]),
       "regional_genus_prior":np.array([.4,.4,.2])}
    with pytest.raises(ValueError,match="ontology"):
        target_species_gain_matrix(code,model,a)
    a["regional_blind"]=np.array([.5,.4,.1])
    a["external_blind"]=np.array([.2,float("nan"),.8])
    with pytest.raises(ValueError,match="ontology"):
        target_species_gain_matrix(code,model,a)


def test_target_private_species_groups_have_distinct_estimator():
    data={}
    # Nine accepted species with 30 repeat occurrences each are insufficient
    # for an independent unique-species source bootstrap, even if n=270.
    for i in range(9):
        data[f"Genus{i} t"]=np.r_[np.full(len(TARGET_GAINS),.1*30),30.]
    out=distinct_target_cluster_summary(data,nboot=99)
    assert out["n_distinct_target_species"]==9
    assert out["n_species_island_occurrences"]==270
    assert not out["species_equal"]["supported_for_inference"]
    assert not out["genus_equal"]["supported_for_inference"]
