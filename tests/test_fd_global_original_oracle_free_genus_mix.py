"""Synthetic true source-donor genus-mix without observed focal island genus oracle."""
import numpy as np
import pandas as pd
import pytest

from fd_global_original_oracle_free_genus_mix import (
    _loss_to_plant_states,
    refit_no_oracle_regional_genus_source,
    score_original_island_without_genus_input,
    group_equal_original_predictive_gains,
    COLUMNS,QUALIFIED,
)


def toy():
    state={}
    # Region-local genus is A-heavy, region-external taxonomy is B-heavy.
    local={f"A local{i}" for i in range(7)}|{"B local1"}
    external={"A outside1","A outside2"}|{
        f"B outside{i}" for i in range(9)}
    for s in local:
        state[s]="blue" if s.startswith("A ") else "red"
    for s in external:
        state[s]="blue" if s.startswith("A ") else "red"
    for i in range(6):state[f"A focal{i}"]="blue"
    for i in range(2):state[f"B focal{i}"]="red"
    return local,external,state,{"blue":0,"red":1}


def test_no_oracle_regional_source_pool_weights_cannot_see_focal_genus():
    local,external,states,idx=toy()
    excluded={"A focal0","B focal0","A nearest_group"}
    model=refit_no_oracle_regional_genus_source(
        local,external,excluded,states,idx)
    assert model["original_focal_island_genus_identity_NOT_USED_IN_MODELLING"]
    assert model["n_distinct_local_direct_source_species"]==8
    assert model["n_distinct_external_direct_source_species"]==11
    assert model["n_distinct_genus_profiles_with_external_direct_congeners"]==2
    # Original genus-conditional estimates are explicitly shrunk toward the
    # outside-region blind prior (ALPHA=5); the source mix need not exceed
    # 50% even when it shifts BLUE probability in the correct direction.
    assert (model["source_genus_weighted_external_flower_profile"][0]>
            model["external_genus_blind_prior"][0])
    assert model["external_genus_blind_prior"][0]<.5
    assert np.isclose(model["source_genus_weighted_external_flower_profile"].sum(),1)
    assert model["local_source_species_fraction_with_external_same_genus_donor"]==1.
    assert "focal" not in " ".join(model.keys())


def test_no_oracle_prediction_all_real_direct_carriers_and_negative_results():
    local,external,state,idx=toy()
    model=refit_no_oracle_regional_genus_source(local,external,
        {"A focal0","A focal1","A nearest_group"},state,idx)
    positive={f"A focal{i}" for i in range(6)}
    p=score_original_island_without_genus_input(positive,positive,state,idx,model)
    assert p["n_original_direct_plant_carriers"]==6
    assert p["genus_mix_vs_external_blind_Brier"]>0
    assert p["genus_mix_vs_external_blind_logloss"]>0
    # Same no-oracle model probabilities: swapping the focal accepted taxon
    # genus names without changing observed TRUE category cannot influence
    # any no-oracle score, while the explicitly oracle-known comparator can.
    state_other=dict(state)
    new_species={f"Z focal{i}" for i in range(6)}
    for sp in new_species:state_other[sp]="blue"
    q=score_original_island_without_genus_input(
        new_species,new_species,state_other,idx,model)
    for key in (
      "genus_mix_vs_external_blind_Brier",
      "genus_mix_vs_external_blind_logloss",
      "genus_mix_vs_local_blind_Brier",
      "genus_mix_vs_local_blind_logloss",
    ):
        assert p[key]==pytest.approx(q[key])
    assert p["oracle_known_genus_vs_genus_mix_Brier"]!=pytest.approx(
        q["oracle_known_genus_vs_genus_mix_Brier"])
    bad={f"B focal{i}" for i in range(2)}
    mixed=positive|bad
    y=score_original_island_without_genus_input(mixed,mixed,state,idx,model)
    assert y["n_original_direct_plant_carriers"]==8
    negative={f"B test{i}" for i in range(5)}
    neg_state=dict(state)
    neg_state.update({sp:"red" for sp in negative})
    z=score_original_island_without_genus_input(
        negative,negative,neg_state,idx,model)
    assert z["genus_mix_vs_external_blind_Brier"]<0
    assert z["genus_mix_vs_external_blind_logloss"]<0


def test_no_oracle_model_falls_back_for_genus_without_outside_source():
    local,external,states,idx=toy()
    local|={"C local1"}
    states["C local1"]="red"
    m=refit_no_oracle_regional_genus_source(
        local,external,{"A focal0"},states,idx)
    assert m["n_distinct_genus_profiles_with_external_direct_congeners"]==2
    assert np.allclose(m["external_conditional_genus_profiles"]["C"],
                       m["external_genus_blind_prior"])
    with pytest.raises(ValueError,match="region-exclusive"):
        refit_no_oracle_regional_genus_source(
            local,local,{"A focal0"},states,idx)
    assert refit_no_oracle_regional_genus_source(local,external,
           set(local),states,idx) is None


def test_original_source_bootstraps_match_proper_losses_and_preserve_bad_signs():
    rows=[]
    for k in range(60):
        row=dict(spatial_block=f"b{k//2}",first_source=f"s{k//3}")
        for col in COLUMNS:
            row[col]=-.03 if col.startswith("genus_mix_vs_local") else .02
        rows.append(row)
    df=pd.DataFrame(rows)
    for group in ("spatial_block","first_source"):
        q=group_equal_original_predictive_gains(
            df,group,seed=20261119,nboot=99)
        assert q["status"]==QUALIFIED
        assert q["genus_mix_vs_external_blind_Brier_group_equal_mean"]==pytest.approx(.02)
        assert q["genus_mix_vs_external_blind_logloss_group_resample_95"]==pytest.approx([.02,.02])
        assert q["genus_mix_vs_local_blind_Brier_group_resample_95"]==pytest.approx([-.03,-.03])
        assert q["original_group_bootstrap_exploratory_not_multiplicity_adjusted"]
    q=group_equal_original_predictive_gains(df.iloc[:9],"spatial_block",
                                              seed=0,nboot=99)
    assert q["status"]=="NOT_IDENTIFIABLE_ORIGINAL_SOURCE_BLOCK_SUPPORT"


def test_original_unobserved_genus_and_category_are_not_imputed():
    local,external,states,idx=toy()
    model=refit_no_oracle_regional_genus_source(
        local,external,{"A focal0"},states,idx)
    with pytest.raises(ValueError,match="not on accepted flora"):
        score_original_island_without_genus_input(
            {"A focal1"},{"A focal0"},states,idx,model)
    with pytest.raises(ValueError,match="category"):
        _loss_to_plant_states(np.array([.5,.5]),[4])
