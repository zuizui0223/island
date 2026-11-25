"""Synthetic-only genus concentration stress tests for numeric flower tube."""
import numpy as np
import pandas as pd
import pytest

from fd.nhm_numeric_genus_concentration import (
    score_genus_sufficient_stats,audit_genus_concentration,
    SCOPES,TIERS,
)


GENERA=("Alpha","Beta","Gamma","Delta","Epsilon")


def one_island_genus_rows(i=0,*,biased=False):
    out=[]
    for scope in SCOPES:
        for tier in TIERS:
            for j,g in enumerate(GENERA):
                n=4 if j==0 else 3
                true_deep=n if j in (0,3) else 0
                pa=.8 if true_deep else .2
                rec=dict(
                    island_id=f"island{i}",regime="tropical",
                    scope=scope,tier=tier,genus=g,
                    spatial_block=f"block{i//2}",
                    first_source=f"src{(i//2)%11}",
                    n_target=n,n_target_deep=true_deep,
                    n_target_original_genus=(n*10 if biased else n),
                )
                for mode in ("geo5","climate10"):
                    rec[mode+"__n_numeric"]=5
                    rec[mode+"__n_original_source_genus"]=(
                        100 if biased else 5
                    )
                    rec[mode+"__genus_pdeep"]=pa
                out.append(rec)
    return pd.DataFrame(out)


def sample(n=40,**kwargs):
    return pd.concat([
        one_island_genus_rows(i,**kwargs) for i in range(n)
    ],ignore_index=True)


def pop(n=40):
    return pd.DataFrame(dict(
        island_id=[f"island{i}" for i in range(n)],
        regime=["tropical"]*n,
    ))


def test_per_island_genus_equal_and_congener_pooled_exact_brier_arithmetic():
    f=one_island_genus_rows()
    q=f.loc[f.scope.eq("all_recorded")&
            f.tier.eq(TIERS[0])]
    d=score_genus_sufficient_stats(q)
    assert d["status"]=="ELIGIBLE_PAIRED_SOURCE_GENUS_CONCENTRATION"
    assert d["n_target"]==16
    assert d["n_genera"]==5
    assert d["largest_genus_carrier_fraction"]==.25
    assert d["effective_genus_number"]>4.
    assert d["documentation_10pp_both"]
    assert d["source_2plus_genus_both"]
    target=np.repeat(q.n_target_deep.to_numpy(float)/
                     q.n_target.to_numpy(float),q.n_target.to_numpy(int))
    pred=np.repeat(q.geo5__genus_pdeep.to_numpy(float),
                   q.n_target.to_numpy(int))
    pool=float(np.average(q.geo5__genus_pdeep,weights=q.geo5__n_numeric))
    assert np.isclose(
        d["geo5__brier_gain_species_equal"],
        np.mean((target-pool)**2-(target-pred)**2)
    )
    assert d["geo5__brier_gain_genus_equal"]>0
    assert np.isclose(d["geo5__brier_gain_species_equal"],
                      d["climate10__brier_gain_species_equal"])


def test_outcome_blind_top_one_and_local_largest_removal_preserve_spatial_groups():
    z=audit_genus_concentration(sample(),pop(),population=40,bootstrap=299)
    assert z["original_global_botanical_floras"]==40
    assert z["n_original_region_scope_corpus_panels"]==5*3*2
    q=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded" and p["tier"]==TIERS[0])
    assert q["n_original_islands_source_numeric_paired_support"]==40
    assert q["n_original_islands_Direct_balance_both_10pp"]==40
    assert q["n_original_islands_effective_three_plus_genera"]==40
    assert q["all_original_eligible__spatial_block"]["status"]=="DESCRIPTIVE_GENUS_CONCENTRATION_SENSITIVITY"
    assert q["remove_top_1_original_carrier_genera__spatial_block"]["n_original_floras"]==40
    assert q["remove_largest_within_each_island__spatial_block"]["n_original_floras"]==40
    assert q["remove_top_5_original_carrier_genera__spatial_block"]["status"]=="NOT_IDENTIFIABLE_GROUP_SUPPORT"
    assert np.isclose(q["all_original_eligible__spatial_block"]["geo5__brier_gain_species_equal__block_mean"],
                      q["all_original_eligible__spatial_block"]["climate10__brier_gain_species_equal__block_mean"])


def test_input_documentation_bias_makes_coverage_sensitivity_fail():
    z=audit_genus_concentration(sample(biased=True),pop(),
                                population=40,bootstrap=299)
    q=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded" and p["tier"]==TIERS[0])
    assert q["n_original_islands_Direct_balance_both_10pp"]==40
    # Both island and source documentation are 10%, hence balanced despite
    # severely low absolute coverage. It is a balance test, not completeness.
    assert q["documentation_balance10__spatial_block"]["n_original_floras"]==40


def test_documentation_discrepancy_not_assumed_zero():
    df=sample()
    df.loc[df.genus.eq("Alpha"),"n_target_original_genus"]=200
    z=audit_genus_concentration(df,pop(),population=40,bootstrap=299)
    q=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded" and p["tier"]==TIERS[0])
    assert q["n_original_islands_Direct_balance_both_10pp"]==0
    assert q["documentation_balance10__spatial_block"]["status"]=="NOT_IDENTIFIABLE_GROUP_SUPPORT"


def test_original_frame_and_private_source_no_invalid_silent_zeros():
    with pytest.raises(ValueError,match="original complete"):
        audit_genus_concentration(sample(),pop(39),population=40,bootstrap=299)
    f=sample()
    with pytest.raises(ValueError,match="duplicate"):
        audit_genus_concentration(
            pd.concat([f,f.iloc[:1]]),pop(),population=40,bootstrap=299)
    f=sample()
    f.loc[0,"geo5__genus_pdeep"]=np.nan
    with pytest.raises(ValueError,match="source numeric"):
        audit_genus_concentration(f,pop(),population=40,bootstrap=299)


def test_source_genus_marginal_score_no_fabricated_target_trait():
    f=one_island_genus_rows()
    x=f.loc[f.scope.eq("all_recorded")&
            f.tier.eq(TIERS[0])].copy()
    x.loc[x.genus.eq("Epsilon"),"n_target_deep"]=4
    with pytest.raises(ValueError,match="Direct"):
        audit_genus_concentration(x,pop(1),population=1,bootstrap=299)
