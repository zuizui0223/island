"""Pure synthetic name-heldout species-direct real floral profile prediction."""
import numpy as np
import pandas as pd
import pytest

from fd.source_only_real_joint_congener_prediction import (
    SCOPES, SECTORS, DONORS, REGIMES,
    predict_joint_by_source_scopes, aggregate_source_congener_prediction
)


def build():
    full={f"Alpha island{i}" for i in range(6)}
    full|={f"Beta island{i}" for i in range(4)}
    full|={"Gamma island0"}
    geo={f"Alpha geo{i}" for i in range(5)}
    geo|={f"Beta geo{i}" for i in range(5)}
    geo|={f"Gamma geo{i}" for i in range(15)}
    climate={f"Alpha climate{i}" for i in range(6)}
    climate|={f"Beta climate{i}" for i in range(6)}
    climate|={f"Gamma climate{i}" for i in range(16)}
    allsp=full|geo|climate
    a=set(x for x in allsp if x.startswith("Alpha"))
    b=set(a)
    # Beta and Gamma are 0 in the true same-species A×B syndrome state.
    profile={"known":frozenset(allsp),
             "A":frozenset(a),"B":frozenset(b),
             "joint":frozenset(a)}
    rosters={
        "all_recorded":full,
        "strict_native":full-{"Gamma island0"},
        "strict_native_nonendemic":full-{"Gamma island0","Beta island3"},
    }
    return full,rosters,geo,climate,profile


def test_source_names_not_target_species_and_genus_is_true_forecast():
    full,rosters,geo,climate,profile=build()
    z=predict_joint_by_source_scopes(
        full_island=full,scope_rosters=rosters,
        geo5_source=geo,climate10_source=climate,profile=profile
    )
    a=z["all_recorded"]["any_congener"]
    assert a["geo5"]["n_target_direct"]==a["climate10"]["n_target_direct"]==11
    assert a["geo5"]["n_target_genera"]==3
    # All recorded genus labels have true different-name source congener donors.
    assert np.isclose(a["geo5"]["joint_observed"],6/11)
    assert np.isclose(a["geo5"]["joint_genus_predicted"],6/11)
    assert a["geo5"]["brier_genus"]==0.
    assert a["geo5"]["brier_factorized"]==0.
    assert a["geo5"]["n_source_genera_with_both_traits_variable"]==0
    assert a["geo5"]["brier_pooled"]>0
    assert a["geo5"]["brier_global"]>0
    assert np.isclose(a["geo5"]["coupling_genus_predicted"],6/11-(6/11)**2)
    assert a["geo5"]["n_target_direct"]==a["climate10"]["n_target_direct"]
    assert z["strict_native"]["at_least_two_congeners"]["geo5"]["n_target_direct"]==10


def test_source_species_focal_name_overlap_rejected():
    full,rosters,geo,climate,profile=build()
    with pytest.raises(ValueError,match="leaked"):
        predict_joint_by_source_scopes(
            full_island=full,scope_rosters=rosters,
            geo5_source=geo|{"Alpha island1"},
            climate10_source=climate,profile=profile
        )


def test_missing_congener_is_unidentified_not_zero_prediction():
    full,rosters,geo,climate,profile=build()
    x=predict_joint_by_source_scopes(
        full_island=full,scope_rosters=rosters,
        geo5_source=geo,climate10_source=set(),
        profile=profile,
    )
    q=x["all_recorded"]["any_congener"]["geo5"]
    assert q["n_target_direct"]==0
    assert np.isnan(q["joint_genus_predicted"])
    assert np.isnan(q["coupling_observed"])


def _frame(n=40,missing=False):
    full,rosters,geo,climate,profile=build()
    z=predict_joint_by_source_scopes(
        full_island=full,scope_rosters=rosters,
        geo5_source=geo,climate10_source=climate,profile=profile)
    out=[]
    for i in range(n):
        for scope in SCOPES:
            for sector in SECTORS:
                for gate in DONORS:
                    info=z[scope][gate]
                    row=dict(
                        island_id=f"synthetic_{i}",regime="tropical",
                        scope=scope,sector=sector,donor_gate=gate,
                        spatial_block=f"b{i//2}",
                        source_signature=f"source{(i//2)%6}",
                        n_target_direct=info["geo5"]["n_target_direct"],
                        n_target_genera=info["geo5"]["n_target_genera"],
                        island_direct_coverage_supported=1.,
                        geo5_direct_coverage_supported=1.,
                        climate10_direct_coverage_supported=1.,
                    )
                    for mode in ("geo5","climate10"):
                        for k,value in info[mode].items():
                            if k.startswith("n_source_matched_total") or k.startswith("target_species_per"):
                                continue
                            row[f"{mode}__{k}"]=value
                        row[f"{mode}__n_source_matched_total"]=25 if mode=="geo5" else 28
                    out.append(row)
    data=pd.DataFrame(out)
    if missing:
        choose=data.island_id.eq("synthetic_0")
        data.loc[choose,"n_target_direct"]=0
        data.loc[choose,"n_target_genera"]=0
        for mode in ("geo5","climate10"):
            data.loc[choose,f"{mode}__n_target_positive"]=0
            for k in ("joint_observed","coupling_observed",
                      "joint_genus_predicted","coupling_genus_predicted",
                      "joint_pooled_predicted","coupling_pooled_predicted",
                      "joint_global_predicted","coupling_global_predicted",
                      "brier_genus","brier_pooled","brier_global"):
                data.loc[choose,f"{mode}__{k}"]=np.nan
    return data


def test_same_original_island_score_and_improvement_over_equal_pool():
    x=aggregate_source_congener_prediction(_frame(),population=40,draws=299)
    assert x["n_region_scope_sector_donor_panels"]==5*3*3*2
    q=next(z for z in x["panels"] if z["region"]=="tropical"
           and z["scope"]=="all_recorded"
           and z["genuine_same_species_joint_sector"]=="actinomorphic_deep"
           and z["congener_donor_requirement"]=="at_least_two_congeners")
    assert q["n_same_islands_with_both_source_congeners"]==40
    assert q["n_coverage_balanced_10pp"]==40
    assert q["median_supported_original_island_genera"]==3
    a=q["all__spatial_block"]
    assert a["status"]=="SOURCE_ONLY_NAME_HELDOUT_CONGENERIC_PREDICTION"
    assert np.isclose(a["geo5__bias_genus__mean"],0.)
    assert a["geo5__brier_gain_vs_pooled__mean"]>0
    assert a["geo5__brier_gain_vs_global__mean"]>0
    assert np.isclose(a["geo5__coupling_residual_congener__mean"],0.)


def test_missing_trait_source_not_imputed_and_original_flora_preserved():
    q=aggregate_source_congener_prediction(_frame(missing=True),population=40,draws=299)
    v=next(z for z in q["panels"] if z["region"]=="tropical"
           and z["scope"]=="all_recorded"
           and z["genuine_same_species_joint_sector"]=="zygomorphic_deep"
           and z["congener_donor_requirement"]=="any_congener")
    assert v["n_same_islands_with_both_source_congeners"]==39
    assert q["n_original_GIFT_botanical_flora_islands"]==40


def test_invalid_duplicate_or_missing_metric_fails_closed():
    x=_frame()
    with pytest.raises(ValueError,match="duplicate"):
        aggregate_source_congener_prediction(
            pd.concat([x,x.iloc[:1]]),population=40,draws=299)
    x=_frame()
    x.loc[0,"geo5__joint_observed"]=np.nan
    with pytest.raises(ValueError,match="missing"):
        aggregate_source_congener_prediction(x,population=40,draws=299)


def test_documentation_imbalance_is_reported_not_treated_as_source_absence():
    x=_frame()
    x["geo5_direct_coverage_supported"]=.1
    report=aggregate_source_congener_prediction(x,population=40,draws=299)
    z=next(k for k in report["panels"] if k["region"]=="tropical"
           and k["scope"]=="all_recorded"
           and k["genuine_same_species_joint_sector"]=="actinomorphic_deep"
           and k["congener_donor_requirement"]=="any_congener")
    assert z["n_same_islands_with_both_source_congeners"]==40
    assert z["n_coverage_balanced_10pp"]==0
    assert z["coverage10__spatial_block"]["status"]=="NOT_IDENTIFIABLE_GEOGRAPHICAL_SOURCE_SUPPORT"


def test_genuine_source_joint_label_linkage_adds_predictive_skill_over_factorized_marginals():
    full,rosters,geo,climate,profile=build()
    observed_A=set(profile["A"])
    observed_B=set(profile["B"])
    for source in (geo,climate):
        alpha=sorted(sp for sp in source if sp.startswith("Alpha"))
        # Independently recorded donor A/B marginals are both variable;
        # 3 of 5 and 3 of 6 sources have both positive states, so
        # P(A and B)>P(A)*P(B) is a genuine different-name donor pattern.
        for sp in alpha[3:]:
            observed_A.remove(sp)
            observed_B.remove(sp)
    prof=dict(profile)
    prof["A"]=frozenset(observed_A)
    prof["B"]=frozenset(observed_B)
    prof["joint"]=prof["A"]&prof["B"]
    results=predict_joint_by_source_scopes(
        full_island=full,scope_rosters=rosters,
        geo5_source=geo,climate10_source=climate,profile=prof,
    )
    for mode in ("geo5","climate10"):
        p=results["all_recorded"]["any_congener"][mode]
        assert p["n_source_genera_with_both_traits_variable"]>=1
        assert p["n_target_carriers_on_variable_source_genera"]>=6
        assert p["joint_genus_predicted"]>p["joint_factorized_predicted"]
        assert p["brier_factorized"]>p["brier_genus"]


def test_no_syndrome_predictive_increment_when_source_genus_margins_are_constant():
    full,rosters,geo,climate,profile=build()
    q=predict_joint_by_source_scopes(
        full_island=full,scope_rosters=rosters,
        geo5_source=geo,climate10_source=climate,profile=profile
    )["all_recorded"]["any_congener"]
    for mode in ("geo5","climate10"):
        assert q[mode]["n_source_genera_with_both_traits_variable"]==0
        assert np.isclose(q[mode]["brier_factorized"],q[mode]["brier_genus"])
        assert np.isclose(
            q[mode]["joint_factorized_predicted"],
            q[mode]["joint_genus_predicted"]
        )
