"""Synthetic NHM numeric tube length source-only genus calibration tests."""
import numpy as np
import pandas as pd
import pytest

from fd.nhm_numeric_source_congener_predict import (
    SCOPES,TIERS,fit_source_only_length_calibrator,
    calibrated_prob,make_genus_heldout_length_calibrator,
    source_predict_same_island,
    aggregate_external_source_numeric_prediction,
)


def training():
    names={f"Source{i//4} training{i}" for i in range(240)}
    values={s: .15+int(s.split("training")[1])/50 for s in names}
    depth={s:("deep" if int(s.split("training")[1])>=120
              else "shallow") for s in names}
    return values,depth,names


def calibration():
    x,y,z=training()
    return fit_source_only_length_calibrator(x,y,z,set())


def test_calibration_is_monotonic_and_no_focal_label_leakage():
    m=calibration()
    assert m["n_calibration_distinct_original_accepted_species"]==240
    assert m["n_original_source_only_deep"]==120
    assert m["log1p_cm_slope"]>0
    assert calibrated_prob(3.5,m)>calibrated_prob(.35,m)
    x,y,z=training()
    notsources=set(s for s in z if int(s.split("training")[1])>=220)
    q=fit_source_only_length_calibrator(x,y,z,notsources)
    assert q["n_calibration_distinct_original_accepted_species"]==220
    assert q["original_island_accepted_species_fully_excluded"]


def test_invalid_source_calibration_and_unmeasured_tube_rejected():
    x,y,z=training()
    with pytest.raises(ValueError,match="insufficient"):
        fit_source_only_length_calibrator(x,y,set(list(z)[:20]),set())
    with pytest.raises(ValueError):
        calibrated_prob(float("nan"),calibration())
    with pytest.raises(ValueError):
        calibrated_prob(-2,calibration())


def example_roster():
    full={"Alpha isle1","Alpha isle2","Alpha isle3",
          "Beta isle1","Beta isle2","Beta isle3"}
    scopes={key:full for key in SCOPES}
    x={s:("deep" if s.startswith("Alpha") else "shallow") for s in full}
    g={f"Alpha geo{i}" for i in range(12)} | {f"Beta geo{i}" for i in range(12)}
    c={f"Alpha clim{i}" for i in range(12)} | {f"Beta clim{i}" for i in range(12)}
    lens={s:(4.0 if s.startswith("Alpha") else .3) for s in g|c}
    return full,scopes,x,g,c,lens


def test_external_numeric_congeners_predict_actual_island_depth_on_heldout_names():
    full,scopes,depth,geo,clim,lens=example_roster()
    q=source_predict_same_island(
        full_island=full,plant_scopes=scopes,
        source_geo5=geo,source_climate10=clim,
        target_depth=depth,length_cm=lens,cal=calibration()
    )
    for scope in SCOPES:
        p=q[scope]
        assert p["n_same_island_direct_target_species"]==6
        assert p["n_original_target_genera"]==2
        assert p["geo5"]["n_distinct_original_source_numeric_congeners"]==24
        assert p["climate10"]["n_distinct_original_source_numeric_congeners"]==24
        assert p["geo5"]["brier_length_genus"]<p["geo5"]["brier_length_pooled"]
        assert p["geo5"]["n_original_source_genera_with_at_least_two_distinct_numeric_donors"]==2


def test_original_target_name_never_used_as_numeric_source():
    full,scopes,depth,geo,clim,lens=example_roster()
    with pytest.raises(ValueError,match="leaked"):
        source_predict_same_island(
            full_island=full,plant_scopes=scopes,
            source_geo5=geo|{"Alpha isle1"},source_climate10=clim,
            target_depth=depth,length_cm=lens,cal=calibration(),
        )


def example_table(n=40,missing=False,coverage=False):
    full,scopes,depth,geo,clim,lens=example_roster()
    m=calibration()
    out=source_predict_same_island(
        full_island=full,plant_scopes=scopes,source_geo5=geo,
        source_climate10=clim,target_depth=depth,
        length_cm=lens,cal=m,
    )
    records=[]
    for i in range(n):
        for scope in SCOPES:
            for tier in TIERS:
                info=out[scope]
                row=dict(island_id=f"isle{i}",region="tropical",
                         scope=scope,tier=tier,
                         spatial_block=f"block{i//2}",
                         first_source=f"source{(i//2)%6}",
                         n_target_species=info["n_same_island_direct_target_species"],
                         n_target_genera=info["n_original_target_genera"],
                         island_documentation_fraction=1.)
                for mode in ("geo5","climate10"):
                    row[mode+"__n_numeric_donors"]=info[mode]["n_distinct_original_source_numeric_congeners"]
                    row[mode+"__source_numeric_coverage"]=(.1 if coverage else 1.)
                    row[mode+"__brier_genus"]=info[mode]["brier_length_genus"]
                    row[mode+"__brier_pooled"]=info[mode]["brier_length_pooled"]
                records.append(row)
    df=pd.DataFrame(records)
    if missing:
        idx=df.island_id.eq("isle0")
        df.loc[idx,"n_target_species"]=0
        df.loc[idx,"n_target_genera"]=0
        for mode in ("geo5","climate10"):
            df.loc[idx,mode+"__n_numeric_donors"]=0
            df.loc[idx,mode+"__brier_genus"]=np.nan
            df.loc[idx,mode+"__brier_pooled"]=np.nan
    return df


def test_aggregate_donor_models_on_same_original_focal_islands():
    q=aggregate_external_source_numeric_prediction(
        example_table(),population=40,bootstrap=299
    )
    assert q["n_original_island_status_source_tier_rows"]==240
    assert q["n_region_status_tier_panels"]==30
    t=next(x for x in q["panels"] if x["region"]=="tropical"
           and x["scope"]=="all_recorded"
           and x["tier"]==TIERS[0])
    assert t["n_same_island_both_source_mode_qualified"]==40
    assert t["n_original_documentation_balanced_both_10pp"]==40
    assert t["all__spatial_block"]["status"]=="EXPLORATORY_EXTERNAL_NUMERIC_SOURCE_PREDICTION"
    assert t["all__spatial_block"]["geo5__brier_gain__equal_group_mean"]>0


def test_documentation_or_source_data_missing_not_imputed():
    x=aggregate_external_source_numeric_prediction(
        example_table(n=40,coverage=True),population=40,bootstrap=299
    )
    z=next(t for t in x["panels"] if t["region"]=="tropical" and
           t["scope"]=="all_recorded" and t["tier"]==TIERS[0])
    assert z["n_same_island_both_source_mode_qualified"]==40
    assert z["n_original_documentation_balanced_both_10pp"]==0
    assert z["coverage10__spatial_block"]["status"]=="NOT_IDENTIFIABLE_ORIGINAL_SPATIAL_SOURCE_SUPPORT"
    y=aggregate_external_source_numeric_prediction(
        example_table(n=40,missing=True),population=40,bootstrap=299
    )
    z=next(t for t in y["panels"] if t["region"]=="tropical" and
           t["scope"]=="all_recorded" and t["tier"]==TIERS[0])
    assert z["n_same_island_both_source_mode_qualified"]==39
    assert y["n_original_botanical_flora_islands"]==40


def test_duplicate_study_original_island_does_not_pass():
    v=example_table()
    with pytest.raises(ValueError,match="duplicate"):
        aggregate_external_source_numeric_prediction(
            pd.concat([v,v.iloc[:1]]),population=40,bootstrap=299)


def test_exact_genus_out_calibration_never_uses_test_genus_original_depth():
    numerical,deep,mainland=training()
    get_cal,meta=make_genus_heldout_length_calibrator(
        numerical,deep,mainland
    )
    assert meta["n_source_accepted_species_with_NHM_numeric_and_Direct_depth"]==240
    assert meta["every_prediction_excludes_its_own_genus_from_all_original_direct_training"]
    a=get_cal("Source0")
    b=get_cal("Source1")
    assert a["heldout_target_genus_never_in_calibration"]
    assert a["n_calibration_distinct_original_accepted_species"]==236
    assert b["n_calibration_distinct_original_accepted_species"]==236
    assert a["log1p_cm_slope"]>0
    assert calibrated_prob(4.0,a)>calibrated_prob(.3,a)


def test_genus_specific_calibrator_used_for_source_species_without_target_genus_label():
    numerical,deep,mainland=training()
    getter,meta=make_genus_heldout_length_calibrator(
        numerical,deep,mainland
    )
    full,scopes,target,geo,climate,lengths=example_roster()
    out=source_predict_same_island(
        full_island=full,plant_scopes=scopes,
        source_geo5=geo,source_climate10=climate,
        target_depth=target,length_cm=lengths,cal=getter,
    )
    p=out["all_recorded"]
    assert p["n_same_island_direct_target_species"]==6
    assert p["geo5"]["brier_length_genus"]<p["geo5"]["brier_length_pooled"]
    assert getter("Alpha")["n_calibration_distinct_original_accepted_species"]==240


def test_no_false_genus_heldout_calibration_when_corolla_ordinal_overlap_absent():
    with pytest.raises(ValueError,match="insufficient"):
        make_genus_heldout_length_calibrator(
            {"Alpha a":1.},
            {"Alpha a":"deep"},
            {"Alpha a"},
        )
