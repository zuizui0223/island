import numpy as np
import pandas as pd
import pytest

from fd.geology_genus_flower_pair import (
    original_geology_join,original_pair_by_geology,
    FOCAL_PAIRS,AXES,
)
from fd.real_joint_syndrome_access import JOINT_SPECS


def frame(population=120):
    rows=[]
    for i in range(population):
        for scope in ("all_recorded","strict_native",
                      "strict_native_nonendemic"):
            for pair in JOINT_SPECS:
                a=i>=population//2
                x=dict(island_id=f"I{i}",regime="northern_midlatitude",
                    scope=scope,sector=pair,
                    spatial_block=f"B{i//2}",
                    first_source=f"S{i//2}",
                    island_total=12,island__n_direct=6,
                    island__n_variable_genera=1)
                for mode in ("geo5","climate10"):
                    x[f"{mode}__total"]=40
                    x[f"{mode}__n_direct"]=17
                    x[f"{mode}__delta_joint"]=.20 if a else .10
                    x[f"{mode}__delta_coupling"]=.09 if a else .07
                    x[f"{mode}__delta_between_genus_coupling"]=.06 if a else .04
                    x[f"{mode}__delta_within_genus_coupling"]=.03
                rows.append(x)
    return pd.DataFrame(rows)


def geology(population=120):
    return pd.DataFrame([
        dict(island_id=f"I{i}",geological_origin=(
            "oceanic" if i>=population//2 else "continental"),
             LGM_connectivity=(
                 "lgm_connected" if i%2 else "lgm_unconnected"),
             area_km2=20+i,
             distance_to_continent_km=float(i))
        for i in range(population)
    ])


def test_global_floras_and_independent_geology_axis_never_narrow_population():
    main=frame()
    geo=geology()
    joined=original_geology_join(main,geo,population=120)
    assert joined.island_id.nunique()==120
    assert len(joined)==120*3*11
    assert joined.geological_origin.isin(("oceanic","continental")).all()
    assert joined.LGM_connectivity.isin(
        ("lgm_connected","lgm_unconnected")).all()
    q=geo.iloc[:-1]
    incomplete=original_geology_join(main,q,population=120)
    assert len(incomplete)==len(main)
    assert incomplete.loc[incomplete.island_id.eq("I119"),
           "geological_origin"].eq("unresolved").all()


def test_origin_class_same_trait_source_modes_cluster_contrast_exact_partition():
    f=original_geology_join(frame(),geology(),population=120)
    z=original_pair_by_geology(f,population=120,draws=199)
    assert z["original_global_botanical_island_floras"]==120
    assert z["n_original_island_status_true_pair_rows"]==120*3*11
    assert z["n_report_panels"]==2*3*3*3
    assert z["n_independent_GIFT_islands_with_oceanic_or_continental_geology"]==120
    q=next(p for p in z["panels"] if p["axis"]=="geological_origin"
           and p["region"]=="all" and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert q["status"]=="DESCRIPTIVE_ORIGIN_CLASS_BOTH_SOURCE_MODES_SUPPORTED"
    assert q["class_support"]["continental"]["n_islands"]==60
    assert q["class_support"]["oceanic"]["n_islands"]==60
    for cluster in ("spatial_block","first_source"):
        d=q[cluster]
        assert d["status"]=="DESCRIPTIVE_GIFT_ORIGIN_SOURCE_RELATIVE_FLORAL_PAIRING"
        assert d["geo5__delta_joint__category_b_minus_a"]==pytest.approx(.1)
        assert d["geo5__delta_coupling__category_b_minus_a"]==pytest.approx(.02)
        assert d["geo5__delta_between_genus_coupling__category_b_minus_a"]==pytest.approx(.02)
        assert d["geo5__delta_within_genus_coupling__category_b_minus_a"]==pytest.approx(0)
        assert d["geo5__delta_joint__bootstrap95"]==pytest.approx([.1,.1])


def test_geology_and_last_glacial_connectivity_are_distinct_analytical_axes():
    f=original_geology_join(frame(),geology(),population=120)
    z=original_pair_by_geology(f,population=120,draws=199)
    o=next(p for p in z["panels"] if p["axis"]=="LGM_connectivity"
           and p["region"]=="all" and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert o["status"]=="DESCRIPTIVE_ORIGIN_CLASS_BOTH_SOURCE_MODES_SUPPORTED"
    assert o["spatial_block"]["geo5__delta_joint__category_b_minus_a"]==pytest.approx(0)
    assert z["identification_limits"][
        "LGM_connectivity_never_substituted_for_geological_origin"]


def test_source_eligibility_insufficient_is_NOT_zero_geological_origin_effect():
    f=original_geology_join(frame(),geology(),population=120)
    f.loc[f.geological_origin.eq("oceanic"),"geo5__n_direct"]=0
    q=original_pair_by_geology(f,population=120,draws=199)
    p=next(p for p in q["panels"] if p["axis"]=="geological_origin"
           and p["region"]=="all" and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert p["class_support"]["oceanic"]["n_islands"]==0
    assert p["status"]=="NOT_IDENTIFIABLE_ORIGIN_CLASS_SOURCE_SUPPORT"
    assert "geo5__delta_joint__category_b_minus_a" not in p["spatial_block"]


def test_unresolved_and_mixed_cannot_become_oceanic():
    geo=geology()
    geo.loc[1,"geological_origin"]="mixed"
    geo.loc[2,"geological_origin"]="unresolved"
    data=original_geology_join(frame(),geo,population=120)
    assert data.loc[data.island_id.eq("I1"),"geological_origin"].eq("mixed").all()
    assert data.loc[data.island_id.eq("I2"),"geological_origin"].eq("unresolved").all()
    z=original_pair_by_geology(data,population=120,draws=199)
    p=next(p for p in z["panels"] if p["axis"]=="geological_origin"
           and p["region"]=="all" and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert p["n_independently_classified_island_floras"]==118


def test_source_frame_duplicates_and_origin_id_conflicts_fail_closed():
    main=frame()
    geo=geology()
    q=pd.concat([geo,geo.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="duplicate"):
        original_geology_join(main,q,population=120)
    joined=original_geology_join(main,geo,population=120)
    doubled=pd.concat([joined,joined.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="global flora"):
        original_pair_by_geology(doubled,population=120,draws=199)
