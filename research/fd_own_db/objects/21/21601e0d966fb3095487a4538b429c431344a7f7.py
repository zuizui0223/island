import pandas as pd
import pytest

from fd.geology_lgm_joint_identifiability import (
    audit_geology_lgm_common_support,CELLS,
)
from fd.real_joint_syndrome_access import JOINT_SPECS


def fake_original(n=160):
    assert n%4==0
    rows=[]
    for i in range(n):
        g,l=CELLS[i//(n//4)]
        ocean=(g=="oceanic")
        unconn=(l=="lgm_unconnected")
        effect=.10 + .03*ocean + .02*unconn + .04*(ocean and unconn)
        between=.025 + .02*ocean + .01*unconn
        within=.005 + .01*(ocean and unconn)
        for scope in ("all_recorded","strict_native","strict_native_nonendemic"):
            for sector in JOINT_SPECS:
                a=dict(
                    island_id=f"I{i}",regime="tropical",
                    scope=scope,sector=sector,
                    spatial_block=f"B{i//2}",first_source=f"S{i//2}",
                    geological_origin=g,LGM_connectivity=l,
                    island_total=15,island__n_direct=8,
                    island__n_variable_genera=1,
                    n_shared_informative_genera_both_source_modes=1,
                )
                for mode in ("geo5","climate10"):
                    a[f"{mode}__n_direct"]=20
                    a[f"{mode}__total"]=40
                    a[f"{mode}__delta_joint"]=effect
                    a[f"{mode}__delta_coupling"]=between+within
                    a[f"{mode}__delta_between_genus_coupling"]=between
                    a[f"{mode}__delta_within_genus_coupling"]=within
                rows.append(a)
    return pd.DataFrame(rows)


def getpanel(report,scope="all_recorded",region="tropical",sector="actinomorphic_deep"):
    return next(p for p in report["panels"] if
        p["scope"]==scope and p["region"]==region and p["sector"]==sector)


def test_4505_population_contract_four_cells_and_source_support():
    f=fake_original()
    z=audit_geology_lgm_common_support(f,population=160,draws=199)
    assert z["n_original_global_botanical_floras"]==160
    assert z["n_original_full_flora_scope_pair_rows"]==160*3*11
    assert z["n_reporting_region_scope_flower_panels"]==3*3*3
    p=getpanel(z)
    assert p["n_original_pair_source_comparable_before_two_history_labels"]==160
    assert p["n_both_geological_and_LGM_history_qualified"]==160
    for g,l in CELLS:
        v=p["original_four_joint_history_cell_support"][f"{g}__{l}"]
        assert v["n_islands"]==40
        assert v["n_geographic_blocks"]==20
        assert v["n_first_mainland_sources"]==20
        assert v["n_same_informative_genus_island_and_both_sources"]==40
    assert p["simultaneous_origin_and_LGM_equal_strata_adjustment"]["n_four_eligible_cells"]==4
    assert z["n_both_history_axes_jointly_supported"]==2*3*3
    assert z["n_history_stratified_panels_supported"]==2*3*3*4


def test_four_cell_partial_adjustment_not_confounded_by_other_history():
    z=audit_geology_lgm_common_support(fake_original(),population=160,draws=199)
    p=getpanel(z)
    for unit in ("spatial_block","first_source"):
        q=p["simultaneous_origin_and_LGM_equal_strata_adjustment"][unit]
        assert q["status"]=="DESCRIPTIVE_FULL_FOUR_CELL_GEOL_LGM_OVERLAP"
        assert q["oceanic_minus_continental_equal_LGM_strata__geo5__delta_joint"]==pytest.approx(.05)
        assert q["LGM_unconnected_minus_connected_equal_geology_strata__geo5__delta_joint"]==pytest.approx(.04)
        assert q["geology_by_LGM_interaction__geo5__delta_joint"]==pytest.approx(.04)
        a=p["within_one_history_class_contrasts"]["geological_origin__within__LGM_connectivity__lgm_unconnected"]
        assert a[unit]["geo5__delta_joint__category_b_minus_a"]==pytest.approx(.07)


def test_missing_one_geology_by_history_cell_does_not_estimate_adjusted_effect():
    f=fake_original()
    missing=f.geological_origin.eq("oceanic")&f.LGM_connectivity.eq("lgm_connected")
    f.loc[missing,"geo5__n_direct"]=0
    q=audit_geology_lgm_common_support(f,population=160,draws=199)
    p=getpanel(q)
    z=p["simultaneous_origin_and_LGM_equal_strata_adjustment"]
    assert z["status"]=="NOT_IDENTIFIABLE_ALL_FOUR_HISTORY_CELLS"
    assert z["n_four_eligible_cells"]==3
    assert "spatial_block" not in z
    c=p["within_one_history_class_contrasts"]
    assert c["geological_origin__within__LGM_connectivity__lgm_connected"]["status"].startswith("NOT_IDENTIFIABLE")
    assert c["geological_origin__within__LGM_connectivity__lgm_unconnected"]["status"]=="DESCRIPTIVE_HISTORY_CONDITIONAL_BOTANICAL_SOURCE_FLORA"


def test_independent_origin_and_GMMC_missingness_stays_unresolved_and_not_zero():
    f=fake_original()
    a=f.island_id.eq("I0")
    f.loc[a,"geological_origin"]="unresolved"
    b=f.island_id.eq("I80")
    f.loc[b,"LGM_connectivity"]="unresolved"
    z=audit_geology_lgm_common_support(f,population=160,draws=199)
    p=getpanel(z)
    assert p["n_both_geological_and_LGM_history_qualified"]==158
    assert p["n_geologically_classified_but_LGM_unknown"]==1
    assert p["n_LGM_classified_but_geology_unknown"]==1
    assert z["limits"]["geological_origin_Roeble_Zizka_independent_not_distance_proxy"]
    assert z["limits"]["original_native_status_nonrandomly_documented"]


def test_complete_original_parent_and_untrusted_geo_input_fail_closed():
    f=fake_original()
    with pytest.raises(ValueError,match="4505-island|complete"):
        audit_geology_lgm_common_support(f.iloc[:-1],population=160,draws=199)
    f=fake_original()
    f.loc[0,"geological_origin"]="distance_imputed_oceanic"
    with pytest.raises(ValueError,match="geology"):
        audit_geology_lgm_common_support(f,population=160,draws=199)
