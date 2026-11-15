"""Synthetic pre-outcome category and geographic/source coverage safeguards."""
import pandas as pd
import pytest
from fd_directional_flower_space_preflight import preflight, CHECK_TRAITS, PAIRS


def synthetic():
    states={t:{} for t in CHECK_TRAITS}
    rows=[]
    geo=[]
    for i in range(12):
        island=f"Island{i}"
        for k in range(5):
            sp=f"Gen{k} sp{i}_{k}"
            rows.append(dict(island_id=island,accepted_species=sp,
                             origin_status="native" if k<4 else "introduced"))
            for t in CHECK_TRAITS:
                states[t][sp]={
                    "floral_form":("tubular" if k<2 else "open_radial"),
                    "floral_symmetry":("zygomorphic" if k<2 else "actinomorphic"),
                    "tube_depth_class":("deep" if k<2 else "shallow"),
                    "flower_size_class":"small",
                    "autonomous_selfing_capacity":"autonomous",
                    "self_incompatibility":"si",
                    "flower_primary_color":"white"}[t]
        geo.append(dict(island_id=island,
             analysis_regime="tropical" if i<6 else "northern_midlatitude",
             spatial_block=f"B{i}",
             island_latitude=65. if i==11 else 20.,
             distance_to_continent_km=100+10*i,
             area_km2=5+i,
             climate_pc1=.1,climate_pc2=.2,climate_pc3=.3,climate_pc4=.4))
    return states,pd.DataFrame(rows),pd.DataFrame(geo)


def test_original_single_trait_and_true_same_plant_pair_support_not_triple_gate():
    s,f,g=synthetic()
    q=preflight(s,f,g)
    assert q["n_original_source_botanical_islands"]==12
    assert q["n_original_islands_with_finite_verified_distance_and_area"]==12
    assert q["original_direct_trait_categories"]["floral_form"]["n_original_unique_Direct_species"]==60
    assert q["original_actual_source_two_trait_combinations"][
        "floral_form__floral_symmetry"]["n_exact_same_original_Direct_species"]==60
    z=next(r for r in q["directional_isolation_trait_and_region_eligibility"]
           if r["scope"]=="all_recorded" and r["endpoint"]=="floral_form"
           and r["region"]=="tropical")
    assert z["n_islands_direct_n5"]==6
    assert z["n_n5_islands_with_verified_isolation_and_geography"]==6
    assert z["n_spatial_blocks_with_isolation_n5"]==6
    strict=next(r for r in q["directional_isolation_trait_and_region_eligibility"]
           if r["scope"]=="original_native" and r["endpoint"]=="floral_form"
           and r["region"]=="tropical")
    assert strict["n_islands_direct_n5"]==0
    nested=next(r for r in q["directional_isolation_trait_and_region_eligibility"]
           if r["scope"]=="all_recorded" and r["endpoint"]=="floral_form"
           and r["region"]=="highlatitude_ge60_nested")
    assert nested["n_islands_direct_n5"]==1
    assert q["claims"]["missing_original_trait_categories_never_imputed_zero"]


def test_unresolved_isolation_retains_original_population_not_false_zero():
    s,f,g=synthetic()
    g.loc[g.island_id.eq("Island1"),"distance_to_continent_km"]=-3.
    a=preflight(s,f,g)
    assert a["n_original_source_botanical_islands"]==12
    assert a["n_original_geographic_isolation_unresolved"]==1
    n=next(x for x in a["directional_isolation_trait_and_region_eligibility"]
           if x["scope"]=="all_recorded" and x["endpoint"]=="floral_form"
           and x["region"]=="tropical")
    assert n["n_islands_direct_n5"]==6
    assert n["n_n5_islands_with_verified_isolation_and_geography"]==5


def test_missing_joint_trait_never_completes_by_cross_species_product():
    s,f,g=synthetic()
    for sp in list(s["floral_symmetry"]):
        if sp.endswith("_0") or sp.endswith("_1") or sp.endswith("_2"):
            del s["floral_symmetry"][sp]
    a=preflight(s,f,g)
    q=a["original_actual_source_two_trait_combinations"][
        "floral_form__floral_symmetry"]
    assert q["n_exact_same_original_Direct_species"]==24
    n=next(x for x in a["directional_isolation_trait_and_region_eligibility"]
           if x["scope"]=="all_recorded" and
           x["endpoint"]=="floral_form__floral_symmetry" and
           x["region"]=="tropical")
    assert n["n_islands_direct_n5"]==0


def test_duplicate_source_or_wrong_geography_rejected():
    s,f,g=synthetic()
    with pytest.raises(ValueError):
        preflight(s,pd.concat([f,f.iloc[:1]]),g)
    with pytest.raises(ValueError):
        preflight(s,f,g.drop(columns=["distance_to_continent_km"]))
