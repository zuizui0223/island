"""Synthetic, source-free tests of genuinely same accepted-species flower pairs."""
import numpy as np
import pandas as pd
import pytest
from fd.real_joint_syndrome_access import (
    JOINT_SPECS, SCOPES, make_joint_species, summarize_one_roster,
    finite_joint_null_counts, same_island_pair_source_frame,
)


def synthetic_states():
    return {
        "floral_symmetry":{
            "plant1":"actinomorphic", "plant2":"zygomorphic",
            "plant3":"actinomorphic", "plant4":"actinomorphic",
        },
        "tube_depth_class":{
            "plant2":"deep", "plant3":"deep", "plant4":"shallow",
            "plant5":"deep",
        },
        "floral_form":{
            "plant2":"tubular", "plant3":"bell_campanulate",
            "plant5":"open_radial",
        },
        "flower_primary_color":{
            "plant1":"red_pink", "plant2":"white", "plant3":"red_pink",
            "plant4":"yellow_orange", "plant5":"white",
        },
    }


def test_real_pair_never_multiplies_marginal_species():
    m=make_joint_species(synthetic_states())
    radial=m["actinomorphic_deep"]
    assert radial["known"]==frozenset({"plant2","plant3","plant4"})
    assert radial["joint"]==frozenset({"plant3"})
    s=summarize_one_roster({"plant1","plant2","plant3","plant4","plant5"},radial)
    assert s["n"]==3 and s["n_A"]==2 and s["n_B"]==2 and s["n_joint"]==1
    assert np.isclose(s["p_joint"],1/3)
    assert np.isclose(s["coupling"],1/3-4/9)
    assert "plant1" not in radial["known"]


def test_colour_and_tube_are_same_species_and_not_effective_guild():
    x=make_joint_species(synthetic_states())
    assert x["red_pink_deep"]["joint"]==frozenset({"plant3"})
    assert x["white_deep"]["joint"]==frozenset({"plant2","plant5"})
    assert x["red_pink_tubular_bell"]["joint"]==frozenset({"plant3"})
    assert x["white_tubular_bell"]["joint"]==frozenset({"plant2"})
    assert len(x)==len(JOINT_SPECS)==11
    assert x["yellow_orange_deep"]["joint"]==frozenset()
    assert x["yellow_orange_tubular_bell"]["joint"]==frozenset()


def test_missing_direct_pair_means_unknown_not_zero():
    m=make_joint_species(synthetic_states())
    z=summarize_one_roster({"plant1"},m["actinomorphic_deep"])
    assert z["n"]==0
    assert np.isnan(z["p_joint"]) and np.isnan(z["coupling"])
    assert z["n_joint"]==0


def test_exact_marginal_null_has_finite_hypergeometric_variance():
    exp,var=finite_joint_null_counts(10,4,5)
    assert np.isclose(exp,2.)
    assert np.isclose(var, 2/3)
    assert finite_joint_null_counts(1,1,1)==(1.,0.)
    with pytest.raises(ValueError):
        finite_joint_null_counts(3,4,2)


def sample(n=40):
    rows=[]
    for i in range(n):
        for scope in SCOPES:
            for sector in JOINT_SPECS:
                rows.append(dict(
                    island_id=f"synthetic_island_{i}",regime="tropical",
                    scope=scope,sector=sector,spatial_block=f"block{i//2}",
                    first_source=f"source{(i//2)%6}",
                    island_total=20,island_direct=10,
                    geo5_total=30,geo5_direct=15,
                    climate10_total=30,climate10_direct=15,
                    island_pair_prevalence=.30,
                    geo5_pair_prevalence=.20,
                    climate10_pair_prevalence=.25,
                    island_pair_coupling=.06,
                    geo5_pair_coupling=.02,
                    climate10_pair_coupling=.01,
                ))
    return pd.DataFrame(rows)


def test_same_island_both_source_modes_separate_pair_prevalence_and_coupling():
    d=same_island_pair_source_frame(sample(),population=40,bootstrap=299)
    assert d["n_real_two_trait_pairs"]==11
    assert d["n_source_mode_provenance_region_panels"]==5*3*11
    z=next(a for a in d["panels"] if a["region"]=="tropical"
           and a["plant_scope"]=="all_recorded"
           and a["documented_joint_sector"]=="actinomorphic_deep")
    assert z["n_both_source_modes_and_joint_carriers"]==40
    assert z["n_original_documentation_balanced_le10pp"]==40
    t=z["all__spatial_block"]
    assert t["status"]=="DESCRIPTIVE_REAL_SAME_SPECIES_PAIR_ONLY"
    assert np.isclose(t["delta_geo5_pair_prevalence_block_equal"],.10)
    assert np.isclose(t["delta_climate10_pair_prevalence_block_equal"],.05)
    assert np.isclose(t["delta_geo5_pair_coupling_block_equal"],.04)
    assert np.isclose(t["delta_climate10_pair_coupling_block_equal"],.05)


def test_missing_source_is_not_zero_and_does_not_remove_original_island():
    x=sample()
    ix=x.island_id.eq("synthetic_island_0")
    x.loc[ix,"climate10_total"]=0
    x.loc[ix,"climate10_direct"]=0
    x.loc[ix,"climate10_pair_prevalence"]=np.nan
    x.loc[ix,"climate10_pair_coupling"]=np.nan
    z=same_island_pair_source_frame(x,population=40,bootstrap=299)
    q=next(a for a in z["panels"] if a["region"]=="tropical"
           and a["plant_scope"]=="all_recorded"
           and a["documented_joint_sector"]=="actinomorphic_deep")
    assert z["n_full_original_botanical_islands"]==40
    assert q["n_both_source_modes_and_joint_carriers"]==39


def test_duplicate_focal_island_and_false_pair_missingness_rejected():
    x=sample()
    with pytest.raises(ValueError,match="duplicate"):
        same_island_pair_source_frame(pd.concat([x,x.iloc[:1]]),
                                      population=40,bootstrap=299)
    x=sample()
    x.loc[0,"island_pair_coupling"]=np.nan
    with pytest.raises(ValueError,match="missingness"):
        same_island_pair_source_frame(x,population=40,bootstrap=299)
