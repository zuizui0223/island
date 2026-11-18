"""Synthetic-only regressions for matched whole-flora botanical genus beta."""
import numpy as np
import pandas as pd
import pytest
from fd_island_genus_disjoint_matched_beta import (
    select_focal_triplets, analyze_matched,
)

def _metadata():
    r=np.array(["tropical"]*5,dtype=object)
    block=np.array(["b0","b1","b2","b3","b4"],dtype=object)
    # Focal i=0. Donor 1 is different genus; donor 2 shares genus,
    # donor 3 shares species, donor 4 has unknown/uncertified genus.
    shared_s=np.eye(5,dtype=int)
    shared_s[0,3]=shared_s[3,0]=1
    shared_g=np.eye(5,dtype=int)
    shared_g[0,2]=shared_g[2,0]=1
    shared_g[0,3]=shared_g[3,0]=1
    km=np.full((5,5),110.,float)
    np.fill_diagonal(km,0)
    km[0,1]=km[1,0]=100.
    km[0,2]=km[2,0]=115.
    return r,block,shared_s,shared_g,km

def test_match_is_independent_of_flower_outcome_and_excludes_unknown():
    r,b,s,g,km=_metadata()
    kwargs=dict(factor=2.,climate_pc1=np.ones(5),
        n_direct=np.full(5,5.),n_full=np.full(5,10.),
        n_genera=np.full(5,2.),area_km2=np.full(5,100.),
        direct_fraction=np.full(5,.5))
    matches,_=select_focal_triplets(r,b,[True]*4+[False],s,g,km,**kwargs)
    assert any(a==0 and j0==1 and j1==2 for a,j0,j1,_ in matches)
    disjoint,_=select_focal_triplets(r,b,[True]*4+[False],s,g,km,
        donor_blocks_distinct=True,**kwargs)
    assert any(a==0 and j0==1 and j1==2 for a,j0,j1,_ in disjoint)
    far=km.copy();far[0,2]=far[2,0]=9999.
    out,_=select_focal_triplets(r,b,[True]*4+[False],s,g,far,**kwargs)
    assert not any(a==0 and j1==2 for a,j0,j1,_ in out)

def test_overlap_match_preserves_species_categories_and_scope():
    names={
        "I0":([f"A a{k}" for k in range(5)]),
        "I1":([f"B b{k}" for k in range(5)]),
        "I2":([f"A c{k}" for k in range(5)]),
        "I3":([f"C d{k}" for k in range(5)])}
    traits={"flower_primary_color":{x:"white" for group in names.values()
                                    for x in group},
            "floral_form":{x:("tube" if x.startswith("A ") else "open")
                           for group in names.values() for x in group}}
    flora=pd.DataFrame([dict(island_id=i,accepted_species=sp,
        origin_status="native") for i,v in names.items() for sp in v])
    geo=pd.DataFrame([dict(island_id=i,analysis_regime="tropical",
        spatial_block=f"b{k}",island_latitude=0.,
        island_longitude=k,area_km2=100.,climate_pc1=0.,
        weigelt_archipelago=f"ar{k}") for k,i in enumerate(names)])
    out=analyze_matched(traits,flora,geo,boot=99)
    assert out["n_source_accepted_island_floras"]==4
    assert out["eligible_direct_islands"]["all_recorded"]==4
    row=next(x for x in out["results"] if
        x["scope"]=="all_recorded" and x["region"]=="GLOBAL_ALL" and
        x["mode"]=="primary" and x["covariate_factor"]==2.)
    assert row["n_matched_focals"]>=1
    assert row["mean_G0_overlap"]>=0 and row["mean_G1_overlap"]>=0
    assert row["source_gift_first_mainland_region_unavailable"]
    assert out["claims"]["matched_outcome_blind"]

def test_invalid_covariates_do_not_manufacture_matches():
    r,b,s,g,km=_metadata()
    matches,_=select_focal_triplets(r,b,[True]*5,s,g,km,
        np.array([5,5,0,5,5]),np.ones(5)*10.,np.ones(5)*2.,
        np.ones(5)*100.,np.ones(5)*.5)
    assert not any(a==0 and j1==2 for a,j0,j1,_ in matches)
