"""Synthetic fixed-island/fixed-placebo cross-corpus checks; no botanical IDs."""
import numpy as np
import pandas as pd
import pytest

from fd_nhm_provider_stratified_colour_placebo import (
    independent_donors,matched_controls,aggregate_comparison)
from fd_nhm_global_source_only_colour import NHM_ALLOWED


def _catalog():
    groups={}
    for provider in ("BHL_text","FloraChina_text"):
        d={}
        for name,block,e,n,g in (
            ("focal","A",("e0","e1"),100,20),
            ("p1","B",("e2","e3"),110,19),
            ("p2","C",("e4","e5"),95,21),
            ("p3","A",("e6","e7"),100,20),
            ("p4","D",("e0","e1"),90,21),
            ("p5","E",("e8","e9"),1000,100),
        ):
            if provider=="FloraChina_text":
                n=int(n*.6)
                g=int(g*.8)
            d[name]=dict(region="northern_midlatitude",block=block,
                         first_source=e[0],source_entities=frozenset(e),
                         n_source_species=n,n_source_genera=g,
                         values=np.array([10,8,6,3,1.]))
        groups[provider]=d
    return groups


def test_common_placebo_choices_are_provider_blind():
    cats=_catalog()
    by={"northern_midlatitude":sorted(cats["BHL_text"])}
    chosen,n=matched_controls("focal",cats,by)
    assert n==2
    assert {p["island"] for p in chosen}=={"p1","p2"}
    assert all(len(p["species_genus_ratios"])==2 for p in chosen)
    assert all(p["source_jaccard"]==0 for p in chosen)
    first=matched_controls("focal",cats,by)
    assert first==(chosen,n)
    assert matched_controls("not-present",cats,by)==([],0)


def test_matching_requirements_hold_in_both_external_sources():
    cats=_catalog()
    cats["FloraChina_text"]["p1"]["n_source_species"]=3
    by={"northern_midlatitude":sorted(cats["BHL_text"])}
    chosen,n=matched_controls("focal",cats,by)
    assert chosen==[]
    assert n==1 # Only p2 still passes


def test_separately_published_labels_never_resurrect_island_species():
    raw=pd.DataFrame([
        dict(taxon="Alpha one",**{"flower colour":"yellow",
             "corolla colour":"","petal colour":"","flower shape":"",
             "flower symmetry":"","source":"BHL","wfo":""}),
        dict(taxon="Beta two",**{"flower colour":"white",
             "corolla colour":"","petal colour":"","flower shape":"",
             "flower symmetry":"","source":"BHL","wfo":""}),
        dict(taxon="Beta two",**{"flower colour":"red",
             "corolla colour":"","petal colour":"","flower shape":"",
             "flower symmetry":"","source":"BHL","wfo":""}),
        dict(taxon="Gamma three",**{"flower colour":"purple",
             "corolla colour":"","petal colour":"","flower shape":"",
             "flower symmetry":"","source":"BHL","wfo":""})
    ])
    donor,meta=independent_donors(raw,
        {"Alpha one","Beta two","Gamma three"},{"Alpha one"})
    assert donor=={"Gamma three":"blue_purple"}
    assert meta["n_species_with_documentation_conflict_by_axis"]["primary"]==1


def test_source_cluster_and_geographic_block_interval_reporting():
    rows=[]
    for i in range(40):
        rows.append(dict(
            block="b"+str(i//2), first_source="s"+str(i//4),
            n_direct=5,source_jaccard=.2,
            contrast_BHL=np.array([.01,.003]),
            contrast_China=np.array([-.005,.001]),
            source_donor_species_ratio_BHL=1.,
            source_donor_species_ratio_China=1.,
            source_donor_genus_ratio_BHL=1.,
            source_donor_genus_ratio_China=1.))
    v=aggregate_comparison(rows,bootstrap=199)
    assert v["status"]=="POSTHOC_SAME_FOCAL_EXTERNAL_NHM_CORPUS_SOURCE_PLACEBO"
    assert v["n_matched_islands"]==40 and v["n_distinct_geographic_blocks"]==20
    assert v["n_distinct_first_GIFT_mainland_source_groups"]==10
    assert v["mean_island_BHL_logloss_gain"]==pytest.approx(.01)
    assert v["mean_island_China_logloss_gain"]==pytest.approx(-.005)
    assert v["mean_island_BHL_minus_China_logloss_advantage"]==pytest.approx(.015)
    assert v["geographic_block_BHL_conditional_logloss_95"]==pytest.approx([.01,.01])
    assert v["first_GIFT_source_China_conditional_logloss_95"]==pytest.approx([-.005,-.005])
    assert v["fraction_focal_floras_both_corpora_favor_real_source"]==0


def test_low_overlap_abstains_not_pseudo_independent_floras():
    rows=[dict(block="b1",first_source="s1",n_direct=5,
               source_jaccard=.2,
               contrast_BHL=np.array([.1,0.]),
               contrast_China=np.array([.1,0.]),
               source_donor_species_ratio_BHL=1.,
               source_donor_species_ratio_China=1.,
               source_donor_genus_ratio_BHL=1.,
               source_donor_genus_ratio_China=1.) for _ in range(4)]
    q=aggregate_comparison(rows)
    assert q["status"]=="INSUFFICIENT_COMMON_EXTERNAL_PROVIDER_SOURCE_SUPPORT"
    assert q["n_matched_islands"]==4
