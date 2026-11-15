"""Falsification tests for exact external-text source-only island flower-colour prediction.

No original source botanical IDs, islands or protected GIFT source names.
"""
import numpy as np
import pandas as pd
import pytest
from fd_nhm_global_source_only_colour import (
    source_probability,proper_flora_losses,
    build_external_source_catalog,placebo_assignment,_bootstrap_equal_group,
)


def test_colour_distribution_is_strictly_proper_and_uniform_in_no_evidence():
    prior=np.array([.2,.3,.5],float)
    assert np.allclose(source_probability(np.zeros(3),prior),prior)
    p=source_probability(np.array([20.,0.,0.]),prior)
    assert p[0]>.2
    assert np.all(p>0)
    assert p.sum()==pytest.approx(1)
    l,b=proper_flora_losses([0,1,2],p)
    assert l>0 and b>=0
    with pytest.raises(ValueError):
        source_probability(np.array([-1.,0,1]),prior)
    with pytest.raises(ValueError):
        proper_flora_losses([3],p)
    with pytest.raises(ValueError):
        proper_flora_losses([],p)


def _example_catalog():
    land=pd.DataFrame([
        dict(entity_ID="R0",work_species="A x"),
        dict(entity_ID="R0",work_species="B y"),
        dict(entity_ID="R1",work_species="C z"),
        dict(entity_ID="R1",work_species="D z"),
        dict(entity_ID="R2",work_species="A u"),
        dict(entity_ID="R2",work_species="B v"),
        dict(entity_ID="R3",work_species="C w"),
        dict(entity_ID="R3",work_species="D w"),
    ])
    labels=dict(zip(land.work_species,["white","red_pink",
                  "white","red_pink","white","red_pink","white","red_pink"]))
    cats={"red_pink":0,"white":1}
    assign=pd.DataFrame([
        dict(island_id=i,source_mode="geo_k5",source_rank=1,entity_ID=e)
        for i,e in (("I0","R0"),("I1","R1"),("I2","R2"),("I3","R3"))])
    geo=pd.DataFrame([
        dict(island_id=i,analysis_regime="tropical",
             spatial_block=f"block{j}")
        for j,i in enumerate(["I0","I1","I2","I3"])])
    return land,labels,cats,assign,geo


def test_source_catalog_counts_actual_taxon_names_only_once_per_near_source():
    land,labels,cats,ass,geo=_example_catalog()
    c=build_external_source_catalog(land,ass,geo,labels,cats)
    assert len(c)==4
    assert all(c[x]["n_source_species"]==2 for x in c)
    assert all(c[x]["n_source_genera"]==2 for x in c)
    assert all(c[x]["values"].sum()==2 for x in c)
    assert all(c[x]["region"]=="tropical" for x in c)
    assert all(c[x]["values"].tolist()==[1.,1.] for x in c)


def test_placebo_choice_is_outcome_blind_and_source_info_matched():
    land,labels,cats,ass,geo=_example_catalog()
    c=build_external_source_catalog(land,ass,geo,labels,cats)
    by={"tropical":sorted(c)}
    choices=placebo_assignment("I0",c,by)
    assert len(choices)==2
    assert len({x["island"] for x in choices})==2
    assert all(x["source_donor_ratio"]==pytest.approx(1) for x in choices)
    assert all(x["source_genus_ratio"]==pytest.approx(1) for x in choices)
    assert all(x["source_entity_jaccard"]==0 for x in choices)
    assert choices==placebo_assignment("I0",c,by)
    assert placebo_assignment("NOTFOUND",c,by)==[]


def test_cluster_equal_is_a_sensitivity_not_many_studies():
    rows=[]
    for i in range(40):
        rows.append(dict(block=f"B{i//2}",
                         first_source=f"S{i//4}",
                         contrast=np.array([.2,-.05])))
    a=_bootstrap_equal_group(rows,"block",seed=7,draws=199)
    b=_bootstrap_equal_group(rows,"first_source",seed=8,draws=199)
    assert a["n_groups"]==20 and b["n_groups"]==10
    assert a["mean_gain_logloss"]==pytest.approx(.2)
    assert b["mean_gain_brier"]==pytest.approx(-.05)
    assert a["blockwise_bootstrap_logloss_95"]==pytest.approx([.2,.2])
    assert b["blockwise_bootstrap_brier_95"]==pytest.approx([-.05,-.05])
    assert _bootstrap_equal_group([],"block") is None
