"""Synthetic privacy and outcome-blind source-placebo diagnostics."""
import numpy as np
import pandas as pd
import pytest

from fd_source_name_holdout_placebo import (
    source_catalog,coverage,source_placebos,predictor_losses,summarize
)


def _catalog():
    # Each source assignment is an independent five-entity native GIFT set,
    # but there are repeat source names; never count as distinct donors.
    source={
        "F":dict(block="B0",source_entities=frozenset(["e0","e1"]),
                 names={"A a","A b","B c"},n_names=3,
                 genera={"A","B"}),
        "P1":dict(block="B1",source_entities=frozenset(["e2","e3"]),
                 names={"A x","A y","B z"},n_names=3,
                 genera={"A","B"}),
        "P2":dict(block="B2",source_entities=frozenset(["e4","e5"]),
                 names={"A u","A v","B w"},n_names=3,
                 genera={"A","B"}),
        "P3":dict(block="B0",source_entities=frozenset(["e6","e7"]),
                 names={"A s","B t"},n_names=2,
                 genera={"A","B"}),
        "P4":dict(block="B4",source_entities=frozenset(["e0","e1"]),
                 names={"A a","B c"},n_names=2,
                 genera={"A","B"}),
    }
    return source


def test_outcome_blind_marginal_coverage_source_matching():
    from collections import Counter
    c=_catalog()
    # The original focal genus frequencies are a pre-outcome fixed oracle.
    target=Counter({"A":3,"B":1})
    assert coverage(target,c["F"]["genera"])==1
    a,n=source_placebos("F",dict(block="B0"),target,c,
                        tolerance=.02,max_ratio=2)
    assert n==2
    assert {x["identity"] for x in a}=={"P1","P2"}
    assert all(x["genus_coverage"]==1 for x in a)
    assert all(x["donor_ratio"]==1 for x in a)
    assert all(x["source_overlap"]==0 for x in a)
    # Changing any recorded flower category is irrelevant to the donor selection.
    b,m=source_placebos("F",dict(block="B0"),target,c,
                        tolerance=.02,max_ratio=2)
    assert a==b and n==m
    with pytest.raises(ValueError):
        source_placebos("UNKNOWN",dict(block="X"),target,c,tolerance=.02)


def test_missing_placebo_controls_abstains_not_promote_positive_association():
    from collections import Counter
    c=_catalog()
    del c["P2"]
    a,n=source_placebos("F",dict(block="B0"),
                        Counter({"A":1}),c,tolerance=.02)
    assert n==1 and a==[]


def test_global_source_probability_and_phenotype_scores_are_finite():
    from collections import defaultdict,Counter
    donor={"A a":0,"A b":0,"B c":1}
    model=defaultdict(Counter)
    for s,k in donor.items():model[s.split()[0]][k]+=1
    prior=np.array([.5,.5])
    pair=predictor_losses(
        ["A target1","A target2"],donor,prior,model,
        {"A a","B c"},{"A target1":0,"A target2":0,**donor})
    assert pair.shape==(2,)
    assert np.isfinite(pair).all()


def test_source_catalog_does_not_invent_or_double_count_species():
    source=pd.DataFrame([
        dict(entity_ID="X",work_species="A a"),
        dict(entity_ID="Y",work_species="A a"),
        dict(entity_ID="Y",work_species="B b")])
    near={"I0":("X","Y"),"I1":("Y",)}
    geo=pd.DataFrame([
        dict(island_id="I0",analysis_regime="northern_midlatitude",
             spatial_block="q0"),
        dict(island_id="I1",analysis_regime="tropical",spatial_block="q1")])
    catalog=source_catalog(source,near,geo)
    assert set(catalog)=={"I0"}
    assert catalog["I0"]["n_names"]==2
    assert catalog["I0"]["genera"]=={"A","B"}


def test_equal_block_resampling_not_215_independent_lineages():
    data=[]
    for i in range(40):
        data.append(dict(
            block=f"B{i//2}",source_cluster=f"S{i//4}",n_carriers=10,
            actual_logloss=1.0,actual_brier=.3,
            placebo_logloss=1.1,placebo_brier=.35,
            coverage_difference=.01,
            donor_ratio=.95,source_overlap=.10))
    q=summarize(data,seed=3,n_boot=199)
    assert q["n_source_matched_islands"]==40
    assert q["n_physical_blocks"]==20
    assert q["n_distinct_first_rank_GIFT_source_clusters"]==10
    assert q["first_GIFT_source_bootstrap_logloss_interval_95"]==pytest.approx([.1,.1])
    assert q["fraction_islands_actual_source_better_than_matched_placebos"]==1
    assert q["mean_actual_minus_placebo_improvement_logloss"]==pytest.approx(.1)
    assert q["block_bootstrap_logloss_interval_95"]==pytest.approx([.1,.1])
    assert summarize(data[:5])["status"].startswith("insufficient_")


def test_no_empty_target_and_outcome_support_cannot_be_scored():
    with pytest.raises(ValueError):coverage({},set())
    with pytest.raises(ValueError):predictor_losses(
        ["A a"],{},np.array([1.,0.]),{},set(),{"A a":0})
