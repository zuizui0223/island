import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from fd.globi_morphology_falsification import (
    effort_stratified_genus_permutation, genus_heldout_guild_breadth
)


def _permutation_example(seed=44, n=90, g=115):
    rng=np.random.default_rng(seed)
    idx=rng.integers(0,g,size=(n,4))
    W=csr_matrix((np.full(n*4, .25),
                  (np.repeat(np.arange(n),4),idx.ravel())),shape=(n,g))
    v=rng.normal(2.5,.4,g)
    refs=rng.integers(3,150,size=g)
    plant=rng.integers(1,12,size=g)
    N=np.column_stack([np.ones(n),rng.normal(size=(n,2))])
    y=np.asarray(W@v).ravel()*.7+N[:,1]*.4+rng.normal(0,.04,n)
    return W,v,refs,plant,y,N


def test_evidence_binned_genus_permutation_detects_synthetic_association():
    r=effort_stratified_genus_permutation(*_permutation_example(),repeats=399)
    assert r["n_islands"]==90
    assert r["n_globi_genera"]==115
    assert r["observed_partial_slope"] > .3
    assert r["two_sided_posthoc_tail_fraction"] < .05
    assert r["biological_phylogenetic_independence_tested"] is False
    assert not {"genus","island_id","spatial_block"} & set(r)


def test_missing_genus_weights_and_rank_deficiency_are_rejected():
    W,v,r,p,y,N=_permutation_example()
    with pytest.raises(ValueError,match="weights"):
        effort_stratified_genus_permutation(W*.8,v,r,p,y,N,repeats=399)
    with pytest.raises(ValueError,match="Rank"):
        effort_stratified_genus_permutation(W,v,r,p,y,np.column_stack([N,N[:,0]]),repeats=399)


def _genera(n=110):
    rng=np.random.default_rng(3)
    morph=rng.uniform(0,1,n)
    return pd.DataFrame({
        "genus":[f"G{i:03d}" for i in range(n)],
        "effective_channel_number":2+.8*morph+.1*rng.normal(size=n),
        "log_refs":np.log1p(rng.integers(3,500,n)),
        "log_plant":np.log1p(rng.integers(1,10,n)),
        "log_n_morph":np.log1p(rng.integers(1,9,n)),
        "colour_red":morph,
        "form_open":rng.uniform(0,1,n),
        "joint_color_form":morph*.7,
    })


def test_one_genus_one_fold_and_source_effort_benchmark():
    frame=_genera()
    out=genus_heldout_guild_breadth(
        frame,
        base_features=["log_refs","log_plant","log_n_morph"],
        marginal_features=["colour_red","form_open"],
        pairing_features=["joint_color_form"],
    )
    assert len(out)==3
    assert set(out.n_heldout_genera)=={110}
    assert out.set_index("model").loc["morph_marginal","RMSE"] < out.set_index("model").loc["evidence_only","RMSE"]
    assert not {"genus","island_id"} & set(out.columns)


def test_duplicate_genus_and_missing_morph_are_not_silently_filled():
    a=_genera()
    b=a.copy()
    b.loc[1,"genus"]=b.loc[0,"genus"]
    with pytest.raises(ValueError,match="one row"):
        genus_heldout_guild_breadth(
            b,base_features=["log_refs"],marginal_features=["colour_red"],
            pairing_features=["joint_color_form"])
    b=a.copy();b.loc[1,"colour_red"]=float("nan")
    with pytest.raises(ValueError,match="Unknown"):
        genus_heldout_guild_breadth(
            b,base_features=["log_refs"],marginal_features=["colour_red"],
            pairing_features=["joint_color_form"])
