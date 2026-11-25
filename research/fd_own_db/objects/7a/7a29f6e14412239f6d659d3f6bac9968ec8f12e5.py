"""Source-verified phylogenetic donor matching cannot invent branch lengths."""
from itertools import product
import numpy as np
import pandas as pd
import pytest

from fd.discrete_functional_volume import TRAITS
from fd.rare_profile_decomposition import RarityReference
from fd.family_conditioned_rare_volume import _compiled
from fd.true_backbone_phylo_source_match import _sample_original_phylogenetic_source


def original_like():
    profiles=list(product(
        ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
        ("open_radial","tubular","brush","composite","zygomorphic"),
        ("actinomorphic","zygomorphic")
    ))
    counts=np.arange(1,101,dtype=float)**1.2
    counts=np.floor(counts/counts.sum()*(2234-100)).astype(int)+1
    for j in range(2234-int(counts.sum())):counts[j%100]+=1
    names=[]
    rows=[]
    for i,(p,n) in enumerate(zip(profiles,counts)):
        for j in range(int(n)):
            names.append(f"G{i%40} species_{i}_{j}")
            rows.append(p)
    assert len(names)==2234
    return pd.DataFrame(rows,index=names,columns=TRAITS)


def test_optimal_real_phylogeny_distance_le_random_geo_donor_and_no_plant_reuse():
    t=original_like()
    ref=RarityReference.from_original(t,null_draws=199,seed=2)
    lookup=_compiled(ref,t)
    names=list(t.index)
    original_island=set(names[0:30])
    original_source=set(names[30:80])
    tip={name:i for i,name in enumerate(names[:100])}
    # A valid symmetric non-negative original tip metric for matching; the
    # real source tree's LCA patristic matrix is tested separately.
    q=np.arange(100)
    D=np.abs(q[:,None]-q[None,:]).astype(float)
    a=_sample_original_phylogenetic_source(
        t,ref,lookup,original_island,original_source,tip,D,
        n=15,repeats=20,seed=9)
    b=_sample_original_phylogenetic_source(
        t,ref,lookup,original_island,original_source,tip,D,
        n=15,repeats=20,seed=9)
    assert a==b
    assert a["d_matched"]<=a["d_random"]+1e-10
    assert a["matched_minus_random_phylo_distance"]<=0
    assert np.isfinite(a["matched_delta_k_adjusted_rare_minus_common"])
    assert np.isfinite(a["random_delta_k_adjusted_rare_minus_common"])
    with pytest.raises(ValueError,match="also observed"):
        _sample_original_phylogenetic_source(
            t,ref,lookup,original_island,original_island,tip,D,
            n=15,repeats=20,seed=9)


def test_no_source_tip_or_unknown_taxon_has_no_biological_effect():
    t=original_like()
    ref=RarityReference.from_original(t,null_draws=199,seed=3)
    lookup=_compiled(ref,t)
    names=list(t.index)
    tips={name:i for i,name in enumerate(names[:35])}
    D=np.abs(np.arange(35)[:,None]-np.arange(35)[None,:]).astype(float)
    with pytest.raises(ValueError,match="Insufficient"):
        _sample_original_phylogenetic_source(t,ref,lookup,set(names[:30]),
           set(names[50:80]),tips,D,n=15,repeats=20)
