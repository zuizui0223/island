"""Synthetic exact categorical union-of-balls Shapley and one-per-genus checks."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from fd.discrete_functional_volume import TRAITS
from fd.rare_profile_decomposition import RarityReference
from fd.genus_shapley_functional_volume import (
    two_group_shapley, _ball_mask, _sample_components, _match_unique_genera
)


def synthetic_original():
    combos=list(product(
        ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
        ("open_radial","tubular","brush","composite","zygomorphic"),
        ("actinomorphic","zygomorphic")
    ))
    assert len(combos)==100
    counts=np.arange(1,101,dtype=float)**1.35
    counts=np.floor(counts/counts.sum()*(2234-100)).astype(int)+1
    while counts.sum()<2234:
        counts[np.argmin(counts)]+=1
    assert counts.sum()==2234
    species=[]
    values=[]
    j=0
    for combo,n in zip(combos,counts):
        for _ in range(int(n)):
            # Repeated names under genera to ensure genus-donor opportunities,
            # but never repeated accepted species identifiers.
            species.append(f"Genus{j%34} species_{j}")
            values.append(combo)
            j+=1
    return pd.DataFrame(values,index=species,columns=TRAITS)


def test_shapley_sum_coverage_exact_even_with_full_overlap():
    t=synthetic_original()
    ref=RarityReference.from_original(t,null_draws=199,seed=3)
    a=list(t.index[:20])
    b=list(t.index[10:40])
    sh=_ball_mask(ref,a,t)
    un=_ball_mask(ref,b,t)
    precompiled={str(n):ref.volume.masks[tuple(values)][1]
                 for n,values in zip(t.index,t[list(TRAITS)].to_numpy(dtype=str))}
    assert sh==_ball_mask(ref,a,t,original_species_mask=precompiled)
    assert un==_ball_mask(ref,b,t,original_species_mask=precompiled)
    r=two_group_shapley(ref,sh,un)
    for rarity in ("rare","common"):
        assert r[f"{rarity}_total"]==pytest.approx(
           r[f"{rarity}_shared"]+r[f"{rarity}_unshared"],abs=1e-12)
        assert 0<=r[f"{rarity}_shared"]<=1
        assert 0<=r[f"{rarity}_unshared"]<=1
    assert r["gap_total"]==pytest.approx(r["gap_shared"]+r["gap_unshared"],abs=1e-12)
    assert two_group_shapley(ref,sh,0)["gap_unshared"]==pytest.approx(0)


def test_genus_decomposition_source_species_disjoint_and_exact():
    t=synthetic_original()
    ref=RarityReference.from_original(t,null_draws=199,seed=6)
    # Allocate names within every genus, not by global even/odd name index.
    # The synthetic species generator cycles even-valued genus indices.
    local={name for name in t.index if int(name.split("species_")[1])//34%2==0}
    source=set(t.index)-local
    r=_sample_components(ref,t,local,source,n=30,repeats=20,seed=2)
    for key in ("rare","common","gap"):
        assert r[key+"_total"]==pytest.approx(
            r[key+"_shared"]+r[key+"_unshared"],abs=1e-12)
    assert r["n_shared_source_genera_in_original_flora"]>=5
    with pytest.raises(ValueError,match="leakage"):
        _sample_components(ref,t,local,local,n=30,repeats=20)


def test_one_source_different_species_per_genus_at_five_and_eight():
    t=synthetic_original()
    ref=RarityReference.from_original(t,null_draws=199,seed=7)
    # Allocate names within every genus, not by global even/odd name index.
    # The synthetic species generator cycles even-valued genus indices.
    local={name for name in t.index if int(name.split("species_")[1])//34%2==0}
    source=set(t.index)-local
    for n in (5,8):
        a=_match_unique_genera(ref,t,local,source,n=n,repeats=20,seed=4)
        b=_match_unique_genera(ref,t,local,source,n=n,repeats=20,seed=4)
        assert a==b
        assert a["n_global_source_shared_genera"]>=n
        assert all(np.isfinite(v) for v in a.values())
    with pytest.raises(ValueError,match="few DISTINCT"):
        _match_unique_genera(ref,t,local,set(),n=5,repeats=20)
