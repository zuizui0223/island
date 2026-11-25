"""Synthetic genus-richness (not same-genus) Direct floral trait space safeguards."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from fd.discrete_functional_volume import TRAITS
from fd.rare_profile_decomposition import RarityReference
from fd.genus_richness_rare_volume import (
    _by_genus,_sample_one_per_genus,_one_island,_one_score
)


def _original_like():
    patterns=list(product(
       ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
       ("open_radial","tubular","brush","composite","zygomorphic"),
       ("actinomorphic","zygomorphic")))
    weights=np.arange(1,101,dtype=float)**1.2
    counts=np.floor(weights/weights.sum()*(2234-100)).astype(int)+1
    for i in range(2234-int(counts.sum())):
        counts[i%100]+=1
    names=[]
    profiles=[]
    for k,(p,n) in enumerate(zip(patterns,counts)):
        for j in range(int(n)):
            names.append(f"Genus{k%38} plant_{k}_{j}")
            profiles.append(p)
    assert len(names)==2234
    return pd.DataFrame(profiles,index=names,columns=TRAITS)


def test_equal_genus_richness_does_not_force_genus_identities_to_match():
    island={f"GenusA{i} sp" for i in range(40)}
    mainland={f"GenusB{i} sp" for i in range(50)}
    local=_sample_one_per_genus(_by_genus(island),n=15,rng=np.random.default_rng(5))
    source=_sample_one_per_genus(_by_genus(mainland),n=15,rng=np.random.default_rng(6))
    assert len({x.split()[0] for x in local})==15
    assert len({x.split()[0] for x in source})==15
    assert {x.split()[0] for x in local}!={x.split()[0] for x in source}
    with pytest.raises(ValueError,match="Too few"):
        _sample_one_per_genus(_by_genus(set(sorted(island)[:4])),
                              n=15,rng=np.random.default_rng(1))


def test_compiled_original_species_bitmap_is_identical_to_direct_categorical_metric():
    x=_original_like()
    ref=RarityReference.from_original(x,null_draws=199)
    bm={name:(tuple(row),ref.volume.masks[tuple(row)][1])
        for name,row in zip(x.index,x.to_numpy(dtype=str))}
    names=list(x.index[:30])
    assert np.array_equal(_one_score(ref,x,names),
                          _one_score(ref,x,names,source_profile_bitmap=bm))


def test_same_n_distinct_genera_with_different_source_species_reproducible():
    x=_original_like()
    ref=RarityReference.from_original(x,null_draws=199)
    local=set(x.index[::2])
    source=set(x.index[1::2])
    bm={name:(tuple(row),ref.volume.masks[tuple(row)][1])
        for name,row in zip(x.index,x.to_numpy(dtype=str))}
    for n in (10,15):
        a=_one_island(ref,x,local,source,n=n,draws=20,seed=37,source_profile_bitmap=bm)
        b=_one_island(ref,x,local,source,n=n,draws=20,seed=37,source_profile_bitmap=bm)
        assert a==b
        assert all(np.isfinite(v) for v in a.values())
    with pytest.raises(ValueError,match="disjoint"):
        _one_island(ref,x,local,local,n=15,draws=20,seed=37)
