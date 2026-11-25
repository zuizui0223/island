"""Strict same-genus source: source-same-name ban and genus slot matching."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from fd.rare_profile_decomposition import RarityReference,genus
from fd.genus_matched_rare_volume import _genus_slots,_simulate


def _fixture():
    profiles=list(product(("white","purple","yellow","red","pink","blue","green","cream","brown","orange"),
                          ("open_radial","tubular","composite_head","brush","zygomorphic"),
                          ("actinomorphic","zygomorphic")))
    counts=1+np.arange(100)*4
    # Nonuniform global species reference with exact frozen 2,234 carriers.
    counts=np.floor(counts/counts.sum()*(2234-100)).astype(int)+1
    for j in range(2234-int(counts.sum())):
        counts[j%100]+=1
    triple=[]
    names=[]
    for i,(p,k) in enumerate(zip(profiles,counts)):
        for j in range(int(k)):
            triple.append(p)
            names.append(f"Genus{i%70} trait_{i}_{j}")
    assert len(names)==2234
    return pd.DataFrame(triple,index=names,columns=(
        "flower_primary_color","floral_form","floral_symmetry"))


def test_taxonomic_opportunity_disjoint_and_matching_is_exact():
    local={"A first","A second","B one","B two","C one"}
    mainland={"A first","A foreign","A another","B new","C new","X other"}
    a,b,slots=_genus_slots(local,mainland)
    assert set(b["A"])=={"A another","A foreign"}
    assert "A first" not in b["A"]
    assert len(slots)==4
    assert len([x for x in slots if x=="A"])==2
    assert len([x for x in slots if x=="B"])==1


def test_model_audit_requires_15_exact_within_source_genus_carriers():
    t=_fixture()
    r=RarityReference.from_original(t,null_draws=199)
    # Synthetic source is a different exact plant species in each genus.
    own=set(t.index[::2])
    other=set(t.index[1::2])
    a,b,slots=_genus_slots(own,other)
    assert len(slots)>=15
    result=_simulate(r,t,own,other,n=15,repeats=20,seed=9)
    assert result["matched"].shape==(5,)
    assert result["unrestricted"].shape==(5,)
    assert np.isfinite(result["matched"]).all()
    assert np.isfinite(result["unrestricted"]).all()
    assert result["n_genus_source_opportunity_slots"]>=15
    # Exact generator reproducibility, without comparing numpy array dictionaries.
    second=_simulate(r,t,own,other,n=15,repeats=20,seed=9)
    assert np.array_equal(result["matched"],second["matched"])
    assert np.array_equal(result["unrestricted"],second["unrestricted"])


def test_an_unmatchable_genus_does_not_become_an_assumed_pollen_loss():
    t=_fixture()
    r=RarityReference.from_original(t,null_draws=199)
    own={str(x) for x in t.index[:15]}
    other={"NeverListed foreign","Unknown source"}
    with pytest.raises(ValueError,match="Insufficient"):
        _simulate(r,t,own,other,n=15,repeats=20,seed=9)
