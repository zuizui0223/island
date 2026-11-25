"""Synthetic true-branch original family/near-genus source matching guards."""
from itertools import product
import numpy as np
import pandas as pd
import pytest

from fd.discrete_functional_volume import TRAITS
from fd.rare_profile_decomposition import RarityReference
from fd.family_conditioned_rare_volume import _compiled
from fd.true_tree_close_family_source import (
    _max_matching,_controlled_sample,_identical_focal_two_nonoverlapping_distance_bands
)


def _original_like():
    profiles=list(product(
        ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
        ("open_radial","tubular","brush","composite","zygomorphic"),
        ("actinomorphic","zygomorphic")))
    counts=np.arange(1,101,dtype=float)**1.2
    counts=np.floor(counts/counts.sum()*(2234-100)).astype(int)+1
    for j in range(2234-int(counts.sum())):counts[j%100]+=1
    species=[]
    traits=[]
    for j,(p,num) in enumerate(zip(profiles,counts)):
        for k in range(int(num)):
            species.append(f"Dummy{j%35} plant_{j}_{k}")
            traits.append(p)
    assert len(species)==2234
    local=[f"Alpha{i} island" for i in range(6)]
    donor=[f"Alpha{i} mainland" for i in range(6)]+[
        f"Other{i} mainland" for i in range(2)]
    names=species[14:]
    source=pd.DataFrame(traits,index=species,columns=TRAITS)
    first=source.iloc[:14].copy()
    first.index=local+donor
    triple=pd.concat([first,source.iloc[14:]])
    fam={name:f"Family{i%3}" for i,name in enumerate(local)}
    fam.update({donor[i]:f"Family{i%3}" for i in range(6)})
    fam[donor[6]]="Family0"
    fam[donor[7]]="Family1"
    assert len(triple)==2234 and len(set(triple.index))==2234
    D=np.full((6,8),300.0)
    for i in range(6):
        D[i,i]=30.
        D[i,6+i%2]=90.
        D[i,(i+3)%6]=80.  # different donor, SAME family, medium distance
    return triple,local,donor,fam,D


def test_true_original_bipartite_maximum_cardinality():
    m=np.array([[1,1,0],[0,1,0],[0,0,1]],dtype=bool)
    p=_max_matching(m)
    assert len(p)==3
    assert sorted(p.tolist())==[0,1,2]
    assert np.array_equal(_max_matching(np.zeros((3,0),bool)),
                          np.array([-1,-1,-1]))


def test_close_phylogenetic_5_real_tree_pairs_and_independent_same_family_source():
    triple,local,donor,family,D=_original_like()
    ref=RarityReference.from_original(triple,null_draws=199,seed=7)
    lookup=_compiled(ref,triple)
    r=_controlled_sample(ref,lookup,donor,local,D,family,
        threshold=50.,same_genus=True,n=5,repeats=20,seed=41)
    assert r["n_matched_unique_pairs_meeting_original_phylo_bound"]==6
    assert r["strict_n5_phylogenetic_source_pairs_available"]
    assert r["matched_tree_distance"]==pytest.approx(30.)
    assert r["family_control_tree_distance"]>=30
    assert np.isfinite(r["phylo_bounded_delta_k_adjusted_rare_minus_common"])
    repeat=_controlled_sample(ref,lookup,donor,local,D,family,
        threshold=50.,same_genus=True,n=5,repeats=20,seed=41)
    assert r==repeat


def test_not_enough_close_real_tree_pairs_is_no_effect_not_zero():
    triple,local,donor,fam,D=_original_like()
    ref=RarityReference.from_original(triple,null_draws=199,seed=3)
    lookup=_compiled(ref,triple)
    D[2:,:]=300
    r=_controlled_sample(ref,lookup,donor,local,D,fam,
       threshold=50.,same_genus=True,n=5,repeats=20,seed=5)
    assert r["n_matched_unique_pairs_meeting_original_phylo_bound"]==2
    assert r["strict_n5_phylogenetic_source_pairs_available"] is False
    assert "phylo_bounded_delta_k_adjusted_rare_minus_common" not in r
    with pytest.raises(ValueError,match="excluded"):
        _controlled_sample(ref,lookup,donor,donor,D,fam,
             threshold=50.,same_genus=True,n=5,repeats=20)

def test_same_five_real_focal_species_compared_across_two_nonoverlapping_source_bands():
    triples,local,donor,fam,D=_original_like()
    ref=RarityReference.from_original(triples,null_draws=199,seed=3)
    lookup=_compiled(ref,triples)
    result=_identical_focal_two_nonoverlapping_distance_bands(
        ref,lookup,local,donor,D,fam,seed=17,repeats=20)
    assert result["n_original_local_true_tip_species_with_two_donor_distance_bands"]==6
    assert result["qualified_exact_focal_taxa_for_two_distance_bands"] is True
    assert result["n_valid_fixed_same_focal_n5_repeats"]==20
    assert result["average_near_original_tree_distance"]==pytest.approx(30.)
    assert 80.<=result["average_mid_original_tree_distance"]<=90.
    again=_identical_focal_two_nonoverlapping_distance_bands(
        ref,lookup,local,donor,D,fam,seed=17,repeats=20)
    assert result==again
    D[2:,:]=300.
    invalid=_identical_focal_two_nonoverlapping_distance_bands(
        ref,lookup,local,donor,D,fam,seed=17,repeats=20)
    assert invalid["qualified_exact_focal_taxa_for_two_distance_bands"] is False
    assert "near_minus_mid_k_adjusted_rare_common" not in invalid
