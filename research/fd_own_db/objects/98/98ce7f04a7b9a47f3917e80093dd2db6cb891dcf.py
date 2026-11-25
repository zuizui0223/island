"""Pure synthetic family identity and genus-switch flower space safeguards."""
from itertools import product
import numpy as np
import pandas as pd
import pytest
from fd.discrete_functional_volume import TRAITS
from fd.rare_profile_decomposition import RarityReference,genus
from fd.family_conditioned_rare_volume import (
    matched_family_source_effect, independent_family_richness_effect,
    _score, _compiled, _mask,
    _shapley_family_scores, _families_with_different_genus_pairs
)
from fd.colxr_family_opportunity import group_by_family


def original_like():
    patterns=list(product(
        ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
        ("open_radial","tubular","brush","composite","zygomorphic"),
        ("actinomorphic","zygomorphic")))
    counts=(np.arange(1,101,dtype=float)**1.3)
    counts=np.floor(counts/counts.sum()*(2234-100)).astype(int)+1
    for j in range(2234-int(counts.sum())):
        counts[j%100]+=1
    names=[]
    rows=[]
    family={}
    k=0
    for p,count in zip(patterns,counts):
        for _ in range(int(count)):
            g=k%70
            n=f"Genus{g} species_{k}"
            names.append(n)
            rows.append(p)
            family[n]=f"Family{g%22}"
            k+=1
    assert len(rows)==2234 and len(set(names))==2234
    return pd.DataFrame(rows,index=names,columns=TRAITS),family


def test_exact_matched_family_allows_genus_or_forces_different_genus():
    t,families=original_like()
    ref=RarityReference.from_original(t,null_draws=199,seed=2)
    look=_compiled(ref,t)
    # Split alternating species WITHIN each genus, not a global even/odd
    # stream (global parity would give each genus to just one study site).
    a={n for n in t.index if int(n.split("species_")[1])//70%2==0}
    b=set(t.index)-a
    src=_families_with_different_genus_pairs(
        group_by_family(a,families),group_by_family(b,families))
    assert len(src)>=8
    for n in (5,8):
        free=matched_family_source_effect(ref,look,a,b,
              family_map=families,n=n,repeats=20,seed=7)
        diff=matched_family_source_effect(ref,look,a,b,
              family_map=families,n=n,repeats=20,seed=7,
              must_change_genus=True)
        assert set(free)==set(diff)
        assert len(free)==5
        assert all(np.isfinite(x) for x in free.values())
        assert all(np.isfinite(x) for x in diff.values())
        assert diff==matched_family_source_effect(ref,look,a,b,
              family_map=families,n=n,repeats=20,seed=7,
              must_change_genus=True)
        free_identity=independent_family_richness_effect(
             ref,look,a,b,family_map=families,
             n=n,repeats=20,seed=17)
        assert set(free_identity)==set(free)
        assert free_identity==independent_family_richness_effect(
             ref,look,a,b,family_map=families,
             n=n,repeats=20,seed=17)
        assert all(np.isfinite(v) for v in free_identity.values())


def test_exact_original_30_species_family_shapley_additivity():
    t,families=original_like()
    ref=RarityReference.from_original(t,null_draws=199)
    look=_compiled(ref,t)
    # Split alternating species WITHIN each genus, not a global even/odd
    # stream (global parity would give each genus to just one study site).
    a={n for n in t.index if int(n.split("species_")[1])//70%2==0}
    b=set(t.index)-a
    score=_shapley_family_scores(ref,look,a,b,families,repeats=20,seed=7)
    for label in ("rare","common","gap"):
        assert score[label+"_total"]==pytest.approx(
            score[label+"_shared"]+score[label+"_unshared"],abs=1e-12)
    assert score["n_shared_candidate_families"]>=8
    assert _mask(sorted(list(a)[:20]),look)>=0


def test_original_family_identity_not_accepted_as_phylogenetic_distances():
    t,families=original_like()
    ref=RarityReference.from_original(t,null_draws=199)
    look=_compiled(ref,t)
    with pytest.raises(ValueError,match="exclude"):
        matched_family_source_effect(ref,look,set(t.index[:30]),
             set(t.index[:30]),family_map=families,repeats=20)
    empty={"OtherPlant a"}
    with pytest.raises(ValueError,match="Insufficient"):
        matched_family_source_effect(ref,look,set(t.index[:30]),
             empty,family_map=families,repeats=20)
