"""True original branch PhyloSor beta and floral support Jaccard are separate metrics."""
from itertools import product
import numpy as np
import pandas as pd
import pytest

from fd.phylo_beta_functional_beta_macro import (
    OriginalTreeSubtree,functional_reference_beta,_weighted_slope,
    _phylo_functional_draw,_family_stratified_flower_trait_shuffle,
    ProfileCountFunctionalBetaNull
)


def original_tiny_tree():
    edge=np.array([[7,5],[7,6],[5,1],[5,2],[6,3],[6,4]],
                  dtype=np.int64)
    lengths=np.array([2.,3.,1.,1.5,2.,1.],dtype=float)
    return OriginalTreeSubtree(edge,lengths,np.array([1,2,3,4]))


def test_original_rooted_phylo_beta_uses_shared_branch_length_not_kde():
    p=original_tiny_tree()
    # [1,3] and [2,4] use both major backbone branches, but
    # exact tips differ. Shared original branch length 2+3=5.
    beta=p.phylosor_beta([1,3],[2,4])
    assert beta==pytest.approx(1-2*5/15.5)
    assert p.phylosor_beta([1,2],[3,4])==pytest.approx(1.)
    assert p.phylosor_beta([1,3],[1,3])==pytest.approx(0.)
    with pytest.raises(ValueError,match="Non-original"):
        p.phylosor_beta([1,5],[2,3])


def test_source_bitset_jaccard_functional_beta_is_not_profile_count():
    assert functional_reference_beta(0b1011,0b1101)==pytest.approx(0.5)
    assert functional_reference_beta(0b1001,0b1001)==pytest.approx(0.)
    assert functional_reference_beta(0b0011,0b1100)==pytest.approx(1.)
    with pytest.raises(ValueError,match="No actually observed"):
        functional_reference_beta(0,0)


def synthetic_original_region():
    rng=np.random.default_rng(13)
    n=45
    x=np.linspace(.13,.78,n)+rng.normal(0,.018,n)
    g=rng.uniform(.1,.8,n)
    f=rng.uniform(.15,.72,n)
    source=rng.uniform(3,8,n)
    area=rng.uniform(2,9,n)
    dist=rng.uniform(1,7,n)
    data=pd.DataFrame(dict(
        phylosor_beta=x,
        functional_support_beta=.38*x+.1*g+.025*f+rng.normal(0,.025,n),
        original_genus_jaccard_turnover=g,
        original_colxr_candidate_family_jaccard_turnover=f,
        log_source_exact_tip_candidates=source,
        log_island_area_km2=area,
        log_distance_to_continent_km=dist,
        spatial_block=[f"block{i%15}" for i in range(n)],
        source_signature=[f"source{i%12}" for i in range(n)],
    ))
    return data


def test_phylo_functional_correlation_cluster_weighted_not_repetition_pseudoreplication():
    z=synthetic_original_region()
    for cluster in ("spatial_block","source_signature"):
        raw=_weighted_slope(z,"phylosor_beta","functional_support_beta",
                            cluster,adjusted=False)
        adjusted=_weighted_slope(z,"phylosor_beta","functional_support_beta",
                                 cluster,adjusted=True)
        assert raw>0
        assert adjusted>0
        assert np.isfinite(adjusted)


def test_not_enough_independent_blocks_is_not_a_zero_slope():
    z=synthetic_original_region()
    z.spatial_block="one_block"
    with pytest.raises(ValueError,match="Insufficient"):
        _weighted_slope(z,"phylosor_beta","functional_support_beta",
                        "spatial_block",adjusted=False)


def test_same_fifteen_tree_tip_flower_triples_not_single_colour_strings():
    from fd.discrete_functional_volume import TRAITS
    from fd.rare_profile_decomposition import RarityReference
    profiles=list(product(
       ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
       ("open_radial","tubular","brush","composite","zygomorphic"),
       ("actinomorphic","zygomorphic")))
    weights=np.arange(1,101,dtype=float)**1.4
    counts=np.floor(weights/weights.sum()*(2234-100)).astype(int)+1
    for i in range(2234-int(counts.sum())):
        counts[i%100]+=1
    rows=[prof for prof,num in zip(profiles,counts) for _ in range(int(num))]
    names=[f"Genus{i%35} species_{i}" for i in range(len(rows))]
    assert len(rows)==2234
    traits=pd.DataFrame(rows,index=names,columns=TRAITS)
    ref=RarityReference.from_original(traits,null_draws=199)
    # Fully rooted binary 32-tip published-tree-like edge graph.
    edges=[]
    clusters=list(range(1,33))
    parent=33
    while len(clusters)>1:
        new=[]
        for i in range(0,len(clusters),2):
            edges.extend([(parent,clusters[i]),(parent,clusters[i+1])])
            new.append(parent)
            parent+=1
        clusters=new
    tree=OriginalTreeSubtree(np.asarray(edges,int),
                             np.ones(len(edges)),np.arange(1,33))
    tip={name:i+1 for i,name in enumerate(names[:32])}
    cat={name:tuple(traits.loc[name].values) for name in names[:32]}
    family={name:"SameSyntheticFamily" for name in names[:32]}
    profile_k=ProfileCountFunctionalBetaNull(
        ref,draws=199,seed=7,nmax=15)
    r=_phylo_functional_draw(ref,tree,names[:15],names[15:30],
                             tip,cat,family,k_profile_null=profile_k)
    assert 0<=r["phylosor_beta"]<=1
    assert 0<=r["functional_support_beta"]<=1
    assert 0<=r["exact_profile_jaccard_beta"]<=1
    assert np.isfinite(r["k_adjusted_rare_minus_common_source_gap"])
    assert np.isfinite(r["functional_beta_excess_given_exact_profile_k"])
    assert np.isfinite(r["functional_beta_excess_given_exact_profile_k_speciesfreq_null"])
    assert 0<=profile_k.expect_uniform[15,15]<=1
    assert profile_k.excess(
        r["functional_support_beta"],
        r["exact_original_island_profile_count_k"],
        r["exact_original_source_profile_count_k"]
    )[0]==pytest.approx(r["functional_beta_excess_given_exact_profile_k"])
    assert r["n_original_shared_families"]==1


def test_family_preserving_null_shuffles_whole_three_trait_states_with_fixed_donors():
    from fd.discrete_functional_volume import TRAITS
    from fd.rare_profile_decomposition import RarityReference
    profiles=list(product(
       ("white","purple","yellow","red","blue","pink","green","orange","brown","cream"),
       ("open_radial","tubular","brush","composite","zygomorphic"),
       ("actinomorphic","zygomorphic")))
    w=np.arange(1,101,dtype=float)**1.4
    counts=np.floor(w/w.sum()*(2234-100)).astype(int)+1
    for i in range(2234-int(counts.sum())):counts[i%100]+=1
    rows=[p for p,c in zip(profiles,counts) for _ in range(int(c))]
    names=[f"Genus{i%70} source_{i}" for i in range(2234)]
    t=pd.DataFrame(rows,index=names,columns=TRAITS)
    ref=RarityReference.from_original(t,null_draws=199)
    subset=names[:1162]
    p={name:tuple(t.loc[name].values) for name in subset}
    fam={name:f"Family{i%40}" for i,name in enumerate(subset)}
    original=[
      dict(scope="all_recorded",regime="tropical",
         spatial_block=f"b{i%10}",source_signature=f"s{i%12}",
         log_source_exact_tip_candidates=5.,log_island_area_km2=3.,
         log_distance_to_continent_km=2.,
         phylosor_beta=float(i)/45.,
         functional_support_beta=.3,
         original_genus_jaccard_turnover=.5,
         original_colxr_candidate_family_jaccard_turnover=.5,
         species_pairs=[(subset[i:i+15],subset[i+400:i+415]) for _ in range(5)])
      for i in range(30)
    ]
    null,audit=_family_stratified_flower_trait_shuffle(
        ref,p,fam,original,seed=14,permutations=99)
    assert len(null)==99 and len(null[0])==30
    assert audit["n_original_source_families_permutable"]==40
    assert audit["n_original_true_tip_species_in_permutable_families"]==1162
    assert audit["preserve_exact_original_within_family_phenotype_counts"] is True
    assert np.isfinite(np.asarray(null)).all()
    assert np.array_equal(np.asarray(null),np.asarray(
        _family_stratified_flower_trait_shuffle(
          ref,p,fam,original,seed=14,permutations=99)[0]))
