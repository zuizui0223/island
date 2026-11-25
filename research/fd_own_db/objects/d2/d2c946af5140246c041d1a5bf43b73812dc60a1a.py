"""Phylogenetic beta vs functional beta of *the same* original island/source floras.

This is a cross-sectional macroecological source-conditional design, not a
test that pollinator loss caused island floral evolution. Exact published
GBOTB tree tips ONLY (no genus/family grafts), frozen original 928-island
high/medium Direct 3-trait florets, and original five-nearest GIFT candidate
mainland native source floras (all focal island names excluded).

Primary phylogenetic beta: true branch-length PhyloSor beta = 1 -
  2 * shared ROOT-TO-TIP SUBTREE branch length /
  (island subtree branch length + source subtree branch length).
Primary functional beta: Jaccard dissimilarity of the original globally
observed 100-profile reference that lies within ONE original category
mismatch from each sampled flora (finite support, NOT 17-D KDE).
Secondary endpoint: rare-minus-common original profile coverage difference
(island-source) corrected for exact number k of recorded profile types.

Both communities are thinned to the SAME N=15 actual published tree-tip
plant species on 20 repetitions per island. Compute ONE mean predictor and
outcome per physical island; repetitions are not biological independent data.
Within-region source-supply, geographic, genus-/family-composition guards,
whole original 10-degree blocks and first-nearest GIFT-source cluster
bootstrap. Correlation != evolutionary process, causality or real pollinators.

Only regional anonymous summaries are returned. Original species, donor,
genus, family, island name, coordinate and tree tip IDs stay in protected CI.
"""
from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
import hashlib
import numpy as np
import pandas as pd

from .discrete_functional_volume import TRAITS, _resample_group
from .rare_profile_decomposition import RarityReference, _seedstr, genus
from .true_backbone_phylo_source_match import _strict_source_taxon_groups

N_RAREFY=15
N_REPEATS=20
N_PERM=99
N_BS=999
REGION_MIN=20
BLOCK_MIN=10
SOURCE_MIN=10


class OriginalTreeSubtree:
    """Path-compressed *original* published rooted tree; no grafts."""
    def __init__(self,edge:np.ndarray,branch_lengths:np.ndarray,
                 original_tips:np.ndarray):
        edge=np.asarray(edge,dtype=np.int64)
        lengths=np.asarray(branch_lengths,dtype=float)
        tips=np.asarray(original_tips,dtype=np.int64)
        if edge.ndim!=2 or edge.shape[1]!=2 or len(lengths)!=len(edge):
            raise ValueError("Original source rooted tree and branches incompatible")
        if np.any(lengths<=0) or not np.isfinite(lengths).all():
            raise ValueError("Published original tree branch lengths missing/invalid")
        nodes=int(edge.max())+1
        if edge.min()<1 or nodes>400000 or len(edge)!=nodes-2:
            raise ValueError("Original GBOTB edge list not single rooted tree")
        parent=np.zeros(nodes,dtype=np.int32)
        weight=np.zeros(nodes,dtype=float)
        parent[edge[:,1]]=edge[:,0]
        weight[edge[:,1]]=lengths
        if np.any(np.bincount(edge[:,1],minlength=nodes)[1:]>1):
            raise ValueError("Source rooted tree has nodes with multiple parents")
        root=np.setdiff1d(np.unique(edge[:,0]),np.unique(edge[:,1]))
        if len(root)!=1:
            raise ValueError("Source root not unique")
        self.root=int(root[0])
        if len(set(tips))!=len(tips) or np.any(tips<1):
            raise ValueError("Original phylogenetic tip species mapping not unique")
        self.original_tip_paths={}
        for tip in tips:
            path=[]
            current=int(tip)
            seen=set()
            while current!=self.root:
                if current in seen or current<1 or current>=nodes or parent[current]==0:
                    raise ValueError("Original plant tree branch path invalid")
                seen.add(current)
                path.append(current)  # each node denotes its incoming edge
                current=int(parent[current])
            self.original_tip_paths[int(tip)]=frozenset(path)
        self.edge_lengths=weight

    def branch_set(self,tips:list[int]|np.ndarray)->set[int]:
        nodes=set()
        for tip in tips:
            if int(tip) not in self.original_tip_paths:
                raise ValueError("Non-original (grafted or absent) phylogenetic tip used")
            nodes.update(self.original_tip_paths[int(tip)])
        return nodes

    def phylosor_beta(self,tips_a:list[int]|np.ndarray,
                      tips_b:list[int]|np.ndarray)->float:
        a=self.branch_set(tips_a)
        b=self.branch_set(tips_b)
        la=float(sum(self.edge_lengths[list(a)]))
        lb=float(sum(self.edge_lengths[list(b)]))
        if la<=0 or lb<=0:
            raise ValueError("Published source flora has zero evolutionary branch length")
        shared=float(sum(self.edge_lengths[list(a&b)]))
        beta=1-2*shared/(la+lb)
        if beta<(-1e-9) or beta>1+1e-9:
            raise ValueError("PhyloSor original dated branch normalized beta invalid")
        return float(np.clip(beta,0,1))


def functional_reference_beta(mask_a:int,mask_b:int)->float:
    union=mask_a|mask_b
    if union==0:
        raise ValueError("No actually observed floral reference profile support")
    return 1-float((mask_a&mask_b).bit_count()/union.bit_count())


class ProfileCountFunctionalBetaNull:
    """Finite categorical geometry benchmark at EXACT (k island, k source).

    Hold the number of unique actually recorded original floral profiles on
    BOTH source floras fixed separately. A random reference replicate
    samples k distinct original 100-profile centers from the frozen global
    reference; then computes one-trait-radius reference-Jaccard beta.
    Both uniform profile centers and species-frequency-weighted distinct
    profile centers are frozen independent reference nulls.
    No genus/family/source phylogenetic history is controlled by this null.
    """
    def __init__(self,reference:RarityReference,*,
                 draws:int=199,seed:int=20261114,nmax:int=15):
        if draws<199 or nmax!=15 or len(reference.volume.profiles)!=100:
            raise ValueError("Exact-k categorical β reference model changed")
        profiles=reference.volume.profiles
        masks=[reference.volume.masks[p][1] for p in profiles]
        frequencies=reference.profile_n_species.astype(float)
        probabilities=frequencies/frequencies.sum()
        self.expect_uniform={}
        self.expect_speciesfreq={}
        rng=np.random.default_rng(seed)
        for ka in range(1,nmax+1):
            for kb in range(1,nmax+1):
                u=[]
                w=[]
                for _ in range(draws):
                    for output,p in ((u,None),(w,probabilities)):
                        first=rng.choice(len(masks),size=ka,replace=False,p=p)
                        second=rng.choice(len(masks),size=kb,replace=False,p=p)
                        x=y=0
                        for j in first:x|=masks[int(j)]
                        for j in second:y|=masks[int(j)]
                        output.append(functional_reference_beta(x,y))
                self.expect_uniform[ka,kb]=float(np.mean(u))
                self.expect_speciesfreq[ka,kb]=float(np.mean(w))
        self.draws=draws

    def excess(self,observed_beta:float,ka:int,kb:int):
        if (ka,kb) not in self.expect_uniform:
            raise ValueError("Original floral distinct profile count k unsupported")
        return (
            float(observed_beta-self.expect_uniform[ka,kb]),
            float(observed_beta-self.expect_speciesfreq[ka,kb])
        )


def _phylo_functional_draw(
    reference:RarityReference,tree:OriginalTreeSubtree,local:list[str],
    donor:list[str],tip_map:dict[str,int],taxon_profiles:dict[str,tuple],
    original_family:dict[str,str],
    k_profile_null:ProfileCountFunctionalBetaNull|None=None,
):
    if len(local)!=N_RAREFY or len(donor)!=N_RAREFY:
        raise ValueError("Both source floras must have exactly 15 original taxa")
    if set(local)&set(donor):
        raise ValueError("Focal accepted island species not disjoint from GIFT source")
    if any(x not in tip_map or x not in taxon_profiles for x in local+donor):
        raise ValueError("One source floral species not an ORIGINAL published dated tree tip")
    pA=[taxon_profiles[x] for x in local]
    pB=[taxon_profiles[x] for x in donor]
    a={z for z in pA}
    b={z for z in pB}
    ma=0
    mb=0
    for t in a:ma|=reference.volume.masks[t][1]
    for t in b:mb|=reference.volume.masks[t][1]
    functional=functional_reference_beta(ma,mb)
    direct_profile_jaccard=1-len(a&b)/len(a|b)
    ea=np.array([(ma&reference.bins[c]).bit_count()/reference.n_bins[c]
                for c in ("rare","middle","common")],dtype=float)
    eb=np.array([(mb&reference.bins[c]).bit_count()/reference.n_bins[c]
                for c in ("rare","middle","common")],dtype=float)
    null_a=reference.expected_uniform[len(a)]
    null_b=reference.expected_uniform[len(b)]
    signed=float((ea[0]-ea[2])-(eb[0]-eb[2]))
    k_signed=signed-float(
        (null_a[0]-null_a[2])-(null_b[0]-null_b[2]))
    phy=tree.phylosor_beta([tip_map[x] for x in local],
                           [tip_map[x] for x in donor])
    gensA={genus(x) for x in local}
    gensB={genus(x) for x in donor}
    famA={original_family[x] for x in local if original_family.get(x)}
    famB={original_family[x] for x in donor if original_family.get(x)}
    k_excess=(k_profile_null.excess(functional,len(a),len(b))
              if k_profile_null is not None else (None,None))
    return dict(
        phylosor_beta=phy,
        functional_support_beta=functional,
        functional_beta_excess_given_exact_profile_k=k_excess[0],
        functional_beta_excess_given_exact_profile_k_speciesfreq_null=k_excess[1],
        exact_original_island_profile_count_k=len(a),
        exact_original_source_profile_count_k=len(b),
        exact_profile_jaccard_beta=float(direct_profile_jaccard),
        signed_rare_minus_common_source_gap=signed,
        k_adjusted_rare_minus_common_source_gap=k_signed,
        original_genus_jaccard_turnover=1-len(gensA&gensB)/len(gensA|gensB),
        original_colxr_candidate_family_jaccard_turnover=(
           1-len(famA&famB)/len(famA|famB) if famA|famB else np.nan),
        n_original_shared_genera=len(gensA&gensB),
        n_original_shared_families=len(famA&famB),
    )


def _weighted_slope(rows:pd.DataFrame,xcol:str,ycol:str,cluster:str,
                    *,adjusted:bool)->float:
    if len(rows)<20 or rows[cluster].nunique()<10:
        raise ValueError("Insufficient independent original clusters for regression")
    d=rows.copy()
    x=d[xcol].to_numpy(float)
    y=d[ycol].to_numpy(float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite original phylogeny/functional response")
    sd=float(np.std(x))
    if sd<1e-12:
        raise ValueError("No physical-island variation in phylogenetic beta")
    z=(x-np.mean(x))/sd
    predictors=[z]
    if adjusted:
        for key in (
          "original_genus_jaccard_turnover",
          "original_colxr_candidate_family_jaccard_turnover",
          "log_source_exact_tip_candidates",
          "log_island_area_km2",
          "log_distance_to_continent_km",
        ):
            values=d[key].to_numpy(float)
            ss=np.std(values)
            if ss<1e-12:
                raise ValueError("No original regional confounder variation")
            predictors.append((values-values.mean())/ss)
    X=np.column_stack([np.ones(len(d)),*predictors])
    # Strict original source regional/regression covariance rank gate. Do not
    # silently drop a confounder to rescue a directional finding.
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Original phylogenetic/richness confounders collinear")
    weights=1/d.groupby(cluster)[cluster].transform("size").to_numpy(float)
    weights=weights/weights.sum()
    root=np.sqrt(weights)
    return float(np.linalg.lstsq(X*root[:,None],y*root,rcond=None)[0][1])


def _boot_beta(rows:pd.DataFrame,ycol:str,*,cluster:str,adjusted:bool,
               draws:int,seed:int):
    """Efficient original-cluster-equal WLS via sufficient statistics.

    Region's originally observed SD scaling is FROZEN during block bootstrap
    so coefficients remain in the same units and resampling cannot recode an
    apparent effect by changing the denominator.
    """
    n_groups=rows[cluster].nunique()
    if len(rows)<20 or n_groups<10:
        return None
    x=rows["phylosor_beta"].to_numpy(float)
    y=rows[ycol].to_numpy(float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return dict(status="NONESTIMABLE_NONFINITE_SOURCE_OUTCOME")
    sd=float(np.std(x))
    if sd<1e-12:
        return dict(status="NONESTIMABLE_CONSTANT_SOURCE_PHYLO_BETA")
    predictors=[(x-float(x.mean()))/sd]
    if adjusted:
        for key in (
            "original_genus_jaccard_turnover",
            "original_colxr_candidate_family_jaccard_turnover",
            "log_source_exact_tip_candidates",
            "log_island_area_km2",
            "log_distance_to_continent_km",
        ):
            v=rows[key].to_numpy(float)
            s=np.std(v)
            if not np.isfinite(v).all() or s<1e-12:
                return dict(status="NONESTIMABLE_ADJUSTMENT_SOURCE_RANK")
            predictors.append((v-v.mean())/s)
    X=np.column_stack((np.ones(len(rows)),*predictors))
    if np.linalg.matrix_rank(X)<X.shape[1]:
        return dict(status="NONESTIMABLE_ADJUSTMENT_SOURCE_RANK")
    levels=list(rows.groupby(cluster,sort=True).indices.values())
    ps=X.shape[1]
    xtx=[]
    xty=[]
    for ids in levels:
        xx=X[ids]
        yy=y[ids]
        xtx.append(xx.T@xx/len(ids))
        xty.append(xx.T@yy/len(ids))
    a=np.asarray(xtx)
    b=np.asarray(xty)
    A=a.sum(axis=0)
    B=b.sum(axis=0)
    if np.linalg.matrix_rank(A)<ps:
        return dict(status="NONESTIMABLE_SOURCE_CLUSTER_RANK")
    observed=float(np.linalg.lstsq(A,B,rcond=None)[0][1])
    rng=np.random.default_rng(seed)
    estimates=[]
    for _ in range(draws):
        selected=rng.integers(len(levels),size=len(levels))
        aa=a[selected].sum(axis=0)
        if np.linalg.matrix_rank(aa)<ps:
            continue
        bb=b[selected].sum(axis=0)
        estimates.append(float(np.linalg.lstsq(aa,bb,rcond=None)[0][1]))
    if len(estimates)<0.8*draws:
        return dict(status="NONESTIMABLE_BOOTSTRAP_PREDICTOR_RANK",
                    n_valid=len(estimates),
                    observed_std_phylo_beta_slope=observed)
    lo,hi=np.quantile(estimates,[.025,.975])
    return dict(status="EXPLORATORY_CLUSTER_CONDITIONAL_ASSOCIATION",
                observed_std_phylo_beta_slope=observed,
                cluster_025=float(lo),cluster_975=float(hi),
                n_valid_block_resamples=len(estimates))


def _family_stratified_flower_trait_shuffle(
    ref:RarityReference,
    true_tip_flower_profiles:dict[str,tuple],
    exact_tip_family:dict[str,str],
    fixed_island_draws:list[dict],
    *,seed:int=20261114,permutations:int=99,
)->list[dict]:
    """POST-HOC original cross-island composition null.

    Every original accepted true-tree-tip plant retains its original
    phylogeny, genus/family, island/GIFT memberships, record opportunities
    and source sampling history. Within EACH original independently accepted
    CoL XR family, shuffle the COMPLETE original 3-trait phenotype tuple
    among its source member true-tree-tip accepted species. Thus exact family
    frequency of colour/form/symmetry triple observations is fixed.
    Do not shuffle individual trait axes or invent new phenotypes.
    """
    if permutations!=99 or len(true_tip_flower_profiles)<1100:
        raise ValueError("Prespecified original family-profile permutation contract changed")
    groups=defaultdict(list)
    for name in sorted(true_tip_flower_profiles):
        if name not in exact_tip_family or not exact_tip_family[name]:
            raise ValueError("Original accepted true tip family unverified")
        groups[exact_tip_family[name]].append(name)
    permutable_species=sum(len(v) for v in groups.values() if len(v)>=2)
    permutable_families=sum(1 for v in groups.values() if len(v)>=2)
    if permutable_species<300 or permutable_families<20:
        raise ValueError("Family-conditioned original flower traits cannot be permuted enough")
    lookup={name:ref.volume.masks[profile][1]
            for name,profile in true_tip_flower_profiles.items()}
    if len(lookup)!=len(true_tip_flower_profiles):
        raise ValueError("Original Direct source tree accepted plant names repeated")
    chosen=[]
    # Family profile permutations happen GLOBAL and ONCE per draw, not
    # independently per island, preserving cross-island source structure.
    rng=np.random.default_rng(seed)
    for iteration in range(permutations):
        trial=lookup.copy()
        for family,names in groups.items():
            if len(names)<2:continue
            permutation=rng.permutation(len(names))
            masks=[lookup[name] for name in names]
            for j,name in enumerate(names):
                trial[name]=masks[int(permutation[j])]
        null_beta=[]
        for r in fixed_island_draws:
            values=[]
            for lhs,rhs in r["species_pairs"]:
                a=0
                b=0
                for x in lhs:a |= trial[x]
                for x in rhs:b |= trial[x]
                values.append(functional_reference_beta(a,b))
            null_beta.append(float(np.mean(values)))
        chosen.append(null_beta)
    # NONE of the original recorded individual species names, original
    # families or island geographic labels are published from this function.
    return chosen,dict(
       n_global_exact_tip_source_taxonomic_families=int(len(groups)),
       n_original_true_tip_species_in_permutable_families=int(permutable_species),
       n_original_source_families_permutable=int(permutable_families),
       n_whole_three_trait_profile_within_family_shuffles=permutations,
       preserve_exact_original_within_family_phenotype_counts=True,
       keep_island_and_GIFT_membership_and_phylogenetic_tree_fixed=True,
       source_family_classification_is_unreviewed_colxr_candidate=True,
       source_profile_permutations_are_posthoc_sensitivity_not_phylogenetic_selection=True,
    )


def _conditional_family_shuffle_comparison(
    fixed_island_draws:list[dict],null_beta:list[list[float]],
    permutation_audit:dict,*,
    min_islands:int=20,min_blocks:int=10,min_sources:int=10,
)->list[dict]:
    if not fixed_island_draws or len(null_beta)!=99:
        raise ValueError("Original source rarefaction and family shuffle result unavailable")
    d=pd.DataFrame([{k:v for k,v in r.items() if k!="species_pairs"}
                    for r in fixed_island_draws])
    if len(d)!=len(null_beta[0]):
        raise ValueError("Family trait null rows do not align with actual original islands")
    null=np.asarray(null_beta,dtype=float)
    if null.shape!=(99,len(d)) or not np.isfinite(null).all():
        raise ValueError("Not all family-shuffled source profiles have valid functional beta")
    result=[]
    for (scope,regime),original in d.groupby(["scope","regime"]):
        inds=original.index.to_numpy()
        physical=original.spatial_block.nunique()
        source=original.source_signature.nunique()
        enough=len(inds)>=min_islands and physical>=min_blocks and source>=min_sources
        row=dict(scope=scope,regime=regime,
            n_original_islands=int(len(original)),n_physical_blocks=int(physical),
            n_first_GIFT_source_signatures=int(source),
            supported_20islands_10blocks_10sources=bool(enough),
            n_within_family_global_flower_triple_null_reassignments=99,
            same_true_tip_species_draws_in_observed_and_all_nulls=True,
        )
        if enough:
            for adjusted in (False,True):
                for block in ("spatial_block","source_signature"):
                    tag=("partial" if adjusted else "raw")+"_"+(
                         "block" if block=="spatial_block" else "source")
                    try:
                        actual=_weighted_slope(original,"phylosor_beta",
                            "functional_support_beta",block,adjusted=adjusted)
                    except ValueError:
                        row[tag]=dict(status="ORIGINAL_SOURCE_TAXONOMIC_COVARIANCE_RANK_INSUFFICIENT")
                        continue
                    permuted=[]
                    for j in range(99):
                        frame=original.copy()
                        frame["functional_support_beta"]=null[j,inds]
                        try:
                            permuted.append(_weighted_slope(frame,"phylosor_beta",
                                 "functional_support_beta",block,adjusted=adjusted))
                        except ValueError:
                            pass
                    if len(permuted)<80:
                        row[tag]=dict(status="FAMILY_TRAIT_NULL_EFFECT_UNIDENTIFIABLE",
                                      n_valid_permutations=len(permuted))
                        continue
                    mean=float(np.mean(permuted))
                    left,right=np.quantile(permuted,[.025,.975])
                    # Two diagnostics: excess over family-preserving null
                    # (one-sided), plus an empirical two-sided comparison
                    # centered on the observed family-preserving null mean.
                    excess=float(actual-mean)
                    p_upper=(1+sum(x>=actual for x in permuted))/(1+len(permuted))
                    p_two=(1+sum(abs(x-mean)>=abs(excess) for x in permuted))/(1+len(permuted))
                    row[tag]=dict(
                        status="POSTHOC_WITHIN_FAMILY_WHOLE_FLORAL_PROFILE_SHUFFLE",
                        observed_original_std_phylogenetic_beta_slope=float(actual),
                        family_trait_shuffle_mean_slope=mean,
                        family_trait_shuffle_025=float(left),
                        family_trait_shuffle_975=float(right),
                        observed_minus_family_shuffle_mean_slope=excess,
                        n_original_family_profile_permutations=len(permuted),
                        empirical_upper_tail_excess_probability=float(p_upper),
                        empirical_two_sided_centered_probability=float(p_two),
                        inference_not_independent_causal_evolution=True)
        result.append(row)
    return result


def macro_original_phylo_functional_beta(
    triple:pd.DataFrame,floras:pd.DataFrame,geo:pd.DataFrame,
    original_source_assignment:pd.DataFrame,
    original_GIFT_native_mainland:pd.DataFrame,
    tree:OriginalTreeSubtree,exact_original_tip_names:list[str],
    exact_tip_node_map:dict[str,int], original_family:dict[str,str],
    *,tree_version:str,repeats:int=20,reference_null_draws:int=199,
    bootstrap:int=999,seed:int=20261114,
)->dict:
    if tree_version not in ("TPL","LCVP","WP"):
        raise ValueError("No original dated GBOTB source version")
    if repeats<20 or reference_null_draws<199 or bootstrap<999:
        raise ValueError("Original phylo/functional beta protocol support weakened")
    if len(triple)!=2234 or len(exact_original_tip_names)<1100:
        raise ValueError("Pinned original Direct/source phylogenetic species cohort changed")
    if len(set(exact_original_tip_names))!=len(exact_original_tip_names):
        raise ValueError("Source phylogenetic exact tip names duplicated")
    original,places,gift_src,gift,universe=_strict_source_taxon_groups(
        triple,floras,geo,original_source_assignment,original_GIFT_native_mainland)
    if len(universe)!=928:
        raise ValueError("Original 928 island Direct macro population changed")
    if not {"island_id","log_island_area_km2",
            "log_distance_to_continent_km"}.issubset(geo):
        raise ValueError("Original island area and geographical isolation needed")
    geographic=geo.set_index("island_id")
    ref=RarityReference.from_original(triple,null_draws=reference_null_draws,seed=seed)
    geometry_k_null=ProfileCountFunctionalBetaNull(
        ref,draws=reference_null_draws,seed=seed+888,nmax=15)
    phen={
        str(name):tuple(row) for name,row in zip(
            triple.index,triple[list(TRAITS)].to_numpy(dtype=str))}
    true=set(exact_original_tip_names)
    if not true.issubset(phen) or not true.issubset(original_family):
        raise ValueError("Only original source verified family/trait true tips may enter")
    if true!=set(exact_tip_node_map) or (
        set(exact_tip_node_map.values())!=set(tree.original_tip_paths)):
        raise ValueError("Dated source original tree tip-node mapping not strictly equivalent")
    per_island={i:v for i,v in original.groupby("island_id")}
    store=[]
    family_null_draws=[]
    for island in sorted(universe):
        original_flora=per_island[island]
        full=set(original_flora.accepted_species)
        pool=(set().union(*(gift.get(x,set()) for x in gift_src[island]))
              & true)-full
        if len(pool)<N_RAREFY:
            continue
        for scope,names in (
            ("all_recorded",full),
            ("source_native",set(original_flora.loc[
                original_flora.origin_status.eq("native"),"accepted_species"]))):
            if len(names)<30:
                continue
            local=sorted(names&true)
            donor=sorted(pool)
            if len(local)<N_RAREFY:
                continue
            random=np.random.default_rng(_seedstr(seed,tree_version,island,scope))
            reps=[]
            first5=[]
            for draw in range(repeats):
                A=random.choice(local,size=N_RAREFY,replace=False).tolist()
                B=random.choice(donor,size=N_RAREFY,replace=False).tolist()
                results=_phylo_functional_draw(
                    ref,tree,A,B,exact_tip_node_map,phen,original_family,
                    k_profile_null=geometry_k_null)
                reps.append(results)
                if draw<5:
                    first5.append((A,B))
            summarized=pd.DataFrame(reps).mean(numeric_only=True).to_dict()
            original_five=pd.DataFrame(reps[:5]).mean(numeric_only=True).to_dict()
            family_null_draws.append(dict(
                scope=scope,regime=str(places.at[island,"analysis_regime"]),
                spatial_block=str(places.at[island,"spatial_block"]),
                source_signature=gift_src[island][0],
                log_source_exact_tip_candidates=float(np.log1p(len(pool))),
                log_island_area_km2=float(
                    geographic.at[island,"log_island_area_km2"]),
                log_distance_to_continent_km=float(
                    geographic.at[island,"log_distance_to_continent_km"]),
                species_pairs=first5,**original_five))
            store.append(dict(
                island_id=island,scope=scope,
                regime=str(places.at[island,"analysis_regime"]),
                spatial_block=str(places.at[island,"spatial_block"]),
                source_signature=gift_src[island][0],
                log_source_exact_tip_candidates=float(np.log1p(len(pool))),
                log_island_area_km2=float(
                    geographic.at[island,"log_island_area_km2"]),
                log_distance_to_continent_km=float(
                    geographic.at[island,"log_distance_to_continent_km"]),
                **summarized))
    df=pd.DataFrame(store)
    reports=[]
    responses=("functional_support_beta",
       "functional_beta_excess_given_exact_profile_k",
       "functional_beta_excess_given_exact_profile_k_speciesfreq_null",
       "exact_profile_jaccard_beta",
       "k_adjusted_rare_minus_common_source_gap")
    for (scope,regime),g in df.groupby(["scope","regime"]):
        n=len(g)
        blocks=g.spatial_block.nunique()
        sources=g.source_signature.nunique()
        pass_gate=n>=REGION_MIN and blocks>=BLOCK_MIN and sources>=SOURCE_MIN
        report=dict(scope=scope,regime=regime,
            n_original_eligible_islands=n,n_physical_spatial_blocks=blocks,
            n_first_nearest_GIFT_source_signatures=sources,
            n_true_tip_focal_species_per_side=N_RAREFY,
            original_flora_repeats_per_island=repeats,
            passed_20islands_10blocks_10source_support=bool(pass_gate),
            median_original_phylo_beta=float(g.phylosor_beta.median()),
            median_original_functional_beta=float(g.functional_support_beta.median()),
            median_original_genus_turnover=float(
                g.original_genus_jaccard_turnover.median()),
            median_original_family_turnover=float(
                g.original_colxr_candidate_family_jaccard_turnover.median()),
        )
        if pass_gate:
            for i,y in enumerate(responses):
                for adj in (False,True):
                    for j,cluster in enumerate(("spatial_block","source_signature")):
                        result=_boot_beta(g,ycol=y,cluster=cluster,
                           adjusted=adj,draws=bootstrap,seed=seed+97*i+31*int(adj)+j)
                        report[f"{y}_{'adjusted' if adj else 'raw'}_"+
                            ("block" if j==0 else "source")]=result
        reports.append(report)
    family_null,diagnostics=_family_stratified_flower_trait_shuffle(
        ref,{x:phen[x] for x in true},original_family,family_null_draws,
        seed=seed+401,permutations=99)
    family_results=_conditional_family_shuffle_comparison(
        family_null_draws,family_null,diagnostics)
    return dict(
        status="ORIGINAL_TRUE_DATED_PHYLO_BETA_VS_CATEGORICAL_FUNCTIONAL_BETA",
        n_original_928_direct_island_population=928,
        n_original_global_direct_flower_species=2234,
        n_original_true_phylogenetic_tip_species=len(true),
        original_tree_backbone_version=tree_version,
        region_support_and_categorical_beta_relationships=reports,
        n_original_private_island_scope_matched_rows=len(df),
        family_preserving_whole_flower_trait_permutation_audit=diagnostics,
        original_family_preserving_flower_trait_null_by_region=family_results,
        original_design={
            "phylogenetic_beta":"one_minus_rooted_original_branch_length_PhyloSor",
            "functional_beta":"one_minus_Jaccard_original_reference_one_trait_neighborhood",
            "secondary":"signed_rare_minus_common_original_k_adjusted_coverage",
            "same_15_true_dated_plant_tips_on_island_and_GIFT_source":True,
            "n_per_island_original_source_random_rarefactions":repeats,
            "n_reference_observed_global_flower_profiles":100,
            "functional_Jaccard_beta_exact_island_source_profile_k_nulls":True,
            "functional_Jaccard_k_null_dissimilarity_reference_draws":reference_null_draws,
            "null_reference_uniform_and_global_direct_species_frequency":True,
            "null_does_not_preserve_source_taxonomy_or_phylogeny":True,
            "no_taxonomically_imputed_phylogenetic_plant_tips":True,
            "remove_all_same_island_species_from_source_GIFT_floras":True,
            "family_genus_source_area_isolation_adjusted_regressions":True,
            "equal_physical_block_and_GIFT_source_bootstraps":bootstrap,
            "one_original_physical_island_one_regression_row":True,
            "conditional_estimates_not_source_independent_ancestries":True,
            "original_plant_trait_assignments_are_species_fixed_not_evolving":True,
            "phylogenetic_source_distance_does_not_prove_ancestor":True,
            "family_preserving_original_whole_triple_profile_permutation_is_posthoc":True,
            "null_based_on_same_five_first_draws_not_main_20_draw_beta_regressions":True,
        },
        no_actual_pollinator_guild_decline_or_evolutionary_trait_change_measured=True,
        original_source_tip_genus_family_island_or_coordinate_ids_not_exported=True,
    )
