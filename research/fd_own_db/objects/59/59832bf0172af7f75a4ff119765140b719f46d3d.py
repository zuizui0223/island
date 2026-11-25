"""Original dated-tree exact-tip island versus candidate mainland floral coupling.

This is phylogenetically bounded *floristic* assembly, NOT ancestral source,
within-species trait evolution, local pollinator switch, or measured pollen
service. Only exact species tips with independently agreed source families
may be paired; inserted genus/family taxa are NOT used.

Source plant species present on a target island must be removed from EVERY
potential GIFT mainland donor region BEFORE matching. CoL family candidate
identity is strictly checked upstream against the original GBIF family.
"""
from __future__ import annotations

from hashlib import sha256
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

from .regional_native_flora_fd import (
    TRAITS, METRICS, exact_rarefied_flora_fd
)
from .regional_native_provenance_volume import cluster_equal_bootstrap

SAMPLE_N=5
REPEATS=24
CONDITIONS=(
    ("same_family_tree_le50",50.,False),
    ("same_family_tree_le100",100.,False),
    ("same_genus_tree_le100",100.,True),
)
SCOPES=("all_recorded","regional_WCVP_compatible","strict_source_native")


def _seed(*tokens:object)->int:
    return int.from_bytes(sha256(":".join(map(str,tokens)).encode()).digest()[:8],"big")


def _genus(name:str)->str:
    v=name.split()
    return v[0] if len(v)>=2 else ""


def _matched_rows(matrix:np.ndarray)->np.ndarray:
    """Original one-to-one source assignment, no trait-dependent edge weights."""
    x=np.asarray(matrix,dtype=bool)
    if x.ndim!=2:
        raise ValueError("Expected botanical island×candidate mainland species graph")
    if not x.size or x.shape[1]==0:
        return np.full(len(x),-1,dtype=np.int32)
    return maximum_bipartite_matching(
        csr_matrix(x.astype(np.int8)),
        perm_type="column").astype(np.int32)


def _score(names:list[str], direct:dict[str,tuple[str,str,str]])->dict:
    if len(set(names))!=SAMPLE_N or not set(names)<=set(direct):
        raise ValueError("Five different exact source Direct plant species required")
    return exact_rarefied_flora_fd(
        np.asarray([direct[x] for x in names],dtype=object),n=SAMPLE_N)


def _original_pool_matrices(
    focal:set[str],donor:set[str], *,
    families:dict[str,str],tip_index:dict[str,int],D:np.ndarray
)->tuple[list[str],list[str],np.ndarray,np.ndarray,np.ndarray]:
    if focal & donor:
        raise ValueError("Identical island species leaked into original mainland donor")
    valid=set(tip_index)&set(families)
    left=sorted(focal&valid)
    right=sorted(donor&valid)
    if not left or not right:
        return left,right,np.zeros((len(left),len(right)),float),np.zeros((len(left),len(right)),bool),np.zeros((len(left),len(right)),bool)
    ix=np.asarray([tip_index[x] for x in left],dtype=np.int32)
    iy=np.asarray([tip_index[x] for x in right],dtype=np.int32)
    distance=D[np.ix_(ix,iy)]
    if not np.isfinite(distance).all() or (distance<0).any():
        raise ValueError("Original published dated tree distance matrix invalid")
    fam=(np.asarray([families[x] for x in left],dtype=object)[:,None] ==
         np.asarray([families[x] for x in right],dtype=object)[None,:])
    genus=(np.asarray([_genus(x) for x in left],dtype=object)[:,None] ==
           np.asarray([_genus(x) for x in right],dtype=object)[None,:])
    return left,right,distance,fam,genus


def sample_true_tip_pairing(
    left:list[str],right:list[str],distance:np.ndarray,
    same_fam:np.ndarray,same_genus:np.ndarray,
    direct:dict[str,tuple[str,str,str]], *,
    threshold:float,require_genus:bool=False,
    repeats:int=REPEATS,seed:int=20261114,
)->dict:
    """Same FIVE focal plants per draw vs original bounded source and family control."""
    if threshold not in (50.,100.) or repeats<20:
        raise ValueError("Original dated tree thresholds or repeated draws drifted")
    if distance.shape!=(len(left),len(right)):
        raise ValueError("Tree matrix row/column source taxon identities not aligned")
    if same_fam.shape!=distance.shape or same_genus.shape!=distance.shape:
        raise ValueError("Source family/genus masks not aligned with branch distances")
    if len(set(left))!=len(left) or len(set(right))!=len(right) or set(left)&set(right):
        raise ValueError("Source species duplicate, or island/source name overlap")
    allowed=same_fam & (distance<=threshold)
    if require_genus:
        allowed &=same_genus
    initial=_matched_rows(allowed)
    eligible=np.flatnonzero(initial>=0)
    if len(eligible)<SAMPLE_N:
        return {"support":int(len(eligible)),"estimable":False}
    rng=np.random.default_rng(seed)
    obs=[]
    for _ in range(repeats):
        picked=rng.choice(eligible,size=SAMPLE_N,replace=False)
        # Permute donor column order to avoid selecting only the first matched
        # donor whenever multiple candidate original published tips qualify.
        order=rng.permutation(len(right))
        assigned=_matched_rows(allowed[np.ix_(picked,order)])
        if (assigned<0).any():
            raise AssertionError("Original source matchability contract violated")
        donor_i=order[assigned]
        ctrl_order=rng.permutation(len(right))
        ctrl=_matched_rows(same_fam[np.ix_(picked,ctrl_order)])
        if (ctrl<0).any():
            raise AssertionError("Family-only control must include bounded source pairs")
        donor_c=ctrl_order[ctrl]
        focal=[left[int(i)] for i in picked]
        bounded=[right[int(i)] for i in donor_i]
        control=[right[int(i)] for i in donor_c]
        fv=_score(focal,direct)
        bv=_score(bounded,direct)
        cv=_score(control,direct)
        z=dict(
            mean_original_bounded_tree_distance=float(np.mean(distance[picked,donor_i])),
            mean_original_family_control_tree_distance=float(np.mean(distance[picked,donor_c])),
            fraction_bounded_pairs_same_genus=float(np.mean(same_genus[picked,donor_i])),
            n_distinct_original_genera_focal=int(len(set(
                np.asarray([_genus(x) for x in focal],dtype=object)))),
        )
        for m in METRICS:
            z["bounded_"+m]=fv[m]-bv[m]
            z["family_control_"+m]=fv[m]-cv[m]
            z["bounded_minus_family_control_"+m]=cv[m]-bv[m]
        obs.append(z)
    mean=pd.DataFrame(obs).mean().to_dict()
    return {"support":int(len(eligible)),"estimable":True,
            **{str(k):float(v) for k,v in mean.items()}}


def sample_same_focal_distance_bands(
    left:list[str],right:list[str],distance:np.ndarray,
    same_fam:np.ndarray,direct:dict[str,tuple[str,str,str]], *,
    repeats:int=REPEATS,seed:int=20261114,
)->dict:
    """Same 5 island plants and different distinct source species in d<=50 vs (50,100].

    Two source windows cannot claim change of local evolutionary traits, as
    both donors and their documented colour/form/symmetry categories differ.
    """
    near=same_fam&(distance<=50)
    mid=same_fam&(distance>50)&(distance<=100)
    support=np.flatnonzero(near.any(axis=1)&mid.any(axis=1))
    if len(support)<SAMPLE_N:
        return {"n_local_true_tips_with_two_distance_bands":int(len(support)),
                "estimable":False}
    rng=np.random.default_rng(seed)
    sample=[]
    for _ in range(repeats*15):
        if len(sample)>=repeats:break
        picked=rng.choice(support,size=SAMPLE_N,replace=False)
        nperm=rng.permutation(len(right))
        mperm=rng.permutation(len(right))
        ni=_matched_rows(near[np.ix_(picked,nperm)])
        mi=_matched_rows(mid[np.ix_(picked,mperm)])
        if (ni<0).any() or (mi<0).any():
            continue
        nindex=nperm[ni];mindex=mperm[mi]
        focal=_score([left[int(i)] for i in picked],direct)
        nval=_score([right[int(i)] for i in nindex],direct)
        mval=_score([right[int(i)] for i in mindex],direct)
        z=dict(
            mean_near_tree_distance=float(distance[picked,nindex].mean()),
            mean_mid_tree_distance=float(distance[picked,mindex].mean()),
        )
        for metric in METRICS:
            z["focal_minus_near_"+metric]=focal[metric]-nval[metric]
            z["focal_minus_mid_"+metric]=focal[metric]-mval[metric]
            z["near_minus_mid_source_sensitivity_"+metric]=mval[metric]-nval[metric]
        sample.append(z)
    if len(sample)<20:
        return {"n_local_true_tips_with_two_distance_bands":int(len(support)),
                "n_completed_draws":len(sample),"estimable":False}
    avg=pd.DataFrame(sample).mean().to_dict()
    return {
        "n_local_true_tips_with_two_distance_bands":int(len(support)),
        "n_completed_draws":len(sample),"estimable":True,
        **{str(k):float(v) for k,v in avg.items()}
    }


def original_true_tip_native_pairing(
    direct:pd.DataFrame,floras:pd.DataFrame,geo:pd.DataFrame,
    nearby:dict[str,tuple[str,...]],giftdonors:dict[str,set[str]],
    families:dict[str,str],name_to_tip:dict[str,int],
    tip_order:list[str],distance:np.ndarray,*,
    tree_version:str,repeats:int=REPEATS,seed:int=20261114,
    cluster_draws:int=999,
)->dict:
    """Compute rights-safe original source region/taxon/dated-tree summaries."""
    if len(direct)!=2234 or list(direct.columns)!=list(TRAITS):
        raise ValueError("Original 2,234 species-direct complete triple not verified")
    if len(floras.island_id.unique())!=928 or floras.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original 928 island flora identity or incidence changed")
    if tree_version not in ("TPL","LCVP","WP") or len(families)<1900:
        raise ValueError("Published source original tree or taxonomy unverified")
    if len(tip_order)!=len(name_to_tip) or set(name_to_tip)!=set(tip_order):
        raise ValueError("Original tip index names mismatch")
    if distance.shape!=(len(tip_order),len(tip_order)):
        raise ValueError("Dated tree distances shape mismatch")
    if cluster_draws<999 or repeats<20:
        raise ValueError("Original output uncertainty or sampling reduced")
    if len(nearby)!=928 or len(set(geo.island_id))!=928:
        raise ValueError("Original geographic/source assignments incomplete")
    g=geo.set_index("island_id")
    reference={
        str(k):tuple(map(str,v))
        for k,v in zip(direct.index,direct.to_numpy(dtype=str))
    }
    support_rows=[]
    result_rows=[]
    for island,raw in floras.groupby("island_id",sort=True):
        island=str(island)
        allnames=set(raw.accepted_species)
        donor=(set().union(*(giftdonors.get(t,set()) for t in nearby[island]))&
               set(reference))-allnames
        if len(nearby[island])!=5 or len(set(nearby[island]))!=5:
            raise ValueError("Original GIFT native five-source IDs missing")
        tiers=(
            ("all_recorded",allnames),
            ("regional_WCVP_compatible",set(raw.loc[
                raw.regional_native_compatible,"accepted_species"])),
            ("strict_source_native",set(raw.loc[
                raw.strict_source_native,"accepted_species"])),
        )
        for scope,own in tiers:
            if len(own)<30:
                continue
            local,source,D,samefam,samegenus=_original_pool_matrices(
                own,donor,families=families,tip_index=name_to_tip,D=distance)
            for label,max_distance,genus_required in CONDITIONS:
                mask=samefam&(D<=max_distance)&(samegenus if genus_required else True)
                match=_matched_rows(mask)
                n_match=int((match>=0).sum())
                row=dict(
                    island_id=island,
                    scope=scope,regime=str(g.at[island,"analysis_regime"]),
                    spatial_block=str(g.at[island,"spatial_block"]),
                    source_signature=nearby[island][0],
                    comparison=label,
                    matched_focal_tip_opportunity=n_match,
                    n_tree_tip_plant_species_locally=len(local),
                    n_source_tip_plant_species=len(source),
                    source_support_at_least_five=bool(n_match>=SAMPLE_N),
                )
                support_rows.append(row)
                if n_match<SAMPLE_N or len(donor)<30:
                    continue
                rr=sample_true_tip_pairing(
                    local,source,D,samefam,samegenus,reference,
                    threshold=max_distance,require_genus=genus_required,
                    repeats=repeats,seed=_seed(seed,tree_version,island,scope,label))
                if not rr["estimable"]:
                    raise AssertionError("Original prequalified true tip source unexpectedly unavailable")
                result_rows.append({**row,**{k:v for k,v in rr.items() if k!="support"}})
            # Identical SAME 5 focal plants, two source tree-distance bands.
            if len(local)>=5 and len(source)>=5:
                pair=sample_same_focal_distance_bands(
                    local,source,D,samefam,reference,
                    repeats=repeats,seed=_seed(seed,tree_version,island,scope,"dual"))
                srow=dict(
                    island_id=island,scope=scope,regime=str(g.at[island,"analysis_regime"]),
                    spatial_block=str(g.at[island,"spatial_block"]),
                    source_signature=nearby[island][0],
                    comparison="two_disjoint_distance_bands",
                    matched_focal_tip_opportunity=pair["n_local_true_tips_with_two_distance_bands"],
                    n_tree_tip_plant_species_locally=len(local),
                    n_source_tip_plant_species=len(source),
                    source_support_at_least_five=bool(pair["estimable"])
                )
                support_rows.append(srow)
                if pair["estimable"] and len(donor)>=30:
                    result_rows.append({**srow,**{
                        k:v for k,v in pair.items()
                        if k not in ("n_local_true_tips_with_two_distance_bands",
                                     "n_completed_draws","estimable")
                    }})
    support=pd.DataFrame(support_rows)
    out=pd.DataFrame(result_rows)
    results=[]
    focus=("colour_architecture_pairing_departure",
           "marginal_rao","joint_rao")
    for keys,sp in support.groupby(
        ["scope","regime","comparison"],sort=True
    ):
        scope,regime,comparison=keys
        eligible=sp.loc[sp.source_support_at_least_five]
        rows=out.loc[
            out.scope.eq(scope)&out.regime.eq(regime)&out.comparison.eq(comparison)]
        # Cases with insufficient original source species might pass taxon
        # overlap but are NOT final eligible. Count tested outputs explicitly.
        n=int(len(rows))
        blocks=int(rows.spatial_block.nunique())
        sources=int(rows.source_signature.nunique())
        exploratory=n>=10 and blocks>=5 and sources>=5
        strong=n>=20 and blocks>=10 and sources>=10
        record=dict(
            scope=scope,regime=regime,tree_version=tree_version,
            comparison=comparison,
            n_original_islands_over_30_direct=len(sp),
            n_original_taxon_matchable_islands_pre_gift_pool=n if n==len(eligible) else int(len(eligible)),
            n_source_and_tree_matchable_islands=n,
            n_physical_spatial_blocks=blocks,
            n_first_GIFT_source_signatures=sources,
            exploratory_support_10island_5block_5source=bool(exploratory),
            strict_support_20island_10block_10source=bool(strong),
        )
        if exploratory:
            if comparison=="two_disjoint_distance_bands":
                cols=[prefix+m for m in focus for prefix in (
                    "focal_minus_near_","focal_minus_mid_",
                    "near_minus_mid_source_sensitivity_")]
            else:
                cols=[prefix+m for m in focus for prefix in (
                    "bounded_","family_control_",
                    "bounded_minus_family_control_")]
                cols+=["mean_original_bounded_tree_distance",
                       "mean_original_family_control_tree_distance",
                       "fraction_bounded_pairs_same_genus"]
            for col in cols:
                for lab,cluster in (
                    ("block","spatial_block"),("source","source_signature")):
                    record[lab+"_"+col]=cluster_equal_bootstrap(
                        rows,col,cluster,draws=cluster_draws,
                        seed=_seed(seed,"bootstrap",tree_version,scope,
                                   regime,comparison,col,lab))
        results.append(record)
    if len(results)<15:
        raise ValueError("Original restricted tree analysis lacks scope/region panels")
    return dict(
        status="ORIGINAL_TRUE_TIP_PHYLO_BOUNDED_NATIVE_PAIRING",
        original_eligible_islands=928,
        n_original_Direct_triple_species=2234,
        n_original_strict_native_eligible=int(
            (floras.groupby("island_id").strict_source_native.sum()>=30).sum()),
        n_original_WCVP_compatible_eligible=int(
            (floras.groupby("island_id").regional_native_compatible.sum()>=30).sum()),
        tree_version=tree_version,
        n_original_true_tip_source_species=len(tip_order),
        original_frozen_seed=seed,
        design=dict(
            matched_number_per_source_and_focal_island=SAMPLE_N,
            phylogenetic_cutoffs_published_original_branch_units=[50.,100.],
            source_gift_candidates_geographical_not_true_ancestral_source=True,
            exact_original_accepted_tip_species_only=True,
            no_taxonomically_inferred_backbone_tip=True,
            family_candidate_unreviewed_but_source_concordant=True,
            original_GIFT_donor_excludes_all_focal_plant_species=True,
            outcome_blind_family_genus_tree_distance_source_selection=True,
            same_focal_five_plant_source_near_vs_mid_distance_bands=True,
            independent_two_numeric_seeds_not_independent_evolution_history=True,
            WCVP_regional_native_does_not_prove_local_island_native=True,
            no_pollination_function_or_historical_guild_replacement_estimated=True,
            no_original_source_species_taxon_island_tip_ids_exported=True,
        ),
        regional_conditional_phylo_pairing=results,
        original_restricted_source_taxa_and_dated_backbone_not_exported=True,
    )
