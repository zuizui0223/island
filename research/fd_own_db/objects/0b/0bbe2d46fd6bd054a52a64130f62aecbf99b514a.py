"""True-backbone phylogenetic source matching for island floral rarity.

Exactly source-verified published GBOTB TPL/LCVP/WP tree tips, NOT inserted
phylogenetic tips or alphabetical genus/family distances. Geographic
candidate donor is the same five-nearest GIFT native mainland regions as
earlier 928-island floral analyses, with ALL focal-island species removed.

For each outcome-blind fixed n15/n10 original phylogeny-tip plant sample,
select different source plant species minimizing real dated PATristic
distance with a one-to-one Hungarian matching. Report independently
sampled source donor plant comparison on exactly SAME physical islands.
No transplant, environmental selection, phylogeny-caused flower evolution,
functional pollinator decline, or true ancestor is established.
"""
from __future__ import annotations

from collections import Counter
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from .discrete_functional_volume import TRAITS,_resample_group
from .family_conditioned_rare_volume import _compiled,_score,METRICS
from .rare_profile_decomposition import RarityReference,_seedstr

N_PRIMARY=15
N_SENSITIVITY=10
MATCH_KINDS=("nearest_phylo","random_geo_source")
PHYLO_RESTRICTED_MIN_ISLANDS=20
PHYLO_RESTRICTED_MIN_BLOCKS=10
PHYLO_RESTRICTED_MIN_DONOR_SOURCES=10


def _sample_original_phylogenetic_source(
    direct:pd.DataFrame,ref:RarityReference,lookup:dict,
    island_names:set[str],donor_names:set[str],tip_index:dict[str,int],
    distances:np.ndarray,*,n:int=15,repeats:int=20,seed:int=20261113,
    original_accepted_family:dict[str,str]|None=None
)->dict:
    if n not in (N_PRIMARY,N_SENSITIVITY) or repeats<10:
        raise ValueError("Frozen original phylogeny matching sampling contract changed")
    if island_names & donor_names:
        raise ValueError("Accepted donor species also observed on focal original island")
    a=sorted(island_names&set(tip_index))
    b=sorted(donor_names&set(tip_index))
    if len(a)<n or len(b)<n:
        raise ValueError("Insufficient ORIGINAL observed dated tree tips on island or donor")
    donor_idx=np.asarray([tip_index[z] for z in b],dtype=np.int32)
    rng=np.random.default_rng(seed)
    rows=[]
    for _ in range(repeats):
        selected=rng.choice(a,size=n,replace=False).tolist()
        i=np.asarray([tip_index[z] for z in selected],dtype=np.int32)
        matrix=distances[np.ix_(i,donor_idx)]
        if not np.isfinite(matrix).all() or matrix.min()<0:
            raise ValueError("Original published tree distances are invalid")
        left,right=linear_sum_assignment(matrix)
        if not np.array_equal(left,np.arange(n)):
            raise AssertionError("Source tree matching not one-to-one")
        donors=[b[int(j)] for j in right]
        random=rng.choice(b,size=n,replace=False).tolist()
        d_match=float(matrix[left,right].mean())
        d_random=float(np.mean(distances[
            i,np.asarray([tip_index[z] for z in random],dtype=np.int32)]))
        if d_match>d_random+1e-7:
            # A minimum-cost matching cannot be farther than an arbitrary
            # random original donor assignment on the SAME source pool.
            raise AssertionError("Minimum-cost dated phylogenetic matching failed")
        orig=_score(ref,selected,lookup)
        matched=_score(ref,donors,lookup)
        arbitrary=_score(ref,random,lookup)
        observed_selected_distances=matrix[left,right]
        family=(original_accepted_family or {})
        row=dict(
            d_matched=d_match,
            d_random=d_random,
            matched_minus_random_phylo_distance=d_match-d_random,
            # POST-HOC descriptive branch-distance support diagnostics.
            # The thresholds below are absolute units of the original
            # dated source tree, NOT proof of close relatives/ancestry.
            share_matched_pairs_le50_tree_length_units=float(
                np.mean(observed_selected_distances<=50.)),
            share_matched_pairs_le100_tree_length_units=float(
                np.mean(observed_selected_distances<=100.)),
            median_matched_pair_distance=float(np.median(observed_selected_distances)),
            maximum_matched_pair_distance=float(np.max(observed_selected_distances)),
            share_matched_pairs_same_original_accepted_genus=float(
                np.mean([p.split()[0]==q.split()[0]
                         for p,q in zip(selected,donors)])),
            share_matched_pairs_same_original_COLXR_candidate_family=float(
                np.mean([bool(family.get(p)) and family.get(p)==family.get(q)
                         for p,q in zip(selected,donors)])),
        )
        for j,k in enumerate(METRICS):
            row["matched_delta_"+k]=float(orig[j]-matched[j])
            row["random_delta_"+k]=float(orig[j]-arbitrary[j])
            row["matched_minus_random_delta_"+k]=float(
                arbitrary[j]-matched[j])
        rows.append(row)
    return pd.DataFrame(rows).mean(numeric_only=True).to_dict()


def _strict_source_taxon_groups(
     direct:pd.DataFrame,floras:pd.DataFrame,geo:pd.DataFrame,
     donors:pd.DataFrame,mainland:pd.DataFrame
):
    for table,cols in [
      (floras,{"island_id","accepted_species","origin_status"}),
      (geo,{"island_id","analysis_regime","spatial_block"}),
      (donors,{"island_id","source_mode","source_rank","entity_ID"}),
      (mainland,{"entity_ID","work_species"})
    ]:
        if not cols.issubset(table):
            raise ValueError("Missing frozen original island or GIFT source")
    x=floras.loc[floras.accepted_species.isin(direct.index)].copy()
    if x.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original source accepted species repeated on island")
    if not x.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Unknown source-native-status is not biological absence")
    c=x.groupby("island_id").accepted_species.nunique()
    eligible=set(c[c>=30].index)
    if len(eligible)!=928:
        raise ValueError("Frozen original 928 Direct island panel drifted")
    x=x.loc[x.island_id.isin(eligible)]
    geom=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates().set_index("island_id")
    if geom.index.has_duplicates or not eligible.issubset(set(geom.index)):
        raise ValueError("Original geography missing or duplicated")
    donors=donors.copy()
    donors["source_rank"]=pd.to_numeric(donors["source_rank"],errors="raise")
    if donors["source_rank"].isna().any():
        raise ValueError("Original GIFT donor source rank missing")
    assignments=donors.loc[
       donors.island_id.isin(eligible)&donors.source_mode.eq("geo_k5")&
       donors.source_rank.between(1,5),["island_id","source_rank","entity_ID"]]
    if assignments.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Candidate native mainland donor source rank duplicate")
    original_source=assignments.sort_values(["island_id","source_rank"]).groupby(
       "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(original_source.get(i,()))!=5 or
       len(set(original_source[i]))!=5 for i in eligible):
        raise ValueError("One of frozen five potential GIFT source IDs is missing")
    dictionary=(mainland.loc[mainland.work_species.isin(direct.index),
       ["entity_ID","work_species"]].drop_duplicates().groupby(
       "entity_ID").work_species.agg(set).to_dict())
    return x,geom,original_source,dictionary,eligible


def original_true_phylogeny_match(
    direct:pd.DataFrame,floras:pd.DataFrame,geo:pd.DataFrame,
    source_assignment:pd.DataFrame,source_native:pd.DataFrame,
    tip_name_to_node_id:dict[str,int],tip_species_order:list[str],
    distance_matrix:np.ndarray,*,
    tree_version:str,seed:int=20261113,
    original_accepted_family:dict[str,str]|None=None,
    rarefactions:int=20,reference_draws:int=399,bootstrap:int=1999
)->dict:
    if rarefactions<15 or reference_draws<199 or bootstrap<999:
        raise ValueError("Original branch-length support null contract cannot be weakened")
    if tree_version not in ("TPL","LCVP","WP"):
        raise ValueError("Unknown original dated phylogenetic backbone")
    if not set(tip_species_order).issubset(direct.index) or len(set(tip_species_order))!=len(tip_species_order):
        raise ValueError("Original tree exact tips include unknown or duplicate Direct plant taxa")
    if not set(tip_name_to_node_id).issuperset(tip_species_order):
        raise ValueError("Measured source tree TIP identity missing")
    if distance_matrix.shape!=(len(tip_species_order),len(tip_species_order)):
        raise ValueError("Original published tree pairwise distance matrix missing")
    if not np.isfinite(distance_matrix).all() or np.any(distance_matrix<0):
        raise ValueError("Invalid original dated branch-length distances")
    tip_index={s:i for i,s in enumerate(tip_species_order)}
    original,loc,near,gift,universe=_strict_source_taxon_groups(
        direct,floras,geo,source_assignment,source_native)
    byisland={i:g for i,g in original.groupby("island_id")}
    ref=RarityReference.from_original(direct,null_draws=reference_draws,seed=seed)
    lookup=_compiled(ref,direct)
    all_stats=[]
    support=[]
    for island in sorted(universe):
        own=byisland[island]
        all_names=set(own.accepted_species)
        donor=(set().union(*(gift.get(entity,set()) for entity in near[island]))&
               set(direct.index))-all_names
        count_donor=len(donor&set(tip_index))
        for scope,local in (
            ("all_recorded",all_names),
            ("source_native",set(own.loc[
                own.origin_status.eq("native"),"accepted_species"]))
        ):
            if len(local)<30:
                continue
            n_local_exact=len(local&set(tip_index))
            common=dict(island_id=island,scope=scope,
               regime=str(loc.at[island,"analysis_regime"]),
               spatial_block=str(loc.at[island,"spatial_block"]),
               source_signature=near[island][0],
               n_direct_original_local=len(local),
               n_exact_tip_original_local=n_local_exact,
               n_exact_tip_source_donor=count_donor)
            for n in (N_PRIMARY,N_SENSITIVITY):
                elig=(n_local_exact>=n and count_donor>=n)
                scope_status=dict(common,n_real_original_tree_tips_per_side=n,
                    sufficient_direct_true_tip_species=bool(elig))
                support.append(scope_status)
                if not elig:continue
                summary=_sample_original_phylogenetic_source(
                    direct,ref,lookup,local,donor,tip_index,distance_matrix,
                    n=n,repeats=rarefactions,
                    seed=_seedstr(seed,tree_version,island,scope,n),
                    original_accepted_family=original_accepted_family)
                all_stats.append(dict(scope_status,**summary))
    summary=pd.DataFrame(all_stats)
    opportunities=pd.DataFrame(support)
    regions=[]
    for (scope,regime,n),q in opportunities.groupby([
       "scope","regime","n_real_original_tree_tips_per_side"]):
        eligible=q.loc[q.sufficient_direct_true_tip_species]
        report=dict(scope=scope,regime=regime,
            n_dated_true_tree_tips_per_side=int(n),
            n_original_30Direct_islands=int(len(q)),
            n_original_tree_tip_eligible_islands=len(eligible),
            n_physical_blocks=int(eligible.spatial_block.nunique()),
            n_first_GIFT_source_signatures=int(eligible.source_signature.nunique()),
            median_n_original_local_exact_tree_tips=float(q.n_exact_tip_original_local.median()),
            median_n_candidate_GIFT_source_exact_tree_tips=float(q.n_exact_tip_source_donor.median()))
        valid=(len(eligible)>=PHYLO_RESTRICTED_MIN_ISLANDS and
           eligible.spatial_block.nunique()>=PHYLO_RESTRICTED_MIN_BLOCKS and
           eligible.source_signature.nunique()>=PHYLO_RESTRICTED_MIN_DONOR_SOURCES)
        report["regional_effect_reportable_20island_10block_10sources"]=bool(valid)
        if valid:
            p=summary.loc[summary.scope.eq(scope)&
                          summary.regime.eq(regime)&
                          summary.n_real_original_tree_tips_per_side.eq(n)].copy()
            if len(p)!=len(eligible):
                raise ValueError("True-tree-tip source sampling dropped a matched island")
            keys=("d_matched","d_random","matched_minus_random_phylo_distance",
                  "share_matched_pairs_le50_tree_length_units",
                  "share_matched_pairs_le100_tree_length_units",
                  "median_matched_pair_distance",
                  "maximum_matched_pair_distance",
                  "share_matched_pairs_same_original_accepted_genus",
                  "share_matched_pairs_same_original_COLXR_candidate_family",
                  "matched_delta_k_adjusted_rare_minus_common",
                  "random_delta_k_adjusted_rare_minus_common",
                  "matched_minus_random_delta_k_adjusted_rare_minus_common",
                  "matched_delta_rare_coverage",
                  "random_delta_rare_coverage")
            for j,key in enumerate(keys):
                report["median_"+key]=float(p[key].median())
                for tier,cluster in (("block","spatial_block"),
                                      ("source","source_signature")):
                    mean,lo,hi=_resample_group(p,key,cluster,
                        seed+n*100+21*j+(tier=="source"),draws=bootstrap)
                    report[tier+"_equal_mean_"+key]=mean
                    report[tier+"_ci025_"+key]=lo
                    report[tier+"_ci975_"+key]=hi
        regions.append(report)
    return {
       "status":"TRUE_BACKBONE_EXACT_TIP_SOURCE_PHYLO_MATCHED_RARE_FD",
       "tree_version":tree_version,
       "n_original_direct_plant_species":len(direct),
       "n_original_strict_tip_species_in_tree":len(tip_species_order),
       "n_original_eligible_full_island_flora":len(universe),
       "n_original_source_matched_island_samples_private":len(summary),
       "region_phylogenetic_tip_matching_support":regions,
       "design":{
          "strict_exact_original_tree_tip_only":True,
          "no_phylogenetic_tip_hypothesis_inserted":True,
          "source_genera_or_families_not_used_to_infer_missing_tree_tips":True,
          "source_species_exact_island_names_removed":True,
          "true_published_branch_length_patristic_distances":True,
          "optimal_unique_nearest_phylogeny_source_species_per_tip":True,
          "random_independent_source_donor_control_same_original_island":True,
          "primary_exact_tree_tip_plants_per_side":N_PRIMARY,
          "sensitivity_exact_tree_tip_plants_per_side":N_SENSITIVITY,
          "sampling_repeats":rarefactions,
          "fixed_profile_count_k_null_draws":reference_draws,
          "block_and_source_conditional_bootstrap_draws":bootstrap,
          "original_original_globally_recorded_100_triple_profile_reference":True,
          "no_ancestral_source_proven_from_geographically_nearest_GIFT":True,
          "source_tree_tip_subset_is_geographically_selected":True,
          "rare_original_categories_may_not_be_in_tree_even_if_in_reference":True,
          "matched_phylo_distance_50_100_thresholds_posthoc_diagnostics_only":True,
          "original_colxr_family_and_genus_pair_support_diagnostics":True,
       },
       "demonstrated_phylogeny_caused_flower_loss_or_pollinator_shift":False,
       "original_tree_tip_island_taxon_or_source_name_exported":False,
    }
