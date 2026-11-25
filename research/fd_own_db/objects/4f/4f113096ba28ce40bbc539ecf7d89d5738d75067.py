"""True-branch-length bounded, SAME-FAMILY source plant comparisons.

This POST-HOC diagnostic addresses an important weakness of the n15/n10
nearest-tree matching: in tropical islands only about 9% of pairs shared a
genus and about 55% a family. Restrict plant pairs to source taxonomically
matching original CoL XR candidate family AND <=50 or <=100 units along the
REAL ORIGINAL GBOTB dated tree. A same-genus+<=100 stricter arm is separate.

One maximum-cardinality bipartite matching is formed using ONLY taxa,
family and dated path lengths BEFORE looking at floral profile states.
Exactly 5 separate island taxa and 5 different source accepted plant taxa
are then randomly sampled 20 times from that outcome-blind pairing. The
control uses exactly the SAME 5 island taxa with independently selected
different mainland accepted plants of the SAME FAMILIES (not necessarily
<= distance cutoff), so it isolates phylogeny threshold sensitivity within
candidate family rather than counting unrelated mainland species.

No observed pollinator loss, actual ancestral source, within-species
selection, genetic sequences newly produced or population floral shifts.
"""
from __future__ import annotations
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

from .discrete_functional_volume import _resample_group
from .family_conditioned_rare_volume import _score,_compiled
from .rare_profile_decomposition import RarityReference,genus,_seedstr
from .true_backbone_phylo_source_match import _strict_source_taxon_groups

CRITERIA=(
    ("same_family_d_le100",100.,False),
    ("same_family_d_le50",50.,False),
    ("same_genus_d_le100",100.,True),
)
SAMPLE_N=5
MIN_REGIONAL_ISLANDS=20
MIN_REGIONAL_BLOCKS=10
MIN_REGIONAL_SOURCES=10


def _max_matching(admissible:np.ndarray)->np.ndarray:
    m=np.asarray(admissible,dtype=bool)
    if m.ndim!=2:
        raise ValueError("A genuine island×source bipartite graph is required")
    if not m.size or m.shape[0]==0 or m.shape[1]==0:
        return np.full(m.shape[0],-1,dtype=np.int32)
    # True support only, no weighted edge score based on floral phenotype.
    return maximum_bipartite_matching(csr_matrix(m.astype(np.int8)),
                                      perm_type="column").astype(np.int32)


def _controlled_sample(
    reference:RarityReference,lookup:dict,source_plant_names:list[str],
    island_plant_names:list[str],distance:np.ndarray,
    original_family:dict[str,str], *,
    threshold:float,same_genus:bool,
    n:int=5,repeats:int=20,seed:int=20261113,
)->dict:
    if n!=5 or repeats<10 or threshold not in (50.,100.):
        raise ValueError("Predeclared close phylogenetic source matching contract changed")
    if len(set(island_plant_names))!=len(island_plant_names) or (
        len(set(source_plant_names))!=len(source_plant_names)):
        raise ValueError("Source accepted plant names not unique")
    if set(island_plant_names)&set(source_plant_names):
        raise ValueError("Focal island species must be excluded from source")
    if distance.shape!=(len(island_plant_names),len(source_plant_names)):
        raise ValueError("Source actual tree distance matrix does not match taxa")
    if not np.isfinite(distance).all() or np.any(distance<0):
        raise ValueError("Nonfinite or negative original tree path lengths")
    a=np.asarray([original_family[z] for z in island_plant_names],dtype=str)
    b=np.asarray([original_family[z] for z in source_plant_names],dtype=str)
    same_fam=a[:,None]==b[None,:]
    if same_genus:
        ag=np.asarray([genus(z) for z in island_plant_names],dtype=str)
        bg=np.asarray([genus(z) for z in source_plant_names],dtype=str)
        valid=same_fam & (ag[:,None]==bg[None,:]) & (distance<=threshold)
    else:
        valid=same_fam & (distance<=threshold)
    matched=_max_matching(valid)
    eligible=np.flatnonzero(matched>=0)
    support=dict(
       n_original_island_exact_tip_plant_names=len(island_plant_names),
       n_original_source_exact_tip_plant_names=len(source_plant_names),
       n_matched_unique_pairs_meeting_original_phylo_bound=len(eligible),
       n_distinct_original_families_in_admissible_pairs=int(
           len(set(a[eligible]))),
       strict_n5_phylogenetic_source_pairs_available=bool(len(eligible)>=n),
    )
    if len(eligible)<n:
        return support
    rng=np.random.default_rng(seed)
    rows=[]
    for _ in range(repeats):
        picked=rng.choice(eligible,size=n,replace=False)
        original=[island_plant_names[int(i)] for i in picked]
        best_idx=matched[picked]
        paired=[source_plant_names[int(j)] for j in best_idx]
        if len(set(paired))!=n:
            raise AssertionError("Unique true-tip source pair matched twice")
        measured=distance[picked,best_idx]
        if np.max(measured)>threshold+1e-8:
            raise AssertionError("Original branch-length threshold violated")
        # Same original island species AND family IDs control, but do not
        # require the donor to fall inside the dated distance cutoff.
        source_perm=rng.permutation(len(source_plant_names))
        candidate=(same_fam[np.ix_(picked,source_perm)])
        randomized=_max_matching(candidate)
        if np.any(randomized<0):
            raise ValueError("Strict family-matched no-distance source absent")
        control=[source_plant_names[int(source_perm[int(j)])]
                 for j in randomized]
        if len(set(control))!=n:
            raise AssertionError("Family-only donor reused in control")
        ctrl_index={sp:j for j,sp in enumerate(source_plant_names)}
        random_d=distance[picked,np.asarray(
            [ctrl_index[z] for z in control],dtype=np.int32)]
        own_score=_score(reference,original,lookup)
        match_score=_score(reference,paired,lookup)
        control_score=_score(reference,control,lookup)
        row=dict(
          matched_tree_distance=float(np.mean(measured)),
          family_control_tree_distance=float(np.mean(random_d)),
          matched_vs_same_family_control_distance=float(np.mean(measured-random_d)),
          n_distinct_families_in_selected_true_tip_pairs=int(len(set(a[picked]))),
        )
        for j,key in enumerate((
             "rare_coverage","common_coverage","raw_rare_minus_common",
             "k_adjusted_rare_minus_common",
             "frequency_k_adjusted_rare_minus_common")):
            row["phylo_bounded_delta_"+key]=float(own_score[j]-match_score[j])
            row["same_family_unbounded_delta_"+key]=float(
                own_score[j]-control_score[j])
            row["bounded_minus_family_delta_"+key]=float(
                control_score[j]-match_score[j])
        rows.append(row)
    support.update({k:float(v) for k,v in pd.DataFrame(rows).mean(
       numeric_only=True).to_dict().items()})
    return support


def _identical_focal_two_nonoverlapping_distance_bands(
    ref,lookup,local,source,dist,family,*,seed:int,repeats:int=20
)->dict:
    """Outcome-blind source taxon gate, same five focal original plant tips.

    Each focal taxon has an original different-species source plant within
    <=50 AND a distinct source plant with distance (50,100] in same family.
    Both donor sets are independently one-to-one. This is posthoc and does
    NOT keep source donor identities fixed between distance bands.
    """
    if dist.shape!=(len(local),len(source)) or set(local)&set(source):
        raise ValueError("Original exact local taxa/source distance contract invalid")
    a=np.asarray([family[x] for x in local],dtype=str)
    b=np.asarray([family[x] for x in source],dtype=str)
    same=a[:,None]==b[None,:]
    near=same&(dist<=50)
    mid=same&(dist>50)&(dist<=100)
    possibilities=np.flatnonzero(near.any(axis=1)&mid.any(axis=1))
    result=dict(
       n_original_local_true_tip_species_with_two_donor_distance_bands=int(len(possibilities)),
       n_valid_fixed_same_focal_n5_repeats=0,
       qualified_exact_focal_taxa_for_two_distance_bands=False)
    if len(possibilities)<5:
        return result
    rng=np.random.default_rng(seed)
    z=[]
    for attempt in range(repeats*10):
        if len(z)>=repeats:
            break
        picked=rng.choice(possibilities,size=5,replace=False)
        perm_a=rng.permutation(len(source))
        perm_b=rng.permutation(len(source))
        an=_max_matching(near[np.ix_(picked,perm_a)])
        bm=_max_matching(mid[np.ix_(picked,perm_b)])
        if np.any(an<0) or np.any(bm<0):
            continue
        own=[local[int(x)] for x in picked]
        source_near=[source[int(perm_a[int(j)])] for j in an]
        source_mid=[source[int(perm_b[int(j)])] for j in bm]
        if len(set(source_near))!=5 or len(set(source_mid))!=5:
            raise AssertionError("Source true tip reused within a distance band")
        d1=dist[picked,perm_a[an]]
        d2=dist[picked,perm_b[bm]]
        if np.any(d1>50) or np.any(d2<=50) or np.any(d2>100):
            raise AssertionError("Actual source original tree distance bands overlapped")
        ownval=_score(ref,own,lookup)
        nearval=_score(ref,source_near,lookup)
        midval=_score(ref,source_mid,lookup)
        z.append(dict(
           local_minus_near_k_adjusted_rare_common=float(ownval[3]-nearval[3]),
           local_minus_mid_k_adjusted_rare_common=float(ownval[3]-midval[3]),
           near_minus_mid_k_adjusted_rare_common=float(midval[3]-nearval[3]),
           near_minus_mid_rare_coverage=float(midval[0]-nearval[0]),
           average_near_original_tree_distance=float(np.mean(d1)),
           average_mid_original_tree_distance=float(np.mean(d2))))
    result["n_valid_fixed_same_focal_n5_repeats"]=len(z)
    result["qualified_exact_focal_taxa_for_two_distance_bands"]=len(z)>=10
    if len(z)>=10:
        result.update({k:float(v) for k,v in pd.DataFrame(z).mean().items()})
    return result


def run_close_phylo_family_source(
    triple:pd.DataFrame,floras:pd.DataFrame,geo:pd.DataFrame,
    sources:pd.DataFrame,mainland:pd.DataFrame,
    strict_tip_species:list[str], true_dated_distance:np.ndarray,
    original_family:dict[str,str], *,
    tree_version:str,reference_draws:int=199,repeat:int=20,
    bootstrap:int=999,seed:int=20261114,
)->dict:
    if tree_version not in ("TPL","LCVP","WP") or repeat<10 or bootstrap<999:
        raise ValueError("Original published dated tree experiment contract changed")
    if len(set(strict_tip_species))!=len(strict_tip_species):
        raise ValueError("Duplicated original true phylogenetic tip")
    if true_dated_distance.shape!=(len(strict_tip_species),len(strict_tip_species)):
        raise ValueError("Actual original tree metric does not cover all true tips")
    original,geo_map,nearest,gift,islands=_strict_source_taxon_groups(
        triple,floras,geo,sources,mainland)
    ref=RarityReference.from_original(
        triple,null_draws=reference_draws,seed=seed)
    lookup=_compiled(ref,triple)
    tip_lookup={z:k for k,z in enumerate(strict_tip_species)}
    eligible_tips=set(tip_lookup)&set(original_family)
    byisland={z:x for z,x in original.groupby("island_id")}
    scored=[]
    fixed_local_bands=[]
    for island in sorted(islands):
        site=byisland[island]
        allplants=set(site.accepted_species)
        donor=(set().union(*(gift.get(z,set()) for z in nearest[island]))&
               set(triple.index))-allplants
        donor=sorted(donor & eligible_tips)
        if len(donor)<SAMPLE_N:
            continue
        donor_arr=np.asarray([tip_lookup[s] for s in donor],dtype=np.int32)
        for scope,names in (
           ("all_recorded",allplants),
           ("source_native",set(site.loc[
               site.origin_status.eq("native"),"accepted_species"]))):
            if len(names)<30:
                continue
            original_names=sorted(names & eligible_tips)
            if len(original_names)<SAMPLE_N:
                continue
            local_idx=np.asarray([tip_lookup[z] for z in original_names],
                                 dtype=np.int32)
            distances=true_dated_distance[np.ix_(local_idx,donor_arr)]
            common=dict(
                original_island_id=island,
                scope=scope,
                regime=str(geo_map.at[island,"analysis_regime"]),
                spatial_block=str(geo_map.at[island,"spatial_block"]),
                first_source_signature=nearest[island][0])
            for i,(criterion,cutoff,require_same_genus) in enumerate(CRITERIA):
                score=_controlled_sample(
                   ref,lookup,donor,original_names,distances,
                   original_family,threshold=cutoff,
                   same_genus=require_same_genus,
                   repeats=repeat,n=SAMPLE_N,
                   seed=_seedstr(seed,island,scope,criterion,tree_version))
                scored.append(dict(common,criterion=criterion,**score))
            dual=_identical_focal_two_nonoverlapping_distance_bands(
                ref,lookup,original_names,donor,distances,original_family,
                repeats=repeat,seed=_seedstr(seed,island,scope,tree_version,"same_focal"))
            fixed_local_bands.append(dict(common,**dual))
    raw=pd.DataFrame.from_records(scored)
    regions=[]
    for (scope,regime,criterion),g in raw.groupby(
         ["scope","regime","criterion"]):
        passed=g.loc[g.strict_n5_phylogenetic_source_pairs_available]
        ok=(
          len(passed)>=MIN_REGIONAL_ISLANDS and
          passed.spatial_block.nunique()>=MIN_REGIONAL_BLOCKS and
          passed.first_source_signature.nunique()>=MIN_REGIONAL_SOURCES)
        report=dict(
          scope=scope,regime=regime,criterion=criterion,
          original_full_30species_flora_with_any_original_true_tips=len(g),
          n_islands_with_5_real_tree_tip_close_pairs=int(len(passed)),
          n_physical_blocks=int(passed.spatial_block.nunique()),
          n_first_GIFT_source_signatures=int(passed.first_source_signature.nunique()),
          n_islands_with_10_or_more_admissible_pairs=int(
              g.n_matched_unique_pairs_meeting_original_phylo_bound.ge(10).sum()),
          n_islands_with_15_or_more_admissible_pairs=int(
              g.n_matched_unique_pairs_meeting_original_phylo_bound.ge(15).sum()),
          median_max_admissible_true_tip_pairs=float(
              g.n_matched_unique_pairs_meeting_original_phylo_bound.median()),
          reportable_under_20_islands_10_blocks_10_sources=bool(ok))
        if ok:
            for j,key in enumerate((
                 "phylo_bounded_delta_k_adjusted_rare_minus_common",
                 "same_family_unbounded_delta_k_adjusted_rare_minus_common",
                 "bounded_minus_family_delta_k_adjusted_rare_minus_common",
                 "phylo_bounded_delta_rare_coverage",
                 "same_family_unbounded_delta_rare_coverage",
                 "matched_tree_distance","family_control_tree_distance",
                 "matched_vs_same_family_control_distance",
                 "n_distinct_families_in_selected_true_tip_pairs"
            )):
                for tier,cluster in (
                    ("block","spatial_block"),
                    ("source","first_source_signature")):
                    mean,low,high=_resample_group(
                        passed,key,cluster,seed+31*j+int(cutoff_from_name(criterion))
                             +(tier=="source"),draws=bootstrap)
                    report[tier+"_mean_"+key]=mean
                    report[tier+"_ci025_"+key]=low
                    report[tier+"_ci975_"+key]=high
        regions.append(report)
    # POST-HOC matched-threshold sensitivity on the EXACT SAME original
    # physical island IDs, same n=5 original source tip pair protocol. This
    # removes geographic-cohort attrition between 50 and 100 cutoffs, but the
    # actual focal plant species and donor pairs CAN still differ between
    # these independent criteria. Never call this a measured phylogenetic
    # evolutionary gradient or causal branch-length effect.
    threshold_paired=[]
    for (scope,regime),q in raw.groupby(["scope","regime"]):
        restrictive=q.loc[
            q.criterion.eq("same_family_d_le50")&
            q.strict_n5_phylogenetic_source_pairs_available].set_index(
                "original_island_id")
        broad=q.loc[
            q.criterion.eq("same_family_d_le100")&
            q.strict_n5_phylogenetic_source_pairs_available].set_index(
                "original_island_id")
        same=restrictive.index.intersection(broad.index)
        if len(same)<10:
            continue
        a=restrictive.loc[same]
        b=broad.loc[same]
        if (a.spatial_block!=b.spatial_block).any() or (
              a.first_source_signature!=b.first_source_signature).any():
            raise ValueError("Paired branch-length thresholds changed original island/source geography")
        good=len(same)>=MIN_REGIONAL_ISLANDS and (
            a.spatial_block.nunique()>=MIN_REGIONAL_BLOCKS) and (
            a.first_source_signature.nunique()>=MIN_REGIONAL_SOURCES)
        item=dict(
            scope=scope,regime=regime,
            n_exact_same_islands_with_50_and_100_source_tip_pairs=len(same),
            n_physical_blocks=int(a.spatial_block.nunique()),
            n_first_nearest_source_signatures=int(a.first_source_signature.nunique()),
            result_reportable_under_20_islands_10_blocks_10_sources=bool(good),
            same_island_but_not_necessarily_same_plant_tip_subset=True,
            cutoff_effect_comparison_is_posthoc_not_causal=True)
        if good:
            df=a[["spatial_block","first_source_signature"]].copy()
            for j,key in enumerate((
                "phylo_bounded_delta_k_adjusted_rare_minus_common",
                "same_family_unbounded_delta_k_adjusted_rare_minus_common",
                "bounded_minus_family_delta_k_adjusted_rare_minus_common",
                "phylo_bounded_delta_rare_coverage",
                "matched_tree_distance")):
                change="d50_minus_d100_"+key
                df[change]=(a[key]-b[key]).to_numpy(float)
                for tier,cluster in (
                    ("block","spatial_block"),
                    ("source","first_source_signature")):
                    mean,lo,hi=_resample_group(df,change,cluster,
                        seed+1200+21*j+(tier=="source"),draws=bootstrap)
                    item[tier+"_mean_"+change]=mean
                    item[tier+"_ci025_"+change]=lo
                    item[tier+"_ci975_"+change]=hi
        threshold_paired.append(item)
    # Further restriction: same IDENTICAL original island plant species
    # under two non-overlapping source distance bands, NOT merely same
    # island as in paired 50 versus 100 above. Do not output a numerical
    # effect if the true original donor/source support is inadequate.
    same_focal=[]
    fixed=pd.DataFrame.from_records(fixed_local_bands)
    for (scope,regime),q in fixed.groupby(["scope","regime"]):
        good=q.loc[q.qualified_exact_focal_taxa_for_two_distance_bands]
        supported=(len(good)>=MIN_REGIONAL_ISLANDS and
                   good.spatial_block.nunique()>=MIN_REGIONAL_BLOCKS and
                   good.first_source_signature.nunique()>=MIN_REGIONAL_SOURCES)
        record=dict(scope=scope,regime=regime,
            n_eligible_source_plant_record_islands=len(q),
            n_islands_with_identical_focal_five_species_two_source_distance_bands=len(good),
            n_physical_blocks=int(good.spatial_block.nunique()),
            n_first_GIFT_source_signatures=int(good.first_source_signature.nunique()),
            median_n_possible_focal_tree_tips_two_source_bands=float(
                q.n_original_local_true_tip_species_with_two_donor_distance_bands.median()),
            original_n5_same_focal_source_bands_adequately_supported=bool(supported))
        if supported:
            for j,key in enumerate((
                "local_minus_near_k_adjusted_rare_common",
                "local_minus_mid_k_adjusted_rare_common",
                "near_minus_mid_k_adjusted_rare_common",
                "near_minus_mid_rare_coverage",
                "average_near_original_tree_distance",
                "average_mid_original_tree_distance")):
                for level,cluster in (
                    ("block","spatial_block"),("source","first_source_signature")):
                    mean,lo,hi=_resample_group(good,key,cluster,
                        seed+2400+17*j+(level=="source"),draws=bootstrap)
                    record[level+"_mean_"+key]=mean
                    record[level+"_ci025_"+key]=lo
                    record[level+"_ci975_"+key]=hi
        same_focal.append(record)
    return dict(
       status="EXPLORATORY_TRUE_BRANCH_BOUNDED_WITHIN_FAMILY_PHYLOGENETIC_SOURCE",
       tree_version=tree_version,
       original_reference_Direct_species=len(triple),
       n_original_source_eligible_full_islands=len(islands),
       n_original_exact_backbone_tips=len(strict_tip_species),
       region_phylogenetic_close_family_source_evidence=regions,
       paired_same_island_50_vs_100_tree_cutoff_sensitivity=threshold_paired,
       identical_focal_five_species_near_vs_mid_source_distance_bands=same_focal,
       design={
         "all_graft_and_taxonomy_inserted_species_excluded":True,
         "families_are_independent_exact_colxr_source_candidates":True,
         "within_family_5_distinct_original_tree_tip_matched_source_species":True,
         "n_required_original_true_species_per_comparison":5,
         "cutoffs_original_branch_length_50_100_are_posthoc":True,
         "same_genus_within_100_support_is_a_distinct_more_restrictive_arm":True,
         "source_selection_does_not_examine_flower_trait_phenotypes":True,
         "source_all_identical_local_plant_species_names_removed":True,
         "comparison_to_unrestricted_phylo_distance_same_family_donor":True,
         "fixed_profile_count_k_null_reference_draws":reference_draws,
         "original_global_Direct_flower_reference_100_triple_profiles":True,
         "numerical_original_phylo_pairing_repetitions":repeat,
         "bootstrap_resampling_geography_and_nearest_GIFT_source":bootstrap,
         "near_branch_length_not_verified_island_ancestral_species":True,
         "same_island_cutoff_pairing_keeps_geography_but_not_identical_plant_species":True,
         "new_strict_nonoverlapping_d50_vs_d50to100_identical_five_focal_plant_species":True,
         "new_strict_nonoverlapping_source_distance_cutoff_is_posthoc":True,
       },
       no_biological_pollen_delivery_or_effective_visitation_measured=True,
       not_phylogeny_or_evolution_causal_island_phenotype_change=True,
       no_original_plant_family_genus_island_source_or_tree_ids_exported=True
    )


def cutoff_from_name(s:str)->int:
    if "le50" in s:return 50
    if "le100" in s:return 100
    raise ValueError("Unknown original dated tree threshold")
