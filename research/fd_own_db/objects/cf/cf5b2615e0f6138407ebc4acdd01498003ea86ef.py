"""Family-conditioned island floral categorical rare-profile support.

Original 2234 Direct species; 100 observed triple categorical profiles.
CoL XR alpha2 strict *candidate* families only (exact accepted source name,
exact source family, exact CoL candidate family, same genus, original source
SHA). Family is NOT a phylogenetic tree or independently reviewed taxonomy.

Two independently predeclared source-disjoint controls:
 (A) identical original 5 family IDs, one different source species on
     each side, allowing the genus to coincide or differ;
 (B) same 5 family IDs, forcing different GENUS identities within
     each family pair (different accepted species always).
8-family sensitivities are fixed before inspecting the response.

A separate exact two-player Shapley bookkeeping decomposition assigns the
original n30 observed rare-category neighborhood coverage to common vs
one-sided family groups. This is NOT a causal decomposition or an independent
phylogenetic randomization. All aggregations are regional; no original
island/site/taxon/family/source IDs can be exported.
"""
from __future__ import annotations
from collections import defaultdict
import numpy as np
import pandas as pd

from .discrete_functional_volume import TRAITS,_resample_group
from .rare_profile_decomposition import RarityReference,genus,_seedstr
from .genus_shapley_functional_volume import two_group_shapley
from .colxr_family_opportunity import group_by_family

PRIMARY=5
SENSITIVITY=8
METRICS=("rare_coverage","common_coverage","rare_minus_common",
         "k_adjusted_rare_minus_common","freq_k_adjusted_rare_minus_common")


def _compiled(reference:RarityReference,triple:pd.DataFrame):
    result={}
    for name,row in zip(triple.index,triple[list(TRAITS)].to_numpy(dtype=str)):
        profile=tuple(row)
        result[str(name)]=(profile,reference.volume.masks[profile][1])
    if len(result)!=len(triple):
        raise ValueError("Original accepted plant source names not unique")
    return result


def _score(ref:RarityReference, names:list[str],lookup:dict)->np.ndarray:
    covered=0
    seen=set()
    for name in names:
        if name not in lookup:
            raise ValueError("Source family plant name not in Direct phenotype data")
        phenotype,bitmap=lookup[name]
        seen.add(phenotype)
        covered |= bitmap
    k=len(seen)
    if k not in ref.expected_uniform:
        raise ValueError("Original profile count k outside frozen null")
    bins=np.asarray([
        (covered & ref.bins[p]).bit_count()/ref.n_bins[p]
        for p in ("rare","middle","common")
    ],dtype=float)
    null=ref.expected_uniform[k]
    weighted=ref.expected_speciesfreq[k]
    return np.asarray((bins[0],bins[2],bins[0]-bins[2],
        bins[0]-bins[2]-null[0]+null[2],
        bins[0]-bins[2]-weighted[0]+weighted[2]),dtype=float)


def _mask(names:list[str],lookup:dict)->int:
    result=0
    for name in names:
        result |= lookup[name][1]
    return result


def _families_with_different_genus_pairs(a:dict,b:dict):
    eligible={}
    for family in sorted(set(a)&set(b)):
        options=[]
        for island_plant in a[family]:
            alternative=[donor_plant for donor_plant in b[family]
                         if genus(donor_plant)!=genus(island_plant)]
            if alternative:
                options.append((island_plant,tuple(alternative)))
        if options:eligible[family]=options
    return eligible


def matched_family_source_effect(ref:RarityReference, lookup:dict,
                                 island_species:set[str], donor_species:set[str], *,
                                 family_map:dict[str,str], n:int=5,
                                 repeats:int=40, seed:int=20261113,
                                 must_change_genus:bool=False)->dict[str,float]:
    if n not in (PRIMARY,SENSITIVITY) or repeats<20:
        raise ValueError("Frozen original n5/n8 family comparison contract changed")
    if island_species&donor_species:
        raise ValueError("Source donor must exclude every original island accepted species")
    a=group_by_family(island_species,family_map)
    b=group_by_family(donor_species,family_map)
    eligible=sorted(set(a)&set(b))
    pairs=_families_with_different_genus_pairs(a,b) if must_change_genus else {}
    if must_change_genus:
        eligible=sorted(pairs)
    if len(eligible)<n:
        raise ValueError("Insufficient supported identical-family different-species source pairs")
    rng=np.random.default_rng(seed)
    out=[]
    for _ in range(repeats):
        chosen=sorted(rng.choice(eligible,size=n,replace=False).tolist())
        if must_change_genus:
            ia=[]
            sb=[]
            for family in chosen:
                options=pairs[family]
                local,available=options[int(rng.integers(len(options)))]
                source=str(rng.choice(available))
                if genus(local)==genus(source):
                    raise AssertionError("Alternative-genus family pair not different")
                ia.append(local)
                sb.append(source)
        else:
            ia=[str(rng.choice(a[family])) for family in chosen]
            sb=[str(rng.choice(b[family])) for family in chosen]
        if (len(set(ia))!=n or len(set(sb))!=n or
                set(ia)&set(sb) or
                [family_map[z] for z in ia]!=[family_map[z] for z in sb]):
            raise AssertionError("Family exact identity or distinct source acceptance failed")
        out.append(_score(ref,ia,lookup)-_score(ref,sb,lookup))
    mean=np.asarray(out,dtype=float).mean(axis=0)
    return dict(zip(METRICS,map(float,mean)))


def independent_family_richness_effect(
    ref:RarityReference, lookup:dict,
    island_species:set[str], donor_species:set[str],*,
    family_map:dict[str,str],n:int=5,
    repeats:int=40,seed:int=20261113
)->dict[str,float]:
    """POST-HOC counter-estimand: equal number of families, not same IDs.

    Explicitly distinct from original same-family identity matching;
    its paired island-level comparison is an exploratory sensitivity, not
    an independent confirmatory permutation test of family selection.
    """
    if n not in (PRIMARY,SENSITIVITY) or repeats<20:
        raise ValueError("Independent original family richness sample contract")
    if island_species&donor_species:
        raise ValueError("Source and island accepted plant species must be disjoint")
    a=group_by_family(island_species,family_map)
    b=group_by_family(donor_species,family_map)
    if len(a)<n or len(b)<n:
        raise ValueError("Insufficient n different botanical families to draw")
    rng=np.random.default_rng(seed)
    draw=[]
    for _ in range(repeats):
        ia=sorted(rng.choice(sorted(a),size=n,replace=False).tolist())
        sb=sorted(rng.choice(sorted(b),size=n,replace=False).tolist())
        names_a=[str(rng.choice(a[fam])) for fam in ia]
        names_b=[str(rng.choice(b[fam])) for fam in sb]
        if len(set(names_a))!=n or len(set(names_b))!=n or set(names_a)&set(names_b):
            raise AssertionError("Family count-only source species sample invalid")
        draw.append(_score(ref,names_a,lookup)-_score(ref,names_b,lookup))
    return dict(zip(METRICS,map(float,np.asarray(draw,float).mean(axis=0))))


def _shapley_family_scores(ref:RarityReference,lookup:dict,
                           local:set[str],source:set[str],families:dict[str,str],
                           *,repeats:int=40,seed:int=20261113):
    if len(local)<30 or len(source)<30:
        raise ValueError("Thirty Direct family-qualified source taxa needed")
    a=group_by_family(local,families)
    b=group_by_family(source,families)
    common=set(a)&set(b)
    rng=np.random.default_rng(seed)
    scores=[]
    for _ in range(repeats):
        lhs=rng.choice(sorted(local),size=30,replace=False).tolist()
        rhs=rng.choice(sorted(source),size=30,replace=False).tolist()
        sides=[]
        for sample in (lhs,rhs):
            shared=[p for p in sample if families[p] in common]
            only=[p for p in sample if families[p] not in common]
            sides.append(two_group_shapley(ref,_mask(shared,lookup),_mask(only,lookup)))
        values={
            key:sides[0][key]-sides[1][key]
            for key in sides[0]
        }
        if abs(values["rare_total"]-values["rare_shared"]-
               values["rare_unshared"])>1e-12:
            raise AssertionError("Family Shapley rare support missing or nonadditive")
        scores.append(values)
    result=pd.DataFrame(scores).mean().to_dict()
    result["n_shared_candidate_families"]=len(common)
    result["n_local_candidate_families"]=len(a)
    result["n_source_candidate_families"]=len(b)
    return {str(k):float(v) for k,v in result.items()}


def original_family_rare_support(
    triple:pd.DataFrame, original_flora:pd.DataFrame,geo:pd.DataFrame,
    assignments:pd.DataFrame,mainland:pd.DataFrame,family_map:dict[str,str], *,
    rarefactions:int=40,reference_draws:int=399,
    bootstraps:int=1999,seed:int=20261113
)->dict:
    if rarefactions<20 or reference_draws<199 or bootstraps<999:
        raise ValueError("Original family morphology numerical support weakened")
    if len(triple)!=2234 or len(family_map)<1900:
        raise ValueError("Original 2234 Direct and sufficiently qualified CoL XR families required")
    for frame,required in (
        (original_flora,{"island_id","accepted_species","origin_status"}),
        (geo,{"island_id","analysis_regime","spatial_block"}),
        (assignments,{"island_id","source_mode","source_rank","entity_ID"}),
        (mainland,{"entity_ID","work_species"})):
        if not required.issubset(frame):
            raise ValueError("Original source botanical family comparator missing")
    flora=original_flora.loc[
       original_flora.accepted_species.isin(triple.index)].copy()
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original source accepted island species duplicates")
    if not flora.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Original source native status unexpectedly changed")
    count=flora.groupby("island_id").accepted_species.nunique()
    eligible=set(count[count.ge(30)].index)
    if len(eligible)!=928:
        raise ValueError("Original broad Direct 928-island floras missing")
    flora=flora.loc[flora.island_id.isin(eligible)]
    geo=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates()
    if geo.island_id.duplicated().any() or not eligible.issubset(set(geo.island_id)):
        raise ValueError("Original geographic block contract not satisfied")
    location=geo.set_index("island_id")
    assignment=assignments.loc[
        assignments.island_id.isin(eligible)&
        assignments.source_mode.eq("geo_k5")&
        assignments.source_rank.between(1,5),
        ["island_id","source_rank","entity_ID"]]
    if assignment.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Conflicting original candidate source slot")
    donors=assignment.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(donors.get(i,()))!=5 or len(set(donors[i]))!=5 for i in eligible):
        raise ValueError("Five original GIFT source assignments missing")
    source=mainland.loc[mainland.work_species.isin(triple.index),
        ["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    byisland={i:q for i,q in flora.groupby("island_id")}
    ref=RarityReference.from_original(triple,null_draws=reference_draws,seed=seed)
    lookup=_compiled(ref,triple)
    # Diagnose whether strict CoL family matching selectively loses rare
    # original flower categories (only 44/2234 species lack exact source
    # family candidates, but those could carry singleton rare profiles).
    mapped_names=set(family_map)&set(triple.index)
    mapped_profile_counts={}
    for species in mapped_names:
        profile=lookup[species][0]
        mapped_profile_counts[profile]=mapped_profile_counts.get(profile,0)+1
    lost_bin={}
    for label,which in (("rare",lambda x:x<=ref.low_quartile_species_count),
                         ("middle",lambda x:ref.low_quartile_species_count<x<ref.high_quartile_species_count),
                         ("common",lambda x:x>=ref.high_quartile_species_count)):
        profiles=[p for p,n in zip(ref.volume.profiles,ref.profile_n_species) if which(int(n))]
        lost_bin[label]=dict(
            n_original_reference_profiles=len(profiles),
            n_original_reference_profiles_with_zero_mapped_carriers=sum(
                mapped_profile_counts.get(p,0)==0 for p in profiles),
            n_original_species_carriers=int(sum(
                int(n) for p,n in zip(ref.volume.profiles,ref.profile_n_species)
                if which(int(n)))),
            n_original_carriers_lacking_strict_family=sum(
                int(n)-int(mapped_profile_counts.get(p,0))
                for p,n in zip(ref.volume.profiles,ref.profile_n_species)
                if which(int(n))),
        )
    rows=[]
    shapley=[]
    for island in sorted(eligible):
        group=byisland[island]
        allnames=set(group.accepted_species)
        raw_source=set().union(*(source.get(j,set()) for j in donors[island]))
        donor=(raw_source & set(triple.index))-allnames
        donor={z for z in donor if z in family_map}
        family_source=group_by_family(donor,family_map)
        for scope,local in (
            ("all_recorded",allnames),
            ("source_native",set(group.loc[
                group.origin_status.eq("native"),"accepted_species"]))):
            if len(local)<30: continue
            local={z for z in local if z in family_map}
            family_island=group_by_family(local,family_map)
            shared=set(family_island)&set(family_source)
            changed=_families_with_different_genus_pairs(
                family_island,family_source)
            info=dict(
                island_id=island,scope=scope,
                regime=str(location.at[island,"analysis_regime"]),
                spatial_block=str(location.at[island,"spatial_block"]),
                source_signature=donors[island][0],
                n_original_local_family_qualified=len(local),
                n_source_family_qualified=len(donor),
                n_shared_families=len(shared),
                n_family_different_genus_options=len(changed))
            if len(local)>=30 and len(donor)>=30:
                sh=_shapley_family_scores(ref,lookup,local,donor,family_map,
                    repeats=rarefactions,
                    seed=_seedstr(seed,island,scope,"family_shapley"))
                shapley.append({**info,**sh})
            for n in (PRIMARY,SENSITIVITY):
                for kind in ("family_richness_only",
                             "same_family","same_family_different_genus"):
                    present=(len(family_island)>=n and len(family_source)>=n
                             if kind=="family_richness_only"
                             else len(shared)>=n if kind=="same_family"
                             else len(changed)>=n)
                    entry=dict(info,condition=kind,n_shared_family_units=n,
                               eligible_original_family_matched=bool(present))
                    if present:
                        if kind=="family_richness_only":
                            value=independent_family_richness_effect(
                                ref,lookup,local,donor,family_map=family_map,n=n,
                                repeats=rarefactions,
                                seed=_seedstr(seed,island,scope,"family_count",n))
                        else:
                            value=matched_family_source_effect(
                                ref,lookup,local,donor,
                                family_map=family_map,n=n,
                                repeats=rarefactions,
                                seed=_seedstr(seed,island,scope,"fam",kind,n),
                                must_change_genus=(
                                   kind=="same_family_different_genus"))
                        entry.update(value)
                    rows.append(entry)
    private=pd.DataFrame(rows)
    shapley_private=pd.DataFrame(shapley)
    match_panels=[]
    for (scope,regime,n,condition),q in private.groupby([
          "scope","regime","n_shared_family_units","condition"]):
        g=q.loc[q.eligible_original_family_matched]
        eligible=len(g)>=10 and g.spatial_block.nunique()>=5 and g.source_signature.nunique()>=5
        result={
            "scope":scope,"regime":regime,"n_exact_shared_families_per_side":int(n),
            "condition":condition,
            "n_original_scope_islands":len(q),
            "n_family_matched_islands":len(g),
            "n_blocks":g.spatial_block.nunique(),
            "n_first_GIFT_source_signatures":g.source_signature.nunique(),
            "reportable_under_10island_5block_5source_gate":bool(eligible)}
        if eligible:
            for j,key in enumerate(METRICS):
                result["median_"+key]=float(g[key].median())
                for level,cluster in (("block","spatial_block"),
                                       ("source","source_signature")):
                    mean,lo,hi=_resample_group(g,key,cluster,
                      seed+n*100+30*j+(level=="source"),draws=bootstraps)
                    result[level+"_mean_"+key]=mean
                    result[level+"_ci025_"+key]=lo
                    result[level+"_ci975_"+key]=hi
        match_panels.append(result)
    # A DIRECT paired, same-island difference of estimands:
    # requiring same FAMILY IDs versus drawing just the same number
    # of family units from either side. This is POST-HOC diagnostic
    # after the original family-matched result was inspected; therefore
    # conditional intervals, not confirmatory family-selection p values.
    paired_panels=[]
    for (scope,regime,n),q in private.groupby([
        "scope","regime","n_shared_family_units"]):
        allowed=q.loc[q.eligible_original_family_matched]
        same=allowed.loc[allowed.condition.eq("same_family")].set_index("island_id")
        free=allowed.loc[allowed.condition.eq("family_richness_only")].set_index("island_id")
        common=same.index.intersection(free.index)
        if len(common)<10:continue
        a=same.loc[common]
        b=free.loc[common]
        if a.spatial_block.nunique()<5 or a.source_signature.nunique()<5:
            continue
        block=a[["spatial_block","source_signature"]].copy()
        result=dict(scope=scope,regime=regime,
            families_per_side=n,
            n_same_island_comparable=len(common),
            n_blocks=a.spatial_block.nunique(),
            n_first_GIFT_source_signatures=a.source_signature.nunique(),
            posthoc_difference_of_family_identity_and_family_count_estimands=True,
            no_causal_family_effect_identified=True)
        for j,metric in enumerate(METRICS):
            key="id_matched_minus_count_only_"+metric
            block[key]=(a[metric]-b[metric]).to_numpy(float)
            for label,cluster in (("block","spatial_block"),
                                  ("source","source_signature")):
                mean,lo,hi=_resample_group(block,key,cluster,
                    seed+910+n*10+j+(label=="source"),draws=bootstraps)
                result[label+"_mean_"+key]=mean
                result[label+"_ci025_"+key]=lo
                result[label+"_ci975_"+key]=hi
        paired_panels.append(result)
    sh_panels=[]
    if not shapley_private.empty:
        for (scope,regime),g in shapley_private.groupby(["scope","regime"]):
            sufficient=len(g)>=10 and g.spatial_block.nunique()>=5 and g.source_signature.nunique()>=5
            if not sufficient:continue
            rec=dict(scope=scope,regime=regime,n_islands=len(g),
                n_spatial_blocks=g.spatial_block.nunique(),
                n_first_GIFT_source_signatures=g.source_signature.nunique(),
                median_island_family_richness=float(g.n_local_candidate_families.median()),
                median_source_family_richness=float(g.n_source_candidate_families.median()),
                median_shared_source_families=float(g.n_shared_candidate_families.median()))
            for j,key in enumerate(("rare_total","rare_shared","rare_unshared",
                                    "gap_total","gap_shared","gap_unshared")):
                rec["median_"+key]=float(g[key].median())
                for label,cluster in (("block","spatial_block"),
                                      ("source","source_signature")):
                    mean,lo,hi=_resample_group(g,key,cluster,seed+600+j*13+
                        (label=="source"),draws=bootstraps)
                    rec[label+"_mean_"+key]=mean
                    rec[label+"_ci025_"+key]=lo
                    rec[label+"_ci975_"+key]=hi
            sh_panels.append(rec)
    return dict(
       status="COL_XR_CANDIDATE_FAMILY_FD_RARITY_SUPPORT_V1",
       n_original_Direct_species=len(triple),
       n_strict_candidate_family_species=len(family_map),
       original_source_family_mapping_ascertainment_by_global_profile_rarity=lost_bin,
       n_original_eligible_island_flora=928,
       n_original_flora_families_source_shapley_eligible=len(shapley),
       family_shapley_rare_category_macro_panels=sh_panels,
       same_family_island_mainland_matched_panels=match_panels,
       posthoc_paired_exact_family_identity_vs_family_count_panels=paired_panels,
       inference={
          "family_labels_are_unreviewed_exact_accepted_COL_XR_candidates":True,
          "family_ranks_not_a_resolved_phylogenetic_tree":True,
          "same_family_genus_may_differ_by_source_and_island":True,
          "family_count_only_keeps_n_family_not_family_identity":True,
          "paired_family_id_vs_count_sensitivity_is_posthoc":True,
          "different_genus_condition_restricts_to_different_genus_in_same_family":True,
          "original_island_species_excluded_from_all_nearest_GIFT_native_source_flora":True,
          "strict_plant_species_source_family_mapping_only":True,
          "equal_30_source_island_species_for_shapley":True,
          "equal_5_family_one_plant_each_source_and_island_primary":True,
          "equal_8_family_one_plant_each_source_and_island_sensitivity":True,
          "reference_100_original_observed_triple_profiles":True,
          "reference_rare_profile_not_necessarily_biological_pollination_function":True,
          "no_pollinator_loss_historical_evolution_or_fitness_tested":True,
          "not_causal_selection_or_phylogenetic_effect":True,
          "bootstrap_whole_original_source_or_spatial_blocks":bootstraps,
          "numerical_rarefaction_repeats":rarefactions
       },
       no_original_species_family_source_island_or_coordinates_exported=True
    )
