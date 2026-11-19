"""Rights-preserving STRICT taxon-name-out-of-mainland floral source predictor.

Retains 4,505 source botanical island floras and original Direct colour+form
states. Prediction treats mainland GIFT lists as *candidate species opportunity*,
never measured independent mainland phenotypes or proven ancestry. All heldout
accepted names are removed from EVERY GIFT mainland entity before model training.

Three posthoc overlapping prospective test cohorts are reported separately,
and a fourth all-island exclusion stress checks whether an actually unobserved
species donor remains. Do not use observed island genus composition as a
prospective colonization predictor; it is a conditional oracle.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from fd_island_functional_beta import strict_traits, membership, geographical_distance

SHRINKAGE=8.0
MIN_TARGETS=30
MIN_BLOCKS=10


def genus(name):
    terms=str(name).strip().split()
    return terms[0] if len(terms)>1 else "UNRESOLVED:"+str(name)


def donor_probability(counts, prior:np.ndarray, shrink:float=SHRINKAGE):
    """Strictly positive shrink-to-prior and genus/focal source hierarchy."""
    if shrink<=0 or prior.ndim!=1 or not np.isclose(prior.sum(),1) or np.any(prior<=0):
        raise ValueError("Finite positive source prior required")
    x=np.asarray(counts,dtype=float)
    if x.shape!=prior.shape or np.any(x<0) or not np.isfinite(x).all():
        raise ValueError("Invalid original genus source evidence")
    p=(x+shrink*prior)/(x.sum()+shrink)
    if not np.isfinite(p).all() or not np.isclose(p.sum(),1):
        raise AssertionError("Invalid donor probability")
    return p


def source_enriched_genus_probability(global_counts, local_counts, prior,
                                       shrink:float=SHRINKAGE):
    """Zero-sum source-reweighting; NEVER double-count source donors.

    Local source donor names are a SUBSET of unique global donor names. Add the
    observed local composition deviation to global counts, preserving the
    original genus donor total exactly. When the source covers all/none of
    the donor genus, the candidate equals the global-genus baseline.
    This is a fixed descriptive recalibration, not an ancestral model.
    """
    g=np.asarray(global_counts,dtype=float)
    l=np.asarray(local_counts,dtype=float)
    if (g.shape!=l.shape or np.any(g<0) or np.any(l<0)
            or np.any(l>g+1e-9) or not np.isfinite(g+l).all()):
        raise ValueError("Local source names must be a subset of global donors")
    n=float(g.sum())
    if n==0:
        return donor_probability(g,prior,shrink)
    alpha=float(l.sum()/n)
    corrected=(1-alpha)*g+l
    if not np.isclose(corrected.sum(),n,atol=1e-9):
        raise AssertionError("Source enrichment changed genus donor sample size")
    return donor_probability(corrected,prior,shrink)


def score(p,truth):
    """Exact categorical full colour×form log-loss and Brier (lower is better)."""
    if not np.isfinite(p).all() or np.any(p<=0) or not np.isclose(p.sum(),1):
        raise ValueError("Invalid prediction simplex")
    if truth<0 or truth>=len(p):raise ValueError("Unknown measured joint phenotype")
    return float(-np.log(p[truth])),float(1+np.sum(p*p)-2*p[truth])


def exact_join_and_source(mainland:pd.DataFrame, assignment:pd.DataFrame,
                          labelled:dict):
    if not {"entity_ID","work_species"}<=set(mainland):
        raise ValueError("Original mainland GIFT entity×name required")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assignment):
        raise ValueError("Original source rank required")
    if mainland[["entity_ID","work_species"]].isna().any().any():
        raise ValueError("Unresolved original mainland identity")
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("Duplicate GIFT native source occurrence")
    if not labelled:raise ValueError("Missing actual Direct profile measurements")
    assigned=assignment.loc[assignment.source_mode.eq("geo_k5")].copy()
    assigned["source_rank"]=pd.to_numeric(assigned.source_rank,errors="raise")
    assigned=assigned.loc[assigned.source_rank.between(1,5)]
    if assigned.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Geographic source rank repeated")
    near=assigned.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(lambda x:tuple(x)).to_dict()
    donor=mainland.loc[mainland.work_species.isin(labelled)].copy()
    donor["work_species"]=donor.work_species.astype(str)
    # A donor name is one species regardless of the number of source regions.
    return donor,near


def predict_one_cohort(cohort:str,island_list:list[str],
                       full_rosters:dict,native_rosters:dict,
                       labelled:dict,donors:pd.DataFrame,
                       nearest:dict,geo:pd.DataFrame,
                       *,test_scope:str) -> dict:
    if test_scope not in ("all_recorded","original_native"):
        raise ValueError("Unknown source-native status scope")
    if not island_list:
        return dict(cohort=cohort,status="no_evaluable_islands")
    measured=full_rosters if test_scope=="all_recorded" else native_rosters
    heldout_all_names=set().union(*(full_rosters[i] for i in island_list))
    if not heldout_all_names:
        raise ValueError("Missing target original full-flora species")
    # This is the CORRECT, stronger holdout over all focal flora names and all
    # source entities, not only the labelled tested species for that trait.
    remaining=donors.loc[~donors.work_species.isin(heldout_all_names)].copy()
    if set(remaining.work_species)&heldout_all_names:
        raise AssertionError("TEST-name leakage in original candidate GIFT sources")
    donor_unique=sorted(set(remaining.work_species))
    x=sorted(set(labelled.values()))
    cats={t:i for i,t in enumerate(x)}
    k=len(x)
    category={sp:cats[value] for sp,value in labelled.items()}
    all_counts=Counter(category[s] for s in donor_unique)
    gcounts=defaultdict(Counter)
    for sp in donor_unique:
        gcounts[genus(sp)][category[sp]]+=1
    n_train=len(donor_unique)
    all_registry=set(donors.work_species)
    train_genus=set(gcounts)
    # Important denominator: distinct actually directly labelled source donors,
    # not repeated entries in 748 GIFT entity lists.
    profile_variable_genera=sum(len(z)>=2 for z in gcounts.values())
    near_donors=remaining.groupby("entity_ID").work_species.agg(set).to_dict()
    geo=geo.set_index("island_id")
    target_total=0
    n_names_in_mainland_raw=0
    n_target_predicted=0
    n_local_supported=0
    n_global_genus_supported=0
    n_local_genus_supported=0
    n_islands_local_supported=0
    n_islands_any_source=0
    n_islands_missing_source=0
    n_islands_without_local_genus=0
    scored_islands=[]
    if n_train:
        prior=(np.array([all_counts.get(i,0) for i in range(k)],float)+.5)
        prior/=prior.sum()
    for island in island_list:
        test=sorted(measured.get(island,set())&set(labelled))
        if not test:
            continue
        n=len(test)
        target_total+=n
        n_names_in_mainland_raw+=sum(s in all_registry for s in test)
        names=nearest.get(island,())
        if not names:
            n_islands_missing_source+=1
        else:
            n_islands_any_source+=1
        source_names=set().union(*(near_donors.get(e,set()) for e in names))
        source_genus=defaultdict(Counter)
        for s in source_names:
            source_genus[genus(s)][category[s]]+=1
        n_local_ref=len(source_names)
        if source_names:n_islands_local_supported+=1
        known_n=sum(genus(s) in train_genus for s in test)
        local_known_n=sum(genus(s) in source_genus for s in test)
        n_global_genus_supported+=known_n
        n_local_genus_supported+=local_known_n
        if not local_known_n:n_islands_without_local_genus+=1
        n_local_supported+=n_local_ref
        if not n_train:
            continue
        local_prior=np.array([sum(v.get(i,0) for v in source_genus.values())
                              for i in range(k)],float)
        local_prior=donor_probability(local_prior,prior)
        losses=np.zeros((4,2),float)
        for sp in test:
            label=category[sp]
            g=genus(sp)
            gt=donor_probability(np.array([gcounts[g].get(i,0) for i in range(k)]),
                                 prior)
            lt=source_enriched_genus_probability(
                np.array([gcounts[g].get(i,0) for i in range(k)]),
                np.array([source_genus[g].get(i,0) for i in range(k)]),
                prior)
            models=(prior,gt,local_prior,lt)
            for j,p in enumerate(models):
                losses[j]+=score(p,label)
        losses/=n
        scored_islands.append(dict(
            island=island,block=geo.at[island,"spatial_block"] if island in geo.index else None,
            n_test=n,n_local_donor=n_local_ref,
            n_global_genus=known_n,n_local_genus=local_known_n,
            losses=losses))
        n_target_predicted+=n
    base=dict(
        cohort=cohort,test_scope=test_scope,
        n_original_focal_islands=len(island_list),
        n_original_focal_full_flora_names=len(heldout_all_names),
        n_test_direct_colour_form_species_carriers=target_total,
        n_test_names_present_in_unheldout_mainland_source=n_names_in_mainland_raw,
        n_unique_source_direct_pair_donor_species_before_holdout=donors.work_species.nunique(),
        n_unique_source_direct_pair_donor_species_after_exact_name_holdout=n_train,
        n_direct_pair_donor_genera_after_holdout=len(gcounts),
        n_donor_genera_with_multiple_joint_profiles=profile_variable_genera,
        n_source_entities_with_some_direct_pair_donor=remaining.entity_ID.nunique(),
        n_islands_with_any_five_nearest_source_assignment=n_islands_any_source,
        n_islands_missing_five_nearest_source_assignment=n_islands_missing_source,
        n_islands_with_nonzero_locally_measured_direct_donors=n_islands_local_supported,
        n_islands_without_source_local_congeneric_donors=n_islands_without_local_genus,
        ratio_test_carriers_with_any_global_same_genus_donor=(
            n_global_genus_supported/target_total if target_total else None),
        ratio_test_carriers_with_same_genus_local_five_source_donor=(
            n_local_genus_supported/target_total if target_total else None),
        mean_local_five_source_trait_known_donor_species=(
            n_local_supported/len(island_list)),
    )
    if not n_train:
        return dict(**base,status="NOT_IDENTIFIABLE_NO_DIRECT_PAIR_GIFT_DONORS_AFTER_HOLDOUT")
    if len(scored_islands)<MIN_TARGETS:
        return dict(**base,status="NOT_IDENTIFIABLE_TOO_FEW_SCORED_ISLANDS")
    # Island and geographic block equal average. No dyad count contributes.
    agg=np.asarray([r["losses"] for r in scored_islands])
    valid=[r for r in scored_islands
           if r["block"] is not None and not pd.isna(r["block"])]
    n_blocks=len({r["block"] for r in valid})
    model_names=["prior","global_genus","source_mixture","local_genus"]
    values=dict()
    for idx,name in enumerate(model_names):
        values[name+"_mean_island_log_loss"]=float(agg[:,idx,0].mean())
        values[name+"_mean_island_brier"]=float(agg[:,idx,1].mean())
        if n_blocks>=MIN_BLOCKS:
            values[name+"_block_equal_log_loss"]=float(np.mean([
                np.mean([r["losses"][idx,0] for r in valid if r["block"]==b])
                for b in sorted({r["block"] for r in valid})]))
    # Comparative source skill on identical named target species and known genera.
    diff=agg[:,1,:]-agg[:,3,:]
    values.update(
        source_local_minus_global_genus_improvement_island_log_loss=float(diff[:,0].mean()),
        source_local_minus_global_genus_improvement_island_brier=float(diff[:,1].mean()),
        fraction_islands_source_local_genus_improves_logloss=float(np.mean(diff[:,0]>0)),
        n_scored_islands=len(scored_islands),
        n_scored_physical_blocks=n_blocks,
        n_evaluated_original_test_names=n_target_predicted,
        # The in-sample selection of cohort is posthoc; no causal p-values.
        status="POSTHOC_OUT_OF_TARGET_NAME_SOURCE_PREDICTION"
               if n_blocks>=MIN_BLOCKS else "INSUFFICIENT_INDEPENDENT_GEOGRAPHY",
    )
    return dict(**base,**values)


def short_distance_block_disjoint_endpoint_union(ids:list[str],full:dict,
                                                  geog:pd.DataFrame,
                                                  *,low=500.,high=2000.) -> tuple[list[str],int,int]:
    """Outcome-blind exact eligibility graph for the original 456 dyads.

    All original accepted plant names, NOT just flower-known/native records,
    define the species-disjoint island-pair condition. The focal genus/flower
    outcome states do not enter this selection.
    """
    if not ids or low<0 or high<=low:
        raise ValueError("Invalid original target island support/distance window")
    required={"island_id","analysis_regime","spatial_block",
              "island_latitude","island_longitude"}
    if not required<=set(geog) or geog.island_id.duplicated().any():
        raise ValueError("Audited physical island geography required")
    ge=geog.set_index("island_id").reindex(ids)
    if len(set(ids))!=len(ids) or len(ids)<2:
        raise ValueError("Duplicate/missing target original islands")
    rosters=[full[i] for i in ids]
    universe=sorted(set().union(*rosters))
    A=membership(rosters,universe)
    same_species=(A@A.T).toarray()
    lat=pd.to_numeric(ge.island_latitude,errors="coerce").to_numpy()
    lon=pd.to_numeric(ge.island_longitude,errors="coerce").to_numpy()
    dist=geographical_distance(lat,lon,0,len(ids))
    blocks=ge.spatial_block.to_numpy()
    regions=ge.analysis_regime.to_numpy()
    I,J=np.triu_indices(len(ids),1)
    valid=(
        (same_species[I,J]==0) & (regions[I]==regions[J]) &
        pd.notna(regions[I]) & pd.notna(regions[J]) &
        pd.notna(blocks[I]) & pd.notna(blocks[J]) &
        (blocks[I]!=blocks[J]) &
        (dist[I,J]>=low) & (dist[I,J]<high) &
        np.isfinite(lat[I]) & np.isfinite(lat[J]) &
        np.isfinite(lon[I]) & np.isfinite(lon[J]) &
        (np.abs(lat[I])<=90) & (np.abs(lat[J])<=90) &
        (np.abs(lon[I])<=180) & (np.abs(lon[J])<=180))
    ii,jj=I[valid],J[valid]
    endpoint=sorted(set(ids[x] for x in np.r_[ii,jj]))
    visited=set(np.r_[blocks[ii],blocks[jj]])
    return endpoint,int(valid.sum()),len(visited)


def analyze(traits,flora,assignment,mainland,geo):
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original flora fields missing")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original accepted species×island duplicated")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original status uncertain")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo):
        raise ValueError("Independent source geography missing")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    names=sorted(full)
    if not {"flower_primary_color","floral_form"}<=set(traits):
        raise ValueError("Both direct flower trait measurements required")
    common=set(traits["flower_primary_color"])&set(traits["floral_form"])
    labelled={s:(traits["flower_primary_color"][s],traits["floral_form"][s])
              for s in common}
    donors,near=exact_join_and_source(mainland,assignment,labelled)
    g=geo.set_index("island_id")
    n5=[s for s in names if len(full[s]&common)>=5]
    native_n5=[s for s in names if len(native.get(s,set())&common)>=5]
    n10=[
        s for s in n5 if (s in g.index and
           g.at[s,"analysis_regime"]=="northern_midlatitude" and
           len(full[s]&common)>=10 and
           len(full[s]&common)/len(full[s])>=.10)]
    if len(n5)!=2324 or len(native_n5)!=327:
        raise ValueError("Original full 2324/native 327 Direct supports changed")
    if not n10:raise ValueError("Original n10 northern-mid opportunity unavailable")
    short_union,n_short_pairs,n_short_blocks=short_distance_block_disjoint_endpoint_union(
        n10,full,geo)
    samples=[
        ("global_all_n5",n5,"all_recorded"),
        ("global_original_native_n5",native_n5,"original_native"),
        ("northmid_n10_coverage10",n10,"all_recorded"),
        ("northmid_500_2000km_exact456_endpoint_union",short_union,"all_recorded"),
    ]
    result=[]
    for title,ids,scope in samples:
        result.append(predict_one_cohort(title,ids,full,native,labelled,
                                        donors,near,geo,test_scope=scope))
    global_all=set().union(*full.values())
    outside=donors.loc[~donors.work_species.isin(global_all)]
    universal=dict(
        n_original_botanical_islands=len(names),
        n_original_short_distance_full_species_disjoint_pairs=n_short_pairs,
        n_original_short_distance_endpoint_islands=len(short_union),
        n_original_short_distance_spatial_blocks=n_short_blocks,
        n_direct_pair_unique_accepted_species=len(labelled),
        n_original_mainland_native_source_rows=len(mainland),
        n_original_distinct_source_regions=mainland.entity_ID.nunique(),
        n_original_mainland_distinct_names=mainland.work_species.nunique(),
        n_source_direct_pair_species_before_global_name_exclusion=donors.work_species.nunique(),
        n_source_direct_pair_species_after_ALL_4505_island_name_exclusion=outside.work_species.nunique(),
        n_source_direct_pair_genera_after_ALL_4505_island_name_exclusion=len({
            genus(s) for s in outside.work_species.unique()}))
    return dict(
        status="ORIGINAL_STRICT_GIFT_SOURCE_NAME_HOLDOUT_PREFLIGHT",
        **universal,
        cohorts=result,
        claims=dict(
            full_4505_botanical_universe_retained=True,
            direct_species_flower_colour_and_form_joint_not_imputed=True,
            same_name_removed_from_ALL_mainland_regions_before_training=True,
            local_five_source_assignments_not_historical_ancestors=True,
            genus_known_is_a_posthoc_oracle_not_colonization_prediction=True,
            prior_genus_local_predictive_model_is_small_posthoc_sensitivity=True,
            no_double_counting_of_local_donors_overlapping_global_genus=True,
            trait_measurement_species_wide_and_originally_island_derived=True,
            missing_genus_donor_is_unknown_not_true_trait_or_colonization_absence=True,
            direct_source_named_categories_not_true_pollinator_guilds=True,
            no_phylogenetic_independence_or_island_local_visitor_function=True,
            no_private_species_island_GIFT_entity_or_source_coordinates_exported=True,
        ),
    )


def sha256(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","output"):
        p.add_argument("--"+k,required=True,type=Path)
    args=p.parse_args()
    files={
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "source":args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assignment":args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
    }
    for name,path in files.items():
        if not path.is_file() or not path.stat().st_size:
            raise ValueError("Missing exact pinned original source: "+name)
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    }
    signatures={k:sha256(v) for k,v in files.items()}
    if any(signatures[k]!=v for k,v in pins.items()):
        raise ValueError("Original protected source checksum changed")
    traits=strict_traits(pd.read_csv(files["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope",
        "evidence_quality"]))
    flora=pd.read_csv(files["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    assignment=pd.read_csv(files["assignment"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(files["source"],dtype=str,usecols=[
        "entity_ID","work_species"])
    geo=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "island_latitude","island_longitude"])
    if flora.island_id.nunique()!=4505:
        raise ValueError("Source original 4505 island flora count changed")
    if mainland.entity_ID.nunique()!=748:
        raise ValueError("Original mainland GIFT region count changed")
    if len(mainland)!=1666186:
        raise ValueError("Original mainland GIFT native source row count changed")
    result=analyze(traits,flora,assignment,mainland,geo)
    if (result["n_original_short_distance_full_species_disjoint_pairs"]!=456
            or result["n_original_short_distance_endpoint_islands"]!=215):
        raise ValueError("Frozen exact 456-pair/215-island source support drift")
    result["pinned_source_sha256"]={
        "direct_traits":signatures["traits"],
        "original_flora":signatures["flora"],
        "gift_mainland_source_file":signatures["source"],
        "gift_source_assignments_file":signatures["assignment"],
    }
    body=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for forbidden in ('"island_id"','"accepted_species"','"work_species"',
                      '"entity_ID"','"source_signature"','"spatial_block"',
                      '"island_latitude"','"island_longitude"'):
        if forbidden in body:
            raise ValueError("Private botanical identity leaked from source")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n")
    print(json.dumps(dict(status=result["status"],
                          n_botanical_flora=result["n_original_botanical_islands"],
                          n_aggregate_cohorts=len(result["cohorts"]),
                          no_true_colonization_or_pollination_inference=True)))


if __name__=="__main__":main()
