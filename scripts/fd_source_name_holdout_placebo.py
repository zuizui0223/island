"""Source-support-matched GIFT placebo on original 456 floral beta dyad endpoints.

Focal island set frozen OUTCOME-BLIND from the original 4,505 botanical floras:
same broad northern-midlat, original >=10 Direct paired flower traits and >=10%
coverage, 500-2000km distance, other physical 10degree block, NO shared full
accepted plant species across each dyad (456 dyads, 215 active island floras).
ALL accepted names present on any of these original 215 focal island floras
are removed from EVERY native GIFT mainland source entity before learning.

Two alternative five-entity source assignments are matched on FOCAL GENUS
coverage and donor source name count, with different physical 10degree source
blocks, and chosen without observing flower phenotype outcomes or score.
This is not mainland-independent phenotype evidence, ancestor inference or
a pollination functional service measurement.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from fd_island_functional_beta import strict_traits
from fd_source_name_heldout_macro_preflight import (
    genus,donor_probability,score,source_enriched_genus_probability,
    exact_join_and_source,short_distance_block_disjoint_endpoint_union,
)

POLICIES=(("strict2pp",.02,2.0),("sensitivity5pp",.05,2.0))
MIN_TARGETS=30
MIN_BLOCKS=10
SHRINK=8.0


def source_catalog(remaining:pd.DataFrame,near:dict,geog:pd.DataFrame,
                   region:str="northern_midlatitude"):
    """All potential source assignments; no outcome traits are inspected."""
    if not {"island_id","analysis_regime","spatial_block"}<=set(geog):
        raise ValueError("Independent source geography required")
    geo=geog.set_index("island_id")
    by_entity=remaining.groupby("entity_ID").work_species.agg(set).to_dict()
    items={}
    for island,entities in near.items():
        if island not in geo.index or geo.at[island,"analysis_regime"]!=region:
            continue
        block=geo.at[island,"spatial_block"]
        if pd.isna(block):continue
        names=set().union(*(by_entity.get(e,set()) for e in entities))
        if not names:continue
        # Genetic group opportunity, not a new independent mainland trait survey.
        groups={genus(n) for n in names}
        items[island]=dict(
            block=block,source_entities=frozenset(entities),
            names=names,n_names=len(names),genera=groups)
    return items


def coverage(target_genus:Counter,source_genera:set)->float:
    n=sum(target_genus.values())
    if n==0:raise ValueError("Empty original focal species")
    return sum(v for g,v in target_genus.items() if g in source_genera)/n


def source_placebos(focal_id,focal:dict,focal_genus:Counter,
                   catalog:dict,*, tolerance:float, max_ratio:float=2.0,
                   n_controls:int=2):
    """Frozen outcome-blind source matching; NEVER inspect trait-state values."""
    if focal_id not in catalog or tolerance<0 or max_ratio<1:
        raise ValueError("Original matched source availability missing")
    actual=catalog[focal_id]
    y=coverage(focal_genus,actual["genera"])
    choices=[]
    for ident,other in catalog.items():
        if ident==focal_id or other["block"]==focal["block"]:
            continue
        if other["source_entities"]==actual["source_entities"]:
            continue
        jac=len(other["source_entities"]&actual["source_entities"])/len(
            other["source_entities"]|actual["source_entities"])
        if jac>=.8:continue
        c=coverage(focal_genus,other["genera"])
        if abs(y-c)>tolerance+1e-12:continue
        ratio=other["n_names"]/actual["n_names"]
        if not (1/max_ratio<=ratio<=max_ratio):continue
        choices.append((
            abs(y-c),abs(np.log(ratio)),jac,str(ident),ident,c,ratio))
    choices.sort()
    selected=[]
    signatures=set()
    for _,_,jac,_,id,c,ratio in choices:
        sig=catalog[id]["source_entities"]
        if sig in signatures:continue
        signatures.add(sig)
        selected.append(dict(identity=id,genus_coverage=c,donor_ratio=ratio,
                             source_overlap=jac))
        if len(selected)==n_controls:break
    return selected if len(selected)==n_controls else [],len(choices)


def predictor_losses(target_names,global_profiles,prior,genus_counts,
                     source_names,category):
    nclass=len(prior)
    local_counts=defaultdict(Counter)
    for sp in source_names:
        local_counts[genus(sp)][category[sp]]+=1
    total=np.zeros(2,float)
    for sp in target_names:
        g=genus(sp)
        global_count=np.array(
            [genus_counts[g].get(i,0) for i in range(nclass)],float)
        local_count=np.array(
            [local_counts[g].get(i,0) for i in range(nclass)],float)
        p=source_enriched_genus_probability(
            global_count,local_count,prior,shrink=SHRINK)
        x=score(p,category[sp])
        total+=x
    return total/len(target_names)


def summarize(rows:list[dict],*,seed:int=20261115,n_boot:int=1999):
    if not rows:return None
    frame=pd.DataFrame(rows)
    blocks=sorted(frame.block.dropna().unique())
    if len(frame)<MIN_TARGETS or len(blocks)<MIN_BLOCKS:
        return dict(status="insufficient_same-cohort_geo_source_support",
                    n_source_matched_islands=len(frame),n_physical_blocks=len(blocks))
    data=np.stack(frame[["actual_logloss","actual_brier",
                         "placebo_logloss","placebo_brier"]].to_numpy())
    advantage=data[:,2:4]-data[:,0:2]
    group={b:advantage[frame.block.to_numpy()==b].mean(0)
           for b in blocks}
    block_mean=np.stack([group[b] for b in blocks])
    rng=np.random.default_rng(seed)
    chosen=rng.integers(len(blocks),size=(n_boot,len(blocks)))
    samples=block_mean[chosen].mean(1)
    trim=np.sort(advantage[:,0])
    d=dict(
        status="POSTHOC_SUPPORT_MATCHED_SOURCE_ASSIGNMENT_PLACEBO",
        n_source_matched_islands=len(frame),
        n_physical_blocks=len(blocks),
        mean_actual_minus_placebo_improvement_logloss=float(advantage[:,0].mean()),
        mean_actual_minus_placebo_improvement_brier=float(advantage[:,1].mean()),
        block_equal_logloss_advantage=float(block_mean[:,0].mean()),
        block_equal_brier_advantage=float(block_mean[:,1].mean()),
        block_bootstrap_logloss_interval_95=[
            float(np.quantile(samples[:,0],.025)),
            float(np.quantile(samples[:,0],.975))],
        block_bootstrap_brier_interval_95=[
            float(np.quantile(samples[:,1],.025)),
            float(np.quantile(samples[:,1],.975))],
        fraction_islands_actual_source_better_than_matched_placebos=float(
            np.mean(advantage[:,0]>0)),
        mean_absolute_focal_vs_placebo_genus_coverage_difference=float(
            frame.coverage_difference.mean()),
        mean_placebo_to_actual_source_donor_species_ratio=float(
            frame.donor_ratio.mean()),
        mean_original_source_entity_jaccard=float(
            frame.source_overlap.mean()),
        n_original_target_species_carrier_occurrences=int(frame.n_carriers.sum()),
        n_source_matched_islands_positive=int(np.sum(advantage[:,0]>0)),
        n_source_matched_islands_nonpositive=int(np.sum(advantage[:,0]<=0)),
        block_equal_leave_one_out_logloss_min=float(min(
            block_mean[np.arange(len(blocks))!=k,0].mean()
            for k in range(len(blocks)))),
        block_equal_leave_one_out_logloss_max=float(max(
            block_mean[np.arange(len(blocks))!=k,0].mean()
            for k in range(len(blocks)))),
    )
    return d


def analyze(traits,flora,assignment,mainland,geo,*,n_boot:int=1999):
    if not {"accepted_species","island_id","origin_status"}<=set(flora):
        raise ValueError("Missing original native island botanical species")
    if flora.duplicated(["island_id","accepted_species"]).any() or (
        not flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original botanical species/status invalid")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("Original global 4505 botanical universe changed")
    direct=({sp:(traits["flower_primary_color"][sp],traits["floral_form"][sp])
             for sp in set(traits["flower_primary_color"])&set(traits["floral_form"])})
    geo=geo.set_index("island_id")
    focal=[
        id for id,species in full.items()
        if id in geo.index and geo.at[id,"analysis_regime"]=="northern_midlatitude"
        and len(species&set(direct))>=10 and
        len(species&set(direct))/len(species)>=.10]
    cohort,n_pairs,nblocks=short_distance_block_disjoint_endpoint_union(
        sorted(focal),full,geo.reset_index())
    if (n_pairs,nblocks,len(cohort))!=(456,29,215):
        raise ValueError("Original exact 456 pair/29 block/215 island cohort drift")
    original_direct,near=exact_join_and_source(mainland,assignment,direct)
    heldout_names=set().union(*(full[id] for id in cohort))
    train=original_direct.loc[
        ~original_direct.work_species.isin(heldout_names)].copy()
    if set(train.work_species)&heldout_names:
        raise AssertionError("Original focal accepted plant name leaked")
    candidates=set(train.work_species)
    if not candidates:raise ValueError("No truly name-disjoint original Direct donor")
    labels=sorted(set(direct.values()))
    category={sp:labels.index(profile) for sp,profile in direct.items()}
    k=len(labels)
    global_count=Counter(category[sp] for sp in candidates)
    prior=np.array([global_count.get(i,0)+.5 for i in range(k)],float)
    prior/=prior.sum()
    genus_counts=defaultdict(Counter)
    for sp in candidates:
        genus_counts[genus(sp)][category[sp]]+=1
    catalog=source_catalog(train,near,geo.reset_index())
    reports=[]
    for title,tol,ratio in POLICIES:
        outcome_blind=[]
        source_matched=[]
        counts=dict(n_focal=len(cohort),n_actual_present_in_catalog=0,
                    n_focal_islands_with_two_source_controls=0,
                    n_focal_islands_without_two_controls=0,
                    n_candidate_alt_source_assignments_sum=0)
        for id in cohort:
            target=sorted(full[id]&set(direct))
            if id not in catalog:
                counts["n_focal_islands_without_two_controls"]+=1
                continue
            counts["n_actual_present_in_catalog"]+=1
            g=Counter(genus(sp) for sp in target)
            controls,n_alt=source_placebos(
                id,dict(block=geo.at[id,"spatial_block"]),g,
                catalog,tolerance=tol,max_ratio=ratio)
            counts["n_candidate_alt_source_assignments_sum"]+=n_alt
            if not controls:
                counts["n_focal_islands_without_two_controls"]+=1
                continue
            counts["n_focal_islands_with_two_source_controls"]+=1
            actual=catalog[id]
            outcome_blind.append((
                id,target,actual,controls,g))
        for id,target,actual,controls,g in outcome_blind:
            obs=predictor_losses(target,direct,prior,genus_counts,
                                  actual["names"],category)
            comp=np.array([
                predictor_losses(target,direct,prior,genus_counts,
                                 catalog[c["identity"]]["names"],category)
                for c in controls])
            same=comp.mean(0)
            cov=coverage(g,actual["genera"])
            source_matched.append(dict(
                block=geo.at[id,"spatial_block"],
                n_carriers=len(target),
                actual_logloss=obs[0],actual_brier=obs[1],
                placebo_logloss=same[0],placebo_brier=same[1],
                coverage_difference=np.mean([
                    abs(cov-c["genus_coverage"]) for c in controls]),
                donor_ratio=np.mean([c["donor_ratio"] for c in controls]),
                source_overlap=np.mean([c["source_overlap"] for c in controls])))
        stats=summarize(source_matched,n_boot=n_boot)
        reports.append(dict(policy=title,max_coverage_abs_difference=tol,
                            source_donor_pool_size_ratio_max=ratio,
                            **counts,**(stats or dict(status="NO_SOURCE_PLACEBO_SUPPORT"))))
    return dict(
        status="ORIGINAL_456PAIR_SOURCE_MATCHED_PLACEBO",
        n_original_islands=len(full),
        n_original_pair_endpoints=len(cohort),n_pair_opportunities=n_pairs,
        n_original_spatial_blocks=nblocks,
        n_out_of_target_cohort_mainland_direct_donors=len(candidates),
        n_external_mainland_candidate_source_islands=len(catalog),
        n_observed_direct_pair_category_states=len(labels),
        n_observed_direct_pair_accepted_species=len(direct),
        comparisons=reports,
        claims=dict(
            outcome_blind_placebo_choice=True,
            global_exclusion_of_all_focal_full_flora_names=True,
            out_of_name_not_out_of_original_island_trait_provider=True,
            donor_species_deduplicated_across_GIFT_source_entities=True,
            placebo_source_evidence_matched_not_ancestral_assignment=True,
            known_target_genus_is_posthoc_oracle=True,
            physical_blocks_not_independent_phylogenetic_histories=True,
            no_island_local_trait_evolution_pollen_fitness_or_guild_shift=True,
            no_private_plant_island_GIFT_entity_or_geography_ids_exported=True))


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda:fp.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for a in ("original-root","geology","output"):
        p.add_argument("--"+a,type=Path,required=True)
    args=p.parse_args()
    paths={
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assignment":args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        "mainland":args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz"}
    expected={"traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
              "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
              "assignment":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
              "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef"}
    for key,path in paths.items():
        if sha256(path)!=expected[key]:
            raise ValueError("Source archive SHA mismatch: "+key)
    s=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    floras=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    assigned=pd.read_csv(paths["assignment"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    gift=pd.read_csv(paths["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    geo=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "island_latitude","island_longitude"])
    if len(gift)!=1666186 or gift.entity_ID.nunique()!=748:
        raise ValueError("GIFT original 748 mainland regions changed")
    q=analyze(s,floras,assigned,gift,geo)
    if q["n_original_pair_endpoints"]!=215 or q["n_pair_opportunities"]!=456:
        raise ValueError("Original eligible 456 island pairs changed")
    if q["n_out_of_target_cohort_mainland_direct_donors"]!=2899:
        raise ValueError("Frozen out-of-215-colour-form-trait GIFT donors changed")
    body=json.dumps(q,indent=2,sort_keys=True,allow_nan=False)
    forbidden=('"island_id"','"accepted_species"','"work_species"',
               '"source_signature"','"entity_ID"','"spatial_block"',
               '"island_latitude"','"island_longitude"','"identity"')
    if any(x in body for x in forbidden):
        raise ValueError("Original plant or source identity escaped")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n")
    print(json.dumps(dict(status=q["status"],n_matched_policies=len(q["comparisons"]),
                          n_original_pair_endpoints=215,source_independence_not_claimed=True)))


if __name__=="__main__":main()
