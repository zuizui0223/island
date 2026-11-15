"""Original 4505 island floral COLOUR composition vs independently published NHM text-source native GIFT donors.

Focal island botanical colour traits are original high/medium single-state
species-direct observations. ALL external training donors are different exact
species names absent from the ENTIRE 4505-island flora. External NHM (CC BY 4)
primary flower colour and separately, organ-expanded descriptive colours are
NOT direct empirical mainland-population phenotype or independent pollen data.

Proper categorical score comparison:
 global independently described GIFT source prior vs five-nearest native
 mainland candidate source text-colour composition vs two outcome-blind
 source-size/genus-richness matched pseudo-source assignments from other
 10degree physical spatial blocks within the same broad regime.

Main macro population remains 4505 botanical islands; do not restrict to
triple-complete 928 or selected posthoc 215 islands. Nothing here proves a
historical colonization source, selection, pollination or trait evolution.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
import pandas as pd
from fd_island_functional_beta import strict_traits
from fd_nhm_external_flora_preflight import (
    nhm_git_blob_sha,nhm_external_categories,NHM_BLOB,NHM_COMMIT,
    EXTERNAL_COLUMNS,accepted_exact_binomial)
from fd_source_name_heldout_macro_preflight import (
    short_distance_block_disjoint_endpoint_union,
)

ALPHA=20.0
MIN_DIRECT=5
MIN_ISLANDS=30
MIN_BLOCKS=10
BOOT=1999
MATCH_MIN_RATIO=.5
MATCH_MAX_RATIO=2.0
SOURCE_JACCARD_MAX=.8
NHM_ALLOWED=frozenset((
    "white","red_pink","yellow_orange","blue_purple",
    "green_brown_inconspicuous",
))


def source_probability(source_counts,global_prior,alpha=ALPHA):
    x=np.asarray(source_counts,dtype=float)
    p=np.asarray(global_prior,dtype=float)
    if (alpha<=0 or x.shape!=p.shape or np.any(x<0) or np.any(p<=0)
            or not np.isclose(p.sum(),1) or not np.isfinite(x+p).all()):
        raise ValueError("Invalid category prior or source opportunity")
    q=(x+alpha*p)/(x.sum()+alpha)
    if not np.isclose(q.sum(),1):raise AssertionError("Source model malformed")
    return q


def proper_flora_losses(known_categories:list[int],p:np.ndarray):
    """Mean over distinct accepted original direct flower-colour species."""
    y=np.asarray(known_categories,int)
    if not len(y) or y.min()<0 or y.max()>=len(p) or np.any(p<=0):
        raise ValueError("Invalid original trait-known target flower colours")
    n=len(y)
    logloss=-float(np.log(p[y]).mean())
    # Full-multiclass strictly proper Brier, 1-2*p_y+sum(p_c²)
    brier=1.0-2.0*float(p[y].mean())+float(np.sum(p*p))
    return np.array([logloss,brier])


def build_external_source_catalog(land:pd.DataFrame,assign:pd.DataFrame,
                                  geo:pd.DataFrame,eligible_states:dict,
                                  cats:dict):
    if not {"entity_ID","work_species"}<=set(land):
        raise ValueError("Original accepted GIFT source identifiers missing")
    if land.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("Native GIFT mainland duplicated entity×species")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assign):
        raise ValueError("GIFT original geographic candidate assignments missing")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo):
        raise ValueError("Geographical independence source missing")
    used=land.loc[land.work_species.isin(eligible_states)]
    by=used.groupby("entity_ID").work_species.agg(set).to_dict()
    rank=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    rank["source_rank"]=pd.to_numeric(rank.source_rank,errors="raise")
    rank=rank.loc[rank.source_rank.between(1,5)]
    if rank.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Geographical source ranks repeated")
    near=rank.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    places=geo.set_index("island_id")
    catalog={}
    for island,entities in near.items():
        if island not in places.index:continue
        region=places.at[island,"analysis_regime"]
        block=places.at[island,"spatial_block"]
        if pd.isna(region) or pd.isna(block):continue
        names=set().union(*(by.get(x,set()) for x in entities))
        if not names:continue
        values=np.bincount([cats[eligible_states[sp]] for sp in names],
                           minlength=len(cats)).astype(float)
        catalog[island]=dict(
            region=str(region),block=str(block),source_entities=frozenset(entities),
            first_source=entities[0],n_source_species=len(names),
            n_source_genera=len({s.split()[0] for s in names}),
            values=values)
    return catalog


def placebo_assignment(focal_id,catalog,by_region,*,size_factor=MATCH_MAX_RATIO,
                       n_controls:int=2):
    """Choose source assignments only by source evidence, NO flower outcomes."""
    actual=catalog.get(focal_id)
    if not actual:return []
    ref_n=actual["n_source_species"]
    ref_g=actual["n_source_genera"]
    if ref_n<1 or ref_g<1: return []
    alternatives=[]
    for ident in by_region.get(actual["region"],[]):
        if ident==focal_id:continue
        other=catalog[ident]
        if other["block"]==actual["block"]:continue
        if other["source_entities"]==actual["source_entities"]:continue
        union=actual["source_entities"]|other["source_entities"]
        share=len(actual["source_entities"]&other["source_entities"])/len(union)
        if share>=SOURCE_JACCARD_MAX:continue
        r=other["n_source_species"]/ref_n
        g=other["n_source_genera"]/ref_g
        if not 1/size_factor<=r<=size_factor:continue
        if not 1/size_factor<=g<=size_factor:continue
        alternatives.append((abs(np.log(r))+abs(np.log(g)),share,str(ident),ident,r,g))
    # Only one source assignment for each distinct SOURCE ENTITY SET.
    chosen=[]
    seen=set()
    for _,jac,_,id,r,g in sorted(alternatives):
        identity=catalog[id]["source_entities"]
        if identity in seen:continue
        seen.add(identity)
        chosen.append(dict(island=id,source_donor_ratio=r,
                           source_genus_ratio=g,source_entity_jaccard=jac))
        if len(chosen)==n_controls:break
    return chosen if len(chosen)==n_controls else []


def _bootstrap_equal_group(rows:list,groupkey:str,seed:int=20261115,draws:int=BOOT):
    if not rows:return None
    groups=defaultdict(list)
    for r in rows:groups[r[groupkey]].append(r["contrast"])
    groups={k:np.asarray(v).mean(0) for k,v in groups.items()}
    if len(groups)<2:return None
    values=np.stack(list(groups.values()))
    mean=values.mean(0)
    rng=np.random.default_rng(seed)
    ix=rng.integers(len(values),size=(draws,len(values)))
    sim=values[ix].mean(1)
    return dict(
        n_groups=len(values),
        mean_gain_logloss=float(mean[0]),
        mean_gain_brier=float(mean[1]),
        blockwise_bootstrap_logloss_95=[
            float(np.quantile(sim[:,0],.025)),
            float(np.quantile(sim[:,0],.975))],
        blockwise_bootstrap_brier_95=[
            float(np.quantile(sim[:,1],.025)),
            float(np.quantile(sim[:,1],.975))],
        leave_one_group_out_logloss_min=float(min(
            values[np.arange(len(values))!=i,0].mean() for i in range(len(values)))),
        leave_one_group_out_logloss_max=float(max(
            values[np.arange(len(values))!=i,0].mean() for i in range(len(values)))),
    )


def analyze_categorical_colour_source(original:dict,floras:pd.DataFrame,
                                      places:pd.DataFrame,land:pd.DataFrame,
                                      assignments:pd.DataFrame,
                                      external:pd.DataFrame, *,
                                      boots:int=BOOT) -> dict:
    if not {"island_id","accepted_species","origin_status"}<=set(floras):
        raise ValueError("Original biological flora missing")
    if floras.duplicated(["island_id","accepted_species"]).any() or (
        not floras.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original botanical accepted source conflict")
    all_plant=floras.groupby("island_id").accepted_species.agg(set).to_dict()
    native=floras.loc[floras.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(all_plant)!=4505:raise ValueError("Global island botanical source changed")
    universe=set().union(*all_plant.values())
    nhm,summary=nhm_external_categories(external)
    direct=original["flower_primary_color"]
    direct_form=original["floral_form"]
    colour_supported=set(direct)&NHM_ALLOWED
    cats={c:i for i,c in enumerate(sorted(NHM_ALLOWED))}
    original_colour_states={s:cats[col] for s,col in direct.items() if col in cats}
    if len(original_colour_states)<1000:raise ValueError("Original colour label ontology changed")
    geo=places.set_index("island_id")
    color_n5=sorted(i for i in all_plant
                    if len(all_plant[i]&set(original_colour_states))>=MIN_DIRECT)
    native_n5=sorted(i for i in all_plant
                    if len(native.get(i,set())&set(original_colour_states))>=MIN_DIRECT)
    paired=set(direct)&set(direct_form)
    candidate_n10=sorted(i for i in all_plant
                         if i in geo.index and geo.at[i,"analysis_regime"]=="northern_midlatitude"
                         and len(all_plant[i]&paired)>=10 and
                         len(all_plant[i]&paired)/len(all_plant[i])>=.10)
    endpoints,n_pairs,n_blocks=short_distance_block_disjoint_endpoint_union(
        candidate_n10,all_plant,places)
    if (len(endpoints),n_pairs,n_blocks)!=(215,456,29):
        raise ValueError("Original 215 focal-flora source-dyad opportunity changed")
    run=[]
    gift_names=set(land.work_species)
    for tier in ("primary","expanded"):
        eligible={
            s:label for s,label in nhm[tier].items()
            if s in gift_names and s not in universe and label in cats}
        if not eligible:continue
        label_freq=Counter(eligible.values())
        global_prior=np.array([label_freq.get(k,0)+.5 for k in sorted(cats)],
                              dtype=float)
        global_prior/=global_prior.sum()
        catalog=build_external_source_catalog(
            land,assignments,places,eligible,cats)
        by_region=defaultdict(list)
        for name,r in catalog.items():
            by_region[r["region"]].append(name)
        for x in by_region.values():x.sort()
        scopes=[
            ("all_recorded_n5",color_n5,all_plant),
            ("original_native_n5",native_n5,native),
            ("northmid_exact456_endpoint215",endpoints,all_plant),
        ]
        reports=[]
        for scope,ids,rosters in scopes:
            rows=[]
            raw=[]
            n_total_test_carriers=0
            n_uncovered=0
            n_with_local=0
            for island in ids:
                names=sorted(rosters.get(island,set())&set(original_colour_states))
                if len(names)<MIN_DIRECT:continue
                n_total_test_carriers+=len(names)
                focal=catalog.get(island)
                if not focal:
                    n_uncovered+=1
                    continue
                n_with_local+=1
                target=[original_colour_states[s] for s in names]
                actual=source_probability(focal["values"],global_prior)
                actual_loss=proper_flora_losses(target,actual)
                prior_loss=proper_flora_losses(target,global_prior)
                raw.append(dict(region=focal["region"],block=focal["block"],
                                first_source=focal["first_source"],
                                global_delta=prior_loss-actual_loss))
                controls=placebo_assignment(island,catalog,by_region)
                if len(controls)!=2:continue
                p_loss=np.mean([
                    proper_flora_losses(target,source_probability(
                        catalog[c["island"]]["values"],global_prior))
                    for c in controls],axis=0)
                rows.append(dict(
                    region=focal["region"],block=focal["block"],
                    first_source=focal["first_source"],
                    contrast=p_loss-actual_loss,
                    n_direct=len(names),n_local_source=focal["n_source_species"],
                    species_ratio=np.mean([c["source_donor_ratio"] for c in controls]),
                    genera_ratio=np.mean([c["source_genus_ratio"] for c in controls]),
                    entity_overlap=np.mean([c["source_entity_jaccard"] for c in controls]),
                ))
            regions=sorted(set(r["region"] for r in raw))
            summaries=[]
            for region in ["GLOBAL_ALL",*regions]:
                r0=raw if region=="GLOBAL_ALL" else [
                    r for r in raw if r["region"]==region]
                r1=rows if region=="GLOBAL_ALL" else [
                    r for r in rows if r["region"]==region]
                clusters={r["block"] for r in r1}
                if len(r1)<MIN_ISLANDS or len(clusters)<MIN_BLOCKS:
                    summaries.append(dict(
                        region=region,status="INSUFFICIENT_MATCHED_SOURCE_SUPPORT",
                        n_islands_with_any_local_source=len(r0),
                        n_islands_with_two_pseudo_source_controls=len(r1),
                        n_distinct_physical_blocks=len(clusters)))
                    continue
                pair=np.stack([r["contrast"] for r in r1])
                glob=np.stack([r["global_delta"] for r in r0])
                summary=dict(region=region,
                    status="EXPLORATORY_GLOBAL_EXTERNAL_SOURCE_TEXT_COLOR_PREDICTION",
                    n_islands_with_any_local_source=len(r0),
                    n_islands_with_two_pseudo_source_controls=len(r1),
                    n_distinct_physical_blocks=len(clusters),
                    n_first_GIFT_source_groups=len({r["first_source"] for r in r1}),
                    island_equal_source_vs_placebo_gain_logloss=float(pair[:,0].mean()),
                    island_equal_source_vs_placebo_gain_brier=float(pair[:,1].mean()),
                    island_equal_source_vs_global_prior_gain_logloss_on_ALL_source_available=float(glob[:,0].mean()),
                    island_equal_source_vs_global_prior_gain_brier_on_ALL_source_available=float(glob[:,1].mean()),
                    fraction_scored_islands_real_source_better_logloss=float(np.mean(pair[:,0]>0)),
                    mean_matched_placebo_to_focal_external_source_species_ratio=float(
                        np.mean([r["species_ratio"] for r in r1])),
                    mean_matched_placebo_to_focal_external_source_genera_ratio=float(
                        np.mean([r["genera_ratio"] for r in r1])),
                    mean_source_entity_Jaccard_between_real_and_matched_placebo=float(
                        np.mean([r["entity_overlap"] for r in r1])),
                    **{key:value for key,value in (
                        ("physical_block_equal",_bootstrap_equal_group(
                            r1,"block",seed=20261115,draws=boots)),
                        ("first_GIFT_source_equal",_bootstrap_equal_group(
                            r1,"first_source",seed=20261116,draws=boots))
                       ) if value is not None})
                summaries.append(summary)
            reports.append(dict(
                scope=scope,
                n_original_source_islands_eligible=len(ids),
                n_original_direct_colour_carrier_rows_in_scope=n_total_test_carriers,
                n_original_islands_without_any_source_donors=n_uncovered,
                n_original_islands_with_external_source_donors=n_with_local,
                summaries=summaries,
            ))
        run.append(dict(
            external_colour_evidence_tier=tier,
            n_distinct_external_primary_colour_source_only_species=len(eligible),
            n_distinct_external_source_only_genera=len({x.split()[0] for x in eligible}),
            n_independent_geographic_island_source_assignments_with_external_colour=len(catalog),
            n_original_4505_botanical_floras_with_near_external_colour_donors=sum(
                x in catalog for x in all_plant),
            n_original_external_colour_category_states=len(cats),
            source_ratio_minimum=MATCH_MIN_RATIO,
            source_ratio_maximum=MATCH_MAX_RATIO,
            source_entity_jaccard_max_exclusive=SOURCE_JACCARD_MAX,
            NHM_independent_original_corpora_provenance_not_available_in_combined=True,
            scope_evidence=reports))
    return dict(
        status="ORIGINAL_NHM_EXTERNAL_GIFT_SOURCE_ONLY_COLOUR_MACRO_PREDICTION",
        n_original_botanical_islands=len(all_plant),
        n_original_exact_source_flora_species=len(universe),
        n_original_direct_colour_species=len(direct),
        n_original_direct_colour_species_in_external_five_class_ontology=len(original_colour_states),
        n_original_colour_n5_island_floras=len(color_n5),
        n_original_native_colour_n5_island_floras=len(native_n5),
        n_original_posthoc_456pair_endpoint_islands=len(endpoints),
        results=run,
        claims=dict(
            training_NHM_text_source_species_absent_from_ALL_original_4505_island_floras=True,
            zero_name_overlap_mainland_external_train_vs_original_island_colour_targets=True,
            main_global_universe_not_replaced_by_complete_case_928=True,
            flower_colour_endpoint_does_not_require_form_or_symmetry_completeness=True,
            external_text_traits_not_mainland_population_measured_phenotypes=True,
            local_source_vs_2_size_matched_placebos_choice_outcome_blind=True,
            botanical_source_plants_not_independent_ancestral_histories=True,
            physical_and_first_source_cluster_conditionals_not_independent_replications=True,
            no_plant_adaptation_pollen_effectiveness_or_guild_shift_claim=True,
            no_original_taxon_island_source_entities_or_coordinates_exported=True,
        ))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for x in ("original-root","nhm-combined","geology","output"):
        p.add_argument("--"+x,type=Path,required=True)
    a=p.parse_args()
    original={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "gift":a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz",
    }
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "gift":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
    }
    for k in original:
        if sha(original[k])!=pins[k]:
            raise ValueError("Original protected GIFT colour baseline SHA drift: "+k)
    if nhm_git_blob_sha(a.nhm_combined)!=NHM_BLOB:
        raise ValueError("Original NHM combined botanical source drift")
    cols=set(pd.read_csv(a.nhm_combined,nrows=0).columns)
    required=set(EXTERNAL_COLUMNS)-{"source"}
    if not required<=cols:raise ValueError("External NHM flower schema mismatch")
    external=pd.read_csv(a.nhm_combined,dtype=str,keep_default_na=False,
                         low_memory=False,usecols=sorted(required))
    if "source" not in external:external["source"]="combined_source_unresolved"
    traits=strict_traits(pd.read_csv(original["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    gift=pd.read_csv(original["gift"],dtype=str,usecols=[
        "entity_ID","work_species"])
    source=pd.read_csv(original["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    places=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "island_latitude","island_longitude"])
    outcome=analyze_categorical_colour_source(traits,flora,places,gift,source,external)
    if (outcome["n_original_botanical_islands"]!=4505 or
        outcome["n_original_posthoc_456pair_endpoint_islands"]!=215):
        raise ValueError("Original source 4505/215 opportunity mismatch")
    outcome["NHM_public_source_metadata"]=dict(
        CC_BY_4_0="https://github.com/NaturalHistoryMuseum/plant-traits",
        immutable_commit=NHM_COMMIT,combined_blob_sha=NHM_BLOB,
        full_source_sha256=sha(a.nhm_combined))
    body=json.dumps(outcome,indent=2,sort_keys=True,allow_nan=False)
    for forbidden in ('"island_id"','"accepted_species"','"work_species"',
                      '"entity_ID"','"source_signature"','"spatial_block"',
                      '"island_latitude"','"island_longitude"'):
        if forbidden in body:
            raise ValueError("Private GIFT source IDs escaped")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(status=outcome["status"],
                          n_main_island_flora=4505,
                          n_original_flower_color_n5_eligible=outcome["n_original_colour_n5_island_floras"],
                          n_source_evidence_tiers=len(outcome["results"]),
                          no_causal_pollination_source_inference=True)))


if __name__=="__main__":main()
