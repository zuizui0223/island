"""Direct original-plant source/strict native-nonendemic/endemic/introduced morphology.

Whole population 4505 original GIFT island botanical floras (do NOT use
n30 three-trait 928 complete-case). Independent original status layers:
origin_status, endemic_status, floristic_status. Only direct species-high/
medium accepted botanical trait states. Mainland pool = candidate GIFT
five-nearest *native* region source species not in focal FULL island flora;
candidate omission is NOT failure to arrive or establish.

Three direct observational comparisons:
(1) strict native_nonendemic island species share minus source candidate-only
    species share without genus conditioning (real genus replacement retained);
(2) original island endemic minus same-island native_nonendemic source
    Direct trait category share, with separately same-island matching GENUS;
(3) original introduced island species minus source candidate-only; this
    does NOT control propagule pressure/cultivation/observation.
These are NOT three causal experimental stages, true adapted traits,
population-specific phenotype or pollen-effective guild evidence.
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

REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
SCOPES=("all_recorded",)
ENDPOINTS=(
    ("zygomorphic","floral_symmetry","zygomorphic"),
    ("tube_or_bell","floral_form","tube_or_bell"),
    ("open_radial","floral_form","open_radial"),
    ("deep_tube","tube_depth_class","deep"),
    ("SC_vs_SI","self_incompatibility","SC"),
    ("autonomous_vs_absent","autonomous_selfing_capacity","autonomous"),
    ("bilateral_tube_or_bell_joint","floral_form__floral_symmetry","joint_special"),
)
MIN_DIRECT_ISLAND=5
MIN_DIRECT_SOURCE=10
MIN_GROUP_ISLANDS=30
MIN_INDEPENDENT_BLOCKS=10
MIN_VARIABLE_GENERA=10
BOOTS=999
BINARY={
 "self_incompatibility":{"SC","SI"},
 "autonomous_selfing_capacity":{"autonomous","absent"}
}


def genus(taxon:str):
    parts=str(taxon).split()
    return parts[0] if len(parts)>=2 else "single_source_name:"+str(taxon)


def status_groups(frame:pd.DataFrame):
    required={"island_id","accepted_species","origin_status",
              "endemic_status","floristic_status"}
    if not required<=set(frame) or frame.duplicated(
            ["island_id","accepted_species"]).any():
        raise ValueError("Original 4505 GIFT accepted flora/states nonunique")
    if not frame.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original native/introduced statuses inconsistent")
    if not frame.endemic_status.isin(("endemic","nonendemic","unresolved")).all():
        raise ValueError("Original separate endemicity status invalid")
    if (frame.floristic_status.eq("endemic") & (
        ~frame.origin_status.eq("native") |
        ~frame.endemic_status.eq("endemic"))).any():
        raise ValueError("Original claimed endemic contradicts indigenous source status")
    if (frame.floristic_status.eq("native_nonendemic") & (
        ~frame.origin_status.eq("native") |
        ~frame.endemic_status.eq("nonendemic"))).any():
        raise ValueError("Original nonendemic source status contradicts indigenous origin")
    d={}
    for island,g in frame.groupby("island_id",sort=True):
        native=set(g.loc[g.origin_status.eq("native"),"accepted_species"])
        e=set(g.loc[
            g.origin_status.eq("native") &
            g.endemic_status.eq("endemic") &
            g.floristic_status.eq("endemic"),"accepted_species"])
        ne=set(g.loc[
            g.origin_status.eq("native") &
            g.endemic_status.eq("nonendemic") &
            g.floristic_status.eq("native_nonendemic"),"accepted_species"])
        intro=set(g.loc[g.origin_status.eq("introduced"),"accepted_species"])
        all_names=set(g.accepted_species)
        if e&ne or native&intro or not (e|ne)<=native:
            raise AssertionError("Endemic/native/introduced source groups overlap")
        d[island]=dict(full=all_names,native=native,endemic=e,
            native_nonendemic=ne,introduced=intro,
            native_endemism_unresolved=native-e-ne,
            unresolved_origin=all_names-native-intro)
    return d


def source_sectors(states:dict):
    necessary={"floral_symmetry","floral_form","tube_depth_class",
               "self_incompatibility","autonomous_selfing_capacity"}
    if not necessary<=set(states):
        raise ValueError("Original individually Direct flora sector ontology missing")
    x={}
    for name,trait,goal in ENDPOINTS:
        if trait=="floral_form__floral_symmetry":
            a,b=states["floral_form"],states["floral_symmetry"]
            x[name]={sp:int(b[sp]=="zygomorphic" and a[sp] in
                {"tubular","bell_campanulate"}) for sp in set(a)&set(b)}
        else:
            d=states[trait]
            if trait in BINARY:d={sp:v for sp,v in d.items()
                                  if v in BINARY[trait]}
            if goal=="tube_or_bell":
                x[name]={sp:int(v in {"tubular","bell_campanulate"})
                         for sp,v in d.items()}
            else:
                x[name]={sp:int(v==goal) for sp,v in d.items()}
        if not x[name]:raise ValueError("Empty source Direct trait state")
    return x


def original_shared_genus(endemic_names:set,nonendemic_names:set,
                          binary_states:dict):
    """Observed same-island same-genus mean differences, not sister comparison."""
    eg=defaultdict(list);ng=defaultdict(list)
    for s in endemic_names&set(binary_states):eg[genus(s)].append(binary_states[s])
    for s in nonendemic_names&set(binary_states):ng[genus(s)].append(binary_states[s])
    common=sorted(set(eg)&set(ng))
    diffs=[];variable=0
    for g in common:
        vals=eg[g]+ng[g]
        variable+=int(len(set(vals))>1)
        diffs.append(float(np.mean(eg[g])-np.mean(ng[g])))
    return dict(n_matched_same_island_genus_units=len(common),
                n_matched_same_island_genus_with_actual_category_variation=variable,
                sum_genus_equal_endemic_minus_native_nonendemic_share=float(sum(diffs)),
                mean_genus_equal_endemic_minus_native_nonendemic_share=(
                    float(np.mean(diffs)) if diffs else None))


def grouped_interval(records:list,cluster:str,*,boots:int=BOOTS,seed:int=20261115):
    if not records:return dict(status="NOT_IDENTIFIABLE_NO_DIRECT_COMPARABLE_FLORAS",
                                n_comparable_floras=0,n_source_clusters=0)
    data=pd.DataFrame(records)
    eligible=data.dropna(subset=[cluster,"difference"])
    n=len(eligible)
    groups=eligible[cluster].nunique()
    summary=dict(
        n_comparable_floras=n,n_source_clusters=groups,
        original_island_equal_difference=float(eligible.difference.mean())
        if n else None,
        n_positive_difference_islands=int((eligible.difference>0).sum()) if n else 0)
    if n<MIN_GROUP_ISLANDS or groups<MIN_INDEPENDENT_BLOCKS:
        return dict(**summary,status="NOT_IDENTIFIABLE_ORIGINAL_SOURCE_CLUSTER_SUPPORT")
    vals=eligible.groupby(cluster,sort=True).difference.mean().to_numpy(float)
    rng=np.random.default_rng(seed)
    draw=vals[rng.integers(groups,size=(boots,groups))].mean(axis=1)
    loo=[float(vals[np.arange(groups)!=i].mean()) for i in range(groups)]
    return dict(**summary,
        status="DESCRIPTIVE_ORIGINAL_BLOCK_EQUAL_FLORISTIC_SECTOR_CONTRAST",
        original_cluster_equal_difference=float(vals.mean()),
        original_cluster_conditional_resample_95=[
            float(np.quantile(draw,.025)),float(np.quantile(draw,.975))],
        original_leave_one_cluster_out_min=min(loo),
        original_leave_one_cluster_out_max=max(loo),
        cluster_draws_not_independent_plant_histories=True)


def analyze(states:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            assignments:pd.DataFrame,mainland:pd.DataFrame,
            *,boots:int=BOOTS):
    clusters={"island_id","analysis_regime","spatial_block",
              "geological_origin"}
    if not clusters<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Original independently mapped GSHHG geography invalid")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assignments):
        raise ValueError("Original five-nearest mainland GIFT entity source invalid")
    if not {"entity_ID","work_species"}<=set(mainland):
        raise ValueError("Native GIFT mainland species opportunity missing")
    d=status_groups(flora)
    if len(d)!=4505:raise ValueError("Main 4505 original island floras changed")
    vocab=source_sectors(states)
    # Exactly one botanical island has unresolved physical GSHHG geography.
    # Preserve it in the 4,505 population and account for it explicitly;
    # do not assign a fictitious mainland location or a regional test group.
    region_geo=geo.set_index("island_id").reindex(sorted(d))
    if int(region_geo.analysis_regime.notna().sum())!=4504:
        raise ValueError("Original 4504 geographically matched botanical floras changed")
    source=assignments.loc[assignments.source_mode.eq("geo_k5")].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise")
    source=source.loc[source.source_rank.between(1,5)]
    if source.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Repeated source potential GIFT region rank")
    near=source.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("Original mainland source duplicated species")
    byentity=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    source_count=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    all_direct=set().union(*(set(d) for d in vocab.values()))
    # Build each island source opportunity once, not by trait or status;
    # every matched source species belongs to GIFT native mainland listings.
    byisland={}
    byisland_mainland_shared={}
    coverage=Counter()
    for id,entry in d.items():
        five=near.get(id,())
        total=(set().union(*(byentity.get(s,set()) for s in five))
               if five else set())
        only=(total-entry["full"])
        coverage["n_islands_with_five_source_regions"]+=int(len(five)==5)
        coverage["n_islands_with_qualifying_source_only_names"]+=int(bool(only))
        coverage["n_candidate_source_island_plant_occurrences"]+=len(only)
        byisland[id]=only&all_direct
        byisland_mainland_shared[id]=total&entry["full"]
    # Direct carrier status counts remain per original source individual island.
    byregion=Counter()
    status_counts=Counter()
    for id,entry in d.items():
        reg=region_geo.at[id,"analysis_regime"]
        region=(str(reg) if pd.notna(reg) and str(reg).strip()
                else "geography_unresolved")
        byregion[region]+=1
        for status in ("native_nonendemic","endemic","introduced",
                       "native_endemism_unresolved","unresolved_origin"):
            status_counts[status+"_source_island_x_plant_rows"]+=len(entry[status])
    analyses=[]
    census=[]
    rows_for={}
    for region in sorted(byregion):
        ids=[i for i in sorted(d) if (
            str(region_geo.at[i,"analysis_regime"])
            if pd.notna(region_geo.at[i,"analysis_regime"])
            else "geography_unresolved")==region]
        if not ids and region!="geography_unresolved":continue
        for endpoint,trait in vocab.items():
            state_names=set(trait)
            records={key:[] for key in (
                "native_nonendemic_minus_mainland_not_recorded",
                "native_nonendemic_GIFT_source_shared_minus_mainland_not_recorded",
                "introduced_minus_mainland_not_recorded",
                "endemic_minus_native_nonendemic_raw",
                "endemic_minus_native_nonendemic_same_island_genus")}
            support=Counter()
            status_species={r:set() for r in (
                "native_nonendemic","native_nonendemic_GIFT_source_shared",
                "endemic","introduced","mainland_not_recorded",
                "native_endemism_unresolved")}
            overall_source_breadth=defaultdict(list)
            informative_distinct_genera=set()
            observed_shared=set()
            for id in ids:
                groups=d[id]
                names={label:groups[label]&state_names for label in (
                    "native_nonendemic","endemic","introduced",
                    "native_endemism_unresolved")}
                names["mainland_not_recorded"]=byisland[id]&state_names
                names["native_nonendemic_GIFT_source_shared"]=(
                    names["native_nonendemic"] & byisland_mainland_shared[id])
                for label,x in names.items():
                    support[label+"_islands_with_any_Direct"]+=int(len(x)>0)
                    support[label+"_islands_Direct_n5"]+=int(len(x)>=5)
                    support[label+"_island_x_plant_Direct_records"]+=len(x)
                    status_species[label].update(x)
                    if x:
                        overall_source_breadth[label].append(float(np.mean([
                            source_count.get(s,0) for s in x])))
                source_only=names["mainland_not_recorded"]
                nnon=names["native_nonendemic"]
                nnon_mainland=names["native_nonendemic_GIFT_source_shared"]
                e=names["endemic"]
                intro=names["introduced"]
                block=region_geo.at[id,"spatial_block"]
                block=str(block) if pd.notna(block) else None
                keysource=near.get(id,())
                cluster=dict(block=block,first_source=keysource[0] if keysource else None)
                for label,one in (
                    ("native_nonendemic_minus_mainland_not_recorded",nnon),
                    ("native_nonendemic_GIFT_source_shared_minus_mainland_not_recorded",
                     nnon_mainland),
                    ("introduced_minus_mainland_not_recorded",intro)):
                    if len(one)>=MIN_DIRECT_ISLAND and len(source_only)>=MIN_DIRECT_SOURCE:
                        records[label].append(dict(**cluster,
                            difference=float(np.mean([trait[s] for s in one])-
                                    np.mean([trait[s] for s in source_only]))))
                if len(e)>=MIN_DIRECT_ISLAND and len(nnon)>=MIN_DIRECT_ISLAND:
                    records["endemic_minus_native_nonendemic_raw"].append(
                        dict(**cluster,difference=float(
                            np.mean([trait[s] for s in e])-
                            np.mean([trait[s] for s in nnon]))))
                m=original_shared_genus(e,nnon,trait)
                support["n_same_island_matched_genera_total"]+=m[
                    "n_matched_same_island_genus_units"]
                support["n_same_island_matched_genera_with_category_variation"]+=m[
                    "n_matched_same_island_genus_with_actual_category_variation"]
                support["n_islands_with_at_least_one_same_island_matched_genus"]+=int(
                    m["n_matched_same_island_genus_units"]>0)
                support["n_islands_with_informative_matched_status_genus"]+=int(
                    m["n_matched_same_island_genus_with_actual_category_variation"]>0)
                if m["n_matched_same_island_genus_units"]:
                    # Each genus equal within island; island equal between islands.
                    records["endemic_minus_native_nonendemic_same_island_genus"].append(
                        dict(**cluster,
                            difference=m[
                                "mean_genus_equal_endemic_minus_native_nonendemic_share"]))
                    for g in set(genus(s) for s in e)&set(genus(s) for s in nnon):
                        if any(s in trait for s in e if genus(s)==g) and (
                            any(s in trait for s in nnon if genus(s)==g)):
                            observed_shared.add(g)
                            trait_values={
                                trait[s] for s in (e|nnon)
                                if s in trait and genus(s)==g}
                            if len(trait_values)>=2:
                                informative_distinct_genera.add(g)
            entry=dict(region=region,trait_sector=endpoint,
                       n_original_region_floras=len(ids),
                       original_source_three_status_groups=[
                           dict(status=r,
                               n_original_distinct_Direct_accepted_species=len(status_species[r]),
                               n_original_island_x_plant_Direct_records=support[
                                   r+"_island_x_plant_Direct_records"],
                               n_original_islands_with_n5_Direct_carriers=support[
                                   r+"_islands_Direct_n5"],
                               mean_candidate_GIFT_mainland_entity_count_per_recorded_species=(
                                   float(np.mean(overall_source_breadth[r]))
                                   if overall_source_breadth[r] else None))
                           for r in status_species],
                       n_original_island_genus_matched_endemic_native_units=support[
                           "n_same_island_matched_genera_total"],
                       n_original_same_island_genus_endemic_native_units_with_ANY_true_category_variation=support[
                           "n_same_island_matched_genera_with_category_variation"],
                       n_original_islands_with_informative_genus_match=support[
                           "n_islands_with_informative_matched_status_genus"],
                       n_original_distinct_shared_endemic_nonendemic_genera=len(observed_shared),
                       n_original_distinct_shared_endemic_nonendemic_genera_with_category_variation=len(informative_distinct_genera))
            census.append(entry)
            for label,rec in records.items():
                block=grouped_interval(rec,"block",boots=boots,
                    seed=int.from_bytes(hashlib.sha256((
                        region+endpoint+label+"block").encode()).digest()[:8],"big"))
                source_grp=grouped_interval(rec,"first_source",boots=boots,
                    seed=int.from_bytes(hashlib.sha256((
                        region+endpoint+label+"source").encode()).digest()[:8],"big"))
                if "same_island_genus" in label and (
                        entry["n_original_distinct_shared_endemic_nonendemic_genera_with_category_variation"]<10):
                    block=dict(status="NOT_IDENTIFIABLE_TOO_FEW_DISTINCT_VARIABLE_STATUS_GENERA",
                               n_comparable_floras=len(rec),
                               n_source_clusters=block.get("n_source_clusters",0),
                               n_original_distinct_variable_genera=entry[
                                   "n_original_distinct_shared_endemic_nonendemic_genera_with_category_variation"])
                    source_grp=dict(status="NOT_IDENTIFIABLE_TOO_FEW_DISTINCT_VARIABLE_STATUS_GENERA",
                               n_comparable_floras=len(rec),
                               n_source_clusters=source_grp.get("n_source_clusters",0),
                               n_original_distinct_variable_genera=entry[
                                   "n_original_distinct_shared_endemic_nonendemic_genera_with_category_variation"])
                analyses.append(dict(
                    region=region,trait_sector=endpoint,contrast=label,
                    n_original_islands_with_both_Direct_groups=len(rec),
                    n_original_genus_informative_matched_source_units=(
                        entry["n_original_same_island_genus_endemic_native_units_with_ANY_true_category_variation"]
                        if "same_island_genus" in label else None),
                    n_original_DISTINCT_status_matched_variable_genera=(
                        entry["n_original_distinct_shared_endemic_nonendemic_genera_with_category_variation"]
                        if "same_island_genus" in label else None),
                    original_physical_block_equal=block,
                    original_first_GIFT_source_equal=source_grp))
    if sum(byregion.values())!=4505 or len(analyses)!=len(byregion)*len(vocab)*5:
        raise AssertionError("Original 4505 source three-stage population not complete")
    return dict(
        status="ORIGINAL_DIRECT_THREE_STAGE_FLORAL_STATUS_SOURCE_COMPARISONS",
        n_original_botanical_island_floras=4505,
        n_original_source_direct_sectors=len(vocab),
        n_original_geographic_regimes_including_unresolved=len(byregion),
        original_region_population=[
            dict(region=k,n_original_botanical_floras=int(v))
            for k,v in sorted(byregion.items())],
        original_endemic_native_introduced_source_census=status_counts,
        original_GIFT_five_source_opportunity_census=coverage,
        original_trait_status_source_opportunity=census,
        original_direct_three_stage_comparisons=analyses,
        claims=dict(
            original_all_4505_botanical_island_floras_preserved=True,
            native_endemic_label_requires_concordant_original_three_status_fields=True,
            no_WCVP_regional_presence_upgraded_to_local_endemism=True,
            introduced_origin_is_original_source_not_confirmed_naturalization=True,
            source_only_species_are_absent_from_recorded_full_flora_not_true_failed_arrivals=True,
            strict_native_nonendemic_GIFT_source_shared_subset_separately_computed=True,
            potential_nearest_five_GIFT_native_mainland_source_not_seed_rain=True,
            establishment_sorting_raw_effect_NOT_genus_conditioned=True,
            endemic_genus_effect_observed_same_island_both_status_not_sisters=True,
            genus_invariant_zero_not_evidence_against_adaptive_evolution=True,
            same_island_matched_endemic_inference_requires_at_least_10_distinct_variable_genera=True,
            biological_source_breadth_proxy_is_GIFT_reporting_not_true_range=True,
            source_species_and_traits_are_literature_wide_not_local_population_phenotypes=True,
            original_any_source_POLLINATOR_service_or_seed_set_NOT_observed=True,
            numerical_resamples_not_evolutionarily_independent_island_histories=True,
            no_original_private_island_taxon_genus_GIFT_source_id_coordinates_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):h.update(piece)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    for key in ("original-root","geology","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    a=parser.parse_args()
    original=dict(
        traits=a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        flora=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        mainland=a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        assign=a.original_root/"fixed/source/result/island_source_assignments.csv.gz")
    signatures=dict(
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        mainland="3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")
    for k,path in original.items():
        if sha(path)!=signatures[k]:
            raise ValueError("Original protected GIFT source file SHA256 changed: "+k)
    states=strict_traits(pd.read_csv(original["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    geography=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block","geological_origin"])
    donors=pd.read_csv(original["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    candidate=pd.read_csv(original["assign"],dtype=str,usecols=[
        "island_id","entity_ID","source_mode","source_rank"])
    if geography.island_id.nunique()!=8264 or donors.entity_ID.nunique()!=748:
        raise ValueError("Independent GSHHG/GIFT source baseline changed")
    answer=analyze(states,flora,geography,candidate,donors)
    if answer["n_original_botanical_island_floras"]!=4505:
        raise ValueError("Original botanical 4505 population drifted")
    output=json.dumps(answer,indent=2,sort_keys=True,allow_nan=False)
    for x in ('"island_id"','"accepted_species"','"work_species"',
              '"entity_ID"','"spatial_block"','"source_signature"',
              '"island_latitude"','"island_longitude"'):
        if x in output:raise ValueError("Private source GIFT flora key exported")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(output+"\n")
    print(json.dumps(dict(
        status=answer["status"],
        n_original_botanical_island_floras=4505,
        n_original_status_response_panels=len(answer[
           "original_direct_three_stage_comparisons"]),
        no_propogule_arrival_or_adaptive_fitness_measured=True)))


if __name__=="__main__":main()
