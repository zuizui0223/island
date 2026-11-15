"""Direct original island endemic versus GIFT source-congeneric phenotype states.

Not a genus predictive model: the signed source contrast is the actually
documented original endemic or native-nonendemic plant's category 0/1
minus the mean actual CATEGORY of OTHER independent accepted species in
the SAME documented genus in that island's five candidate GIFT native
mainland source checklists. Exclude ALL plant names in the focal island's
FULL recorded botanical flora from source candidate lists.

True geographic source lists are NOT verified ancestors, arrivals or
biologically independent mainland morphometry. Original species-wide Direct
labels are NOT island-population phenotype evolution. The whole 4505
botanical source population and original native/endemic/intro status are
kept; never require three floral traits on every species.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    status_groups,source_sectors,grouped_interval,genus
)

REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
SCOPES=("endemic","native_nonendemic")
MIN_FOCAL_PLANTS_PER_ISLAND=3
N_BOOT=999


def direct_plant_vs_actual_source_congeners(focal:set[str],donors:set[str],
                                             phenotype:dict[str,int]):
    """Actual 0/1 source congener state distance for EACH original plant name."""
    if focal&donors:
        raise ValueError("Focal plant name appears in independent source candidate")
    if not focal<=set(phenotype) or not donors<=set(phenotype):
        raise ValueError("Original high/medium species-Direct labels missing")
    by_genus=defaultdict(list)
    for sp in donors:
        by_genus[genus(sp)].append(phenotype[sp])
    rows=[]
    for plant in sorted(focal):
        reference=by_genus.get(genus(plant),[])
        if not reference:continue
        x=phenotype[plant]
        rows.append(dict(
            _private_accepted_plant_name=plant,
            _private_documented_genus=genus(plant),
            real_source_congener_n=len(reference),
            real_source_congener_has_variable_categories=len(set(reference))>=2,
            signed_plant_minus_actual_congener_category=float(x-np.mean(reference)),
            absolute_plant_congener_discordance=float(np.mean(
                [int(x!=v) for v in reference])),
        ))
    return rows


def analyze(states:dict,original:pd.DataFrame,geo:pd.DataFrame,
            assignment:pd.DataFrame,mainland:pd.DataFrame,
            *,n_boot=N_BOOT):
    # Independently original native/endemic/floristic status agreed.
    grouped=status_groups(original)
    if len(grouped)!=4505:raise ValueError("Full 4505 source island floras changed")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or (
        geo.island_id.duplicated().any()):
        raise ValueError("Original GSHHG island source region/block invalid")
    if not {"island_id","entity_ID","source_mode","source_rank"}<=set(assignment):
        raise ValueError("Original island candidate geographic source rank invalid")
    if not {"entity_ID","work_species"}<=set(mainland) or mainland.duplicated(
            ["entity_ID","work_species"]).any():
        raise ValueError("Original GIFT native source accepted plant names duplicate")
    g=geo.set_index("island_id").reindex(sorted(grouped))
    if g.analysis_regime.notna().sum()!=4504:
        raise ValueError("Original one geography-unresolved island status changed")
    region_of={
        i:(str(g.at[i,"analysis_regime"])
           if pd.notna(g.at[i,"analysis_regime"]) else "geography_unresolved")
        for i in grouped}
    region_counts=Counter(region_of.values())
    if sum(region_counts.values())!=4505 or set(region_counts)!=set(REGIONS):
        raise ValueError("Original GSHHG region universe no longer complete")
    src=assignment.loc[assignment.source_mode.eq("geo_k5")].copy()
    src.source_rank=pd.to_numeric(src.source_rank,errors="raise")
    src=src.loc[src.source_rank.between(1,5)]
    if src.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original GIFT nearest-mainland rank")
    bysource=src.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    flora_mainland=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    labels=source_sectors(states)
    direct_names=set().union(*(set(p) for p in labels.values()))
    candidates={}
    for i in grouped:
        sources=bysource.get(i,())
        candidates[i]=(set().union(*(flora_mainland.get(x,set()) for x in sources))
                       -grouped[i]["full"])&direct_names
    panels=[]
    for region in REGIONS:
        island_ids=[i for i in sorted(grouped) if region_of[i]==region]
        for sector,phenotype in labels.items():
            all_labelled=set(phenotype)
            collected={}
            for status in SCOPES:
                collected[status]=[]
            comparisons=[]
            visited=defaultdict(set)
            source_control=Counter()
            for i in island_ids:
                donor=candidates[i]&all_labelled
                bystatus={}
                source_focal_names=grouped[i]
                for role in SCOPES:
                    focal=source_focal_names[role]&all_labelled
                    source_control[role+"_n_direct_species_island_rows"]+=len(focal)
                    source_control[role+"_n_islands_with_any_direct_carriers"]+=int(bool(focal))
                    # Distinct accepted species not source congener opportunity per se.
                    records=direct_plant_vs_actual_source_congeners(
                        focal,donor,phenotype)
                    bystatus[role]=records
                    source_control[role+"_n_island_by_plant_congener_rows"]+=len(records)
                    for r in records:
                        visited[role].add(r["_private_accepted_plant_name"])
                        visited[role+"_genus"].add(r["_private_documented_genus"])
                    if len(records)>=MIN_FOCAL_PLANTS_PER_ISLAND:
                        block=g.at[i,"spatial_block"]
                        near=bysource.get(i,())
                        collected[role].append(dict(
                            spatial_block=str(block) if pd.notna(block) else None,
                            first_source=near[0] if near else None,
                            difference=float(np.mean([
                                x["signed_plant_minus_actual_congener_category"]
                                for x in records])),
                            discordance=float(np.mean([
                                x["absolute_plant_congener_discordance"]
                                for x in records])),
                            n_focal_Direct_plant_congener_rows=len(records),
                            n_direct_original_trait_plant_names=len(focal),
                            n_variable_source_congener_records=sum(
                                r["real_source_congener_has_variable_categories"]
                                for r in records),
                            mean_n_original_congeners=float(np.mean([
                                r["real_source_congener_n"] for r in records])),
                            island_identity_private=i))
                e=collected["endemic"][-1] if (
                    collected["endemic"] and
                    collected["endemic"][-1]["island_identity_private"]==i) else None
                n=collected["native_nonendemic"][-1] if (
                    collected["native_nonendemic"] and
                    collected["native_nonendemic"][-1]["island_identity_private"]==i) else None
                if e is not None and n is not None:
                    comparisons.append(dict(
                        spatial_block=e["spatial_block"],first_source=e["first_source"],
                        difference=e["difference"]-n["difference"],
                        n_endemic_source_congener_plant_rows=e["n_focal_Direct_plant_congener_rows"],
                        n_native_nonendemic_source_congener_plant_rows=n[
                            "n_focal_Direct_plant_congener_rows"]))
            report=dict(
                region=region,sector=sector,
                n_original_botanical_floras_in_region=len(island_ids),
                source_only_candidate_mainland_species_not_ancestors=True)
            for status,records in collected.items():
                z=dict(
                    n_original_islands_with_3plus_signed_name_distinct_source_congeners=len(records),
                    n_original_island_x_plant_original_Direct_records=source_control[
                        status+"_n_direct_species_island_rows"],
                    n_original_island_x_plant_name_disjoint_source_congener_qualified_records=source_control[
                        status+"_n_island_by_plant_congener_rows"],
                    n_original_distinct_focal_accepted_source_congener_qualified_species=len(
                        visited[status]),
                    n_original_distinct_focal_source_qualified_genera=len(
                        visited[status+"_genus"]),
                    n_qualified_focal_islands_with_at_least_one_variable_source_genus=sum(
                        int(r["n_variable_source_congener_records"]>0) for r in records),
                    mean_original_focal_direct_species_carriers_with_source_congeners_per_qualifying_island=(
                        float(np.mean([r["n_focal_Direct_plant_congener_rows"] for r in records]))
                        if records else None))
                if records:
                    for response,key in (
                        ("difference","signed_actual_focal_minus_source_congener"),
                        ("discordance","actual_categorical_mismatch_to_source_congeners")):
                        for group,cluster in (("spatial_block","physical_10degree_block"),
                                              ("first_source","first_GIFT_candidate_mainland_source")):
                            input_data=[dict(**{group:r[group]},difference=r[response])
                                        for r in records]
                            z[key+"__"+cluster]=grouped_interval(
                                input_data,group,boots=n_boot,
                                seed=int.from_bytes(hashlib.sha256(
                                    (region+sector+status+response+cluster).encode()
                                ).digest()[:8],"big"))
                report[status+"_versus_name_distinct_same_genus_GIFT_candidate_sources"]=z
            for grouping,label in (("spatial_block","physical_10degree_block"),
                                   ("first_source","first_GIFT_candidate_mainland_source")):
                report["endemic_minus_native_nonendemic_congener_centered__"+label]=grouped_interval(
                    comparisons,grouping,boots=n_boot,
                    seed=int.from_bytes(hashlib.sha256(
                        (region+sector+"endemic_vs_native_direct"+label).encode()
                    ).digest()[:8],"big"))
            report["n_original_islands_with_endemic_AND_native_nonendemic_3plus_source_congener_plants"]=len(comparisons)
            panels.append(report)
    if len(panels)!=len(REGIONS)*len(labels):
        raise AssertionError("Original source trait geography omitted")
    return dict(
        status="ORIGINAL_4505_ENDEMIC_SAME_GENUS_CANDIDATE_MAINLAND_STATE_CONTRAST",
        n_original_botanical_island_floras=4505,
        n_original_Direct_traits_sectors=len(labels),
        original_complete_flora_GSHHG_geography=[
            dict(region=region,n_original_botanical_island_floras=int(region_counts[region]))
            for region in REGIONS],
        original_signed_endemic_and_native_same_genus_mainland_source_panels=panels,
        claims=dict(
            origin_endemic_nonendemic_from_three_concordant_archived_labels=True,
            original_focal_entire_island_accepted_plant_names_excluded_from_donors=True,
            comparisons_use_actual_source_congener_species_direct_state_not_predictive_accuracy=True,
            actual_signed_focal_minus_same_genus_source_direct_trait_not_ancestor_state=True,
            same_source_genus_donor_presence_not_independent_mainland_population_phenotype=True,
            genus_category_invariant_zero_is_not_no_local_adaptation=True,
            raw_4505_global_flora_not_n30_triple_complete_restricted=True,
            no_phylogenetic_independent_sisters_history_pollen_fitness_or_guild_loss=True,
            original_private_taxon_island_genus_source_entity_coordinates_not_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):
            h.update(data)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    paths={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "mainland":a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    expected={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for k,path in paths.items():
        if sha(path)!=expected[k]:
            raise ValueError("Original accepted plant Direct source SHA mismatch "+k)
    states=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    candidates=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(paths["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    if geo.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original GSHHG geography and GIFT mainland entity drift")
    out=analyze(states,flora,geo,candidates,mainland)
    if out["n_original_botanical_island_floras"]!=4505:
        raise ValueError("Original 4505 source botanical population lost")
    blob=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    for bad in ('"island_id"','"accepted_species"','"entity_ID"',
                '"work_species"','"spatial_block"','"source_signature"',
                '"island_latitude"','"island_longitude"'):
        if bad in blob:
            raise ValueError("Rights-limited original source botanical data exported")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(blob+"\n")
    print(json.dumps(dict(status=out["status"],
                          n_original_floras=4505,
                          n_source_trait_region_panels=len(
                              out["original_signed_endemic_and_native_same_genus_mainland_source_panels"]),
                          no_ancestral_adaptation_claim=True)))


if __name__=="__main__":main()
