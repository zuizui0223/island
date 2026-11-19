"""Direct original island floristic composition vs native GIFT source breadth.

Outcomes are ACTUAL category-frequency differences, not prediction gains.
Within each original island, strict native-nonendemic accepted species that
are themselves GIFT mainland source listed are compared to distinct GIFT
native mainland candidates NOT recorded in the island's FULL original
flora. Match plant names on GIFT source-listing breadth (count of original
native mainland GIFT region entities), not botanical genus, to preserve
whole-flora genus sorting. This is a documentation/geographic range PROXY,
not seed dispersal, propagule pressure or proven establishment.

Also compare introduced vs strict native-nonendemic category fractions on
the same original island; introductions may be horticultural and their
true naturalized population status cannot be inferred here.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    status_groups,source_sectors,grouped_interval,
)

BINS=("source_once","source_2_3","source_4_10","source_11_30","source_31plus")
REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
MIN_SOURCE_MATCHED_EFFECTIVE=5
BOOTS=999


def listing_breadth_bin(n:int)->str:
    if not isinstance(n,(int,np.integer)) or n<1:
        raise ValueError("Original GIFT native candidate must be source listed")
    if n==1:return "source_once"
    if n<=3:return "source_2_3"
    if n<=10:return "source_4_10"
    if n<=30:return "source_11_30"
    return "source_31plus"


def breadth_matched_sector_difference(
        native:list[str]|set[str],external:list[str]|set[str],
        state:dict[str,int],listed_regions:dict[str,int]):
    """Equal within-bin contrast; min(counts) weights without replacement.

    Actual recorded candidate plants are never re-labelled if the other
    origin/breadth group lacks adequate independent named source coverage.
    """
    if not native or not external:
        return None
    if set(native)&set(external):
        raise ValueError("Recorded island and unrecorded source plants overlap")
    if not set(native)|set(external)<=set(state):
        raise ValueError("Not independently direct recorded on both sets")
    pools={r:{k:[] for k in BINS} for r in ("native","mainland_only")}
    for role,sp in (("native",s) for s in native):
        pools[role][listing_breadth_bin(listed_regions[sp])].append(state[sp])
    for role,sp in (("mainland_only",s) for s in external):
        pools[role][listing_breadth_bin(listed_regions[sp])].append(state[sp])
    sums=0.
    m=0
    used=0
    for name in BINS:
        one,two=pools["native"][name],pools["mainland_only"][name]
        w=min(len(one),len(two))
        if w:
            m+=w
            used+=1
            sums+=w*(float(np.mean(one))-float(np.mean(two)))
    return dict(
        n_original_breadth_overlapping_bin_support=used,
        n_effectively_matched_source_listing_breadth_carriers=m,
        fraction_native_carriers_from_overlapping_source_breadth_bins=(
            sum(len(pools["native"][name]) for name in BINS
                if pools["native"][name] and pools["mainland_only"][name])/len(native)),
        fraction_mainland_candidate_carriers_from_overlap_breadth_bins=(
            sum(len(pools["mainland_only"][name]) for name in BINS
                if pools["native"][name] and pools["mainland_only"][name])/len(external)),
        source_breadth_matched_native_minus_candidate_sector_share=(
            sums/m if m else None),
        same_island_original_unmatched_native_minus_candidate_sector_share=(
            float(np.mean([state[s] for s in native])-
                  np.mean([state[s] for s in external]))),
    )


def analyze(state:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            source_assign:pd.DataFrame,mainland:pd.DataFrame,
            *,nboot=BOOTS):
    if not {"island_id","accepted_species","origin_status",
            "endemic_status","floristic_status"}<=set(flora):
        raise ValueError("Original 3 independent botanical statuses required")
    full=status_groups(flora)
    if len(full)!=4505:
        raise ValueError("4505 original accepted source botanical island floras changed")
    if not {"island_id","spatial_block","analysis_regime"}<=set(geo) or (
        geo.island_id.duplicated().any()):
        raise ValueError("Independent original GSHHG geography invalid")
    if not {"island_id","source_rank","source_mode","entity_ID"}<=set(source_assign):
        raise ValueError("Original native GIFT mainland source ranking missing")
    if not {"entity_ID","work_species"}<=set(mainland) or mainland.duplicated(
            ["entity_ID","work_species"]).any():
        raise ValueError("Original GIFT source region×plant membership invalid")
    g=geo.set_index("island_id").reindex(sorted(full))
    if g.analysis_regime.notna().sum()!=4504:
        raise ValueError("Original one geography unresolved GIFT flora changed")
    a=source_assign.loc[source_assign.source_mode.eq("geo_k5")].copy()
    a.source_rank=pd.to_numeric(a.source_rank,errors="raise")
    a=a.loc[a.source_rank.between(1,5)]
    if a.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original candidate mainland GIFT source rank duplicate")
    byisland=a.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    byentity=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    count=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    measured=source_sectors(state)
    regions={
        i:(str(g.at[i,"analysis_regime"]) if pd.notna(g.at[i,"analysis_regime"])
           else "geography_unresolved") for i in full
    }
    categories=Counter(regions.values())
    if set(categories)!=set(REGIONS) or sum(categories.values())!=4505:
        raise ValueError("Original named independent 4505 GSHHG geography drift")
    census=[]
    results=[]
    for region in REGIONS:
        ids=[i for i in sorted(full) if regions[i]==region]
        for trait,scored in measured.items():
            labels=set(scored)
            native_matched_records=[]
            introduced_records=[]
            support=Counter()
            for i in ids:
                q=full[i]
                nearest=byisland.get(i,())
                source=set().union(*(byentity.get(e,set()) for e in nearest))
                outside=(source-q["full"])&labels
                native_shared=(source&q["native_nonendemic"])&labels
                introduced=q["introduced"]&labels
                native_original=q["native_nonendemic"]&labels
                support["n_original_source_candidate_flora_with_5_nearest"]+=int(
                    len(nearest)==5)
                support["n_strict_native_shared_direct_carrier_records"]+=len(native_shared)
                support["n_candidate_mainland_only_direct_carrier_records"]+=len(outside)
                support["n_original_introduced_direct_carrier_records"]+=len(introduced)
                support["n_original_strict_native_nonendemic_direct_carrier_records"]+=len(native_original)
                block=g.at[i,"spatial_block"]
                block=str(block) if pd.notna(block) else None
                first=nearest[0] if nearest else None
                cluster=dict(spatial_block=block,first_source=first)
                if len(native_shared)>=5 and len(outside)>=10:
                    support["n_original_islands_with_raw_GIFT_listed_native_and_source_eligible"]+=1
                    z=breadth_matched_sector_difference(
                        native_shared,outside,scored,count)
                    if z and z["n_effectively_matched_source_listing_breadth_carriers"]>=MIN_SOURCE_MATCHED_EFFECTIVE:
                        native_matched_records.append(dict(
                            **cluster,
                            breadth_overlap_weight=z["n_effectively_matched_source_listing_breadth_carriers"],
                            n_original_shared_breadth_bins=z["n_original_breadth_overlapping_bin_support"],
                            original_native_carrier_overlap_fraction=z[
                                "fraction_native_carriers_from_overlapping_source_breadth_bins"],
                            original_mainland_carrier_overlap_fraction=z[
                                "fraction_mainland_candidate_carriers_from_overlap_breadth_bins"],
                            raw=z["same_island_original_unmatched_native_minus_candidate_sector_share"],
                            matched=z["source_breadth_matched_native_minus_candidate_sector_share"],
                            delta=z["source_breadth_matched_native_minus_candidate_sector_share"]-
                                  z["same_island_original_unmatched_native_minus_candidate_sector_share"]))
                if len(introduced)>=5 and len(native_original)>=5:
                    introduced_records.append(dict(
                        **cluster,
                        raw=float(np.mean([scored[x] for x in introduced])-
                                  np.mean([scored[x] for x in native_original]))))
            res=dict(region=region,trait_sector=trait,
                n_original_region_botanical_floras=len(ids),
                n_strict_native_shared_direct_carrier_records=support[
                    "n_strict_native_shared_direct_carrier_records"],
                n_candidate_mainland_only_direct_carrier_records=support[
                    "n_candidate_mainland_only_direct_carrier_records"],
                n_original_introduced_direct_carrier_records=support[
                    "n_original_introduced_direct_carrier_records"],
                n_original_strict_native_nonendemic_direct_carrier_records=support[
                    "n_original_strict_native_nonendemic_direct_carrier_records"],
                n_original_islands_with_raw_GIFT_listed_native_and_source_eligible=support[
                    "n_original_islands_with_raw_GIFT_listed_native_and_source_eligible"],
                n_original_breadth_matched_islands=len(native_matched_records),
                n_original_same_island_introduced_vs_nonendemic_floras=len(
                    introduced_records))
            for comp,records,fields in (
                ("original_source_breadth_matched_native_vs_source",native_matched_records,
                 ("raw","matched","delta")),
                ("original_introduced_vs_strict_native_nonendemic",introduced_records,
                 ("raw",))):
                output=dict()
                if records:
                    table=pd.DataFrame(records)
                    if comp=="original_source_breadth_matched_native_vs_source":
                        output.update(
                            mean_effective_source_listing_breadth_match_count=float(
                                table.breadth_overlap_weight.mean()),
                            n_original_islands_with_more_than_one_common_breadth_bin=int(
                                (table.n_original_shared_breadth_bins>=2).sum()),
                            mean_native_carrier_fraction_with_candidate_breadth_support=float(
                                table.original_native_carrier_overlap_fraction.mean()),
                            mean_mainland_candidate_carrier_fraction_with_native_breadth_support=float(
                                table.original_mainland_carrier_overlap_fraction.mean()))
                    for key in fields:
                        for grouping,label in (("spatial_block","physical_10degree_block"),
                                               ("first_source","first_GIFT_source")):
                            data=[dict(**{grouping:r[grouping]},difference=r[key])
                                  for r in records]
                            output[key+"__"+label]=grouped_interval(
                                data,grouping,boots=nboot,
                                seed=int.from_bytes(hashlib.sha256(
                                    (region+trait+comp+key+label).encode()
                                ).digest()[:8],"big"))
                res[comp]=output or dict(
                    status="NOT_IDENTIFIABLE_NO_DIRECT_SOURCE_OVERLAP")
            census.append(res)
    return dict(
        status="ORIGINAL_4505_DIRECT_NATIVE_SOURCE_BREADTH_MATCHING",
        n_original_botanical_floras=4505,
        n_original_independently_observed_GIFT_native_source_entities=mainland.entity_ID.nunique(),
        n_original_real_Direct_morphology_and_reproductive_sectors=len(measured),
        source_breadth_cutoffs=list(BINS),
        original_botanical_GSHHG_population=[
            dict(region=r,n_original_floras=int(categories[r])) for r in REGIONS],
        original_floristic_breadth_and_introduction_contrasts=census,
        claims=dict(
            original_entire_botanical_island_flora_retained=True,
            strictly_native_nonendemic_requires_original_all_three_status_agreement=True,
            GIFT_breadth_counts_original_DISTINCT_listed_mainland_regions_not_range_area=True,
            no_genus_matching_in_first_stage_candidate_flora_sorting_contrast=True,
            both_native_and_source_names_MUST_be_listed_in_GIFT_candidate_mainland_source=True,
            source_only_candidate_absent_from_full_recorded_island_flora_not_proven_to_fail_arrival=True,
            within_island_breadth_overlap_weighted_not_replacement_or_new_species=True,
            introduced_group_not_independent_confirmed_naturalized_population=True,
            no_ignored_original_genus_same_trait_ambiguous_label_imputation=True,
            no_independent_phylogeny_pollinator_visits_fitness_or_selection_claim=True,
            original_private_botanical_island_taxon_GIFT_entity_ids_never_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    paths={
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "mainland":a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    expected={
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for key,filename in paths.items():
        if sha(filename)!=expected[key]:
            raise ValueError("Original GIFT full-4505 source provenance SHA mismatch: "+key)
    states=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assignment=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(paths["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    if geo.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original GSHHG and GIFT candidate region reference changed")
    result=analyze(states,flora,geo,assignment,mainland)
    if result["n_original_botanical_floras"]!=4505 or (
         sum(r["n_original_floras"] for r in result["original_botanical_GSHHG_population"])!=4505):
        raise ValueError("Original all islands lost from the result universe")
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for forbidden in ('"island_id"','"accepted_species"','"work_species"',
                      '"entity_ID"','"spatial_block"','"source_signature"',
                      '"island_latitude"','"island_longitude"'):
        if forbidden in body:
            raise ValueError("Original accepted botany/source ID escaped protected CI")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(status=result["status"],n_source_floras=4505,
                          n_sector_panels=len(result[
                              "original_floristic_breadth_and_introduction_contrasts"]),
                          no_independent_colonization_or_adaptation_inference=True)))


if __name__=="__main__":main()
