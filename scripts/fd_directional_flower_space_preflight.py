"""Original 4505 island source Direct flower-category directional isolation support audit.

No ecological hypotheses are tested by this preflight; count recorded
categories and exact region×isolation observation opportunities. This freezes
the original label dictionary BEFORE running directional source gradients.
Original island IDs, species names and coordinates remain in CI.
"""
from __future__ import annotations
import argparse
import json
import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits

CHECK_TRAITS=(
    "floral_form","floral_symmetry","tube_depth_class",
    "flower_size_class","autonomous_selfing_capacity","self_incompatibility",
    "flower_primary_color")
PAIRS=(("floral_form","floral_symmetry"),
       ("floral_symmetry","tube_depth_class"),
       ("floral_form","tube_depth_class"))
MIN_DIRECT=5


def preflight(states:dict,flora:pd.DataFrame,geo:pd.DataFrame)->dict:
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Incomplete original accepted flora")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original botanical flora duplicates")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Invalid independent original native status")
    fields={"island_id","analysis_regime","spatial_block",
            "island_latitude","distance_to_continent_km",
            "area_km2","climate_pc1","climate_pc2","climate_pc3","climate_pc4"}
    if not fields<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Independent GSHHG geography incomplete")
    if not set(CHECK_TRAITS)<=set(states):
        raise ValueError("Missing original real Direct trait axes")
    allmap=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    if not allmap:raise ValueError("No original GIFT botanical floras")
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    islands=sorted(allmap)
    g=geo.set_index("island_id").reindex(islands)
    distance=pd.to_numeric(g.distance_to_continent_km,errors="coerce")
    latitude=pd.to_numeric(g.island_latitude,errors="coerce")
    area=pd.to_numeric(g.area_km2,errors="coerce")
    region=g.analysis_regime.fillna("geography_unresolved").astype(str)
    block=g.spatial_block
    valid=(np.isfinite(distance)&(distance>=0)&
           np.isfinite(latitude)&(latitude.abs()<=90)&
           np.isfinite(area)&(area>0)&block.notna()&
           region.ne("geography_unresolved"))
    trait_names={}
    for t in CHECK_TRAITS:
        counts=Counter(states[t].values())
        trait_names[t]=dict(
            n_original_unique_Direct_species=len(states[t]),
            original_Direct_category_counts=[
                dict(category=str(name),n_unique_Direct_source_species=int(n))
                for name,n in sorted(counts.items(),key=lambda kv:str(kv[0]))])
    additional={}
    for a,b in PAIRS:
        names=set(states[a])&set(states[b])
        pair_ct=Counter((states[a][s],states[b][s]) for s in names)
        additional[a+"__"+b]=dict(
            n_exact_same_original_Direct_species=len(names),
            n_observed_joint_category_states=len(pair_ct),
            original_observed_category_combinations=[
                dict(first=str(u),second=str(v),n_unique_Direct_source_species=int(n))
                for (u,v),n in sorted(pair_ct.items(),key=lambda kv:repr(kv[0]))])
    support=[]
    for scope,mapping in (("all_recorded",allmap),("original_native",native)):
        for endpoint,known in (
            *[(t,set(states[t])) for t in CHECK_TRAITS],
            *[(a+"__"+b,set(states[a])&set(states[b])) for a,b in PAIRS]):
            eligible=pd.Series({
                island:len(mapping.get(island,set())&known) for island in islands
            })
            covered=(eligible>=MIN_DIRECT)
            for label,selector in (
                ("ALL_4505",pd.Series(True,index=islands)),
                ("northern_midlatitude",region.eq("northern_midlatitude")),
                ("tropical",region.eq("tropical")),
                ("southern_extratropical",region.eq("southern_extratropical")),
                ("highlatitude_ge60_nested",latitude.abs()>=60)):
                pre=covered&selector
                usable=pre&valid
                subset=usable[usable].index
                support.append(dict(
                    scope=scope,endpoint=endpoint,region=label,
                    n_original_botanical_islands=int(selector.sum()),
                    n_islands_direct_n5=int(pre.sum()),
                    n_n5_islands_with_verified_isolation_and_geography=int(usable.sum()),
                    n_spatial_blocks_with_isolation_n5=int(block.reindex(subset).nunique()),
                    n_n10_islands_with_verified_isolation=int((usable&(eligible>=10)).sum()),
                    n_n5_islands_direct_coverage_ge10pct=int(sum(
                        eligible.loc[subset] / np.array([len(mapping.get(i,set()))
                                                       for i in subset],float)>=.10)),
                    median_mainland_isolation_km=(
                        float(distance.loc[subset].median()) if len(subset) else None),
                    min_mainland_isolation_km=(
                        float(distance.loc[subset].min()) if len(subset) else None),
                    max_mainland_isolation_km=(
                        float(distance.loc[subset].max()) if len(subset) else None),
                ))
    return dict(
        status="ORIGINAL_DIRECTIONAL_FLORAL_CATEGORY_PREOUTCOME_SOURCE_CENSUS",
        n_original_source_botanical_islands=len(islands),
        n_original_islands_matched_geography=int(g.analysis_regime.notna().sum()),
        n_original_islands_with_finite_verified_distance_and_area=int(valid.sum()),
        n_original_geographic_isolation_unresolved=int((~valid).sum()),
        original_direct_trait_categories=trait_names,
        original_actual_source_two_trait_combinations=additional,
        directional_isolation_trait_and_region_eligibility=support,
        claims=dict(
            source_island_data_4505_populational_frame=True,
            source_category_states_recorded_not_island_population_morphometry=True,
            highlatitude_is_overlapping_nested_sensitivity=True,
            target_region_cutoffs_and_traits_fixed_before_gradient=True,
            missing_original_trait_categories_never_imputed_zero=True,
            no_isolation_slope_or_pollination_adaptation_inference_in_preflight=True,
            no_private_original_plant_island_GIFT_source_or_coordinates_exported=True))


def sha(path):
    digest=hashlib.sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda:fp.read(1<<20),b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","output"):
        p.add_argument("--"+k,type=Path,required=True)
    args=p.parse_args()
    paths=dict(
        traits=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        flora=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz")
    expected=dict(
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a")
    for name,path in paths.items():
        if sha(path)!=expected[name]:
            raise ValueError("Protected source changed: "+name)
    states=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(args.geology,dtype=str)
    if flora.island_id.nunique()!=4505 or geo.island_id.nunique()!=8264:
        raise ValueError("Original GIFT 4505 / independently mapped 8264 changed")
    result=preflight(states,flora,geo)
    if result["n_original_source_botanical_islands"]!=4505:
        raise ValueError("Source original botanical island count drift")
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for token in ('"island_id"','"accepted_species"','"work_species"',
                  '"entity_ID"','"spatial_block"','"island_latitude"',
                  '"island_longitude"','"source_signature"'):
        if token in body:raise ValueError("Original private identity exported")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n")
    print(json.dumps(dict(status=result["status"],
                          n_original_islands=result["n_original_source_botanical_islands"],
                          n_trait_source_vocabulary=len(CHECK_TRAITS),
                          n_paired_source_vocabulary=len(PAIRS),
                          no_biological_isolation_slope_claim=True)))


if __name__=="__main__":main()
