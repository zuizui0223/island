"""Original 4505-island floral reproductive assurance channels as SEPARATE traits.

Self-compatibility (SC vs SI) and documented autonomous selfing capacity
(autonomous vs absent) are not interchangeable. Exclude variable/delayed
category records as explicitly uncertain from strict binary contrasts.
Original high/medium species-direct trait states are globally fixed per
accepted species. Compare island-equal regional isolation slopes and
in-region genus-preserving expected category share, not causal pollen
assurance, selection, or measured plant pollination rates.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_directional_flower_isolation_gradients import (
    weighted_cluster_fit,holm,MIN_DIRECT
)

BINARIES={
    "self_incompatibility":("SC","SI"),
    "autonomous_selfing_capacity":("autonomous","absent")
}
REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "highlatitude_ge60_nested")
SCOPES=("all_recorded","original_native")


def direct_binary_states(states):
    if not set(BINARIES)<=set(states):
        raise ValueError("Original reproductive Direct trait names missing")
    results={}
    audits=[]
    for trait,(yes,no) in BINARIES.items():
        observed=Counter(states[trait].values())
        if yes not in observed or no not in observed:
            raise ValueError("Original binary selfing/SC ontology changed")
        # Every other single-state value is not equated to 'no'.
        clean={sp:int(val==yes) for sp,val in states[trait].items()
               if val in (yes,no)}
        results[trait]=clean
        audits.append(dict(
            trait=trait,yes_category=yes,no_category=no,
            n_original_Direct_category_yes=observed[yes],
            n_original_Direct_category_no=observed[no],
            n_other_direct_status_excluded=sum(n for k,n in observed.items()
                                               if k not in (yes,no)),
            n_distinct_strict_binary_direct_source_species=len(clean)))
    return results,audits


def expected_from_genus(focal_rosters:list[set],regional_rosters:list[set],
                        states:dict[str,int]):
    ref=set().union(*regional_rosters)&set(states)
    g=defaultdict(list)
    for name in ref:
        name_genus=name.split()[0] if len(name.split())>=2 else "fixed:"+name
        g[name_genus].append(states[name])
    p={name:np.mean(x) for name,x in g.items()}
    expected=[]
    for ro in focal_rosters:
        vals=[p[s.split()[0] if len(s.split())>=2 else "fixed:"+s]
              for s in sorted(ro)]
        expected.append(float(np.mean(vals)) if vals else float("nan"))
    return np.asarray(expected)


def analyze(states,flora,geo,assign,*,n_boot=999):
    need={"island_id","accepted_species","origin_status"}
    if not need<=set(flora) or flora.duplicated(
            ["island_id","accepted_species"]).any():
        raise ValueError("Original botanical flora incomplete")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original reproductive plants have unknown nativeness tokens")
    g_need={"island_id","analysis_regime","spatial_block",
            "distance_to_continent_km","area_km2","island_latitude"}
    if not g_need<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Original physical island geography invalid")
    n_original=flora.island_id.nunique()
    if n_original!=4505:raise ValueError("Original all-4505 botanical frame changed")
    B,audit=direct_binary_states(states)
    species_both=set(B["self_incompatibility"])&set(
        B["autonomous_selfing_capacity"])
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    islands=sorted(full)
    g=geo.set_index("island_id").reindex(islands)
    s=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    s.source_rank=pd.to_numeric(s.source_rank,errors="raise")
    if s.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original GIFT mainland source rank")
    first=s.loc[s.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    rows=[]
    shared=[]
    for scope,rosters in (("all_recorded",full),("original_native",native)):
        for region in REGIONS:
            if region=="highlatitude_ge60_nested":
                choices=[i for i in islands if pd.notna(g.at[i,"island_latitude"])
                         and abs(float(g.at[i,"island_latitude"]))>=60]
            else:
                choices=[i for i in islands if g.at[i,"analysis_regime"]==region]
            regional=[rosters.get(i,set()) for i in choices]
            outcome_frames={}
            for trait,label in BINARIES.items():
                t=B[trait]
                fields=[]
                own=[]
                own_ids=[]
                for island in choices:
                    orig=rosters.get(island,set())
                    direct=orig&set(t)
                    if len(direct)<MIN_DIRECT:continue
                    ge=g.loc[island]
                    try:
                        distance=float(ge.distance_to_continent_km)
                        area=float(ge.area_km2)
                        lat=float(ge.island_latitude)
                    except (TypeError,ValueError):continue
                    if (not np.isfinite([distance,area,lat]).all() or
                            distance<0 or area<=0 or abs(lat)>90 or
                            pd.isna(ge.spatial_block)):
                        continue
                    fields.append(dict(
                        response=float(np.mean([t[x] for x in direct])),
                        log2_isolation_plus1=math.log2(1+distance),
                        log1p_area=math.log1p(area),
                        direct_coverage=len(direct)/len(orig),
                        log_direct_n=math.log(len(direct)),
                        absolute_latitude=abs(lat),
                        spatial_block=str(ge.spatial_block),
                        source_cluster=first.get(island),
                        n_binary_Direct=len(direct),
                        n_scope_original_flora_species=len(orig)))
                    own.append(direct)
                    own_ids.append(island)
                outcome=pd.DataFrame(fields,index=own_ids)
                if len(outcome):
                    null=expected_from_genus(own,regional,t)
                    outcome["in_region_genus_relative_response"]=(
                        outcome.response.to_numpy()-null)
                outcome_frames[trait]=outcome
                report=dict(
                    scope=scope,region=region,trait=trait,
                    sector_binary_yes=label[0],sector_binary_no=label[1],
                    n_original_geographic_region_botanical_floras=len(choices),
                    n_original_islands_with_n5_strict_binary_Direct=len(fields),
                    response_refers_to_documented_Direct_strict_binary_carriers_only=True)
                if len(fields)<30:
                    report["status"]="INSUFFICIENT_SOURCE_DIRECT_BINARY_SUPPORT"
                    rows.append(report);continue
                for column,nickname in (("response","observed"),
                                        ("in_region_genus_relative_response","genus_relative")):
                    frame=outcome.copy()
                    frame["response"]=outcome[column]
                    for group,weight in (("spatial_block","physical_block"),
                                         ("source_cluster","first_GIFT_source")):
                        m=weighted_cluster_fit(frame,cluster=group,n_boot=n_boot,
                            seed=int.from_bytes(hashlib.sha256(
                                (scope+region+trait+nickname+weight).encode()).digest()[:8],"big"))
                        for k,v in m.items():report[nickname+"__"+weight+"__"+k]=v
                report["status"]=report.get("observed__physical_block__status")
                rows.append(report)
            # Two different Direct species rosters are NOT the same-plant test.
            both_available=set(outcome_frames["self_incompatibility"].index)&set(
                outcome_frames["autonomous_selfing_capacity"].index)
            shared.append(dict(
                scope=scope,region=region,
                n_islands_with_n5_BOTH_separate_reproductive_trait_rosters=len(both_available),
                n_unique_original_source_species_with_BOTH_strict_binary_labels=len(species_both),
                not_same_species_or_independent_trait_observation_by_island=True))
    return dict(
        status="ORIGINAL_SEPARATE_SELF_COMPATIBILITY_AND_AUTONOMOUS_SELFING_ISOLATION",
        n_original_botanical_islands=n_original,
        n_original_distinct_source_species_with_BOTH_reproductive_categories=len(species_both),
        source_original_strict_binary_ontology=audit,
        original_regional_reproductive_shares_and_genus_isolation=rows,
        original_same_island_separate_trait_documentation_opportunity=shared,
        claims=dict(
            morphology_and_reproduction_are_separate_functional_axes=True,
            self_compatibility_not_equated_to_autonomous_self_fertilization=True,
            delayed_and_variable_categories_excluded_not_recoded_as_absent=True,
            original_direct_species_wide_not_population_reproductive_fitness=True,
            original_island_population_4505_not_forced_three_trait_complete=True,
            source_region_genus_conditioning_is_in_sample_not_adaptive_selection=True,
            no_island_species_GIFT_source_ids_coordinates_or_phylogeny_leaked=True))


def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):
            h.update(z)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","output"):
        p.add_argument("--"+k,required=True,type=Path)
    p.add_argument("--n-boot",type=int,default=999)
    a=p.parse_args()
    if a.n_boot<99:raise ValueError("Original source bootstrap cannot be weakened")
    f={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for k in f:
        if sha(f[k])!=pins[k]:
            raise ValueError("Protected original reproductive source changed")
    original=strict_traits(pd.read_csv(f["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(f["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "distance_to_continent_km","area_km2","island_latitude"])
    source=pd.read_csv(f["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if geo.island_id.nunique()!=8264:raise ValueError("Original geography frame changed")
    result=analyze(original,flora,geo,source,n_boot=a.n_boot)
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for term in ('"island_id"','"accepted_species"','"entity_ID"',
                 '"spatial_block"','"source_signature"',
                 '"island_latitude"','"island_longitude"'):
        if term in body:raise ValueError("Rights-restricted botanical evidence leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(status=result["status"],n_original_islands=4505,
                          n_original_source_reproductive_rows=len(result[
                              "original_regional_reproductive_shares_and_genus_isolation"]),
                          no_pollination_reproductive_causation=True)))


if __name__=="__main__":main()
