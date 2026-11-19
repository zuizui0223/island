"""Original GIFT 4505 floral syndrome-access *morphology* joint-pair audit.

These independently coded morphological conjunctions are not observed
effective guild labels. Every qualifying profile is an ACTUAL accepted
species with both trait axes directly documented, never marginal products
of traits recorded on different plants.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd.real_joint_syndrome_access import (
    JOINT_SPECS, SCOPES, REGIMES, make_joint_species, summarize_one_roster,
    same_island_pair_source_frame,
)

SOURCE_HASHES={
    "direct_only_trait_ledger.csv.gz":
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def file_sha(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def run(
    ledger:pd.DataFrame, flora:pd.DataFrame, geography:pd.DataFrame,
    assignments:pd.DataFrame, mainland:pd.DataFrame,
)->dict:
    full=status_groups(flora)
    if len(full)!=4505:
        raise ValueError("Original 4505 GIFT floras not available")
    if geography.island_id.duplicated().any() or geography.island_id.nunique()<8200:
        raise ValueError("original independently geocoded GSHHG islands missing")
    geo=geography.set_index("island_id").reindex(sorted(full))
    if int(geo.analysis_regime.notna().sum())!=4504:
        raise ValueError("original geography reconciliation drift")
    ranks=assignments.loc[assignments.source_mode.isin(MODES)].copy()
    if ranks.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("GIFT source mode assignment duplicated")
    ranks.source_rank=pd.to_numeric(ranks.source_rank,errors="raise").astype(int)
    for mode,(_,k) in MODES.items():
        if not ranks.loc[ranks.source_mode.eq(mode),"source_rank"].between(1,k).all():
            raise ValueError("GIFT mainland source ranking changed")
    if {mode:int(ranks.source_mode.eq(mode).sum()) for mode in MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("source geography/climate assignment original snapshot changed")
    choice=ranks.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("mainland GIFT entity species names duplicated")
    region_species=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(region_species)!=748:
        raise ValueError("original 748 mainland native GIFT regions changed")

    @lru_cache(maxsize=25)
    def pooled_source(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(region_species.get(x,set()) for x in ids))

    states=strict_traits(ledger)
    profiles=make_joint_species(states)
    if len(profiles)!=len(JOINT_SPECS):
        raise AssertionError("real pair trait sector ontology drift")
    original_color_counts=Counter(states["flower_primary_color"].values())
    original_species_pair_counts={
        k:len(v["known"]) for k,v in profiles.items()
    }
    raw=[]
    climate_count=0
    for island in sorted(full):
        flora_i=full[island]
        all_island_names=flora_i["full"]
        src={}
        for mode,(_,n) in MODES.items():
            ids=choice.get((island,mode),())
            if ids and len(ids)!=n:
                raise ValueError("incomplete mainland GIFT regional source set")
            if ids:
                src[mode]=set(pooled_source(ids))-all_island_names
            else:
                src[mode]=set()
        if choice.get((island,"geo50_climate10")):
            climate_count+=1
        regime=(str(geo.at[island,"analysis_regime"])
                if pd.notna(geo.at[island,"analysis_regime"])
                else "geography_unresolved")
        if regime not in REGIMES:
            raise ValueError("unexpected GSHHG region")
        block=geo.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        top=choice.get((island,"geo_k5"),())
        first_source=top[0] if len(top)==5 else None
        scopes={
            "all_recorded":flora_i["full"],
            "strict_native":flora_i["native"],
            "strict_native_nonendemic":flora_i["native_nonendemic"],
        }
        for sector,profile in profiles.items():
            source_stats={
                mode:summarize_one_roster(src[mode],profile)
                for mode in MODES
            }
            for scope,names in scopes.items():
                island_stats=summarize_one_roster(names,profile)
                z=dict(
                    island_id=island,regime=regime,scope=scope,sector=sector,
                    spatial_block=block,first_source=first_source,
                    island_total=len(names),
                    island_direct=island_stats["n"],
                    island_pair_prevalence=island_stats["p_joint"],
                    island_pair_coupling=island_stats["coupling"],
                )
                for mode,(label,_) in MODES.items():
                    q=source_stats[mode]
                    z[label+"_total"]=len(src[mode])
                    z[label+"_direct"]=q["n"]
                    z[label+"_pair_prevalence"]=q["p_joint"]
                    z[label+"_pair_coupling"]=q["coupling"]
                raw.append(z)
    if climate_count!=4485:
        raise ValueError("original source-covered 4485 climate10 islands drift")
    if len(raw)!=4505*len(SCOPES)*len(profiles):
        raise AssertionError("real same-species pair original full island frame reduced")
    results=same_island_pair_source_frame(
        pd.DataFrame(raw),population=4505,bootstrap=999,
    )
    if results["n_source_mode_provenance_region_panels"]!=len(REGIMES)*len(SCOPES)*len(JOINT_SPECS):
        raise AssertionError("floral syndrome joint support panel coverage drift")
    results.update(
        n_original_748_mainland_source_entities=len(region_species),
        n_original_geographic5_source_assignment_rows=41325,
        n_original_climate10_source_assignment_rows=82190,
        n_original_climate10_covered_botanical_floras=climate_count,
        original_unique_Direct_accepted_species_with_both_trait_axes=original_species_pair_counts,
        original_single_accepted_species_Direct_colour_category_counts=dict(sorted(original_color_counts.items())),
    )
    return results


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--geology-crosswalk",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    ledger_path=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger_path,flora_path):
        if file_sha(path)!=SOURCE_HASHES[path.name]:
            raise ValueError("original GIFT botanical species Direct source SHA changed")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality",
    ])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status",
        "floristic_status",
    ])
    geography=pd.read_csv(args.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
    ])
    assigned=pd.read_csv(
        args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"],
    )
    mainland=pd.read_csv(
        args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"],
    )
    output=run(ledger,flora,geography,assigned,mainland)
    output["original_botanical_flora_sha256"]=file_sha(flora_path)
    output["original_Direct_species_trait_ledger_sha256"]=file_sha(ledger_path)
    raw=json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)
    forbidden=(
        "island_id","accepted_species","source_signature","spatial_block",
        "work_species","entity_ID","latitude","longitude","source_genus",
        "source_family","observer_id","GBIF_id","photo_id",
    )
    if any('"'+key+'"' in raw for key in forbidden):
        raise ValueError("original protected source IDs leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(raw+"\n",encoding="utf-8")
    qualified=sum(z["all__spatial_block"]["status"]==
        "DESCRIPTIVE_REAL_SAME_SPECIES_PAIR_ONLY" for z in output["panels"])
    print(json.dumps(dict(
        status=output["status"],global_GIFT_floras=4505,
        true_same_species_two_trait_combinations=len(JOINT_SPECS),
        n_provenance_region_pair_panels=len(output["panels"]),
        n_with_geographic_support=qualified,
        historical_bird_butterfly_effect_inferred=False,
    ),sort_keys=True))


if __name__=="__main__":
    main()
