"""Source-protected 4505 GIFT island REAL flower-pair source-only congener prediction.

Original accepted island names stay in protected CI. Source candidate accepted
species with ANY focal island flora listing are excluded from both source
modes before genus donor matching. Report only aggregate eligibility, source
group/spatial uncertainty and same-species joint phenotype forecast score.
"""
from __future__ import annotations
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd.real_joint_syndrome_access import make_joint_species
from fd.source_genus_same_island import genus
from fd.source_only_real_joint_congener_prediction import (
    SECTORS,SCOPES,REGIMES,DONORS,
    predict_joint_by_source_scopes,aggregate_source_congener_prediction,
)

SHA={
  "direct_only_trait_ledger.csv.gz":
    "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
  "chapter1_status_flora.csv.gz":
    "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
SOURCE_MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def digest(path:Path)->str:
    d=sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            d.update(block)
    return d.hexdigest()


def original_source_only_prediction_frame(
    direct_ledger:pd.DataFrame, flora:pd.DataFrame,
    geography:pd.DataFrame,assignments:pd.DataFrame,
    mainland:pd.DataFrame
)->pd.DataFrame:
    plants=status_groups(flora)
    if len(plants)!=4505:
        raise ValueError("original 4505 GIFT island flora frame changed")
    if geography.island_id.duplicated().any() or geography.island_id.nunique()<8200:
        raise ValueError("original independently sourced physical GSHHG islands missing")
    G=geography.set_index("island_id").reindex(sorted(plants))
    if G.analysis_regime.notna().sum()!=4504:
        raise ValueError("original 4504 georeferenced and one unknown island identity changed")
    source=assignments.loc[assignments.source_mode.isin(SOURCE_MODES)].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise").astype(int)
    if source.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("duplicate GIFT candidate mainland source assignment")
    for mode,(_,ranks) in SOURCE_MODES.items():
        if not source.loc[source.source_mode.eq(mode),"source_rank"].between(1,ranks).all():
            raise ValueError("invalid original source country rank")
    if {k:int(source.source_mode.eq(k).sum()) for k in SOURCE_MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("source assignment cardinalities changed")
    nearest=source.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("original accepted mainland botanical source names duplicated")
    region_plants=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(region_plants)!=748:
        raise ValueError("original 748 GIFT potential mainland source regions not present")
    @lru_cache(maxsize=16)
    def pooled(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(region_plants.get(r,set()) for r in ids))

    profiles=make_joint_species(strict_traits(direct_ledger))
    if not set(SECTORS)<=set(profiles):
        raise ValueError("actual same species symmetry and deep-tube characters unavailable")
    original_pair_sp={
        name:len(profiles[name]["known"]) for name in SECTORS
    }
    rows=[]
    with_climate=0
    for island in sorted(plants):
        item=plants[island]
        full=item["full"]
        sources={}
        for mode,(label,k) in SOURCE_MODES.items():
            keys=nearest.get((island,mode),())
            if keys and len(keys)!=k:
                raise ValueError("incomplete original selected nearest source ranks")
            sources[label]=set(pooled(keys))-full if keys else set()
        if nearest.get((island,"geo50_climate10")):
            with_climate+=1
        region=(str(G.at[island,"analysis_regime"])
                if pd.notna(G.at[island,"analysis_regime"])
                else "geography_unresolved")
        if region not in REGIMES:
            raise ValueError("unexpected original island GSHHG region")
        block=G.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        first=nearest.get((island,"geo_k5"),())
        donor_sig=first[0] if len(first)==5 else None
        plant_scopes={
            "all_recorded":full,
            "strict_native":item["native"],
            "strict_native_nonendemic":item["native_nonendemic"],
        }
        for sector in SECTORS:
            profile=profiles[sector]
            scores=predict_joint_by_source_scopes(
                full_island=full,scope_rosters=plant_scopes,
                geo5_source=sources["geo5"],
                climate10_source=sources["climate10"],
                profile=profile,
            )
            complete_by_source={
                mode:Counter(genus(sp) for sp in (source&profile["known"]))
                for mode,source in sources.items()
            }
            all_by_source={
                mode:Counter(genus(sp) for sp in source)
                for mode,source in sources.items()
            }
            for scope,names in plant_scopes.items():
                for gate in DONORS:
                    minimum=1 if gate=="any_congener" else 2
                    valid_genus={
                        g for g,n in complete_by_source["geo5"].items()
                        if n>=minimum
                    }&{
                        g for g,n in complete_by_source["climate10"].items()
                        if n>=minimum
                    }
                    observed_target=[
                        sp for sp in names&profile["known"]
                        if genus(sp) in valid_genus
                    ]
                    matched_genus={genus(sp) for sp in observed_target}
                    target_total=sum(genus(sp) in matched_genus for sp in names)
                    score=scores[scope][gate]
                    n_target=len(observed_target)
                    if n_target!=score["geo5"]["n_target_direct"] or (
                            n_target!=score["climate10"]["n_target_direct"]):
                        raise AssertionError("focal plant species roster changed by source mode")
                    row=dict(
                        island_id=island,regime=region,scope=scope,sector=sector,
                        donor_gate=gate,spatial_block=block,
                        source_signature=donor_sig,
                        n_target_direct=n_target,
                        n_target_genera=len(matched_genus),
                        island_direct_coverage_supported=(
                            n_target/target_total if target_total else np.nan
                        ),
                    )
                    for mode in ("geo5","climate10"):
                        pred=score[mode]
                        source_total=sum(
                            all_by_source[mode][g] for g in matched_genus
                        )
                        source_direct=sum(
                            complete_by_source[mode][g] for g in matched_genus
                        )
                        if source_direct!=pred["n_source_matched_direct"]:
                            raise AssertionError("congener real source Direct carrier count drifted")
                        for key,value in pred.items():
                            if key in ("n_source_matched_total",
                                       "target_species_per_source_genus_support"):
                                continue
                            row[mode+"__"+key]=value
                        row[mode+"__n_source_matched_total"]=source_total
                        row[mode+"_direct_coverage_supported"]=(
                            source_direct/source_total if source_total else np.nan
                        )
                    if (row["geo5__n_target_positive"]!=
                            row["climate10__n_target_positive"]):
                        raise AssertionError("same island target phenotypes differ by source mode")
                    rows.append(row)
    if with_climate!=4485:
        raise ValueError("original GIFT climate10 coverage changed from 4485")
    if len(rows)!=4505*len(SCOPES)*len(SECTORS)*len(DONORS):
        raise AssertionError("original all 4505 GIFT island rosters not retained")
    output=pd.DataFrame.from_records(rows)
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    ledger=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger,flora_path):
        if digest(path)!=SHA[path.name]:
            raise ValueError("original GIFT flora and flower Direct trait source HASH drifted")
    direct=pd.read_csv(ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"
    ])
    geography=pd.read_csv(args.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"
    ])
    source=pd.read_csv(
        args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"]
    )
    mainland=pd.read_csv(
        args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"]
    )
    original=original_source_only_prediction_frame(
        direct,flora,geography,source,mainland
    )
    result=aggregate_source_congener_prediction(
        original,population=4505,draws=999
    )
    result.update(
        n_original_GIFT_native_mainland_source_regions=748,
        n_original_selected_source_assigned_geo5=41325,
        n_original_selected_source_assigned_climate10=82190,
        original_direct_trait_ledger_sha256=digest(ledger),
        original_botanical_flora_sha256=digest(flora_path),
        original_accepted_botanical_source_ids_remain_private=True,
    )
    raw=json.dumps(result,sort_keys=True,ensure_ascii=False,indent=2,allow_nan=False)
    protected=(
        "island_id","accepted_species","work_species","entity_ID",
        "source_signature","spatial_block","source_region_id","source_genus",
        "source_family","latitude","longitude","taxon_id","GBIF_id","photo_id"
    )
    if any('"'+k+'"' in raw for k in protected):
        raise ValueError("original GIFT botanist identity/source rights violated")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "n_original_botanical_islands":result["n_original_GIFT_botanical_flora_islands"],
        "n_source_only_congener_joint_prediction_panels":len(result["panels"]),
        "n_panels_with_10_or_more_blocks":sum(
            q["all__spatial_block"]["status"]==
            "SOURCE_ONLY_NAME_HELDOUT_CONGENERIC_PREDICTION" for q in result["panels"]
        ),
        "no_local_evolution_or_effective_pollinator_service_claim":True,
    },sort_keys=True))


if __name__=="__main__":
    main()
