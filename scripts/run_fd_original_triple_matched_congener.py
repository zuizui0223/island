"""Original 4,505 GIFT flora, SAME target and SAME mainland donor plant three-trait check.

Private accepted plant/island/source names remain in original island Actions.
Only anonymized regional conditional estimates leave the protected runner.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd.real_joint_syndrome_access import make_joint_species
from fd.triple_matched_congener import (
    SCOPES, REGIMES, triple_matched_island, aggregate_triple,
)
from fd.source_loo_flower_placebo import summarize_source_loo_placebo

EXPECTED_SHA = {
    "direct_only_trait_ledger.csv.gz":
       "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
       "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
SOURCE_MODES = {"geo_k5":("geo5",5),
                "geo50_climate10":("climate10",10)}


def digest(path: Path) -> str:
    q=sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            q.update(b)
    return q.hexdigest()


def real_triple_frame(
    direct:pd.DataFrame, flora:pd.DataFrame, geography:pd.DataFrame,
    assignments:pd.DataFrame, mainland:pd.DataFrame,
) -> pd.DataFrame:
    plants=status_groups(flora)
    if len(plants)!=4505:
        raise ValueError("original 4505 botanical islands changed")
    if geography.island_id.duplicated().any() or geography.island_id.nunique()<8200:
        raise ValueError("original GSHHG physical island geography unavailable")
    g=geography.set_index("island_id").reindex(sorted(plants))
    if g.analysis_regime.notna().sum()!=4504:
        raise ValueError("original geographic island match changed")
    source=assignments.loc[
        assignments.source_mode.isin(SOURCE_MODES)
    ].copy()
    source.source_rank=pd.to_numeric(
        source.source_rank,errors="raise"
    ).astype(int)
    if source.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("GIFT original mainland region ranks duplicated")
    for mode,(_,ranks) in SOURCE_MODES.items():
        if not source.loc[source.source_mode.eq(mode),"source_rank"].between(1,ranks).all():
            raise ValueError("original mainland source rank out of bounds")
    if {k:int(source.source_mode.eq(k).sum()) for k in SOURCE_MODES} != {
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("original candidate mainland source assignments changed")
    nearest=source.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("mainland source accepted species duplicated")
    by_source=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(by_source)!=748:
        raise ValueError("original 748 mainland GIFT botanical regions changed")

    @lru_cache(maxsize=32)
    def pooled(region_ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(
            *(by_source.get(k,set()) for k in region_ids)
        )

    profiles=make_joint_species(strict_traits(direct))
    rows=[]
    for island in sorted(plants):
        item=plants[island]
        full=item["full"]
        candidate={}
        for original_mode,(label,nsource) in SOURCE_MODES.items():
            region_ids=nearest.get((island,original_mode),())
            if region_ids and len(region_ids)!=nsource:
                raise ValueError("original candidate source rank list truncated")
            candidate[label]=set(pooled(region_ids))-full if region_ids else set()
        original_regime=g.at[island,"analysis_regime"]
        region=(str(original_regime) if pd.notna(original_regime)
                else "geography_unresolved")
        if region not in REGIMES:
            raise ValueError("invalid original GSHHG geographical regime")
        block=g.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        first=nearest.get((island,"geo_k5"),())
        first=first[0] if len(first)==5 else None
        statuses={"all_recorded":full,
                  "strict_native":item["native"],
                  "strict_native_nonendemic":item["native_nonendemic"]}
        for res in triple_matched_island(
            full_island=full,scope_rosters=statuses,
            geo5_source=candidate["geo5"],
            climate10_source=candidate["climate10"],
            profiles=profiles,
        ):
            row=dict(island_id=island,regime=region,
                     spatial_block=block,source_signature=first)
            row.update(res)
            rows.append(row)
    if len(rows)!=4505*len(SCOPES)*2*2:
        raise AssertionError("source-protected original population lost")
    return pd.DataFrame(rows)


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--geology-crosswalk",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    root=a.original_root
    direct_path=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (direct_path,flora_path):
        if digest(path)!=EXPECTED_SHA[path.name]:
            raise ValueError("original floral Direct or GIFT flora archive SHA drift")
    direct=pd.read_csv(direct_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality",
    ])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status",
    ])
    geography=pd.read_csv(a.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
    ])
    assignments=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"],
    )
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"],
    )
    frame=real_triple_frame(direct,flora,geography,assignments,mainland)
    result=aggregate_triple(frame,population=4505,draws=999)
    result["source_loo_negative_control"]=(
        summarize_source_loo_placebo(frame,population=4505,draws=999)
    )
    result.update(
        original_ledger_sha256=digest(direct_path),
        original_flora_sha256=digest(flora_path),
        n_original_mainland_regions=748,
    )
    payload=json.dumps(result,ensure_ascii=False,sort_keys=True,
                       indent=2,allow_nan=False)
    # The output must be aggregate-only. No private botanical identity or
    # precise locality/source record can leave protected original Actions.
    forbidden=("island_id","accepted_species","work_species","entity_ID",
               "source_signature","spatial_block","source_genus","latitude",
               "longitude","photo_id","observer_id")
    if any('"'+key+'"' in payload for key in forbidden):
        raise ValueError("private original GIFT botanical identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(payload+"\n",encoding="utf8")
    print(json.dumps(dict(status=result["status"],
        n_original_botanical_islands=result["n_original_botanical_islands"],
        n_panels=result["n_panels"],
        n_source_loo_panels=result["source_loo_negative_control"]["n_reporting_panels"],
        n_source_loo_identifiable=sum(
            p["status"]=="DESCRIPTIVE_REAL_ISLAND_VS_SOURCE_LOO"
            for p in result["source_loo_negative_control"]["panels"]
        ),
        identifiable_panels=sum(
            p["status"]=="DESCRIPTIVE_TRIPLE_DIRECT_MATCHED_SOURCE"
            for p in result["panels"]
        )),sort_keys=True))


if __name__=="__main__":
    main()
