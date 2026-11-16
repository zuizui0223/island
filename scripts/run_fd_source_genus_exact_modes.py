"""Exact signed GIFT mainland candidate genus-accounting on fixed 4,505 islands.

Reuses the independently synthetic-tested 4-term per-genus decomposition
but does NOT infer that genus composition caused island establishment or
that shared-genus differences are independent evolutionary adaptations.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    status_groups, source_sectors,
)
from fd_direct_source_genus_accounting import exact_genus_bin_decomposition
from fd.source_genus_exact_partition_modes import (
    REGIMES, SCOPES, MORPHOLOGY, TERMS,
    summarize_exact_genus_partition,
)

SHA = {
    "direct_only_trait_ledger.csv.gz":
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def digest(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):
            h.update(z)
    return h.hexdigest()


def build_original_private_decomposition(
    states:dict,flora:pd.DataFrame,geo:pd.DataFrame,
    assignment:pd.DataFrame,mainland:pd.DataFrame,
)->pd.DataFrame:
    island=status_groups(flora)
    if len(island)!=4505:raise ValueError("original GIFT 4505 islands missing")
    if geo.island_id.duplicated().any() or geo.island_id.nunique()<8200:
        raise ValueError("original GSHHG geography invalid")
    georef=geo.set_index("island_id").reindex(sorted(island))
    if georef.analysis_regime.notna().sum()!=4504:
        raise ValueError("one original unresolved island source geography changed")
    A=assignment.loc[assignment.source_mode.isin(MODES)].copy()
    A.source_rank=pd.to_numeric(A.source_rank,errors="raise").astype(int)
    if A.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("duplicate GIFT source mode rank")
    if {mode:int(A.source_mode.eq(mode).sum()) for mode in MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("source assignment mode population drift")
    for mode,(_,n) in MODES.items():
        if not A.loc[A.source_mode.eq(mode),"source_rank"].between(1,n).all():
            raise ValueError("invalid original mainland source rank")
    source=A.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("duplicate original mainland source plant record")
    members=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(members)!=748:
        raise ValueError("748 mainland native GIFT source regions changed")
    @lru_cache(maxsize=15)
    def pooled(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(members.get(x,set()) for x in ids))

    sectors=source_sectors(states)
    if not set(MORPHOLOGY)<=set(sectors):
        raise ValueError("Direct sector original symmetry/form data missing")
    rows=[]
    n_climate=0
    for id in sorted(island):
        F=island[id]
        full=F["full"]
        all_sources={}
        for mode,(_,k) in MODES.items():
            ranks=source.get((id,mode),())
            if ranks and len(ranks)!=k:
                raise ValueError("partial mainland native source ranking")
            all_sources[mode]=pooled(ranks) if ranks else frozenset()
        if all_sources["geo50_climate10"]:
            n_climate+=1
        region=str(georef.at[id,"analysis_regime"]) if pd.notna(
            georef.at[id,"analysis_regime"]) else "geography_unresolved"
        if region not in REGIMES:raise ValueError("GSHHG original region unknown")
        block=georef.at[id,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        r0=source.get((id,"geo_k5"),())
        sig=r0[0] if len(r0)==5 else None
        scopes={
            "all_recorded":full,
            "strict_native":F["native"],
            "strict_native_nonendemic":F["native_nonendemic"],
            "native_nonendemic_source_listed_both":
                F["native_nonendemic"] &
                all_sources["geo_k5"] & all_sources["geo50_climate10"],
        }
        for sector in MORPHOLOGY:
            scored=sectors[sector]
            known=set(scored)
            candidates={
                name:(set(all_sources[name])-full)&known
                for name in MODES
            }
            for scope,plants in scopes.items():
                target=plants&known
                n=len(target)
                row=dict(
                    island_id=id,regime=region,scope=scope,sector=sector,
                    spatial_block=block,source_signature=sig,
                    n_scope_total=len(plants),n_scope_direct=n,
                )
                for mode,(label,_) in MODES.items():
                    donors=candidates[mode]
                    row[label+"__n_source_total"]=len(set(all_sources[mode])-full)
                    row[label+"__n_source_direct"]=len(donors)
                    for field in ("observed_difference",*TERMS):
                        row[label+"__"+field]=np.nan
                    row[label+"__n_shared_genera"]=0
                    row[label+"__n_native_carriers_from_shared_genera"]=0
                    if n<5 or len(donors)<10:
                        continue
                    part=exact_genus_bin_decomposition(target,donors,scored)
                    row[label+"__observed_difference"]=part[
                        "original_direct_binary_sector_difference"
                    ]
                    for t in TERMS:
                        row[label+"__"+t]=part[t]
                    row[label+"__n_shared_genera"]=part["n_both_status_shared_genera"]
                    row[label+"__n_native_carriers_from_shared_genera"]=part[
                        "n_native_source_carriers_from_shared_genera"]
                rows.append(row)
    if n_climate!=4485:
        raise ValueError("original 4485 climate-source-covered islands drifted")
    if len(rows)!=4505*len(SCOPES)*len(MORPHOLOGY):
        raise AssertionError("original all-island genus analysis inadvertently restricted")
    return pd.DataFrame(rows)


def run(root:Path,geo_path:Path)->dict:
    ledgerpath=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    florapath=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for p in (ledgerpath,florapath):
        if digest(p)!=SHA[p.name]:
            raise ValueError("original accepted botanical source SHA mismatch")
    ledger=pd.read_csv(ledgerpath,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(florapath,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"
    ])
    geography=pd.read_csv(geo_path,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"
    ])
    assigned=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"]
    )
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"]
    )
    states=strict_traits(ledger)
    private=build_original_private_decomposition(
        states,flora,geography,assigned,mainland
    )
    res=summarize_exact_genus_partition(private,population=4505,bootstrap=999)
    if res["n_original_botanical_floras"]!=4505 or res["n_region_scope_sector_panels"]!=100:
        raise AssertionError("protected original 4505 genus accounting cohort invalid")
    return dict(
        **res,real_original_GIFT_source_entities=748,
        original_source_assignment_counts={"geo5":41325,"climate10":82190},
        original_accepted_flora_sha256=digest(florapath),
        original_Direct_trait_ledger_sha256=digest(ledgerpath),
        original_taxa_and_GIFT_entities_kept_private=True,
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--geology-crosswalk",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    obj=run(a.original_root,a.geology_crosswalk)
    txt=json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)
    excluded=("island_id","accepted_species","entity_ID","work_species",
              "source_signature","spatial_block","source_genus",
              "source_family","taxon_id","GBIF_id","latitude","longitude")
    if any('"'+k+'"' in txt for k in excluded):
        raise ValueError("source rights violation")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(txt+"\n",encoding="utf-8")
    supported=sum(
        r["all__spatial_block"]["status"] ==
        "EXACT_SIGNED_GENUS_ACCOUNTING_DESCRIPTIVE_NOT_CAUSAL"
        for r in obj["panels"]
    )
    print(json.dumps(dict(
        status=obj["status"],
        n_original_botanical_floras=4505,
        regional_scope_sector_panels=obj["n_region_scope_sector_panels"],
        supported_original_exact_source_genus_partition_panels=supported,
        not_causal_or_phylogenetically_independent=True,
    ),sort_keys=True))


if __name__=="__main__":
    main()
