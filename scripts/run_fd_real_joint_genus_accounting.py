"""Original protected botanical genus decomposition of real same-species flower pairs.

A reported pair always belongs to an accepted species with BOTH trait axes
genuinely Direct-labelled. Uses GIFT source-only distinct accepted mainland
candidate species, not historical ancestors or effective pollinator guilds.
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
from fd_direct_three_stage_floristic_comparison import status_groups
from fd.real_joint_syndrome_access import make_joint_species
from fd.real_joint_genus_accounting import (
    TARGETS, SCOPES, REGIMES, PARTS, summarize_real_joint_genus,
)
from run_fd_source_genus_exact_modes import (
    source_genus_counts,exact_genus_from_source_counts,
)

SHA={
    "direct_only_trait_ledger.csv.gz":
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def digest(p:Path)->str:
    h=sha256()
    with p.open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):
            h.update(z)
    return h.hexdigest()


def original_joint_genus_rows(
    ledger:pd.DataFrame,fl:pd.DataFrame,geo:pd.DataFrame,
    assigned:pd.DataFrame,mainland:pd.DataFrame,
)->pd.DataFrame:
    islands=status_groups(fl)
    if len(islands)!=4505:
        raise ValueError("source original 4505 botanical island study frame reduced")
    if geo.island_id.duplicated().any() or geo.island_id.nunique()<8200:
        raise ValueError("original independently sourced physical island geography invalid")
    g=geo.set_index("island_id").reindex(sorted(islands))
    if g.analysis_regime.notna().sum()!=4504:
        raise ValueError("original one geographically unresolved GIFT island lost")
    A=assigned.loc[assigned.source_mode.isin(MODES)].copy()
    A.source_rank=pd.to_numeric(A.source_rank,errors="raise").astype(int)
    if A.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("original GIFT source ranks duplicated")
    if {k:int(A.source_mode.eq(k).sum()) for k in MODES}!={"geo_k5":41325,"geo50_climate10":82190}:
        raise ValueError("independent geographic/climate source assignment changed")
    byisland=A.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True).entity_ID.agg(tuple).to_dict()
    byentity=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(byentity)!=748 or mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("original mainland accepted GIFT source identities changed")

    @lru_cache(maxsize=15)
    def pooled(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(byentity.get(z,set()) for z in ids))

    pair=make_joint_species(strict_traits(ledger))
    states={}
    for name in TARGETS:
        profile=pair[name]
        states[name]={
            sp:int(sp in profile["joint"]) for sp in profile["known"]
        }
    rows=[]
    n_climate=0
    for island in sorted(islands):
        item=islands[island]
        full=item["full"]
        source={}
        for mode,(_,n) in MODES.items():
            ranking=byisland.get((island,mode),())
            if ranking and len(ranking)!=n:
                raise ValueError("original ranked GIFT source missing entities")
            source[mode]=set(pooled(ranking))-full if ranking else set()
        if byisland.get((island,"geo50_climate10")):
            n_climate+=1
        region=(str(g.at[island,"analysis_regime"])
                if pd.notna(g.at[island,"analysis_regime"])
                else "geography_unresolved")
        if region not in REGIMES:
            raise ValueError("original source regime unknown")
        block=g.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        src5=byisland.get((island,"geo_k5"),())
        sig=src5[0] if len(src5)==5 else None
        status={
            "all_recorded":full,
            "strict_native":item["native"],
            "strict_native_nonendemic":item["native_nonendemic"],
        }
        for sector in TARGETS:
            direct=states[sector]
            available=set(direct)
            donors={mode:source[mode]&available for mode in MODES}
            dgenus={
                mode:source_genus_counts(donors[mode],direct)
                for mode in MODES
            }
            for scope,names in status.items():
                target=names&available
                row=dict(
                    island_id=island,region=region,scope=scope,sector=sector,
                    spatial_block=block,first_source=sig,
                    n_island_total=len(names),n_island_direct=len(target),
                )
                for mode,(label,_) in MODES.items():
                    row[label+"__n_mainland_total"]=len(source[mode])
                    row[label+"__n_mainland_direct"]=len(donors[mode])
                    row[label+"__n_shared_genera"]=0
                    for term in ("observed_difference",*PARTS):
                        row[label+"__"+term]=np.nan
                    if len(target)<5 or len(donors[mode])<10:
                        continue
                    z=exact_genus_from_source_counts(target,dgenus[mode],direct)
                    row[label+"__observed_difference"]=z[
                        "original_direct_binary_sector_difference"]
                    for term in PARTS:
                        row[label+"__"+term]=z[term]
                    row[label+"__n_shared_genera"]=z["n_both_status_shared_genera"]
                rows.append(row)
    if n_climate!=4485:
        raise ValueError("only original 4485/4505 climate mainland source regions qualify")
    if len(rows)!=4505*len(SCOPES)*len(TARGETS):
        raise AssertionError("original real pair phenotype study frame collapsed")
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-root",type=Path,required=True)
    ap.add_argument("--geology-crosswalk",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    ledger=a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for f in (ledger,flora):
        if digest(f)!=SHA[f.name]:
            raise ValueError("original accepted plant botanical SHA source drift")
    traits=pd.read_csv(ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    fl=pd.read_csv(flora,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    g=pd.read_csv(a.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    sources=pd.read_csv(
        a.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    frame=original_joint_genus_rows(traits,fl,g,sources,mainland)
    obj=summarize_real_joint_genus(frame,population=4505,draws=999)
    obj["original_direct_ledger_SHA256"]=digest(ledger)
    obj["original_accepted_flora_SHA256"]=digest(flora)
    obj["mainland_GIFT_native_source_entities"]=748
    text=json.dumps(obj,ensure_ascii=False,sort_keys=True,
                    indent=2,allow_nan=False)
    prohibited=("island_id","accepted_species","work_species","entity_ID",
                "spatial_block","source_signature","source_genus","latitude",
                "longitude","taxon_id","GBIF_id","source_family")
    if any('"'+k+'"' in text for k in prohibited):
        raise ValueError("private original accepted island/plant/source evidence leak")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=obj["status"],
        n_original_global_GIFT_botanical_floras=4505,
        n_actual_same_species_structural_joints=len(TARGETS),
        n_original_region_scope_sector_panels=len(obj["panels"]),
        n_geographically_supported=sum(
            p["all__spatial_block"]["status"]==
            "EXACT_REAL_PAIRED_FLORA_GENUS_ACCOUNTING_NOT_CAUSAL"
            for p in obj["panels"]
        ),
    ),sort_keys=True))


if __name__=="__main__":
    main()
