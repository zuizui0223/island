"""Original GIFT 4505 island REAL species-Direct 3- and 4-axis syndrome-compatible motifs.

Accepted species, original island IDs, GIFT candidate source regions, and
site records never leave protected CI. Anonymous aggregate report only.
A floral motif is NOT a validated guild assignment or effective pollen service.
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
from fd.syndrome_ensemble_source import (
    MOTIFS, SCOPES, REGIMES, make_motifs,
    make_island_motif_rows, aggregate_motifs,
)

SOURCE_HASHES={
 "direct_only_trait_ledger.csv.gz":
  "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
 "chapter1_status_flora.csv.gz":
  "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def file_hash(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def original_4505_motif_frame(
    ledger:pd.DataFrame,flora:pd.DataFrame,geo:pd.DataFrame,
    assigned:pd.DataFrame,mainland:pd.DataFrame,
)->pd.DataFrame:
    plants=status_groups(flora)
    if len(plants)!=4505:
        raise ValueError("original GIFT 4505 island botanical population missing")
    if geo.island_id.duplicated().any() or geo.island_id.nunique()<8200:
        raise ValueError("original authenticated GSHHG physical island crosswalk unavailable")
    G=geo.set_index("island_id").reindex(sorted(plants))
    if int(G.analysis_regime.notna().sum())!=4504:
        raise ValueError("original 4504 mapped botanical islands drifted")

    source=assigned.loc[assigned.source_mode.isin(MODES)].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise").astype(int)
    if source.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("GIFT candidate source region rank duplicated")
    if {k:int(source.source_mode.eq(k).sum()) for k in MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("GIFT original source mode assignments drifted")
    for mode,(_,limit) in MODES.items():
        if not source.loc[source.source_mode.eq(mode),"source_rank"].between(1,limit).all():
            raise ValueError("original potential mainland source rank invalid")

    nearest=source.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("duplicated accepted mainland species name")
    per_region=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(per_region)!=748:
        raise ValueError("original 748 native GIFT mainland source region archive changed")

    @lru_cache(maxsize=40)
    def pooled(keys:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(per_region.get(r,set()) for r in keys))

    motifs=make_motifs(strict_traits(ledger))
    if len(motifs)!=len(MOTIFS):
        raise ValueError("full genuine 3/4-axis motif ledger unavailable")

    rows=[]
    n_climate=0
    for island in sorted(plants):
        item=plants[island]
        full=item["full"]
        source_rosters={}
        for original_mode,(label,nregion) in MODES.items():
            region_ids=nearest.get((island,original_mode),())
            if region_ids and len(region_ids)!=nregion:
                raise ValueError("original GIFT candidate mainland pool truncated")
            source_rosters[label]=(set(pooled(region_ids))-full
                                    if region_ids else set())
        if nearest.get((island,"geo50_climate10")):
            n_climate+=1
        reg=G.at[island,"analysis_regime"]
        reg=str(reg) if pd.notna(reg) else "geography_unresolved"
        if reg not in REGIMES:
            raise ValueError("source-original regional stratum unrecognized")
        block=G.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        ids=nearest.get((island,"geo_k5"),())
        first=ids[0] if len(ids)==5 else None
        scopes=dict(all_recorded=full,strict_native=item["native"],
                    strict_native_nonendemic=item["native_nonendemic"])
        rows.extend(make_island_motif_rows(
            island_id=island,region=reg,spatial_block=block,first_source=first,
            scope_rosters=scopes,sources=source_rosters,motifs=motifs))
    if n_climate!=4485:
        raise ValueError("original GIFT climate10 source availability changed")
    if len(rows)!=4505*len(SCOPES)*len(MOTIFS):
        raise AssertionError("4505 island original motif-denominator contracted")
    return pd.DataFrame.from_records(rows)


def main()->None:
    cli=argparse.ArgumentParser()
    cli.add_argument("--original-root",type=Path,required=True)
    cli.add_argument("--geology-crosswalk",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    a=cli.parse_args()
    root=a.original_root
    ledger_path=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger_path,flora_path):
        if file_hash(path)!=SOURCE_HASHES[path.name]:
            raise ValueError("original botanical accepted plant/Direct trait SHA mismatch")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    geo=pd.read_csv(a.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assigned=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    rows=original_4505_motif_frame(ledger,flora,geo,assigned,mainland)
    out=aggregate_motifs(rows,population=4505,draws=999)
    out["original_flora_sha256"]=file_hash(flora_path)
    out["original_Direct_phenotype_sha256"]=file_hash(ledger_path)
    out["n_original_mainland_native_GIFT_regions"]=748
    raw=json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)
    forbidden=("island_id","accepted_species","work_species","entity_ID",
               "spatial_block","first_source","source_genus",
               "latitude","longitude","observer_id","GBIF_id","photo_id")
    if any('"'+k+'"' in raw for k in forbidden):
        raise ValueError("protected original island/botanical/source identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=out["status"],
        original_islands=out["n_original_GIFT_botanical_islands"],
        n_motifs=out["n_motifs"],
        n_panels=out["n_reporting_panels"],
        source_supported=sum(p["all__spatial_block"]["status"]==
            "EXPLORATORY_REAL_SAME_PLANT_MOTIF_SOURCE_COMPOSITION" for p in out["panels"])
    ),sort_keys=True))


if __name__=="__main__":
    main()
