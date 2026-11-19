"""Source-protected real 4,505-island GIFT genus-conditional floral coupling.

The analysis asks whether previously observed excess co-occurrence of actual
same-species floral symmetry×depth can be reconstructed by fixed botanic genus
representation and observed marginal state frequencies WITHIN each genus.
Within-genus A/B label exchange reference is exact, not Monte Carlo.
No original accepted names, island identifiers, region IDs or coordinates are
exported from source CI.
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
from fd.real_joint_within_genus_coupling import (
    PRIORITY,SCOPES,REGIMES,
    one_roster_analytic_genus_pair,
    summarize_original_genus_joint_coupling,
)

HASHES={
 "direct_only_trait_ledger.csv.gz":
 "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
 "chapter1_status_flora.csv.gz":
 "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def sha(path:Path)->str:
    d=sha256()
    with path.open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):
            d.update(z)
    return d.hexdigest()


def original_global_genus_coupling_rows(
    ledger:pd.DataFrame,flora:pd.DataFrame,
    geodata:pd.DataFrame,assign:pd.DataFrame,
    mainland:pd.DataFrame,
)->pd.DataFrame:
    original=status_groups(flora)
    if len(original)!=4505:
        raise ValueError("original fixed 4505 GIFT botanical floras altered")
    if geodata.island_id.duplicated().any() or geodata.island_id.nunique()<8200:
        raise ValueError("original physically resolved GSHHG island crosswalk missing")
    geo=geodata.set_index("island_id").reindex(sorted(original))
    if int(geo.analysis_regime.notna().sum())!=4504:
        raise ValueError("one original unresolved source geography altered")
    source=assign.loc[assign.source_mode.isin(MODES)].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise").astype(int)
    if source.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("duplicate original GIFT mainland rank")
    for mode,(_,rank) in MODES.items():
        if not source.loc[source.source_mode.eq(mode),"source_rank"].between(1,rank).all():
            raise ValueError("invalid GIFT nearest source region rank")
    expected={"geo_k5":41325,"geo50_climate10":82190}
    if {mode:int(source.source_mode.eq(mode).sum()) for mode in MODES}!=expected:
        raise ValueError("original GIFT full source rank assignment drifted")
    rankmap=source.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("original accepted source entity×species duplicated")
    donors=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(donors)!=748:
        raise ValueError("original 748 native GIFT candidate source regions missing")

    @lru_cache(maxsize=20)
    def region_union(regions:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(donors.get(r,set()) for r in regions))

    profiles=make_joint_species(strict_traits(ledger))
    if not set(PRIORITY)<=set(profiles):
        raise ValueError("true original same-plant symmetry/depth state missing")
    rows=[]
    n_climate=0
    for i in sorted(original):
        q=original[i]
        whole=q["full"]
        candidate={}
        for mode,(_,nranks) in MODES.items():
            ranks=rankmap.get((i,mode),())
            if ranks and len(ranks)!=nranks:
                raise ValueError("incomplete original regional source ranking")
            candidate[mode]=(set(region_union(ranks))-whole) if ranks else set()
        if rankmap.get((i,"geo50_climate10")):n_climate+=1
        regime=(str(geo.at[i,"analysis_regime"])
                if pd.notna(geo.at[i,"analysis_regime"]) else "geography_unresolved")
        if regime not in REGIMES:
            raise ValueError("unexpected original island analysis regime")
        block=geo.at[i,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        nearest=rankmap.get((i,"geo_k5"),())
        source_signature=nearest[0] if len(nearest)==5 else None
        botanical={
            "all_recorded":whole,
            "strict_native":q["native"],
            "strict_native_nonendemic":q["native_nonendemic"],
        }
        if set(botanical)!=set(SCOPES):
            raise AssertionError("original native and nonendemic controls lost")
        for sector in PRIORITY:
            profile=profiles[sector]
            mainland_stats={
                label:one_roster_analytic_genus_pair(candidate[mode],profile)
                for mode,(label,_) in MODES.items()
            }
            for scope,plants in botanical.items():
                local=one_roster_analytic_genus_pair(plants,profile)
                record=dict(
                    island_id=i,regime=regime,scope=scope,sector=sector,
                    spatial_block=block,source_signature=source_signature,
                    island_total=len(plants),
                    geo5_total=len(candidate["geo_k5"]),
                    climate10_total=len(candidate["geo50_climate10"]),
                )
                for prefix,stats in [
                    ("island",local),
                    ("geo5",mainland_stats["geo5"]),
                    ("climate10",mainland_stats["climate10"]),
                ]:
                    for k,val in stats.items():
                        record[prefix+"_"+k]=val
                rows.append(record)
    if n_climate!=4485:
        raise ValueError("original climate source source covered islands not 4485")
    if len(rows)!=4505*len(SCOPES)*len(PRIORITY):
        raise AssertionError("original complete island roster narrowed to trait subset")
    return pd.DataFrame.from_records(rows)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--geology-crosswalk",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    ledger_path=(a.original_root/
        "fixed/source/current/evidence/direct_only_trait_ledger.csv.gz")
    flora_path=(a.original_root/
        "fixed/canonical/input/chapter1_status_flora.csv.gz")
    for file in (ledger_path,flora_path):
        if sha(file)!=HASHES[file.name]:
            raise ValueError("original published botanical Direct source SHA mismatch")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality",
    ])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status",
    ])
    geodata=pd.read_csv(a.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
    ])
    assigned=pd.read_csv(
        a.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"],
    )
    mainland=pd.read_csv(
        a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"],
    )
    original_rows=original_global_genus_coupling_rows(
        ledger,flora,geodata,assigned,mainland,
    )
    obj=summarize_original_genus_joint_coupling(
        original_rows,population=4505,boots=999,
    )
    obj.update(
        original_source_748_GIFT_native_mainland_regions=748,
        original_geo5_source_assignment_rows=41325,
        original_climate10_source_assignment_rows=82190,
        original_Direct_trait_ledger_sha256=sha(ledger_path),
        original_GIFT_botanical_flora_sha256=sha(flora_path),
        original_source_island_taxon_identifiers_not_exported=True,
    )
    raw=json.dumps(obj,sort_keys=True,ensure_ascii=False,indent=2,allow_nan=False)
    prohibited=(
        "island_id","accepted_species","work_species","entity_ID",
        "spatial_block","source_signature","source_region_id",
        "latitude","longitude","source_genus","source_family",
        "photo_id","GBIF_id","taxon_id",
    )
    if any('"'+k+'"' in raw for k in prohibited):
        raise ValueError("original island/plant/GIFT source identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=obj["status"],original_islands=obj["n_original_GIFT_botanical_islands"],
        n_report_panels=obj["n_region_scope_structural_pair_panels"],
        n_block_supported_panels=sum(
            q["all_direct_qualified__spatial_block"]["status"]==
            "DESCRIPTIVE_GENUS_CONDITIONAL_COUPLING_NOT_CAUSAL"
            for q in obj["panels"]
        ),
        no_effective_pollinator_guild_inferred=True,
    ),sort_keys=True))


if __name__=="__main__":
    main()
