"""Original 4,505 GIFT floras, source-genus flower coupling by geological history.

Separate REAL source-verified geologically oceanic vs continental from
Weigelt LGM mainland connection; never classifies missing source geology by
distance/isolation. Only anonymous region/trait/source/qualifier results.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from run_fd_original_genus_flower_pair_covariance import (
    original_4505_genus_pair_frame, file_hash, SOURCE_HASHES,
)
from fd.geology_genus_flower_pair import (
    original_geology_join,original_pair_by_geology,
)


def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-root",type=Path,required=True)
    ap.add_argument("--geology-crosswalk",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    root=args.original_root
    ledger_path=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger_path,flora_path):
        if file_hash(path)!=SOURCE_HASHES[path.name]:
            raise ValueError("original real accepted species GIFT/Direct source archive drifted")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    geo=pd.read_csv(args.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "geological_origin","LGM_connectivity","area_km2",
        "distance_to_continent_km"])
    if len(geo)<8200 or geo.island_id.duplicated().any():
        raise ValueError("independent GSHHG geography and geology unverifiable")
    assigned=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    original=original_4505_genus_pair_frame(
        ledger,flora,geo,assigned,mainland)
    full=original_geology_join(original,geo,population=4505)
    report=original_pair_by_geology(full,population=4505,draws=999)
    n_geo=report["n_independent_GIFT_islands_with_oceanic_or_continental_geology"]
    n_glacial=report["n_original_islands_with_LGM_mainland_connectivity"]
    if not 400<=n_geo<=450:
        raise ValueError("original independent geology island source support changed")
    if n_glacial<1000:
        raise ValueError("Weigelt mainland past geography source no longer matches GIFT")
    report["original_flora_sha256"]=file_hash(flora_path)
    report["original_Direct_traits_sha256"]=file_hash(ledger_path)
    report["original_independent_geography_crosswalk_sha256"]=file_hash(
        args.geology_crosswalk)
    report["n_original_potential_mainland_GIFT_sources"]=748
    report["geology_context"]={
        "geological_origin_label_source":"Roeble_2024_Zizka_2022_independent_join",
        "LGM_connectivity_source":"Weigelt_2013_GMMC_distinct",
        "original_GIFT_island_floras_n":4505,
        "source_archipelago_and_study_independence_not_heldout":True,
        "matched_oceanic_continental_disjoint_pair_causality_not_claimed":True,
    }
    raw=json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)
    forbidden=(
        "island_id","accepted_species","work_species","entity_ID",
        "spatial_block","first_source","source_genus",
        "latitude","longitude","weigelt_record_id","geological_record_id",
        "original_species","archipelago","observer_id"
    )
    if any('"'+x+'"' in raw for x in forbidden):
        raise ValueError("protected GIFT accepted plant/source/location leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(raw+"\n",encoding="utf8")
    print(json.dumps(dict(
        status=report["status"],
        global_islands=report["original_global_botanical_island_floras"],
        geological_origin_classified_floras=n_geo,
        LGM_connected_or_unconnected_floras=n_glacial,
        n_panels=report["n_report_panels"],
        supported_spatial_block_panels=sum(
            p["spatial_block"]["status"]==
            "DESCRIPTIVE_GIFT_ORIGIN_SOURCE_RELATIVE_FLORAL_PAIRING"
            for p in report["panels"]
        ),
        supported_first_source_panels=sum(
            p["first_source"]["status"]==
            "DESCRIPTIVE_GIFT_ORIGIN_SOURCE_RELATIVE_FLORAL_PAIRING"
            for p in report["panels"]
        ),
    ),sort_keys=True))


if __name__=="__main__":
    main()
