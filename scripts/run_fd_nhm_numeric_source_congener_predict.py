"""Original GIFT 4505 independent NHM numeric corolla tube length prediction.

Train deep/nondeep calibration ONLY on GIFT native mainland accepted species
with BOTH NHM numerical cm and original Direct ordinal depth, after excluding
EVERY original island flora plant accepted name from the calibration. Use
different source-only congener accepted names in two region-mode sources to
predict the original island deep category. No original island label enters fit.

NHM source provider corpora are separately pinned; the combined.csv release
does not create a sixth independent botanical source.
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
from fd.nhm_numeric_source_congener_predict import (
    REGIMES,SCOPES,TIERS,
    make_genus_heldout_length_calibrator,source_predict_same_island,
    aggregate_external_source_numeric_prediction,
)
from run_fd_independent_nhm_numeric_tube_length import (
    parse_nhm_numeric_file,combine_external_provider_lengths,
)
from fd_nhm_original_provider_audit import PROVIDERS,BLOB_PINS
from fd_nhm_external_flora_preflight import nhm_git_blob_sha

PIN={
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


def fit_original_and_predict(
    ledger:pd.DataFrame,flora:pd.DataFrame,
    geography:pd.DataFrame,assigned:pd.DataFrame,
    mainland:pd.DataFrame,
    numeric_provider:dict[str,dict]
)->dict:
    original=status_groups(flora)
    if len(original)!=4505:
        raise ValueError("original 4505 accepted botanical island flora unavailable")
    if geography.island_id.duplicated().any() or geography.island_id.nunique()<8200:
        raise ValueError("GSHHG physical island geographic crosswalk incomplete")
    geo=geography.set_index("island_id").reindex(sorted(original))
    if int(geo.analysis_regime.notna().sum())!=4504:
        raise ValueError("original 4504 geocoded one unresolved island split changed")
    if mainland.duplicated(["entity_ID","work_species"]).any() or (
        mainland.entity_ID.nunique()!=748):
        raise ValueError("original GIFT native mainland 748 source entities drifted")
    gift=set(mainland.work_species)
    all_plant_names=set().union(*(a["full"] for a in original.values()))
    traits=strict_traits(ledger)
    original_depth=traits["tube_depth_class"]
    external,two,prov=combine_external_provider_lengths(numeric_provider)
    # The original strictly OUTSIDE-ALL-ISLANDS named-species
    # calibration failed with precisely zero numeric×ordinal overlap.
    # Use exact leave-one-BOTANICAL-GENUS-out calibration, with no target
    # genus categorical labels included in its numeric-to-depth mapping.
    cal,calibration_meta=make_genus_heldout_length_calibrator(
        external,original_depth,gift
    )
    source_modes=assigned.loc[assigned.source_mode.isin(MODES)].copy()
    source_modes.source_rank=pd.to_numeric(
        source_modes.source_rank,errors="raise").astype(int)
    if source_modes.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("duplicate GIFT candidate source island×mode×rank")
    if {m:int(source_modes.source_mode.eq(m).sum()) for m in MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("original GIFT island mainland candidate rank cardinality drifted")
    for mode,(_,limit) in MODES.items():
        if not source_modes.loc[
            source_modes.source_mode.eq(mode),"source_rank"
        ].between(1,limit).all():
            raise ValueError("original regional source ordering invalid")
    byisland=source_modes.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    members=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    @lru_cache(maxsize=20)
    def pooled(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(members.get(name,set()) for name in ids))

    records=[]
    nclimate=0
    for island in sorted(original):
        q=original[island]
        full=q["full"]
        sources={}
        for mode,(label,n_rank) in MODES.items():
            ids=byisland.get((island,mode),())
            if ids and len(ids)!=n_rank:
                raise ValueError("incomplete source region shortlist")
            sources[label]=set(pooled(ids))-full if ids else set()
            if sources[label]&full:
                raise AssertionError("focal plant name leaked into mainland numeric predictor")
        if byisland.get((island,"geo50_climate10")):
            nclimate+=1
        region=(str(geo.at[island,"analysis_regime"])
                if pd.notna(geo.at[island,"analysis_regime"])
                else "geography_unresolved")
        if region not in REGIMES:raise ValueError("unknown original GSHHG region")
        block=geo.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        near=byisland.get((island,"geo_k5"),())
        first=near[0] if len(near)==5 else None
        scopes={
            "all_recorded":full,
            "strict_native":q["native"],
            "strict_native_nonendemic":q["native_nonendemic"],
        }
        for tier,lengths in zip(TIERS,(external,two)):
            predictions=source_predict_same_island(
                full_island=full,plant_scopes=scopes,
                source_geo5=sources["geo5"],
                source_climate10=sources["climate10"],
                target_depth=original_depth,length_cm=lengths,cal=cal,
            )
            for scope,estimate in predictions.items():
                row=dict(
                    island_id=island,region=region,scope=scope,tier=tier,
                    spatial_block=block,first_source=first,
                    n_target_species=estimate["n_same_island_direct_target_species"],
                    n_target_genera=estimate["n_original_target_genera"],
                    island_documentation_fraction=estimate["island_documentation_fraction"],
                )
                for mode in ("geo5","climate10"):
                    z=estimate[mode]
                    row[mode+"__n_numeric_donors"]=z["n_distinct_original_source_numeric_congeners"]
                    row[mode+"__source_numeric_coverage"]=z["source_direct_numeric_coverage"]
                    row[mode+"__brier_genus"]=z["brier_length_genus"]
                    row[mode+"__brier_pooled"]=z["brier_length_pooled"]
                records.append(row)
    if nclimate!=4485:
        raise ValueError("original 4485/4505 climate source assigned botanicals changed")
    if len(records)!=4505*len(SCOPES)*len(TIERS):
        raise AssertionError("original 4505 island status×numeric-source frame collapsed")
    report=aggregate_external_source_numeric_prediction(
        pd.DataFrame.from_records(records),population=4505,bootstrap=999
    )
    report.update(
        calibration_excluding_each_target_botanical_genus=calibration_meta,
        original_nhm_numeric_source_provider_summary=prov,
        n_external_numeric_original_NHM_species=len(external),
        n_external_numeric_two_corpus_consistent_species=len(two),
        n_original_GIFT_native_mainland_source_entities=748,
        n_original_climate_nearest_10_covered_islands=nclimate,
        original_direct_pair_columns_not_required=True,
    )
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--nhm-directory",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    ledger_path=a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ledger_path,flora_path):
        if digest(path)!=PIN[path.name]:
            raise ValueError("original Direct floral trait/GIFT island roster SHA changed")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
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
        a.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    numeric={}
    for provider,name in PROVIDERS.items():
        path=a.nhm_directory/name
        if not path.is_file() or nhm_git_blob_sha(path)!=BLOB_PINS[provider]:
            raise ValueError("original independent NHM corpus file git hash mismatch")
        numeric[provider],_=parse_nhm_numeric_file(path,provider)
    report=fit_original_and_predict(
        ledger,flora,geography,assignments,mainland,numeric
    )
    report.update(
        original_GIFT_direct_ordinal_ledger_SHA256=digest(ledger_path),
        original_GIFT_botanical_flora_SHA256=digest(flora_path),
        original_NHM_source_git_blobs=BLOB_PINS,
        original_names_islands_GIFT_regions_never_released=True,
    )
    raw=json.dumps(report,indent=2,ensure_ascii=False,
                   sort_keys=True,allow_nan=False)
    private=("island_id","accepted_species","work_species","entity_ID",
             "spatial_block","source_signature","source_genus","source_family",
             "latitude","longitude","GBIF_id","photo_id","taxon_id")
    if any('"'+key+'"' in raw for key in private):
        raise ValueError("original botanical taxa/island/source private identifiers leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps({
        "status":report["status"],
        "original_islands":report["n_original_botanical_flora_islands"],
        "calibration_species":report["calibration_excluding_each_target_botanical_genus"][
            "n_source_accepted_species_with_NHM_numeric_and_Direct_depth"],
        "eligible_original_source_panels":sum(
            x["all__spatial_block"]["status"]==
            "EXPLORATORY_EXTERNAL_NUMERIC_SOURCE_PREDICTION"
            for x in report["panels"]
        )
    },sort_keys=True))


if __name__=="__main__":
    main()
