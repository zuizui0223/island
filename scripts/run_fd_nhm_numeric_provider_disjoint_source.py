"""Original GIFT 4,505-island five-NHM-named-corpus heldout source test.

An original NHM corpus supplies NUMERIC tube lengths of accepted GIFT
source-only congener species, and the other FOUR supply the *calibration*
lengths. In each target genus, all Direct ordinal tube-depth labels of that
genus are excluded from numerical-length/deep-category calibration. A
released corpus is NOT a unique underlying publication or locality.
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
from fd.nhm_numeric_source_congener_predict import (
    botanical_genus,make_genus_heldout_length_calibrator,
    source_predict_same_island,
)
from fd.nhm_numeric_provider_holdout import (
    SCOPES,REGIMES,PROVIDER_NAMES,summarize_provider_disjoint,
)
from run_fd_independent_nhm_numeric_tube_length import (
    combine_external_provider_lengths,parse_nhm_numeric_file,
)
from fd_nhm_original_provider_audit import PROVIDERS,BLOB_PINS
from fd_nhm_external_flora_preflight import nhm_git_blob_sha

SHA={
 "direct_only_trait_ledger.csv.gz":
 "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
 "chapter1_status_flora.csv.gz":
 "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
SOURCE_MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def hashfile(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def original_provider_model_and_source_support(
    ledger:pd.DataFrame,flora:pd.DataFrame,geo:pd.DataFrame,
    assign:pd.DataFrame,mainland:pd.DataFrame,
    nhm_by_corpus:dict[str,dict[str,float]]
)->dict:
    islands=status_groups(flora)
    if len(islands)!=4505:
        raise ValueError("original GIFT botanical 4505 islands not available")
    if geo.island_id.duplicated().any() or geo.island_id.nunique()<8200:
        raise ValueError("original independently published GSHHG source invalid")
    geography=geo.set_index("island_id").reindex(sorted(islands))
    if geography.analysis_regime.notna().sum()!=4504:
        raise ValueError("original one geographically unresolved botanical flora lost")
    if mainland.duplicated(["entity_ID","work_species"]).any() or mainland.entity_ID.nunique()!=748:
        raise ValueError("original mainland 748 accepted source region sets changed")
    gift=set(mainland.work_species)
    depth=strict_traits(ledger)["tube_depth_class"]
    if set(nhm_by_corpus)!=set(PROVIDER_NAMES):
        raise ValueError("exact five original NHM text corpora required, combined is not sixth")
    summary=[]
    models={}
    for provider in PROVIDER_NAMES:
        other={k:v for k,v in nhm_by_corpus.items() if k!=provider}
        other_lengths=combine_external_provider_lengths(other)[0]
        train=set(other_lengths)&gift&set(depth)
        n_deep=sum(depth[sp]=="deep" for sp in train)
        record=dict(
            donor_numeric_named_corpus=provider,
            original_numeric_NHM_source_corpora_used_for_calibration=[
                q for q in PROVIDER_NAMES if q!=provider
            ],
            n_original_numeric_donor_accepted_species=len(nhm_by_corpus[provider]),
            n_original_donor_GIFT_mainland_accepted_species=len(
                set(nhm_by_corpus[provider])&gift
            ),
            n_different_corpus_train_numeric_Direct_accepted_species=len(train),
            n_other_corpus_source_deep=n_deep,
            n_other_corpus_source_nondeep=len(train)-n_deep,
            donor_same_corpus_excluded_from_calibration=True,
            original_target_genus_all_Direct_depth_excluded_per_prediction=True,
        )
        try:
            underlying,model_meta=make_genus_heldout_length_calibrator(
                other_lengths,depth,gift
            )
        except ValueError:
            models[provider]=None
            record["calibration_status"]="NOT_IDENTIFIABLE_OTHER_FOUR_CORPORA_DEPTH_CLASSES"
            record["n_target_genera_failing_specific_calibration"]=None
        else:
            unsupported=set()
            @lru_cache(maxsize=None)
            def safe_calibrate(g,base=underlying,unsupported=unsupported):
                try:return base(g)
                except ValueError:
                    unsupported.add(g)
                    return None
            models[provider]=(safe_calibrate,unsupported)
            record["calibration_status"]="OTHER_FOUR_CORPORA_LENGTH_CALIBRATION_POSSIBLE"
            record["original_fitted_train_scope"]=model_meta
        summary.append(record)

    src=assign.loc[assign.source_mode.isin(SOURCE_MODES)].copy()
    src.source_rank=pd.to_numeric(src.source_rank,errors="raise").astype(int)
    if src.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("original GIFT source ranking duplicate")
    if {k:int(src.source_mode.eq(k).sum()) for k in SOURCE_MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("original GIFT source ranks drifted")
    for mode,(_,ranks) in SOURCE_MODES.items():
        if not src.loc[src.source_mode.eq(mode),"source_rank"].between(1,ranks).all():
            raise ValueError("original mainland source rank range invalid")
    byisland=src.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    entities=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    @lru_cache(maxsize=20)
    def regional_roster(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(entities.get(v,set()) for v in ids))

    records=[]
    nclimate=0
    for island in sorted(islands):
        info=islands[island]
        full=info["full"]
        sources={}
        for mode,(short,n_ranks) in SOURCE_MODES.items():
            codes=byisland.get((island,mode),())
            if codes and len(codes)!=n_ranks:
                raise ValueError("incomplete original GIFT botanical source shortlist")
            sources[short]=set(regional_roster(codes))-full if codes else set()
            if sources[short]&full:
                raise AssertionError("original focal accepted island name leaked")
        if byisland.get((island,"geo50_climate10")):
            nclimate+=1
        regime=geography.at[island,"analysis_regime"]
        regime=str(regime) if pd.notna(regime) else "geography_unresolved"
        if regime not in REGIMES:
            raise ValueError("original GSHHG island region drift")
        block=geography.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        first=byisland.get((island,"geo_k5"),())
        source_signature=first[0] if len(first)==5 else None
        scopes={
            "all_recorded":full,
            "strict_native":info["native"],
            "strict_native_nonendemic":info["native_nonendemic"],
        }
        for provider in PROVIDER_NAMES:
            obj=models[provider]
            if obj is None:
                predictions=None
            else:
                prediction_calibrator,_=obj
                predictions=source_predict_same_island(
                    full_island=full,plant_scopes=scopes,
                    source_geo5=sources["geo5"],
                    source_climate10=sources["climate10"],
                    target_depth=depth,
                    length_cm=nhm_by_corpus[provider],
                    cal=prediction_calibrator,
                )
            for scope in SCOPES:
                rec=dict(
                    island_id=island,region=regime,
                    scope=scope,provider=provider,
                    spatial_block=block,first_source=source_signature,
                    calibration_supported=obj is not None,
                )
                if predictions is None:
                    rec.update(n_target=0,n_target_genera=0,
                               island_documentation_fraction=np.nan)
                    for mode in ("geo5","climate10"):
                        rec[mode+"__n_numeric_donors"]=0
                        rec[mode+"__source_coverage"]=np.nan
                        rec[mode+"__brier_genus"]=np.nan
                        rec[mode+"__brier_pooled"]=np.nan
                else:
                    q=predictions[scope]
                    rec.update(
                        n_target=q["n_same_island_direct_target_species"],
                        n_target_genera=q["n_original_target_genera"],
                        island_documentation_fraction=q["island_documentation_fraction"],
                    )
                    for mode in ("geo5","climate10"):
                        z=q[mode]
                        rec[mode+"__n_numeric_donors"]=z[
                            "n_distinct_original_source_numeric_congeners"
                        ]
                        rec[mode+"__source_coverage"]=z[
                            "source_direct_numeric_coverage"
                        ]
                        rec[mode+"__brier_genus"]=z["brier_length_genus"]
                        rec[mode+"__brier_pooled"]=z["brier_length_pooled"]
                records.append(rec)
    if nclimate!=4485:
        raise ValueError("original climate-nearest10 source assignment lost")
    if len(records)!=4505*len(SCOPES)*len(PROVIDER_NAMES):
        raise AssertionError("full original GIFT study population contracted")
    for row in summary:
        label=row["donor_numeric_named_corpus"]
        item=models[label]
        row["n_target_genera_failing_specific_calibration"]=(
            len(item[1]) if item is not None else None
        )
    result=summarize_provider_disjoint(
        pd.DataFrame.from_records(records),population=4505,draws=999
    )
    result.update(
        n_original_748_mainland_source_regions=748,
        original_geo5_ranked_source_rows=41325,
        original_climate10_ranked_source_rows=82190,
        n_climate10_covered_original_botanical_floras=nclimate,
        original_source_corpus_calibration_opportunity=summary,
        original_five_corpus_NHM_valid_length_species_counts={
            key:len(value) for key,value in nhm_by_corpus.items()
        },
    )
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--nhm-directory",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    ledger_path=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for f in (ledger_path,flora_path):
        if hashfile(f)!=SHA[f.name]:
            raise ValueError("original botanical Direct trait / GIFT island source hash mismatch")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"
    ])
    geography=pd.read_csv(args.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"
    ])
    assign=pd.read_csv(
        args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"]
    )
    mainland=pd.read_csv(
        args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"]
    )
    raw={}
    for provider,path in PROVIDERS.items():
        file=args.nhm_directory/path
        if not file.is_file() or nhm_git_blob_sha(file)!=BLOB_PINS[provider]:
            raise ValueError("original named NHM corpus Git hash changed")
        raw[provider],_=parse_nhm_numeric_file(file,provider)
    report=original_provider_model_and_source_support(
        ledger,flora,geography,assign,mainland,raw
    )
    report.update(
        original_GIFT_Direct_species_ledger_sha256=hashfile(ledger_path),
        original_GIFT_4505_botanical_flora_sha256=hashfile(flora_path),
        pinned_NHM_five_named_original_git_blobs=BLOB_PINS,
        all_original_island_taxa_and_source_identifiers_remain_private=True,
    )
    out=json.dumps(report,sort_keys=True,ensure_ascii=False,indent=2,
                   allow_nan=False)
    prohibited=("island_id","accepted_species","work_species",
                "entity_ID","source_region_id","source_signature",
                "spatial_block","source_genus","source_family","genus",
                "latitude","longitude","photo_id","GBIF_id","taxon_id")
    if any('"'+x+'"' in out for x in prohibited):
        raise ValueError("original botanical/visitor/source identity leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(out+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=report["status"],n_original_islands=report["n_original_GIFT_island_floras"],
        n_region_status_provider_panels=report["n_regional_status_source_provider_panels"],
        n_original_named_NHM_provider_source_supported=sum(
            z["calibration_status"]=="OTHER_FOUR_CORPORA_LENGTH_CALIBRATION_POSSIBLE"
            for z in report["original_source_corpus_calibration_opportunity"]
        ),
        n_geographically_supported_panels=sum(
            p["all__spatial_block"]["status"]==
            "EXPLORATORY_ORIGINAL_PROVIDER_DISJOINT_NUMERIC_SOURCE_PREDICTION"
            for p in report["panels"]
        ),
    ),sort_keys=True))


if __name__=="__main__":
    main()
