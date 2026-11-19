"""Original 4,505 GIFT island source versus independently published NHM corolla tube length.

NHM numeric *corolla tube min/max length [cm]* was previously ignored by the
fd colour/form/symmetry scripts. This audits actual species-level overlap
against original categorical Direct tube-depth and source-region opportunity.
NHM text-extracted lengths are NOT independent field measurements of island
populations, not the same construct as an ordinal 'deep' state, and not a
bird/butterfly pollination assignment. Original taxon/island IDs stay in CI.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import json
import re

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd_nhm_external_flora_preflight import (
    accepted_exact_binomial,nhm_git_blob_sha, NHM_BLOB,
)
from fd_nhm_original_provider_audit import PROVIDERS,BLOB_PINS
from fd.source_genus_same_island import genus

LENGTH_COLS=("corolla tube min. length [cm]",
             "corolla tube max. length [cm]")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}
TIER=("any_valid_NHM_numeric_tube_length",
      "two_source_corpora_consistent_tube_length")
PINNED_SHA={
    "direct_only_trait_ledger.csv.gz":
      "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
      "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}


def sha(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for s in iter(lambda:f.read(1<<20),b""):
            h.update(s)
    return h.hexdigest()


def strict_single_length(vmin:str,vmax:str)->float|None:
    """Reject unknown ranges, impossible min>max and nonfinite/OCR junk."""
    def parse(v):
        if v is None:return None
        t=str(v).strip()
        # NHM original values are e.g. '0.5 cm'; strip ONLY the
        # explicitly recorded centimetre suffix, never infer other units.
        t=re.sub(r"\s*cm\s*$","",t,flags=re.I)
        if not re.fullmatch(r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",t):
            return None
        z=float(t)
        return z if 0<z<=100 and np.isfinite(z) else None
    x,y=parse(vmin),parse(vmax)
    if x is not None and y is not None:
        if x>y:return None
        return (x+y)/2
    return x if x is not None else y


def parse_nhm_numeric_file(path:Path,provider:str)->tuple[dict,dict]:
    columns=set(pd.read_csv(path,nrows=0).columns)
    if "taxon" not in columns:raise ValueError("NHM provider lacks plant taxon field")
    available=[v for v in LENGTH_COLS if v in columns]
    if not available:
        return {},dict(provider=provider,n_numeric_cols=0,
                       n_raw_rows_with_numeric_tube_length=0,
                       n_unique_exact_names_with_numeric_length=0,
                       n_invalid_numeric_range_rows=0)
    use=["taxon",*available]
    raw=pd.read_csv(path,dtype=str,keep_default_na=False,
                    low_memory=False,usecols=use)
    for c in LENGTH_COLS:
        if c not in raw:raw[c]=""
    states=defaultdict(list)
    valid_rows=0;invalid=0
    for name,minv,maxv in raw[["taxon",*LENGTH_COLS]].itertuples(
        index=False,name=None
    ):
        sp=accepted_exact_binomial(name)
        if sp is None:continue
        z=strict_single_length(minv,maxv)
        if z is None:
            if str(minv).strip() or str(maxv).strip():invalid+=1
            continue
        valid_rows+=1
        states[sp].append(z)
    out={s:float(np.median(a)) for s,a in states.items()}
    return out,dict(
        provider=provider,n_numeric_cols=len(available),
        n_raw_rows_with_numeric_tube_length=valid_rows,
        n_unique_exact_names_with_numeric_length=len(out),
        n_invalid_numeric_range_rows=invalid,
        n_names_with_multiple_numeric_entries=sum(len(a)>1 for a in states.values()),
    )


def combine_external_provider_lengths(raw:dict[str,dict])->tuple[dict,dict,dict]:
    byname=defaultdict(dict)
    for provider,source in raw.items():
        for sp,value in source.items():
            if not 0<value<=100 or not np.isfinite(value):
                raise ValueError("source numeric corolla tube out of original units")
            byname[sp][provider]=float(value)
    pooled={}
    two={}
    conflicts=0
    for sp,records in byname.items():
        numbers=np.array(list(records.values()),float)
        pooled[sp]=float(np.median(numbers))
        if len(records)>=2:
            ratio=float(numbers.max()/numbers.min())
            if ratio<=1.5:two[sp]=pooled[sp]
            else:conflicts+=1
    return pooled,two,dict(
        n_any_original_NHM_accepted_names=len(pooled),
        n_original_names_recorded_by_two_consistent_corpora=len(two),
        n_multi_corpus_gt_1p5_max_to_min_conflicts=conflicts,
        n_external_source_corpora=len(raw),
    )


def study(ledger:pd.DataFrame,flora:pd.DataFrame,geo:pd.DataFrame,
          mainland:pd.DataFrame,assign:pd.DataFrame,
          numeric_provider:dict[str,dict],provider_info:list)->dict:
    original=status_groups(flora)
    if len(original)!=4505:raise ValueError("original 4505 GIFT flora universe narrowed")
    if geo.island_id.duplicated().any() or geo.island_id.nunique()<8200:
        raise ValueError("original independent GSHHG island geography incomplete")
    G=geo.set_index("island_id").reindex(sorted(original))
    if int(G.analysis_regime.notna().sum())!=4504:
        raise ValueError("GSHHG source geography changed")
    if mainland.duplicated(["entity_ID","work_species"]).any() or (
        mainland.entity_ID.nunique()!=748):
        raise ValueError("original 748 mainland GIFT region memberships changed")
    traits=strict_traits(ledger)
    ordinal=traits["tube_depth_class"]
    ext,two,info=combine_external_provider_lengths(numeric_provider)
    accepted=set(mainland.work_species)
    all_island=set().union(*(q["full"] for q in original.values()))
    n_intersection=set(ext)&set(ordinal)
    known_depths=defaultdict(list)
    for s in n_intersection:
        known_depths[ordinal[s]].append(ext[s])
    direct_comparison=[
        dict(
            original_Direct_ordinal_depth_category=k,
            n_independently_released_NHM_numeric_lengths=len(values),
            observed_original_NHM_tube_length_cm_median=float(np.median(values)),
            q25_q75_cm=[float(v) for v in np.quantile(values,[.25,.75])],
        ) for k,values in sorted(known_depths.items())
    ]
    tiers={TIER[0]:ext,TIER[1]:two}
    external_source_audit=[]
    for key,vals in tiers.items():
        ids=(set(vals)&accepted)
        external_source_audit.append(dict(
            independent_NHM_numeric_provider_tier=key,
            n_original_GIFT_mainland_species_with_numeric_length=len(ids),
            n_original_GIFT_mainland_species_not_on_ANY_GIFT_island=len(ids-all_island),
            n_original_GIFT_mainland_species_in_any_original_island=len(ids&all_island),
            n_original_GIFT_mainland_species_also_Direct_tube_depth=len(ids&set(ordinal)),
            n_original_numeric_source_only_mainland_genera=len({genus(s) for s in ids}),
            original_length_units="cm",
        ))

    original_assign=assign.loc[assign.source_mode.isin(MODES)].copy()
    original_assign.source_rank=pd.to_numeric(original_assign.source_rank,errors="raise").astype(int)
    if original_assign.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("original source regional rankings duplicated")
    if {mode:int(original_assign.source_mode.eq(mode).sum()) for mode in MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("original GIFT source name opportunity changed")
    byisland=original_assign.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    for mode,(_,max_rank) in MODES.items():
        if not original_assign.loc[
            original_assign.source_mode.eq(mode),"source_rank"
        ].between(1,max_rank).all():
            raise ValueError("incomplete original GIFT candidate source rank")

    eligible=(set(ext)|set(two))&accepted
    region_names=mainland.loc[mainland.work_species.isin(eligible)].groupby(
        "entity_ID"
    ).work_species.agg(set).to_dict()

    counters=defaultdict(Counter)
    blocks=defaultdict(set)
    firstsources=defaultdict(set)
    for island,q in original.items():
        reg=G.at[island,"analysis_regime"]
        region=str(reg) if pd.notna(reg) else "geography_unresolved"
        if region not in REGIMES:
            raise ValueError("unknown GSHHG original island regime")
        block=G.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        candidate={}
        for mode,(label,nranks) in MODES.items():
            source=byisland.get((island,mode),())
            if source and len(source)!=nranks:
                raise ValueError("partial original source region ranks")
            # Exclude ONLY focal accepted species, unlike the deliberately
            # strict global island-species-excluded sensitivity.
            candidate[label]=set().union(*(region_names.get(v,set()) for v in source))-q["full"]
            if candidate[label]&q["full"]:
                raise AssertionError("focal botanical accepted name leaked into source")
        for scope,names in (
            ("all_recorded",q["full"]),("strict_native",q["native"]),
            ("strict_native_nonendemic",q["native_nonendemic"]),
        ):
            known_target=names&set(ordinal)
            nt=len(known_target)
            for tier,values in tiers.items():
                key=(region,scope,tier)
                c=counters[key]
                c["n_original_region_islands"]+=1
                if nt>=5:c["n_islands_at_least_five_original_Direct_depth_species"]+=1
                donors={mode:candidate[mode]&set(values)
                        for mode in ("geo5","climate10")}
                for mode in ("geo5","climate10"):
                    count=len(donors[mode])
                    for n in (1,5,10):
                        if count>=n:c[f"n_source_{mode}_NHM_length_species_ge_{n}"]+=1
                genera_geo={genus(s) for s in donors["geo5"]}
                genera_clim={genus(s) for s in donors["climate10"]}
                target_with_congener={
                    sp for sp in known_target
                    if genus(sp) in (genera_geo & genera_clim)
                }
                if nt>=5 and len(donors["geo5"])>=10 and (
                    len(donors["climate10"])>=10):
                    c["n_source_mode_eligible_original_floras"]+=1
                if (len(target_with_congener)>=5 and
                    len({genus(x) for x in target_with_congener})>=2 and
                    len(donors["geo5"])>=10 and
                    len(donors["climate10"])>=10):
                    c["n_supported_same_island_source_congener_comparisons"]+=1
                    if block is not None:blocks[key].add(block)
                    source=byisland.get((island,"geo_k5"),())
                    if source:firstsources[key].add(source[0])
    summary=[]
    for r in REGIMES:
        for scope in SCOPES:
            for tier in TIER:
                k=(r,scope,tier)
                c=counters[k]
                n=c["n_supported_same_island_source_congener_comparisons"]
                nblocks=len(blocks[k]);ngroups=len(firstsources[k])
                summary.append(dict(
                    region=r,scope=scope,NHM_external_numeric_tier=tier,
                    **{name:int(c[name]) for name in (
                        "n_original_region_islands",
                        "n_islands_at_least_five_original_Direct_depth_species",
                        "n_source_geo5_NHM_length_species_ge_1",
                        "n_source_geo5_NHM_length_species_ge_5",
                        "n_source_geo5_NHM_length_species_ge_10",
                        "n_source_climate10_NHM_length_species_ge_1",
                        "n_source_climate10_NHM_length_species_ge_5",
                        "n_source_climate10_NHM_length_species_ge_10",
                        "n_source_mode_eligible_original_floras",
                        "n_supported_same_island_source_congener_comparisons",
                    )},
                    n_original_independent_10deg_spatial_blocks=nblocks,
                    n_original_first_source_clusters=ngroups,
                    status=("SOURCE_ONLY_EXTERNAL_NUMERIC_DEPTH_OPPORTUNITY"
                            if n>=30 and nblocks>=10 and ngroups>=5
                            else "NOT_IDENTIFIABLE_FOR_CLUSTERED_EXTERNAL_SOURCE_PREDICTION")
                ))
    return dict(
        status="ORIGINAL_4505_INDEPENDENT_NHM_NUMERIC_COROLLA_TUBE_LENGTH_PREFLIGHT",
        n_original_botanical_flora_islands=4505,
        n_original_mainland_GIFT_source_entities=748,
        original_five_NHM_raw_provider_tables=provider_info,
        original_external_species_numeric_summary=info,
        original_external_numeric_depth_and_Direct_ordinal_overlap=direct_comparison,
        original_source_accepted_name_opportunity=external_source_audit,
        n_island_scope_provider_region_panels=len(summary),
        island_source_numeric_depth_support=summary,
        n_analytically_eligible_region_scopes=sum(
            r["status"]=="SOURCE_ONLY_EXTERNAL_NUMERIC_DEPTH_OPPORTUNITY" for r in summary
        ),
        scientific_limits=dict(
            independently_published_NHM_ocr_numeric_length_not_same_as_focal_local_population_flower=True,
            tube_length_cm_not_equipollent_to_original_ordinal_depth_without_calibration=True,
            source_species_absent_from_focal_original_island_checklist=True,
            original_all_4505_study_frame_kept=True,
            requires_plant_specific_independent_visitors_to_label_syndrome=True,
            no_source_biological_arrival_establishment_or_island_trait_evolution_inferred=True,
            no_new_effect_or_causal_fit_until_source_overlap_audited=True,
            independent_released_text_corpus_not_independent_original_studies=True,
            combined_numeric_NHM_release_lacks_original_study_identification=True,
            combined_release_excluded_from_two_original_provider_voting=True,
            constituent_five_corpora_and_combined_numeric_release_separate=True,
        ),
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--geology-crosswalk",required=True,type=Path)
    p.add_argument("--nhm-directory",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    ledger_path=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for file in (ledger_path,flora_path):
        if sha(file)!=PINNED_SHA[file.name]:
            raise ValueError("original frozen GIFT accepted plant source SHA drift")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    geo=pd.read_csv(args.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assignments=pd.read_csv(
        args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    numeric_provider={};metadata=[]
    for label,name in PROVIDERS.items():
        file=args.nhm_directory/name
        if not file.is_file() or nhm_git_blob_sha(file)!=BLOB_PINS[label]:
            raise ValueError("NHM independent released source git blob changed")
        numeric_provider[label],meta=parse_nhm_numeric_file(file,label)
        metadata.append(meta)
    # All five original constituent CSVs contain explicitly labelled
    # 'x cm' corolla-tube lengths. The separately released combined.csv
    # reproduces/aggregates them; NEVER add it as a sixth independent
    # 'provider' or manufacture a two-provider consensus from duplicates.
    combined=args.nhm_directory/"combined.csv"
    if not combined.is_file() or nhm_git_blob_sha(combined)!=NHM_BLOB:
        raise ValueError("NHM combined numeric corolla tube Git blob changed")
    combined_values,meta=parse_nhm_numeric_file(
        combined,"NHM_combined_opaque_numeric_release")
    metadata.append(meta)
    original_union,two_original,_=combine_external_provider_lengths(
        numeric_provider
    )
    paired_overlap=set(original_union)&set(combined_values)
    numeric_reconciliation=dict(
        n_original_constituent_provider_union=len(original_union),
        n_combined_report_accepted_species=len(combined_values),
        n_exact_name_overlap=len(paired_overlap),
        n_original_constituent_only=len(set(original_union)-set(combined_values)),
        n_combined_only=len(set(combined_values)-set(original_union)),
        n_original_real_two_provider_consensus_species=len(two_original),
        n_overlap_length_disagreement_gt50pct=sum(
            max(original_union[sp],combined_values[sp])/
            min(original_union[sp],combined_values[sp])>1.5
            for sp in paired_overlap),
        combined_counts_only_a_checksum_not_an_extra_source_provider=True,
    )
    output=study(ledger,flora,geo,mainland,assignments,numeric_provider,metadata)
    output["combined_vs_five_independent_corpora_numeric_receipt"]=numeric_reconciliation
    output["original_Direct_ledger_SHA256"]=sha(ledger_path)
    output["original_GIFT_4505_flora_SHA256"]=sha(flora_path)
    output["NHM_pinned_git_corpus_blobs"]={
        **BLOB_PINS,"NHM_combined_opaque_numeric_release":NHM_BLOB
    }
    output["combined_numeric_release_not_original_provider_independence"]=True
    output["real_numeric_source_corpora_for_consensus"]=len(numeric_provider)
    output["not_six_corpora_after_combined_release_deduplication"]=True
    raw=json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)
    forbidden=("island_id","accepted_species","work_species","entity_ID",
               "spatial_block","source_signature","source_region_id",
               "source_genus","source_family","latitude","longitude",
               "GBIF_id","photo_id","taxon_id")
    if any('"'+k+'"' in raw for k in forbidden):
        raise ValueError("original plant or mainland source identifiers leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=output["status"],original_GIFT_floras=4505,
        n_NHM_provider_numeric_names=output["original_external_species_numeric_summary"][
            "n_any_original_NHM_accepted_names"],
        eligible_original_region_status_panels=output["n_analytically_eligible_region_scopes"],
    ),sort_keys=True))


if __name__=="__main__":
    main()
