"""Original-data: WCVP native floral pairing with exact genus/family controls.

Protected original island flora, GIFT donors, WCVP candidate nativeness and
CoL exact family candidate source blobs are hash-verified in the workflow.
Only coarse region×provenance×taxonomic-depth aggregate is released.

This does NOT infer realized pollination, historical partner loss, within-
species evolutionary flower phenotype change, nor that WCVP TDWG-L3 native
range proves local island nativeness.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.regional_native_flora_fd import TRAITS
from fd.regional_native_provenance_volume import regional_native_candidates
from fd.colxr_family_opportunity import strict_candidate_family_map
from fd.native_taxon_matched_pairing import compare_native_source_taxon_pairing


EXPECTED_COLXR_SHA="079d9c27b0ca3aa51f0c7310f43ad4901b07a319d6063eab5a7d62b291bdd9e3"


def verify_file(path: Path, expected: str):
    if not path.is_file():
        raise ValueError("Missing frozen protected taxonomic original source file")
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1<<20), b""):
            h.update(part)
    if h.hexdigest()!=expected:
        raise ValueError("Protected original source SHA failed")


def run(root:Path, wcvp_status:Path, island_tdwg:Path, colxr_crosswalk:Path,
        *, seed:int)->dict:
    verify_file(colxr_crosswalk,EXPECTED_COLXR_SHA)
    df=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    original=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status"])
    eligible,gates=strict_direct_928_flora(df,original)
    if gates["n_original_eligible_islands"]!=928:
        raise ValueError("Original eligible flora changed from 928 islands")
    from fd.regional_native_flora_fd import TRAITS
    rec=df.loc[df.trait_name.isin(TRAITS)&
               df.evidence_scope.eq("species_direct")&
               df.evidence_quality.isin(("high","medium"))&
               df.normalized_value.notna()&
               ~df.normalized_value.astype(str).str.contains("|",regex=False)]
    if rec.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Duplicate original species Direct trait records")
    direct=rec.pivot(index="accepted_species",columns="trait_name",
                    values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(direct)!=2234:
        raise ValueError("Original complete species Direct reference changed")
    ranges=pd.read_csv(wcvp_status,dtype=str,
        usecols=["accepted_species","native_l3_codes"]).drop_duplicates()
    tdwg=pd.read_csv(island_tdwg,dtype=str,
        usecols=["island_id","tdwg_l3_code","tdwg_match_status"])
    if len(tdwg)!=8265 or len(ranges.accepted_species.unique())!=55137:
        raise ValueError("Original independent WCVP 2026 source status drift")
    flora=regional_native_candidates(eligible,tdwg,ranges)
    if flora.island_id.nunique()!=928:
        raise ValueError("Original WCVP candidate source mapping dropped islands")
    native=flora.groupby("island_id")[
        ["strict_source_native","regional_native_compatible"]].sum()
    if int(native.strict_source_native.ge(30).sum())!=87:
        raise ValueError("Original 87 strict-source-native islands changed")
    if int(native.regional_native_compatible.ge(30).sum())!=357:
        raise ValueError("Original 357 WCVP compatible floral island cohort changed")
    crosswalk=pd.read_csv(colxr_crosswalk,dtype=str)
    families,family_audit=strict_candidate_family_map(crosswalk,set(direct.index))
    if len(families)<1900 or family_audit["n_candidate_families"]<50:
        raise ValueError("Unreviewed original CoL family candidate evidence inadequate")
    geography=pd.read_csv(
        root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        dtype=str,usecols=["island_id","analysis_regime","spatial_block"])
    geography=geography.loc[geography.island_id.isin(flora.island_id)].drop_duplicates("island_id")
    if len(geography)!=928 or geography.spatial_block.isna().any():
        raise ValueError("Original geography completeness requirement failed")
    src=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    src=src.loc[src.island_id.isin(flora.island_id)&
                src.source_mode.eq("geo_k5")].copy()
    src["source_rank"]=pd.to_numeric(src.source_rank,errors="raise")
    src=src.loc[src.source_rank.between(1,5)]
    if src.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original source candidates are not unique")
    near=src.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(near.get(island,()))!=5 or len(set(near[island]))!=5
           for island in flora.island_id.unique()):
        raise ValueError("Missing original 5 GIFT native regional source IDs")
    donors=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    donors=donors.loc[donors.work_species.isin(direct.index)].drop_duplicates()
    by_source=donors.groupby("entity_ID").work_species.agg(set).to_dict()
    output=compare_native_source_taxon_pairing(
        direct,flora,geography,near,by_source,families,
        repeats=40,seed=seed,boots=1999)
    output["source_qualification"]={
        "original_CoLXR_taxon_candidate_not_human_reviewed":True,
        "original_CoLXR_candidate_names":len(families),
        "original_CoLXR_candidate_families":family_audit["n_candidate_families"],
        "original_CoLXR_crosswalk_SHA256":EXPECTED_COLXR_SHA,
        "original_source_native_status_is_not_reconstructed_ancestry":True,
        "original_WCVP_native_TDWG_is_not_local_island_nativeness":True,
    }
    output["seed"]=int(seed)
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--wcvp-status",type=Path,required=True)
    p.add_argument("--island-tdwg",type=Path,required=True)
    p.add_argument("--colxr-crosswalk",type=Path,required=True)
    p.add_argument("--seed",type=int,default=20261114)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    report=run(a.original_root,a.wcvp_status,a.island_tdwg,a.colxr_crosswalk,
               seed=a.seed)
    if report["n_original_island_floras"]!=928:
        raise ValueError("Islands were silently lost")
    text=json.dumps(report,sort_keys=True,indent=2,ensure_ascii=False,
                    allow_nan=False)
    confidential=(
        "island_id","accepted_species","entity_ID","work_species",
        "spatial_block","source_signature","native_l3_codes",
        "tdwg_l3_code","source_genus","source_family","family_map",
        "candidate_family","submitted_name","longitude","latitude")
    if any('"'+col+'"' in text for col in confidential):
        raise ValueError("Protected original plant, island or source ID escaped")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":report["status"],
        "seed":a.seed,
        "n_report_panels":len(report["regional_taxon_matched_panels"]),
        "reported_scope_rank_n":[
            (x["scope"],x["regime"],x["level"],x["n_matched_distinct_groups"],
             x["n_taxon_matchable_source_comparable_islands"],
             x["meets_pre_outcome_minimum_reporting_gate"])
            for x in report["regional_taxon_matched_panels"]],
    },ensure_ascii=False))


if __name__=="__main__":
    main()
