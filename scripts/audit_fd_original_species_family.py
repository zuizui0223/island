"""Detect a primary source-backed plant family map for original fd Direct taxa.

Runs only in the rights-protected island CI with the frozen progressive
checkpoint. Emits ONLY anonymous file-level eligibility statistics and
taxonomic conflict counts. No original taxon names, genus IDs, site IDs,
coordinates or original taxon×family mapping leave the protected job.

This is a provenance/positivity audit; never infer family from genus string.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import pandas as pd


def read_headers(path:Path):
    try:
        if path.name.endswith(".tsv") or path.name.endswith(".tsv.gz"):
            sep="\t"
        else:
            sep=","
        cols=list(pd.read_csv(path,nrows=0,sep=sep,low_memory=False).columns)
        return cols,sep
    except (UnicodeError,ValueError,EOFError,OSError):
        return [],","


def audit(root:Path)->dict:
    anchor=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    ledger=pd.read_csv(anchor,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"],dtype=str)
    x=ledger.loc[
        ledger.trait_name.isin(("flower_primary_color","floral_form","floral_symmetry"))&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.str.contains("|",regex=False)
    ].copy()
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original source direct taxon/trait duplication")
    xt=x.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    xt=xt.dropna(subset=[
       "flower_primary_color","floral_form","floral_symmetry"])
    if len(xt)!=2234:
        raise ValueError("Original frozen 2234 Direct species universe drifted")
    names=set(xt.index)
    candidates=[]
    all_paths=sorted(p for p in root.rglob("*") if p.is_file() and
        (p.name.endswith(".csv") or p.name.endswith(".csv.gz") or
         p.name.endswith(".tsv") or p.name.endswith(".tsv.gz")))
    for path in all_paths:
        header,sep=read_headers(path)
        key=("accepted_species" if "accepted_species" in header else None)
        fam=("family" if "family" in header else
            "accepted_family" if "accepted_family" in header else
            "family_name" if "family_name" in header else None)
        if not key or not fam:continue
        # Require taxon-family mapping with native master metadata, not
        # pollinator/visitor research rows; all are evaluated for consistency
        # but only the full/unique taxonomic mapping can be chosen.
        try:
            frame=pd.read_csv(path,usecols=[key,fam],sep=sep,dtype=str)
        except Exception:
            continue
        frame=frame.loc[
            frame[key].isin(names)&frame[fam].notna()].copy()
        frame[fam]=frame[fam].str.strip()
        frame=frame.loc[frame[fam].ne("")&~frame[fam].str.lower().isin(
            ("unknown","unresolved","na","n/a","none"))]
        if frame.empty:continue
        conflict=(frame.groupby(key)[fam].nunique()>1)
        bad=set(conflict[conflict].index)
        frame=frame.loc[~frame[key].isin(bad)]
        collapsed=frame.drop_duplicates(key)
        cov=int(collapsed[key].nunique())
        rec=dict(
            source_basename=path.name,
            # Only source basename, NEVER private original local path.
            n_recorded_source_files_with_species_family_columns=1,
            n_2234_original_Direct_species_family_qualified=cov,
            n_taxa_with_conflicting_family_labels=len(bad),
            n_matched_original_source_rows=len(frame),
            n_distinct_families=int(collapsed[fam].nunique()),
            n_distinct_genera=(int(collapsed[key].str.split().str[0].nunique())),
            source_mapping_is_unvalidated_for_family_phylogenetic_tree=True,
        )
        candidates.append(rec)
    candidates.sort(key=lambda p:(
       -p["n_2234_original_Direct_species_family_qualified"],
       p["n_taxa_with_conflicting_family_labels"],p["source_basename"]))
    chosen=candidates[0] if candidates else None
    return {
      "status":("SOURCE_FAMILY_MAPPING_CANDIDATE_FOUND_NOT_YET_AUTHORITATIVE"
                if chosen else "NO_SOURCE_FAMILY_MAPPING_AVAILABLE"),
      "n_original_Direct_complete_species":2234,
      "n_original_source_csv_tsv_files_scanned":len(all_paths),
      "n_species_family_table_candidates":len(candidates),
      "largest_source_candidates":candidates[:8],
      "best_candidate_basename_only":chosen["source_basename"] if chosen else None,
      "n_Direct_species_with_source_family_if_best":(
         chosen["n_2234_original_Direct_species_family_qualified"]
         if chosen else 0),
      "source_family_mapping_not_yet_used_to_fit_biological_effect":True,
      "source_family_names_and_taxon_map_never_exported":True,
      "genus_names_not_used_to_guess_family":True,
      "phylogenetic_tree_not_reconstructed":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=audit(a.root)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    txt=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for key in ('"accepted_species"','"island_id"','"taxon"','"genus"',
                '"family"','"source_entity"','"GIFT_entity_ID"'):
        if key in txt:
            raise ValueError("Private source identity leaked to family-map preflight")
    a.output.write_text(txt+"\n")
    print(json.dumps({
       "status":result["status"],
       "n_files_scanned":result["n_original_source_csv_tsv_files_scanned"],
       "n_candidates":result["n_species_family_table_candidates"],
       "n_Direct_taxa_with_verified_family_candidate":result["n_Direct_species_with_source_family_if_best"],
       "best_source_basename":result["best_candidate_basename_only"],
    },sort_keys=True))
if __name__=="__main__":
    main()
