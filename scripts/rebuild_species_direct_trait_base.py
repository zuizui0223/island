#!/usr/bin/env python3
"""Rebuild the v2 trait base from source-backed species-level records only."""
from __future__ import annotations
import argparse, csv, gzip, json
from collections import defaultdict
from pathlib import Path

TARGET_TRAITS=("flower_primary_color","floral_form","floral_symmetry","pollination_functional_guild","pollen_vector_mode","sex_system","mating_system","self_incompatibility","autonomous_selfing_capacity")
BIEN_MAP={"flower color":"flower_primary_color","flower pollination syndrome":"pollen_vector_mode","whole plant sexual system":"sex_system"}
AUSTRAITS_MAP={"flower_colour":"flower_primary_color","perianth_colour":"flower_primary_color","flower_perianth_symmetry":"floral_symmetry","flower_structural_sex_type":"sex_system"}
FIELDS=["accepted_species","trait_name","candidate_value","source_name","source_trait","source_url","source_record_id","source_excerpt","evidence_scope","evidence_status"]

def norm(value: object) -> str: return " ".join(str(value or "").strip().split())

def master_names(path: Path) -> list[str]:
    with path.open(newline="",encoding="utf-8-sig") as handle:
        reader=csv.DictReader(handle)
        if not reader.fieldnames or "accepted_species" not in reader.fieldnames: raise ValueError("accepted_species missing from master")
        return sorted({norm(row.get("accepted_species")) for row in reader if norm(row.get("accepted_species"))})

def iter_csv(path: Path):
    if path.suffix==".gz": handle=gzip.open(path,"rt",encoding="utf-8",newline="")
    else: handle=path.open("r",encoding="utf-8",newline="")
    with handle: yield from csv.DictReader(handle)

def add(evidence,by_source,master,row,source,target,value,raw_trait):
    species=norm(row.get("accepted_species") or row.get("scrubbed_species_binomial") or row.get("taxon_name"))
    if not species or species not in master or target not in TARGET_TRAITS or not value: return
    by_source[source].add(species)
    record={
        "accepted_species":species,"trait_name":target,"candidate_value":value,"source_name":source,
        "source_trait":raw_trait,"source_url":norm(row.get("source_url") or row.get("url_source")),
        "source_record_id":norm(row.get("source_record_id") or row.get("id")),
        "source_excerpt":norm(row.get("source_excerpt")),"evidence_scope":"species_direct",
        "evidence_status":norm(row.get("evidence_status")) or "source_backed_unreviewed",
    }
    key=(species,target,value,source,raw_trait,record["source_url"],record["source_excerpt"])
    evidence[key]=record

def collect(root: Path, master: set[str]):
    evidence={}; by_source=defaultdict(set)
    for path in root.rglob("master_matched_trait_records*.csv.gz"):
        for row in iter_csv(path):
            raw=norm(row.get("bien_query_trait") or row.get("trait_name")); target=BIEN_MAP.get(raw.lower()); value=norm(row.get("trait_value"))
            if target: add(evidence,by_source,master,row,"BIEN",target,value,raw)
    for path in root.rglob("austraits_all.csv"):
        for row in iter_csv(path):
            raw=norm(row.get("trait_name")); target=AUSTRAITS_MAP.get(raw); value=norm(row.get("trait_value"))
            if target: add(evidence,by_source,master,row,"AusTraits",target,value,raw)
    for path in root.rglob("efloras_species_direct_candidates.csv"):
        for row in iter_csv(path):
            target=norm(row.get("trait_name")); value=norm(row.get("candidate_value")); raw=norm(row.get("source_trait"))
            add(evidence,by_source,master,row,"eFloras",target,value,raw)
    return list(evidence.values()),by_source

def main() -> int:
    parser=argparse.ArgumentParser(); parser.add_argument("--artifact-root",type=Path,required=True); parser.add_argument("--master-csv",type=Path,required=True); parser.add_argument("--output-dir",type=Path,required=True); args=parser.parse_args()
    names=master_names(args.master_csv); master=set(names); rows,by_source=collect(args.artifact_root,master); args.output_dir.mkdir(parents=True,exist_ok=True)
    with gzip.open(args.output_dir/"species_direct_evidence.csv.gz","wt",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=FIELDS); writer.writeheader(); writer.writerows(sorted(rows,key=lambda r:(r["accepted_species"],r["trait_name"],r["source_name"],r["candidate_value"])))
    cells={(r["accepted_species"],r["trait_name"]) for r in rows}; covered_species={s for s,_ in cells}; by_trait={trait:len({s for s,t in cells if t==trait}) for trait in TARGET_TRAITS}
    with gzip.open(args.output_dir/"species_trait_gap_ledger.csv.gz","wt",encoding="utf-8",newline="") as handle:
        writer=csv.writer(handle); writer.writerow(["accepted_species","trait_name","direct_evidence_present","acquisition_status"])
        for species in names:
            for trait in TARGET_TRAITS:
                present=(species,trait) in cells; writer.writerow([species,trait,str(present).lower(),"covered" if present else "pending_direct_search"])
    report={"master_species":len(names),"target_traits":len(TARGET_TRAITS),"direct_evidence_rows":len(rows),"direct_species_trait_cells":len(cells),"species_with_any_mapped_target_trait":len(covered_species),"species_with_any_mapped_target_trait_rate":len(covered_species)/len(names) if names else 0,"species_by_source_for_mapped_targets":{k:len(v) for k,v in sorted(by_source.items())},"species_direct_by_target_trait":by_trait,"policy":{"taxonomic_resolution":"species_direct_only","inference_included":False,"rule_extracted_flora_evidence":"unreviewed with source excerpt retained","missing_cells":"pending_direct_search; never biological absence"}}
    (args.output_dir/"species_direct_coverage.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8"); print(json.dumps(report,indent=2,sort_keys=True)); return 0
if __name__=="__main__": raise SystemExit(main())
