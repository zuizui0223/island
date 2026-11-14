"""Protected original oceanic-vs-continental SAME-GENUS endemic source comparison.

Original trait, endemic/native status, candidate GIFT mainland flora and
geological crosswalk stay inside original source Actions. Only trait×region
×geological origin aggregate comparisons are exported, with no island or
taxonomic identifiers. Source congener != ancestor; no adaptation claim.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd

from fd.geological_endemic_shared_genus import calculate_endemic_shared_genus_geology

def run(original:Path,geology:Path,*,seed:int)->dict:
    ledger=pd.read_csv(
        original/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    raw=pd.read_csv(
        original/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status",
                           "endemic_status","floristic_status"])
    geo=pd.read_csv(geology,dtype=str,usecols=[
        "island_id","geological_origin","analysis_regime","spatial_block"])
    typed=geo.loc[geo.geological_origin.isin(("oceanic","continental"))]
    if (typed.island_id.nunique()<430 or
        typed.geological_origin.value_counts().get("oceanic",0)<200 or
        typed.geological_origin.value_counts().get("continental",0)<150):
        raise ValueError("Original independently typed geological island count drifted")
    source=pd.read_csv(
        original/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        original/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    out=calculate_endemic_shared_genus_geology(
        ledger,raw,geo,source,mainland,seed=seed,draws=1999)
    if out["original_geological_islands_with_valid_labels"]!=len(typed):
        raise ValueError("Original geology direct trait source sample drifted")
    out["source_execution_seed"]=seed
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--seed",type=int,default=20261114)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.original_root,a.geology_crosswalk,seed=a.seed)
    text=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    private=(
        "island_id","accepted_species","origin_status","endemic_status",
        "floristic_status","entity_ID","work_species","spatial_block",
        "source_signature","family","genus","native_l3_codes","longitude","latitude")
    if any('"'+v+'"' in text for v in private):
        raise ValueError("Original plant or island source identity in public output")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":out["status"],
        "seed":a.seed,
        "n_geology_typed":out["original_geological_islands_with_valid_labels"],
        "n_single_trait_source_qualified_rows":out["n_original_distinct_endemic_source_qualified_trait_island_plant_rows"],
        "n_original_regional_source_opportunity":len(out["original_endemic_species_and_independent_genus_opportunity"]),
        "n_geological_shared_genus_trait_tests":len(out["original_shared_genus_geological_comparisons"]),
        "n_reportable_geological_effects":sum(
            x["qualified_to_compare_source_endemic_discordance_within_shared_genera"]
            for x in out["original_shared_genus_geological_comparisons"])
    },sort_keys=True))

if __name__=="__main__":
    main()
