"""Protected original 8k-island all-available-direct-trait floral diversity.

This replaces the 928 complete-three-trait main-sample restriction with
PER-TRAIT data availability. Single-trait and TWO-trait endpoints each use
their own eligible directly observed accepted species; no missingness
imputation and no invented zero. Original full global island botanical
floras, species IDs, GIFT source IDs, and source geometry stay private.

The historical 928>=30-Direct flower-triple subset is check ONLY and
does not gate eligible global islands.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from fd.global_partial_trait_flora_fd import original_global_partial_traits


def run(root:Path,geology:Path,*,seed:int=20261114)->dict:
    ledger=pd.read_csv(
        root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    flora=pd.read_csv(
        root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status",
                           "endemic_status","floristic_status"])
    if flora.island_id.nunique()<4000:
        raise ValueError("Original all-island GIFT flora not loaded")
    geo=pd.read_csv(geology,dtype=str,
        usecols=["island_id","geological_origin","analysis_regime","spatial_block"])
    if (geo.island_id.nunique()<8200 or
        geo.loc[geo.geological_origin.eq("continental")].shape[0]!=176 or
        geo.loc[geo.geological_origin.eq("oceanic")].shape[0]!=258):
        raise ValueError("Original independent geological island crosswalk changed")
    assignment=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    out=original_global_partial_traits(
        ledger,flora,geo,assignment,mainland,
        seed=seed,bootstrap_draws=999)
    if out["n_original_three_trait_complete30_islands_retained_for_SENSITIVITY"]!=928:
        raise ValueError("Original old 928 sensitivity population changed")
    if out["n_original_distinct_islands_with_geological_origin_typed"]<400:
        raise ValueError("Geological class not linkable to original broad botanical floras")
    out["original_numeric_seed"]=seed
    out["original_source_private_rights_restricted_all_island_floras"]=True
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--seed",type=int,default=20261114)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.original_root,a.geology_crosswalk,seed=a.seed)
    body=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    private=(
        "island_id","accepted_species","entity_ID","work_species",
        "source_signature","spatial_block","source_genus","source_family",
        "latitude","longitude","taxon_id","photographer_id","trait_value",
        "geological_record_id","native_l3_codes","tdwg_l3_code")
    if any('"'+k+'"' in body for k in private):
        raise ValueError("Rights-restricted original island/plant source identity escaped")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps({
        "status":out["status"],
        "n_full_original_global_floras":out["n_original_global_distinct_island_floras"],
        "n_original_origin_typed_floras":out["n_original_distinct_islands_with_geological_origin_typed"],
        "n_historic_three_trait_complete_sensitivity":out[
            "n_original_three_trait_complete30_islands_retained_for_SENSITIVITY"],
        "n_global_coverage_trait_panels":len(out["marginal_trait_and_pairwise_coverage"]),
        "n_oceanic_continental_coverage_contrasts":len(
            out["oceanic_minus_continental_conditional_traits"]),
        "n_supported_origin_contrasts":sum(
            x["conditional_sampling_support"]
            for x in out["oceanic_minus_continental_conditional_traits"]),
        "original_multi_trait_complete_NOT_main_criterion":True,
    },sort_keys=True))


if __name__=="__main__":main()
