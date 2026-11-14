"""Rights-safe original 928-island genus membership Shapley functional breadth.

Original Direct plant traits and GIFT native continental source lists are read
only in protected source job. Every original plant taxon, source, island,
genus and coordinate label is withheld from publicly uploaded aggregate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.discrete_functional_volume import TRAITS
from fd.genus_shapley_functional_volume import original_genus_shapley_macro


def _complete_original_direct(ledger:pd.DataFrame)->pd.DataFrame:
    required={"accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"}
    if not required.issubset(ledger):
        raise ValueError("Original Direct typed floral evidence columns missing")
    d=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if d.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Source Direct accepted species/trait conflicting")
    traits=d.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    traits=traits.dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(traits)!=2234:
        raise ValueError("Global original Direct triple complete species cohort changed")
    return traits


def _strict_json(obj):
    if isinstance(obj,dict):
        return {str(k):_strict_json(v) for k,v in obj.items()}
    if isinstance(obj,list):
        return [_strict_json(v) for v in obj]
    if isinstance(obj,float) and not math.isfinite(obj):
        return None # unidentifiable is NOT a zero effect
    return obj


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-progressive-root",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    a.add_argument("--rarefactions",type=int,default=40)
    a.add_argument("--reference-null-draws",type=int,default=399)
    a.add_argument("--seed",type=int,default=20261113)
    x=a.parse_args()
    root=x.original_progressive_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    qualifying,receipt=strict_direct_928_flora(ledger,flora)
    if receipt["n_original_eligible_islands"]!=928:
        raise ValueError("Original source 928 eligible islands changed")
    direct=_complete_original_direct(ledger)
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    assign=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    native=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    scores=original_genus_shapley_macro(direct,qualifying,geo,assign,native,
        repeats=x.rarefactions,reference_draws=x.reference_null_draws,
        bootstraps=1999,seed=x.seed)
    assert scores["n_original_eligible_Direct_islands"]==928
    assert scores["n_original_reference_Direct_species"]==2234
    assert scores["n_original_reference_triple_profiles"]==100
    assert scores["method"]["exact_additive_shapley_coverage_proven"] is True
    assert scores["method"]["strict_same_plant_name_source_exclusion"] is True
    assert scores["no_original_island_taxon_species_genus_source_identity_released"] is True
    assert scores["no_actual_pollinator_service_or_historical_decline_tested"] is True
    scores["execution"]={
      "source_original_rights_restricted_private":True,
      "n_repeats":x.rarefactions,
      "n_reference_null_draws":x.reference_null_draws,
      "fixed_seed":x.seed
    }
    scores=_strict_json(scores)
    data=json.dumps(scores,indent=2,sort_keys=True,allow_nan=False)
    forbidden=['"island_id"','"accepted_species"','"work_species"',
       '"entity_ID"','"genus"','"family"','"source_signature"',
       '"spatial_block"','"longitude"','"latitude"']
    if any(k in data for k in forbidden):
        raise ValueError("Original restricted source identifiers present in public aggregate")
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(data+"\n")
    print(json.dumps({
       "status":scores["status"],
       "n_original_macro_islands":scores["n_original_eligible_Direct_islands"],
       "n_source_disjoint_flora_scored":scores["n_original_island_scope_rows_decomposed"],
       "n_regional_shapley_panels":len(scores["all_source_disjoint_flora_genus_shapley_regions"]),
       "matched5_and8_eligibility":[{
           "scope":z["scope"],"region":z["regime"],
           "n_genus_shared_islands":z["n_islands_matching_genus_support"],
           "n_blocks":z["n_spatial_blocks"],
           "n_nearest_GIFT_sources":z["n_first_GIFT_sources"],
           "n_genus":z["matched_distinct_genera_per_island"],
           "supported":z["reportable_under_10islands_5blocks_5sources"]
       } for z in scores["restricted_distinct_genus_source_pairing_panels"]]
    },sort_keys=True))
if __name__=="__main__":
    main()
