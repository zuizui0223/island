"""Protected source-replay for 928-island rare versus common floral metric volume.

Never push original species, genus, GIFT source entity, provider and island
identities. Only anonymous blocks, region aggregates, source-eligibility.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.rare_profile_decomposition import original_rare_category_coverage
from fd.discrete_functional_volume import TRAITS


def direct_triples(direct:pd.DataFrame)->pd.DataFrame:
    required={"accepted_species","trait_name","normalized_value",
              "evidence_scope","evidence_quality"}
    if not required.issubset(direct):
        raise ValueError("Original Direct ledger columns changed")
    sample=direct.loc[
       direct.trait_name.isin(TRAITS) &
       direct.evidence_scope.eq("species_direct") &
       direct.evidence_quality.isin(("high","medium")) &
       direct.normalized_value.notna() &
       ~direct.normalized_value.astype(str).str.contains("|",regex=False)
    ].copy()
    if sample.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original source duplicate species-trait evidence")
    result=sample.pivot(index="accepted_species",columns="trait_name",
        values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(result)!=2234:
        raise ValueError("Original global Direct 2234 species drifted")
    return result


def main()->None:
    arg=argparse.ArgumentParser()
    arg.add_argument("--progressive-root",type=Path,required=True)
    arg.add_argument("--output",type=Path,required=True)
    arg.add_argument("--seed",type=int,default=20261113)
    arg.add_argument("--rarefactions",type=int,default=40)
    arg.add_argument("--null-draws",type=int,default=399)
    args=arg.parse_args()
    base=args.progressive_root
    direct=pd.read_csv(base/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    flora=pd.read_csv(base/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    eligible,gate=strict_direct_928_flora(direct,flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original 928-island flower trio cohort drifted")
    triples=direct_triples(direct)
    geo=pd.read_csv(base/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    source=pd.read_csv(base/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(base/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    result=original_rare_category_coverage(
       triples,eligible,geo,source,mainland,
       n=30,repeats=args.rarefactions,
       null_draws=args.null_draws,bootstraps=1999,seed=args.seed)
    assert result["n_original_reference_complete_Direct_species"]==2234
    assert result["n_original_reference_triple_profiles"]==100
    assert result["n_original_928_island_flora"]==928
    assert not result["genus_matched_null_effect_estimated"]
    assert result["source_species_island_taxa_genera_coordinates_not_released"]
    assert not result["measured_actual_pollinator_function_or_adaptation"]
    text=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for forbidden in ('"island_id"','"accepted_species"','"entity_ID"',
                       '"work_species"','"spatial_block"','"genus"','"family"'):
        if forbidden in text:
            raise ValueError("Protected original source metadata leaked to aggregate")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(text+"\n")
    print(json.dumps({
        "status":result["status"],
        "n_original_islands":result["n_original_928_island_flora"],
        "original_profiles":result["n_original_reference_triple_profiles"],
        "species_frequency_bin_cutoffs":{
            "low":result["rarity_bins"]["rare_profile_global_direct_species_frequency_max"],
            "high":result["rarity_bins"]["common_profile_global_direct_species_frequency_min"]},
        "by_region":[{
           "scope":r["scope"],"region":r["regime"],"n_islands":r["n_islands"],
           "n_disjoint_source":r["n_source_disjoint_eligible"],
           "genus_30_support":r["n_source_disjoint_genus_matching_max_30plus"],
           "genus_15_support":r["n_source_disjoint_genus_matching_max_15plus"],
        }for r in result["rare_vs_common_island_source_aggregate"]]
    },sort_keys=True))
if __name__=="__main__":
    main()
