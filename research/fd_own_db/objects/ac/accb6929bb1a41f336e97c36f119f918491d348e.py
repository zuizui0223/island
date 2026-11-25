"""Run the licensed-original fd pollinator opportunity vs floral FD audit.

Reads original frozen sources locally. Writes anonymous JSON only. No private
island, taxon, GIFT region, provider or coordinate rows may be exported.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from fd.pollinator_functional_diversity import audit_original_direct_channel_opportunity

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--progressive-root",required=True,type=Path)
    a.add_argument("--qualified-channel-csv",required=True,type=Path)
    a.add_argument("--output",required=True,type=Path)
    a.add_argument("--rarefactions",type=int,default=8)
    z=a.parse_args()
    root=z.progressive_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"])
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"])
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    channels=pd.read_csv(z.qualified_channel_csv,
        usecols=["island_id","channel_id","channel_state","source_state"])
    sources=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"])
    gift=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        usecols=["entity_ID","work_species"])
    result=audit_original_direct_channel_opportunity(
        ledger,flora,geo,channels,sources,gift,repeats=z.rarefactions)
    assert result["n_original_qualifying_islands"]==928
    assert result["n_direct_triple_species"]>=2000
    assert not result["raw_island_taxon_or_source_identifiers_exported"]
    assert not result["real_local_visitor_shift_measured"]
    assert result["sampling"]["equal_species_per_island_or_source"]==30
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("Rights-safe FD functional and pollinator opportunity aggregates written")
    print(json.dumps({
       "status":result["status"],
       "n_islands":result["n_original_qualifying_islands"],
       "n_source_comparable_by_region":
         {r["regime"]:r["n_source_comparable"] for r in
          result["floral_functional_assemblage_by_region"]
          if r["scope"]=="all_recorded"},
       "n_guild_region_contrasts_admissible":
         result["n_qualified_channel_region_contrasts_passed"],
       "n_Bombus_candidate_cohorts":len(
         result["source_relative_bombus_nondetected_bird_detected_candidates"]),
    },sort_keys=True))
if __name__=="__main__":
    main()
