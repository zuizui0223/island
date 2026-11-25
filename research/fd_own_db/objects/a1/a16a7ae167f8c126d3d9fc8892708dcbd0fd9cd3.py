"""Rights-safe original-source Direct reproductive-portfolio versus floral FD audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from fd.reproductive_assurance_functional_diversity import audit_assurance_flora


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--progressive-root",type=Path,required=True)
    p.add_argument("--qualified-channel-csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    root=a.progressive_root
    ledger=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
       usecols=["accepted_species","trait_name","normalized_value",
                "evidence_scope","evidence_quality"])
    flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
       usecols=["island_id","accepted_species","origin_status"])
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
       usecols=["island_id","analysis_regime","spatial_block"])
    sources=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
       usecols=["island_id","source_mode","source_rank","entity_ID"])
    gift=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
       usecols=["entity_ID","work_species"])
    states=pd.read_csv(a.qualified_channel_csv,
       usecols=["island_id","channel_id","channel_state","source_state"])
    result=audit_assurance_flora(ledger,flora,geo,sources,gift,states)
    assert result["n_original_direct_triple_flora_islands"]==928
    assert result["n_original_focal_Bombus_nondetection_bird_detection_candidates"]==21
    assert not result["no_actual_pollinator_loss_or_shift_measured"] is False
    assert not result["no_species_island_or_source_level_data_exported"] is False
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
       "status":result["status"],
       "n_original_islands":result["n_original_direct_triple_flora_islands"],
       "n_native_sources":result["n_species_with_direct_assurance_by_trait"],
       "n_region_trait_outputs":len(result["reproductive_assurance_functional_assemblage_by_region"]),
       "n_focal_assurance_support_classes":len(result["candidate_potential_pollinator_shift_assurance_support"]),
    },sort_keys=True))
if __name__=="__main__":
    main()
