"""Read only the signed anonymous protected GIFT genus environment/source result.

Actual GIFT 4505 accepted botanical island flora, original four GSHHG
regions and two nested original provenance scopes. No source island, plant,
genus, coordinate, first candidate mainland ID or taxon labels exported.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from fd.environment_genus_source_audit import audit_geo_source_genus_prediction


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    src=json.loads(a.source.read_text(encoding="utf-8"))
    audit=audit_geo_source_genus_prediction(src)
    raw=json.dumps(audit,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)
    for forbidden in ("island_id","accepted_species","entity_ID","source_genus",
                      "genus","genus_name","latitude","longitude",
                      "spatial_block","source_cluster","source_species_names"):
        if '"'+forbidden+'"' in raw:
            raise ValueError("Original GIFT protected botanical/geographic identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    for x in audit["original_full_geographic_composition_cells"]:
        print("ORIGINAL_GIFT_GENUS_ENVIRONMENT="+json.dumps({
            "scope":x["scope"],"region":x["region"],
            "n_source_focal_islands":x["n_original_islands_scored_without_local_genus_oracle"],
            "source_models":{k:v["status"] for k,v in x["original_source_models"].items()},
            "source_paired_increments":{
                k:v["status"] for k,v in x["original_source_paired_predictor_increments"].items()
            },
            "source_similarity_direct_paired_Brier":x["original_source_paired_predictor_increments"][
                "source_similarity_increment_over_environment"][
                    "proper_paired_difference_source_groups"]["Brier__physical_10degree_block"],
        },sort_keys=True))
    print("ANONYMOUS_GIFT_ORIGINAL_GENUS_SOURCE_AUDIT="+json.dumps({
        "status":audit["status"],
        "original_island_parent":audit["n_original_global_botanical_floras"],
        "n_source_four_region_and_two_scope_panels":len(
            audit["original_full_geographic_composition_cells"]),
        "n_paired_source_incremental_sensitivity_families":len(
            audit["original_paired_model_direct_source_increment_summary"]),
        "not_actual_mainland_ancestral_colonization":True,
    },sort_keys=True))
if __name__=="__main__":
    main()
