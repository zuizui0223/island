"""Source-protected recap: original 4505 GIFT regional genus transfer generality.

Consumes original signed anonymous 2026-11-15 source artifact directly.
A new analysis is a synthesis of original data; it is NOT new trait/animal
records, independent botanical phylogeny, colonization or local adaptation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from fd.source_regional_transfer_replication import audit_regional_transfer

EXPECTED_EIGHT={
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","tube_depth_class","inflorescence_display",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
}

def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-anonymous",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    source=json.loads(args.original_anonymous.read_text(encoding="utf-8"))
    r=audit_regional_transfer(source)
    if r["n_original_botanical_global_island_floras"]!=4505:
        raise ValueError("Original GIFT global parent changed")
    if r["n_original_source_comparison_panels"]!=100:
        raise ValueError("The original 100 region×trait×status panels were not replayed")
    if (set(r["all_recorded_replicated_trait_endpoints"])!=EXPECTED_EIGHT
            or r["n_native_only_endpoints_replicating_across_all_three_geographies"]!=0):
        raise ValueError("Original regional source score/interval replication changed")
    local={
        (x["region"],x["endpoint"]):x for x in r[
            "original_same_focal_local_vs_external_source_transfer"]}
    north=local[("northern_midlatitude","floral_form")]
    trop=local[("tropical","floral_form")]
    if (not north["local_training_beats_external_in_both_grouped_intervals"]
        or not trop["external_training_beats_local_in_both_grouped_intervals"]):
        raise ValueError("Original locality-genus conditional floral form sign reversal not confirmed")
    if r["strict_approach"]["genus_conditioning_is_oracle_not_prospective_dispersal_or_colonization_prediction"] is not True:
        raise ValueError("Unsupported prospective source claim")
    body=json.dumps(r,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)
    for term in ("accepted_species","island_id","source_signature",
                 "source_region_id","spatial_block","entity_ID",
                 "source_genus","work_species","latitude","longitude",
                 "reference_id","provider_id"):
        if '"'+term+'"' in body:
            raise ValueError("Private source original GIFT island/taxon/location leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print("ORIGINAL_REGION_GENUS_TRANSFER_REPLICATION="+json.dumps(dict(
        status=r["status"],
        n_original_global_flora_islands=4505,
        n_original_source_comparison_panels=100,
        n_all_recorded_three_region_replicated_endpoints=len(EXPECTED_EIGHT),
        repeated_across_regions=sorted(EXPECTED_EIGHT),
        n_original_native_three_region_replicated_endpoints=0,
        floral_form_north_region_local_minus_external_block_equal=(
            north["local_vs_external_genus_Brier_gain_physical_block_equal"]),
        floral_form_tropical_region_local_minus_external_block_equal=(
            trop["local_vs_external_genus_Brier_gain_physical_block_equal"]),
        no_evolutionary_adaptation_or_effective_pollination_inferred=True,
    ),sort_keys=True))
if __name__=="__main__":
    main()
