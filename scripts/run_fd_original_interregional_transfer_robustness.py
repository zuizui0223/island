"""Archive-replayed original source 4505 island interregional genus trait transfer audit.

Two original source runs are the SAME original GIFT botanical dataset; their
agreement is a reproducibility check, not independent ecological replication.
Inputs may contain only original protected anonymous, signed output, not
taxa/island/genus coordinates or source mainland raw data.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from fd.interregional_transfer_robustness import (
    analyze_original_global_regionality,
)


def run(current:Path,previous:Path)->dict:
    latest=json.loads(current.read_text(encoding="utf-8"))
    older=json.loads(previous.read_text(encoding="utf-8"))
    r=analyze_original_global_regionality(latest,older)
    if r["source_original_four_region_and_missing_denominator"][
        "source_qualified_count_by_status"].get(
         "ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON")!=63:
        raise ValueError("Original 63 qualified source regional panels drifted")
    counts={
        (x["scope"],x["region"]):x for x in r["source_robust_13_endpoint_panel_counts"]
    }
    if (counts[("all_recorded","northern_midlatitude")][
           "n_robust_local_better_both_proper_losses"]<8 or
        counts[("all_recorded","tropical")][
           "n_robust_cross_region_better_both_proper_losses"]<3):
        raise ValueError("Preexisting regional source-transfer result drift")
    reversal={
       (x["scope"],x["trait_endpoint"]):x
       for x in r["original_north_to_tropical_double_metric_sign_reversal"]
    }
    if not all(("all_recorded",ep) in reversal for ep in (
        "floral_form","flower_primary_color__floral_symmetry",
        "self_incompatibility")):
        raise ValueError("Original cross-region source donor directional reversal cannot replicate")
    r["fixed_source_original_artifact_runs"]=dict(
        newer=43700248327,older=43700017352)
    r["original_github_artifact_archives_SHA256"]=dict(
        newer="4b91489394109a75cf600bc1248ca1d1013e10ee525b071a7abc790b95163413",
        older="aedfdc3409072752bd71a98a9f68de0a65f678d836799f1c0492d9efab4dbdb9")
    return r


def main():
    p=argparse.ArgumentParser()
    for arg in ("source-latest","source-older","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    r=run(a.source_latest,a.source_older)
    body=json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)
    for key in (
        "island_id","accepted_species","source_entity_ID","entity_ID",
        "genus","plant","latitude","longitude","source_taxon","exact_location",
        "original_source_island_id","source_region_id","spatial_block",
        "source_species_names","pollinator_species"):
        if '"'+key+'"' in body:
            raise ValueError("Source private original plant/island genus data leak")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    q=r["source_robust_13_endpoint_panel_counts"]
    for x in q:
        print("ORIGINAL_REGION_TRAIT_TRANSFER="+json.dumps(x,sort_keys=True))
    for x in r["original_north_to_tropical_double_metric_sign_reversal"]:
        print("ORIGINAL_CROSS_REGION_REVERSAL="+json.dumps(x,sort_keys=True))
    print("ANONYMOUS_REGIONAL_TRANSFER_VERIFY="+json.dumps(dict(
        status=r["status"],original_GIFT_floras=4505,
        original_13_Direct_trait_endpoints=13,
        n_original_qualified_regional_panels=63,
        n_original_matched_older_source_cells=r[
            "independent_execution_same_original_data_replay"][
                "n_original_matched_cells"],
        no_evolution_or_pollinator_causal_claim=True,
        source_genus_oracle_conditional=True)))
if __name__=="__main__":
    main()
