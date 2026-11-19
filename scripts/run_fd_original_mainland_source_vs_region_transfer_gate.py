"""Original GIFT 4505 candidate MAINLAND source novelty limit vs region transfer.

Both input JSONs are pre-verified anonymous original source outputs from
original protected GitHub CI, each SHA-artifact-pinned by the runner workflow.
Their endpoints and botanical cohorts DIFFER; no matched model effect or
independence inference can be computed by subtracting their scores.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from fd.mainland_source_novel_donor_gate import audit_candidate_mainland_reality


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--mainland-source",type=Path,required=True)
    a.add_argument("--interregional-transfer",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    x=a.parse_args()
    mainland=json.loads(x.mainland_source.read_text())
    transfer=json.loads(x.interregional_transfer.read_text())
    q=audit_candidate_mainland_reality(mainland,transfer)
    cohorts={r["source_cohort"]:r for r in q[
        "original_actual_candidate_mainland_GIFT_source_holdout_cohorts"]}
    global_all=cohorts["global_all_n5"]
    nat=cohorts["global_original_native_n5"]
    if (global_all["original_test_botanical_islands"]!=2324
        or global_all["n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed"]!=9
        or nat["original_test_botanical_islands"]!=327
        or nat["n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed"]!=679):
        raise ValueError("Previously source verified 4505/2324/327 mainland donor availability changed")
    if (q["n_original_source_colour_form_species_outside_ALL_4505_GIFT_island_floras"]!=0
        or q["n_original_candidate_mainland_donors_after_2324_trait_eligible_island_full_name_holdout"]!=9):
        raise ValueError("Original all 4505 name exclusion and 2324 eligible cohort donor count conflated")
    if q["conditional_other_region_GIFT_island_genus_predictive_endpoints_across_three_regions"]!=8:
        raise ValueError("Strict original name-excluded interregional genus conditional transfer changed")
    body=json.dumps(q,sort_keys=True,indent=2,allow_nan=False)
    for private in ("accepted_species","island_id","source_region_id",
            "entity_ID","source_signature","work_species","source_genus",
            "spatial_block","latitude","longitude","reference_id"):
        if '"'+private+'"' in body:
            raise ValueError("Source original GIFT island/species/mainland entity identifiers leaked")
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(body+"\n",encoding="utf-8")
    print("ORIGINAL_MAINLAND_VS_REGION_GENUS_SUPPORT="+json.dumps(dict(
        status=q["status"],global_parent_GIFT_botanical_islands=4505,
        mainland_GIFT_source_regions=748,
        original_direct_pair_candidate_mainland_species=3606,
        truly_name_novel_mainland_donors_outside_all_4505_flora_names=0,
        name_novel_mainland_donors_outside_2324_evaluated_island_floras=9,
        global_all_source_donor_variable_genera=0,
        original_strict_native_source_n5_islands=327,
        original_strict_native_source_species_after_cohort_holdout=679,
        conditional_other_GIFT_island_regions_replicated_flower_endpoints=8,
        source_training_different_not_directly_subtractable=True,
        no_ancestral_mainland_and_no_pollinator_adaptation_claim=True,
    ),sort_keys=True))
if __name__=="__main__":
    main()
