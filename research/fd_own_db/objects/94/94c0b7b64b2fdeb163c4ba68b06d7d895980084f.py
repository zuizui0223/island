"""Source original GIFT 4,505 no-focal-genus predictive portability audit.

Takes ONLY the original signed anonymous 13-endpoint output of
zuizui0223/island/scripts/fd_global_original_oracle_free_genus_mix.py.
Both proper categorical Brier and log-loss gains are required in both
original GSHHG physical 10-degree block and first geo_k5 source-group
equal resampling summaries. A local region empirical flower trait
frequency source baseline is checked SEPARATELY from external source
genus-blind category prior. No p-value, ancestral mainland phenotype,
effective pollinator or colonization claim is estimated.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from fd.oracle_free_assembly_evidence import audit_no_oracle_source_flora


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-anonymous",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    origin=json.loads(a.original_anonymous.read_text(encoding="utf-8"))
    result=audit_no_oracle_source_flora(origin)
    body=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False,
                    allow_nan=False)
    for private in (
        "island_id","accepted_species","entity_ID","spatial_block",
        "source_cluster","source_genus","source_name","latitude","longitude",
        "genus_name","source_region_entity","original_taxon_id"):
        if '"'+private+'"' in body:
            raise ValueError("Private protected GIFT floristic/source identifier leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    for row in result["source_region_scope_mode_oracle_free_prediction_counts"]:
        print("ORACLE_FREE_GIFT_SOURCE_PREDICTION="+json.dumps(row,sort_keys=True))
    for test,by_mode in result["three_region_supported_direct_endpoints_by_holdout"].items():
        print("ORACLE_FREE_THREE_REGION_SOURCE="+json.dumps(dict(
            study_scope_and_predictor=test,source_modes=by_mode,
            strict_all_3_mode_intersection=result["all_three_geography_holdout_intersection"][test],
            original_full_GIFT_island_flora_parent=4505),sort_keys=True))
    print("ANONYMOUS_ORACLE_FREE_4505_GIFT_CONCLUSION="+json.dumps(dict(
        source_status=result["status"],
        original_global_botanical_floras=4505,
        n_original_direct_trait_endpoints=13,
        n_original_represented_analytical_support_panels=len(
            result["source_original_13_endpoint_region_scope_training_mode_results"]),
        not_a_prospective_island_genus_colonization_prediction=True,
        source_trait_assignments_are_from_other_GIFT_ISLAND_floras=True,
        no_pollinator_guild_selection_or_ancestral_mainland_evolution_inferred=True),sort_keys=True))
if __name__=="__main__":
    main()
