"""Audit external-region SAME-genus versus genus-blind source transfer after
genuinely deleting original co-block/co-first-source whole-island floras FROM
TRAINING, not merely resampling errors.

Inputs: only source-original anonymous 4,505 GIFT island trait summaries.
Both proper Brier score and log-loss must support positive genus-vs-blind
gain in BOTH geographically independent GSHHG block and first candidate
GIFT source-group resampling summaries. Confidence bounds are exploratory,
not family-wise adjusted or phylogenetic independent replications.

The input conditional model is GIVEN focal island genus composition; training
plants are other ORIGINAL GIFT island floras, NOT separately phenotyped
ancestral mainland floral source, and no pollinator/pollen fitness is observed.
"""
from __future__ import annotations
from collections import Counter
import math
from pathlib import Path
import json
import argparse

MODES=("block_out","first_source_out","block_and_first_source_out")
REGIONS=("northern_midlatitude","tropical",
         "southern_extratropical","northern_high_latitude")
MAIN_REGIONS=REGIONS[:3]
SCOPES=("all_recorded","original_native")
FIELDS=(
  "strict_other_region_genus_vs_blind_Brier",
  "strict_other_region_genus_vs_blind_logloss",
)
PREVIOUS=(
  "baseline_on_identical_carriers_other_region_genus_vs_blind_Brier",
  "baseline_on_identical_carriers_other_region_genus_vs_blind_logloss",
)
GROUPS=("physical_10degree_block","first_candidate_source_group")
QUALIFIED="ORIGINAL_STRICT_TRAINING_EXCLUSION_BOTH_GROUPS_QUALIFIED"
GROUP_QUALIFIED="ORIGINAL_STRICT_TRAINING_EXCLUSION_PAIRED_FOCAL_NAME_HOLDOUT"


def _finite(v):
    return type(v) in (int,float) and math.isfinite(v)


def _metric(group:dict|None,key:str)->dict|None:
    if not isinstance(group,dict) or group.get("status")!=GROUP_QUALIFIED:
        return None
    mid=group.get(key+"_group_equal_mean")
    ci=group.get(key+"_group_resample_95")
    clusters=group.get("n_distinct_original_source_clusters")
    if (not _finite(mid) or not isinstance(ci,list) or len(ci)!=2
        or not all(_finite(x) for x in ci) or ci[0]>ci[1]
        or not _finite(clusters) or clusters<10):
        return None
    return dict(group_equal_mean=float(mid),lower_95=float(ci[0]),
                upper_95=float(ci[1]),n_original_clusters=int(clusters))


def _signed_group_pair(row:dict,key:str)->dict:
    groups={g:_metric(row.get(g),key) for g in GROUPS}
    valid=all(groups.values())
    if valid and all(v["lower_95"]>0 for v in groups.values()):
        status="external_genus_better_than_blind"
    elif valid and all(v["upper_95"]<0 for v in groups.values()):
        status="external_genus_worse_than_blind"
    else:status="insufficient_or_ambiguous"
    return dict(direction=status,physical_10degree_block=groups[GROUPS[0]],
                first_candidate_source_group=groups[GROUPS[1]])


def _validate_source(x:dict)->dict:
    if (x.get("status")!=
        "ORIGINAL_GIFT_4505_GLOBAL_GEO_SOURCE_TRAINING_EXCLUDED_REGION_TRAIT_TRANSFER"
        or x.get("n_original_botanical_GIFT_island_floras")!=4505
        or x.get("n_original_floral_direct_endpoint_families")!=13
        or x.get("n_actual_trait_endpoint_families_tested")!=13
        or tuple(x.get("original_source_training_exclusion_modes",()))!=MODES):
        raise ValueError("Original 4,505 global botanical source or 13 Direct traits changed")
    den=x.get("original_source_training_eligibility_and_scoring")
    scored=x.get("original_source_training_scored_results")
    if not isinstance(den,list) or len(den)!=312 or not isinstance(scored,list):
        raise ValueError("Original full source 312-panel support census missing")
    keys=set()
    for row in den:
        key=(row.get("original_scope"),row.get("original_trait_endpoint"),
             row.get("original_source_region"),
             row.get("source_training_exclusion_mode"))
        if key in keys:
            raise ValueError("Duplicated original flora trait/status/region/mode")
        keys.add(key)
    src=x.get("scientific_limitations",{})
    for flag in (
        "focal_accepted_full_flora_names_excluded_from_all_original_source_donor_profiles",
        "SAME_physical_block_OTHER_island_species_names_removed_from_local_and_external_training",
        "SAME_original_first_mainland_source_OTHER_island_species_names_removed_from_local_and_external_training",
        "rescore_baseline_on_IDENTICAL_strict_source_supported_accepted_species",
        "original_training_species_borrowed_from_GIFT_island_floras_NOT_independently_phenotyped_MAINLAND",
        "original_plant_genus_composition_given_as_ORACLE_not_predicted",
    ):
        if src.get(flag) is not True:
            raise ValueError("Original source rigorous training exclusion contract missing "+flag)
    by_mode=Counter(r["source_training_exclusion_mode"] for r in den)
    if by_mode!={q:104 for q in MODES}:
        raise ValueError("Source population has unbalanced 13x4x2 inclusion modes")
    qualified=sum(r.get("status")==QUALIFIED for r in den)
    if qualified!=178:
        raise ValueError("Original 178 source-qualified panels changed")
    names={r["original_trait_endpoint"] for r in den}
    if len(names)!=13:
        raise ValueError("Original 13 Direct single and true-pair endpoints changed")
    regions=set(r["original_source_region"] for r in den)
    if regions!=set(REGIONS):
        raise ValueError("Four original GSHHG regions changed")
    if set(r["original_scope"] for r in den)!=set(SCOPES):
        raise ValueError("Original direct status scopes changed")
    if set(keys)!={(scope,trait,reg,mode) for scope in SCOPES
           for trait in names for reg in REGIONS for mode in MODES}:
        raise ValueError("Original full 312 source support key identity missing")
    seen=set()
    for row in scored:
        key=(row["original_scope"],row["original_trait_endpoint"],
             row["original_source_region"],row["source_training_exclusion_mode"])
        if key not in keys or key in seen:
            raise ValueError("Source scored subset not unique original support")
        seen.add(key)
        if row.get("status")==QUALIFIED and (
            not _finite(row.get("n_original_original_GIFT_islands_with_strict_common_focal_scores"))
            or row["n_original_original_GIFT_islands_with_strict_common_focal_scores"]<30):
            raise ValueError("Qualified original source panel lacks n>=30 focal islands")
    return dict(
        original_full_parent_island_floras=4505,
        original_direct_endpoints=13,
        original_region_status_trait_mode_eligibility_panels=len(den),
        n_original_source_qualified_region_status_modes=qualified,
        n_original_regression_scored_region_status_modes=len(scored),
        original_names_heldout_source_group_keys=len(keys),
        source_training_eligibility_by_mode=dict(by_mode),
        exact_matched_target_species_after_training_deletion=True,
    )


def audit_genus_blind_dual_scores(x:dict)->dict:
    prov=_validate_source(x)
    results=[]
    scored=x["original_source_training_scored_results"]
    d={}
    for row in scored:
        key=(row["original_scope"],row["original_trait_endpoint"],
             row["original_source_region"],row["source_training_exclusion_mode"])
        both={key_:_signed_group_pair(row,key_) for key_ in FIELDS}
        before={key_:_signed_group_pair(row,key_) for key_ in PREVIOUS}
        b=both[FIELDS[0]]["direction"]
        l=both[FIELDS[1]]["direction"]
        if b==l and b in (
            "external_genus_better_than_blind",
            "external_genus_worse_than_blind"):
            dual=b
        else:dual="not_robust_across_two_proper_losses_and_two_geographic_cluster_methods"
        if row.get("status")!=QUALIFIED and dual!="not_robust_across_two_proper_losses_and_two_geographic_cluster_methods":
            raise AssertionError("Original group-insufficient source panel cannot be robust")
        v=dict(
            original_scope=key[0],
            original_direct_trait_endpoint=key[1],
            original_GSHHG_region=key[2],
            training_group_exclusion=key[3],
            original_source_support_status=row.get("status"),
            n_original_matched_direct_island_floras=int(row.get(
                "n_original_original_GIFT_islands_with_strict_common_focal_scores",0)),
            strict_other_region_genus_vs_blind_dual_loss_direction=dual,
            strict_other_region_genus_vs_blind_Brier=both[FIELDS[0]],
            strict_other_region_genus_vs_blind_logloss=both[FIELDS[1]],
            matched_original_baseline_genus_vs_blind_Brier=before[PREVIOUS[0]],
            matched_original_baseline_genus_vs_blind_logloss=before[PREVIOUS[1]],
            source_only_same_true_direct_focal_species_compared_before_and_after=True,
        )
        results.append(v)
        d[key]=v
    summary=[]
    for mode in MODES:
        for scope in SCOPES:
            for region in REGIONS:
                rows=[r for r in results if
                    r["training_group_exclusion"]==mode and
                    r["original_scope"]==scope and
                    r["original_GSHHG_region"]==region]
                summary.append(dict(mode=mode,scope=scope,region=region,
                    n_of_13_scored_original_source_endpoints=len(rows),
                    n_original_source_qualified_direct_endpoints=sum(
                        r["original_source_support_status"]==QUALIFIED for r in rows),
                    n_both_proper_scores_and_both_groupings_positive=sum(
                        r["strict_other_region_genus_vs_blind_dual_loss_direction"]==
                        "external_genus_better_than_blind" for r in rows),
                    n_both_proper_scores_and_both_groupings_negative=sum(
                        r["strict_other_region_genus_vs_blind_dual_loss_direction"]==
                        "external_genus_worse_than_blind" for r in rows),
                    n_scored_without_robust_dual_source_direction=sum(
                        r["strict_other_region_genus_vs_blind_dual_loss_direction"]==
                        "not_robust_across_two_proper_losses_and_two_geographic_cluster_methods"
                        for r in rows),
                    not_independent_endpoint_evolution_replicates=True))
    triple=[]
    traits=sorted({row["original_direct_trait_endpoint"] for row in results})
    for mode in MODES:
        for scope in SCOPES:
            for endpoint in traits:
                vals=[d.get((scope,endpoint,region,mode)) for region in MAIN_REGIONS]
                if not all(vals):continue
                if all(v["strict_other_region_genus_vs_blind_dual_loss_direction"]==
                       "external_genus_better_than_blind" for v in vals):
                    triple.append(dict(
                        mode=mode,scope=scope,original_direct_endpoint=endpoint,
                        n_independent_broad_geographical_regions=3,
                        source_original_matched_islands_by_region={
                            region:int(v["n_original_matched_direct_island_floras"])
                            for region,v in zip(MAIN_REGIONS,vals)},
                        source_genus_oracle_transfer_better_than_genus_blind_only=True))
    return dict(
        status="ORIGINAL_GIFT_4505_STRICT_TRAINING_EXCLUDED_CROSS_REGION_GENUS_OVER_BLIND_DOUBLE_PROPER_SCORE",
        original_global_botanical_source=prov,
        full_region_status_trait_mode_scored_results=results,
        source_geo_region_status_mode_dual_proper_score_counts=summary,
        original_three_region_strict_double_proper_score_transfers=triple,
        scientific_claim_ceiling=dict(
            source_exclusion_already_refitted_BOTH_genus_and_blind_prior=True,
            proper_Brier_and_logloss_both_required_positive_in_both_original_groupings=True,
            positive_resamples_are_exploratory_not_multiplicity_corrected=True,
            bootstrap_clustering_is_not_independent_genetic_evolution=True,
            actual_focal_island_genus_roster_PROVIDED_as_oracle=True,
            actual_train_source_plant_states_from_other_GIFT_island_floras_NOT_mainland_ancestors=True,
            original_anthropogenic_or_pollinator_mechanism_NOT_measured=True,
            original_4505_botanical_parent_not_restricted_to_matched_source_subset=True,
            no_original_accepted_species_island_genus_candidate_mainland_entity_ID_exposed=True,
        ),
    )


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    x=json.loads(args.source.read_text())
    out=audit_genus_blind_dual_scores(x)
    body=json.dumps(out,indent=2,sort_keys=True,allow_nan=False)
    for forbidden in (
        "accepted_species","island_id","genus","entity_ID",
        "latitude","longitude","spatial_block","source_cluster",
        "original_taxon_name","source_species_names",
    ):
        if '"'+forbidden+'"' in body:
            raise ValueError("Private original GIFT source identifier in output")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    for row in out["source_geo_region_status_mode_dual_proper_score_counts"]:
        print("ORIGINAL_STRICT_GENUS_VS_BLIND="+json.dumps(row,sort_keys=True))
    for row in out["original_three_region_strict_double_proper_score_transfers"]:
        print("ORIGINAL_STRICT_THREE_REGION_GENUS_SIGNAL="+json.dumps(row,sort_keys=True))
    print("ANONYMOUS_STRICT_GLOBAL_GENUS_PORTABILITY="+json.dumps(dict(
        original_parent_botanical_floras=4505,
        original_Direct_endpoints=13,
        source_original_geographical_panels=312,
        n_original_qualified_panels=178,
        n_three_region_strict_dual_proper_score_transfers=len(
            out["original_three_region_strict_double_proper_score_transfers"]),
        no_mainland_phylogenetic_or_pollinator_adaptation_claim=True),
        sort_keys=True))
if __name__=="__main__":
    main()
