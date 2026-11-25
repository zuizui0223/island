"""Across-region GENERALIZABILITY audit of original GIFT genus-bound floral states.

Re-score ONLY anonymous original protected source artifact from GIFT 4,505
source-verified floras: strictly distinct plant species names, genus-balanced
local-vs-other-region donors and held-out full focal accepted plant names.

A credible regional source-transfer signal must be positive on BOTH
proper Brier and log losses, and in BOTH 10-degree geographic-block and
first mainland-candidate-source-cluster bootstrap intervals. Only compare
three geographically separate broad nonoverlapping regimes; 13 botanical
trait endpoints are correlated descriptions, NOT 13 independent experiments.
Prediction assumes actual focal botanical genus known and does not model
independent mainland ancestry, within-plant evolution or pollinator service.
"""
from __future__ import annotations
from collections import defaultdict

TRAITS=(
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","tube_depth_class","inflorescence_display",
    "self_incompatibility","autonomous_selfing_capacity","mating_system",
    "cleistogamy","flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
)
REGIONS=(
    "northern_midlatitude","tropical","southern_extratropical",
)
SCOPES=("all_recorded","original_native")
GROUPS=("physical_10degree_block","first_GIFT_mainland_source")
LOSSES=(
    "external_genus_vs_blind_proper_Brier_gain_per_comparable_plant",
    "external_genus_vs_blind_proper_logloss_gain_per_comparable_plant",
)
LOCALITY="local_minus_external_proper_Brier_gain_per_comparable_plant"
MIN_ISLANDS=30
MIN_BLOCKS=10
MIN_SOURCES=10


def _stats(v:dict,key:str)->dict:
    """Numerical keys must be independently sampled original group summaries."""
    if not isinstance(v,dict) or v.get("status")!=(
            "ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON"):
        return dict(eligible=False,mean=None,lo=None,hi=None)
    n=v.get("n_clusters",0)
    mean=v.get(key+"_group_equal_mean")
    ci=v.get(key+"_group_resample_95")
    if (not isinstance(n,int) or n<1 or
        not isinstance(ci,list) or len(ci)!=2 or
        any(q is None or not isinstance(q,(int,float)) for q in ci) or
        mean is None or not isinstance(mean,(int,float))):
        return dict(eligible=False,mean=None,lo=None,hi=None)
    import math
    if not all(map(math.isfinite,(mean,ci[0],ci[1]))) or ci[0]>ci[1]:
        raise ValueError("Original group bootstrap uncertainty not valid")
    return dict(eligible=True,mean=float(mean),
        lo=float(ci[0]),hi=float(ci[1]),n_clusters=int(n))


def _region_panel(row:dict)->dict:
    """Never count positive intervals on under-supported source groups."""
    n=row.get("n_same_focal_islands_matched_carrier_n5")
    if not isinstance(n,(int,float)) or n<MIN_ISLANDS:
        valid=False
    else:
        valid=row.get("status")=="ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON"
    all_scores={}
    valid=valid and all(isinstance(row.get(g),dict) for g in GROUPS)
    if valid:
        valid=(
            row[GROUPS[0]].get("n_clusters",0)>=MIN_BLOCKS and
            row[GROUPS[1]].get("n_clusters",0)>=MIN_SOURCES
        )
    for metric in (*LOSSES,LOCALITY):
        z={g:_stats(row.get(g),metric) for g in GROUPS}
        eligible=bool(valid and all(q["eligible"] for q in z.values()))
        all_scores[metric]=dict(
            jointly_supported=eligible,
            positive_in_all_grouped_intervals=bool(
                eligible and all(q["lo"]>0 for q in z.values())),
            negative_in_all_grouped_intervals=bool(
                eligible and all(q["hi"]<0 for q in z.values())),
            physical_block_equal_mean=z[GROUPS[0]]["mean"],
            physical_block_resample_95=(
                [z[GROUPS[0]]["lo"],z[GROUPS[0]]["hi"]]
                if z[GROUPS[0]]["eligible"] else None),
            first_GIFT_mainland_source_equal_mean=z[GROUPS[1]]["mean"],
            first_GIFT_mainland_source_resample_95=(
                [z[GROUPS[1]]["lo"],z[GROUPS[1]]["hi"]]
                if z[GROUPS[1]]["eligible"] else None))
    both=all(
        all_scores[m]["positive_in_all_grouped_intervals"] for m in LOSSES)
    return dict(
        endpoint=row["trait_endpoint"],
        region=row["region"],
        scope=row["scope"],
        original_n_matched_GIFT_islands=(int(n) if n is not None else None),
        n_original_GIFT_blocks=(
            row.get(GROUPS[0],{}).get("n_clusters")),
        n_original_source_groups=(
            row.get(GROUPS[1],{}).get("n_clusters")),
        original_same_genus_source_comparison_supported=bool(valid),
        transfer_from_other_regions_beats_genus_blind_prior_on_both_proper_losses_and_both_groupings=both,
        source_grouped_scores=all_scores,
    )


def audit_regional_transfer(original:dict)->dict:
    """Source-audited replication across geography, not causal evolution."""
    if (original.get("status")!="ORIGINAL_4505_REGION_EXCLUSIVE_EQUAL_GENUS_DONOR_TRANSFER"
        or original.get("n_original_botanical_island_floras")!=4505
        or original.get("n_original_accepted_source_trait_endpoints")!=13):
        raise ValueError("Not a source-verified all-4505 original GIFT transfer")
    classes=original.get("original_complete_botanical_region_universe")
    if not isinstance(classes,list) or sum(
            x["n_original_botanical_floras"] for x in classes)!=4505:
        raise ValueError("Original full GIFT 4505 regional population not preserved")
    coverage=original.get("source_original_n5_endpoint_coverage")
    if (not isinstance(coverage,dict) or
        coverage.get("all_recorded__flower_primary_color__floral_form")!=2324 or
        coverage.get("original_native__flower_primary_color__floral_form")!=327):
        raise ValueError("Original GIFT trait-wise rather than triple complete-case mask drift")
    claims=original.get("claims",{})
    needed=(
        "entire_4505_original_botanical_population_preserved",
        "target_original_all_accepted_plant_names_removed_from_both_donor_sets",
        "names_occurring_in_target_and_other_regions_removed_from_both_source_sets",
        "local_and_external_same_original_genus_and_exact_equal_donor_species_count",
        "donor_species_selected_by_original_name_hash_not_by_trait_state",
    )
    if any(claims.get(k) is not True for k in needed):
        raise ValueError("Original source excluded-donor and genus-oracle contract broken")
    data=original.get("original_equal_genus_donor_transfer_results")
    if not isinstance(data,list) or not data:
        raise ValueError("Original all-endpoint interregional scores absent")
    seen={}
    for row in data:
        if not isinstance(row,dict):
            raise ValueError("Original regional transfer row malformed")
        k=(row.get("scope"),row.get("trait_endpoint"),row.get("region"))
        if (k in seen or k[0] not in SCOPES or k[1] not in TRAITS or
                k[2] not in (*REGIONS,"northern_high_latitude")):
            raise ValueError("Region/trait/scope source panel duplicated or unknown")
        seen[k]=row
    checks=[]
    panels=[]
    for scope in SCOPES:
        for name in TRAITS:
            cross=[]
            for region in REGIONS:
                original_row=seen.get((scope,name,region))
                if original_row is None:
                    cross.append(dict(endpoint=name,region=region,scope=scope,
                        original_same_genus_source_comparison_supported=False,
                        transfer_from_other_regions_beats_genus_blind_prior_on_both_proper_losses_and_both_groupings=False))
                    continue
                cross.append(_region_panel(original_row))
            panels.extend(cross)
            good=[x for x in cross if x[
                "transfer_from_other_regions_beats_genus_blind_prior_on_both_proper_losses_and_both_groupings"]]
            supported=[x for x in cross if x[
                "original_same_genus_source_comparison_supported"]]
            checks.append(dict(scope=scope,endpoint=name,
                n_required_different_broad_geographic_regions=len(REGIONS),
                n_regions_with_source_qualified_same_genus_name_exclusive_comparison=len(supported),
                n_regions_with_external_genus_better_than_blind_in_all_proper_loss_and_grouped_sensitivities=len(good),
                three_region_external_genus_information_replication=(
                    len(good)==3),
                no_global_conclusion_from_two_of_three_regions=True))
    # Regional locality-specific difference uses SAME matched focal islands
    # and same number of regional/external donors within each genus, with
    # proper Brier response. Never interpret as actual source ancestry.
    locality=[]
    for name in ("floral_form","flower_primary_color__floral_form",
                 "floral_symmetry","tube_depth_class"):
        for region in REGIONS:
            r=seen.get(("all_recorded",name,region))
            if r is None:
                continue
            x=_region_panel(r)
            q=x["source_grouped_scores"][LOCALITY]
            locality.append(dict(
                endpoint=name,region=region,
                n_original_matched_islands=x["original_n_matched_GIFT_islands"],
                n_original_independent_GIFT_blocks=x["n_original_GIFT_blocks"],
                n_original_first_source_groups=x["n_original_source_groups"],
                local_vs_external_genus_Brier_gain_physical_block_equal=q["physical_block_equal_mean"],
                local_vs_external_genus_Brier_gain_physical_block_resample_95=q["physical_block_resample_95"],
                local_vs_external_genus_Brier_gain_first_source_equal=q["first_GIFT_mainland_source_equal_mean"],
                local_vs_external_genus_Brier_gain_first_source_resample_95=q["first_GIFT_mainland_source_resample_95"],
                local_training_beats_external_in_both_grouped_intervals=q[
                    "positive_in_all_grouped_intervals"],
                external_training_beats_local_in_both_grouped_intervals=q[
                    "negative_in_all_grouped_intervals"],
            ))
    cross=[x for x in checks if x["scope"]=="all_recorded"]
    replicated=[x["endpoint"] for x in cross if x[
        "three_region_external_genus_information_replication"]]
    native=[x for x in checks if x["scope"]=="original_native"]
    return dict(
        status="ORIGINAL_4505_GIFT_THREE_REGION_TWO_LOSS_TWO_GROUP_GENUS_TRANSFER_REPLICATION",
        n_original_botanical_global_island_floras=4505,
        n_original_direct_endpoints=13,
        n_original_source_comparison_panels=len(data),
        n_distinct_test_regions=len(REGIONS),
        source_exclusive_region_genus_transfer_panels=panels,
        original_13_endpoints_cross_region_evidence=checks,
        n_all_recorded_endpoints_replicating_across_all_three_geographies=len(replicated),
        all_recorded_replicated_trait_endpoints=replicated,
        n_native_only_endpoints_replicating_across_all_three_geographies=sum(
            x["three_region_external_genus_information_replication"]
            for x in native),
        original_same_focal_local_vs_external_source_transfer=locality,
        strict_approach=dict(
            original_all_4505_source_GIFT_island_floras_retained=True,
            exclude_full_focal_accepted_species_names_from_source_donors=True,
            exclude_shared_region_source_species_names_from_both_training_pools=True,
            condition_prediction_on_known_focal_plant_genus_composition=True,
            equal_donor_species_per_genus_both_sources=True,
            require_at_least_30_matched_islands_10_physical_blocks_and_10_first_source_groups=True,
            both_Brier_and_proper_log_loss_10degree_and_source_cluster_intervals_positive=True,
            original_bootstrap_intervals_exploratory_not_familywise_significance=True,
            original_shared_species_Direct_trait_corpus_not_independent_trait_measurement=True,
            source_training_other_regions_island_floras_not_proven_ancestral_mainland_species=True,
            three_region_traits_correlated_not_independent_evolutionary_replication=True,
            genus_conditioning_is_oracle_not_prospective_dispersal_or_colonization_prediction=True,
            no_independent_phylogenetic_cross_genus_selection_or_pollen_effect_identified=True,
            no_private_taxa_islands_source_entity_or_exact_location_exported=True,
        ))
