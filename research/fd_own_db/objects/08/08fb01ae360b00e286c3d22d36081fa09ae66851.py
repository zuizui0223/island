"""Synthetic original 4505 island interregional source, prior and no-oracle claims."""
import copy
import pytest
from fd.source_regional_transfer_replication import (
    audit_regional_transfer,LOSSES,LOCALITY,GROUPS,REGIONS,
)


def _source():
    q=dict(
        status="ORIGINAL_4505_REGION_EXCLUSIVE_EQUAL_GENUS_DONOR_TRANSFER",
        n_original_botanical_island_floras=4505,
        n_original_accepted_source_trait_endpoints=13,
        original_complete_botanical_region_universe=[
            dict(region="northern_midlatitude",n_original_botanical_floras=2205),
            dict(region="tropical",n_original_botanical_floras=1650),
            dict(region="southern_extratropical",n_original_botanical_floras=500),
            dict(region="northern_high_latitude",n_original_botanical_floras=149),
            dict(region="geography_unresolved",n_original_botanical_floras=1)],
        source_original_n5_endpoint_coverage={
            "all_recorded__flower_primary_color__floral_form":2324,
            "original_native__flower_primary_color__floral_form":327,
        },
        claims={
            "entire_4505_original_botanical_population_preserved":True,
            "target_original_all_accepted_plant_names_removed_from_both_donor_sets":True,
            "names_occurring_in_target_and_other_regions_removed_from_both_source_sets":True,
            "local_and_external_same_original_genus_and_exact_equal_donor_species_count":True,
            "donor_species_selected_by_original_name_hash_not_by_trait_state":True,
        },original_equal_genus_donor_transfer_results=[])
    for i,region in enumerate(REGIONS):
        z=dict(trait_endpoint="floral_form",scope="all_recorded",region=region,
             n_same_focal_islands_matched_carrier_n5=120,
             status="ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON")
        for group in GROUPS:
            t=dict(status=z["status"],n_clusters=20)
            for score in (*LOSSES,LOCALITY):
                t[score+"_group_equal_mean"]=.12
                if score==LOCALITY:
                    t[score+"_group_equal_mean"]=.05 if i==0 else -.05
                    t[score+"_group_resample_95"]=(
                        [.03,.07] if i==0 else [-.08,-.02])
                else:
                    t[score+"_group_resample_95"]=[.08,.17]
            z[group]=t
        q["original_equal_genus_donor_transfer_results"].append(z)
    return q


def test_three_original_broad_regions_and_both_proper_losses_replicate():
    original=_source()
    q=audit_regional_transfer(original)
    assert q["n_original_botanical_global_island_floras"]==4505
    assert q["n_all_recorded_endpoints_replicating_across_all_three_geographies"]==1
    assert q["all_recorded_replicated_trait_endpoints"]==["floral_form"]
    assert q["n_native_only_endpoints_replicating_across_all_three_geographies"]==0
    z=[x for x in q["source_exclusive_region_genus_transfer_panels"]
       if x["endpoint"]=="floral_form" and x["scope"]=="all_recorded"]
    assert len(z)==3
    assert all(x["transfer_from_other_regions_beats_genus_blind_prior_on_both_proper_losses_and_both_groupings"] for x in z)
    local=[x for x in q["original_same_focal_local_vs_external_source_transfer"]
           if x["endpoint"]=="floral_form"]
    assert local[0]["local_training_beats_external_in_both_grouped_intervals"]
    assert local[1]["external_training_beats_local_in_both_grouped_intervals"]
    assert q["strict_approach"]["genus_conditioning_is_oracle_not_prospective_dispersal_or_colonization_prediction"]


def test_wrong_geographic_independence_and_brier_interval_fails():
    source=_source()
    source["original_equal_genus_donor_transfer_results"][0][GROUPS[1]][
       LOSSES[0]+"_group_resample_95"]=[-.02,.07]
    q=audit_regional_transfer(source)
    assert q["n_all_recorded_endpoints_replicating_across_all_three_geographies"]==0
    source=_source()
    source["original_equal_genus_donor_transfer_results"][0][GROUPS[1]]["n_clusters"]=9
    q=audit_regional_transfer(source)
    assert q["n_all_recorded_endpoints_replicating_across_all_three_geographies"]==0


def test_holdout_source_origins_never_broadened_and_duplicates_rejected():
    source=_source()
    source["claims"]["target_original_all_accepted_plant_names_removed_from_both_donor_sets"]=False
    with pytest.raises(ValueError,match="excluded-donor"):
        audit_regional_transfer(source)
    source=_source()
    source["original_equal_genus_donor_transfer_results"].append(
        copy.deepcopy(source["original_equal_genus_donor_transfer_results"][0]))
    with pytest.raises(ValueError,match="duplicated"):
        audit_regional_transfer(source)
    source=_source()
    source["n_original_botanical_island_floras"]=928
    with pytest.raises(ValueError,match="4505"):
        audit_regional_transfer(source)
