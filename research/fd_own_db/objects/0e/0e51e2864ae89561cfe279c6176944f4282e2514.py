"""Original full-GIFT flora must not be preselected by source photo/guild effect.

The source-same-cohort and matched-block diagnostic only analyzes eligible
subsets transparently; all 4,505 islands remain in the original frame.
Synthetic inputs test that an association may vanish on the exact same
islands, and that paired geographic source blocks are needed for any
conditional interpretation.
"""
import numpy as np
import pandas as pd

from fd.global_guild_source_same_cohort_audit import (
    ENDPOINTS,PRIMARY_SPECS,
    group_status_contrast,matched_block_delta,
    source_cohort_and_guild_robustness,
)


def _source(n=4001):
    all=[]
    for i in range(n):
        island=f"isl{i:05}"
        for scope in ("all_recorded","original_native"):
            for trait in ENDPOINTS:
                is_pair="__" in trait
                eligible=i<120
                source=i<72
                # exact raw data everywhere in first 120, source-data
                # availability only first 72. The raw source-matched
                # stratum has a different mean than the whole raw frame.
                raw=.45 +(.12 if i%2 else -.12)
                src_rel=(.015 if i%2 else .015)
                all.append(dict(
                    island_id=island,regime="tropical" if i<1000 else "northern_midlatitude",
                    geological_origin="oceanic" if i<120 else "unresolved",
                    original_scope=scope,outcome=trait,
                    spatial_block=f"block{i//2%12}",
                    source_signature=f"source{i//2%8}",
                    n_recorded_direct_species=8 if eligible else 0,
                    observed_direct_coverage=.4 if eligible else 0,
                    observed_finite5=raw if eligible else np.nan,
                    source_relative_finite5=src_rel if source else np.nan,
                    observed_genus_equal_finite5=raw if eligible and not is_pair else np.nan,
                    source_relative_genus_equal_finite5=src_rel if source and not is_pair else np.nan,
                    observed_two_trait_pairing=raw/20 if eligible and is_pair else np.nan,
                    source_relative_pairing=src_rel/20 if source and is_pair else np.nan,
                ))
    states=[]
    for i in range(120):
        for channel in ("bombus","flower_visiting_birds","non_bombus_bees"):
            states.append(dict(island_id=f"isl{i:05}",
                channel_id=channel,
                channel_state="retained" if i%2 else "disrupted",
                source_state="available"))
    return pd.DataFrame(all),pd.DataFrame(states)


def test_raw_and_reference_are_same_physical_islands_with_block_mixing():
    df,states=_source()
    result=source_cohort_and_guild_robustness(df,states,seed=20261114,draws=999)
    assert result["n_original_individual_botanical_islands"]==4001
    assert result["n_source_same_island_panels_with_support"]>0
    p=next(x for x in result["original_all_cohort_and_same_GIFT_source_analyses"]
           if x["region"]=="GLOBAL_ALL"
           and x["original_source_flora_status"]=="all_recorded"
           and x["GBIF_guild"]=="bombus"
           and x["trait_or_pair"]=="flower_primary_color__floral_form"
           and x["original_raw_metric"]=="observed_finite5")
    assert p["full_raw"]["n_GBIF_detected"]==60
    assert p["exact_same_island_source"]["n_GBIF_detected"]==36
    assert p["full_raw"]["n_qualified_nondetection"]==60
    assert p["exact_same_island_source"]["n_qualified_nondetection"]==36
    assert p["physical_block_within_full_raw"]["n_mixed_clusters"]==12
    assert p["physical_block_within_raw"]["n_mixed_clusters"]==12
    assert p["physical_block_within_raw"]["mean"]>0
    assert np.isclose(p["physical_block_within_source_relative"]["mean"],0)
    assert result["n_original_primary_eight_slot_testable_with_mixed_blocks"]==8
    assert all(x["BH_q_with_untestable_slots_conservatively_p1"] is not None
               for x in result["source_fixed_posthoc_primary_eight_slot_family"])
    assert result["source_limits"]["original_global_4505_flora_denominator_not_cut_to_50_islands"]


def test_within_cluster_inadequate_when_states_in_different_clusters():
    x=pd.DataFrame([
        dict(spatial_block=f"b{i}",channel_state="retained",
             observed_finite5=.4) for i in range(8)
    ] + [
        dict(spatial_block=f"d{i}",channel_state="disrupted",
             observed_finite5=.3) for i in range(8)
    ])
    d=matched_block_delta(x,"observed_finite5","spatial_block",
                          bootstrap=999,seed=1)
    assert d["status"]=="MIXED_CLUSTER_POWER_TOO_LOW"
    assert d["n_mixed_clusters"]==0


def test_block_sign_exact_diagnostic_not_causal():
    x=pd.DataFrame([
        dict(spatial_block=f"b{i}",channel_state=s,
             v=1. if s=="retained" else 0.)
        for i in range(5) for s in ("retained","disrupted")
    ])
    z=matched_block_delta(x,"v","spatial_block",bootstrap=999,seed=3)
    assert z["n_mixed_clusters"]==5
    assert z["two_sided_cluster_sign_diagnostic_p"]==2/32
    assert z["sign_method"].startswith("EXACT")
    assert z["status"]=="DESCRIPTIVE_NOT_CAUSAL"
