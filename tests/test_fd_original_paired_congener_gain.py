"""Source-safe synthetic tests for ORIGINAL 4505-population paired aggregator."""
from __future__ import annotations
import numpy as np
import pandas as pd
from run_fd_original_paired_congener_gain import paired_aggregate,SECTORS,MODES,SCOPES,GATES

def _frame(*, n_supported=40):
    rows=[]
    for i in range(4505):
        region="tropical" if i<40 else "northern_midlatitude"
        for scope in SCOPES:
            for gate in GATES:
                for sector in (*SECTORS,"zygomorphic_deep"):
                    eligible=(i<n_supported and sector in SECTORS)
                    row=dict(island_id=f"island{i}",regime=region,
                        scope=scope,donor_gate=gate,sector=sector,
                        spatial_block=f"block{i//2}",source_signature=f"source{i//2}")
                    for mode in MODES:
                        row.update({
                            mode+"__n_target_direct":10 if eligible else 0,
                            mode+"__n_source_matched_direct":20 if eligible else 0,
                            mode+"__brier_pooled":(
                                (.3 if sector==SECTORS[0] else .1)
                                if eligible else np.nan),
                            mode+"__brier_genus":(
                                (.1 if sector==SECTORS[0] else .3)
                                if eligible else np.nan),
                        })
                    rows.append(row)
    return pd.DataFrame(rows)

def test_original_full_population_matched_results():
    z=paired_aggregate(_frame(),draws=199)
    assert z["n_original_botanical_islands"]==4505
    assert z["n_panels"]==2*len(SCOPES)*len(GATES)*len(MODES)
    t=[x for x in z["panels"] if x["region"]=="tropical"]
    assert all(x["n_common_source_supported_islands"]==40 for x in t)
    assert all(x["status"]=="DESCRIPTIVE_SOURCE_ONLY_COMMON_ISLAND_PAIRED" for x in t)
    assert all(abs(x["block_equal_delta_gain_radial_minus_tubular"]-.4)<1e-12 for x in t)
    north=[x for x in z["panels"] if x["region"]=="northern_midlatitude"]
    assert all(x["n_common_source_supported_islands"]==0 for x in north)
    assert all(x["status"]=="NOT_IDENTIFIABLE_COMMON_ISLAND_SUPPORT" for x in north)

def test_no_pair_support_is_not_evidence_of_zero_effect():
    z=paired_aggregate(_frame(n_supported=0),draws=199)
    assert all(x["status"]=="NOT_IDENTIFIABLE_COMMON_ISLAND_SUPPORT" for x in z["panels"])
    assert all("block_equal_delta_gain_radial_minus_tubular" not in x for x in z["panels"])
