"""Synthetic controls for ecological island-level GloPL vs floral diversity association.

The independent pollen supplementation experiment is never treated as a
plant-specific visit, pollinator abundance/decline, or causal mediation.
"""
import numpy as np
import pandas as pd
import pytest

from fd.glopl_floral_fd_association import audit_crosssectional_pollen_fd


def _frames(n=32):
    rng=np.random.default_rng(1002)
    site=[]
    exp=[]
    for j in range(n):
        region=("northern_midlatitude" if j<14 else
                "tropical" if j<27 else "southern_extratropical")
        area=np.log1p(8+j*5)
        dist=np.log1p(j*11+1)
        q=float(.15*np.sin(j*1.713)+rng.normal(0,.018))
        qi=float(q*.65+rng.normal(0,.014))
        f=float(.03*np.cos(j*2.17)+rng.normal(0,.02))
        island=f"fake-{j:02d}"
        site.append(dict(
            island_id=island,analysis_regime=region,
            spatial_block=f"block{j%22}",
            log_distance_to_continent_km=dist,
            log_island_area_km2=area,
            source_delta_marginal_rao=q,
            source_delta_joint_rao=qi,
            source_delta_open_radial_share=f))
        exp.append(dict(island_id=island,study_key=f"fake-study{j}",
            PL_Effect_Size_mean=.15-1.2*q+.04*rng.normal()))
    return pd.DataFrame(site),pd.DataFrame(exp)


def test_analytic_island_response_is_observational_not_plant_specific():
    a,b=_frames()
    r=audit_crosssectional_pollen_fd(a,b,block_draws=999)
    assert r["n_islands_with_source_Q_and_finite_publication_effects"]==32
    assert r["n_independent_original_publications"]==32
    assert r["n_spatial_blocks"]>=15
    assert r["temporal_decline_or_effective_visitor_identity_tested"] is False
    assert r["ecological_mediation_or_adaptation_identified"] is False
    assert r["original_source_taxa_sites_island_ids_exported"] is False
    assert len(r["endpoints"])==3
    assert r["endpoints"][0]["original_sd_scale_slope"]<0
    assert r["endpoints"][0]["status"]=="POSTHOC_CONDITIONAL_ASSOCIATION_NOT_CAUSAL"


def test_insufficient_blocks_suppress_biological_beta():
    a,b=_frames(32)
    a["spatial_block"]="single-block"
    r=audit_crosssectional_pollen_fd(a,b,block_draws=999)
    assert r["status"]=="ECOLOGICAL_ASSOCIATION_NOT_IDENTIFIABLE"
    assert r["endpoints"]==[]


def test_duplicate_publication_island_rejected():
    a,b=_frames()
    with pytest.raises(ValueError,match="Multiple experiment"):
        audit_crosssectional_pollen_fd(a,pd.concat([b,b.iloc[[0]]],ignore_index=True),
          block_draws=999)
