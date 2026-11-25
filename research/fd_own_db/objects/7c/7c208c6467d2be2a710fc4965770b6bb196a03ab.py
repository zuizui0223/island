"""4,505 original botanical floras are retained when GBIF guild source is absent."""
import numpy as np
import pandas as pd
import pytest

from fd.global_partial_pollinator_guild_bridge import (
    CHANNELS,TARGETS,
    normalize_guild_evidence,
    source_guild_partial_traits,
)
from fd.global_partial_trait_flora_fd import TRAITS,PAIR_TRAITS


def _fixture(n=4001):
    endpoints=list(TRAITS)+["__".join(pair) for pair in PAIR_TRAITS]
    rows=[]
    for i in range(n):
        eligible=i<600
        for scope in ("all_recorded","original_native"):
            for trait in endpoints:
                is_pair="__" in trait
                rows.append(dict(
                    island_id=f"island{i:04}",regime="tropical" if i<1000 else "northern_midlatitude",
                    geological_origin="oceanic" if i<500 else "continental" if i<1000 else "unresolved",
                    original_scope=scope,outcome=trait,
                    outcome_type="pair" if is_pair else "single",
                    spatial_block=f"block{i%29}",
                    source_signature=f"GIFTsource{i%11}",
                    n_original_flora_species=20,
                    n_recorded_direct_species=8 if eligible else 0,
                    observed_direct_coverage=0.4 if eligible else 0.0,
                    observed_finite5=.32+(i%7)*.01 if eligible else np.nan,
                    observed_genus_equal_finite5=(
                        .34+(i%7)*.01 if eligible and not is_pair else np.nan),
                    observed_two_trait_pairing=(
                        .001*(i%7) if eligible and is_pair else np.nan),
                    source_relative_finite5=(
                        -.05+.01*(i%7) if eligible else np.nan),
                    source_relative_genus_equal_finite5=(
                        -.04+.01*(i%7) if eligible and not is_pair else np.nan),
                    source_relative_pairing=(
                        -.003+.001*(i%7) if eligible and is_pair else np.nan),
                ))
    states=[]
    for i in range(600):
        state="retained" if i%2 else "disrupted"
        for guild in CHANNELS:
            states.append(dict(
                island_id=f"island{i:04}",channel_id=guild,
                channel_state=state,source_state="available"))
    return pd.DataFrame(rows),pd.DataFrame(states)


def test_missing_gbif_rows_remain_unknown_not_pollinator_absence():
    private,channels=_fixture()
    out=source_guild_partial_traits(private,channels,seed=10,draws=999)
    assert out["n_original_global_botanical_flora_islands"]==4001
    assert out["n_original_botanical_islands_without_ANY_guild_source_row"]==3401
    assert out["n_original_plant_scope_trait_endpoints"]==4001*2*13
    assert out["n_original_plant_scope_trait_GBIF_channel_rows"]==4001*2*13*5
    assert out["n_exploratory_GBIF_channel_vs_qualified_nondetection_endpoint_contrasts"]>0
    p=next(x for x in out["regional_and_geological_pollinator_guild_opportunity"]
           if x["region"]=="GLOBAL_ALL" and x["geological_origin"]=="all" and
           x["original_flora_status"]=="all_recorded" and
           x["trait_or_pair"]=="flower_primary_color" and
           x["pollinator_guild"]=="bombus")
    assert p["n_ALL_islands_in_this_original_scope"]==4001
    assert p["n_islands_GBIF_guild_row_MISSING_not_absent"]==3401
    assert p["n_guild_qualified_adequate_NONdetection"]==300
    assert p["n_guild_GBIF_observation_retained"]==300
    assert len(p["exploratory_guild_floral_trait_associations"])>=1
    assert out["inference_limits"]["source_disrupted_means_ADEQUATE_nondetection_not_historical_loss"]


def test_gbif_qualified_nondetection_without_available_source_fails():
    _,channels=_fixture()
    channels.loc[channels.channel_state.eq("disrupted"),"source_state"]="unresolved"
    with pytest.raises(ValueError,match="adequate non-detection"):
        normalize_guild_evidence(channels,{"island0000","island0001"})


def test_missing_source_channel_only_is_canonical_unresolved():
    x=pd.DataFrame([dict(
        island_id="A",channel_id="bombus",channel_state="retained",
        source_state="available")])
    z=normalize_guild_evidence(x,{"A","B"})
    assert len(z)==10
    assert z.loc[z.island_id.eq("B"),"channel_state"].eq("unresolved").all()
    assert not z.loc[z.island_id.eq("B"),"source_channel_row_present"].any()
    assert z.loc[z.island_id.eq("A")&z.channel_id.ne("bombus"),
                 "channel_state"].eq("unresolved").all()


def test_gbif_repeated_or_unverified_source_fails():
    x=pd.DataFrame([
        dict(island_id="A",channel_id="bombus",
             channel_state="retained",source_state="available")
    ])
    with pytest.raises(ValueError):
        normalize_guild_evidence(pd.concat([x,x]),{"A"})
    x.loc[0,"channel_state"]="structurally_absent"
    with pytest.raises(ValueError):
        normalize_guild_evidence(x,{"A"})
