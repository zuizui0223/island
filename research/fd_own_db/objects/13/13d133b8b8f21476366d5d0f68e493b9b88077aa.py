import numpy as np
import pandas as pd
import pytest

from fd.pollinator_assemblage_bridge import (
    CHANNELS, adjusted_globi_breadth, channel_support,
)


def fixture(n=60):
    meta = pd.DataFrame({
        "island_id": [f"secret_island_{i}" for i in range(n)],
        "analysis_regime": ["north" if i < n//2 else "tropics" for i in range(n)],
        "spatial_block": [f"block_{i//3}" for i in range(n)],
    })
    rows = []
    for i in range(n):
        for ch in CHANNELS:
            state = ("disrupted" if ch == "bombus" and i % 3 == 0 else
                     "unresolved" if ch == "diptera" and i % 5 == 0 else
                     "retained")
            rows.append({
                "island_id": f"secret_island_{i}",
                "channel_id": ch, "channel_state": state,
                "source_state": "unresolved" if state == "unresolved" else "available",
                "channel_record_count": 1 if state == "retained" else 0,
                "background_record_count": 50 if state == "disrupted" else 10,
            })
    return meta, pd.DataFrame(rows)


def test_frozen_states_are_not_inferred_absences_and_grouped_only():
    meta, poll = fixture()
    result, receipt = channel_support(meta, poll)
    assert receipt["n_islands"] == 60
    assert receipt["n_channel_rows"] == 300
    assert receipt["n_islands_any_confirmed_disrupted"] == 20
    assert receipt["all_unresolved_is_not_absent"]
    bombus = result.loc[result.channel.eq("bombus")]
    assert (bombus.n_disrupted == 10).all()
    assert (bombus.n_retained == 20).all()
    assert not {"island_id", "source_region_id", "spatial_block"} & set(result.columns)


def test_unresolved_cannot_turn_into_disruption():
    meta, poll = fixture()
    mask = poll.channel_state.eq("unresolved")
    poll.loc[mask, "channel_state"] = "disrupted"
    with pytest.raises(ValueError, match="certified without source"):
        channel_support(meta, poll)


def test_missing_channel_is_fail_closed():
    meta, poll = fixture()
    with pytest.raises(ValueError, match="five-channel"):
        channel_support(meta, poll.drop(index=0))


def test_retained_requires_detected_record():
    meta, poll = fixture()
    ix = poll.index[poll.channel_state.eq("retained")][0]
    poll.loc[ix, "channel_record_count"] = 0
    with pytest.raises(ValueError, match="inconsistent"):
        channel_support(meta, poll)


def synthetic_breadth(n=78):
    rng = np.random.default_rng(42)
    refs = rng.normal(3, .4, n)
    cov = {v:rng.normal(0,1,n) for v in (
        "rao_marginal", "log_completeness", "trait_globi_fraction",
        "mean_log_genus_refs", "log_island_area_km2",
        "log_distance_to_continent_km",
        "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"
    )}
    return pd.DataFrame({
        "island_id":[f"private{i}" for i in range(n)],
        "analysis_regime":["north" if i%2==0 else "tropical" for i in range(n)],
        "spatial_block":[f"block{i//3}" for i in range(n)],
        "mean_source_genus_effective_channels": refs,
        "color_arch_excess": .2 + .3*refs + .01*rng.normal(size=n),
        **cov
    })


def test_regression_only_aggregate_and_spatial_cluster_uncertainty():
    z = adjusted_globi_breadth(synthetic_breadth())
    assert z["n_islands"] == 78
    assert z["n_spatial_blocks"] == 26
    assert .25 < z["partial_guild_breadth_slope"] < .35
    assert z["cluster_normal_ci_low"] < z["partial_guild_breadth_slope"]
    assert z["cluster_normal_ci_high"] > z["partial_guild_breadth_slope"]
    assert not {"island_id", "spatial_block", "genus"} & set(z)


def test_repeated_island_and_insufficient_independence_rejected():
    x = synthetic_breadth()
    with pytest.raises(ValueError, match="repeated island"):
        adjusted_globi_breadth(pd.concat([x,x.iloc[[0]]],ignore_index=True))
    x["spatial_block"] = "same_block"
    with pytest.raises(ValueError, match="Insufficient"):
        adjusted_globi_breadth(x)
