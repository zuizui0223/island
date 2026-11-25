"""Synthetic tests: no original protected island or source identities."""
import pandas as pd
import pytest

from fd.source_guild_support import audit_source_guild_support


def sample(n=40, mixed=True):
    data = []
    for i in range(n):
        state = "retained" if i % 2 == 0 else "disrupted"
        data.append(dict(
            island_id=f"synthetic{i}", regime="tropical", guild="bombus",
            channel_state=state, source_state="available",
            spatial_block=f"block{i // 2}",
            source_signature=f"source{(i // 2) % 6}" if mixed else f"source{int(i % 2)}",
            sector="zygomorphic", source_sector_share=.1 + (i // 2) / 100,
            n_direct_carriers=8, n_source_direct_carriers=20,
        ))
    return pd.DataFrame(data)


def test_eligible_when_each_state_has_blocks_and_mixed_sources():
    r = audit_source_guild_support(sample())
    x = r["strata"][0]
    assert x["status"] == "ELIGIBLE_FOR_EXPLORATORY_ASSOCIATION"
    assert x["mixed_source_signatures"] >= 3
    assert r["no_guild_effect_estimated"]


def test_unresolved_is_not_absence():
    df = sample()
    df.loc[df.channel_state.eq("disrupted"), "channel_state"] = "unresolved"
    x = audit_source_guild_support(df)["strata"][0]
    assert x["n_adequately_nondetected"] == 0
    assert x["n_unresolved"] == 20
    assert x["status"] == "NOT_IDENTIFIABLE"


def test_source_signature_separation_fails():
    x = audit_source_guild_support(sample(mixed=False))["strata"][0]
    assert x["mixed_source_signatures"] == 0
    assert "INSUFFICIENT_MIXED_SOURCES" in x["failed_gates"]


def test_tiny_nondetection_group_fails():
    df = sample()
    df.loc[df.index[1:38:2], "channel_state"] = "unresolved"
    x = audit_source_guild_support(df)["strata"][0]
    assert x["n_adequately_nondetected"] == 1
    assert "INSUFFICIENT_BOTH_STATES" in x["failed_gates"]


def test_source_share_no_overlap_fails():
    df = sample()
    df.loc[df.channel_state.eq("retained"), "source_sector_share"] = .1
    df.loc[df.channel_state.eq("disrupted"), "source_sector_share"] = .9
    x = audit_source_guild_support(df)["strata"][0]
    assert "SOURCE_SHARE_NO_COMMON_SUPPORT" in x["failed_gates"]


def test_duplicate_island_guild_sector_rejected():
    df = sample()
    with pytest.raises(ValueError):
        audit_source_guild_support(pd.concat([df, df.iloc[:1]]))


def test_missing_channel_is_not_silently_reclassified():
    df = sample()
    df.loc[0, "channel_state"] = None
    with pytest.raises(ValueError):
        audit_source_guild_support(df)


def test_unqualified_disrupted_state_rejected():
    df = sample()
    df.loc[1, "source_state"] = "unresolved"
    with pytest.raises(ValueError, match="qualified non-detection"):
        audit_source_guild_support(df)


def test_structural_absence_not_mistaken_for_qualified_nondetection():
    df = sample()
    ix = df.channel_state.eq("disrupted")
    df.loc[ix, "channel_state"] = "structurally_absent"
    df.loc[ix, "source_state"] = "structurally_absent"
    x = audit_source_guild_support(df)["strata"][0]
    assert x["n_adequately_nondetected"] == 0
    assert x["status"] == "NOT_IDENTIFIABLE"


def test_mainland_source_documentation_gate():
    df = sample()
    df.loc[df.channel_state.eq("disrupted"), "n_source_direct_carriers"] = 3
    x = audit_source_guild_support(df)["strata"][0]
    assert x["n_adequately_nondetected"] == 0
    assert x["status"] == "NOT_IDENTIFIABLE"


def test_both_states_in_separate_blocks_are_not_an_interaction_support():
    df = sample(n=40)
    df["spatial_block"] = [f"state_{z}_{i}" for i, z in enumerate(df.channel_state)]
    x = audit_source_guild_support(df)["strata"][0]
    assert x["retained_blocks"] == 20
    assert x["nondetected_blocks"] == 20
    assert x["mixed_source_signatures"] >= 3
    assert x["mixed_spatial_blocks"] == 0
    assert x["status"] == "NOT_IDENTIFIABLE"
    assert "INSUFFICIENT_WITHIN_BLOCK_OVERLAP" in x["failed_gates"]
