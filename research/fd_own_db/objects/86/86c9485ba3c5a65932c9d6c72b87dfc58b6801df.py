"""Prospective source×guild interaction support gate; no effect fitting.

All inputs are already-authorized, original protected-source island-level rows.
Only non-identifying aggregate counts and eligibility statuses are returned.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

REQUIRED = {
    "island_id", "regime", "guild", "channel_state", "source_state",
    "spatial_block", "source_signature", "sector", "source_sector_share",
    "n_direct_carriers", "n_source_direct_carriers",
}
RETAINED = "retained"
NONDETECTED = "disrupted"  # Original NEE qualified non-detection, NOT historical loss
UNRESOLVED = "unresolved"


def audit_source_guild_support(rows: pd.DataFrame) -> dict:
    """Return one eligibility decision per regime×guild×sector, without effects.

    A source signature must contain both retained and adequately nondetected
    *islands*; blocks must be distinct within each state. Source-sector
    variation must be nondegenerate within both states. No imputation of
    unresolved channels and no global/regional pooling.
    """
    missing = REQUIRED.difference(rows.columns)
    if missing:
        raise ValueError(f"missing fields: {sorted(missing)}")
    if rows.duplicated(["island_id", "guild", "sector"]).any():
        raise ValueError("duplicate island×guild×sector records")
    allowed = {RETAINED, NONDETECTED, UNRESOLVED, "structurally_absent"}
    if not rows.channel_state.isin(allowed).all():
        raise ValueError("invalid channel_state; do not recode missing as absence")
    if not rows.source_state.isin(["available", "unresolved", "structurally_absent"]).all():
        raise ValueError("invalid source_state")
    if (rows.channel_state.eq(NONDETECTED) & ~rows.source_state.eq("available")).any():
        raise ValueError("qualified non-detection requires actual source availability")
    if (rows.channel_state.eq("structurally_absent") & ~rows.source_state.eq("structurally_absent")).any():
        raise ValueError("source structurally absent state inconsistent")
    if (pd.to_numeric(rows.n_direct_carriers, errors="coerce").dropna() < 0).any():
        raise ValueError("negative carrier counts")
    if (pd.to_numeric(rows.n_source_direct_carriers, errors="coerce").dropna() < 0).any():
        raise ValueError("negative mainland source carrier counts")
    share = pd.to_numeric(rows.source_sector_share, errors="coerce")
    if ((share.dropna() < 0) | (share.dropna() > 1)).any():
        raise ValueError("source sector shares must be probabilities")
    if rows[["island_id", "regime", "guild", "sector"]].isna().any().any():
        raise ValueError("missing island/regime/guild/sector identity")

    out = []
    for (regime, guild, sector), sub in rows.groupby(
        ["regime", "guild", "sector"], sort=True
    ):
        # Only source-positive guilds can be interpreted as potential loss.
        valid = sub.loc[sub.source_state.eq("available")].copy()
        valid["source_sector_share"] = pd.to_numeric(
            valid.source_sector_share, errors="coerce"
        )
        valid["n_direct_carriers"] = pd.to_numeric(
            valid.n_direct_carriers, errors="coerce"
        )
        valid["n_source_direct_carriers"] = pd.to_numeric(
            valid.n_source_direct_carriers, errors="coerce"
        )
        valid = valid.loc[
            valid.n_direct_carriers.ge(5)
            & valid.n_source_direct_carriers.ge(10)
            & valid.source_sector_share.notna()
            & np.isfinite(valid.source_sector_share)
        ]
        a = valid.loc[valid.channel_state.eq(RETAINED)]
        b = valid.loc[valid.channel_state.eq(NONDETECTED)]
        source_states = valid.loc[
            valid.channel_state.isin([RETAINED, NONDETECTED])
        ].groupby("source_signature", dropna=True).channel_state.nunique()
        both_sources = int(source_states.eq(2).sum())
        block_states = valid.loc[
            valid.channel_state.isin([RETAINED, NONDETECTED])
        ].groupby("spatial_block", dropna=True).channel_state.nunique()
        both_blocks = int(block_states.eq(2).sum())
        n_a, n_b = len(a), len(b)
        blocks_a = int(a.spatial_block.nunique())
        blocks_b = int(b.spatial_block.nunique())
        varied_a = a.source_sector_share.nunique() >= 2
        varied_b = b.source_sector_share.nunique() >= 2
        overlap = (
            n_a > 0 and n_b > 0 and
            max(a.source_sector_share.min(), b.source_sector_share.min())
            < min(a.source_sector_share.max(), b.source_sector_share.max())
        )
        reasons = []
        if min(n_a, n_b) < 10:
            reasons.append("INSUFFICIENT_BOTH_STATES")
        if min(blocks_a, blocks_b) < 10:
            reasons.append("INSUFFICIENT_BLOCKS")
        if both_sources < 3:
            reasons.append("INSUFFICIENT_MIXED_SOURCES")
        if both_blocks < 5:
            reasons.append("INSUFFICIENT_WITHIN_BLOCK_OVERLAP")
        if not (varied_a and varied_b and overlap):
            reasons.append("SOURCE_SHARE_NO_COMMON_SUPPORT")
        out.append(dict(
            regime=str(regime), guild=str(guild), sector=str(sector),
            n_source_available=int(sub.source_state.eq("available").sum()),
            n_retained=n_a, n_adequately_nondetected=n_b,
            n_unresolved=int(valid.channel_state.eq(UNRESOLVED).sum()),
            source_pool_condition="GIFT_source_only_distinct_Direct_species_n_ge_10",
            retained_blocks=blocks_a, nondetected_blocks=blocks_b,
            mixed_source_signatures=both_sources,
            mixed_spatial_blocks=both_blocks,
            status="ELIGIBLE_FOR_EXPLORATORY_ASSOCIATION" if not reasons
                   else "NOT_IDENTIFIABLE",
            failed_gates=reasons,
        ))
    return {"n_unique_islands_in_input": int(rows.island_id.nunique()),
            "n_strata": len(out), "strata": out,
            "no_guild_effect_estimated": True}
