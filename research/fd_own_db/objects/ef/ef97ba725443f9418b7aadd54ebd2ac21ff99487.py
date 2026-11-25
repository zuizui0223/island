"""Candidate pollinator-channel turnover audit for Direct floral island cohorts.

This scores only documented island-level guild observation opportunities.
It NEVER calls a local species x guild interaction or an evolutionary
pollinator shift merely because a guild was detected/adequately undetected.
"""
from __future__ import annotations

import pandas as pd

GUILDS = ("bombus", "non_bombus_bees", "lepidoptera",
          "diptera", "flower_visiting_birds")
STATES = {"retained", "disrupted", "unresolved"}


def audit_shift_opportunity(
    focal_islands: pd.DataFrame,
    channel_evidence: pd.DataFrame,
    *, min_trait_complete_species: int = 30,
) -> tuple[pd.DataFrame, dict]:
    """Report only anonymous pair-of-guild support counts and evidence limits.

    `disrupted` means source-available plus adequate local nondetection,
    not extinction or historically observed pollinator loss. `retained` is
    GBIF animal occurrence, NOT the animal visiting any particular flower.
    """
    if min_trait_complete_species < 30:
        raise ValueError("Requires at least 30 Direct complete floral species")
    need = {"island_id", "analysis_regime", "spatial_block",
            "n_direct_complete_species"}
    if not need.issubset(focal_islands):
        raise ValueError("Missing original trait-complete island cohort")
    if not {"island_id", "channel_id", "channel_state", "source_state"}.issubset(channel_evidence):
        raise ValueError("Missing original qualified pollinator-source evidence")
    x = focal_islands[list(need)].copy()
    if (x.isna().any().any() or x.island_id.duplicated().any()
            or (x.n_direct_complete_species < min_trait_complete_species).any()):
        raise ValueError("Invalid island cohort / species count / spatial identities")
    q = channel_evidence[
        ["island_id", "channel_id", "channel_state", "source_state"]
    ].copy()
    if q.isna().any().any() or q.duplicated(["island_id", "channel_id"]).any():
        raise ValueError("Missing or duplicate source-channel evidence")
    if not q.channel_id.isin(GUILDS).all() or not q.channel_state.isin(STATES).all():
        raise ValueError("Unexpected pollinator channel or status")
    if (q.channel_state.eq("disrupted") & q.source_state.ne("available")).any():
        raise ValueError("Adequate nondetection requires established source opportunity")
    if set(q.island_id) != set(x.island_id) or len(q) != len(x) * len(GUILDS):
        raise ValueError("Island/cohort channel evidence does not match exactly")
    wide = q.pivot(index="island_id", columns="channel_id",
                   values="channel_state")
    if wide.isna().any().any() or set(wide.columns) != set(GUILDS):
        raise ValueError("Missing one or more observation channels")
    wide = wide.join(x.set_index("island_id")[["analysis_regime", "spatial_block"]])
    rows = []
    for absent in GUILDS:
        for candidate in GUILDS:
            if absent == candidate:
                continue
            considered = wide[absent].eq("disrupted")
            observed = wide[candidate].eq("retained")
            both = considered & observed
            alternative = [g for g in GUILDS if g not in (absent, candidate)]
            n_other = wide.loc[both, alternative].eq("retained").sum(axis=1)
            rows.append({
                "adequately_undetected_guild": absent,
                "co_recorded_candidate_guild": candidate,
                "n_islands_candidate_co_record": int(both.sum()),
                "n_spatial_blocks": int(wide.loc[both, "spatial_block"].nunique()),
                "n_regions": int(wide.loc[both, "analysis_regime"].nunique()),
                "n_with_at_least_one_other_retained_guild": int((n_other > 0).sum()),
                "n_with_all_other_guilds_retained": int((n_other == len(alternative)).sum()),
                "n_islands_both_guilds_retained": int(
                    (wide[absent].eq("retained") & observed).sum()),
                "n_islands_candidate_guild_adequately_undetected": int(
                    (wide[absent].eq("retained") &
                     wide[candidate].eq("disrupted")).sum()),
                "actual_species_level_pollinator_shift_certified": False,
            })
    group = pd.DataFrame(rows)
    receipt = {
        "n_direct_complete_morphology_islands": int(len(x)),
        "n_exact_island_pollinator_channel_rows": int(len(q)),
        "n_islands_with_at_least_one_adequately_undetected_guild":
            int(wide[list(GUILDS)].eq("disrupted").any(axis=1).sum()),
        "n_islands_with_both_retained_and_adequately_undetected_guilds":
            int((wide[list(GUILDS)].eq("disrupted").any(axis=1) &
                 wide[list(GUILDS)].eq("retained").any(axis=1)).sum()),
        "n_confirmed_same_plant_species_mainland_island_visitor_links_in_input": 0,
        "n_certified_evolutionary_pollinator_shifts": 0,
        "note": "No plant-specific interaction, pollen delivery or colonization ancestor links in these guild occurrence tables. Both zero counts refer to absent input evidence, NOT absence of real shifts.",
    }
    return group, receipt
