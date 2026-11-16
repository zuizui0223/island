"""Protected original GIFT-4505 floral-sector supply × qualified guild support.

All island, taxon, source-region and GBIF IDs remain inside the original
Actions workspace. Export aggregate regional sector/guild coverage only.
This is a *support audit*, not a pollinator effect or historical loss fit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups, source_sectors
from fd.global_partial_pollinator_guild_bridge import (
    CHANNELS, normalize_guild_evidence,
)
from fd.source_guild_support import audit_source_guild_support

# Morphology only; SC/SI and autonomous seed set have separate reproductive scales.
MORPHOLOGY = (
    "zygomorphic", "tube_or_bell", "open_radial", "deep_tube",
    "bilateral_tube_or_bell_joint",
)
ORIGINAL_SHA = {
    "direct_only_trait_ledger.csv.gz":
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
REGIONS = (
    "northern_midlatitude", "tropical", "southern_extratropical",
    "northern_high_latitude", "geography_unresolved",
)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def build_rows(
    ledger: pd.DataFrame, flora: pd.DataFrame, geography: pd.DataFrame,
    assignments: pd.DataFrame, mainland: pd.DataFrame, channels: pd.DataFrame,
) -> pd.DataFrame:
    """Construct every original island×guild×morphology sector, no complete cases."""
    byisland = status_groups(flora)
    if len(byisland) != 4505:
        raise ValueError("original 4,505 botanical floras drifted")
    if geography.island_id.duplicated().any() or geography.island_id.nunique() < 8200:
        raise ValueError("independent GSHHG crosswalk incomplete")
    geo = geography.set_index("island_id").reindex(sorted(byisland))
    if int(geo.analysis_regime.notna().sum()) != 4504:
        raise ValueError("original geographic source reconciliation drifted")
    if assignments.duplicated(["island_id", "source_mode", "source_rank"]).any():
        raise ValueError("duplicate original source ranking")
    a = assignments.loc[assignments.source_mode.eq("geo_k5")].copy()
    a["source_rank"] = pd.to_numeric(a.source_rank, errors="raise")
    a = a.loc[a.source_rank.between(1, 5)]
    nearest = a.sort_values(["island_id", "source_rank"]).groupby(
        "island_id"
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID", "work_species"]).any():
        raise ValueError("original mainland species source memberships duplicated")
    bysource = mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(bysource) != 748:
        raise ValueError("original 748 GIFT mainland source-region universe drifted")
    states = source_sectors(strict_traits(ledger))
    if set(MORPHOLOGY) - set(states):
        raise ValueError("source original floral-sector ontology incomplete")
    guilds = normalize_guild_evidence(channels, set(byisland))
    if len(guilds) != 4505 * len(CHANNELS):
        raise AssertionError("five guilds not available for all 4,505 island IDs")
    guild_lookup = {
        island: group.set_index("channel_id")[["channel_state", "source_state"]].to_dict("index")
        for island, group in guilds.groupby("island_id", sort=False)
    }

    @lru_cache(maxsize=4505)
    def source_roster(five: tuple[str, ...]) -> frozenset[str]:
        return frozenset().union(*(bysource.get(region, set()) for region in five))

    direct_sets = {k: set(states[k]) for k in MORPHOLOGY}
    records = []
    for i in sorted(byisland):
        q = byisland[i]
        five = nearest.get(i, ())
        has_five = len(five) == 5
        candidates = (set(source_roster(five)) - q["full"]) if has_five else set()
        region = str(geo.at[i, "analysis_regime"]) if pd.notna(
            geo.at[i, "analysis_regime"]) else "geography_unresolved"
        if region not in REGIONS:
            raise ValueError("unexpected geography regime")
        block = geo.at[i, "spatial_block"]
        block = str(block) if pd.notna(block) else None
        source_signature = five[0] if has_five else None
        for sector in MORPHOLOGY:
            category = states[sector]
            documented = direct_sets[sector]
            island_names = q["full"] & documented
            source_names = candidates & documented
            n_island = len(island_names)
            n_source = len(source_names)
            island_share = (
                float(sum(category[x] for x in island_names) / n_island)
                if n_island else np.nan
            )
            source_share = (
                float(sum(category[x] for x in source_names) / n_source)
                if n_source else np.nan
            )
            for guild in CHANNELS:
                z = guild_lookup[i][guild]
                records.append(dict(
                    island_id=i, regime=region, guild=guild,
                    channel_state=z["channel_state"],
                    source_state=z["source_state"],
                    spatial_block=block, source_signature=source_signature,
                    sector=sector, source_sector_share=source_share,
                    island_sector_share=island_share,
                    n_direct_carriers=n_island,
                    n_source_direct_carriers=n_source,
                ))
    out = pd.DataFrame.from_records(records)
    if len(out) != 4505 * len(CHANNELS) * len(MORPHOLOGY):
        raise AssertionError("global source-sector and guild grid unexpectedly shrank")
    if out.island_id.nunique() != 4505:
        raise AssertionError("full original botanical flora frame was lost")
    return out


def summarize_source_only(rows: pd.DataFrame) -> list[dict]:
    """Conditional documented island-minus-source proportions, no guild labels."""
    # Use exactly one copy per island/sector; the five guild rows are duplicates.
    core = rows.loc[rows.guild.eq(CHANNELS[0])].copy()
    valid = core.loc[
        core.n_direct_carriers.ge(5)
        & core.n_source_direct_carriers.ge(10)
        & core.island_sector_share.notna()
        & core.source_sector_share.notna()
    ].copy()
    valid["island_minus_candidate_source"] = (
        valid.island_sector_share - valid.source_sector_share
    )
    out = []
    for (regime, sector), frame in valid.groupby(["regime", "sector"], sort=True):
        blocks = frame.loc[frame.spatial_block.notna()].groupby(
            "spatial_block", sort=True
        ).island_minus_candidate_source.mean().to_numpy(float)
        status = ("DESCRIPTIVE_CLUSTER_INTERVAL" if
                  len(frame) >= 30 and len(blocks) >= 10
                  else "INSUFFICIENT_GEOGRAPHIC_SUPPORT")
        row = dict(
            regime=str(regime), sector=str(sector),
            n_original_islands_with_island5_and_source10=int(len(frame)),
            n_geographic_blocks=int(len(blocks)),
            status=status,
            no_guild_conditioning_or_historical_colonization_inference=True,
        )
        if status == "DESCRIPTIVE_CLUSTER_INTERVAL":
            seed = int.from_bytes(hashlib.sha256(
                (str(regime) + "|" + str(sector)).encode()
            ).digest()[:8], "big")
            rng = np.random.default_rng(seed)
            resampled = blocks[
                rng.integers(len(blocks), size=(999, len(blocks)))
            ].mean(axis=1)
            row.update(
                island_minus_candidate_source_fraction_block_mean=float(blocks.mean()),
                original_999_block_resample_95=[
                    float(x) for x in np.quantile(resampled, [.025, .975])
                ],
            )
        out.append(row)
    return out


def run(root: Path, geo_path: Path, channels_path: Path) -> dict:
    ledger_path = root / "fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path = root / "fixed/canonical/input/chapter1_status_flora.csv.gz"
    if digest(ledger_path) != ORIGINAL_SHA[ledger_path.name]:
        raise ValueError("original strictly Direct ledger hash changed")
    if digest(flora_path) != ORIGINAL_SHA[flora_path.name]:
        raise ValueError("original botanical 4,505 source hash changed")
    ledger = pd.read_csv(ledger_path, dtype=str,
                         usecols=["accepted_species", "trait_name",
                                  "normalized_value", "evidence_scope", "evidence_quality"])
    flora = pd.read_csv(flora_path, dtype=str,
                        usecols=["island_id", "accepted_species", "origin_status",
                                 "endemic_status", "floristic_status"])
    geo = pd.read_csv(geo_path, dtype=str,
                      usecols=["island_id", "analysis_regime", "spatial_block"])
    assignments = pd.read_csv(
        root / "fixed/source/result/island_source_assignments.csv.gz", dtype=str,
        usecols=["island_id", "source_mode", "source_rank", "entity_ID"])
    mainland = pd.read_csv(
        root / "fixed/source/source/gift_native_mainland_flora.csv.gz", dtype=str,
        usecols=["entity_ID", "work_species"])
    channel = pd.read_csv(
        channels_path, dtype=str,
        usecols=["island_id", "channel_id", "channel_state", "source_state"])
    rows = build_rows(ledger, flora, geo, assignments, mainland, channel)
    support = audit_source_guild_support(rows)
    if support["n_unique_islands_in_input"] != 4505:
        raise AssertionError("full source-GIFT island frame invalid")
    if not all("mixed_spatial_blocks" in row for row in support["strata"]):
        raise AssertionError("independent within-physical-block overlap guard absent")
    source_only = summarize_source_only(rows)
    return {
        "status": "SOURCE_GUILD_FUNCTIONAL_SECTOR_POSITIVITY_ORIGINAL_4505",
        "original_global_botanical_flora_population": 4505,
        "original_source_mode": "geo_k5",
        "original_mainland_source_region_entities": 748,
        "n_original_guilds": len(CHANNELS),
        "n_real_morphological_sectors": len(MORPHOLOGY),
        "n_original_joint_trait_states_imputed": 0,
        "n_source_island_guild_sector_rows_private": int(len(rows)),
        "guild_support": support,
        "unconditioned_observed_sector_source_comparison": source_only,
        "provenance": {
            "direct_ledger_sha256": digest(ledger_path),
            "botanical_flora_sha256": digest(flora_path),
            "qualified_channel_sha256": digest(channels_path),
        },
        "interpretation": (
            "Qualified GBIF nondetection is not guild extinction. "
            "Source GIFT candidates are not observed seed arrivals. "
            "All statuses include unresolved nativeness. "
            "Morphology is species-wide. No historical pollination shift, "
            "within-species floral evolution or pollen service estimated."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--original-root", type=Path, required=True)
    p.add_argument("--geology-crosswalk", type=Path, required=True)
    p.add_argument("--channel-evidence-csv", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = run(args.original_root, args.geology_crosswalk,
                 args.channel_evidence_csv)
    raw = json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2,
                     allow_nan=False)
    forbidden = (
        "island_id", "accepted_species", "entity_ID", "work_species",
        "source_signature", "spatial_block", "source_region_id",
        "latitude", "longitude", "GBIF_id", "taxon_id",
        "original_photo_id", "source_genus", "source_family",
    )
    if any('"' + x + '"' in raw for x in forbidden):
        raise ValueError("restricted original island/source evidence leaked")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "original_global_botanical_flora_population": 4505,
        "source_guild_sector_strata": result["guild_support"]["n_strata"],
        "n_strata_eligible_for_interaction": sum(
            x["status"] == "ELIGIBLE_FOR_EXPLORATORY_ASSOCIATION"
            for x in result["guild_support"]["strata"]
        ),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
