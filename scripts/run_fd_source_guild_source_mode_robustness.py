"""Original source-protected same-island floral structure: geo5 vs climate10.

Three *fixed* plant member scopes (all, known-native, known-native nonendemic),
five genuine Direct morphological sectors, original 4,505 GIFT island frame.
Both source pools are potential GIFT regional checklists, not seed arrivals.
Exports only anonymous region/sector/coverage/guild aggregates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    source_sectors, status_groups,
)
from fd.source_mode_robustness import (
    MORPHOLOGY, SCOPES, REGIONS, _valid, same_island_source_mode_summary,
)
from fd.global_partial_pollinator_guild_bridge import (
    CHANNELS, normalize_guild_evidence,
)
from fd.source_guild_support import audit_source_guild_support
from fd.global_guild_source_same_cohort_audit import matched_block_delta

ORIGINAL_SHA = {
    "direct_only_trait_ledger.csv.gz":
    "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
    "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODES = {"geo_k5": 5, "geo50_climate10": 10}


def _sha(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


def build_mode_frame(
    ledger: pd.DataFrame,
    flora: pd.DataFrame,
    geography: pd.DataFrame,
    assignments: pd.DataFrame,
    mainland: pd.DataFrame,
) -> pd.DataFrame:
    frame = status_groups(flora)
    if len(frame) != 4505:
        raise ValueError("original GIFT-4505 flora source scope changed")
    if geography.island_id.duplicated().any() or geography.island_id.nunique() < 8200:
        raise ValueError("original GSHHG geometry or block source invalid")
    geo = geography.set_index("island_id").reindex(sorted(frame))
    if int(geo.analysis_regime.notna().sum()) != 4504:
        raise ValueError("original geographic reconciliation changed")
    src = assignments.loc[assignments.source_mode.isin(MODES)].copy()
    if src.duplicated(["island_id", "source_mode", "source_rank"]).any():
        raise ValueError("original native GIFT source ranks duplicated")
    src.source_rank = pd.to_numeric(src.source_rank, errors="raise").astype(int)
    for mode, n_rank in MODES.items():
        found = src.loc[src.source_mode.eq(mode)]
        if not found.source_rank.between(1, n_rank).all():
            raise ValueError("original source rank range changed")
    if {m: int(src.source_mode.eq(m).sum()) for m in MODES} != {
        "geo_k5": 41325, "geo50_climate10": 82190
    }:
        raise ValueError("source-ranked original GIFT geography opportunity drift")
    bysource = src.sort_values("source_rank").groupby(
        ["island_id", "source_mode"], sort=True
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID", "work_species"]).any():
        raise ValueError("original mainland checklist region×taxon duplicate")
    members = mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(members) != 748:
        raise ValueError("original 748 GIFT mainland source entities changed")

    # Small cache only: a large cache of huge candidate species pools would
    # needlessly hold gigabytes of protected taxon names in CI memory.
    @lru_cache(maxsize=30)
    def source_members(entities: tuple[str, ...]) -> frozenset[str]:
        return frozenset().union(*(members.get(e, set()) for e in entities))

    scored = source_sectors(strict_traits(ledger))
    state = {sector: scored[sector] for sector in MORPHOLOGY}
    documented = {sector: set(scored[sector]) for sector in MORPHOLOGY}
    records = []
    n_climate = 0
    for island in sorted(frame):
        item = frame[island]
        sources = {}
        for mode, rank_count in MODES.items():
            ids = bysource.get((island, mode), ())
            if ids and len(ids) != rank_count:
                raise ValueError("incomplete original ranked mainland source list")
            if ids:
                sources[mode] = set(source_members(ids)) - item["full"]
        if "geo50_climate10" in sources:
            n_climate += 1
        geo5 = bysource.get((island, "geo_k5"), ())
        regime = geo.at[island, "analysis_regime"]
        regime = str(regime) if pd.notna(regime) else "geography_unresolved"
        if regime not in REGIONS:
            raise ValueError("original GSHHG regime invalid")
        block = geo.at[island, "spatial_block"]
        block = str(block) if pd.notna(block) else None
        scopes = {
            "all_recorded": item["full"], "strict_native": item["native"],
            "strict_native_nonendemic": item["native_nonendemic"],
        }
        for sector in MORPHOLOGY:
            direct = documented[sector]
            binary = state[sector]
            counts = {}
            for mode in MODES:
                candidates = sources.get(mode, set())
                names = candidates & direct
                num = len(names)
                counts[mode] = (
                    len(candidates), num,
                    float(sum(binary[s] for s in names) / num)
                    if num else np.nan,
                )
            for scope in SCOPES:
                plants = scopes[scope]
                observed = plants & direct
                n_island = len(observed)
                records.append(dict(
                    island_id=island, regime=regime, scope=scope, sector=sector,
                    spatial_block=block,
                    source_signature=geo5[0] if len(geo5) == 5 else None,
                    n_island_total=len(plants), n_island_direct=n_island,
                    island_sector_share=float(
                        sum(binary[s] for s in observed) / n_island
                    ) if n_island else np.nan,
                    n_geo5_total=counts["geo_k5"][0],
                    n_geo5_direct=counts["geo_k5"][1],
                    geo5_sector_share=counts["geo_k5"][2],
                    n_climate10_total=counts["geo50_climate10"][0],
                    n_climate10_direct=counts["geo50_climate10"][1],
                    climate10_sector_share=counts["geo50_climate10"][2],
                    climate10_source_assigned=("geo50_climate10" in sources),
                ))
    if n_climate != 4485:
        raise ValueError("original 4485/4505 climate source assignments drift")
    result = pd.DataFrame.from_records(records)
    if len(result) != 4505 * len(SCOPES) * len(MORPHOLOGY):
        raise ValueError("original phenotype scope not complete")
    return result


def guild_paired_mode_sensitivity(
    fixed: pd.DataFrame, channels: pd.DataFrame,
) -> dict:
    island_guild = normalize_guild_evidence(channels, set(fixed.island_id))
    results = []
    summaries = []
    for scope in SCOPES:
        sub = fixed.loc[fixed.scope.eq(scope)].copy()
        q = sub.merge(island_guild, on="island_id",
                      how="left", validate="many_to_many")
        if len(q) != 4505 * len(MORPHOLOGY) * len(CHANNELS):
            raise ValueError("missing original island×guild×scope coverage")
        # Preserve all 4505 islands. Only genuinely paired Direct trait/source
        # supported rows qualify; unsupported cells become *missing*, not zero.
        qualify = _valid(q)
        q["source_sector_share"] = q.geo5_sector_share.where(qualify, np.nan)
        q["n_direct_carriers"] = q.n_island_direct
        q["n_source_direct_carriers"] = q[
            ["n_geo5_direct", "n_climate10_direct"]
        ].min(axis=1)
        q = q.rename(columns={"channel_id": "guild"})
        g = audit_source_guild_support(q)
        if len(g["strata"]) != len(REGIONS) * len(CHANNELS) * len(MORPHOLOGY):
            raise ValueError("original 125 morphology×guild×regime support slots changed")
        for record in g["strata"]:
            record["scope"] = scope
            record["both_mainland_source_definitions_required"] = True
            summaries.append(record)
            if record["status"] != "ELIGIBLE_FOR_EXPLORATORY_ASSOCIATION":
                continue
            f = q.loc[
                q.regime.eq(record["regime"]) &
                q.guild.eq(record["guild"]) &
                q.sector.eq(record["sector"]) &
                q.source_state.eq("available") &
                q.channel_state.isin(["retained", "disrupted"]) &
                qualify
            ].copy()
            if len(f) != record["n_retained"] + record["n_adequately_nondetected"]:
                raise AssertionError("original discovery cohort changed after source mode")
            re = f.loc[f.channel_state.eq("retained"), "climate10_sector_share"]
            ne = f.loc[f.channel_state.eq("disrupted"), "climate10_sector_share"]
            common = (
                re.nunique() >= 2 and ne.nunique() >= 2
                and max(re.min(), ne.min()) < min(re.max(), ne.max())
            )
            out = dict(
                scope=scope, region=record["regime"], guild=record["guild"],
                morphological_sector=record["sector"],
                n_retained=int(len(re)), n_adequately_nondetected=int(len(ne)),
                n_physically_mixed_blocks=int(record["mixed_spatial_blocks"]),
                actual_climate_source_fraction_common_support=bool(common),
                source_mode_comparison_uses_same_exact_focal_islands=True,
            )
            if not common:
                out["status"] = "NOT_IDENTIFIABLE_CLIMATE_SOURCE_SHARE_OVERLAP"
                results.append(out)
                continue
            f["island_minus_geo5_source"] = (
                f.island_sector_share - f.geo5_sector_share
            )
            f["island_minus_climate10_source"] = (
                f.island_sector_share - f.climate10_sector_share
            )
            f["climate10_minus_geo5_source_relative"] = (
                f.island_minus_climate10_source -
                f.island_minus_geo5_source
            )
            out["status"] = "DESCRIPTIVE_WITHIN_BLOCK_NOT_CAUSAL"
            for var in (
                "island_minus_geo5_source", "island_minus_climate10_source",
                "climate10_minus_geo5_source_relative",
            ):
                seed = int.from_bytes(hashlib.sha256(
                    f"{scope}|{record['regime']}|{record['guild']}|"
                    f"{record['sector']}|{var}".encode()
                ).digest()[:8], "big")
                out[var + "_retained_minus_nondetected"] = matched_block_delta(
                    f, var, "spatial_block", bootstrap=1999, seed=seed
                )
            results.append(out)
    return {
        "n_preserved_original_scope_guild_sector_regional_slots": len(summaries),
        "n_all_gates_passed_across_all_plant_scopes": int(sum(
            r["status"] == "ELIGIBLE_FOR_EXPLORATORY_ASSOCIATION"
            for r in summaries
        )),
        "n_real_both_source_guild_contrasts_with_climate_overlap": int(sum(
            r["status"] == "DESCRIPTIVE_WITHIN_BLOCK_NOT_CAUSAL"
            for r in results
        )),
        "gates": summaries,
        "source_mode_guild_robustness_comparisons": results,
        "no_measured_pollinator_decline_or_pollen_service": True,
        "different_plant_status_scopes_are_nested_not_independent": True,
    }


def run(root: Path, geo_path: Path, channels: Path) -> dict:
    ledger_path = root / "fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path = root / "fixed/canonical/input/chapter1_status_flora.csv.gz"
    for f in (ledger_path, flora_path):
        if _sha(f) != ORIGINAL_SHA[f.name]:
            raise ValueError("original botanical source SHA drift")
    ledger = pd.read_csv(ledger_path, dtype=str,
                         usecols=["accepted_species", "trait_name",
                                  "normalized_value", "evidence_scope",
                                  "evidence_quality"])
    flora = pd.read_csv(flora_path, dtype=str,
                        usecols=["island_id", "accepted_species", "origin_status",
                                 "endemic_status", "floristic_status"])
    geo = pd.read_csv(geo_path, dtype=str,
                      usecols=["island_id", "analysis_regime", "spatial_block"])
    assign = pd.read_csv(
        root / "fixed/source/result/island_source_assignments.csv.gz",
        dtype=str, usecols=["island_id", "source_mode", "source_rank", "entity_ID"])
    mainland = pd.read_csv(
        root / "fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str, usecols=["entity_ID", "work_species"])
    channel = pd.read_csv(channels, dtype=str, usecols=[
        "island_id", "channel_id", "channel_state", "source_state"])
    private = build_mode_frame(ledger, flora, geo, assign, mainland)
    paired = same_island_source_mode_summary(
        private, expected_islands=4505, draws=999
    )
    guilds = guild_paired_mode_sensitivity(private, channel)
    if paired["n_original_islands_with_climate10_source_assignment"] != 4485:
        raise AssertionError("original source climate10 roster changed")
    return {
        "status": "ORIGINAL_4505_SAME_ISLAND_GEO5_VS_CLIMATE10_FLORAL_SECTORS",
        "n_original_botanical_islands": 4505,
        "n_original_mainland_candidate_source_entities": 748,
        "n_geographic5_source_assignment_rows": 41325,
        "n_climate10_source_assignment_rows": 82190,
        "original_direct_trait_ledger_SHA256": _sha(ledger_path),
        "original_GIFT_flora_SHA256": _sha(flora_path),
        "independent_guild_source_SHA256": _sha(channels),
        "same_island_source_definitions": paired,
        "qualified_guild_robustness_on_exact_paired_sources": guilds,
        "no_original_sensitive_island_species_or_source_IDs_exported": True,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original-root", required=True, type=Path)
    ap.add_argument("--geology-crosswalk", required=True, type=Path)
    ap.add_argument("--channel-evidence-csv", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    result = run(args.original_root, args.geology_crosswalk,
                 args.channel_evidence_csv)
    raw = json.dumps(result, ensure_ascii=False, indent=2,
                     sort_keys=True, allow_nan=False)
    protected = (
        "island_id", "accepted_species", "source_signature", "spatial_block",
        "entity_ID", "work_species", "source_region_id", "latitude",
        "longitude", "taxon_id", "GBIF_id", "source_genus", "source_family",
    )
    if any('"' + q + '"' in raw for q in protected):
        raise ValueError("original sensitive source fields leaked")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "original_islands": result["n_original_botanical_islands"],
        "source_paired_panels": result["same_island_source_definitions"][
            "n_scope_sector_regional_panels"],
        "eligible_guild_slots": result[
            "qualified_guild_robustness_on_exact_paired_sources"
        ]["n_all_gates_passed_across_all_plant_scopes"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
