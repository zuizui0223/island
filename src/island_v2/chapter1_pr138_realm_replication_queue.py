from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

CORE_SYNDROMES = [
    "bird_like",
    "butterfly_like",
    "generalized_accessible",
    "large_bee_like",
    "selfing_syndrome",
]


def build_realm_replication_queue(
    status_flora: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    species_syndrome: pd.DataFrame,
    observation_support: pd.DataFrame,
    *,
    target_realm: str,
    min_flora_species: int = 20,
    max_flora_species: int = 400,
    min_species_per_syndrome: int = 10,
    min_syndromes_meeting_support: int = 4,
) -> pd.DataFrame:
    """Build an outcome-blind queue for a realm that is short of the replication gate.

    Ranking uses only flora size, syndrome-score availability, geography, and spatial
    support. Syndrome values, distance slopes, p-values, and preliminary effect signs
    are never used.
    """
    required_status = {
        "island_id",
        "accepted_species",
        "origin_status",
        "status_resolved",
    }
    required_realm = {"island_id", "biogeographic_realm"}
    required_species = {"accepted_species", "syndrome"}
    required_observation = {
        "island_id",
        "flora_recorded",
        "n_flora_species_recorded",
        "distance_to_continent_km",
        "area_km2",
        "spatial_block",
    }
    for name, frame, required in [
        ("status_flora", status_flora, required_status),
        ("realm_assignment", realm_assignment, required_realm),
        ("species_syndrome", species_syndrome, required_species),
        ("observation_support", observation_support, required_observation),
    ]:
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{name} missing required columns: {sorted(missing)}")

    realm_islands = set(
        realm_assignment.loc[
            realm_assignment["biogeographic_realm"].eq(target_realm), "island_id"
        ].astype(str)
    )
    if not realm_islands:
        return pd.DataFrame()

    status = status_flora.copy()
    status["island_id"] = status["island_id"].astype(str)
    resolved_native = set(
        status.loc[
            status["island_id"].isin(realm_islands)
            & status["status_resolved"].fillna(False).astype(bool)
            & status["origin_status"].eq("native"),
            "island_id",
        ]
    )

    unresolved = status.loc[
        status["island_id"].isin(realm_islands)
        & ~status["island_id"].isin(resolved_native),
        ["island_id", "accepted_species"],
    ].drop_duplicates()

    syndrome_availability = species_syndrome[
        ["accepted_species", "syndrome"]
    ].drop_duplicates()
    coverage = (
        unresolved.merge(syndrome_availability, on="accepted_species", how="inner")
        .groupby(["island_id", "syndrome"])["accepted_species"]
        .nunique()
        .unstack(fill_value=0)
    )
    for syndrome in CORE_SYNDROMES:
        if syndrome not in coverage.columns:
            coverage[syndrome] = 0
    coverage = coverage[CORE_SYNDROMES].copy()
    coverage["min_syndrome_species"] = coverage.min(axis=1)
    coverage["n_syndromes_meeting_support"] = (
        coverage[CORE_SYNDROMES] >= min_species_per_syndrome
    ).sum(axis=1)
    coverage = coverage.reset_index()

    realm_cols = ["island_id", "biogeographic_realm"]
    for optional in ["island_latitude", "island_longitude"]:
        if optional in realm_assignment.columns:
            realm_cols.append(optional)

    queue = (
        realm_assignment.loc[
            realm_assignment["biogeographic_realm"].eq(target_realm), realm_cols
        ]
        .assign(island_id=lambda d: d["island_id"].astype(str))
        .merge(observation_support, on="island_id", how="left", suffixes=("_realm", ""))
        .merge(coverage, on="island_id", how="left")
    )
    queue = queue.loc[
        queue["flora_recorded"].fillna(False).astype(bool)
        & ~queue["island_id"].isin(resolved_native)
    ].copy()

    numeric_coverage = CORE_SYNDROMES + [
        "min_syndrome_species",
        "n_syndromes_meeting_support",
    ]
    for column in numeric_coverage:
        queue[column] = pd.to_numeric(queue[column], errors="coerce").fillna(0).astype(int)

    queue = queue.loc[
        queue["n_flora_species_recorded"].between(min_flora_species, max_flora_species)
        & queue["n_syndromes_meeting_support"].ge(min_syndromes_meeting_support)
    ].copy()

    queue["distance_bin"] = pd.cut(
        pd.to_numeric(queue["distance_to_continent_km"], errors="coerce"),
        [-1, 5, 25, 100, 500, float("inf")],
        labels=["0-5", "5-25", "25-100", "100-500", "500+"],
    ).astype("string")

    queue["priority_basis"] = (
        "outcome_blind:syndrome_availability+flora_size+distance_bin+spatial_block"
    )
    queue = queue.sort_values(
        [
            "n_syndromes_meeting_support",
            "min_syndrome_species",
            "n_flora_species_recorded",
            "island_id",
        ],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    queue.insert(0, "queue_rank", range(1, len(queue) + 1))
    return queue


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--status-flora", required=True)
    parser.add_argument("--realm-assignment", required=True)
    parser.add_argument("--species-syndrome", required=True)
    parser.add_argument("--observation-support", required=True)
    parser.add_argument("--target-realm", default="Nearctic")
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--manifest-json", required=True)
    parser.add_argument("--min-flora-species", type=int, default=20)
    parser.add_argument("--max-flora-species", type=int, default=400)
    parser.add_argument("--min-species-per-syndrome", type=int, default=10)
    parser.add_argument("--min-syndromes-meeting-support", type=int, default=4)
    args = parser.parse_args()

    queue = build_realm_replication_queue(
        pd.read_csv(args.status_flora),
        pd.read_csv(args.realm_assignment),
        pd.read_csv(args.species_syndrome),
        pd.read_csv(args.observation_support),
        target_realm=args.target_realm,
        min_flora_species=args.min_flora_species,
        max_flora_species=args.max_flora_species,
        min_species_per_syndrome=args.min_species_per_syndrome,
        min_syndromes_meeting_support=args.min_syndromes_meeting_support,
    )
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    queue.to_csv(output, index=False)

    manifest = {
        "contract": "chapter1_pr138_realm_replication_queue_v1",
        "target_realm": args.target_realm,
        "n_candidates": int(len(queue)),
        "ranking_uses_outcome_values": False,
        "ranking_uses_effect_estimates": False,
        "ranking_uses_p_values": False,
        "core_syndromes": CORE_SYNDROMES,
        "min_flora_species": args.min_flora_species,
        "max_flora_species": args.max_flora_species,
        "min_species_per_syndrome": args.min_species_per_syndrome,
        "min_syndromes_meeting_support": args.min_syndromes_meeting_support,
    }
    manifest_path = Path(args.manifest_json)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
