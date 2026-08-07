"""Create a 106,295-species acquisition outcome ledger and next-source queue."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd


def now_utc() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def joined_unique(series: pd.Series) -> str:
    return "|".join(sorted({str(value) for value in series if str(value)}))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--wfo-name-matches", type=Path, required=True)
    parser.add_argument("--descriptions", type=Path, nargs="+", required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    species = coverage[["accepted_species"]].drop_duplicates()
    if len(species) != 106_295:
        raise ValueError(f"fixed universe mismatch: {len(species)}")
    master = pd.read_csv(args.master, dtype=str).fillna("")
    master = master.loc[
        master["accepted_species"].isin(set(species["accepted_species"]))
    ].drop_duplicates("accepted_species")
    ledger = species.merge(master, on="accepted_species", how="left", validate="1:1")

    strict_name_matches = pd.read_csv(args.wfo_name_matches, dtype=str).fillna("")
    strict_name_matches = strict_name_matches.loc[
        strict_name_matches["family_match"].eq("true")
        & strict_name_matches["rank_match"].eq("true")
        & strict_name_matches["wfo_taxonomic_status"].isin(["Accepted", "Synonym"])
    ]
    name_summary = (
        strict_name_matches.groupby("accepted_species", sort=False)
        .agg(
            strict_wfo_ids=("wfo_id", joined_unique),
            wfo_name_statuses=("wfo_taxonomic_status", joined_unique),
            strict_wfo_name_rows=("wfo_id", "size"),
        )
        .reset_index()
    )
    ledger = ledger.merge(name_summary, on="accepted_species", how="left")

    descriptions = pd.concat(
        [pd.read_csv(path, dtype=str).fillna("") for path in args.descriptions],
        ignore_index=True,
    )
    description_summary = (
        descriptions.groupby("accepted_species", sort=False)
        .agg(
            acquired_description_rows=("description_sha256", "size"),
            acquired_description_providers=("provider", "nunique"),
            acquired_description_languages=("language", joined_unique),
            acquired_description_types=("description_type", joined_unique),
            acquired_wfo_ids=("wfo_id", joined_unique),
        )
        .reset_index()
    )
    ledger = ledger.merge(description_summary, on="accepted_species", how="left")

    candidates = pd.read_csv(args.candidates, dtype=str).fillna("")
    candidate_summary = (
        candidates.groupby("accepted_species", sort=False)
        .agg(
            candidate_rows=("candidate_id", "size"),
            candidate_traits=("trait_name", joined_unique),
            candidate_axes=("axis", joined_unique),
            candidate_providers=("provider", "nunique"),
        )
        .reset_index()
    )
    ledger = ledger.merge(candidate_summary, on="accepted_species", how="left")

    formal = (
        coverage.assign(_resolved=coverage["quality"].ne(""))
        .groupby("accepted_species", sort=False)
        .agg(
            formal_resolved_axes=("_resolved", "sum"),
            formal_high_axes=("quality", lambda x: int((x == "high").sum())),
            formal_medium_axes=("quality", lambda x: int((x == "medium").sum())),
            formal_low_axes=("quality", lambda x: int((x == "low").sum())),
        )
        .reset_index()
    )
    ledger = ledger.merge(formal, on="accepted_species", how="left")

    integer_columns = [
        "strict_wfo_name_rows",
        "acquired_description_rows",
        "acquired_description_providers",
        "candidate_rows",
        "candidate_providers",
        "formal_resolved_axes",
        "formal_high_axes",
        "formal_medium_axes",
        "formal_low_axes",
    ]
    for column in integer_columns:
        ledger[column] = (
            pd.to_numeric(ledger[column], errors="coerce").fillna(0).astype(int)
        )
    for column in (
        "strict_wfo_ids",
        "wfo_name_statuses",
        "acquired_description_languages",
        "acquired_description_types",
        "acquired_wfo_ids",
        "candidate_traits",
        "candidate_axes",
    ):
        ledger[column] = ledger[column].fillna("")

    ledger["unresolved_formal_axes"] = 3 - ledger["formal_resolved_axes"]
    ledger["acquisition_outcome"] = "strict_wfo_identity_unresolved"
    identity = ledger["strict_wfo_name_rows"].gt(0)
    has_description = ledger["acquired_description_rows"].gt(0)
    has_candidate = ledger["candidate_rows"].gt(0)
    ledger.loc[identity, "acquisition_outcome"] = (
        "strict_wfo_identity_resolved_no_eligible_description"
    )
    ledger.loc[has_description, "acquisition_outcome"] = (
        "description_acquired_no_supported_trait_candidate"
    )
    ledger.loc[has_candidate, "acquisition_outcome"] = (
        "trait_candidate_extracted_review_pending"
    )
    ledger["recommended_next_lane"] = (
        "strict_name_reconciliation_with_second_backbone"
    )
    ledger.loc[identity, "recommended_next_lane"] = (
        "source_first_open_web_or_current_portal_description"
    )
    ledger.loc[has_description, "recommended_next_lane"] = (
        "source_grounded_multilingual_llm_extraction"
    )
    ledger.loc[has_candidate, "recommended_next_lane"] = (
        "quote_context_conflict_and_cultivar_review"
    )
    ledger["_n_islands"] = pd.to_numeric(
        ledger["n_islands"], errors="coerce"
    ).fillna(0)
    stage_priority = {
        "description_acquired_no_supported_trait_candidate": 0,
        "trait_candidate_extracted_review_pending": 1,
        "strict_wfo_identity_resolved_no_eligible_description": 2,
        "strict_wfo_identity_unresolved": 3,
    }
    ledger["_stage_priority"] = ledger["acquisition_outcome"].map(stage_priority)
    ledger = ledger.sort_values(
        [
            "_stage_priority",
            "unresolved_formal_axes",
            "_n_islands",
            "accepted_species",
        ],
        ascending=[True, False, False, True],
    ).drop(columns=["_stage_priority", "_n_islands"])
    ledger.insert(0, "acquisition_priority_rank", range(1, len(ledger) + 1))

    summary = {
        "contract": "wfo_full_universe_acquisition_status_v1",
        "created_at_utc": now_utc(),
        "fixed_universe_species": 106_295,
        "fixed_species_axis_denominator": 318_885,
        "outcomes": {
            key: int(value)
            for key, value in ledger["acquisition_outcome"].value_counts().items()
        },
        "strict_wfo_identity_species": int(identity.sum()),
        "species_with_acquired_descriptions": int(has_description.sum()),
        "species_with_trait_candidates": int(has_candidate.sum()),
        "species_without_acquired_descriptions": int((~has_description).sum()),
        "species_without_trait_candidates": int((~has_candidate).sum()),
        "formal_promotion_status": "unchanged_review_pending",
        "policy": {
            "family_inference": False,
            "global_fallback": False,
            "missing_as_absence": False,
            "candidate_rows_autoaccepted": 0,
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    ledger_path = args.output / "full_universe_acquisition_status.csv.gz"
    ledger.to_csv(ledger_path, index=False, compression="gzip")
    (args.output / "run_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (args.output / "file_manifest.sha256").write_text(
        "\n".join(
            [
                f"{sha256_file(ledger_path)}  {ledger_path.name}",
                f"{sha256_file(args.output / 'run_manifest.json')}  run_manifest.json",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
