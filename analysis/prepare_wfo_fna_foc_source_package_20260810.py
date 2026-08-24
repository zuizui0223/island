"""Build the audited FNA, Flora of China, and Flora d'Afrique package.

The archive descriptions are acquired separately and stored as common-schema
evidence.  This step applies the final quote-completeness gates, deterministically
replaces filtered audit rows inside each provider/trait stratum, and combines the
result with the preceding WFO high-yield package.  The combined audit is used so
that production approval remains trait-specific and is based on 320 reviewed
rows rather than treating each flora as an independent source-quality claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.prepare_wfo_high_yield_source_package_20260810 import (
    safe_colour_row,
    safe_display_row,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

SOURCE_RUN_ID = "wfo-high-yield-plus-fna-foc-africa-20260810"
SOURCE_ARTIFACT = "wfo-high-yield-plus-fna-foc-africa-source-package-20260810"
AUDIT_REPLACEMENT_SEED = "wfo-fna-foc-audit-v2-replacement"
PACKAGE_TRAITS = {
    "floral_form",
    "flower_primary_color",
    "flower_size_class",
    "inflorescence_display",
}
AUDIT_PROVENANCE_COLUMNS = {
    "accepted_correct",
    "cultivar_status",
    "reviewer",
    "reviewed_at_utc",
    "audit_reason",
}


def filter_incremental_evidence(evidence: pd.DataFrame) -> pd.DataFrame:
    """Apply fail-closed colour and display gates to extracted evidence."""

    frame = evidence.copy().fillna("")
    missing = set(EVIDENCE_COLUMNS).difference(frame.columns)
    if missing:
        raise ValueError(f"extracted evidence is missing columns: {sorted(missing)}")
    frame = frame.loc[frame["trait_name"].isin(PACKAGE_TRAITS)].copy()
    colour = frame["trait_name"].eq("flower_primary_color")
    display = frame["trait_name"].eq("inflorescence_display")
    frame = frame.loc[
        (~colour | frame.apply(safe_colour_row, axis=1))
        & (~display | frame.apply(safe_display_row, axis=1))
    ].copy()
    frame = frame.drop_duplicates("source_record_id")
    frame = frame.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"],
        keep="first",
    )
    return frame.sort_values(
        ["trait_name", "source_provider", "accepted_species", "source_record_id"],
        kind="stable",
    ).reset_index(drop=True)


def deterministic_refill_audit(
    evidence: pd.DataFrame,
    seed_audit: pd.DataFrame,
) -> pd.DataFrame:
    """Preserve each provider/trait quota while replacing filtered rows."""

    seed = seed_audit.copy().fillna("")
    if seed["candidate_id"].duplicated().any():
        raise ValueError("audit seed contains duplicate candidate IDs")
    evidence_ids = set(evidence["source_record_id"])
    kept = seed.loc[seed["candidate_id"].isin(evidence_ids)].copy()
    used_seed_ids = set(seed["candidate_id"])
    parts = [kept]
    for (trait, provider), quota in seed.groupby(
        ["trait_name", "provider"], sort=True
    ).size().items():
        have = len(
            kept.loc[
                kept["trait_name"].eq(trait) & kept["provider"].eq(provider)
            ]
        )
        need = int(quota) - have
        if need <= 0:
            continue
        pool = evidence.loc[
            evidence["trait_name"].eq(trait)
            & evidence["source_provider"].eq(provider)
            & ~evidence["source_record_id"].isin(used_seed_ids)
        ].copy()
        pool["_audit_order"] = pool["source_record_id"].map(
            lambda value, trait_name=trait, provider_name=provider: hashlib.sha256(
                f"{AUDIT_REPLACEMENT_SEED}|{trait_name}|{provider_name}|{value}".encode()
            ).hexdigest()
        )
        replacement = pool.sort_values("_audit_order", kind="stable").head(need)
        if len(replacement) != need:
            raise ValueError(f"not enough replacements for {provider}/{trait}")
        replacement = replacement.drop(columns="_audit_order").rename(
            columns={
                "source_provider": "provider",
                "source_record_id": "candidate_id",
                "source_excerpt": "exact_supporting_quote",
            }
        )
        for column in AUDIT_PROVENANCE_COLUMNS:
            replacement[column] = ""
        parts.append(replacement.reindex(columns=seed.columns))
    audit = pd.concat(parts, ignore_index=True).sort_values(
        ["trait_name", "provider", "accepted_species", "candidate_id"],
        kind="stable",
    )
    if len(audit) != len(seed) or audit["candidate_id"].duplicated().any():
        raise ValueError("refilled audit does not preserve seed size and uniqueness")
    return audit.reset_index(drop=True)


def _validate_manual_audit(queue: pd.DataFrame, manual: pd.DataFrame) -> None:
    if set(queue["candidate_id"]) != set(manual["candidate_id"]):
        raise ValueError("manual audit IDs do not exactly match the deterministic queue")
    required = AUDIT_PROVENANCE_COLUMNS | {"candidate_id", "trait_name"}
    missing = required.difference(manual.columns)
    if missing:
        raise ValueError(f"manual audit is missing columns: {sorted(missing)}")
    if manual[list(AUDIT_PROVENANCE_COLUMNS)].fillna("").eq("").any().any():
        raise ValueError("manual audit has incomplete decisions or review provenance")


def run(args: argparse.Namespace) -> dict[str, Any]:
    incremental = filter_incremental_evidence(
        pd.read_csv(args.extracted_evidence, dtype=str).fillna("")
    )
    queue = deterministic_refill_audit(
        incremental,
        pd.read_csv(args.audit_seed, dtype=str).fillna(""),
    )
    manual = pd.read_csv(args.manual_audit, dtype=str).fillna("")
    _validate_manual_audit(queue, manual)
    first_evidence = pd.read_csv(args.first_evidence, dtype=str).fillna("")
    first_audit = pd.read_csv(args.first_audit, dtype=str).fillna("")
    combined_evidence = pd.concat(
        [first_evidence.reindex(columns=EVIDENCE_COLUMNS), incremental],
        ignore_index=True,
    ).drop_duplicates("source_record_id")
    combined_audit = pd.concat([first_audit, manual], ignore_index=True)
    if combined_audit["candidate_id"].duplicated().any():
        raise ValueError("combined audit contains duplicate candidate IDs")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    incremental.to_csv(
        args.output_dir / "wfo_fna_foc_africa_incremental_evidence.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    queue.to_csv(args.output_dir / "wfo_fna_foc_africa_audit_queue_120.csv", index=False)
    combined_evidence.to_csv(
        args.output_dir / "wfo_combined_high_yield_evidence.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    combined_audit.to_csv(
        args.output_dir / "wfo_combined_high_yield_manual_audit_320.csv", index=False
    )
    summary = {
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "incremental_candidate_rows": len(incremental),
        "incremental_species": int(incremental["accepted_species"].nunique()),
        "incremental_species_trait": int(
            incremental[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "combined_candidate_rows": len(combined_evidence),
        "combined_reviewed": len(combined_audit),
        "incremental_reviewed": len(manual),
        "incremental_accepted_correct": int(
            manual["accepted_correct"].str.casefold().isin({"true", "1", "yes"}).sum()
        ),
    }
    (args.output_dir / "wfo_combined_high_yield_source_package_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extracted-evidence", type=Path, required=True)
    parser.add_argument("--audit-seed", type=Path, required=True)
    parser.add_argument("--manual-audit", type=Path, required=True)
    parser.add_argument("--first-evidence", type=Path, required=True)
    parser.add_argument("--first-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), sort_keys=True))
