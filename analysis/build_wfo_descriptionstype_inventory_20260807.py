"""Acquire strict species descriptions from WFO's consolidated description export."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from build_wfo_bulk_description_inventory_20260807 import (
    build_wfo_identity_map,
    fixed_universe,
    normalized_description,
    sha256_file,
)


ELIGIBLE_TYPES = {"general", "morphology", "diagnostic"}


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def now_utc() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def write_file_manifest(output_dir: Path) -> None:
    rows = []
    for path in sorted(output_dir.iterdir()):
        if not path.is_file() or path.name == "file_manifest.sha256":
            continue
        rows.append(f"{sha256_file(path)}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text(
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-csv", type=Path, required=True)
    parser.add_argument("--backbone-zip", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--prior-descriptions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    universe = fixed_universe(args.coverage, args.master)
    _, identity_map, ambiguous, backbone_stats = build_wfo_identity_map(
        args.backbone_zip,
        universe,
    )
    prior = pd.read_csv(
        args.prior_descriptions,
        usecols=["wfo_id", "description_sha256"],
        dtype=str,
    ).fillna("")
    prior_keys = set(
        prior[["wfo_id", "description_sha256"]].itertuples(index=False, name=None)
    )

    type_counts: Counter[str] = Counter()
    rank_counts: Counter[str] = Counter()
    authority_rows: Counter[tuple[str, str]] = Counter()
    authority_species: dict[tuple[str, str], set[str]] = {}
    total_rows = 0
    eligible_species_rows = 0
    strict_identity_rows = 0
    duplicate_prior = 0
    duplicate_within_export = 0
    seen: set[tuple[str, str]] = set()
    records: list[dict[str, str]] = []

    with args.export_csv.open(
        encoding="utf-8-sig",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        expected = [
            "identifier",
            "type",
            "en",
            "description",
            "WFO-ID",
            "taxonRank",
            "accessRights",
            "rights",
            "rightsHolder",
            "authority_id",
            "identifier",
            "title",
        ]
        if header != expected:
            raise ValueError(f"unexpected DescriptionsType header: {header}")
        for row in reader:
            total_rows += 1
            if len(row) < 12:
                continue
            (
                source_record,
                description_type,
                language,
                raw_description,
                wfo_id,
                taxon_rank,
                access_rights,
                rights,
                rights_holder,
                authority_id,
                authority_identifier,
                authority_title,
            ) = row[:12]
            type_counts[description_type or "(blank)"] += 1
            rank_counts[taxon_rank or "(blank)"] += 1
            if taxon_rank.casefold() != "species":
                continue
            if description_type.casefold() not in ELIGIBLE_TYPES:
                continue
            eligible_species_rows += 1
            wfo_id = wfo_id.casefold().strip()
            identity = identity_map.get(wfo_id)
            if identity is None:
                continue
            strict_identity_rows += 1
            description = normalized_description(raw_description)
            if not description:
                continue
            description_sha256 = hashlib.sha256(
                description.encode("utf-8")
            ).hexdigest()
            content_key = (wfo_id, description_sha256)
            if content_key in prior_keys:
                duplicate_prior += 1
                continue
            if content_key in seen:
                duplicate_within_export += 1
                continue
            seen.add(content_key)
            authority_key = (
                authority_id.strip() or "unknown",
                authority_title.strip() or "unknown_not_recorded",
            )
            authority_rows[authority_key] += 1
            authority_species.setdefault(authority_key, set()).add(
                identity["accepted_species"]
            )
            provider = f"wfo_authority_{authority_key[0]}"
            records.append(
                {
                    "accepted_species": identity["accepted_species"],
                    "wfo_id": wfo_id,
                    "matched_wfo_name": identity["matched_wfo_name"],
                    "wfo_name_status": identity["wfo_name_status"],
                    "name_match_method": identity["name_match_method"],
                    "provider": provider,
                    "provider_title": authority_key[1],
                    "authority_id": authority_key[0],
                    "authority_identifier": authority_identifier,
                    "source_record": source_record,
                    "description_type": description_type,
                    "language": language,
                    "description": description,
                    "description_sha256": description_sha256,
                    "source_url": f"https://list.worldfloraonline.org/{wfo_id}",
                    "archive_url": (
                        "https://files.worldfloraonline.org/Files/Temp/"
                        "NLP_ML_Descriptions_Type/DescriptionsType.csv"
                    ),
                    "license": access_rights,
                    "rights": rights,
                    "rights_holder": rights_holder,
                    "retrieved_at_utc": now_utc(),
                    "source_lineage": (
                        f"provider_treatment:wfo_authority_{authority_key[0]}:{wfo_id}"
                    ),
                    "review_status": "source_description_acquired_unreviewed",
                }
            )

    descriptions = pd.DataFrame(records)
    provider_rows = []
    for key, count in authority_rows.most_common():
        provider_rows.append(
            {
                "authority_id": key[0],
                "authority_title": key[1],
                "matched_new_description_rows": count,
                "matched_master_species": len(authority_species[key]),
            }
        )
    providers = pd.DataFrame(provider_rows)
    args.output.mkdir(parents=True, exist_ok=True)
    descriptions.to_csv(
        args.output / "matched_new_master_species_descriptions.csv.gz",
        index=False,
        compression="gzip",
    )
    providers.to_csv(args.output / "provider_coverage.csv", index=False)
    ambiguous.to_csv(args.output / "ambiguous_wfo_identity_map.csv", index=False)

    summary = {
        "contract": "wfo_consolidated_description_export_inventory_v1",
        "created_at_utc": now_utc(),
        "formal_promotion_status": "not_promoted_acquisition_stage_only",
        "source_export": {
            "file": args.export_csv.name,
            "bytes": args.export_csv.stat().st_size,
            "sha256": sha256_file(args.export_csv),
            "url": (
                "https://files.worldfloraonline.org/Files/Temp/"
                "NLP_ML_Descriptions_Type/DescriptionsType.csv"
            ),
            "search_api_queries": 0,
            "cost_usd": 0.0,
        },
        "fixed_universe_species": 106_295,
        "backbone": backbone_stats,
        "export_rows": total_rows,
        "taxon_rank_counts": dict(rank_counts),
        "description_type_counts": dict(type_counts),
        "eligible_species_description_rows": eligible_species_rows,
        "strict_identity_rows_before_content_deduplication": strict_identity_rows,
        "duplicates_of_prior_archives_removed": duplicate_prior,
        "duplicates_within_export_removed": duplicate_within_export,
        "new_descriptions": {
            "rows": len(descriptions),
            "species": int(
                descriptions["accepted_species"].nunique()
                if not descriptions.empty
                else 0
            ),
            "wfo_ids": int(
                descriptions["wfo_id"].nunique() if not descriptions.empty else 0
            ),
            "authorities": int(
                descriptions["authority_id"].nunique()
                if not descriptions.empty
                else 0
            ),
        },
        "evidence_contract": {
            "accepted_rows": "none_at_this_stage",
            "identity": "species-rank exact accepted or exact synonym, family-consistent",
            "lineage": "original WFO authority/resource plus WFO taxon identifier",
            "prior_content_deduplication": "WFO ID plus normalized full-description SHA-256",
            "next_gate": "trait extraction, exact quote validation, conflict and cultivar review",
        },
    }
    (args.output / "run_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_file_manifest(args.output)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
