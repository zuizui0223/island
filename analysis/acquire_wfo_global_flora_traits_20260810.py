"""Acquire floral morphology from six official WFO-hosted flora archives.

This source-scale wave reprocesses complete pinned provider archives with the
audited shared morphology extractor. It excludes every species x trait already
resolved in the preceding formal artifact and writes a deterministic holdout
queue; it does not promote evidence or perform genus inference.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable, Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.acquire_flora_of_australia_traits_20260809 import (
    AXIS,
    _cultivar_or_hybrid_treatment,
    extract_description,
)
from analysis.acquire_wfo_kew_africa_traits_20260810 import (
    AUDIT_COLUMNS,
    BACKBONE_BYTES,
    BACKBONE_MEMBER,
    BACKBONE_SHA256,
    EXPECTED_UNIVERSE_SPECIES,
    _fixed_universe,
    _lineage_columns,
    _resolve_identity_map,
    _sha,
    _text,
    _validate_file,
    _write_manifest,
)
from analysis.build_wfo_bulk_description_inventory_20260807 import (
    dict_rows_from_zip,
    normalized_description,
    remap_rows,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import load_rules

SOURCE_RUN_ID = "wfo-global-six-floras-20260810"
SOURCE_ARTIFACT = "wfo-global-six-official-dwca-20260810"
ACCEPTANCE_CONTRACT = "official_wfo_dwca_exact_species_description_quote_v2"
AUDIT_SEED = "wfo-global-six-floras-holdout-audit-20260810-v1"
AUDIT_SIZE = 200
AUDIT_QUOTA = {
    "floral_form": 30,
    "floral_symmetry": 10,
    "flower_primary_color": 70,
    "flower_size_class": 30,
    "inflorescence_display": 40,
    "tube_depth_class": 20,
}

SOURCE_META = {
    "brazilian_flora_2020": {
        "archive": "brazilian_flora_2020.zip",
        "bytes": 13_266_709,
        "sha256": "79455efc837678d0812c4f83247c9f024ab4d8dae43fe72f5c40b2302d50c923",
        "member": "description.txt",
        "citation": "Brazilian Flora 2020, official WFO Darwin Core Archive",
    },
    "eflora_thailand": {
        "archive": "eflora_thailand.zip",
        "bytes": 2_208_210,
        "sha256": "26a135ce728fc6ebddb2576127113858f3222acbfbf00196aee4fae4736ce64e",
        "member": "description.txt",
        "citation": "e-Flora Thailand, official WFO Darwin Core Archive",
    },
    "flora_of_panama": {
        "archive": "flora_panama.zip",
        "bytes": 3_686_955,
        "sha256": "b5b7d67e6038be4aee9ff3e9ffcfb9764b314bb025fa87aa887bfe96951b9859",
        "member": "description.txt",
        "citation": "Flora of Panama, Missouri Botanical Garden, official WFO DwC-A",
    },
    "nybg_floraneotropica": {
        "archive": "nybg_floraneotropica.zip",
        "bytes": 4_711_792,
        "sha256": "29fd3502f3b75b46f07f80a3cd8a248d570a998f25811668d6ca451f6b075fe7",
        "member": "description.txt",
        "citation": "Flora Neotropica, New York Botanical Garden, official WFO DwC-A",
    },
    "nybg_memoirs": {
        "archive": "nybg_memoirs.zip",
        "bytes": 6_256_344,
        "sha256": "168015ec1f7fd17df00e1b5c07b4171453366451bc12e69630227f42a99ef93e",
        "member": "description.txt",
        "citation": "Memoirs of the New York Botanical Garden, official WFO DwC-A",
    },
    "rhododendron_monograph": {
        "archive": "rhododendron.zip",
        "bytes": 302_618,
        "sha256": "f6a2cee64e3d4f17acf6cafb667781680df9aa64e908687b1b1e951a499d66c0",
        "member": "RhododendronDescriptions.txt",
        "citation": "Rhododendron monograph, Royal Botanic Garden Edinburgh, official WFO DwC-A",
    },
}

ELIGIBLE_TYPES = {
    "brazilian_flora_2020": {"morphology"},
    "eflora_thailand": {"general"},
    "flora_of_panama": {"general"},
    "nybg_floraneotropica": {"general"},
    "nybg_memoirs": {"general"},
    "rhododendron_monograph": {"morphology"},
}


def _provider_rows(base: Path) -> Iterable[tuple[str, Iterator[dict[str, str]]]]:
    yield (
        "brazilian_flora_2020",
        dict_rows_from_zip(base / "brazilian_flora_2020.zip", "description.txt"),
    )
    yield (
        "eflora_thailand",
        dict_rows_from_zip(base / "eflora_thailand.zip", "description.txt"),
    )
    dwca_fields = [
        "taxonID",
        "description",
        "type",
        "source",
        "language",
        "created",
        "contributor",
        "audience",
        "rights",
        "license",
        "rightsHolder",
        "creator",
    ]
    yield (
        "flora_of_panama",
        dict_rows_from_zip(
            base / "flora_panama.zip",
            "description.txt",
            header=False,
            fieldnames=dwca_fields,
        ),
    )
    for provider, filename in (
        ("nybg_floraneotropica", "nybg_floraneotropica.zip"),
        ("nybg_memoirs", "nybg_memoirs.zip"),
    ):
        yield (
            provider,
            dict_rows_from_zip(base / filename, "description.txt", delimiter=","),
        )
    yield (
        "rhododendron_monograph",
        remap_rows(
            dict_rows_from_zip(
                base / "rhododendron.zip",
                "RhododendronDescriptions.txt",
            ),
            defaults={"language": "en"},
        ),
    )


def _description_inventory(
    archive_dir: Path,
    identity_map: dict[str, dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, str]] = []
    provider_rows: list[dict[str, Any]] = []
    for provider, rows in _provider_rows(archive_dir):
        total = eligible = matched = cultivar_excluded = 0
        seen: set[tuple[str, str, str]] = set()
        for row in rows:
            total += 1
            if _text(row.get("type")).casefold() not in ELIGIBLE_TYPES[provider]:
                continue
            eligible += 1
            wfo_id = _text(row.get("taxonID")).casefold()
            identity = identity_map.get(wfo_id)
            if identity is None:
                continue
            description = normalized_description(_text(row.get("description")))
            species = identity["accepted_species"]
            if not description:
                continue
            if _cultivar_or_hybrid_treatment(species, description):
                cultivar_excluded += 1
                continue
            fingerprint = _sha(description)
            duplicate_key = (species, wfo_id, fingerprint)
            if duplicate_key in seen:
                continue
            seen.add(duplicate_key)
            matched += 1
            records.append(
                {
                    "accepted_species": species,
                    "wfo_id": wfo_id,
                    "matched_wfo_name": identity["matched_wfo_name"],
                    "name_match_method": identity["name_match_method"],
                    "provider": provider,
                    "source_record": _text(row.get("source")),
                    "description": description,
                    "description_sha256": fingerprint,
                }
            )
        provider_rows.append(
            {
                "provider": provider,
                "description_rows_total": total,
                "description_rows_trait_eligible": eligible,
                "matched_description_rows": matched,
                "cultivar_or_hybrid_rows_excluded": cultivar_excluded,
            }
        )
    return pd.DataFrame(records), pd.DataFrame(provider_rows)


def build_evidence(
    inventory: pd.DataFrame,
    current_direct_pairs: set[tuple[str, str]],
    rules: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for record in inventory.sort_values(
        ["provider", "accepted_species", "wfo_id", "description_sha256"],
        kind="stable",
    ).itertuples(index=False):
        for trait, item in extract_description(record.description, rules).items():
            if (record.accepted_species, trait) in current_direct_pairs:
                continue
            quote = _text(item["quote"])
            if not quote or quote not in record.description:
                continue
            value = _text(item["value"])
            record_hash = _sha(
                f"{record.provider}|{record.accepted_species}|{trait}|{value}|{quote}"
            )[:24]
            meta = SOURCE_META[record.provider]
            rows.append(
                {
                    "accepted_species": record.accepted_species,
                    "axis": AXIS[trait],
                    "trait_name": trait,
                    "normalized_value": value,
                    "quality": "high",
                    "source_group": "latest_public_web",
                    "source_provider": record.provider,
                    "source_url": f"https://list.worldfloraonline.org/{record.wfo_id}",
                    "source_record_id": f"wfo-global-{record_hash}",
                    "source_citation": meta["citation"],
                    "source_excerpt": quote,
                    "evidence_scope": "species_direct",
                    "name_match_method": record.name_match_method,
                    "source_lineage": record.source_lineage,
                    "lineage_method": record.lineage_method,
                    "source_run_id": SOURCE_RUN_ID,
                    "source_artifact": SOURCE_ARTIFACT,
                    "source_file": f"{meta['archive']}::{meta['member']}",
                    "acceptance_contract": ACCEPTANCE_CONTRACT,
                }
            )
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
    if evidence.empty:
        return evidence
    return (
        evidence.drop_duplicates("source_record_id")
        .sort_values(
            ["trait_name", "accepted_species", "source_lineage", "source_record_id"],
            kind="stable",
        )
        .reset_index(drop=True)
    )


def deterministic_audit_sample(evidence: pd.DataFrame) -> pd.DataFrame:
    samples: list[pd.DataFrame] = []
    for trait, quota in AUDIT_QUOTA.items():
        candidates = evidence.loc[evidence["trait_name"].eq(trait)]
        quota = min(quota, len(candidates))
        seed = int(_sha(f"{AUDIT_SEED}|{trait}")[:8], 16)
        samples.append(candidates.sample(n=quota, random_state=seed))
    audit = pd.concat(samples, ignore_index=True)
    if len(audit) < AUDIT_SIZE:
        selected = set(audit["source_record_id"])
        remainder = evidence.loc[~evidence["source_record_id"].isin(selected)].copy()
        remainder["_audit_order"] = remainder["source_record_id"].map(
            lambda value: _sha(f"{AUDIT_SEED}|remainder|{value}")
        )
        needed = AUDIT_SIZE - len(audit)
        if len(remainder) < needed:
            raise ValueError(f"not enough total evidence rows for {AUDIT_SIZE}-row audit")
        audit = pd.concat(
            [
                audit,
                remainder.sort_values("_audit_order", kind="stable")
                .head(needed)
                .drop(columns="_audit_order"),
            ],
            ignore_index=True,
        )
    audit = audit.sort_values(
        ["trait_name", "accepted_species", "source_record_id"], kind="stable"
    ).rename(
        columns={
            "source_provider": "provider",
            "source_record_id": "candidate_id",
            "source_excerpt": "exact_supporting_quote",
        }
    )
    for column in (
        "accepted_correct",
        "cultivar_status",
        "reviewer",
        "reviewed_at_utc",
        "audit_reason",
    ):
        audit[column] = ""
    return audit.reindex(columns=AUDIT_COLUMNS)


def acquire(
    *,
    backbone_zip: Path,
    archive_dir: Path,
    master_csv: Path,
    strict_coverage_csv: Path,
    current_direct_ledger_csv: Path,
    measurement_config: Path,
    output_dir: Path,
) -> dict[str, Any]:
    _validate_file(backbone_zip, BACKBONE_BYTES, BACKBONE_SHA256)
    for meta in SOURCE_META.values():
        _validate_file(archive_dir / meta["archive"], meta["bytes"], meta["sha256"])
    universe = _fixed_universe(master_csv, strict_coverage_csv)
    identity_map, identity_stats = _resolve_identity_map(backbone_zip, universe)
    inventory, provider_stats = _description_inventory(archive_dir, identity_map)
    inventory = _lineage_columns(inventory)

    direct = pd.read_csv(current_direct_ledger_csv, dtype=str).fillna("")
    if "resolution_status" in direct:
        direct = direct.loc[direct["resolution_status"].eq("resolved")]
    current_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=True))
    evidence = build_evidence(inventory, current_pairs, load_rules(measurement_config))
    audit = deterministic_audit_sample(evidence)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output_dir / "wfo_global_six_evidence.csv.gz", index=False, compression="gzip")
    inventory.to_csv(
        output_dir / "wfo_global_six_treatment_inventory.csv.gz",
        index=False,
        compression="gzip",
    )
    provider_stats.to_csv(output_dir / "wfo_global_six_provider_coverage.csv", index=False)
    audit.to_csv(output_dir / "wfo_global_six_audit_queue_200.csv", index=False)

    summary: dict[str, Any] = {
        "contract": "wfo_global_six_source_package_v1",
        "generated_at_utc": datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        "fixed_universe_species": EXPECTED_UNIVERSE_SPECIES,
        "backbone": {
            "member": BACKBONE_MEMBER,
            "bytes": BACKBONE_BYTES,
            "sha256": BACKBONE_SHA256,
        },
        "source_archives": SOURCE_META,
        "identity": identity_stats,
        "matched_description_rows": len(inventory),
        "matched_species": int(inventory["accepted_species"].nunique()),
        "cross_provider_content_rows": int(
            inventory["lineage_method"]
            .eq("accepted_species_normalized_content_fingerprint_cross_provider")
            .sum()
        ),
        "incremental_evidence_rows": len(evidence),
        "incremental_species_trait": int(
            evidence.drop_duplicates(["accepted_species", "trait_name"]).shape[0]
        ),
        "incremental_species_axis": int(
            evidence.drop_duplicates(["accepted_species", "axis"]).shape[0]
        ),
        "incremental_species": int(evidence["accepted_species"].nunique()),
        "trait_rows": {
            str(key): int(value)
            for key, value in evidence["trait_name"].value_counts().sort_index().items()
        },
        "audit_rows": len(audit),
        "audit_seed": AUDIT_SEED,
        "audit_quota": AUDIT_QUOTA,
        "source_run_id": SOURCE_RUN_ID,
        "source_artifact": SOURCE_ARTIFACT,
        "quality": "high",
        "promotion_performed": False,
        "strict_exclusions": [
            "non_species_and_family_conflicting_names",
            "ambiguous_wfo_identity_clusters",
            "non_morphology_description_types",
            "cultivar_and_hybrid_treatments",
            "existing_direct_species_trait_pairs",
            "family_and_global_inference",
        ],
    }
    (output_dir / "wfo_global_six_source_package_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_manifest(output_dir)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbone-zip", type=Path, required=True)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=Path("data/v2/staging/gbif/collected/island_taxa.csv"),
    )
    parser.add_argument("--strict-coverage-csv", type=Path, required=True)
    parser.add_argument("--current-direct-ledger-csv", type=Path, required=True)
    parser.add_argument(
        "--measurement-config",
        type=Path,
        default=Path("config/measurement_classification.yml"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = acquire(
        backbone_zip=args.backbone_zip,
        archive_dir=args.archive_dir,
        master_csv=args.master_csv,
        strict_coverage_csv=args.strict_coverage_csv,
        current_direct_ledger_csv=args.current_direct_ledger_csv,
        measurement_config=args.measurement_config,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
