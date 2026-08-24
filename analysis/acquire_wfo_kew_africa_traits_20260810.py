"""Acquire floral morphology from three official Kew/WFO African floras.

The command is source-scale rather than query-scale: it verifies a fixed WFO
Plant List backbone and the official FTEA, Flora Zambesiaca, and FWTA Darwin
Core Archives, reconciles accepted names and synonyms against the fixed
106,295-species universe, extracts six visible floral traits, and writes a
deterministic 200-row audit queue.  Promotion remains a separate, review-gated
step in ``open-web-review-promote.yml``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import zipfile
from collections import defaultdict
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.acquire_flora_of_australia_traits_20260809 import (
    AXIS,
    _cultivar_or_hybrid_treatment,
    extract_description,
)
from analysis.build_wfo_bulk_description_inventory_20260807 import (
    ELIGIBLE_DESCRIPTION_TYPES,
    normalized_description,
    provider_rows_round3,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.trait_measurements import load_rules

EXPECTED_UNIVERSE_SPECIES = 106_295
BACKBONE_MEMBER = "classification.csv"
BACKBONE_BYTES = 121_660_019
BACKBONE_SHA256 = "ccfc7e3ca85c8b80aa96eb52d4e50e07e1811d28916b544280d8a4f1a043ba72"
SOURCE_RUN_ID = "wfo-kew-africa-three-floras-20260810"
SOURCE_ARTIFACT = "wfo-kew-africa-official-dwca-20260810"
ACCEPTANCE_CONTRACT = "official_wfo_dwca_exact_species_description_quote_v2"
AUDIT_SEED = "wfo-kew-africa-holdout-audit-20260810-v1"
AUDIT_SIZE = 200
AUDIT_QUOTA = {
    "floral_form": 25,
    "floral_symmetry": 5,
    "flower_primary_color": 80,
    "flower_size_class": 15,
    "inflorescence_display": 50,
    "tube_depth_class": 25,
}

SOURCE_META = {
    "kew_flora_tropical_east_africa": {
        "archive": "kew_ftea.zip",
        "bytes": 5_401_533,
        "sha256": "14e95ad3742225793a667f24e40c166db118dae4917269bb30f78178954f42a6",
        "citation": (
            "Flora of Tropical East Africa, Royal Botanic Gardens, Kew, "
            "official WFO DwC-A"
        ),
    },
    "kew_flora_zambesiaca_current": {
        "archive": "kew_fz.zip",
        "bytes": 3_960_398,
        "sha256": "5aec2c09934091db2f0fda65314b8b232ff4f291ac4846d8688e6900428c86e2",
        "citation": (
            "Flora Zambesiaca, Royal Botanic Gardens, Kew, official WFO DwC-A"
        ),
    },
    "kew_flora_west_tropical_africa": {
        "archive": "kew_fwta_wfo.zip",
        "bytes": 888_526,
        "sha256": "243cd67ae56196ff61b940949eda102c96a834fd41253b06a370ae802baf843c",
        "citation": (
            "Flora of West Tropical Africa, Royal Botanic Gardens, Kew, "
            "official WFO DwC-A"
        ),
    },
}

AUDIT_COLUMNS = [
    "accepted_species",
    "axis",
    "trait_name",
    "normalized_value",
    "quality",
    "source_group",
    "provider",
    "source_url",
    "candidate_id",
    "source_citation",
    "exact_supporting_quote",
    "evidence_scope",
    "name_match_method",
    "source_lineage",
    "lineage_method",
    "source_run_id",
    "source_artifact",
    "source_file",
    "acceptance_contract",
    "accepted_correct",
    "cultivar_status",
    "reviewer",
    "reviewed_at_utc",
    "audit_reason",
]


def _text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_file(path: Path, expected_bytes: int, expected_sha256: str) -> None:
    if path.stat().st_size != expected_bytes:
        raise ValueError(
            f"byte-size drift for {path.name}: {path.stat().st_size} != {expected_bytes}"
        )
    observed = _file_sha256(path)
    if observed != expected_sha256:
        raise ValueError(f"SHA-256 drift for {path.name}: {observed}")


def _backbone_rows(path: Path) -> Iterator[dict[str, str]]:
    with zipfile.ZipFile(path) as archive, archive.open(BACKBONE_MEMBER) as raw:
        handle = io.TextIOWrapper(
            raw,
            encoding="utf-8-sig",
            errors="replace",
            newline="",
        )
        yield from csv.DictReader(handle, delimiter="\t")


def _fixed_universe(master_csv: Path, strict_coverage_csv: Path) -> pd.DataFrame:
    strict = pd.read_csv(
        strict_coverage_csv,
        dtype=str,
        usecols=["accepted_species"],
    ).fillna("")
    names = set(strict["accepted_species"])
    if len(names) != EXPECTED_UNIVERSE_SPECIES:
        raise ValueError("strict coverage does not contain the fixed 106,295 species")
    master = pd.read_csv(
        master_csv,
        dtype=str,
        usecols=["accepted_species", "family"],
    ).fillna("")
    master = master.loc[master["accepted_species"].isin(names)].drop_duplicates(
        "accepted_species"
    )
    if len(master) != EXPECTED_UNIVERSE_SPECIES or master["family"].eq("").any():
        raise ValueError("master taxonomy does not match the fixed strict universe")
    return master


def _resolve_identity_map(
    backbone_zip: Path,
    universe: pd.DataFrame,
) -> tuple[dict[str, dict[str, str]], dict[str, int]]:
    """Resolve exact fixed names, then their WFO synonym clusters, fail closed."""

    master = universe.set_index("accepted_species")["family"].to_dict()
    candidates: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    exact_rows = 0
    for row in _backbone_rows(backbone_zip):
        name = _text(row.get("scientificName"))
        if name not in master:
            continue
        if _text(row.get("taxonRank")).casefold() != "species":
            continue
        if _text(row.get("family")).casefold() != master[name].casefold():
            continue
        status = _text(row.get("taxonomicStatus")).casefold()
        if status not in {"accepted", "synonym"}:
            continue
        exact_rows += 1
        taxon_id = _text(row.get("taxonID")).casefold()
        accepted_id = _text(row.get("acceptedNameUsageID")).casefold()
        method = (
            "wfo_june_2026_exact_accepted_species_family_consistent"
            if status == "accepted"
            else "wfo_june_2026_exact_synonym_species_family_consistent"
        )
        item = {
            "accepted_species": name,
            "matched_wfo_name": name,
            "name_match_method": method,
        }
        if taxon_id:
            candidates[taxon_id].append(item)
        if status == "synonym" and accepted_id:
            candidates[accepted_id].append(
                {
                    **item,
                    "name_match_method": "wfo_accepted_usage_from_exact_fixed_synonym",
                }
            )

    direct_map = _unique_identity_map(candidates)
    accepted_clusters = {
        key: value
        for key, value in direct_map.items()
        if key.startswith("wfo-")
    }
    synonym_rows = 0
    for row in _backbone_rows(backbone_zip):
        if _text(row.get("taxonomicStatus")).casefold() != "synonym":
            continue
        accepted_id = _text(row.get("acceptedNameUsageID")).casefold()
        identity = accepted_clusters.get(accepted_id)
        if identity is None:
            continue
        if _text(row.get("taxonRank")).casefold() != "species":
            continue
        family = master[identity["accepted_species"]]
        if _text(row.get("family")).casefold() != family.casefold():
            continue
        taxon_id = _text(row.get("taxonID")).casefold()
        if not taxon_id:
            continue
        synonym_rows += 1
        candidates[taxon_id].append(
            {
                "accepted_species": identity["accepted_species"],
                "matched_wfo_name": _text(row.get("scientificName")),
                "name_match_method": "wfo_synonym_cluster_from_exact_fixed_name",
            }
        )

    identity_map = _unique_identity_map(candidates)
    stats = {
        "exact_fixed_backbone_rows": exact_rows,
        "strict_synonym_cluster_rows": synonym_rows,
        "resolved_wfo_ids": len(identity_map),
        "resolved_fixed_species": len(
            {value["accepted_species"] for value in identity_map.values()}
        ),
        "ambiguous_wfo_ids_excluded": len(candidates) - len(identity_map),
    }
    return identity_map, stats


def _unique_identity_map(
    candidates: dict[str, list[dict[str, str]]],
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for wfo_id, rows in candidates.items():
        if len({row["accepted_species"] for row in rows}) != 1:
            continue
        result[wfo_id] = min(
            rows,
            key=lambda row: (row["name_match_method"], row["matched_wfo_name"]),
        )
    return result


def _description_inventory(
    archive_dir: Path,
    identity_map: dict[str, dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, str]] = []
    provider_rows: list[dict[str, Any]] = []
    for provider, archive_url, rows in provider_rows_round3(archive_dir):
        total = eligible = matched = cultivar_excluded = 0
        seen: set[tuple[str, str, str]] = set()
        for row in rows:
            total += 1
            if _text(row.get("type")).casefold() not in ELIGIBLE_DESCRIPTION_TYPES[provider]:
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
                    "archive_url": archive_url,
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


def _lineage_columns(inventory: pd.DataFrame) -> pd.DataFrame:
    cross_provider = (
        inventory.groupby(["accepted_species", "description_sha256"])["provider"]
        .nunique()
        .gt(1)
        .to_dict()
    )
    inventory = inventory.copy()
    inventory["source_lineage"] = inventory.apply(
        lambda row: (
            "cross_provider_content:"
            f"{_sha(row['accepted_species'] + '|' + row['description_sha256'])[:32]}"
            if cross_provider.get(
                (row["accepted_species"], row["description_sha256"]), False
            )
            else f"provider_treatment:{row['provider']}:{row['accepted_species']}"
        ),
        axis=1,
    )
    inventory["lineage_method"] = inventory.apply(
        lambda row: (
            "accepted_species_normalized_content_fingerprint_cross_provider"
            if cross_provider.get(
                (row["accepted_species"], row["description_sha256"]), False
            )
            else "provider_accepted_species_treatment"
        ),
        axis=1,
    )
    return inventory


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
        extracted = extract_description(record.description, rules)
        for trait, item in extracted.items():
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
                    "source_record_id": f"wfo-africa-{record_hash}",
                    "source_citation": meta["citation"],
                    "source_excerpt": quote,
                    "evidence_scope": "species_direct",
                    "name_match_method": record.name_match_method,
                    "source_lineage": record.source_lineage,
                    "lineage_method": record.lineage_method,
                    "source_run_id": SOURCE_RUN_ID,
                    "source_artifact": SOURCE_ARTIFACT,
                    "source_file": f"{meta['archive']}::description.txt",
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
        ["trait_name", "accepted_species", "source_record_id"],
        kind="stable",
    )
    audit = audit.rename(
        columns={
            "source_provider": "provider",
            "source_record_id": "candidate_id",
            "source_excerpt": "exact_supporting_quote",
        }
    )
    audit["accepted_correct"] = ""
    audit["cultivar_status"] = ""
    audit["reviewer"] = ""
    audit["reviewed_at_utc"] = ""
    audit["audit_reason"] = ""
    return audit.reindex(columns=AUDIT_COLUMNS)


def _write_manifest(output_dir: Path) -> None:
    lines = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "file_manifest.sha256":
            lines.append(f"{_file_sha256(path)}  {path.name}")
    (output_dir / "file_manifest.sha256").write_text(
        "\n".join(lines) + "\n",
        encoding="ascii",
    )


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
    evidence_path = output_dir / "wfo_kew_africa_evidence.csv.gz"
    inventory_path = output_dir / "wfo_kew_africa_treatment_inventory.csv.gz"
    provider_path = output_dir / "wfo_kew_africa_provider_coverage.csv"
    audit_path = output_dir / "wfo_kew_africa_audit_queue_200.csv"
    evidence.to_csv(evidence_path, index=False, compression="gzip")
    inventory.to_csv(inventory_path, index=False, compression="gzip")
    provider_stats.to_csv(provider_path, index=False)
    audit.to_csv(audit_path, index=False)

    summary: dict[str, Any] = {
        "contract": "wfo_kew_africa_source_package_v2",
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
    summary_path = output_dir / "wfo_kew_africa_source_package_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
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
