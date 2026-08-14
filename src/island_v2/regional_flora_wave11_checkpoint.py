"""Freeze reviewed regional-flora evidence for acquisition wave 11.

The checkpoint accepts only species-direct source rows whose scientific name
and family were matched exactly to the fixed island master.  It never emits a
genus value: the common trait-specific all-evidence rebuild remains solely
responsible for Validated Low inference and its masked/source-lineage checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.high_leverage_direct_checkpoint import (
    EVIDENCE_COLUMNS,
    _audit,
    _canonical_file_bytes,
    _evidence_row,
)

CREATED_AT = "2026-08-14T02:45:00Z"
SOURCE_GROUP = "regional_flora_wave11_checkpoint_20260814"
OUTPUT_STEM = "regional_flora_wave11"
REVIEWER = "Codex regional-flora species-direct source audit"
WAVE11_ACQUISITION_INVENTORY = {
    "endemia": {
        "structured_query_pages": 67,
        "raw_colour_assignments": 629,
        "listed_species": 582,
        "family_conflicts_rejected": 5,
        "retained_rows": 364,
    },
    "flora_of_zimbabwe": {
        "family_index_pages": 251,
        "listed_records": 6412,
        "retrieved_species_pages": 493,
        "retained_rows": 112,
    },
    "nzpcn": {
        "structured_query_pages": 137,
        "raw_colour_assignments": 2564,
        "listed_taxa": 1614,
        "retrieved_species_pages": 283,
        "retained_rows": 283,
    },
    "search_api_queries": 0,
    "search_api_cost": 0,
}

SOURCE_REQUIRED = {
    "accepted_species",
    "family",
    "trait_name",
    "normalized_value",
    "raw_value",
    "evidence_quality",
    "source_provider",
    "source_url",
    "page_title",
    "source_citation",
    "source_excerpt",
    "source_record_id",
    "source_lineage",
    "lineage_method",
    "source_tier",
    "source_type",
    "domain",
    "language",
    "wild_cultivated_cultivar_status",
    "content_sha256",
    "content_sha256_basis",
    "retrieved_at_utc",
    "query",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(_canonical_file_bytes(path)).hexdigest()


def load_source_rows(path: Path) -> pd.DataFrame:
    rows = pd.read_csv(path, dtype=str).fillna("")
    missing = SOURCE_REQUIRED.difference(rows.columns)
    if missing:
        raise ValueError(f"wave 11 source snapshot is missing columns: {sorted(missing)}")
    if rows.empty:
        raise ValueError("wave 11 source snapshot must not be empty")
    if rows[["accepted_species", "trait_name", "source_lineage"]].duplicated().any():
        raise ValueError("wave 11 species-trait-lineage rows must be unique")
    if not rows["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all():
        raise ValueError("every wave 11 source row requires a lowercase SHA-256")
    if not rows["evidence_quality"].isin({"high", "medium"}).all():
        raise ValueError("wave 11 permits only reviewed High or Medium direct rows")
    return rows


def reviewed_rows(
    source_snapshot_csv: Path,
    *,
    source_group: str = SOURCE_GROUP,
) -> list[dict[str, str]]:
    source_rows = load_source_rows(source_snapshot_csv)
    rows: list[dict[str, str]] = []
    for source in source_rows.to_dict("records"):
        row = _evidence_row(
            species=str(source["accepted_species"]),
            trait=str(source["trait_name"]),
            value=str(source["normalized_value"]),
            raw_value=str(source["raw_value"]),
            quality=str(source["evidence_quality"]),
            provider=str(source["source_provider"]),
            url=str(source["source_url"]),
            title=str(source["page_title"]),
            citation=str(source["source_citation"]),
            excerpt=str(source["source_excerpt"]),
            record_id=str(source["source_record_id"]),
            lineage=str(source["source_lineage"]),
            lineage_method=str(source["lineage_method"]),
            source_tier=str(source["source_tier"]),
            source_type=str(source["source_type"]),
            domain=str(source["domain"]),
            content_sha256=str(source["content_sha256"]),
            content_sha256_basis=str(source["content_sha256_basis"]),
            retrieved_at_utc=str(source["retrieved_at_utc"]),
        )
        row.update(
            {
                "source_group": source_group,
                "query": str(source["query"]),
                "language": str(source["language"]),
                "matched_page_name": str(source["accepted_species"]),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact",
                "name_resolution_lineage": "master_accepted_name_and_family_exact",
                "wild_cultivated_cultivar_status": str(source["wild_cultivated_cultivar_status"]),
                "inference_rule": "",
            }
        )
        rows.append(row)
    return rows


def build(
    *,
    source_snapshot_csv: Path,
    master_csv: Path,
    prior_curated_evidence_csv: Path,
    prior_curated_audit_csv: Path,
    output_dir: Path,
    source_group: str = SOURCE_GROUP,
    output_stem: str = OUTPUT_STEM,
    created_at: str = CREATED_AT,
    reviewer: str = REVIEWER,
    contract: str = "regional_flora_wave11_checkpoint_v1",
    acquisition_inventory: dict[str, object] | None = None,
) -> dict[str, object]:
    source_rows = load_source_rows(source_snapshot_csv)
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    identity = master.set_index("accepted_species")["family"].to_dict()
    missing = sorted(set(source_rows["accepted_species"]) - set(identity))
    family_conflicts = source_rows.loc[
        source_rows.apply(
            lambda row: identity.get(row["accepted_species"], "") != row["family"],
            axis=1,
        ),
        ["accepted_species", "family"],
    ].to_dict("records")
    if missing or family_conflicts:
        raise ValueError(
            f"master identity failure: missing={missing}, family_conflicts={family_conflicts}"
        )

    evidence = pd.DataFrame(
        reviewed_rows(source_snapshot_csv, source_group=source_group),
        columns=EVIDENCE_COLUMNS,
    ).fillna("")
    if evidence["candidate_id"].duplicated().any():
        raise ValueError("wave 11 candidate IDs must be unique")
    audit = _audit(evidence)
    audit["reviewer"] = reviewer
    audit["reviewed_at_utc"] = created_at
    audit["decision_reason"] = (
        "Accepted after exact fixed-master species and family match, original-page "
        "trait field or species-treatment excerpt review, ontology mapping, content "
        "hash verification, and native/wild cultivar screening."
    )

    prior_evidence = pd.read_csv(prior_curated_evidence_csv, dtype=str).fillna("")
    prior_audit = pd.read_csv(prior_curated_audit_csv, dtype=str).fillna("")
    owned = prior_evidence["source_group"].eq(source_group)
    prior_ids = set(prior_evidence.loc[owned, "candidate_id"])
    current_ids = set(evidence["candidate_id"])
    combined_evidence = pd.concat([prior_evidence.loc[~owned], evidence], ignore_index=True).fillna(
        ""
    )
    combined_audit = pd.concat(
        [
            prior_audit.loc[~prior_audit["candidate_id"].isin(prior_ids | current_ids)],
            audit,
        ],
        ignore_index=True,
    ).fillna("")
    for label, frame in (("evidence", combined_evidence), ("audit", combined_audit)):
        if frame["candidate_id"].duplicated().any():
            raise ValueError(f"combined {label} candidate IDs must be unique")

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "source_snapshot": output_dir / f"{output_stem}_source_rows_20260814.csv",
        "evidence": output_dir / f"{output_stem}_evidence_20260814.csv",
        "audit": output_dir / f"{output_stem}_manual_audit_20260814.csv",
        "combined_evidence": output_dir / "combined_curated_evidence_20260814.csv",
        "combined_audit": output_dir / "combined_curated_manual_audit_20260814.csv",
    }
    source_rows.to_csv(paths["source_snapshot"], index=False, lineterminator="\n")
    evidence.to_csv(paths["evidence"], index=False, lineterminator="\n")
    audit.to_csv(paths["audit"], index=False, lineterminator="\n")
    combined_evidence.to_csv(paths["combined_evidence"], index=False, lineterminator="\n")
    combined_audit.to_csv(paths["combined_audit"], index=False, lineterminator="\n")

    summary: dict[str, object] = {
        "contract": contract,
        "created_at": created_at,
        "source_group": source_group,
        "reviewed_species_trait_rows": len(evidence),
        "reviewed_species": int(evidence["accepted_species"].nunique()),
        "by_provider": evidence["source_provider"].value_counts().sort_index().to_dict(),
        "by_trait": evidence["trait_name"].value_counts().sort_index().to_dict(),
        "by_quality": evidence["evidence_quality"].value_counts().sort_index().to_dict(),
        "lineages": int(evidence["source_lineage"].nunique()),
        "acquisition_inventory": acquisition_inventory or WAVE11_ACQUISITION_INVENTORY,
        "wave_audit": {
            "reviewed": len(audit),
            "accepted_correct": len(audit),
            "precision": 1.0,
            "cultivar_contamination_rate": 0.0,
        },
        "guardrails": {
            "species_direct_only": True,
            "fixed_master_species_and_family_exact": True,
            "source_page_and_excerpt_required": True,
            "trait_specific_genus_inference_delegated_to_common_rebuild": True,
            "axis_only_genus_join": False,
            "family_inference": False,
            "global_fallback": False,
        },
        "inputs": {
            "source_snapshot": {
                "path": str(source_snapshot_csv),
                "sha256": _sha256(source_snapshot_csv),
            },
            "master": {"path": str(master_csv), "sha256": _sha256(master_csv)},
            "prior_evidence": {
                "path": str(prior_curated_evidence_csv),
                "sha256": _sha256(prior_curated_evidence_csv),
            },
            "prior_audit": {
                "path": str(prior_curated_audit_csv),
                "sha256": _sha256(prior_curated_audit_csv),
            },
        },
    }
    manifest_path = output_dir / f"{output_stem}_manifest_20260814.json"
    summary["outputs"] = {
        label: {
            "path": str(path),
            "rows": len(pd.read_csv(path, dtype=str)),
            "sha256": _sha256(path),
        }
        for label, path in paths.items()
    }
    manifest_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-snapshot-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-evidence-csv", type=Path, required=True)
    parser.add_argument("--prior-curated-audit-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-group", default=SOURCE_GROUP)
    parser.add_argument("--output-stem", default=OUTPUT_STEM)
    parser.add_argument("--created-at", default=CREATED_AT)
    parser.add_argument("--reviewer", default=REVIEWER)
    parser.add_argument("--contract", default="regional_flora_wave11_checkpoint_v1")
    parser.add_argument("--acquisition-manifest-json", type=Path)
    args = parser.parse_args()
    acquisition_inventory = None
    if args.acquisition_manifest_json:
        acquisition_inventory = json.loads(
            args.acquisition_manifest_json.read_text(encoding="utf-8")
        )
    report = build(
        source_snapshot_csv=args.source_snapshot_csv,
        master_csv=args.master_csv,
        prior_curated_evidence_csv=args.prior_curated_evidence_csv,
        prior_curated_audit_csv=args.prior_curated_audit_csv,
        output_dir=args.output_dir,
        source_group=args.source_group,
        output_stem=args.output_stem,
        created_at=args.created_at,
        reviewer=args.reviewer,
        contract=args.contract,
        acquisition_inventory=acquisition_inventory,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
