"""Promote audited open-Web evidence through the shared all-evidence pipeline."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.all_evidence_trait_audit import load_ontology, normalise_state_set
from island_v2.open_web_common import (
    production_candidate_ids,
    promoted_to_common_evidence,
    rebuild_with_common_all_evidence,
    reviewed_audit_metrics,
    reviewed_source_package_evidence,
)
from island_v2.open_web_network_pilot import build_priority_queue
from island_v2.open_web_search import (
    build_query_tasks,
    load_config,
    load_strict_name_variants,
    load_vernacular_names,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

_TRUE_VALUES = {"1", "true", "yes", "y"}
_DIRECT_NAME_MATCHES = {
    "accepted_name_exact",
    "exact_accepted_name",
    "exact_synonym",
    "synonym_exact",
}
_DIRECT_SCOPES = {"species_direct", "synonym_direct"}
_MIN_SOURCE_PACKAGE_NOVELTY = 0.50


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bool(value: object) -> bool:
    return str(value).strip().casefold() in _TRUE_VALUES


def validate_individually_reviewed_evidence(
    candidates: pd.DataFrame,
    audit: pd.DataFrame,
    master: pd.DataFrame,
    ontology_yaml: Path,
) -> pd.DataFrame:
    """Fail closed on small, line-by-line reviewed species-direct batches.

    The domain x trait production gate is for automated inventory/search lanes.
    A separately committed curated batch may bypass that aggregate gate only
    when every candidate has an explicit, correct manual review and complete
    source-backed provenance.
    """

    candidate_required = {
        "candidate_id",
        "accepted_species",
        "trait_name",
        "normalized_value",
        "evidence_quality",
        "evidence_scope",
        "source_url",
        "source_citation",
        "source_excerpt",
        "source_lineage",
        "name_match_method",
        "source_tier",
        "content_sha256",
        "content_sha256_basis",
        "retrieved_at_utc",
        "inference_rule",
    }
    audit_required = {
        "candidate_id",
        "decision",
        "species_identity_correct",
        "value_correct",
        "provenance_complete",
        "cultivar_contamination",
        "reviewer",
        "reviewed_at_utc",
        "supporting_excerpt",
        "normalized_excerpt_sha256",
        "content_fingerprint",
    }
    missing_candidates = candidate_required.difference(candidates.columns)
    missing_audit = audit_required.difference(audit.columns)
    if missing_candidates:
        raise ValueError(
            f"individually reviewed evidence is missing columns: {sorted(missing_candidates)}"
        )
    if missing_audit:
        raise ValueError(f"individual audit is missing columns: {sorted(missing_audit)}")
    if candidates.empty:
        return candidates.copy()
    if candidates["candidate_id"].duplicated().any() or audit["candidate_id"].duplicated().any():
        raise ValueError("individual candidate_id values must be unique")
    candidate_ids = set(candidates["candidate_id"].astype(str))
    audit_ids = set(audit["candidate_id"].astype(str))
    if candidate_ids != audit_ids:
        raise ValueError("individual audit must cover every candidate exactly once")

    reviewed = candidates.merge(
        audit,
        on="candidate_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_audit"),
    ).fillna("")
    correct = (
        reviewed["decision"].str.casefold().eq("accept")
        & reviewed["species_identity_correct"].map(_bool)
        & reviewed["value_correct"].map(_bool)
        & reviewed["provenance_complete"].map(_bool)
        & ~reviewed["cultivar_contamination"].map(_bool)
        & reviewed["reviewer"].astype(str).str.strip().ne("")
        & reviewed["reviewed_at_utc"].astype(str).str.strip().ne("")
    )
    if not bool(correct.all()):
        rejected = reviewed.loc[~correct, "candidate_id"].astype(str).tolist()
        raise ValueError(f"individually curated evidence must be accepted and correct: {rejected}")

    accepted_names = set(master["accepted_species"].astype(str))
    ontology = load_ontology(ontology_yaml)
    failures: list[str] = []
    for row in reviewed.to_dict("records"):
        candidate_id = str(row["candidate_id"])
        trait = str(row["trait_name"])
        valid, invalid = normalise_state_set(trait, row["normalized_value"], ontology)
        excerpt = " ".join(str(row["source_excerpt"]).split())
        normalized_excerpt = excerpt.casefold()
        expected_excerpt_sha256 = hashlib.sha256(normalized_excerpt.encode("utf-8")).hexdigest()
        expected_fingerprint = hashlib.sha256(
            (str(row["source_lineage"]) + "\n" + normalized_excerpt).encode("utf-8")
        ).hexdigest()
        expected_candidate_id = hashlib.sha256(
            "|".join(
                [
                    str(row["accepted_species"]),
                    trait,
                    str(row["normalized_value"]),
                    str(row["source_lineage"]),
                ]
            ).encode("utf-8")
        ).hexdigest()[:24]
        checks = {
            "unstable_candidate_id": candidate_id != expected_candidate_id,
            "species_not_in_master": str(row["accepted_species"]) not in accepted_names,
            "invalid_direct_scope": str(row["evidence_scope"]) not in _DIRECT_SCOPES,
            "invalid_name_match": str(row["name_match_method"]) not in _DIRECT_NAME_MATCHES,
            "invalid_source_tier": str(row["source_tier"]).upper() not in {"A", "B"},
            "invalid_quality": str(row["evidence_quality"]).casefold() not in {"high", "medium"},
            "high_without_tier_a": (
                str(row["evidence_quality"]).casefold() == "high"
                and str(row["source_tier"]).upper() != "A"
            ),
            "invalid_ontology_value": bool(invalid) or not valid,
            "invalid_source_url": not str(row["source_url"]).startswith(("http://", "https://")),
            "missing_source_citation": not str(row["source_citation"]).strip(),
            "missing_exact_excerpt": not excerpt,
            "audit_excerpt_mismatch": excerpt != " ".join(str(row["supporting_excerpt"]).split()),
            "invalid_excerpt_sha256": str(row["normalized_excerpt_sha256"])
            != expected_excerpt_sha256,
            "invalid_content_fingerprint": str(row["content_fingerprint"]) != expected_fingerprint,
            "missing_source_lineage": not str(row["source_lineage"]).strip(),
            "invalid_content_sha256": re.fullmatch(
                r"[0-9a-f]{64}", str(row["content_sha256"]).strip().casefold()
            )
            is None,
            "missing_content_sha256_basis": not str(row["content_sha256_basis"]).strip(),
            "invalid_retrieval_date": re.fullmatch(
                r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z",
                str(row["retrieved_at_utc"]).strip(),
            )
            is None,
            "contains_inference_rule": bool(str(row["inference_rule"]).strip()),
        }
        failures.extend(f"{candidate_id}:{reason}" for reason, failed in checks.items() if failed)
    duplicate_keys = reviewed.duplicated(
        ["accepted_species", "trait_name", "normalized_value", "source_lineage"]
    )
    if duplicate_keys.any():
        failures.extend(
            f"{value}:duplicate_species_trait_value_lineage"
            for value in reviewed.loc[duplicate_keys, "candidate_id"].astype(str)
        )
    if failures:
        raise ValueError("individual evidence validation failed: " + "; ".join(failures))
    return candidates.copy().fillna("")


def _load_baseline(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage = pd.read_csv(
        root / "integrated_species_axis_coverage.csv.gz",
        dtype=str,
    ).fillna("")
    evidence = pd.read_csv(
        root / "integrated_evidence_lineage.csv.gz",
        dtype=str,
    ).fillna("")
    if (
        len(coverage) != 318_885
        or coverage["accepted_species"].nunique() != 106_295
        or coverage.duplicated(["accepted_species", "axis"]).any()
    ):
        raise ValueError("baseline does not satisfy the fixed 106295 x 3 contract")
    return coverage, evidence


def combine_public_web_ledgers(
    prior: pd.DataFrame,
    promoted: pd.DataFrame,
    *,
    prior_run_id: str = "",
    prior_artifact: str = "",
) -> pd.DataFrame:
    """Append reviewed evidence without replacing the prior formal Web ledger."""

    old = prior.copy().fillna("")
    new = promoted.copy().fillna("")
    for frame in (old, new):
        for column in (
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_lineage",
            "source_url",
            "source_run_id",
            "source_artifact",
        ):
            if column not in frame:
                frame[column] = ""
    if not old.empty:
        old["source_run_id"] = old["source_run_id"].where(
            old["source_run_id"].ne(""),
            prior_run_id,
        )
        old["source_artifact"] = old["source_artifact"].where(
            old["source_artifact"].ne(""),
            prior_artifact,
        )
    combined = pd.concat([old, new], ignore_index=True, sort=False).fillna("")
    if combined.empty:
        return combined
    combined["_dedup_lineage"] = combined["source_lineage"].where(
        combined["source_lineage"].ne(""),
        "url:" + combined["source_url"].str.rstrip("/").str.casefold(),
    )
    combined = combined.drop_duplicates(
        [
            "accepted_species",
            "trait_name",
            "_dedup_lineage",
            "normalized_value",
        ],
        keep="last",
    )
    return combined.drop(columns="_dedup_lineage").reset_index(drop=True)


def _report_groups(reviewed: pd.DataFrame, column: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for value, group in reviewed.groupby(column, sort=True):
        accepted_correct = int(group["_accepted_correct"].sum())
        rows.append(
            {
                column: str(value) or "unknown_not_recorded",
                "reviewed": len(group),
                "accepted_correct": accepted_correct,
                "precision": accepted_correct / len(group),
                "cultivar_contamination_rate": float(group["_cultivar"].mean()),
            }
        )
    return rows


def finalize_review(
    *,
    baseline_dir: Path,
    candidates_csv: Path,
    audit_csv: Path,
    master_csv: Path,
    output_dir: Path,
    ontology_yaml: Path = Path("config/trait_ontology.yml"),
    source_run_id: str = "",
    source_artifact: str = "",
    prior_public_web_csv: Path | None = None,
    prior_public_web_run_id: str = "",
    prior_public_web_artifact: str = "",
    synonym_csv: Path | None = None,
    vernacular_csv: Path | None = None,
    acquisition_config_yaml: Path = Path("config/open_web_trait_acquisition.yml"),
    curated_evidence_csv: Path | None = None,
    curated_audit_csv: Path | None = None,
    curated_source_run_id: str = "",
    curated_source_artifact: str = "",
    source_package_evidence_csv: Path | None = None,
    source_package_audit_csv: Path | None = None,
    source_package_run_id: str = "",
    source_package_artifact: str = "",
) -> dict[str, Any]:
    """Apply the audit gate and run PR #131's shared all-evidence functions."""

    coverage, evidence = _load_baseline(baseline_dir)
    candidates = pd.read_csv(candidates_csv, dtype=str).fillna("")
    audit = pd.read_csv(audit_csv, dtype=str).fillna("")
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    scope_table, audit_summary = reviewed_audit_metrics(audit)

    reviewed = audit.loc[audit["decision"].str.casefold().isin({"accept", "reject"})].copy()
    reviewed["_accepted_correct"] = (
        reviewed["decision"].str.casefold().eq("accept")
        & reviewed["species_identity_correct"]
        .astype(str)
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
        & reviewed["value_correct"].astype(str).str.casefold().isin({"1", "true", "yes", "y"})
        & reviewed["provenance_complete"].astype(str).str.casefold().isin({"1", "true", "yes", "y"})
        & ~reviewed["cultivar_contamination"]
        .astype(str)
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
    )
    reviewed["_cultivar"] = (
        reviewed["cultivar_contamination"]
        .astype(str)
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
    )
    approved_ids = (
        production_candidate_ids(candidates, audit, scope_table)
        if audit_summary["pilot_gate_passed"]
        else set()
    )
    promoted = candidates.loc[candidates["candidate_id"].isin(approved_ids)].copy()
    if not promoted.empty:
        promoted["review_status"] = "accepted_manual_audit"
        promoted["source_run_id"] = source_run_id
        promoted["source_artifact"] = source_artifact
        promoted["source_file"] = str(candidates_csv)

    formal = promoted.copy()
    if not formal.empty:
        formal["evidence_quality"] = formal.get("confidence", "medium")
        formal["source_provider"] = formal.get("domain", "")
        formal["source_record_id"] = formal.get("candidate_id", "")
        formal["source_citation"] = formal.get("source_citation", formal.get("page_title", ""))
        formal["source_excerpt"] = formal.get(
            "source_excerpt",
            formal.get("supporting_excerpt", ""),
        )
        formal["inference_rule"] = ""
    curated_formal = pd.DataFrame()
    curated_audit_summary: dict[str, Any] = {
        "reviewed": 0,
        "accepted_correct": 0,
        "precision": None,
        "cultivar_contamination_rate": None,
    }
    if (curated_evidence_csv is None) != (curated_audit_csv is None):
        raise ValueError("curated_evidence_csv and curated_audit_csv must be supplied together")
    if curated_evidence_csv is not None and curated_audit_csv is not None:
        curated_candidates = pd.read_csv(curated_evidence_csv, dtype=str).fillna("")
        curated_audit = pd.read_csv(curated_audit_csv, dtype=str).fillna("")
        curated_formal = validate_individually_reviewed_evidence(
            curated_candidates,
            curated_audit,
            master,
            ontology_yaml,
        )
        curated_formal["review_status"] = "accepted_individual_manual_audit"
        curated_formal["evidence_status"] = "accepted_manual_audit"
        curated_formal["source_run_id"] = curated_source_run_id
        curated_formal["source_artifact"] = curated_source_artifact or str(curated_evidence_csv)
        curated_formal["source_file"] = str(curated_evidence_csv)
        curated_audit_summary = {
            "reviewed": len(curated_audit),
            "accepted_correct": len(curated_formal),
            "precision": len(curated_formal) / len(curated_audit),
            "cultivar_contamination_rate": 0.0,
        }
    web_incremental_formal = pd.concat(
        [formal, curated_formal],
        ignore_index=True,
        sort=False,
    ).fillna("")
    prior_formal = (
        pd.read_csv(prior_public_web_csv, dtype=str).fillna("")
        if prior_public_web_csv is not None
        else pd.DataFrame()
    )
    source_package_formal = pd.DataFrame()
    source_package_scopes = pd.DataFrame()
    source_package_summary: dict[str, Any] = {
        "reviewed": 0,
        "accepted_correct": 0,
        "precision": None,
        "cultivar_contamination_rate": None,
        "package_gate_passed": False,
        "approved_traits": [],
        "candidate_evidence_rows": 0,
        "selected_evidence_rows": 0,
        "unique_species_trait": 0,
        "novel_species_trait": 0,
        "novelty_rate": None,
        "minimum_novelty_rate": _MIN_SOURCE_PACKAGE_NOVELTY,
    }
    if (source_package_evidence_csv is None) != (source_package_audit_csv is None):
        raise ValueError(
            "source_package_evidence_csv and source_package_audit_csv must be supplied together"
        )
    if source_package_evidence_csv is not None and source_package_audit_csv is not None:
        source_package = pd.read_csv(source_package_evidence_csv, dtype=str).fillna("")
        source_package_audit = pd.read_csv(source_package_audit_csv, dtype=str).fillna("")
        (
            source_package_formal,
            source_package_scopes,
            source_package_summary,
        ) = reviewed_source_package_evidence(source_package, source_package_audit)
        if not source_package_summary["package_gate_passed"]:
            raise ValueError("source package did not pass its review gate")
        if source_package_run_id:
            source_package_formal["source_run_id"] = source_package_formal[
                "source_run_id"
            ].where(source_package_formal["source_run_id"].ne(""), source_package_run_id)
        if source_package_artifact:
            source_package_formal["source_artifact"] = source_package_formal[
                "source_artifact"
            ].where(
                source_package_formal["source_artifact"].ne(""),
                source_package_artifact,
            )
        source_package_formal["source_file"] = source_package_formal["source_file"].where(
            source_package_formal["source_file"].ne(""),
            str(source_package_evidence_csv),
        )

        current_direct = pd.concat(
            [
                evidence,
                promoted_to_common_evidence(prior_formal),
                promoted_to_common_evidence(web_incremental_formal),
            ],
            ignore_index=True,
            sort=False,
        ).fillna("")
        current_direct = current_direct.loc[
            current_direct["quality"].str.casefold().isin({"high", "medium"})
            & current_direct["evidence_scope"].str.casefold().isin(_DIRECT_SCOPES)
        ]
        current_pairs = set(
            zip(
                current_direct["accepted_species"],
                current_direct["trait_name"],
                strict=True,
            )
        )
        package_pairs = set(
            zip(
                source_package_formal["accepted_species"],
                source_package_formal["trait_name"],
                strict=True,
            )
        )
        novel_pairs = package_pairs.difference(current_pairs)
        novelty_rate = len(novel_pairs) / len(package_pairs) if package_pairs else 0.0
        source_package_summary.update(
            {
                "unique_species_trait": len(package_pairs),
                "novel_species_trait": len(novel_pairs),
                "novelty_rate": novelty_rate,
                "minimum_novelty_rate": _MIN_SOURCE_PACKAGE_NOVELTY,
            }
        )
        if novelty_rate < _MIN_SOURCE_PACKAGE_NOVELTY:
            raise ValueError(
                f"source-package novelty {novelty_rate:.6f} is below "
                f"{_MIN_SOURCE_PACKAGE_NOVELTY:.2f}"
            )

    incremental_formal = pd.concat(
        [web_incremental_formal, source_package_formal],
        ignore_index=True,
        sort=False,
    ).fillna("")
    combined_formal = combine_public_web_ledgers(
        prior_formal,
        incremental_formal,
        prior_run_id=prior_public_web_run_id,
        prior_artifact=prior_public_web_artifact,
    )

    common_report, tables = rebuild_with_common_all_evidence(
        coverage=coverage,
        evidence=evidence,
        promoted=incremental_formal,
        prior_promoted=prior_formal,
        master=master,
        ontology_yaml=ontology_yaml,
    )
    strict_for_queue = tables["strict_coverage"].rename(columns={"quality": "after_quality"})
    species_queue = build_priority_queue(
        strict_for_queue,
        tables["lineages"],
        master,
        species_limit=1000,
        traits_per_species=2,
    )
    search_tasks = build_query_tasks(
        species_queue.to_dict("records"),
        config=load_config(acquisition_config_yaml),
        synonyms=load_strict_name_variants(synonym_csv),
        vernacular_names=load_vernacular_names(vernacular_csv),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    scope_table.to_csv(
        output_dir / "domain_trait_audit_precision.csv",
        index=False,
    )
    reviewed.drop(columns=["_accepted_correct", "_cultivar"]).to_csv(
        output_dir / "reviewed_audit_sample.csv",
        index=False,
    )
    incremental_formal.to_csv(
        output_dir / "accepted_open_web_evidence.csv.gz",
        index=False,
    )
    source_package_formal.to_csv(
        output_dir / "accepted_source_package_evidence.csv.gz",
        index=False,
    )
    source_package_scopes.to_csv(
        output_dir / "source_package_trait_audit_precision.csv",
        index=False,
    )

    combined_formal.to_csv(
        output_dir / "broad_web_medium_evidence.csv.gz",
        index=False,
    )

    table_files = {
        "lineages": "reviewed_source_lineages.csv.gz",
        "lineage_duplicates": "source_lineage_duplicates.csv.gz",
        "direct_cells": "direct_species_trait_ledger.csv.gz",
        "direct_conflicts": "source_lineage_conflict_table.csv.gz",
        "rules": "trait_specific_genus_rule_audit.csv.gz",
        "prior_rules": "prior_trait_specific_genus_rule_audit.csv.gz",
        "rebuilt_low": "rebuilt_all_evidence_validated_low.csv.gz",
        "prior_rebuilt_low": "prior_rebuilt_all_evidence_validated_low.csv.gz",
        "pilot_low_comparison": "pilot_validated_low_change.csv.gz",
        "old_low_comparison": "old_low_comparison.csv.gz",
        "strict_coverage": "strict_species_axis_coverage.csv.gz",
        "next_queue": "prioritized_unresolved_search_queue.csv.gz",
    }
    for key, filename in table_files.items():
        tables[key].to_csv(output_dir / filename, index=False)
    species_queue.to_csv(
        output_dir / "prioritized_search_queue.csv.gz",
        index=False,
    )
    pd.DataFrame([asdict(task) for task in search_tasks]).to_csv(
        output_dir / "search_tasks.csv", index=False
    )

    report = {
        "contract": "reviewed_open_web_formal_ledger_v1",
        "completed_at_utc": _now(),
        "source_run_id": source_run_id,
        "source_artifact": source_artifact,
        "search_queries": 0,
        "search_cost_usd": 0.0,
        "credential_free_inventory": True,
        "audit": {
            **audit_summary,
            "by_domain": _report_groups(reviewed, "domain"),
            "by_language": _report_groups(reviewed, "language"),
            "by_trait": _report_groups(reviewed, "trait_name"),
        },
        "production_approved_domain_trait_scopes": (
            scope_table.loc[
                scope_table["production_approved"].astype(bool),
                ["domain", "trait_name"],
            ].to_dict("records")
            if not scope_table.empty
            else []
        ),
        "promoted_reviewed_evidence_rows": len(promoted),
        "individually_reviewed_evidence_rows": len(curated_formal),
        "total_incremental_reviewed_evidence_rows": len(incremental_formal),
        "individual_manual_audit": curated_audit_summary,
        "source_package": {
            **source_package_summary,
            "source_run_id": source_package_run_id,
            "source_artifact": source_package_artifact,
        },
        "prior_formal_public_web_evidence_rows": len(prior_formal),
        "combined_formal_public_web_evidence_rows": len(combined_formal),
        "next_search_queue": {
            "species_trait_rows": len(species_queue),
            "search_tasks": len(search_tasks),
            "strict_synonym_snapshot_used": synonym_csv is not None,
            "vernacular_snapshot_used": vernacular_csv is not None,
        },
        "shared_all_evidence": common_report,
    }
    (output_dir / "open_web_formalization_summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "broad_web_medium_pilot_summary.json").write_text(
        json.dumps(
            {
                "contract": "broad_web_medium_pilot_summary_v1",
                "run_id": source_run_id,
                "artifact": source_artifact,
                "accepted": {
                    "evidence_rows": len(combined_formal),
                },
                "accepted_evidence_rows": len(combined_formal),
                "new_promoted_evidence_rows": len(promoted),
                "prior_evidence_rows": len(prior_formal),
                "manual_reviewed": True,
                "audit_precision": audit_summary["precision"],
                "cultivar_contamination_rate": audit_summary["cultivar_contamination_rate"],
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    files = sorted(path for path in output_dir.rglob("*") if path.is_file())
    manifest = {
        "contract": "reviewed_open_web_source_run_manifest_v1",
        "source_run_id": source_run_id,
        "source_artifact": source_artifact,
        "common_validated_low_implementation": ("island_v2.all_evidence_trait_audit"),
        "prior_public_web": {
            "run_id": prior_public_web_run_id,
            "artifact": prior_public_web_artifact,
            "source_file": (str(prior_public_web_csv) if prior_public_web_csv is not None else ""),
        },
        "curated_inputs": [
            {
                "path": str(path),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in (curated_evidence_csv, curated_audit_csv)
            if path is not None
        ],
        "source_package": {
            "run_id": source_package_run_id,
            "artifact": source_package_artifact,
            "inputs": [
                {
                    "path": str(path),
                    "sha256": _sha256(path),
                    "size_bytes": path.stat().st_size,
                }
                for path in (source_package_evidence_csv, source_package_audit_csv)
                if path is not None
            ],
        },
        "files": [
            {
                "path": str(path.relative_to(output_dir)),
                "sha256": _sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in files
        ],
    }
    (output_dir / "source_run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("apply-audit")
def apply_audit_command(
    baseline_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    audit_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    master_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    ontology_yaml: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ] = Path("config/trait_ontology.yml"),
    source_run_id: str = "",
    source_artifact: str = "",
    prior_public_web_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    prior_public_web_run_id: str = "",
    prior_public_web_artifact: str = "",
    synonym_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    vernacular_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    acquisition_config_yaml: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ] = Path("config/open_web_trait_acquisition.yml"),
    curated_evidence_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    curated_audit_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    curated_source_run_id: str = "",
    curated_source_artifact: str = "",
    source_package_evidence_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    source_package_audit_csv: Annotated[
        Path | None,
        typer.Option(exists=True, dir_okay=False),
    ] = None,
    source_package_run_id: str = "",
    source_package_artifact: str = "",
) -> None:
    report = finalize_review(
        baseline_dir=baseline_dir,
        candidates_csv=candidates_csv,
        audit_csv=audit_csv,
        master_csv=master_csv,
        output_dir=output_dir,
        ontology_yaml=ontology_yaml,
        source_run_id=source_run_id,
        source_artifact=source_artifact,
        prior_public_web_csv=prior_public_web_csv,
        prior_public_web_run_id=prior_public_web_run_id,
        prior_public_web_artifact=prior_public_web_artifact,
        synonym_csv=synonym_csv,
        vernacular_csv=vernacular_csv,
        acquisition_config_yaml=acquisition_config_yaml,
        curated_evidence_csv=curated_evidence_csv,
        curated_audit_csv=curated_audit_csv,
        curated_source_run_id=curated_source_run_id,
        curated_source_artifact=curated_source_artifact,
        source_package_evidence_csv=source_package_evidence_csv,
        source_package_audit_csv=source_package_audit_csv,
        source_package_run_id=source_package_run_id,
        source_package_artifact=source_package_artifact,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
