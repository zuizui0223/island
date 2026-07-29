"""Promote audited open-Web evidence through the shared all-evidence pipeline."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.open_web_common import (
    accepted_review_ids,
    rebuild_with_common_all_evidence,
    reviewed_audit_metrics,
)
from island_v2.open_web_network_pilot import build_priority_queue
from island_v2.open_web_search import (
    build_query_tasks,
    load_config,
    load_strict_name_variants,
    load_vernacular_names,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
        accepted_review_ids(audit, scope_table) if audit_summary["pilot_gate_passed"] else set()
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
    prior_formal = (
        pd.read_csv(prior_public_web_csv, dtype=str).fillna("")
        if prior_public_web_csv is not None
        else pd.DataFrame()
    )
    combined_formal = combine_public_web_ledgers(
        prior_formal,
        formal,
        prior_run_id=prior_public_web_run_id,
        prior_artifact=prior_public_web_artifact,
    )

    common_report, tables = rebuild_with_common_all_evidence(
        coverage=coverage,
        evidence=evidence,
        promoted=formal,
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
    promoted.to_csv(output_dir / "accepted_open_web_evidence.csv.gz", index=False)

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
            "source_file": (
                str(prior_public_web_csv) if prior_public_web_csv is not None else ""
            ),
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
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
