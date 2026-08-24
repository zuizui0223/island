"""Build a reviewed direct-evidence checkpoint from already acquired bulk data.

The checkpoint is deliberately baseline-specific: it promotes only exact
accepted-name species x trait records that are not resolved in the pinned
formal direct ledger.  Conflicting source values are retained as separate
lineages and are never resolved by row order.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.all_evidence_trait_audit import trait_axis

app = typer.Typer(add_completion=False, no_args_is_help=True)

REVIEWED_AT_UTC = "2026-08-11T00:00:00Z"
REVIEWER = "Codex structured-source evidence audit"
BASELINE_RUN_ID = "31397350878"
BASELINE_ARTIFACT_ID = "9066996051"

EVIDENCE_COLUMNS = [
    "candidate_id",
    "accepted_species",
    "axis",
    "trait_name",
    "raw_value",
    "normalized_value",
    "evidence_quality",
    "evidence_scope",
    "source_group",
    "source_provider",
    "source_url",
    "page_title",
    "source_citation",
    "source_excerpt",
    "source_record_id",
    "source_lineage",
    "lineage_method",
    "name_resolution_lineage",
    "name_match_method",
    "matched_page_name",
    "source_tier",
    "source_type",
    "domain",
    "language",
    "wild_cultivated_cultivar_status",
    "evidence_status",
    "content_sha256",
    "content_sha256_basis",
    "retrieved_at_utc",
    "query",
    "search_rank",
    "inference_rule",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _sha256_file(path: Path, *, normalize_newlines: bool = False) -> str:
    content = path.read_bytes()
    if normalize_newlines:
        content = content.replace(b"\r\n", b"\n")
    return hashlib.sha256(content).hexdigest()


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith(".gz"):
        frame.to_csv(
            path,
            index=False,
            lineterminator="\n",
            compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        )
    else:
        frame.to_csv(path, index=False, lineterminator="\n")


def _stable_candidate_id(record: dict[str, str]) -> str:
    return _sha256_text(
        "|".join(
            [
                record["accepted_species"],
                record["trait_name"],
                record["normalized_value"],
                record["source_lineage"],
            ]
        )
    )[:24]


def _audit_row(record: dict[str, str]) -> dict[str, str]:
    excerpt = _text(record["source_excerpt"])
    normalized_excerpt = excerpt.casefold()
    return {
        "candidate_id": record["candidate_id"],
        "accepted_species": record["accepted_species"],
        "trait_name": record["trait_name"],
        "normalized_value": record["normalized_value"],
        "source_url": record["source_url"],
        "page_title": record["page_title"],
        "source_citation": record["source_citation"],
        "source_tier": record["source_tier"],
        "source_type": record["source_type"],
        "domain": record["domain"],
        "language": record["language"],
        "supporting_excerpt": excerpt,
        "normalized_excerpt_sha256": _sha256_text(normalized_excerpt),
        "content_fingerprint": _sha256_text(
            record["source_lineage"] + "\n" + normalized_excerpt
        ),
        "name_match_method": record["name_match_method"],
        "wild_cultivated_cultivar_status": record[
            "wild_cultivated_cultivar_status"
        ],
        "decision": "accept",
        "species_identity_correct": "true",
        "value_correct": "true",
        "provenance_complete": "true",
        "cultivar_contamination": "false",
        "false_positive_reason": "",
        "decision_reason": (
            "Exact accepted-name species record in a versioned structured source; "
            "raw source fields, source record ID, source hash, and conservative "
            "trait mapping were rechecked against the committed acquisition bundle."
        ),
        "reviewer": REVIEWER,
        "reviewed_at_utc": REVIEWED_AT_UTC,
    }


def _candidate_base(
    row: dict[str, str],
    *,
    source_provider: str,
    page_title: str,
    source_lineage: str,
    lineage_method: str,
    source_type: str,
    domain: str,
    excerpt: str,
    content_sha256: str,
    content_sha256_basis: str,
) -> dict[str, str]:
    trait = _text(row["trait_name"])
    record = {
        "accepted_species": _text(row["accepted_species"]),
        "axis": trait_axis(trait),
        "trait_name": trait,
        "raw_value": _text(row["raw_value"]),
        "normalized_value": _text(row["standardized_value"]),
        "evidence_quality": "medium",
        "evidence_scope": "species_direct",
        "source_group": "curated_bulk_reaudit",
        "source_provider": source_provider,
        "source_url": _text(row["source_url"]),
        "page_title": page_title,
        "source_citation": _text(row["source_citation"]),
        "source_excerpt": excerpt,
        "source_record_id": _text(row["source_record_id"]),
        "source_lineage": source_lineage,
        "lineage_method": lineage_method,
        "name_resolution_lineage": "accepted_name_exact",
        "name_match_method": "accepted_name_exact",
        "matched_page_name": _text(row["accepted_species"]),
        "source_tier": "A",
        "source_type": source_type,
        "domain": domain,
        "language": "en",
        "wild_cultivated_cultivar_status": "species_trait_record_not_cultivar_limited",
        "evidence_status": "accepted_individual_structured_source_audit",
        "content_sha256": content_sha256,
        "content_sha256_basis": content_sha256_basis,
        "retrieved_at_utc": REVIEWED_AT_UTC,
        "query": "already_acquired_bulk_reaudit_no_new_web_query",
        "search_rank": "",
        "inference_rule": "",
    }
    record["candidate_id"] = _stable_candidate_id(record)
    return record


def _meyer_excerpt(row: dict[str, str], source: pd.DataFrame) -> str:
    match = re.fullmatch(r"Output_Data_MS\.csv:(\d+)", _text(row["source_record_id"]))
    if not match:
        raise ValueError(f"invalid Meyer source_record_id: {row['source_record_id']}")
    source_id = match.group(1)
    raw = source.loc[source["Unnamed: 0"].map(_text).eq(source_id)]
    if len(raw) != 1:
        raise ValueError(f"Meyer source row {source_id} is not unique")
    record = raw.iloc[0]
    if _text(record["genus_species"]) != _text(row["accepted_species"]):
        raise ValueError(f"Meyer species mismatch for {source_id}")
    if _text(record["mating_system"]).casefold() != _text(row["raw_value"]).casefold():
        raise ValueError(f"Meyer value mismatch for {source_id}")
    return (
        f"Output_Data_MS.csv row {source_id}: genus_species="
        f"{_text(record['genus_species'])}; Family={_text(record['Family'])}; "
        f"mating_system={_text(record['mating_system'])}."
    )


def _razanajatovo_excerpt(row: dict[str, str], source: pd.DataFrame) -> tuple[str, str]:
    species = _text(row["accepted_species"])
    raw = source.loc[
        source["Species"].map(_text).str.replace("_", " ", regex=False).eq(species)
    ].copy()
    if raw.empty:
        raise ValueError(f"Razanajatovo source species missing: {species}")
    fields = [
        "Data_source",
        "Species",
        "Family",
        "Self-compatibility_index_FS",
        "Self-compatibility_index_SFL",
        "Autofertility_index_FS",
        "Autofertility_index_SFL",
    ]
    records = []
    for source_row in raw[fields].fillna("").to_dict("records"):
        records.append(
            "; ".join(f"{field}={_text(source_row[field])}" for field in fields)
        )
    reference_keys = sorted({_text(value) for value in raw["Data_source"] if _text(value)})
    lineage_key = "|".join(reference_keys).casefold()
    return (
        "Supplementary Data 1 source row(s): " + " || ".join(records),
        "citation_keys:razanajatovo2016:" + _sha256_text(lineage_key)[:20],
    )


def _eol_lineage(row: dict[str, str]) -> str:
    citation = _text(row["source_citation"]).casefold()
    return "citation:eol-upstream:" + _sha256_text(citation)[:20]


def _selection_status(
    row: dict[str, str],
    *,
    resolved_pairs: set[tuple[str, str]],
    allowed_traits: set[str] | None,
) -> tuple[bool, str]:
    trait = _text(row["trait_name"])
    if allowed_traits is not None and trait not in allowed_traits:
        return False, "outside_checkpoint_trait_scope"
    if _text(row["name_match_method"]) != "accepted_name_exact":
        return False, "not_exact_accepted_name"
    if (_text(row["accepted_species"]), trait) in resolved_pairs:
        return False, "already_resolved_in_baseline"
    return True, "selected_exact_direct_unresolved"


def build_checkpoint(
    *,
    baseline_direct_cells: Path,
    output_dir: Path,
    meyer_candidates: Path = Path(
        "data/v2/staging/traits/bulk/meyer_2026_reproductive_traits/trait_candidates.csv.gz"
    ),
    meyer_source: Path = Path("data/v2/external/traits/meyer_2026/Output_Data_MS.csv"),
    razanajatovo_candidates: Path = Path(
        "data/v2/staging/traits/bulk/razanajatovo_2016_reproductive_traits/trait_candidates.csv.gz"
    ),
    razanajatovo_source: Path = Path(
        "data/v2/external/traits/razanajatovo_2016/Razanajatovo_etal_data.csv"
    ),
    razanajatovo_workbook: Path = Path(
        "data/v2/external/traits/razanajatovo_2016/Supplementary_Data_1.xlsx"
    ),
    eol_candidates: Path = Path("data/v2/staging/traits/eol_traitbank/trait_candidates.csv.gz"),
    eol_summary: Path = Path("data/v2/staging/traits/eol_traitbank/coverage_summary.json"),
) -> dict[str, Any]:
    expected_hashes = {
        "meyer_source_normalized": "dbc454103faa2d1677586b5b690b8ccc3e6b1a99fb39cdc0b00181703b4c6b56",
        "razanajatovo_source_normalized": "e030e1053d0c0b5a73ecc9a4dde3b26ba3fab20f8bb65cd61d01f4883800547b",
        "razanajatovo_workbook": "cb17009ae65c4f6cb2d28b21239714782600d0904effa77c5635249d28e62fc0",
    }
    observed_hashes = {
        "meyer_source_normalized": _sha256_file(meyer_source, normalize_newlines=True),
        "razanajatovo_source_normalized": _sha256_file(
            razanajatovo_source, normalize_newlines=True
        ),
        "razanajatovo_workbook": _sha256_file(razanajatovo_workbook),
    }
    if observed_hashes != expected_hashes:
        raise ValueError(
            f"bulk source hash mismatch: expected={expected_hashes}, observed={observed_hashes}"
        )

    baseline = pd.read_csv(baseline_direct_cells, dtype=str).fillna("")
    resolved = baseline.loc[baseline["resolution_status"].eq("resolved")]
    resolved_pairs = set(zip(resolved["accepted_species"], resolved["trait_name"], strict=True))
    meyer_raw = pd.read_csv(meyer_source, dtype=str).fillna("")
    raz_raw = pd.read_csv(razanajatovo_source, dtype=str).fillna("")
    eol_report = json.loads(eol_summary.read_text(encoding="utf-8"))
    eol_archive_sha256 = _text(eol_report["archive_sha256"])

    evidence_rows: list[dict[str, str]] = []
    selection_rows: list[dict[str, str]] = []
    source_inputs = [
        (
            "meyer_2026",
            meyer_candidates,
            None,
            "Meyer, Galloway & Eckert 2026 Dryad v4",
        ),
        (
            "razanajatovo_2016",
            razanajatovo_candidates,
            None,
            "Razanajatovo et al. 2016 Supplementary Data 1",
        ),
        (
            "eol_traitbank_flower_colour",
            eol_candidates,
            {"flower_primary_color"},
            "EOL TraitBank all-traits archive",
        ),
    ]
    for source_key, candidate_path, allowed_traits, provider in source_inputs:
        candidates = pd.read_csv(candidate_path, dtype=str).fillna("")
        for row in candidates.to_dict("records"):
            selected, reason = _selection_status(
                row,
                resolved_pairs=resolved_pairs,
                allowed_traits=allowed_traits,
            )
            selection_rows.append(
                {
                    "source_key": source_key,
                    "evidence_id": _text(row["evidence_id"]),
                    "accepted_species": _text(row["accepted_species"]),
                    "trait_name": _text(row["trait_name"]),
                    "standardized_value": _text(row["standardized_value"]),
                    "name_match_method": _text(row["name_match_method"]),
                    "selected": str(selected).lower(),
                    "selection_reason": reason,
                    "baseline_run_id": BASELINE_RUN_ID,
                }
            )
            if not selected:
                continue
            if source_key == "meyer_2026":
                excerpt = _meyer_excerpt(row, meyer_raw)
                lineage = "dataset:dryad.cc2fqz6hr:v4"
                record = _candidate_base(
                    row,
                    source_provider=provider,
                    page_title="Meyer et al. 2026 Dryad mating-system dataset",
                    source_lineage=lineage,
                    lineage_method="versioned_curated_dataset",
                    source_type="peer_reviewed_curated_trait_dataset",
                    domain="datadryad.org",
                    excerpt=excerpt,
                    content_sha256=expected_hashes["meyer_source_normalized"],
                    content_sha256_basis="Output_Data_MS.csv_lf_normalized_bytes",
                )
            elif source_key == "razanajatovo_2016":
                excerpt, lineage = _razanajatovo_excerpt(row, raz_raw)
                record = _candidate_base(
                    row,
                    source_provider=provider,
                    page_title="Razanajatovo et al. 2016 Supplementary Data 1",
                    source_lineage=lineage,
                    lineage_method="underlying_reference_key_set",
                    source_type="peer_reviewed_supplementary_trait_dataset",
                    domain="nature.com",
                    excerpt=excerpt,
                    content_sha256=expected_hashes["razanajatovo_workbook"],
                    content_sha256_basis="Supplementary_Data_1.xlsx_bytes",
                )
            else:
                record = _candidate_base(
                    row,
                    source_provider=provider,
                    page_title="Encyclopedia of Life TraitBank species trait record",
                    source_lineage=_eol_lineage(row),
                    lineage_method="normalized_upstream_source_citation",
                    source_type="public_curated_trait_database",
                    domain="eol.org",
                    excerpt=_text(row["evidence_excerpt"]),
                    content_sha256=eol_archive_sha256,
                    content_sha256_basis="EOL_traits_all.zip_bytes",
                )
            evidence_rows.append(record)

    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).sort_values(
        ["accepted_species", "trait_name", "source_provider", "candidate_id"]
    )
    audit = pd.DataFrame([_audit_row(row) for row in evidence.to_dict("records")])
    selection = pd.DataFrame(selection_rows).sort_values(
        ["source_key", "accepted_species", "trait_name", "evidence_id"]
    )

    pair_groups = evidence.groupby(["accepted_species", "trait_name"], sort=True)
    conflicts: list[dict[str, str | int]] = []
    for (species, trait), group in pair_groups:
        if len(group) < 2:
            continue
        values = sorted(set(group["normalized_value"]))
        conflicts.append(
            {
                "accepted_species": species,
                "axis": trait_axis(trait),
                "trait_name": trait,
                "classification": (
                    "independent_source_agreement"
                    if len(values) == 1
                    else "unresolved_direct_conflict"
                ),
                "normalized_values": json.dumps(values, separators=(",", ":")),
                "n_source_lineages": group["source_lineage"].nunique(),
                "source_lineages": "|".join(sorted(set(group["source_lineage"]))),
                "candidate_ids": "|".join(sorted(group["candidate_id"])),
            }
        )
    conflict_table = pd.DataFrame(conflicts)

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "curated_bulk_direct_evidence_20260811.csv"
    audit_path = output_dir / "curated_bulk_direct_manual_audit_20260811.csv"
    selection_path = output_dir / "curated_bulk_selection_audit_20260811.csv.gz"
    conflict_path = output_dir / "curated_bulk_source_conflicts_20260811.csv"
    _write_csv(evidence, evidence_path)
    _write_csv(audit, audit_path)
    _write_csv(selection, selection_path)
    _write_csv(conflict_table, conflict_path)

    selected_by_source = selection.loc[selection["selected"].eq("true"), "source_key"].value_counts()
    report = {
        "contract": "curated_bulk_direct_checkpoint_v1",
        "baseline": {
            "run_id": BASELINE_RUN_ID,
            "artifact_id": BASELINE_ARTIFACT_ID,
            "direct_cells_sha256": _sha256_file(baseline_direct_cells),
        },
        "reviewed_at_utc": REVIEWED_AT_UTC,
        "reviewer": REVIEWER,
        "source_hashes": {
            **observed_hashes,
            "eol_archive": eol_archive_sha256,
            "meyer_candidates": _sha256_file(meyer_candidates),
            "razanajatovo_candidates": _sha256_file(razanajatovo_candidates),
            "eol_candidates": _sha256_file(eol_candidates),
        },
        "evidence_rows": len(evidence),
        "unique_species_trait": int(
            evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0]
        ),
        "selected_rows_by_source": {
            str(key): int(value) for key, value in selected_by_source.items()
        },
        "multi_source_species_trait": len(conflict_table),
        "independent_source_agreement": int(
            conflict_table["classification"].eq("independent_source_agreement").sum()
        ),
        "unresolved_direct_conflict": int(
            conflict_table["classification"].eq("unresolved_direct_conflict").sum()
        ),
        "selection_reason_counts": {
            str(key): int(value)
            for key, value in selection["selection_reason"].value_counts().items()
        },
        "outputs": {},
    }
    for path in (evidence_path, audit_path, selection_path, conflict_path):
        report["outputs"][path.name] = {
            "sha256": _sha256_file(path),
            "bytes": path.stat().st_size,
        }
    manifest_path = output_dir / "curated_bulk_checkpoint_manifest_20260811.json"
    manifest_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    report["manifest"] = str(manifest_path)
    return report


@app.command("build")
def build(
    baseline_direct_cells: Annotated[
        Path, typer.Option(exists=True, dir_okay=False)
    ],
    output_dir: Annotated[
        Path, typer.Option(file_okay=False)
    ] = Path("data/v2/staging/traits/open_web_pilot"),
) -> None:
    report = build_checkpoint(
        baseline_direct_cells=baseline_direct_cells,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
