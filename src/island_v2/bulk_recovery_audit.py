"""Recover evidence the pipeline already fetched and then discarded.

Two independent lanes, measured against the rejection ledgers the pipeline
itself wrote. Neither lane invents evidence and neither promotes anything: each
separates "rejected for a good reason" from "rejected by a gap in our own
matching or extraction", and emits a queue for the second group.

**Name matching.** Bulk joins reject a large share of source rows as
``reject_not_exact_synonym`` even where GBIF returned an EXACT match on an
ACCEPTED name. Those rows all carry a GBIF usage key, so the accepted name is
recoverable by key lookup rather than by string matching. This module classifies
which rejections qualify and emits the key queue; the lookup itself needs
network and runs in CI.

**Extraction.** Source pages were fetched, quoted, and then rejected by the
context gate. Every rejected candidate retains its exact supporting quote, so
the classification is reproducible offline. Reading the quotes rather than the
reason labels shows three distinct populations: measurements that belong to a
non-target organ (correctly rejected -- a style is not a flower), descriptive
templates that enumerate every state instead of asserting one (correctly
rejected), and plain statements in a language the English-centric rules never
handled (recoverable).

Nothing here writes accepted evidence. The output is an audited queue.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

NAME_MATCH_REQUIRED = {"match_status", "gbif_match_type", "gbif_status", "reason"}
EXTRACTION_REQUIRED = {"exact_supporting_quote", "context_gate_reason", "trait_name"}

# name-match classes
RECOVERABLE_KEY = "recoverable_by_gbif_usage_key"
BLOCKED_GBIF = "blocked_gbif_not_confident"
BLOCKED_REASON = "blocked_rejection_reason_is_final"
BLOCKED_NO_KEY = "blocked_no_usage_key"

# extraction classes
CORRECT_ORGAN = "correctly_rejected_non_target_organ"
CORRECT_TEMPLATE = "correctly_rejected_enumeration_template"
RECOVERABLE_LANG = "recoverable_non_english_statement"
NEEDS_REVIEW = "needs_manual_review"


@app.callback()
def main() -> None:
    """Audit and queue evidence discarded by name matching or by extraction."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned recovery configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"name_match", "extraction", "language_aliases"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "config must contain name_match, extraction and language_aliases"
        )
    return config


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).strip().split())


def normalize_language(raw: object, aliases: dict[str, Any]) -> str:
    """Map a free-form language tag onto a code, or "unknown".

    The acquired ledgers carry en, English, EN, eng, fre, es, Spanish and empty
    in the same column. An unresolvable tag becomes ``unknown`` and never
    defaults to English, because defaulting would silently re-apply the
    English-only rules that caused the loss in the first place.
    """
    text = _text(raw).lower()
    if not text:
        return "unknown"
    for code, variants in aliases.items():
        if text == code or text in {str(v).lower() for v in variants}:
            return code
    return "unknown"


def classify_name_match(audit: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Split rejected bulk rows into key-resolvable and finally-rejected."""
    missing = NAME_MATCH_REQUIRED - set(audit.columns)
    if missing:
        raise typer.BadParameter(f"name-match audit is missing columns: {sorted(missing)}")

    rules = config["name_match"]
    allowed_type = {str(v).upper() for v in rules["required_gbif_match_type"]}
    allowed_status = {str(v).upper() for v in rules["required_gbif_status"]}
    final_reasons = {str(v) for v in (rules.get("never_recover_reasons") or ())}

    frame = audit.loc[audit["match_status"].map(_text).eq("unmatched")].copy()
    if frame.empty:
        return frame.assign(recovery_class="", gbif_key="")

    key_columns = [c for c in ("gbif_accepted_usage_key", "gbif_usage_key") if c in frame.columns]

    classes: list[str] = []
    keys: list[str] = []
    for _, row in frame.iterrows():
        key = ""
        for column in key_columns:
            candidate = _text(row.get(column))
            if candidate:
                key = candidate
                break
        keys.append(key)

        if _text(row.get("reason")) in final_reasons:
            classes.append(BLOCKED_REASON)
        elif (
            _text(row.get("gbif_match_type")).upper() not in allowed_type
            or _text(row.get("gbif_status")).upper() not in allowed_status
        ):
            classes.append(BLOCKED_GBIF)
        elif not key:
            classes.append(BLOCKED_NO_KEY)
        else:
            classes.append(RECOVERABLE_KEY)

    frame["recovery_class"] = classes
    frame["gbif_key"] = keys
    frame["recovered_tier"] = rules.get("recovered_tier", "")
    return frame


def name_match_queue(classified: pd.DataFrame) -> pd.DataFrame:
    """The distinct GBIF keys that need one lookup each to finish recovery."""
    if classified.empty:
        return pd.DataFrame(columns=["gbif_key", "n_rows", "example_source_name"])
    recoverable = classified.loc[classified["recovery_class"].eq(RECOVERABLE_KEY)]
    if recoverable.empty:
        return pd.DataFrame(columns=["gbif_key", "n_rows", "example_source_name"])

    name_column = (
        "source_scientific_name" if "source_scientific_name" in recoverable.columns else None
    )
    rows: list[dict[str, Any]] = []
    for key, group in recoverable.groupby("gbif_key", sort=True):
        rows.append(
            {
                "gbif_key": key,
                "n_rows": int(len(group)),
                "example_source_name": (
                    _text(group.iloc[0][name_column]) if name_column else ""
                ),
            }
        )
    return pd.DataFrame(rows)


def resolve_name_match_queue(
    queue: pd.DataFrame,
    master_species: set[str],
    fetch: "Callable[[str], dict[str, Any] | None]",
) -> pd.DataFrame:
    """Turn resolved GBIF usage keys into master-species hits.

    ``fetch`` takes a usage key and returns the GBIF species record, or None
    when the key cannot be resolved. It is injected so the decision logic is
    testable without network; the CI workflow supplies the real client.

    A key is only accepted when GBIF reports an ACCEPTED species-rank record
    whose canonical name is present in the frozen island master. Anything else
    -- unresolved, synonym chains that end nowhere, non-species ranks, names
    outside the island universe -- is recorded with its reason and dropped.
    """
    columns = [
        "gbif_key",
        "n_rows",
        "example_source_name",
        "resolved_accepted_species",
        "resolved_rank",
        "resolved_status",
        "on_island_master",
        "outcome",
    ]
    if queue.empty:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, Any]] = []
    for _, item in queue.iterrows():
        key = _text(item["gbif_key"])
        record = fetch(key)
        resolved = _text((record or {}).get("species") or (record or {}).get("canonicalName"))
        rank = _text((record or {}).get("rank")).upper()
        status = _text((record or {}).get("taxonomicStatus")).upper()

        if record is None:
            outcome = "unresolved_key"
        elif rank != "SPECIES":
            outcome = "rejected_non_species_rank"
        elif status != "ACCEPTED":
            outcome = "rejected_not_accepted"
        elif not resolved:
            outcome = "rejected_no_canonical_name"
        elif resolved not in master_species:
            outcome = "rejected_off_island_master"
        else:
            outcome = "recovered"

        rows.append(
            {
                "gbif_key": key,
                "n_rows": int(item["n_rows"]),
                "example_source_name": _text(item.get("example_source_name")),
                "resolved_accepted_species": resolved,
                "resolved_rank": rank,
                "resolved_status": status,
                # True only for rows that cleared every gate. A rejected record
                # whose name happens to appear in the master must not read as a
                # master hit, or a downstream filter on this column would pull
                # rejected rows back in.
                "on_island_master": outcome == "recovered",
                "outcome": outcome,
            }
        )
    return pd.DataFrame(rows, columns=columns)


def classify_extraction(rejected: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Classify rejected candidates by what their retained quote actually says."""
    missing = EXTRACTION_REQUIRED - set(rejected.columns)
    if missing:
        raise typer.BadParameter(f"rejection ledger is missing columns: {sorted(missing)}")

    rules = config["extraction"]
    organs = [str(v).lower() for v in (rules.get("non_target_organ_terms") or ())]
    markers = [str(v) for v in (rules.get("enumeration_markers") or ())]
    min_hits = int(rules.get("min_enumeration_hits", 2))
    multilingual = rules.get("multilingual_terms") or {}
    aliases = config["language_aliases"]

    frame = rejected.copy()
    frame["quote"] = frame["exact_supporting_quote"].map(_text)
    frame["language_code"] = (
        frame["language"].map(lambda v: normalize_language(v, aliases))
        if "language" in frame.columns
        else "unknown"
    )

    classes: list[str] = []
    matched_terms: list[str] = []
    for _, row in frame.iterrows():
        quote = row["quote"].lower()
        if not quote:
            classes.append(NEEDS_REVIEW)
            matched_terms.append("")
            continue

        if any(term in quote for term in organs):
            classes.append(CORRECT_ORGAN)
            matched_terms.append("")
            continue

        if sum(quote.count(marker.lower()) for marker in markers) >= min_hits:
            classes.append(CORRECT_TEMPLATE)
            matched_terms.append("")
            continue

        # Language-conditioned recovery, but the term must actually be present.
        # A declared language alone never recovers a row.
        hit = ""
        for terms in multilingual.values():
            for term, canonical in terms.items():
                if term.lower() in quote:
                    hit = f"{term}->{canonical}"
                    break
            if hit:
                break

        if hit:
            classes.append(RECOVERABLE_LANG)
            matched_terms.append(hit)
        else:
            classes.append(NEEDS_REVIEW)
            matched_terms.append("")

    frame["recovery_class"] = classes
    frame["matched_multilingual_term"] = matched_terms
    return frame


def _summary(frame: pd.DataFrame, label: str) -> dict[str, Any]:
    counts = frame["recovery_class"].value_counts().to_dict() if len(frame) else {}
    total = int(len(frame))
    return {
        "lane": label,
        "n_rejected_rows": total,
        "classes": {str(k): int(v) for k, v in counts.items()},
        "class_share": (
            {str(k): round(int(v) / total, 4) for k, v in counts.items()} if total else {}
        ),
    }


@app.command("run")
def run(
    output_dir: Path = typer.Option(..., "--output-dir"),
    config_path: Path = typer.Option(Path("config/bulk_recovery.yml"), "--config", exists=True),
    name_match_audit: list[Path] = typer.Option([], "--name-match-audit"),
    rejected_candidates: list[Path] = typer.Option([], "--rejected-candidates"),
) -> None:
    """Classify both rejection lanes and write the recovery queues."""
    config = load_config(config_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {"version": "1.0", "lanes": []}

    if name_match_audit:
        frames = [pd.read_csv(path) for path in name_match_audit]
        classified = classify_name_match(pd.concat(frames, ignore_index=True), config)
        classified.to_csv(output_dir / "name_match_recovery_audit.csv.gz", index=False)
        queue = name_match_queue(classified)
        queue.to_csv(output_dir / "name_match_recovery_queue.csv", index=False)
        lane = _summary(classified, "name_match")
        lane["n_distinct_gbif_keys_to_resolve"] = int(len(queue))
        summary["lanes"].append(lane)

    if rejected_candidates:
        frames = [pd.read_csv(path) for path in rejected_candidates]
        classified = classify_extraction(pd.concat(frames, ignore_index=True), config)
        classified.to_csv(output_dir / "extraction_recovery_audit.csv.gz", index=False)
        recoverable = classified.loc[classified["recovery_class"].eq(RECOVERABLE_LANG)]
        recoverable.to_csv(output_dir / "extraction_recovery_queue.csv.gz", index=False)
        lane = _summary(classified, "extraction")
        lane["n_recoverable_rows"] = int(len(recoverable))
        lane["language_codes"] = {
            str(k): int(v) for k, v in classified["language_code"].value_counts().items()
        }
        summary["lanes"].append(lane)

    if not summary["lanes"]:
        raise typer.BadParameter(
            "give at least one --name-match-audit or --rejected-candidates input"
        )

    summary["interpretation"] = (
        "Counts are rejected rows re-examined against the pipeline's own ledgers. "
        "Recoverable rows are a queue for review, never accepted evidence. Rows "
        "classified as correctly rejected stay rejected."
    )
    (output_dir / "bulk_recovery_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
