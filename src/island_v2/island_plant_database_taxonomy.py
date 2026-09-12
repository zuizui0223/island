"""Versioned taxonomy-normalization layer for Island Plant Database 2.0.

This module never mutates the frozen alpha1 candidate-flora tables. It matches
provisional species-name identities against GBIF's current v2 matcher using the
Catalogue of Life Extended Release (COL XR) and emits an auditable crosswalk.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

COL_XR_CHECKLIST_KEY = "7ddf754f-d193-4cc9-b351-99906754a03b"
MATCH_URL = "https://api.gbif.org/v2/species/match"
METADATA_URL = "https://api.gbif.org/v2/species/match/metadata"
SOURCE_LICENSE = "CC-BY-4.0"
ALPHA1_TAXA_SUMMARY_SHA256 = "c0586264d9a877b88c26e866264f2c77b20423cd93380a8eebe269582d53352e"

OUTPUT_COLUMNS = [
    "submitted_name",
    "submitted_genus",
    "submitted_family",
    "matched_usage_key",
    "matched_usage_name",
    "matched_canonical_name",
    "matched_status",
    "matched_rank",
    "synonym",
    "candidate_accepted_key",
    "candidate_accepted_name",
    "candidate_accepted_authorship",
    "candidate_accepted_status",
    "candidate_accepted_rank",
    "accepted_target_basis",
    "candidate_genus",
    "candidate_family",
    "candidate_kingdom",
    "match_type",
    "confidence",
    "processing_flags",
    "issues",
    "automatic_resolution_candidate",
    "resolution_status",
    "review_status",
    "checklist_key",
    "source_license",
]


@app.callback()
def main() -> None:
    """Build non-destructive taxonomy crosswalks for Database 2.0."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _pipe(values: object) -> str:
    if not isinstance(values, list):
        return ""
    return "|".join(sorted({_text(value) for value in values if _text(value)}))


def load_source_taxa(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"accepted_species", "genus", "family"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"alpha1 taxonomy source missing columns: {missing}")
    frame = frame.copy()
    frame["accepted_species"] = frame["accepted_species"].map(_text)
    frame["genus"] = frame["genus"].map(_text)
    frame["family"] = frame["family"].map(_text)
    frame = frame.loc[frame["accepted_species"].ne("")]
    if frame["accepted_species"].duplicated().any():
        raise ValueError("alpha1 taxonomy source contains duplicate accepted_species")
    return frame.sort_values("accepted_species", kind="stable").reset_index(drop=True)


def deterministic_pilot(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    if n < 1:
        raise ValueError("pilot size must be positive")
    if n >= len(frame):
        return frame.copy()
    positions = [(i * len(frame)) // n for i in range(n)]
    return frame.iloc[positions].reset_index(drop=True)


def query_for_row(row: pd.Series) -> dict[str, Any]:
    query: dict[str, Any] = {
        "scientificName": _text(row["accepted_species"]),
        "taxonRank": "SPECIES",
        "kingdom": "Plantae",
        "strict": True,
        "verbose": False,
    }
    if _text(row.get("genus", "")):
        query["genus"] = _text(row["genus"])
    if _text(row.get("family", "")):
        query["family"] = _text(row["family"])
    return query


def _classification_value(result: dict[str, Any], rank: str) -> str:
    for item in result.get("classification") or []:
        if str(item.get("rank", "")).upper() == rank.upper():
            return _text(item.get("name"))
    return ""


def normalize_match(
    source_row: pd.Series,
    result: dict[str, Any],
    checklist_key: str = COL_XR_CHECKLIST_KEY,
) -> dict[str, Any]:
    usage = result.get("usage") or {}
    diagnostics = result.get("diagnostics") or {}
    synonym = bool(result.get("synonym", False))
    accepted_usage = result.get("acceptedUsage") or {}
    accepted = accepted_usage or (usage if not synonym else {})
    match_type = _text(diagnostics.get("matchType")).upper() or "NONE"
    confidence_raw = diagnostics.get("confidence")
    confidence = int(confidence_raw) if confidence_raw is not None else 0
    flags = diagnostics.get("processingFlags") or []
    issues = diagnostics.get("issues") or []

    target_rank = _text(accepted.get("rank")).upper()
    target_status = _text(accepted.get("status")).upper()
    matched_status = _text(usage.get("status")).upper()
    # In production v2 responses, acceptedUsage for an exact synonym can omit
    # its own status. The acceptedUsage field is itself the API's pointer to the
    # accepted concept, so its presence is sufficient only when the matched
    # usage is explicitly SYNONYM. We retain the raw blank target status.
    target_is_accepted = target_status == "ACCEPTED" or bool(
        synonym and accepted_usage and matched_status == "SYNONYM"
    )
    automatic = bool(
        usage
        and accepted
        and match_type == "EXACT"
        and confidence >= 95
        and not flags
        and not issues
        and target_rank == "SPECIES"
        and target_is_accepted
        and _classification_value(result, "KINGDOM").casefold() == "plantae"
    )

    if not usage or match_type == "NONE":
        resolution_status = "unmatched"
        review_status = "needs_taxonomy_review"
    elif automatic and synonym:
        resolution_status = "exact_synonym_to_accepted_candidate"
        review_status = "automatic_candidate_not_promoted"
    elif automatic:
        resolution_status = "exact_accepted_candidate"
        review_status = "automatic_candidate_not_promoted"
    else:
        resolution_status = "review_required"
        review_status = "needs_taxonomy_review"

    return {
        "submitted_name": _text(source_row["accepted_species"]),
        "submitted_genus": _text(source_row.get("genus", "")),
        "submitted_family": _text(source_row.get("family", "")),
        "matched_usage_key": _text(usage.get("key")),
        "matched_usage_name": _text(usage.get("name")),
        "matched_canonical_name": _text(usage.get("canonicalName")),
        "matched_status": _text(usage.get("status")),
        "matched_rank": _text(usage.get("rank")),
        "synonym": str(synonym).lower(),
        "candidate_accepted_key": _text(accepted.get("key")),
        "candidate_accepted_name": _text(accepted.get("canonicalName") or accepted.get("name")),
        "candidate_accepted_authorship": _text(accepted.get("authorship")),
        "candidate_accepted_status": _text(accepted.get("status")),
        "candidate_accepted_rank": _text(accepted.get("rank")),
        "accepted_target_basis": "acceptedUsage" if accepted_usage else ("usage" if accepted else ""),
        "candidate_genus": _classification_value(result, "GENUS"),
        "candidate_family": _classification_value(result, "FAMILY"),
        "candidate_kingdom": _classification_value(result, "KINGDOM"),
        "match_type": match_type,
        "confidence": confidence,
        "processing_flags": _pipe(flags),
        "issues": _pipe(issues),
        "automatic_resolution_candidate": str(automatic).lower(),
        "resolution_status": resolution_status,
        "review_status": review_status,
        "checklist_key": checklist_key,
        "source_license": SOURCE_LICENSE,
    }


def _request_with_retry(
    client: httpx.Client,
    method: str,
    url: str,
    *,
    params: dict[str, str],
    json_body: object | None = None,
    attempts: int = 4,
) -> Any:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = client.request(method, url, params=params, json=json_body)
            response.raise_for_status()
            return response.json()
        except (httpx.HTTPError, ValueError) as exc:
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"GBIF taxonomy request failed after {attempts} attempts: {last_error}")


def fetch_metadata(client: httpx.Client, checklist_key: str) -> dict[str, Any]:
    payload = _request_with_retry(
        client,
        "GET",
        METADATA_URL,
        params={"checklistKey": checklist_key},
    )
    if not isinstance(payload, dict):
        raise RuntimeError("GBIF taxonomy metadata response is not an object")
    return payload


def match_rows(
    source: pd.DataFrame,
    *,
    checklist_key: str = COL_XR_CHECKLIST_KEY,
    batch_size: int = 1000,
    timeout_seconds: float = 90.0,
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    if not 1 <= batch_size <= 1000:
        raise ValueError("GBIF v2 batch_size must be between 1 and 1000")
    normalized: list[dict[str, Any]] = []
    raw_batches: list[dict[str, Any]] = []
    headers = {"User-Agent": "island-plant-database/2.0 taxonomy-normalization"}
    with httpx.Client(timeout=timeout_seconds, headers=headers, follow_redirects=True) as client:
        metadata = fetch_metadata(client, checklist_key)
        for start in range(0, len(source), batch_size):
            batch = source.iloc[start : start + batch_size]
            queries = [query_for_row(row) for _, row in batch.iterrows()]
            results = _request_with_retry(
                client,
                "POST",
                MATCH_URL,
                params={"checklistKey": checklist_key},
                json_body=queries,
            )
            if not isinstance(results, list) or len(results) != len(queries):
                raise RuntimeError(
                    f"GBIF batch response cardinality mismatch: {len(queries)} queries, "
                    f"{len(results) if isinstance(results, list) else 'non-list'} results"
                )
            raw_batches.append({"start": start, "queries": queries, "results": results})
            normalized.extend(
                normalize_match(row, result, checklist_key)
                for (_, row), result in zip(batch.iterrows(), results, strict=True)
            )
    return pd.DataFrame(normalized, columns=OUTPUT_COLUMNS), raw_batches, metadata


def summarize(frame: pd.DataFrame) -> dict[str, Any]:
    return {
        "n_input": int(len(frame)),
        "n_automatic_resolution_candidates": int(
            frame["automatic_resolution_candidate"].astype(str).eq("true").sum()
        ),
        "n_review_required": int(frame["resolution_status"].eq("review_required").sum()),
        "n_unmatched": int(frame["resolution_status"].eq("unmatched").sum()),
        "match_type_counts": frame["match_type"].value_counts(dropna=False).sort_index().to_dict(),
        "resolution_status_counts": (
            frame["resolution_status"].value_counts(dropna=False).sort_index().to_dict()
        ),
    }


@app.command("pilot")
def pilot(
    source_taxa: Path = typer.Option(..., "--source-taxa", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
    pilot_size: int = typer.Option(512, "--pilot-size", min=1, max=5000),
    batch_size: int = typer.Option(512, "--batch-size", min=1, max=1000),
    expected_source_sha256: str = typer.Option(ALPHA1_TAXA_SUMMARY_SHA256, "--expected-source-sha256"),
) -> None:
    actual_sha = sha256_file(source_taxa)
    if actual_sha != expected_source_sha256:
        raise typer.BadParameter(
            f"taxonomy source SHA mismatch: expected {expected_source_sha256}, got {actual_sha}"
        )
    source = load_source_taxa(source_taxa)
    sample = deterministic_pilot(source, pilot_size)
    crosswalk, raw_batches, metadata = match_rows(sample, batch_size=batch_size)

    output_dir.mkdir(parents=True, exist_ok=True)
    crosswalk.to_csv(output_dir / "taxonomy_crosswalk_pilot.csv", index=False)
    (output_dir / "gbif_colxr_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    raw_dir = output_dir / "raw_batches"
    raw_dir.mkdir(exist_ok=True)
    for index, payload in enumerate(raw_batches):
        (raw_dir / f"batch_{index:04d}.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )

    summary = summarize(crosswalk)
    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "target_version": "2.0.0-alpha2-taxonomy",
        "stage": "colxr_taxonomy_pilot",
        "source_database_version": "2.0.0-alpha1",
        "source_taxa_summary_sha256": actual_sha,
        "source_taxa_total": int(len(source)),
        "pilot_sampling": "alphabetically_stable_evenly_spaced",
        "pilot_size": int(len(sample)),
        "batch_size": int(batch_size),
        "checklist_key": COL_XR_CHECKLIST_KEY,
        "taxonomy": "Catalogue of Life Extended Release",
        "source_license": SOURCE_LICENSE,
        "promotion_policy": (
            "No alpha1 taxon is mutated. Only EXACT, confidence>=95, flag-free, "
            "species-rank Plantae matches to accepted concepts are automatic-resolution candidates; "
            "for a matched SYNONYM, the v2 acceptedUsage field is the accepted target even when "
            "its optional status field is absent. All other matches remain review-required or unmatched."
        ),
        "summary": summary,
    }
    (output_dir / "TAXONOMY_PILOT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    app()
