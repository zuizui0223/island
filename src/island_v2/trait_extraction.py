"""Run broad web-search trait extraction for a manually selected taxon batch.

This module writes only reviewable staging artifacts. It never edits the curated
trait table and it never commits data back to GitHub. The corresponding Actions
workflow uploads the output as an artifact for review.

A large global batch must be resilient to one failed API request: raw responses,
completed candidate rows, and a machine-readable failure record are written after
every completed sub-batch and again before a failing run exits.
"""

from __future__ import annotations

import hashlib
import json
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml
from openai import OpenAI

from island_v2.schemas import TraitEvidenceCandidate, TraitExtractionResponse

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_TAXON_COLUMNS = {"accepted_species", "genus", "family"}
RETRIABLE_STATUS_CODES = {408, 409, 429, 500, 502, 503, 504}
SUPPORTED_FORMATS = {
    "date-time",
    "time",
    "date",
    "duration",
    "email",
    "hostname",
    "ipv4",
    "ipv6",
    "uuid",
}
UNSUPPORTED_SCHEMA_KEYS = {"default", "examples", "title", "$schema"}


@app.callback()
def main() -> None:
    """Search the web and write reviewable trait-evidence candidates."""


def read_text(path: Path) -> str:
    if not path.exists():
        raise typer.BadParameter(f"File does not exist: {path}")
    return path.read_text(encoding="utf-8")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def strict_json_schema(model: type[TraitExtractionResponse]) -> dict[str, Any]:
    """Normalize a Pydantic schema to the supported strict-outputs subset.

    The Responses API requires every object property to be required and every
    object to set ``additionalProperties: false``. Pydantic also emits metadata
    such as ``default`` and ``title``; those are not biological content and are
    removed so API compatibility does not depend on SDK-specific schema details.
    """
    schema = model.model_json_schema()

    def visit(node: Any) -> None:
        if isinstance(node, dict):
            for key in UNSUPPORTED_SCHEMA_KEYS:
                node.pop(key, None)
            if node.get("format") not in SUPPORTED_FORMATS:
                node.pop("format", None)
            properties = node.get("properties")
            if isinstance(properties, dict):
                node["additionalProperties"] = False
                node["required"] = list(properties)
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)

    visit(schema)
    return schema


def candidate_columns() -> list[str]:
    """Return stable CSV columns even when a request fails before any row exists."""
    return [
        "run_id",
        "accepted_species",
        "taxon_id",
        "batch_notes",
        "no_evidence_found",
        *TraitEvidenceCandidate.model_fields.keys(),
        "retrieved_source_urls_json",
        "input_taxon_context_json",
        "input_trait_candidate_status",
        "input_release_gate",
        "review_status",
    ]


def load_taxa(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    """Load selected taxa while preserving any review-gate metadata columns."""
    taxa = pd.read_csv(path, dtype=str).fillna("")
    missing = REQUIRED_TAXON_COLUMNS.difference(taxa.columns)
    if missing:
        raise typer.BadParameter(
            f"Taxon input is missing required columns: {', '.join(sorted(missing))}"
        )
    taxa["accepted_species"] = taxa["accepted_species"].astype(str).str.strip()
    taxa = taxa.loc[taxa["accepted_species"].ne("")].drop_duplicates(subset=["accepted_species"])
    if limit is not None:
        taxa = taxa.head(limit)
    if taxa.empty:
        raise typer.BadParameter("Taxon input has no species to process.")
    return taxa.to_dict(orient="records")


def taxon_context_by_species(taxa: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Retain input provenance beside each extracted candidate, without changing traits.

    The trait extractor can receive a general taxon list or a gated pilot list.
    Additional input columns such as ``trait_candidate_status`` and ``release_gate``
    are copied into candidate rows as provenance. They are never used in the web
    search prompt and never influence extracted trait values.
    """
    contexts: dict[str, dict[str, str]] = {}
    for taxon in taxa:
        species = str(taxon.get("accepted_species", "")).strip()
        if not species:
            continue
        contexts[species] = {
            str(key): str(value)
            for key, value in taxon.items()
            if key not in {"accepted_species", "genus", "family"} and str(value).strip()
        }
    return contexts


def build_user_request(taxa: list[dict[str, str]], ontology: dict[str, Any]) -> str:
    """Build one transparent task request for a manually selected taxon batch."""
    trait_names = list(ontology["traits"])
    taxon_lines = "\n".join(
        f"- accepted_species: {row['accepted_species']}; genus: {row['genus']}; family: {row['family']}"
        for row in taxa
    )
    return f"""Search the live web for the following accepted plant taxa and return a v2.1 structured trait-candidate response.

Taxa:\n{taxon_lines}

Allowed canonical trait names:\n{', '.join(trait_names)}

For every trait candidate, use a web-retrieved source URL when possible. Broad retrieval is expected: scholarly sources, floras, herbarium or botanic-garden pages, curated databases, and credible specialist web pages are all eligible. When direct species information is sparse, genus/family inference is allowed only as a declared hierarchical candidate with supporting taxa and a stated transfer rule.

Do not return a complete matrix by force. Return unresolved candidates only after searching broadly. Every candidate must require human review."""


def response_to_jsonable(response: Any) -> dict[str, Any]:
    """Serialize SDK response objects without depending on one SDK object layout."""
    if hasattr(response, "model_dump"):
        return response.model_dump(mode="json")
    if hasattr(response, "to_dict"):
        return response.to_dict()
    return json.loads(json.dumps(response, default=str))


def extract_url_citations(node: Any) -> list[dict[str, str]]:
    """Extract URL citations and search-source links from an API response payload."""
    found: dict[str, dict[str, str]] = {}

    def add(url: Any, title: Any = "") -> None:
        if isinstance(url, str) and url.startswith(("https://", "http://")):
            found.setdefault(url, {"url": url, "title": str(title or "")})

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            if value.get("type") == "url_citation":
                add(value.get("url"), value.get("title"))
            if "url" in value and any(key in value for key in ("title", "name", "source")):
                add(value.get("url"), value.get("title") or value.get("name"))
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    walk(node)
    return list(found.values())


def flatten_candidates(
    parsed: TraitExtractionResponse,
    run_id: str,
    response_sources: list[dict[str, str]],
    taxon_contexts: dict[str, dict[str, str]] | None = None,
) -> list[dict[str, Any]]:
    """Flatten trait candidates while retaining taxon-selection provenance."""
    contexts = taxon_contexts or {}
    rows: list[dict[str, Any]] = []
    for species in parsed.species:
        context = contexts.get(species.accepted_species, {})
        for candidate in species.trait_candidates:
            row = {
                "run_id": run_id,
                "accepted_species": species.accepted_species,
                "taxon_id": species.taxon_id,
                "batch_notes": species.batch_notes,
                "no_evidence_found": species.no_evidence_found,
                **candidate.model_dump(mode="json"),
                "retrieved_source_urls_json": json.dumps(response_sources, ensure_ascii=False),
                "input_taxon_context_json": json.dumps(context, ensure_ascii=False, sort_keys=True),
                "input_trait_candidate_status": context.get("trait_candidate_status", ""),
                "input_release_gate": context.get("release_gate", ""),
                "review_status": "pending",
            }
            row["supporting_taxa"] = json.dumps(row["supporting_taxa"], ensure_ascii=False)
            rows.append(row)
    return rows


def safe_error_payload(error: Exception) -> dict[str, Any]:
    """Keep useful API diagnostics without recording environment secrets."""
    message = str(error).replace("\n", " ").strip()
    if len(message) > 1200:
        message = f"{message[:1200]}…"
    return {
        "exception_type": type(error).__name__,
        "message": message,
        "status_code": getattr(error, "status_code", None),
        "request_id": getattr(error, "request_id", None),
    }


def is_retriable_error(error: Exception) -> bool:
    """Return whether an API/network error can reasonably succeed on retry."""
    status_code = getattr(error, "status_code", None)
    if status_code in RETRIABLE_STATUS_CODES:
        return True
    return type(error).__name__ in {"APIConnectionError", "APITimeoutError", "RateLimitError"}


def create_response_with_retry(
    client: OpenAI,
    request: dict[str, Any],
    max_retries: int,
) -> tuple[Any, list[dict[str, Any]]]:
    """Call Responses API with bounded retries and an auditable attempt history."""
    attempt_history: list[dict[str, Any]] = []
    for attempt in range(1, max_retries + 2):
        try:
            response = client.responses.create(**request)
            return response, attempt_history
        except Exception as error:  # API and transport errors differ across SDK releases.
            detail = safe_error_payload(error)
            detail["attempt"] = attempt
            detail["retriable"] = is_retriable_error(error)
            attempt_history.append(detail)
            if not detail["retriable"] or attempt > max_retries:
                raise RuntimeError(json.dumps(detail, ensure_ascii=False)) from error
            time.sleep(min(2 ** (attempt - 1), 8))
    raise RuntimeError("Unreachable retry state")


def response_request(
    *,
    model: str,
    prompt: str,
    user_request: str,
    schema: dict[str, Any],
) -> dict[str, Any]:
    """Build one supported, bounded Responses API request.

    The default web-search returned-token budget is intentionally used. Unlimited
    tool output is reserved for deep-research tasks and can inflate latency and
    failure exposure in high-throughput trait batches.
    """
    return {
        "model": model,
        "reasoning": {"effort": "medium"},
        "tools": [{"type": "web_search", "external_web_access": True}],
        "tool_choice": "required",
        "include": ["web_search_call.action.sources"],
        "text": {
            "format": {
                "type": "json_schema",
                "name": "island_v2_trait_candidates",
                "strict": True,
                "schema": schema,
            }
        },
        "input": [
            {"role": "system", "content": prompt},
            {"role": "user", "content": user_request},
        ],
    }


def write_run_artifacts(
    *,
    output_dir: Path,
    run_id: str,
    rows: list[dict[str, Any]],
    manifest: dict[str, Any],
) -> tuple[Path, Path]:
    """Write repeatable snapshots so completed batches survive a later failure."""
    candidate_path = output_dir / f"trait_candidates_{run_id}.csv"
    pd.DataFrame(rows, columns=candidate_columns()).to_csv(candidate_path, index=False)
    manifest_path = output_dir / f"manifest_{run_id}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return candidate_path, manifest_path


@app.command()
def run(
    taxa_csv: Path = typer.Option(..., exists=True, help="CSV with accepted_species, genus, family."),
    output_dir: Path = typer.Option(..., help="Directory for reviewable run artifacts."),
    prompt_path: Path = typer.Option(
        Path("prompts/trait_evidence_extraction_v2.md"), help="Versioned LLM instruction file."
    ),
    ontology_path: Path = typer.Option(
        Path("config/trait_ontology.yml"), help="Versioned trait ontology."
    ),
    model: str = typer.Option("gpt-5.5", help="Responses API model."),
    batch_size: int = typer.Option(10, min=1, max=25, help="Species per web-search request."),
    limit: int | None = typer.Option(None, min=1, help="Optional cap for a pilot run."),
    preflight_taxa: int = typer.Option(
        1,
        min=0,
        max=1,
        help="Run the first actual taxon as a one-taxon compatibility preflight before full chunks.",
    ),
    max_retries: int = typer.Option(2, min=0, max=5, help="Retries for transient API failures."),
) -> None:
    """Search the web and write LLM trait candidates as reviewable artifacts."""
    prompt = read_text(prompt_path)
    ontology_text = read_text(ontology_path)
    ontology = yaml.safe_load(ontology_text)
    taxa = load_taxa(taxa_csv, limit=limit)
    contexts = taxon_context_by_species(taxa)
    schema = strict_json_schema(TraitExtractionResponse)

    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "raw_responses"
    raw_dir.mkdir(exist_ok=True)
    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    client = OpenAI()
    all_rows: list[dict[str, Any]] = []
    manifests: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    batches: list[tuple[str, list[dict[str, str]]]] = []
    if preflight_taxa and taxa:
        batches.append(("preflight", taxa[:preflight_taxa]))
        remaining_taxa = taxa[preflight_taxa:]
    else:
        remaining_taxa = taxa
    for start in range(0, len(remaining_taxa), batch_size):
        batches.append(("extraction", remaining_taxa[start : start + batch_size]))

    base_manifest: dict[str, Any] = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "model": model,
        "prompt_path": str(prompt_path),
        "prompt_sha256": sha256_text(prompt),
        "ontology_path": str(ontology_path),
        "ontology_sha256": sha256_text(ontology_text),
        "structured_output_schema_sha256": sha256_text(json.dumps(schema, sort_keys=True)),
        "taxa_csv": str(taxa_csv),
        "taxa_count": len(taxa),
        "batch_size": batch_size,
        "preflight_taxa": preflight_taxa,
        "max_retries": max_retries,
        "input_context_columns": sorted({key for context in contexts.values() for key in context}),
        "curation_rule": "All rows are candidates; review_status=pending is required before curation.",
    }

    for batch_index, (batch_kind, batch) in enumerate(batches, start=1):
        request = response_request(
            model=model,
            prompt=prompt,
            user_request=build_user_request(batch, ontology),
            schema=schema,
        )
        entry: dict[str, Any] = {
            "batch_index": batch_index,
            "batch_kind": batch_kind,
            "species": [item["accepted_species"] for item in batch],
        }
        try:
            response, retry_history = create_response_with_retry(client, request, max_retries)
            raw_payload = response_to_jsonable(response)
            raw_path = raw_dir / f"{run_id}_batch_{batch_index:03d}.json"
            raw_path.write_text(json.dumps(raw_payload, indent=2, ensure_ascii=False), encoding="utf-8")
            parsed = TraitExtractionResponse.model_validate_json(response.output_text)
            sources = extract_url_citations(raw_payload)
            batch_rows = flatten_candidates(parsed, run_id, sources, contexts)
            all_rows.extend(batch_rows)
            entry.update(
                {
                    "status": "completed",
                    "raw_response": str(raw_path),
                    "retrieved_source_count": len(sources),
                    "candidate_rows": len(batch_rows),
                    "retry_history": retry_history,
                }
            )
            manifests.append(entry)
        except Exception as error:
            failure = safe_error_payload(error)
            failure.update(
                {
                    "batch_index": batch_index,
                    "batch_kind": batch_kind,
                    "species": entry["species"],
                }
            )
            failures.append(failure)
            entry.update({"status": "failed", "failure": failure})
            manifests.append(entry)

        status = "complete" if not failures else "partial_failure"
        snapshot = {
            **base_manifest,
            "status": status,
            "n_completed_batches": sum(item["status"] == "completed" for item in manifests),
            "n_failed_batches": len(failures),
            "candidate_row_count": len(all_rows),
            "batches": manifests,
            "failures": failures,
        }
        candidate_path, manifest_path = write_run_artifacts(
            output_dir=output_dir,
            run_id=run_id,
            rows=all_rows,
            manifest=snapshot,
        )
        if failures:
            typer.echo(
                f"Partial trait extraction failure after batch {batch_index}; "
                f"wrote {len(all_rows)} candidates and diagnostics to {manifest_path}",
                err=True,
            )
            raise typer.Exit(code=1)

    typer.echo(f"Wrote {len(all_rows)} candidates to {candidate_path}")
    typer.echo(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    app()
