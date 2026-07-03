"""Run broad web-search trait extraction for a manually selected taxon batch.

This module writes only reviewable staging artifacts. It never edits the curated
trait table and it never commits data back to GitHub. The corresponding Actions
workflow uploads the output as an artifact for review.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml
from openai import OpenAI

from island_v2.schemas import TraitExtractionResponse

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_TAXON_COLUMNS = {"accepted_species", "genus", "family"}


def read_text(path: Path) -> str:
    if not path.exists():
        raise typer.BadParameter(f"File does not exist: {path}")
    return path.read_text(encoding="utf-8")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def strict_json_schema(model: type[TraitExtractionResponse]) -> dict[str, Any]:
    """Convert a Pydantic schema to the strict object form required by the API."""

    schema = model.model_json_schema()

    def visit(node: Any) -> None:
        if isinstance(node, dict):
            properties = node.get("properties")
            if isinstance(properties, dict):
                node["additionalProperties"] = False
                node["required"] = list(properties)
                for child in properties.values():
                    visit(child)
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)

    visit(schema)
    return schema


def load_taxa(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    taxa = pd.read_csv(path, dtype=str).fillna("")
    missing = REQUIRED_TAXON_COLUMNS.difference(taxa.columns)
    if missing:
        raise typer.BadParameter(
            f"Taxon input is missing required columns: {', '.join(sorted(missing))}"
        )
    taxa = taxa.drop_duplicates(subset=["accepted_species"])
    if limit is not None:
        taxa = taxa.head(limit)
    if taxa.empty:
        raise typer.BadParameter("Taxon input has no species to process.")
    return taxa.to_dict(orient="records")


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
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for species in parsed.species:
        for candidate in species.trait_candidates:
            row = {
                "run_id": run_id,
                "accepted_species": species.accepted_species,
                "taxon_id": species.taxon_id,
                "batch_notes": species.batch_notes,
                "no_evidence_found": species.no_evidence_found,
                **candidate.model_dump(mode="json"),
                "retrieved_source_urls_json": json.dumps(response_sources, ensure_ascii=False),
                "review_status": "pending",
            }
            row["supporting_taxa"] = json.dumps(row["supporting_taxa"], ensure_ascii=False)
            rows.append(row)
    return rows


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
) -> None:
    """Search the web and write LLM trait candidates as reviewable artifacts."""

    prompt = read_text(prompt_path)
    ontology_text = read_text(ontology_path)
    ontology = yaml.safe_load(ontology_text)
    taxa = load_taxa(taxa_csv, limit=limit)

    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = output_dir / "raw_responses"
    raw_dir.mkdir(exist_ok=True)
    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    client = OpenAI()
    all_rows: list[dict[str, Any]] = []
    manifests: list[dict[str, Any]] = []

    for batch_index, start in enumerate(range(0, len(taxa), batch_size), start=1):
        batch = taxa[start : start + batch_size]
        response = client.responses.create(
            model=model,
            reasoning={"effort": "medium"},
            tools=[
                {
                    "type": "web_search",
                    "external_web_access": True,
                    "return_token_budget": "unlimited",
                }
            ],
            tool_choice="required",
            include=["web_search_call.action.sources"],
            text={
                "format": {
                    "type": "json_schema",
                    "name": "island_v2_trait_candidates",
                    "strict": True,
                    "schema": strict_json_schema(TraitExtractionResponse),
                }
            },
            input=[
                {"role": "system", "content": prompt},
                {"role": "user", "content": build_user_request(batch, ontology)},
            ],
        )
        raw_payload = response_to_jsonable(response)
        raw_path = raw_dir / f"{run_id}_batch_{batch_index:03d}.json"
        raw_path.write_text(json.dumps(raw_payload, indent=2, ensure_ascii=False), encoding="utf-8")

        parsed = TraitExtractionResponse.model_validate_json(response.output_text)
        sources = extract_url_citations(raw_payload)
        all_rows.extend(flatten_candidates(parsed, run_id, sources))
        manifests.append(
            {
                "batch_index": batch_index,
                "species": [item["accepted_species"] for item in batch],
                "raw_response": str(raw_path),
                "retrieved_source_count": len(sources),
            }
        )

    candidate_path = output_dir / f"trait_candidates_{run_id}.csv"
    pd.DataFrame(all_rows).to_csv(candidate_path, index=False)
    manifest = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "model": model,
        "prompt_path": str(prompt_path),
        "prompt_sha256": sha256_text(prompt),
        "ontology_path": str(ontology_path),
        "ontology_sha256": sha256_text(ontology_text),
        "taxa_csv": str(taxa_csv),
        "taxa_count": len(taxa),
        "batch_size": batch_size,
        "candidate_csv": str(candidate_path),
        "batches": manifests,
        "curation_rule": "All rows are candidates; review_status=pending is required before curation.",
    }
    manifest_path = output_dir / f"manifest_{run_id}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    typer.echo(f"Wrote {len(all_rows)} candidates to {candidate_path}")
    typer.echo(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    app()
