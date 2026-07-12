"""High-relevance search-result evidence for plant traits.

This experiment prioritizes exact scientific-name relevance over raw search-result count.
It queries DuckDuckGo HTML, retains only results whose title/snippet explicitly contains
the searched binomial, and freezes the search snippet as an auditable evidence lead.
Fetched page text is preferred when available; the snippet remains usable when anti-bot,
JavaScript or TLS prevents page retrieval.
"""

from __future__ import annotations

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote_plus

import httpx
import pandas as pd
import typer

from island_v2.trait_dynamic_discovery import (
    USER_AGENT,
    _clean,
    _parse_duck_html,
    fetch_source,
    gbif_name_bundle,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

QUERY_TEMPLATES = (
    '"{name}" flowers',
    '"{name}" flower colour color corolla',
    '"{name}" floral morphology',
    '"{name}" pollination breeding system',
    '"{name}" "self-compatible"',
    '"{name}" "self-incompatible"',
    '"{name}" "self-fertile"',
    '"{name}" "self-sterile"',
    '"{name}" autogamy selfing',
    '"{name}" outcrossing',
    '"{name}" "breeding system"',
    '"{name}" "mating system"',
)

OUTPUT_COLUMNS = [
    "accepted_species",
    "searched_name",
    "query_family",
    "search_query",
    "search_engine",
    "search_rank",
    "search_result_title",
    "search_result_snippet",
    "source_text",
    "source_url",
    "source_citation",
    "source_type",
    "evidence_scope",
    "content_sha256",
    "evidence_carrier",
]


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _binomial(name: str) -> str:
    parts = _clean(name).lower().split()
    return " ".join(parts[:2]) if len(parts) >= 2 else _clean(name).lower()


def _relevant(text: str, name: str) -> bool:
    return _binomial(name) in _clean(text).lower()


def _query_family(template: str) -> str:
    low = template.lower()
    if any(
        token in low
        for token in (
            "pollination",
            "self-compatible",
            "self-incompatible",
            "self-fertile",
            "self-sterile",
            "autogamy",
            "selfing",
            "outcrossing",
            "breeding system",
            "mating system",
        )
    ):
        return "reproductive_ecology"
    return "floral_morphology"


def _source_type(url: str, carrier: str) -> str:
    if carrier == "search_result_snippet":
        return "search_result_snippet"
    low = url.lower()
    if low.endswith(".pdf"):
        return "paper_or_flora_pdf"
    if any(token in low for token in ("flora", "efloras", "kew.org", "powo", "nzpcn")):
        return "flora_or_monograph"
    return "public_web_page"


def _search_duck(client: httpx.Client, query: str, limit: int) -> list[dict[str, str]]:
    response = client.get(
        f"https://html.duckduckgo.com/html/?q={quote_plus(query)}",
        timeout=25,
        follow_redirects=True,
    )
    response.raise_for_status()
    return _parse_duck_html(response.text)[:limit]


def _one_species(
    species: str,
    max_synonyms: int,
    max_results: int,
    max_sources: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    search_rows: list[dict[str, str]] = []
    evidence_rows: list[dict[str, str]] = []
    seen_urls: set[str] = set()
    with httpx.Client(headers=headers) as client:
        names = gbif_name_bundle(client, species, max_synonyms=max_synonyms)
        for searched_name in names:
            for template in QUERY_TEMPLATES:
                query = template.format(name=searched_name)
                try:
                    results = _search_duck(client, query, max_results)
                except Exception:
                    continue
                for rank, result in enumerate(results, start=1):
                    title = _clean(result.get("title"))
                    snippet = _clean(result.get("snippet"))
                    url = _clean(result.get("url"))
                    combined = f"{title} {snippet}"
                    if not url or not _relevant(combined, searched_name):
                        continue
                    base = {
                        "accepted_species": species,
                        "searched_name": searched_name,
                        "query_family": _query_family(template),
                        "search_query": query,
                        "search_engine": "duckduckgo_html",
                        "search_rank": str(rank),
                        "search_result_title": title,
                        "search_result_snippet": snippet,
                        "source_url": url,
                        "source_citation": title or url,
                        "evidence_scope": "species_direct_search_result",
                    }
                    search_rows.append(base.copy())
                    if snippet:
                        evidence_rows.append(
                            {
                                **base,
                                "source_text": snippet,
                                "source_type": _source_type(url, "search_result_snippet"),
                                "content_sha256": _sha(snippet),
                                "evidence_carrier": "search_result_snippet",
                            }
                        )
                    if url in seen_urls or len(seen_urls) >= max_sources:
                        continue
                    seen_urls.add(url)
                    page_text = fetch_source(client, url)
                    if page_text and _relevant(page_text, searched_name):
                        evidence_rows.append(
                            {
                                **base,
                                "source_text": page_text,
                                "source_type": _source_type(url, "fetched_page"),
                                "content_sha256": _sha(page_text),
                                "evidence_carrier": "fetched_page",
                                "evidence_scope": "species_direct",
                            }
                        )
    return search_rows, evidence_rows


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    species_column: str = typer.Option("accepted_species"),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_synonyms: int = typer.Option(1, min=0, max=10),
    max_results_per_query: int = typer.Option(8, min=1, max=20),
    max_sources_per_species: int = typer.Option(8, min=1, max=30),
    workers: int = typer.Option(8, min=1, max=32),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    if species_column not in table.columns:
        raise typer.BadParameter(f"missing species column: {species_column}")
    species: list[str] = []
    seen: set[str] = set()
    for value in table[species_column]:
        name = _clean(value)
        if name and name not in seen:
            seen.add(name)
            species.append(name)
        if len(species) >= max_taxa:
            break

    all_search: list[dict[str, str]] = []
    all_evidence: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _one_species,
                name,
                max_synonyms,
                max_results_per_query,
                max_sources_per_species,
            ): name
            for name in species
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                search_rows, evidence_rows = future.result()
                all_search.extend(search_rows)
                all_evidence.extend(evidence_rows)
            except Exception as exc:
                errors.append({"accepted_species": name, "error": str(exc)})

    search = pd.DataFrame(all_search)
    evidence = pd.DataFrame(all_evidence, columns=OUTPUT_COLUMNS)
    if len(evidence):
        evidence = evidence.drop_duplicates(
            ["accepted_species", "source_url", "content_sha256", "evidence_carrier"]
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    search.to_csv(output_dir / "relevant_search_results.csv", index=False)
    evidence.to_csv(output_dir / "source_texts.csv", index=False)
    pd.DataFrame(errors, columns=["accepted_species", "error"]).to_csv(
        output_dir / "search_errors.csv", index=False
    )
    summary = {
        "n_species_requested": len(species),
        "n_species_with_relevant_search_results": int(search["accepted_species"].nunique()) if len(search) else 0,
        "n_species_with_evidence_text": int(evidence["accepted_species"].nunique()) if len(evidence) else 0,
        "n_relevant_search_results": len(search),
        "n_evidence_rows": len(evidence),
        "carrier_counts": evidence["evidence_carrier"].value_counts().to_dict() if len(evidence) else {},
        "n_errors": len(errors),
    }
    (output_dir / "relevant_search_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
