"""Recover reproductive traits from full pages discovered by exact-phrase web search.

Search-result snippets often identify the right source without exposing the decisive sentence.
This pilot follows the highest-ranked source-agnostic results, freezes extractable HTML/PDF text,
and searches the full text for species/genus-scoped reproductive evidence.
"""

from __future__ import annotations

import io
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote_plus, urlparse

import httpx
import pandas as pd
import typer
from pypdf import PdfReader

from island_v2.reproductive_web_recovery_pilot import (
    MATING,
    SI,
    USER_AGENT,
    clean,
    genus,
    parse_bing,
    parse_duck,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

QUERY_GROUPS = {
    "self_incompatibility": (
        '"{name}" "self-compatible"',
        '"{name}" "self-incompatible"',
        '"{name}" "self compatibility"',
        '"{name}" "self incompatibility"',
        '"{name}" "self-fertile"',
        '"{name}" "self-sterile"',
    ),
    "mating_system": (
        '"{name}" autogamy autogamous selfing',
        '"{name}" outcrossing outcrossed',
        '"{name}" "mixed mating"',
        '"{name}" "breeding system"',
        '"{name}" "mating system"',
        '"{name}" cleistogamy cleistogamous',
    ),
}

OUTPUT_COLUMNS = [
    "accepted_species",
    "trait_group",
    "value",
    "inference_status",
    "confidence",
    "query_taxon",
    "search_scope",
    "search_engine",
    "search_query",
    "search_rank",
    "source_url",
    "source_title",
    "evidence_quote",
    "source_type",
]


class TextParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []
        self.hidden = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() in {"script", "style", "noscript", "svg"}:
            self.hidden += 1

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() in {"script", "style", "noscript", "svg"} and self.hidden:
            self.hidden -= 1

    def handle_data(self, data: str) -> None:
        if not self.hidden:
            self.parts.append(data)

    def text(self) -> str:
        return clean(" ".join(self.parts))


def fetch_text(client: httpx.Client, url: str, max_pdf_pages: int = 60) -> tuple[str, str]:
    try:
        response = client.get(url, timeout=30, follow_redirects=True)
        response.raise_for_status()
    except Exception:
        return "", ""
    content_type = response.headers.get("content-type", "").lower()
    if "pdf" in content_type or urlparse(str(response.url)).path.lower().endswith(".pdf"):
        try:
            reader = PdfReader(io.BytesIO(response.content))
            text = clean(" ".join((page.extract_text() or "") for page in reader.pages[:max_pdf_pages]))
            return text, "pdf"
        except Exception:
            return "", "pdf"
    if "html" in content_type or not content_type:
        parser = TextParser()
        try:
            parser.feed(response.text)
        except Exception:
            return "", "html"
        return parser.text(), "html"
    return "", content_type


def search_results(client: httpx.Client, query: str, limit: int) -> list[tuple[str, int, dict[str, str]]]:
    encoded = quote_plus(query)
    rows: list[tuple[str, int, dict[str, str]]] = []
    try:
        response = client.get(f"https://www.bing.com/search?format=rss&q={encoded}", timeout=25)
        response.raise_for_status()
        rows.extend(("bing_rss", rank, row) for rank, row in enumerate(parse_bing(response.text)[:limit], 1))
    except Exception:
        pass
    try:
        response = client.get(f"https://html.duckduckgo.com/html/?q={encoded}", timeout=25)
        response.raise_for_status()
        rows.extend(("duckduckgo_html", rank, row) for rank, row in enumerate(parse_duck(response.text)[:limit], 1))
    except Exception:
        pass
    return rows


def taxon_regex(query_taxon: str, scope: str) -> re.Pattern[str]:
    if scope == "species":
        parts = query_taxon.split()
        if len(parts) >= 2:
            return re.compile(rf"\b{re.escape(parts[0])}\s+{re.escape(parts[1])}\b", re.I)
    return re.compile(rf"\b{re.escape(genus(query_taxon))}\b", re.I)


def sentences(text: str) -> list[str]:
    return [part.strip() for part in re.split(r"(?<=[.!?;])\s+|\n+", text) if part.strip()]


def extract_evidence(group: str, text: str, query_taxon: str, scope: str) -> list[tuple[str, str]]:
    name_pattern = taxon_regex(query_taxon, scope)
    parts = sentences(text)
    hits: list[tuple[str, str]] = []
    patterns = SI if group == "self_incompatibility" else MATING
    for index, sentence in enumerate(parts):
        low = sentence.casefold()
        values = [
            value
            for value, regexes in patterns.items()
            if any(re.search(regex, low) for regex in regexes)
        ]
        values = list(dict.fromkeys(values))
        if len(values) != 1:
            continue
        # Species/genus can be in the same sentence or the immediately preceding sentence,
        # which covers headings and taxon-led descriptions without treating a whole page as scoped.
        context = sentence
        if index > 0:
            context = f"{parts[index - 1]} {sentence}"
        if not name_pattern.search(context):
            continue
        hits.append((values[0], clean(context)))
    return hits


def one_species(species: str, results_per_query: int, max_pages: int) -> list[dict[str, str]]:
    plan = [(species, "species")]
    genus_name = genus(species)
    if genus_name:
        plan.append((genus_name, "genus"))
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    rows: list[dict[str, str]] = []
    seen_urls: set[str] = set()
    with httpx.Client(headers=headers, follow_redirects=True) as client:
        for query_taxon, scope in plan:
            for group, templates in QUERY_GROUPS.items():
                pages_used = 0
                for template in templates:
                    query = template.format(name=query_taxon)
                    for engine, rank, result in search_results(client, query, results_per_query):
                        url = clean(result.get("url"))
                        if not url or url in seen_urls or pages_used >= max_pages:
                            continue
                        seen_urls.add(url)
                        pages_used += 1
                        text, source_type = fetch_text(client, url)
                        if not text:
                            continue
                        for value, quote in extract_evidence(group, text, query_taxon, scope):
                            out_value = value
                            if group == "self_incompatibility" and scope == "genus":
                                out_value = "likely_SI" if value == "SI" else "likely_SC"
                            rows.append(
                                {
                                    "accepted_species": species,
                                    "trait_group": group,
                                    "value": out_value,
                                    "inference_status": "reported" if scope == "species" else "likely",
                                    "confidence": "high" if scope == "species" else "low",
                                    "query_taxon": query_taxon,
                                    "search_scope": scope,
                                    "search_engine": engine,
                                    "search_query": query,
                                    "search_rank": str(rank),
                                    "source_url": url,
                                    "source_title": clean(result.get("title")),
                                    "evidence_quote": quote,
                                    "source_type": source_type,
                                }
                            )
    return rows


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=200),
    results_per_query: int = typer.Option(5, min=1, max=10),
    max_pages_per_scope_group: int = typer.Option(8, min=1, max=20),
    workers: int = typer.Option(12, min=1, max=24),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    species = [value for value in table["accepted_species"].astype(str).str.strip() if value][:max_taxa]
    rows: list[dict[str, str]] = []
    errors = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(one_species, name, results_per_query, max_pages_per_scope_group): name
            for name in species
        }
        for future in as_completed(futures):
            try:
                rows.extend(future.result())
            except Exception:
                errors += 1
    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(
            ["accepted_species", "trait_group", "value", "inference_status", "source_url", "evidence_quote"]
        )
    coverage = {
        group: int(frame.loc[frame["trait_group"].eq(group), "accepted_species"].nunique())
        if len(frame)
        else 0
        for group in QUERY_GROUPS
    }
    report = {
        "n_species": len(species),
        "coverage": coverage,
        "n_rows": len(frame),
        "n_errors": errors,
        "policy": "full-page species-scoped explicit evidence = reported; genus-scoped explicit evidence = likely",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "reproductive_page_recovery.csv", index=False)
    (output_dir / "reproductive_page_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
