"""Dynamic, query-driven public web discovery for plant trait evidence.

Unlike fixed-site lookup, this lane searches the public web with trait-aware queries,
expands accepted names with GBIF synonyms, freezes search-result provenance, fetches
candidate pages, and writes species-direct source text for downstream deterministic
or LLM extraction.

No trait value is decided here. The output is an evidence corpus.
"""

from __future__ import annotations

import hashlib
import io
import json
import re
import time
from collections.abc import Iterable
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import parse_qs, quote_plus, unquote, urlparse
from xml.etree import ElementTree

import httpx
import pandas as pd
import typer
from pypdf import PdfReader

app = typer.Typer(add_completion=False, no_args_is_help=True)

USER_AGENT = "island-floral-v2/0.1 reproducible trait evidence discovery"
SEARCH_ENDPOINTS = (
    "https://www.bing.com/search?format=rss&q={query}",
    "https://html.duckduckgo.com/html/?q={query}",
)
QUERY_FAMILIES = {
    "floral_morphology": (
        '"{name}" flower color morphology corolla',
        '"{name}" flora flower description',
    ),
    "reproductive_ecology": (
        '"{name}" pollination breeding system self compatibility',
        '"{name}" reproductive biology pollinator',
    ),
}
SOURCE_COLUMNS = [
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
]
SEARCH_COLUMNS = [
    "accepted_species",
    "searched_name",
    "query_family",
    "search_query",
    "search_engine",
    "search_rank",
    "result_title",
    "result_url",
    "result_snippet",
]


class _TextParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []
        self._hidden = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() in {"script", "style", "noscript", "svg"}:
            self._hidden += 1

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() in {"script", "style", "noscript", "svg"} and self._hidden:
            self._hidden -= 1

    def handle_data(self, data: str) -> None:
        if not self._hidden:
            self.parts.append(data)

    def text(self) -> str:
        return _clean(" ".join(self.parts))


class _DuckParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.results: list[dict[str, str]] = []
        self._in_link = False
        self._in_snippet = False
        self._href = ""
        self._title_parts: list[str] = []
        self._snippet_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attrs_dict = dict(attrs)
        classes = set((attrs_dict.get("class") or "").split())
        if tag == "a" and "result__a" in classes:
            self._in_link = True
            self._href = attrs_dict.get("href") or ""
            self._title_parts = []
        if "result__snippet" in classes:
            self._in_snippet = True
            self._snippet_parts = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "a" and self._in_link:
            self._in_link = False
            url = _unwrap_duck_url(self._href)
            if url:
                self.results.append(
                    {"title": _clean(" ".join(self._title_parts)), "url": url, "snippet": ""}
                )
        if self._in_snippet and tag in {"a", "div", "span"}:
            self._in_snippet = False
            if self.results:
                self.results[-1]["snippet"] = _clean(" ".join(self._snippet_parts))

    def handle_data(self, data: str) -> None:
        if self._in_link:
            self._title_parts.append(data)
        if self._in_snippet:
            self._snippet_parts.append(data)


def _clean(value: object) -> str:
    return " ".join(str(value or "").replace("\x00", " ").split())


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _unwrap_duck_url(url: str) -> str:
    if not url:
        return ""
    parsed = urlparse(url)
    if "duckduckgo.com" in parsed.netloc and parsed.path.startswith("/l/"):
        target = parse_qs(parsed.query).get("uddg", [""])[0]
        return unquote(target)
    if url.startswith("//"):
        return "https:" + url
    return url


def _parse_bing_rss(text: str) -> list[dict[str, str]]:
    try:
        root = ElementTree.fromstring(text)
    except ElementTree.ParseError:
        return []
    rows = []
    for item in root.findall(".//item"):
        rows.append(
            {
                "title": _clean(item.findtext("title")),
                "url": _clean(item.findtext("link")),
                "snippet": _clean(item.findtext("description")),
            }
        )
    return [row for row in rows if row["url"]]


def _parse_duck_html(text: str) -> list[dict[str, str]]:
    parser = _DuckParser()
    parser.feed(text)
    return parser.results


def _source_type(url: str) -> str:
    host = urlparse(url).netloc.lower()
    path = urlparse(url).path.lower()
    if path.endswith(".pdf"):
        return "paper_or_flora_pdf"
    if any(token in host for token in ("flora", "efloras", "powo", "kew", "nzpcn")):
        return "flora_or_monograph"
    if any(token in host for token in ("edu", "ac.", "gov", "org")):
        return "institutional_or_literature_web"
    if "wikipedia" in host:
        return "encyclopedia"
    return "public_web_page"


def _name_tokens(name: str) -> tuple[str, str]:
    parts = _clean(name).lower().split()
    if len(parts) < 2:
        return (parts[0] if parts else "", "")
    return parts[0], parts[1]


def _mentions_name(text: str, names: Iterable[str]) -> bool:
    low = text.lower()
    for name in names:
        genus, epithet = _name_tokens(name)
        if genus and epithet and (f"{genus} {epithet}" in low or (genus in low and epithet in low)):
            return True
    return False


def _extract_pdf(content: bytes, max_pages: int = 40) -> str:
    try:
        reader = PdfReader(io.BytesIO(content))
    except Exception:
        return ""
    parts = []
    for page in reader.pages[:max_pages]:
        try:
            parts.append(page.extract_text() or "")
        except Exception:
            continue
    return _clean(" ".join(parts))


def _extract_html(text: str) -> str:
    parser = _TextParser()
    parser.feed(text)
    return parser.text()


def gbif_name_bundle(client: httpx.Client, accepted_name: str, max_synonyms: int = 10) -> list[str]:
    names = [accepted_name]
    try:
        match = client.get(
            "https://api.gbif.org/v1/species/match",
            params={"name": accepted_name},
            timeout=20,
        ).json()
        key = match.get("usageKey") or match.get("acceptedUsageKey")
        canonical = _clean(match.get("canonicalName"))
        if canonical and canonical not in names:
            names.append(canonical)
        if key:
            response = client.get(
                f"https://api.gbif.org/v1/species/{key}/synonyms",
                params={"limit": max_synonyms},
                timeout=20,
            )
            data = response.json()
            rows = data.get("results", data if isinstance(data, list) else [])
            for row in rows[:max_synonyms]:
                synonym = _clean(row.get("canonicalName") or row.get("scientificName"))
                if synonym and synonym not in names:
                    names.append(synonym)
    except Exception:
        pass
    return names[: max_synonyms + 1]


def search_web(client: httpx.Client, query: str, max_results: int = 8) -> tuple[str, list[dict[str, str]]]:
    encoded = quote_plus(query)
    for template in SEARCH_ENDPOINTS:
        url = template.format(query=encoded)
        engine = "bing_rss" if "bing.com" in url else "duckduckgo_html"
        try:
            response = client.get(url, timeout=25, follow_redirects=True)
            response.raise_for_status()
            results = (
                _parse_bing_rss(response.text)
                if engine == "bing_rss"
                else _parse_duck_html(response.text)
            )
            if results:
                return engine, results[:max_results]
        except Exception:
            continue
    return "none", []


def fetch_source(client: httpx.Client, url: str, max_chars: int = 60000) -> str:
    try:
        response = client.get(url, timeout=30, follow_redirects=True)
        response.raise_for_status()
    except Exception:
        return ""
    content_type = response.headers.get("content-type", "").lower()
    if "pdf" in content_type or urlparse(str(response.url)).path.lower().endswith(".pdf"):
        text = _extract_pdf(response.content)
    else:
        text = _extract_html(response.text)
    return text[:max_chars]


def discover_species(
    client: httpx.Client,
    accepted_species: str,
    *,
    max_synonyms: int = 5,
    max_results_per_query: int = 6,
    max_pages_per_species: int = 12,
    delay_seconds: float = 0.15,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    names = gbif_name_bundle(client, accepted_species, max_synonyms=max_synonyms)
    search_rows: list[dict[str, str]] = []
    candidate_rows: list[dict[str, str]] = []
    seen_urls: set[str] = set()

    for searched_name in names:
        for family, templates in QUERY_FAMILIES.items():
            for template in templates:
                query = template.format(name=searched_name)
                engine, results = search_web(client, query, max_results=max_results_per_query)
                for rank, result in enumerate(results, start=1):
                    url = _clean(result.get("url"))
                    if not url:
                        continue
                    search_rows.append(
                        {
                            "accepted_species": accepted_species,
                            "searched_name": searched_name,
                            "query_family": family,
                            "search_query": query,
                            "search_engine": engine,
                            "search_rank": str(rank),
                            "result_title": _clean(result.get("title")),
                            "result_url": url,
                            "result_snippet": _clean(result.get("snippet")),
                        }
                    )
                    if url in seen_urls or len(seen_urls) >= max_pages_per_species:
                        continue
                    seen_urls.add(url)
                    text = fetch_source(client, url)
                    time.sleep(delay_seconds)
                    if len(text) < 120 or not _mentions_name(text, names):
                        continue
                    candidate_rows.append(
                        {
                            "accepted_species": accepted_species,
                            "searched_name": searched_name,
                            "query_family": family,
                            "search_query": query,
                            "search_engine": engine,
                            "search_rank": str(rank),
                            "search_result_title": _clean(result.get("title")),
                            "search_result_snippet": _clean(result.get("snippet")),
                            "source_text": text,
                            "source_url": url,
                            "source_citation": _clean(result.get("title")) or url,
                            "source_type": _source_type(url),
                            "evidence_scope": "species_direct",
                            "content_sha256": _sha256(text),
                        }
                    )
    return search_rows, candidate_rows


def run_discovery(
    species: Iterable[str],
    *,
    max_synonyms: int = 5,
    max_results_per_query: int = 6,
    max_pages_per_species: int = 12,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    all_search: list[dict[str, str]] = []
    all_sources: list[dict[str, str]] = []
    with httpx.Client(headers=headers) as client:
        for accepted_species in species:
            search_rows, source_rows = discover_species(
                client,
                accepted_species,
                max_synonyms=max_synonyms,
                max_results_per_query=max_results_per_query,
                max_pages_per_species=max_pages_per_species,
            )
            all_search.extend(search_rows)
            all_sources.extend(source_rows)
    search_frame = pd.DataFrame(all_search, columns=SEARCH_COLUMNS).drop_duplicates()
    source_frame = pd.DataFrame(all_sources, columns=SOURCE_COLUMNS).drop_duplicates(
        ["accepted_species", "source_url", "content_sha256"]
    )
    summary = {
        "n_species_requested": len(list(species)) if not isinstance(species, list) else len(species),
        "n_species_with_search_results": int(search_frame["accepted_species"].nunique()) if len(search_frame) else 0,
        "n_species_with_fetched_sources": int(source_frame["accepted_species"].nunique()) if len(source_frame) else 0,
        "n_search_results": len(search_frame),
        "n_fetched_sources": len(source_frame),
        "source_type_counts": source_frame["source_type"].value_counts().to_dict() if len(source_frame) else {},
    }
    return search_frame, source_frame, summary


@app.command("discover")
def discover_command(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    species_column: str = typer.Option("accepted_species"),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_synonyms: int = typer.Option(5, min=0, max=20),
    max_results_per_query: int = typer.Option(6, min=1, max=20),
    max_pages_per_species: int = typer.Option(12, min=1, max=50),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    if species_column not in table.columns:
        raise typer.BadParameter(f"missing species column: {species_column}")
    species = []
    seen = set()
    for value in table[species_column]:
        name = _clean(value)
        if name and name not in seen:
            seen.add(name)
            species.append(name)
        if len(species) >= max_taxa:
            break
    search_frame, source_frame, summary = run_discovery(
        species,
        max_synonyms=max_synonyms,
        max_results_per_query=max_results_per_query,
        max_pages_per_species=max_pages_per_species,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    search_frame.to_csv(output_dir / "search_results.csv", index=False)
    source_frame.to_csv(output_dir / "source_texts.csv", index=False)
    (output_dir / "dynamic_discovery_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
