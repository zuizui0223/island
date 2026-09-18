"""Metadata-only discovery for the prospective H4 temporal validation cohort.

This script intentionally retrieves bibliographic metadata only. It does not request
abstracts, full text, treatment outcomes, pollen-limitation effect sizes, or result
direction. Screening and trait-overlap preflight happen before any outcome extraction.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

QUERIES = (
    "plant pollen limitation",
    "plant pollen supplementation",
    "plant hand pollination fruit set",
    "wild plant supplemental pollination",
    "plant pollen addition reproductive success",
)
START = "2016-01-01"
END = "2026-09-18"


def _get_json(url: str) -> dict:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "island-h4-prospective-metadata/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
        return json.loads(response.read().decode("utf-8"))


def _doi(value: object) -> str:
    text = str(value or "").strip().casefold()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if text.startswith(prefix):
            text = text[len(prefix) :]
    return text


def discover_openalex(max_per_query: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for query in QUERIES:
        url = (
            "https://api.openalex.org/works?"
            + urllib.parse.urlencode(
                {
                    "search": query,
                    "filter": (
                        f"from_publication_date:{START},"
                        f"to_publication_date:{END},type:article"
                    ),
                    "per-page": min(max_per_query, 100),
                }
            )
        )
        payload = _get_json(url)
        for item in payload.get("results", [])[:max_per_query]:
            rows.append(
                {
                    "source_database": "OpenAlex",
                    "source_record_id": str(item.get("id", "")),
                    "doi": _doi(item.get("doi")),
                    "title": str(item.get("title", "") or ""),
                    "publication_date": str(item.get("publication_date", "") or ""),
                    "publication_year": str(item.get("publication_year", "") or ""),
                    "query_family": query,
                }
            )
        time.sleep(0.2)
    return rows


def _crossref_date(item: dict) -> str:
    for key in ("published-print", "published-online", "published", "issued"):
        parts = item.get(key, {}).get("date-parts", [])
        if parts and parts[0]:
            values = list(parts[0])
            year = int(values[0])
            month = int(values[1]) if len(values) > 1 else 1
            day = int(values[2]) if len(values) > 2 else 1
            return f"{year:04d}-{month:02d}-{day:02d}"
    return ""


def discover_crossref(max_per_query: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for query in QUERIES:
        params = urllib.parse.urlencode(
            {
                "query.bibliographic": query,
                "filter": (
                    f"from-pub-date:{START},until-pub-date:{END},"
                    "type:journal-article"
                ),
                "rows": min(max_per_query, 100),
                "select": "DOI,title,published,published-print,published-online,issued",
            }
        )
        payload = _get_json("https://api.crossref.org/works?" + params)
        for item in payload.get("message", {}).get("items", [])[:max_per_query]:
            title_values = item.get("title") or [""]
            date = _crossref_date(item)
            rows.append(
                {
                    "source_database": "Crossref",
                    "source_record_id": _doi(item.get("DOI")),
                    "doi": _doi(item.get("DOI")),
                    "title": str(title_values[0] if title_values else ""),
                    "publication_date": date,
                    "publication_year": date[:4] if date else "",
                    "query_family": query,
                }
            )
        time.sleep(0.2)
    return rows


def deduplicate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    chosen: dict[str, dict[str, str]] = {}
    for row in rows:
        date = row["publication_date"].strip()
        if not date or date < START or date > END:
            continue
        doi = row["doi"].strip()
        key = (
            f"doi:{doi}"
            if doi
            else "title:" + row["title"].strip().casefold() + "|" + row["publication_year"]
        )
        current = chosen.get(key)
        if current is None:
            chosen[key] = dict(row)
            continue
        sources = sorted(
            {
                token
                for value in (
                    current.get("source_database", ""),
                    row.get("source_database", ""),
                )
                for token in str(value).split("|")
                if token
            }
        )
        current["source_database"] = "|".join(sources)
        queries = sorted(
            {
                token
                for value in (
                    current.get("query_family", ""),
                    row.get("query_family", ""),
                )
                for token in str(value).split("|")
                if token
            }
        )
        current["query_family"] = "|".join(queries)
    return sorted(
        chosen.values(),
        key=lambda row: (row["publication_date"], row["doi"], row["title"]),
    )


@app.command("run")
def run(
    output_csv: Path = typer.Option(...),
    max_per_query: int = typer.Option(100, min=1, max=100),
) -> None:
    rows = deduplicate(
        discover_openalex(max_per_query) + discover_crossref(max_per_query)
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_database",
        "source_record_id",
        "doi",
        "title",
        "publication_date",
        "publication_year",
        "query_family",
    ]
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    typer.echo(json.dumps({"n_metadata_candidates": len(rows)}, indent=2))


if __name__ == "__main__":
    app()
