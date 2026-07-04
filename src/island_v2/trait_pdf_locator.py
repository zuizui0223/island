"""Locate trait-relevant pages in public PDFs without retaining article text.

This is a source-navigation step only. It accepts OpenAlex/Crossref/Unpaywall
leads that declare an open PDF candidate, downloads one bounded PDF at a time to
memory, emits matched keyword page numbers, and discards the bytes immediately.
It never assigns a trait value or changes taxon, establishment, applicability, or
analysis-inclusion status.
"""

from __future__ import annotations

import io
import json
from collections.abc import Callable, Iterable
from typing import Any

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

MAX_PDF_BYTES = 25 * 1024 * 1024
PDF_KEYWORDS = (
    "flower", "floral", "colour", "color", "corolla", "inflorescence", "nectar",
    "pollination", "pollinator", "pollen transfer", "breeding system",
    "self-incompatibility", "self incompatibility", "selfing", "mating system",
    "autonomous selfing", "herkogamy", "dichogamy", "cleistogamy",
)
LOCATOR_COLUMNS = [
    "lead_row", "query_taxon", "island_id", "doi", "source_api", "source_url",
    "pdf_access_status", "http_status", "content_type", "pdf_bytes", "page_count",
    "candidate_pages_json", "locator_status", "locator_note", "do_not_infer",
]
DO_NOT_INFER = (
    "Candidate pages locate keywords in a public PDF only. They do not verify a trait value, "
    "taxon scope, island establishment, Bombus applicability, or analysis inclusion."
)

PdfFetcher = Callable[[str], tuple[int, str, bytes]]
PdfScanner = Callable[[bytes], tuple[int, list[dict[str, Any]]]]


@app.callback()
def main() -> None:
    """Locate candidate trait passages in public open-access PDFs."""


def _text(value: object) -> str:
    return str(value or "").strip()


def public_pdf_candidates(leads: pd.DataFrame) -> pd.DataFrame:
    """Retain only explicitly open-access candidate PDF rows."""
    required = {"candidate_pdf_url", "is_open_access"}
    missing = required.difference(leads.columns)
    if missing:
        raise typer.BadParameter(f"Lead table missing columns: {', '.join(sorted(missing))}")
    table = leads.copy().fillna("")
    open_access = table["is_open_access"].astype(str).str.strip().str.lower().isin({"true", "yes", "1"})
    urls = table["candidate_pdf_url"].astype(str).str.strip().str.startswith(("https://", "http://"))
    return table.loc[open_access & urls].reset_index(names="lead_row")


def fetch_public_pdf(url: str) -> tuple[int, str, bytes]:
    """Fetch one bounded public candidate PDF into memory; callers never persist it."""
    with httpx.Client(timeout=60.0, follow_redirects=True) as client:
        with client.stream("GET", url, headers={"Accept": "application/pdf,*/*"}) as response:
            response.raise_for_status()
            content_type = _text(response.headers.get("content-type")).lower()
            declared = response.headers.get("content-length")
            if declared and int(declared) > MAX_PDF_BYTES:
                raise ValueError(f"response exceeds {MAX_PDF_BYTES} byte limit")
            chunks: list[bytes] = []
            size = 0
            for chunk in response.iter_bytes():
                size += len(chunk)
                if size > MAX_PDF_BYTES:
                    raise ValueError(f"response exceeds {MAX_PDF_BYTES} byte limit")
                chunks.append(chunk)
    return response.status_code, content_type, b"".join(chunks)


def scan_pdf_pages(pdf_bytes: bytes) -> tuple[int, list[dict[str, Any]]]:
    """Return only keyword hit locations; never return PDF text."""
    from pypdf import PdfReader

    reader = PdfReader(io.BytesIO(pdf_bytes))
    hits: list[dict[str, Any]] = []
    for number, page in enumerate(reader.pages, start=1):
        text = (page.extract_text() or "").casefold()
        matched = [keyword for keyword in PDF_KEYWORDS if keyword.casefold() in text]
        if matched:
            hits.append({"page": number, "matched_terms": matched, "score": len(matched)})
    hits.sort(key=lambda row: (-int(row["score"]), int(row["page"])))
    return len(reader.pages), hits[:8]


def _base(row: pd.Series) -> dict[str, str]:
    return {
        "lead_row": str(row.get("lead_row", "")),
        "query_taxon": _text(row.get("query_taxon")),
        "island_id": _text(row.get("island_id")),
        "doi": _text(row.get("doi")),
        "source_api": _text(row.get("source_api")),
        "source_url": _text(row.get("candidate_pdf_url")),
        "pdf_access_status": "declared_public_pdf_candidate",
        "http_status": "",
        "content_type": "",
        "pdf_bytes": "",
        "page_count": "",
        "candidate_pages_json": "[]",
        "locator_status": "not_run",
        "locator_note": "No PDF attempt made.",
        "do_not_infer": DO_NOT_INFER,
    }


def locate_candidates(
    candidates: Iterable[pd.Series],
    *,
    max_pdfs: int,
    fetcher: PdfFetcher = fetch_public_pdf,
    scanner: PdfScanner = scan_pdf_pages,
) -> pd.DataFrame:
    """Download and scan at most ``max_pdfs`` candidates, recording pages only."""
    rows: list[dict[str, str]] = []
    for row in list(candidates)[:max_pdfs]:
        result = _base(row)
        try:
            status, content_type, payload = fetcher(result["source_url"])
            if status >= 400:
                raise RuntimeError(f"HTTP {status}")
            if not payload.startswith(b"%PDF"):
                raise RuntimeError("response is not a PDF")
            page_count, matches = scanner(payload)
            result.update(
                {
                    "pdf_access_status": "public_pdf_recovered",
                    "http_status": str(status),
                    "content_type": content_type,
                    "pdf_bytes": str(len(payload)),
                    "page_count": str(page_count),
                    "candidate_pages_json": json.dumps(matches, ensure_ascii=False),
                    "locator_status": "candidate_pages_located" if matches else "no_keyword_pages_located",
                    "locator_note": "PDF bytes were parsed in memory and discarded; only page locations are retained.",
                }
            )
        except Exception as error:  # noqa: BLE001 - one inaccessible source must not abort batch
            result.update(
                {
                    "pdf_access_status": "public_pdf_not_recovered",
                    "locator_status": "access_or_parse_failed",
                    "locator_note": f"{type(error).__name__}: {error}",
                }
            )
        rows.append(result)
    return pd.DataFrame(rows, columns=LOCATOR_COLUMNS)


@app.command("locate-pages")
def locate_pages(
    leads_csv: str = typer.Option(..., help="trait_source_leads.csv from source discovery."),
    output_csv: str = typer.Option(..., help="Destination trait_pdf_page_locators.csv."),
    max_pdfs: int = typer.Option(20, min=1, max=200, help="Bounded number of public PDF candidates."),
) -> None:
    """Write candidate-page locators from public PDFs without storing articles."""
    leads = pd.read_csv(leads_csv, dtype=str).fillna("")
    candidates = public_pdf_candidates(leads)
    result = locate_candidates((row for _, row in candidates.iterrows()), max_pdfs=max_pdfs)
    destination = pd.io.common.stringify_path(output_csv)
    parent = pd.io.common.get_handle(destination, "w", encoding="utf-8", is_text=True)
    parent.close()
    # Re-open through pathlib only after parent directories exist.
    from pathlib import Path

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(path, index=False)
    typer.echo(f"Wrote {len(result)} public-PDF candidate-page locator rows.")


if __name__ == "__main__":
    app()
