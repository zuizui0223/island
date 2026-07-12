"""Recover low-coverage plant traits from OpenAlex titles and abstracts.

This pilot complements web snippets. Species-level title/abstract matches are treated as
reported candidates when wording is explicit; genus-only matches are low-confidence
likely inferences. The full title/abstract and OpenAlex work id are retained.
"""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import httpx
import pandas as pd
import typer

from island_v2.trait_source_discovery import OPENALEX_WORKS_URL, openalex_abstract, score_taxon_relevance

app = typer.Typer(add_completion=False, no_args_is_help=True)

QUERY_GROUPS = {
    "flower_color": "flower floral color colour corolla petal morphology",
    "pollination": "pollination pollinator floral visitor bee fly bird wind",
    "self_incompatibility": "self-incompatibility self-compatible self-incompatible self-compatibility",
    "mating_system": "mating system selfing outcrossing autogamy mixed mating breeding system",
}

COLOR_PATTERNS = {
    "white": ("white flower", "white flowers", "white corolla", "cream flower"),
    "yellow_orange": ("yellow flower", "yellow flowers", "orange flower", "orange flowers"),
    "red_pink": ("red flower", "red flowers", "pink flower", "pink flowers"),
    "blue_purple": ("blue flower", "blue flowers", "purple flower", "purple flowers", "violet flower"),
}

GUILD_PATTERNS = {
    "bumblebees": ("bumblebee", "bumblebees", "bombus"),
    "bees": ("pollinated by bees", "bee pollination", "bee-pollinated", "bee pollinators"),
    "flies": ("pollinated by flies", "fly pollination", "fly-pollinated"),
    "butterflies": ("butterfly pollination", "pollinated by butterflies"),
    "moths": ("moth pollination", "pollinated by moths", "hawkmoth pollination"),
    "birds": ("bird pollination", "bird-pollinated", "pollinated by birds", "hummingbird pollination"),
    "wind": ("wind pollination", "wind-pollinated", "pollinated by wind", "anemophil"),
    "mixed": ("generalist pollination", "generalist pollinators", "multiple pollinator guilds"),
}

SI_PATTERNS = {
    "SI": (r"\bself[- ]incompatib", r"\bself[- ]sterile\b"),
    "SC": (r"\bself[- ]compatib", r"\bself[- ]fertile\b"),
}

MATING_PATTERNS = {
    "obligate_outcrossing": (r"\bobligate(?:ly)? outcross", r"\bstrict(?:ly)? outcross", r"\bdioecious\b"),
    "mixed_mating": (r"\bmixed[- ]mating\b", r"\bfacultative selfing\b", r"\bselfing and outcrossing\b", r"\bpredominantly outcrossing\b"),
    "mainly_selfing": (r"\bpredominantly selfing\b", r"\bmainly selfing\b", r"\bhigh selfing rate\b", r"\bautogamous\b", r"\bself[- ]pollinating\b", r"\bcleistogamous\b"),
    "obligate_selfing": (r"\bobligate(?:ly)? selfing\b", r"\bstrict(?:ly)? selfing\b"),
}

OUTPUT_COLUMNS = [
    "accepted_species", "trait_group", "value", "inference_status", "confidence",
    "query_taxon", "taxon_relevance_score", "taxon_match_kind", "title", "abstract",
    "openalex_id", "doi", "publication_year", "source_type", "evidence_type",
]


def _genus(species: str) -> str:
    parts = str(species).strip().split()
    return parts[0] if parts else ""


def _extract(group: str, text: str) -> str:
    low = " ".join(str(text).casefold().split())
    if group == "flower_color":
        hits = [value for value, pats in COLOR_PATTERNS.items() if any(pat in low for pat in pats)]
    elif group == "pollination":
        hits = [value for value, pats in GUILD_PATTERNS.items() if any(pat in low for pat in pats)]
    elif group == "self_incompatibility":
        hits = [value for value, pats in SI_PATTERNS.items() if any(re.search(pat, low) for pat in pats)]
    elif group == "mating_system":
        hits = [value for value, pats in MATING_PATTERNS.items() if any(re.search(pat, low) for pat in pats)]
    else:
        hits = []
    hits = list(dict.fromkeys(hits))
    if len(hits) == 1:
        return hits[0]
    if group in {"flower_color", "pollination"} and len(hits) > 1:
        return "multicolored_variable" if group == "flower_color" else "mixed"
    return ""


def _search_taxon(client: httpx.Client, species: str, query_taxon: str, scope: str, per_page: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for group, terms in QUERY_GROUPS.items():
        query = f'"{query_taxon}" {terms}'
        try:
            response = client.get(OPENALEX_WORKS_URL, params={"search": query, "per_page": per_page}, timeout=30)
            response.raise_for_status()
            works = response.json().get("results") or []
        except Exception:
            continue
        for work in works:
            title = str(work.get("display_name") or work.get("title") or "").strip()
            abstract = openalex_abstract(work)
            score, match_kind = score_taxon_relevance(species if scope == "species" else f"{query_taxon} placeholder", title, abstract)
            if scope == "species" and score < 2:
                continue
            if scope == "genus" and query_taxon.casefold() not in f"{title} {abstract}".casefold():
                continue
            value = _extract(group, f"{title}. {abstract}")
            if not value:
                continue
            if group == "self_incompatibility" and scope == "genus":
                value = "likely_SI" if value == "SI" else "likely_SC"
            rows.append(
                {
                    "accepted_species": species,
                    "trait_group": group,
                    "value": value,
                    "inference_status": "reported" if scope == "species" else "likely",
                    "confidence": "high" if scope == "species" and group in {"self_incompatibility", "mating_system"} else ("medium" if scope == "species" else "low"),
                    "query_taxon": query_taxon,
                    "taxon_relevance_score": str(score),
                    "taxon_match_kind": match_kind,
                    "title": title,
                    "abstract": abstract,
                    "openalex_id": str(work.get("id") or ""),
                    "doi": str(work.get("doi") or ""),
                    "publication_year": str(work.get("publication_year") or ""),
                    "source_type": "openalex_title_abstract",
                    "evidence_type": "field_study" if work.get("doi") else "review",
                }
            )
    return rows


def _one_species(species: str, per_page: int) -> list[dict[str, str]]:
    headers = {"User-Agent": "island-floral-v2/0.1 literature trait recovery"}
    rows: list[dict[str, str]] = []
    with httpx.Client(headers=headers) as client:
        rows.extend(_search_taxon(client, species, species, "species", per_page))
        genus = _genus(species)
        if genus:
            rows.extend(_search_taxon(client, species, genus, "genus", per_page))
    return rows


def recover_literature_traits(species_table: pd.DataFrame, max_taxa: int = 25, per_page: int = 8, workers: int = 12) -> tuple[pd.DataFrame, dict[str, object]]:
    species = [value for value in species_table["accepted_species"].astype(str).str.strip() if value][:max_taxa]
    rows: list[dict[str, str]] = []
    errors = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_one_species, sp, per_page): sp for sp in species}
        for future in as_completed(futures):
            try:
                rows.extend(future.result())
            except Exception:
                errors += 1
    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(["accepted_species", "trait_group", "value", "inference_status", "openalex_id"])
    coverage = {
        group: int(frame.loc[frame["trait_group"].eq(group), "accepted_species"].nunique()) if len(frame) else 0
        for group in QUERY_GROUPS
    }
    report: dict[str, object] = {
        "n_species": len(species),
        "coverage": coverage,
        "n_rows": len(frame),
        "n_errors": errors,
        "policy": "explicit species-level title/abstract wording is reported; genus-only wording is low-confidence likely",
    }
    return frame, report


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    per_page: int = typer.Option(8, min=1, max=50),
    workers: int = typer.Option(12, min=1, max=32),
) -> None:
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = recover_literature_traits(species, max_taxa=max_taxa, per_page=per_page, workers=workers)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "literature_trait_recovery.csv", index=False)
    (output_dir / "literature_trait_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
