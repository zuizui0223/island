"""Isolated OpenAlex pilot for SC/SI, mating system, pollination and flower colour.

Species-level explicit title/abstract wording is treated as reported evidence. Genus-only
matches are retained as low-confidence likely inference. This pilot is intentionally small
and bounded to the validation species list.
"""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)
OPENALEX_WORKS_URL = "https://api.openalex.org/works"

QUERY_GROUPS = {
    "flower_color": "flower floral color colour corolla petal morphology",
    "pollination": "pollination pollinator floral visitor bee fly bird wind",
    "self_incompatibility": "self-incompatibility self-compatible self-incompatible self-compatibility",
    "mating_system": "mating system selfing outcrossing autogamy mixed mating breeding system",
}

COLOR = {
    "white": ("white flower", "white flowers", "white corolla", "cream flower"),
    "yellow_orange": ("yellow flower", "yellow flowers", "orange flower", "orange flowers"),
    "red_pink": ("red flower", "red flowers", "pink flower", "pink flowers"),
    "blue_purple": ("blue flower", "blue flowers", "purple flower", "purple flowers", "violet flower"),
}

GUILD = {
    "bumblebees": ("bumblebee", "bumblebees", "bombus"),
    "bees": ("pollinated by bees", "bee pollination", "bee-pollinated", "bee pollinators"),
    "flies": ("pollinated by flies", "fly pollination", "fly-pollinated"),
    "butterflies": ("butterfly pollination", "pollinated by butterflies"),
    "moths": ("moth pollination", "pollinated by moths", "hawkmoth pollination"),
    "birds": ("bird pollination", "bird-pollinated", "pollinated by birds", "hummingbird pollination"),
    "wind": ("wind pollination", "wind-pollinated", "pollinated by wind", "anemophil"),
    "mixed": ("generalist pollination", "generalist pollinators", "multiple pollinator guilds"),
}

SI = {
    "SI": (r"\bself[- ]incompatib", r"\bself[- ]sterile\b"),
    "SC": (r"\bself[- ]compatib", r"\bself[- ]fertile\b"),
}

MATING = {
    "obligate_outcrossing": (r"\bobligate(?:ly)? outcross", r"\bstrict(?:ly)? outcross", r"\bdioecious\b"),
    "mixed_mating": (
        r"\bmixed[- ]mating\b",
        r"\bfacultative selfing\b",
        r"\bselfing and outcrossing\b",
        r"\bpredominantly outcrossing\b",
        r"\bmainly outcrossing\b",
    ),
    "mainly_selfing": (
        r"\bpredominantly selfing\b",
        r"\bmainly selfing\b",
        r"\bhigh selfing rate\b",
        r"\bautogamous\b",
        r"\bself[- ]pollinating\b",
        r"\bcleistogamous\b",
    ),
    "obligate_selfing": (r"\bobligate(?:ly)? selfing\b", r"\bstrict(?:ly)? selfing\b"),
}

OUTPUT_COLUMNS = [
    "accepted_species",
    "trait_group",
    "value",
    "inference_status",
    "confidence",
    "query_taxon",
    "search_scope",
    "title",
    "abstract",
    "openalex_id",
    "doi",
    "publication_year",
]


def _abstract(work: dict[str, object]) -> str:
    inverted = work.get("abstract_inverted_index")
    if not isinstance(inverted, dict):
        return ""
    positions: dict[int, str] = {}
    for word, indices in inverted.items():
        if not isinstance(indices, list):
            continue
        for index in indices:
            positions[int(index)] = str(word)
    return " ".join(positions[index] for index in sorted(positions))


def _genus(species: str) -> str:
    parts = species.split()
    return parts[0] if parts else ""


def _taxon_present(text: str, taxon: str, scope: str) -> bool:
    low = text.casefold()
    if scope == "species":
        return taxon.casefold() in low
    return _genus(taxon).casefold() in low


def _extract(group: str, text: str) -> str:
    low = " ".join(text.casefold().split())
    if group == "flower_color":
        hits = [value for value, terms in COLOR.items() if any(term in low for term in terms)]
    elif group == "pollination":
        hits = [value for value, terms in GUILD.items() if any(term in low for term in terms)]
    elif group == "self_incompatibility":
        hits = [value for value, patterns in SI.items() if any(re.search(pattern, low) for pattern in patterns)]
    else:
        hits = [value for value, patterns in MATING.items() if any(re.search(pattern, low) for pattern in patterns)]
    hits = list(dict.fromkeys(hits))
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1 and group == "flower_color":
        return "multicolored_variable"
    if len(hits) > 1 and group == "pollination":
        return "mixed"
    return ""


def _search_scope(
    client: httpx.Client,
    species: str,
    query_taxon: str,
    scope: str,
    per_page: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for group, terms in QUERY_GROUPS.items():
        query = f'"{query_taxon}" {terms}'
        try:
            response = client.get(
                OPENALEX_WORKS_URL,
                params={"search": query, "per_page": per_page},
                timeout=30,
            )
            response.raise_for_status()
            works = response.json().get("results") or []
        except Exception:
            continue
        for work in works:
            title = str(work.get("display_name") or work.get("title") or "").strip()
            abstract = _abstract(work)
            text = f"{title}. {abstract}"
            if not _taxon_present(text, query_taxon, scope):
                continue
            value = _extract(group, text)
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
                    "search_scope": scope,
                    "title": title,
                    "abstract": abstract,
                    "openalex_id": str(work.get("id") or ""),
                    "doi": str(work.get("doi") or ""),
                    "publication_year": str(work.get("publication_year") or ""),
                }
            )
    return rows


def _one_species(species: str, per_page: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    headers = {"User-Agent": "island-floral-v2/0.1 reproductive-literature-pilot"}
    with httpx.Client(headers=headers) as client:
        rows.extend(_search_scope(client, species, species, "species", per_page))
        genus = _genus(species)
        if genus:
            rows.extend(_search_scope(client, species, genus, "genus", per_page))
    return rows


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    per_page: int = typer.Option(8, min=1, max=25),
    workers: int = typer.Option(12, min=1, max=32),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    species = [value for value in table["accepted_species"].astype(str).str.strip() if value][:max_taxa]
    rows: list[dict[str, str]] = []
    errors = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_one_species, name, per_page): name for name in species}
        for future in as_completed(futures):
            try:
                rows.extend(future.result())
            except Exception:
                errors += 1
    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(
            ["accepted_species", "trait_group", "value", "inference_status", "openalex_id"]
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
        "policy": "species explicit title/abstract wording = reported; genus-only wording = likely",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "reproductive_literature_recovery.csv", index=False)
    (output_dir / "reproductive_literature_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
