"""Trait-targeted web search for low-coverage reproductive and floral fields.

The generic search lane is good at broad morphology but weak for colour, pollination,
self-incompatibility and mating system. This module issues narrow exact-name queries for
those traits and, when species-level evidence is unavailable, optionally falls back to
same-genus evidence as a low-confidence ``likely`` inference.
"""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote_plus

import httpx
import pandas as pd
import typer

from island_v2.trait_dynamic_discovery import USER_AGENT, _clean, _parse_duck_html, gbif_name_bundle

app = typer.Typer(add_completion=False, no_args_is_help=True)

QUERY_GROUPS = {
    "flower_color": (
        '"{name}" "flower color"',
        '"{name}" "flower colour"',
        '"{name}" flowers corolla petals',
        '"{name}" blooms flower',
    ),
    "pollination": (
        '"{name}" "pollinated by"',
        '"{name}" pollinator pollination',
        '"{name}" bee pollination',
        '"{name}" fly bird moth butterfly pollination',
        '"{name}" "wind pollinated"',
    ),
    "self_incompatibility": (
        '"{name}" "self-compatible"',
        '"{name}" "self-incompatible"',
        '"{name}" "self compatibility"',
        '"{name}" "self incompatibility"',
        '"{name}" "self-fertile" OR "self-sterile"',
    ),
    "mating_system": (
        '"{name}" autogamy selfing',
        '"{name}" outcrossing outcrossed',
        '"{name}" "mixed mating"',
        '"{name}" "mating system"',
        '"{name}" "breeding system"',
    ),
}

COLOR_PATTERNS = {
    "white": ("white", "whitish", "cream", "ivory"),
    "yellow_orange": ("yellow", "orange", "golden"),
    "red_pink": ("red", "pink", "rose", "crimson", "scarlet"),
    "blue_purple": ("blue", "purple", "violet", "lavender"),
    "green_brown_inconspicuous": ("greenish", "brownish", "brown", "green"),
}

GUILD_PATTERNS = {
    "bumblebees": ("bumblebee", "bumblebees"),
    "bees": ("bee-pollinated", "bee pollinated", "pollinated by bees", "bees are the main pollinators"),
    "flies": ("pollinated by flies", "fly-pollinated", "fly pollinated"),
    "butterflies": ("pollinated by butterflies", "butterfly-pollinated"),
    "moths": ("pollinated by moths", "moth-pollinated", "hawkmoth"),
    "birds": ("pollinated by birds", "bird-pollinated", "hummingbird-pollinated"),
    "wind": ("wind-pollinated", "wind pollinated", "pollinated by wind", "anemophil"),
}

SI_PATTERNS = {
    "SI": (r"\bself[- ]incompatible\b", r"\bself[- ]sterile\b", r"\bself incompatibility\b"),
    "SC": (r"\bself[- ]compatible\b", r"\bself[- ]fertile\b", r"\bself compatibility\b"),
}

MATING_PATTERNS = {
    "obligate_outcrossing": (r"\bobligate(?:ly)? outcross", r"\bstrict(?:ly)? outcross"),
    "mixed_mating": (r"\bmixed[- ]mating\b", r"\bboth selfing and outcrossing\b"),
    "mainly_selfing": (r"\bpredominantly selfing\b", r"\bmainly selfing\b", r"\bprimarily selfing\b"),
    "obligate_selfing": (r"\bobligate(?:ly)? selfing\b", r"\bstrict(?:ly)? selfing\b"),
}

OUTPUT_COLUMNS = [
    "accepted_species", "trait_group", "value", "inference_status", "confidence",
    "searched_name", "search_scope", "search_query", "search_rank", "result_title",
    "result_snippet", "result_url", "evidence_type",
]


def _binomial(name: str) -> str:
    parts = _clean(name).split()
    return " ".join(parts[:2]) if len(parts) >= 2 else _clean(name)


def _genus(name: str) -> str:
    parts = _clean(name).split()
    return parts[0] if parts else ""


def _search(client: httpx.Client, query: str, limit: int) -> list[dict[str, str]]:
    response = client.get(
        f"https://html.duckduckgo.com/html/?q={quote_plus(query)}",
        timeout=25,
        follow_redirects=True,
    )
    response.raise_for_status()
    return _parse_duck_html(response.text)[:limit]


def _contains_name(text: str, name: str, scope: str) -> bool:
    low = _clean(text).casefold()
    needle = _binomial(name).casefold() if scope == "species" else _genus(name).casefold()
    return bool(needle) and needle in low


def _extract_value(group: str, text: str) -> str:
    low = _clean(text).casefold()
    if group == "flower_color":
        hits = [value for value, terms in COLOR_PATTERNS.items() if any(re.search(rf"(?<!\w){re.escape(t)}(?!\w)", low) for t in terms)]
        hits = list(dict.fromkeys(hits))
        return hits[0] if len(hits) == 1 else ("multicolored_variable" if len(hits) > 1 else "")
    if group == "pollination":
        hits = [value for value, terms in GUILD_PATTERNS.items() if any(term in low for term in terms)]
        hits = list(dict.fromkeys(hits))
        return hits[0] if len(hits) == 1 else ("mixed" if len(hits) > 1 else "")
    if group == "self_incompatibility":
        hits = [value for value, pats in SI_PATTERNS.items() if any(re.search(pat, low) for pat in pats)]
        hits = list(dict.fromkeys(hits))
        return hits[0] if len(hits) == 1 else ""
    if group == "mating_system":
        hits = [value for value, pats in MATING_PATTERNS.items() if any(re.search(pat, low) for pat in pats)]
        hits = list(dict.fromkeys(hits))
        return hits[0] if len(hits) == 1 else ""
    return ""


def _one_species(species: str, max_synonyms: int, max_results: int) -> list[dict[str, str]]:
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    rows: list[dict[str, str]] = []
    with httpx.Client(headers=headers) as client:
        names = gbif_name_bundle(client, species, max_synonyms=max_synonyms)
        search_plan: list[tuple[str, str]] = [(name, "species") for name in names]
        genus = _genus(species)
        if genus:
            search_plan.append((genus, "genus"))
        seen: set[tuple[str, str, str, str]] = set()
        for searched_name, scope in search_plan:
            for group, templates in QUERY_GROUPS.items():
                for template in templates:
                    query = template.format(name=searched_name)
                    try:
                        results = _search(client, query, max_results)
                    except Exception:
                        continue
                    for rank, result in enumerate(results, start=1):
                        title = _clean(result.get("title"))
                        snippet = _clean(result.get("snippet"))
                        url = _clean(result.get("url"))
                        combined = f"{title}. {snippet}".strip()
                        if not url or not _contains_name(combined, searched_name, scope):
                            continue
                        value = _extract_value(group, combined)
                        if not value:
                            continue
                        key = (group, value, scope, url)
                        if key in seen:
                            continue
                        seen.add(key)
                        if group == "self_incompatibility" and scope == "genus":
                            value = "likely_SI" if value == "SI" else "likely_SC"
                        rows.append({
                            "accepted_species": species,
                            "trait_group": group,
                            "value": value,
                            "inference_status": "reported" if scope == "species" else "likely",
                            "confidence": "high" if scope == "species" and group in {"self_incompatibility", "mating_system"} else ("medium" if scope == "species" else "low"),
                            "searched_name": searched_name,
                            "search_scope": scope,
                            "search_query": query,
                            "search_rank": str(rank),
                            "result_title": title,
                            "result_snippet": snippet,
                            "result_url": url,
                            "evidence_type": "inference" if scope == "genus" else "field_study" if ("journal" in url or "doi.org" in url) else "inference",
                        })
    return rows


def collect_targeted(
    species_table: pd.DataFrame,
    *,
    max_taxa: int = 25,
    max_synonyms: int = 2,
    max_results_per_query: int = 6,
    workers: int = 8,
) -> tuple[pd.DataFrame, dict[str, object]]:
    species = [name for name in species_table["accepted_species"].astype(str).str.strip() if name][:max_taxa]
    all_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_one_species, sp, max_synonyms, max_results_per_query): sp for sp in species}
        for future in as_completed(futures):
            sp = futures[future]
            try:
                all_rows.extend(future.result())
            except Exception as exc:
                errors.append({"accepted_species": sp, "error": str(exc)})
    frame = pd.DataFrame(all_rows, columns=OUTPUT_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(["accepted_species", "trait_group", "value", "search_scope", "result_url"])
    coverage = {
        group: int(frame.loc[frame["trait_group"].eq(group), "accepted_species"].nunique()) if len(frame) else 0
        for group in QUERY_GROUPS
    }
    report: dict[str, object] = {
        "n_species": len(species),
        "coverage": coverage,
        "n_rows": len(frame),
        "n_errors": len(errors),
        "policy": "species-level explicit snippets are reported; same-genus fallback is low-confidence likely inference",
    }
    return frame, report


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_synonyms: int = typer.Option(2, min=0, max=10),
    max_results_per_query: int = typer.Option(6, min=1, max=20),
    workers: int = typer.Option(8, min=1, max=32),
) -> None:
    species = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = collect_targeted(
        species,
        max_taxa=max_taxa,
        max_synonyms=max_synonyms,
        max_results_per_query=max_results_per_query,
        workers=workers,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "targeted_trait_recovery.csv", index=False)
    (output_dir / "targeted_trait_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
