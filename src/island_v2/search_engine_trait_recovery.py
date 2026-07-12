"""Scalable search-engine-only plant trait recovery.

Designed for large taxon sets where literature APIs and full-page retrieval are too slow.
The primary pass uses only exact scientific-name search-result snippets. A second fallback
pass is restricted to unresolved fields and may use the genus name as low-confidence
``likely`` evidence. All queries and snippets are retained for audit and resumability.
"""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import parse_qs, quote_plus, unquote, urlparse
from xml.etree import ElementTree

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

USER_AGENT = "island-floral-v2/0.1 search-engine-trait-recovery"

# Short first-pass queries work better than keyword-stuffed queries in search snippets.
PRIMARY_QUERIES = {
    "phenotype": '"{name}" flower',
    "reproduction": '"{name}" pollination',
}

# Only unresolved fields trigger these short exact-phrase variants. The expensive
# variants are gap-only, so well-described taxa still finish after two searches.
FALLBACK_QUERY_VARIANTS = {
    "flower_color": (
        '"{name}" flower color',
        '"{name}" flower colour',
        '"{name}" flowers are',
    ),
    "flower_shape": (
        '"{name}" flower shape',
        '"{name}" floral morphology',
        '"{name}" corolla',
    ),
    "pollination_guild": (
        '"{name}" pollinated by',
        '"{name}" pollinator',
        '"{name}" pollination biology',
    ),
    "mating_system": (
        '"{name}" mating system',
        '"{name}" breeding system',
        '"{name}" selfing outcrossing',
        '"{name}" autogamy',
    ),
    "self_incompatibility": (
        '"{name}" self-compatible',
        '"{name}" self-incompatible',
        '"{name}" self-fertile',
        '"{name}" self-sterile',
    ),
}

COLOR_TERMS = {
    "white": ("white", "whitish", "cream", "ivory"),
    "yellow_orange": ("yellow", "orange", "golden"),
    "red_pink": ("red", "pink", "rose", "crimson", "scarlet"),
    "blue_purple": ("blue", "purple", "violet", "lavender"),
    "green_brown_inconspicuous": ("greenish", "brownish", "brown", "green flowers"),
}
FLOWER_CONTEXT = ("flower", "flowers", "floral", "corolla", "petal", "petals", "bloom", "blooms")

SHAPE_PATTERNS = {
    "bell_campanulate": (r"\bcampanulat", r"\bbell[- ]shaped\b"),
    "tubular": (r"\btubular\b", r"\btube[- ]shaped\b", r"\blong[- ]tubed\b"),
    "funnel_trumpet": (r"\bfunnel[- ]shaped\b", r"\btrumpet[- ]shaped\b", r"\binfundibuliform\b"),
    "open_radial": (r"\brotate\b", r"\bopen flower", r"\bflat flower"),
    "composite_head": (r"\bcapitulum\b", r"\bcapitula\b", r"\bflower head\b"),
    "reduced_wind": (r"\bspikelet", r"\bgrass flower", r"\bsedge flower"),
    "zygomorphic": (r"\bzygomorph", r"\bbilateral(?:ly)? symmetric\b"),
    "actinomorphic": (r"\bactinomorph", r"\bradial(?:ly)? symmetric\b"),
}

GUILD_PATTERNS = {
    "bumblebees": (r"\bbumblebee", r"\bBombus\b"),
    "bees": (r"\bpollinated by bees\b", r"\bbee[- ]pollinat", r"\bbees? (?:are )?(?:the )?(?:main )?pollinators?\b"),
    "flies": (r"\bpollinated by flies\b", r"\bfly[- ]pollinat"),
    "butterflies": (r"\bpollinated by butterflies\b", r"\bbutterfly[- ]pollinat"),
    "moths": (r"\bpollinated by moths\b", r"\bmoth[- ]pollinat", r"\bhawkmoth"),
    "birds": (r"\bpollinated by birds\b", r"\bbird[- ]pollinat", r"\bhummingbird[- ]pollinat"),
    "wind": (r"\bwind[- ]pollinat", r"\bpollinated by wind\b", r"\banemophil"),
    "self": (r"\bautonomous selfing\b", r"\bself[- ]pollinat(?:ing|ed)\b", r"\bautogamous\b"),
    "mixed": (r"\bmixed pollination\b", r"\bgeneralist pollinat", r"\binsect[- ]pollinat", r"\bentomophil"),
}

SI_PATTERNS = {
    "SI": (r"\bself[- ]incompatib", r"\bself[- ]sterile\b"),
    "SC": (r"\bself[- ]compatib", r"\bself[- ]fertile\b"),
}

MATING_PATTERNS = {
    "obligate_outcrossing": (
        r"\bobligate(?:ly)? outcross",
        r"\bstrict(?:ly)? outcross",
        r"\bself-incompatible and outcross",
    ),
    "mixed_mating": (
        r"\bmixed[- ]mating\b",
        r"\bselfing and outcrossing\b",
        r"\bboth selfing and outcrossing\b",
        r"\bfacultative selfing\b",
    ),
    "mainly_selfing": (
        r"\bpredominantly selfing\b",
        r"\bmainly selfing\b",
        r"\bprimarily selfing\b",
        r"\bhigh selfing rate\b",
        r"\bautogamous\b",
    ),
    "obligate_selfing": (r"\bobligate(?:ly)? selfing\b", r"\bstrict(?:ly)? selfing\b"),
}

TRAITS = (
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "mating_system",
    "self_incompatibility",
)

EVIDENCE_COLUMNS = [
    "accepted_species",
    "trait",
    "value",
    "inference_status",
    "confidence",
    "query_taxon",
    "search_scope",
    "query_family",
    "search_engine",
    "search_query",
    "search_rank",
    "result_title",
    "result_snippet",
    "result_url",
]

TABLE_COLUMNS = [
    "species",
    "flower_color",
    "flower_shape",
    "pollination_guild",
    "pollination_notes",
    "mating_system",
    "self_incompatibility",
    "evidence_type",
    "confidence",
]


class DuckParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.results: list[dict[str, str]] = []
        self.in_link = False
        self.in_snippet = False
        self.href = ""
        self.title_parts: list[str] = []
        self.snippet_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        classes = set((values.get("class") or "").split())
        if tag == "a" and "result__a" in classes:
            self.in_link = True
            self.href = values.get("href") or ""
            self.title_parts = []
        if "result__snippet" in classes:
            self.in_snippet = True
            self.snippet_parts = []

    def handle_endtag(self, tag: str) -> None:
        if tag == "a" and self.in_link:
            self.in_link = False
            url = unwrap_duck_url(self.href)
            if url:
                self.results.append({"title": clean(" ".join(self.title_parts)), "url": url, "snippet": ""})
        if self.in_snippet and tag in {"a", "div", "span"}:
            self.in_snippet = False
            if self.results:
                self.results[-1]["snippet"] = clean(" ".join(self.snippet_parts))

    def handle_data(self, data: str) -> None:
        if self.in_link:
            self.title_parts.append(data)
        if self.in_snippet:
            self.snippet_parts.append(data)


def clean(value: object) -> str:
    return " ".join(str(value or "").replace("\x00", " ").split())


def unwrap_duck_url(url: str) -> str:
    if not url:
        return ""
    parsed = urlparse(url)
    if "duckduckgo.com" in parsed.netloc and parsed.path.startswith("/l/"):
        return unquote(parse_qs(parsed.query).get("uddg", [""])[0])
    if url.startswith("//"):
        return "https:" + url
    return url


def parse_duck(text: str) -> list[dict[str, str]]:
    parser = DuckParser()
    parser.feed(text)
    return parser.results


def parse_bing(text: str) -> list[dict[str, str]]:
    try:
        root = ElementTree.fromstring(text)
    except ElementTree.ParseError:
        return []
    rows: list[dict[str, str]] = []
    for item in root.findall(".//item"):
        rows.append(
            {
                "title": clean(item.findtext("title")),
                "url": clean(item.findtext("link")),
                "snippet": clean(item.findtext("description")),
            }
        )
    return [row for row in rows if row["url"]]


def genus(species: str) -> str:
    parts = clean(species).split()
    return parts[0] if parts else ""


def taxon_present(text: str, taxon: str, scope: str) -> bool:
    low = clean(text).casefold()
    needle = clean(taxon).casefold()
    if scope == "species":
        return needle in low
    return bool(genus(taxon)) and genus(taxon).casefold() in low


def extract_traits(text: str) -> dict[str, str]:
    low = clean(text).casefold()
    values: dict[str, str] = {}

    if any(term in low for term in FLOWER_CONTEXT):
        hits = [
            value
            for value, terms in COLOR_TERMS.items()
            if any(re.search(rf"(?<!\w){re.escape(term)}(?!\w)", low) for term in terms)
        ]
        hits = list(dict.fromkeys(hits))
        if len(hits) == 1:
            values["flower_color"] = hits[0]
        elif len(hits) > 1:
            values["flower_color"] = "multicolored_variable"

        shape_hits = [
            value
            for value, patterns in SHAPE_PATTERNS.items()
            if any(re.search(pattern, low, flags=re.IGNORECASE) for pattern in patterns)
        ]
        shape_hits = list(dict.fromkeys(shape_hits))
        if len(shape_hits) == 1:
            values["flower_shape"] = shape_hits[0]

    guild_hits = [
        value
        for value, patterns in GUILD_PATTERNS.items()
        if any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)
    ]
    guild_hits = list(dict.fromkeys(guild_hits))
    if len(guild_hits) == 1:
        values["pollination_guild"] = guild_hits[0]
    elif len(guild_hits) > 1:
        values["pollination_guild"] = "mixed"

    si_hits = [
        value
        for value, patterns in SI_PATTERNS.items()
        if any(re.search(pattern, low) for pattern in patterns)
    ]
    si_hits = list(dict.fromkeys(si_hits))
    if len(si_hits) == 1:
        values["self_incompatibility"] = si_hits[0]

    mating_hits = [
        value
        for value, patterns in MATING_PATTERNS.items()
        if any(re.search(pattern, low) for pattern in patterns)
    ]
    mating_hits = list(dict.fromkeys(mating_hits))
    if len(mating_hits) == 1:
        values["mating_system"] = mating_hits[0]

    return values


def search_bing(client: httpx.Client, query: str, limit: int) -> list[tuple[str, dict[str, str]]]:
    encoded = quote_plus(query)
    response = client.get(f"https://www.bing.com/search?format=rss&q={encoded}", timeout=20)
    response.raise_for_status()
    return [("bing_rss", row) for row in parse_bing(response.text)[:limit]]


def search_duck(client: httpx.Client, query: str, limit: int) -> list[tuple[str, dict[str, str]]]:
    encoded = quote_plus(query)
    response = client.get(f"https://html.duckduckgo.com/html/?q={encoded}", timeout=20)
    response.raise_for_status()
    return [("duckduckgo_html", row) for row in parse_duck(response.text)[:limit]]


def search_with_fallback(client: httpx.Client, query: str, limit: int) -> list[tuple[str, dict[str, str]]]:
    try:
        rows = search_bing(client, query, limit)
        if rows:
            return rows
    except Exception:
        pass
    try:
        return search_duck(client, query, limit)
    except Exception:
        return []


def collect_query(
    client: httpx.Client,
    *,
    species: str,
    query_taxon: str,
    scope: str,
    query_family: str,
    query: str,
    limit: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rank, (engine, result) in enumerate(search_with_fallback(client, query, limit), start=1):
        title = clean(result.get("title"))
        snippet = clean(result.get("snippet"))
        url = clean(result.get("url"))
        combined = f"{title}. {snippet}".strip()
        if not url:
            continue
        explicit_taxon = taxon_present(combined, query_taxon, scope)
        # Search engines often omit the queried scientific name from the rendered
        # snippet. Exact-name query provenance is therefore retained as low-confidence
        # likely evidence instead of being discarded outright.
        status = "reported" if scope == "species" and explicit_taxon else "likely"
        confidence = "medium" if status == "reported" else "low"
        target_trait = ""
        if query_family.startswith("gap_"):
            target_trait = query_family.removeprefix("gap_")
        elif query_family.startswith("genus_"):
            target_trait = query_family.removeprefix("genus_")
        for trait, value in extract_traits(combined).items():
            if target_trait and trait != target_trait:
                continue
            if status == "likely" and trait == "self_incompatibility":
                value = "likely_SI" if value == "SI" else "likely_SC"
            rows.append(
                {
                    "accepted_species": species,
                    "trait": trait,
                    "value": value,
                    "inference_status": status,
                    "confidence": confidence,
                    "query_taxon": query_taxon,
                    "search_scope": scope,
                    "query_family": query_family,
                    "search_engine": engine,
                    "search_query": query,
                    "search_rank": str(rank),
                    "result_title": title,
                    "result_snippet": snippet,
                    "result_url": url,
                }
            )
    return rows


def best_values(evidence: pd.DataFrame, species: str) -> dict[str, str]:
    subset = evidence.loc[evidence["accepted_species"].eq(species)].copy()
    values: dict[str, str] = {}
    rank_status = {"reported": 0, "likely": 1}
    rank_conf = {"high": 0, "medium": 1, "low": 2}
    for trait in TRAITS:
        current = subset.loc[subset["trait"].eq(trait)].copy()
        if current.empty:
            continue
        current["_s"] = current["inference_status"].map(rank_status).fillna(9)
        current["_c"] = current["confidence"].map(rank_conf).fillna(9)
        current = current.sort_values(["_s", "_c", "search_rank"], kind="stable")
        values[trait] = clean(current.iloc[0]["value"])
    return values


def one_species(species: str, result_limit: int, use_genus_fallback: bool) -> tuple[list[dict[str, str]], dict[str, str]]:
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    evidence: list[dict[str, str]] = []
    with httpx.Client(headers=headers, follow_redirects=True) as client:
        for family, template in PRIMARY_QUERIES.items():
            query = template.format(name=species)
            evidence.extend(
                collect_query(
                    client,
                    species=species,
                    query_taxon=species,
                    scope="species",
                    query_family=family,
                    query=query,
                    limit=result_limit,
                )
            )

        frame = pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS)
        resolved = set(frame["trait"]) if not frame.empty else set()

        # Short exact species-level variants only for unresolved fields. Stop as soon
        # as a variant resolves the target trait.
        for trait in TRAITS:
            if trait in resolved:
                continue
            for template in FALLBACK_QUERY_VARIANTS[trait]:
                query = template.format(name=species)
                evidence.extend(
                    collect_query(
                        client,
                        species=species,
                        query_taxon=species,
                        scope="species",
                        query_family=f"gap_{trait}",
                        query=query,
                        limit=result_limit,
                    )
                )
                if any(row["trait"] == trait for row in evidence):
                    break

        frame = pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS)
        resolved = set(frame["trait"]) if not frame.empty else set()

        if use_genus_fallback:
            genus_name = genus(species)
            for trait in TRAITS:
                if trait in resolved or not genus_name:
                    continue
                for template in FALLBACK_QUERY_VARIANTS[trait]:
                    query = template.format(name=genus_name)
                    evidence.extend(
                        collect_query(
                            client,
                            species=species,
                            query_taxon=genus_name,
                            scope="genus",
                            query_family=f"genus_{trait}",
                            query=query,
                            limit=result_limit,
                        )
                    )
                    if any(row["trait"] == trait for row in evidence):
                        break

    frame = pd.DataFrame(evidence, columns=EVIDENCE_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(
            ["accepted_species", "trait", "value", "inference_status", "result_url"]
        )
    values = best_values(frame, species) if not frame.empty else {}
    return frame.to_dict("records"), values


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=100000),
    results_per_query: int = typer.Option(8, min=1, max=20),
    workers: int = typer.Option(24, min=1, max=64),
    use_genus_fallback: bool = typer.Option(True, "--genus-fallback/--no-genus-fallback"),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    species_list = [
        clean(value)
        for value in table["accepted_species"].astype(str)
        if clean(value)
    ][:max_taxa]

    all_evidence: list[dict[str, str]] = []
    compact_rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(one_species, species, results_per_query, use_genus_fallback): species
            for species in species_list
        }
        for future in as_completed(futures):
            species = futures[future]
            try:
                rows, values = future.result()
                all_evidence.extend(rows)
                compact_rows.append(
                    {
                        "species": species,
                        "flower_color": values.get("flower_color", "unknown"),
                        "flower_shape": values.get("flower_shape", "unknown"),
                        "pollination_guild": values.get("pollination_guild", "unknown"),
                        "pollination_notes": "search-engine snippet evidence retained in detailed output",
                        "mating_system": values.get("mating_system", "unknown"),
                        "self_incompatibility": values.get("self_incompatibility", "unknown"),
                        "evidence_type": "inference",
                        "confidence": "medium" if any(values.values()) else "low",
                    }
                )
            except Exception as exc:
                errors.append({"species": species, "error": str(exc)})
                compact_rows.append(
                    {
                        "species": species,
                        "flower_color": "unknown",
                        "flower_shape": "unknown",
                        "pollination_guild": "unknown",
                        "pollination_notes": "unknown",
                        "mating_system": "unknown",
                        "self_incompatibility": "unknown",
                        "evidence_type": "inference",
                        "confidence": "low",
                    }
                )

    evidence = pd.DataFrame(all_evidence, columns=EVIDENCE_COLUMNS)
    compact = pd.DataFrame(compact_rows, columns=TABLE_COLUMNS)
    order = {name: i for i, name in enumerate(species_list)}
    compact["_order"] = compact["species"].map(order)
    compact = compact.sort_values("_order").drop(columns="_order")

    coverage = {
        trait: int(compact[trait].ne("unknown").sum())
        for trait in TRAITS
    }
    report = {
        "n_species": len(species_list),
        "coverage": coverage,
        "n_species_with_any_non_unknown_trait": int(
            compact[list(TRAITS)].ne("unknown").any(axis=1).sum()
        ),
        "n_evidence_rows": len(evidence),
        "n_errors": len(errors),
        "query_policy": (
            "2 short exact-name queries per species; unresolved fields trigger short exact-phrase variants; "
            "query-bound snippets without repeated taxon names are retained as low-confidence likely evidence; "
            "remaining gaps optionally trigger genus-level likely fallback; Bing RSS primary, DuckDuckGo fallback"
        ),
        "page_fetching": False,
        "literature_api": False,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output_dir / "search_engine_trait_evidence.csv", index=False)
    compact.to_csv(output_dir / "search_engine_trait_table.csv", index=False)
    pd.DataFrame(errors, columns=["species", "error"]).to_csv(
        output_dir / "search_engine_trait_errors.csv", index=False
    )
    (output_dir / "search_engine_trait_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
