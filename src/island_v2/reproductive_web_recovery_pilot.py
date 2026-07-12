"""Isolated search-engine pilot for low-coverage plant reproductive traits."""

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
USER_AGENT = "island-floral-v2/0.1 reproductive-web-pilot"

QUERY_GROUPS = {
    "flower_color": (
        '"{name}" "flower color"',
        '"{name}" "flower colour"',
        '"{name}" flowers corolla petals',
    ),
    "pollination": (
        '"{name}" "pollinated by"',
        '"{name}" pollinator pollination',
        '"{name}" "wind pollinated"',
    ),
    "self_incompatibility": (
        '"{name}" "self-compatible"',
        '"{name}" "self-incompatible"',
        '"{name}" "self compatibility"',
        '"{name}" "self incompatibility"',
    ),
    "mating_system": (
        '"{name}" autogamy autogamous selfing',
        '"{name}" outcrossing outcrossed',
        '"{name}" "mixed mating"',
        '"{name}" "breeding system"',
    ),
}

COLOR = {
    "white": ("white", "whitish", "cream", "ivory"),
    "yellow_orange": ("yellow", "orange", "golden"),
    "red_pink": ("red", "pink", "rose", "crimson", "scarlet"),
    "blue_purple": ("blue", "purple", "violet", "lavender"),
    "green_brown_inconspicuous": ("greenish", "brownish", "brown"),
}
FLOWER_CONTEXT = ("flower", "flowers", "floral", "corolla", "petal", "petals", "bloom", "blooms")
GUILD = {
    "bumblebees": ("bumblebee", "bumblebees", "bombus"),
    "bees": ("pollinated by bees", "bee pollination", "bee-pollinated", "bee pollinators", "visited by bees"),
    "flies": ("pollinated by flies", "fly pollination", "fly-pollinated", "visited by flies"),
    "butterflies": ("pollinated by butterflies", "butterfly pollination"),
    "moths": ("pollinated by moths", "moth pollination", "hawkmoth"),
    "birds": ("pollinated by birds", "bird pollination", "bird-pollinated", "hummingbird"),
    "wind": ("wind-pollinated", "wind pollinated", "pollinated by wind", "anemophil"),
    "mixed": ("insect-pollinated", "insect pollinated", "entomophil", "generalist pollinat", "attracts pollinators"),
}
SI = {
    "SI": (r"\bself[- ]incompatib", r"\bself[- ]sterile\b"),
    "SC": (r"\bself[- ]compatib", r"\bself[- ]fertile\b"),
}
MATING = {
    "obligate_outcrossing": (r"\bobligate(?:ly)? outcross", r"\bstrict(?:ly)? outcross", r"\bdioecious\b"),
    "mixed_mating": (r"\bmixed[- ]mating\b", r"\bfacultative selfing\b", r"\bselfing and outcrossing\b", r"\bpredominantly outcrossing\b", r"\bmainly outcrossing\b"),
    "mainly_selfing": (r"\bpredominantly selfing\b", r"\bmainly selfing\b", r"\bhigh selfing rate\b", r"\bautogamous\b", r"\bself[- ]pollinating\b", r"\bcleistogamous\b"),
    "obligate_selfing": (r"\bobligate(?:ly)? selfing\b", r"\bstrict(?:ly)? selfing\b"),
}

OUTPUT_COLUMNS = [
    "accepted_species", "trait_group", "value", "inference_status", "confidence",
    "query_taxon", "search_scope", "search_engine", "search_query", "search_rank",
    "result_title", "result_snippet", "result_url",
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
    return " ".join(str(value or "").split())


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
    rows = []
    for item in root.findall(".//item"):
        rows.append({
            "title": clean(item.findtext("title")),
            "url": clean(item.findtext("link")),
            "snippet": clean(item.findtext("description")),
        })
    return [row for row in rows if row["url"]]


def genus(species: str) -> str:
    parts = species.split()
    return parts[0] if parts else ""


def taxon_present(text: str, taxon: str, scope: str) -> bool:
    low = text.casefold()
    return (taxon.casefold() in low) if scope == "species" else (genus(taxon).casefold() in low)


def extract_value(group: str, text: str) -> str:
    low = clean(text).casefold()
    if group == "flower_color":
        if not any(term in low for term in FLOWER_CONTEXT):
            return ""
        hits = [value for value, terms in COLOR.items() if any(re.search(rf"(?<!\w){re.escape(term)}(?!\w)", low) for term in terms)]
    elif group == "pollination":
        hits = [value for value, terms in GUILD.items() if any(term in low for term in terms)]
    elif group == "self_incompatibility":
        hits = [value for value, patterns in SI.items() if any(re.search(pattern, low) for pattern in patterns)]
    else:
        hits = [value for value, patterns in MATING.items() if any(re.search(pattern, low) for pattern in patterns)]
    hits = list(dict.fromkeys(hits))
    if len(hits) == 1:
        return hits[0]
    if group == "flower_color" and len(hits) > 1:
        return "multicolored_variable"
    if group == "pollination" and len(hits) > 1:
        return "mixed"
    return ""


def search(client: httpx.Client, query: str, limit: int) -> list[tuple[str, dict[str, str]]]:
    rows: list[tuple[str, dict[str, str]]] = []
    encoded = quote_plus(query)
    try:
        response = client.get(f"https://www.bing.com/search?format=rss&q={encoded}", timeout=25)
        response.raise_for_status()
        rows.extend(("bing_rss", row) for row in parse_bing(response.text)[:limit])
    except Exception:
        pass
    try:
        response = client.get(f"https://html.duckduckgo.com/html/?q={encoded}", timeout=25)
        response.raise_for_status()
        rows.extend(("duckduckgo_html", row) for row in parse_duck(response.text)[:limit])
    except Exception:
        pass
    return rows


def one_species(species: str, limit: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    plan = [(species, "species")]
    genus_name = genus(species)
    if genus_name:
        plan.append((genus_name, "genus"))
    headers = {"User-Agent": USER_AGENT, "Accept-Language": "en;q=0.9,*;q=0.5"}
    seen: set[tuple[str, str, str, str]] = set()
    with httpx.Client(headers=headers, follow_redirects=True) as client:
        for query_taxon, scope in plan:
            for group, templates in QUERY_GROUPS.items():
                for template in templates:
                    query = template.format(name=query_taxon)
                    for engine, result in search(client, query, limit):
                        title = clean(result.get("title"))
                        snippet = clean(result.get("snippet"))
                        url = clean(result.get("url"))
                        text = f"{title}. {snippet}"
                        if not url or not taxon_present(text, query_taxon, scope):
                            continue
                        value = extract_value(group, text)
                        if not value:
                            continue
                        if group == "self_incompatibility" and scope == "genus":
                            value = "likely_SI" if value == "SI" else "likely_SC"
                        key = (group, value, scope, url)
                        if key in seen:
                            continue
                        seen.add(key)
                        rows.append({
                            "accepted_species": species,
                            "trait_group": group,
                            "value": value,
                            "inference_status": "reported" if scope == "species" else "likely",
                            "confidence": "high" if scope == "species" and group in {"self_incompatibility", "mating_system"} else ("medium" if scope == "species" else "low"),
                            "query_taxon": query_taxon,
                            "search_scope": scope,
                            "search_engine": engine,
                            "search_query": query,
                            "search_rank": "",
                            "result_title": title,
                            "result_snippet": snippet,
                            "result_url": url,
                        })
    return rows


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(25, min=1, max=500),
    results_per_query: int = typer.Option(6, min=1, max=20),
    workers: int = typer.Option(16, min=1, max=32),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    species = [value for value in table["accepted_species"].astype(str).str.strip() if value][:max_taxa]
    rows: list[dict[str, str]] = []
    errors = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one_species, name, results_per_query): name for name in species}
        for future in as_completed(futures):
            try:
                rows.extend(future.result())
            except Exception:
                errors += 1
    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(["accepted_species", "trait_group", "value", "inference_status", "result_url"])
    coverage = {
        group: int(frame.loc[frame["trait_group"].eq(group), "accepted_species"].nunique()) if len(frame) else 0
        for group in QUERY_GROUPS
    }
    report = {
        "n_species": len(species),
        "coverage": coverage,
        "n_rows": len(frame),
        "n_errors": errors,
        "policy": "species explicit snippet = reported; genus-only snippet = low-confidence likely",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "reproductive_web_recovery.csv", index=False)
    (output_dir / "reproductive_web_recovery_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
