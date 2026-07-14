#!/usr/bin/env python3
"""Harvest species treatments and floral traits from any eFloras flora_id.

The traversal follows the eFloras taxonomic hierarchy rather than generic page
crawling. Every extracted value retains the source URL and original excerpt.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import random
import re
import time
from collections import deque
from pathlib import Path
from urllib.parse import parse_qs, urljoin, urlparse

import requests
from bs4 import BeautifulSoup

BASE = "http://www.efloras.org/"
BINOMIAL = re.compile(r"\b([A-Z][a-zA-Z-]+\s+[a-z][a-zA-Z.-]+)\b")
FLORAL = re.compile(r"\b(flowers?|corolla|petals?|perianth|tepals?|calyx|sepals?)\b", re.I)
SPLIT = re.compile(r"(?<=[.;!?])\s+|\n+")
COLORS = ("white", "whitish", "cream", "yellow", "yellowish", "orange", "red", "reddish", "pink", "purple", "violet", "blue", "green", "brown", "black")
FORMS = {
    "campanulate": "campanulate", "bell-shaped": "campanulate",
    "tubular": "tubular", "funnelform": "funnelform",
    "funnel-shaped": "funnelform", "rotate": "rotate",
    "salverform": "salverform", "urceolate": "urceolate",
    "bilabiate": "bilabiate", "two-lipped": "bilabiate",
    "2-lipped": "bilabiate", "papilionaceous": "papilionaceous",
    "spurred": "spurred",
}
SYMMETRY = {
    "actinomorphic": "radial", "radially symmetric": "radial",
    "radially symmetrical": "radial", "zygomorphic": "bilateral",
    "bilaterally symmetric": "bilateral", "bilaterally symmetrical": "bilateral",
}
SEX = ("bisexual", "hermaphroditic", "unisexual", "monoecious", "dioecious", "andromonoecious", "gynomonoecious")


def clean(value: str) -> str:
    return " ".join(html.unescape(value or "").split())


def query_id(url: str, keys=("taxon_id", "start_taxon_id")) -> str | None:
    query = parse_qs(urlparse(url).query)
    for key in keys:
        values = query.get(key)
        if values and values[0].isdigit():
            return values[0]
    return None


class Client:
    def __init__(self, delay_min: float, delay_max: float):
        self.delay_min = delay_min
        self.delay_max = delay_max
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "island-v2-efloras-harvester/1.1 (+https://github.com/zuizui0223/island)"

    def get(self, url: str, attempts: int = 6) -> requests.Response:
        last = None
        for attempt in range(attempts):
            try:
                response = self.session.get(url, timeout=120)
                if response.status_code == 200:
                    time.sleep(random.uniform(self.delay_min, self.delay_max))
                    return response
                last = RuntimeError(f"HTTP {response.status_code}")
            except Exception as exc:
                last = exc
            if attempt + 1 < attempts:
                time.sleep(min(60, 2 ** attempt))
        raise RuntimeError(f"GET failed: {url}: {last}")


def load_master(path: Path | None, column: str) -> set[str] | None:
    if path is None:
        return None
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or column not in reader.fieldnames:
            raise ValueError(f"Missing master column {column!r}")
        return {clean(row[column]) for row in reader if row.get(column)}


def listing_links(client: Client, url: str, flora_id: int, max_pages: int) -> list[tuple[str, str]]:
    found: dict[str, str] = {}
    for page in range(1, max_pages + 1):
        sep = "&" if "?" in url else "?"
        response = client.get(f"{url}{sep}page={page}")
        soup = BeautifulSoup(response.text, "html.parser")
        before = len(found)
        for link in soup.find_all("a", href=True):
            absolute = urljoin(response.url, link["href"])
            q = parse_qs(urlparse(absolute).query)
            linked_flora = q.get("flora_id", [str(flora_id)])[0]
            if linked_flora != str(flora_id):
                continue
            taxon_id = query_id(absolute)
            label = clean(link.get_text(" ", strip=True))
            if taxon_id and label and label.lower() not in {"next", "previous", "first", "last", "back"}:
                found.setdefault(taxon_id, label)
        if len(found) == before:
            break
    return list(found.items())


def treatment(client: Client, flora_id: int, taxon_id: str) -> tuple[str, str, str]:
    response = client.get(f"{BASE}florataxon.aspx?flora_id={flora_id}&taxon_id={taxon_id}")
    soup = BeautifulSoup(response.text, "html.parser")
    node = soup.find(id="lblTaxonDesc")
    description = clean(node.get_text(" ", strip=True)) if node else ""
    return response.text, response.url, description


def candidates(flora_id: int, flora_name: str, species: str, taxon_id: str, text: str, url: str) -> list[dict[str, str]]:
    rows = []
    seen = set()
    for sentence in SPLIT.split(text):
        sentence = clean(sentence)
        if not 12 <= len(sentence) <= 1500:
            continue
        low = sentence.lower()
        hits = []
        if FLORAL.search(sentence):
            color_hits = sorted({word[:-3] if word.endswith("ish") else word for word in COLORS if re.search(rf"\b{re.escape(word)}\b", low)})
            if color_hits:
                hits.append(("flower_primary_color", "|".join(color_hits), "flower-colour description"))
            form_hits = sorted({value for phrase, value in FORMS.items() if phrase in low})
            if form_hits:
                hits.append(("floral_form", "|".join(form_hits), "corolla/form description"))
            symmetry_hits = sorted({value for phrase, value in SYMMETRY.items() if phrase in low})
            if symmetry_hits:
                hits.append(("floral_symmetry", "|".join(symmetry_hits), "symmetry description"))
        sex_hits = sorted({word for word in SEX if re.search(rf"\b{word}\b", low)})
        if sex_hits:
            hits.append(("sex_system", "|".join(sex_hits), "sexual-system description"))
        if re.search(r"\bself[- ]incompat", low):
            hits.append(("self_incompatibility", "self_incompatible", "self-incompatibility description"))
        elif re.search(r"\bself[- ]compat", low):
            hits.append(("self_incompatibility", "self_compatible", "self-compatibility description"))
        for trait, value, rule in hits:
            key = (trait, value, sentence)
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "accepted_species": species,
                "trait_name": trait,
                "candidate_value": value,
                "source_name": flora_name,
                "flora_id": str(flora_id),
                "source_trait": rule,
                "source_url": url,
                "source_record_id": f"efloras:{flora_id}:{taxon_id}",
                "source_excerpt": sentence,
                "evidence_scope": "species_direct",
                "evidence_status": "source_backed_rule_extracted_unreviewed",
            })
    return rows


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--flora-id", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("data/v2/staging/traits/bulk/efloras"))
    parser.add_argument("--master-csv", type=Path)
    parser.add_argument("--master-name-column", default="accepted_species")
    parser.add_argument("--max-pages", type=int, default=20)
    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--max-species", type=int, default=0)
    parser.add_argument("--delay-min", type=float, default=1.5)
    parser.add_argument("--delay-max", type=float, default=3.5)
    args = parser.parse_args()

    client = Client(args.delay_min, args.delay_max)
    master = load_master(args.master_csv, args.master_name_column)
    out = args.output_dir / f"flora-{args.flora_id}"
    html_dir = out / "html"
    html_dir.mkdir(parents=True, exist_ok=True)

    landing = client.get(f"{BASE}flora_page.aspx?flora_id={args.flora_id}")
    soup = BeautifulSoup(landing.text, "html.parser")
    flora_name = clean(soup.title.get_text(" ", strip=True) if soup.title else f"eFloras flora {args.flora_id}")

    queue = deque()
    for link in soup.find_all("a", href=True):
        absolute = urljoin(landing.url, link["href"])
        if "volume_page.aspx" in absolute or "browse.aspx" in absolute:
            if parse_qs(urlparse(absolute).query).get("flora_id", [str(args.flora_id)])[0] == str(args.flora_id):
                queue.append((absolute, 0))
    queue.append((f"{BASE}browse.aspx?flora_id={args.flora_id}", 0))

    visited = set()
    discovered: dict[str, str] = {}
    errors = []
    while queue:
        url, depth = queue.popleft()
        canonical = re.sub(r"&page=\d+$", "", url)
        if canonical in visited or depth > args.max_depth:
            continue
        visited.add(canonical)
        try:
            links = listing_links(client, canonical, args.flora_id, args.max_pages)
        except Exception as exc:
            errors.append({"stage": "browse", "url": canonical, "error": str(exc)})
            continue
        for taxon_id, label in links:
            match = BINOMIAL.search(re.sub(r"^\s*\d+[a-zA-Z]?\.\s*", "", label))
            if match:
                discovered.setdefault(taxon_id, match.group(1))
            elif depth < args.max_depth:
                queue.append((f"{BASE}browse.aspx?flora_id={args.flora_id}&start_taxon_id={taxon_id}", depth + 1))

    selected = [(taxon_id, species) for taxon_id, species in sorted(discovered.items()) if master is None or species in master]
    if args.max_species > 0:
        selected = selected[:args.max_species]

    treatment_rows = []
    candidate_rows = []
    for number, (taxon_id, species) in enumerate(selected, 1):
        try:
            raw, final, description = treatment(client, args.flora_id, taxon_id)
        except Exception as exc:
            errors.append({"stage": "species", "taxon_id": taxon_id, "species": species, "error": str(exc)})
            continue
        if len(description) < 20:
            continue
        (html_dir / f"{taxon_id}.html").write_text(raw, encoding="utf-8")
        extracted = candidates(args.flora_id, flora_name, species, taxon_id, description, final)
        candidate_rows.extend(extracted)
        treatment_rows.append({
            "accepted_species": species,
            "flora_id": args.flora_id,
            "flora_name": flora_name,
            "taxon_id": taxon_id,
            "source_url": final,
            "description_chars": len(description),
            "trait_candidate_rows": len(extracted),
        })
        if number % 25 == 0:
            print(f"{number}/{len(selected)} treatments; candidates={len(candidate_rows)}", flush=True)

    write_csv(out / "species_treatments.csv", treatment_rows, ["accepted_species", "flora_id", "flora_name", "taxon_id", "source_url", "description_chars", "trait_candidate_rows"])
    write_csv(out / "trait_candidates.csv", candidate_rows, ["accepted_species", "trait_name", "candidate_value", "source_name", "flora_id", "source_trait", "source_url", "source_record_id", "source_excerpt", "evidence_scope", "evidence_status"])

    report = {
        "flora_id": args.flora_id,
        "flora_name": flora_name,
        "browse_pages_visited": len(visited),
        "species_discovered": len(discovered),
        "master_species_selected": len(selected),
        "species_treatments_saved": len(treatment_rows),
        "species_with_trait_candidates": len({row["accepted_species"] for row in candidate_rows}),
        "trait_candidate_rows": len(candidate_rows),
        "trait_candidate_cells": len({(row["accepted_species"], row["trait_name"]) for row in candidate_rows}),
        "species_by_trait": {trait: len({row["accepted_species"] for row in candidate_rows if row["trait_name"] == trait}) for trait in sorted({row["trait_name"] for row in candidate_rows})},
        "errors": errors[:200],
    }
    (out / "coverage_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
