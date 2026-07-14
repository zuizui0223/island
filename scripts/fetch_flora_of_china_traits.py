#!/usr/bin/env python3
"""Harvest species descriptions and direct floral-trait candidates from Flora of China.

This follows the proven eFloras hierarchy used by Bamba et al.:
volume -> family -> genus -> species taxon IDs -> species treatment text.
It avoids generic site crawling and stores every source excerpt and URL.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import re
import time
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

BASE = "http://www.efloras.org/"
FLORA_ID = 2
MASTER = Path("data/v2/staging/gbif/collected/island_taxa.csv")
OUT_ROOT = Path("data/v2/staging/traits/bulk/efloras")
UA = {"User-Agent": "island-v2-flora-of-china-harvester/1.0 (+https://github.com/zuizui0223/island)"}

COLOR_WORDS = (
    "white", "whitish", "cream", "creamy", "yellow", "yellowish", "orange",
    "red", "reddish", "pink", "pinkish", "purple", "purplish", "violet",
    "blue", "bluish", "green", "greenish", "brown", "black",
)
FORM_MAP = {
    "campanulate": "campanulate", "bell-shaped": "campanulate",
    "tubular": "tubular", "funnelform": "funnelform", "funnel-shaped": "funnelform",
    "rotate": "rotate", "salverform": "salverform", "urceolate": "urceolate",
    "bilabiate": "bilabiate", "2-lipped": "bilabiate", "two-lipped": "bilabiate",
    "papilionaceous": "papilionaceous", "spurred": "spurred",
}
SYMMETRY_MAP = {
    "actinomorphic": "radial", "radially symmetric": "radial", "radially symmetrical": "radial",
    "zygomorphic": "bilateral", "bilaterally symmetric": "bilateral", "bilaterally symmetrical": "bilateral",
}
SEX_WORDS = ("bisexual", "hermaphroditic", "unisexual", "monoecious", "dioecious", "andromonoecious", "gynomonoecious")
FLORAL_CONTEXT = re.compile(r"\b(flowers?|corolla|petals?|perianth|tepals?|calyx)\b", re.I)
SENTENCE_SPLIT = re.compile(r"(?<=[.;!?])\s+|\n+")


class Parser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[tuple[str, str]] = []
        self._href = ""
        self._anchor: list[str] = []
        self._capture_desc = False
        self._desc_depth = 0
        self.desc: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        attrs = dict(attrs)
        if tag.lower() == "a":
            self._href = attrs.get("href", "")
            self._anchor = []
        if attrs.get("id") == "lblTaxonDesc":
            self._capture_desc = True
            self._desc_depth = 1
        elif self._capture_desc:
            self._desc_depth += 1
        if self._capture_desc and tag.lower() in {"p", "br", "div", "li", "tr"}:
            self.desc.append("\n")

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "a" and self._href:
            self.links.append((self._href, " ".join("".join(self._anchor).split())))
            self._href = ""
            self._anchor = []
        if self._capture_desc:
            self._desc_depth -= 1
            if self._desc_depth <= 0:
                self._capture_desc = False

    def handle_data(self, data: str) -> None:
        if self._href:
            self._anchor.append(data)
        if self._capture_desc:
            self.desc.append(data)


def get(url: str, attempts: int = 6) -> tuple[str, str]:
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=180) as response:
                return response.read().decode("utf-8", "replace"), response.geturl()
        except Exception as exc:
            last = exc
            if attempt + 1 < attempts:
                time.sleep(min(60, 2 ** attempt))
    raise RuntimeError(f"GET failed: {url}: {last}")


def master_names() -> set[str]:
    with MASTER.open(encoding="utf-8-sig", newline="") as handle:
        return {
            " ".join((row.get("accepted_species") or "").split())
            for row in csv.DictReader(handle)
            if row.get("accepted_species")
        }


def stable_bucket(value: str, count: int) -> int:
    return int.from_bytes(hashlib.sha1(value.encode()).digest()[:8], "big") % count


def parse_links(raw: str) -> list[tuple[str, str]]:
    parser = Parser()
    parser.feed(raw)
    return parser.links


def taxon_links(url: str, max_pages: int = 20) -> list[tuple[str, str]]:
    """Return unique taxon_id/name links from a volume or browse listing."""
    found: dict[str, str] = {}
    for page in range(1, max_pages + 1):
        sep = "&" if "?" in url else "?"
        raw, final = get(f"{url}{sep}page={page}")
        before = len(found)
        for href, label in parse_links(raw):
            absolute = urllib.parse.urljoin(final, href)
            match = re.search(r"(?:taxon_id|start_taxon_id)=(\d+)", absolute)
            name = " ".join(html.unescape(label).split())
            if match and name and name.lower() not in {"next", "previous", "first", "last"}:
                found.setdefault(match.group(1), name)
        if len(found) == before:
            break
    return [(taxon_id, name) for taxon_id, name in found.items()]


def species_description(taxon_id: str) -> tuple[str, str, str]:
    url = f"{BASE}florataxon.aspx?flora_id={FLORA_ID}&taxon_id={taxon_id}"
    raw, final = get(url)
    parser = Parser()
    parser.feed(raw)
    text = " ".join(html.unescape(" ".join(parser.desc)).split())
    return raw, final, text


def extract_candidates(species: str, text: str, url: str, record_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for sentence in SENTENCE_SPLIT.split(text):
        sentence = " ".join(sentence.split())
        if len(sentence) < 12 or len(sentence) > 1200:
            continue
        low = sentence.lower()
        hits: list[tuple[str, str, str]] = []
        if FLORAL_CONTEXT.search(sentence):
            colors = sorted({w.rstrip("ish") if w.endswith("ish") else w for w in COLOR_WORDS if re.search(rf"\b{re.escape(w)}\b", low)})
            if colors:
                hits.append(("flower_primary_color", "|".join(colors), "flower-colour description"))
            forms = sorted({value for phrase, value in FORM_MAP.items() if phrase in low})
            if forms:
                hits.append(("floral_form", "|".join(forms), "corolla/form description"))
            sym = sorted({value for phrase, value in SYMMETRY_MAP.items() if phrase in low})
            if sym:
                hits.append(("floral_symmetry", "|".join(sym), "symmetry description"))
        sex = sorted({word for word in SEX_WORDS if re.search(rf"\b{word}\b", low)})
        if sex:
            hits.append(("sex_system", "|".join(sex), "sexual-system description"))
        if "self-incompat" in low or "self incompatible" in low:
            hits.append(("self_incompatibility", "self_incompatible", "self-incompatibility description"))
        elif "self-compat" in low or "self compatible" in low:
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
                "source_name": "Flora of China (eFloras)",
                "source_trait": rule,
                "source_url": url,
                "source_record_id": record_id,
                "source_excerpt": sentence,
                "evidence_scope": "species_direct",
                "evidence_status": "source_backed_rule_extracted_unreviewed",
            })
    return rows


def main() -> None:
    shard_index = int(os.environ.get("EFLORAS_SHARD_INDEX", "0"))
    shard_count = max(1, int(os.environ.get("EFLORAS_SHARD_COUNT", "1")))
    max_species = int(os.environ.get("EFLORAS_MAX_SPECIES", "0"))
    out = OUT_ROOT / f"shard-{shard_index}"
    html_dir = out / "html"
    html_dir.mkdir(parents=True, exist_ok=True)
    master = master_names()

    top_raw, top_final = get(f"{BASE}flora_page.aspx?flora_id={FLORA_ID}")
    volumes: list[str] = []
    for href, _ in parse_links(top_raw):
        absolute = urllib.parse.urljoin(top_final, href)
        if "volume_page.aspx" in absolute and "flora_id=2" in absolute and absolute not in volumes:
            volumes.append(absolute)
    volumes = volumes[1:]  # volume 1 is introductory
    assigned_volumes = [v for i, v in enumerate(volumes) if i % shard_count == shard_index]

    candidates: list[dict[str, str]] = []
    treatments: list[dict[str, object]] = []
    discovered_species: dict[str, str] = {}
    errors: list[dict[str, str]] = []

    for volume_url in assigned_volumes:
        try:
            families = taxon_links(volume_url)
        except Exception as exc:
            errors.append({"stage": "volume", "url": volume_url, "error": str(exc)})
            continue
        for family_id, family_name in families:
            try:
                genera = taxon_links(f"{BASE}browse.aspx?flora_id={FLORA_ID}&start_taxon_id={family_id}")
            except Exception as exc:
                errors.append({"stage": "family", "url": family_id, "error": str(exc)})
                continue
            genera = [(i, n) for i, n in genera if len(n.split()) == 1]
            for genus_id, genus_name in genera:
                try:
                    species_links = taxon_links(f"{BASE}browse.aspx?flora_id={FLORA_ID}&start_taxon_id={genus_id}")
                except Exception as exc:
                    errors.append({"stage": "genus", "url": genus_id, "error": str(exc)})
                    continue
                for species_id, label in species_links:
                    cleaned = re.sub(r"^\d+\.\s*", "", label).strip()
                    match = re.search(r"\b([A-Z][a-zA-Z-]+\s+[a-z][a-zA-Z.-]+)\b", cleaned)
                    if match:
                        discovered_species.setdefault(species_id, match.group(1))

    selected = [
        (taxon_id, species)
        for taxon_id, species in sorted(discovered_species.items())
        if species in master and stable_bucket(species, shard_count) == shard_index
    ]
    if max_species > 0:
        selected = selected[:max_species]

    for number, (taxon_id, species) in enumerate(selected, 1):
        try:
            raw, final, description = species_description(taxon_id)
        except Exception as exc:
            errors.append({"stage": "species", "url": taxon_id, "error": str(exc)})
            continue
        if len(description) < 20:
            continue
        filename = f"{taxon_id}.html"
        (html_dir / filename).write_text(raw, encoding="utf-8")
        rows = extract_candidates(species, description, final, f"foc:{taxon_id}")
        candidates.extend(rows)
        treatments.append({
            "accepted_species": species,
            "taxon_id": taxon_id,
            "source_url": final,
            "description_chars": len(description),
            "trait_candidate_rows": len(rows),
        })
        if number % 25 == 0:
            time.sleep(1)

    fields = ["accepted_species", "trait_name", "candidate_value", "source_name", "source_trait", "source_url", "source_record_id", "source_excerpt", "evidence_scope", "evidence_status"]
    with (out / "efloras_species_direct_candidates.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(candidates)
    with (out / "species_treatments.csv").open("w", encoding="utf-8", newline="") as handle:
        fields2 = ["accepted_species", "taxon_id", "source_url", "description_chars", "trait_candidate_rows"]
        writer = csv.DictWriter(handle, fieldnames=fields2)
        writer.writeheader()
        writer.writerows(treatments)

    cells = {(row["accepted_species"], row["trait_name"]) for row in candidates}
    by_trait: dict[str, int] = {}
    for trait in sorted({row["trait_name"] for row in candidates}):
        by_trait[trait] = len({row["accepted_species"] for row in candidates if row["trait_name"] == trait})
    report = {
        "source": "Flora of China via eFloras",
        "method": "volume-family-genus-species taxon-ID traversal",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "volumes_assigned": len(assigned_volumes),
        "species_discovered": len(discovered_species),
        "master_species_selected": len(selected),
        "species_treatments_saved": len(treatments),
        "species_with_trait_candidates": len({row["accepted_species"] for row in candidates}),
        "trait_candidate_rows": len(candidates),
        "trait_candidate_cells": len(cells),
        "species_by_trait": by_trait,
        "errors": errors[:200],
    }
    (out / "coverage_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
