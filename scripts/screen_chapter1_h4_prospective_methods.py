"""Outcome-blind Methods screening for prospective H4 temporal replication.

Only Europe PMC bibliographic lite records and XML sections whose headings match the
frozen Methods-section contract are inspected. Results/Discussion text, abstracts,
reproductive-output values, and pollen-limitation effects are never written.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_h4_prospective_methods_preflight_v1"

BINOMIAL_RE = re.compile(r"\b([A-Z][a-z]{2,})\s+([a-z][a-z-]{2,})\b")
COORD_RE = re.compile(
    r"(?P<lat>\d{1,2}(?:\.\d+)?)\s*(?:°|º)?\s*(?P<ns>[NS])"
    r"[^\d]{0,20}"
    r"(?P<lon>\d{1,3}(?:\.\d+)?)\s*(?:°|º)?\s*(?P<ew>[EW])",
    flags=re.IGNORECASE,
)
COORD_RE_REVERSED = re.compile(
    r"(?P<lon>\d{1,3}(?:\.\d+)?)\s*(?:°|º)?\s*(?P<ew>[EW])"
    r"[^\d]{0,20}"
    r"(?P<lat>\d{1,2}(?:\.\d+)?)\s*(?:°|º)?\s*(?P<ns>[NS])",
    flags=re.IGNORECASE,
)


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected H4 Methods preflight contract")
    return config


def _get(url: str) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "island-h4-methods-preflight/1.0"},
    )
    for attempt in range(7):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 6:
                raise
            retry_after = exc.headers.get("Retry-After")
            delay = (
                float(retry_after)
                if retry_after and retry_after.isdigit()
                else min(30.0, 2.0 ** attempt)
            )
            time.sleep(delay)
    raise RuntimeError("unreachable Methods retry loop")


def _get_json(url: str) -> dict[str, Any]:
    return json.loads(_get(url).decode("utf-8"))


def _normalise_doi(value: object) -> str:
    text = str(value or "").strip().casefold()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if text.startswith(prefix):
            return text[len(prefix) :]
    return text


def discover_europe_pmc(config: dict[str, Any]) -> list[dict[str, str]]:
    source = config["source"]["Europe_PMC"]
    terms = [f'"{term}"' for term in source["query_terms"]]
    query = (
        "(" + " OR ".join(terms) + ")"
        + f" AND FIRST_PDATE:[{source['date_start']} TO {source['date_end']}]"
        + " AND OPEN_ACCESS:Y"
    )
    cursor = "*"
    rows: list[dict[str, str]] = []
    while True:
        params = urllib.parse.urlencode(
            {
                "query": query,
                "format": "json",
                "resultType": "lite",
                "pageSize": 1000,
                "cursorMark": cursor,
            }
        )
        payload = _get_json(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + params
        )
        result_list = payload.get("resultList", {}).get("result", [])
        if not result_list:
            break
        for item in result_list:
            pmcid = str(item.get("pmcid", "") or "").strip()
            if not pmcid:
                continue
            rows.append(
                {
                    "doi": _normalise_doi(item.get("doi")),
                    "pmcid": pmcid,
                    "publication_date": str(
                        item.get("firstPublicationDate")
                        or item.get("journalInfo", {}).get("printPublicationDate")
                        or ""
                    ),
                    "title": str(item.get("title", "") or ""),
                }
            )
        next_cursor = str(payload.get("nextCursorMark", "") or "")
        if not next_cursor or next_cursor == cursor:
            break
        cursor = next_cursor
        time.sleep(0.15)
    dedup: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row["pmcid"]
        dedup[key] = row
    return sorted(dedup.values(), key=lambda row: (row["publication_date"], row["pmcid"]))


def _section_title(sec: ET.Element) -> str:
    title_node = sec.find("title")
    return " ".join(title_node.itertext()).strip() if title_node is not None else ""


def _collect_allowed_section_text(
    node: ET.Element,
    *,
    exclude: list[re.Pattern[str]],
) -> list[str]:
    if node.tag == "sec":
        title = _section_title(node)
        if title and any(pattern.search(title) for pattern in exclude):
            return []

    pieces: list[str] = []
    if node.text and node.text.strip():
        pieces.append(node.text.strip())
    for child in node:
        if child.tag == "sec":
            pieces.extend(_collect_allowed_section_text(child, exclude=exclude))
        else:
            pieces.extend(
                part.strip()
                for part in child.itertext()
                if part and part.strip()
            )
        if child.tail and child.tail.strip():
            pieces.append(child.tail.strip())
    return pieces


def _section_texts(root: ET.Element, config: dict[str, Any]) -> list[str]:
    include = [
        re.compile(pattern, re.IGNORECASE)
        for pattern in config["methods_section_titles"]["include_regex"]
    ]
    exclude = [
        re.compile(pattern, re.IGNORECASE)
        for pattern in config["methods_section_titles"]["exclude_regex"]
    ]

    selected: list[str] = []

    def walk(node: ET.Element) -> None:
        for sec in node.findall("sec"):
            title = _section_title(sec)
            if title and any(pattern.search(title) for pattern in exclude):
                continue
            if title and any(pattern.search(title) for pattern in include):
                text = " ".join(
                    _collect_allowed_section_text(sec, exclude=exclude)
                ).strip()
                if text:
                    selected.append(text)
                continue
            walk(sec)

    body = root.find(".//body")
    walk(body if body is not None else root)
    return selected


def _species_lookup(trait_states: pd.DataFrame) -> dict[str, str]:
    if "accepted_species" not in trait_states.columns:
        raise typer.BadParameter("trait-state table missing accepted_species")
    lookup: dict[str, str] = {}
    for value in trait_states["accepted_species"].fillna("").astype(str):
        tokens = value.strip().split()
        if len(tokens) != 2:
            continue
        key = f"{tokens[0]} {tokens[1]}".casefold()
        if key in lookup and lookup[key] != value:
            continue
        lookup[key] = value
    return lookup


def _matched_species(text: str, lookup: dict[str, str]) -> list[str]:
    found: set[str] = set()
    for match in BINOMIAL_RE.finditer(text):
        key = f"{match.group(1)} {match.group(2)}".casefold()
        accepted = lookup.get(key)
        if accepted:
            found.add(accepted)
    return sorted(found)


def _coordinates(text: str) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    seen: set[tuple[float, float]] = set()
    for pattern in (COORD_RE, COORD_RE_REVERSED):
        for match in pattern.finditer(text):
            lat = float(match.group("lat"))
            lon = float(match.group("lon"))
            if match.group("ns").upper() == "S":
                lat *= -1
            if match.group("ew").upper() == "W":
                lon *= -1
            if not (-90 <= lat <= 90 and -180 <= lon <= 180):
                continue
            key = (round(lat, 6), round(lon, 6))
            if key in seen:
                continue
            seen.add(key)
            rows.append({"lat": key[0], "lon": key[1]})
    return rows


def _contains_any(text: str, values: list[str]) -> bool:
    lowered = text.casefold()
    return any(str(value).casefold() in lowered for value in values)


def screen_record(
    record: dict[str, str],
    *,
    lookup: dict[str, str],
    config: dict[str, Any],
) -> dict[str, Any]:
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/"
        + urllib.parse.quote(record["pmcid"])
        + "/fullTextXML"
    )
    try:
        raw = _get(url)
        root = ET.fromstring(raw)
    except Exception as exc:  # noqa: BLE001
        return {
            **record,
            "methods_section_sha256": "",
            "matched_species": "",
            "coordinate_pairs_json": "[]",
            "has_supplementation_term": False,
            "has_natural_control_term": False,
            "has_field_term": False,
            "exclusion_flags": "",
            "candidate_status": "xml_unavailable",
            "screen_error": type(exc).__name__,
        }

    sections = _section_texts(root, config)
    if not sections:
        return {
            **record,
            "methods_section_sha256": "",
            "matched_species": "",
            "coordinate_pairs_json": "[]",
            "has_supplementation_term": False,
            "has_natural_control_term": False,
            "has_field_term": False,
            "exclusion_flags": "",
            "candidate_status": "methods_section_unresolved",
            "screen_error": "",
        }

    text = "\n".join(sections)
    methods_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    species = _matched_species(text, lookup)
    coordinates = _coordinates(text)
    flags = config["design_flags"]
    has_supplement = _contains_any(text, flags["supplementation_terms"])
    has_control = _contains_any(text, flags["natural_control_terms"])
    has_field = _contains_any(text, flags["field_terms"])
    exclusion = [
        term
        for term in flags["exclusion_flags"]
        if str(term).casefold() in text.casefold()
    ]
    candidate = bool(
        species
        and has_supplement
        and has_control
        and has_field
        and coordinates
        and not exclusion
    )
    status = "automatic_methods_candidate" if candidate else "unresolved_or_excluded_by_conservative_screen"
    return {
        **record,
        "methods_section_sha256": methods_sha,
        "matched_species": "|".join(species),
        "coordinate_pairs_json": json.dumps(coordinates, separators=(",", ":")),
        "has_supplementation_term": has_supplement,
        "has_natural_control_term": has_control,
        "has_field_term": has_field,
        "exclusion_flags": "|".join(sorted(set(exclusion))),
        "candidate_status": status,
        "screen_error": "",
    }


@app.command("run")
def run(
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_csv: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    traits = pd.read_csv(trait_states_csv)
    lookup = _species_lookup(traits)
    records = discover_europe_pmc(config)
    rows: list[dict[str, Any]] = []
    for index, record in enumerate(records, start=1):
        rows.append(screen_record(record, lookup=lookup, config=config))
        if index % 25 == 0:
            typer.echo(f"screened {index}/{len(records)}")
        time.sleep(0.20)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "doi",
        "pmcid",
        "publication_date",
        "title",
        "methods_section_sha256",
        "matched_species",
        "coordinate_pairs_json",
        "has_supplementation_term",
        "has_natural_control_term",
        "has_field_term",
        "exclusion_flags",
        "candidate_status",
        "screen_error",
    ]
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "n_europe_pmc_oa_records": len(records),
        "n_automatic_methods_candidates": sum(
            row["candidate_status"] == "automatic_methods_candidate" for row in rows
        ),
        "n_methods_section_unresolved": sum(
            row["candidate_status"] == "methods_section_unresolved" for row in rows
        ),
        "n_xml_unavailable": sum(row["candidate_status"] == "xml_unavailable" for row in rows),
        "methods_text_materialized": False,
        "results_text_materialized": False,
        "effect_sizes_materialized": False,
    }
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
