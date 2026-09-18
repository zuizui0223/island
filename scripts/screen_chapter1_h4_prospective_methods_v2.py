"""Outcome-blind Methods v2 screen for the prospective H4 validation cohort.

V2 fixes two design issues found in the v1 outcome-blind pilot:
1. only records in the already frozen OpenAlex/Crossref sampling frame are eligible;
2. target species must come from the article title or a target-species Methods section,
   never from arbitrary binomials elsewhere in Methods.

No abstract, Results/Discussion text, reproductive output, or pollen-limitation effect is
materialized or used.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen_chapter1_h4_prospective_methods as v1  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_h4_prospective_methods_preflight_v2"


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected H4 Methods v2 contract")
    return config


def _compile(values: list[str]) -> list[re.Pattern[str]]:
    return [re.compile(value, re.IGNORECASE) for value in values]


def _section_title(sec: ET.Element) -> str:
    title = sec.find("title")
    return " ".join(title.itertext()).strip() if title is not None else ""


def _safe_section_text(
    sec: ET.Element,
    *,
    hard_exclude: list[re.Pattern[str]],
) -> str:
    """Return section text while dropping nested excluded sections entirely."""

    pieces: list[str] = []

    def walk(node: ET.Element) -> None:
        if node.tag == "sec":
            title = _section_title(node)
            if title and any(pattern.search(title) for pattern in hard_exclude):
                return
        if node.text and node.text.strip():
            pieces.append(node.text.strip())
        for child in node:
            if child.tag == "sec":
                walk(child)
            else:
                pieces.extend(
                    part.strip()
                    for part in child.itertext()
                    if part and part.strip()
                )
            if child.tail and child.tail.strip():
                pieces.append(child.tail.strip())

    walk(sec)
    return " ".join(pieces).strip()


def _selected_sections(
    root: ET.Element,
    config: dict[str, Any],
) -> list[tuple[str, str]]:
    policy = config["section_policy"]
    include = _compile([str(x) for x in policy["methods_include_regex"]])
    hard_exclude = _compile([str(x) for x in policy["hard_exclude_regex"]])
    selected: list[tuple[str, str]] = []

    def walk(parent: ET.Element) -> None:
        for sec in parent.findall("sec"):
            title = _section_title(sec)
            if title and any(pattern.search(title) for pattern in hard_exclude):
                continue
            if title and any(pattern.search(title) for pattern in include):
                text = _safe_section_text(sec, hard_exclude=hard_exclude)
                if text:
                    selected.append((title, text))
                # Keep scanning nested allowed sections so a specific
                # "Study species" / "Study system" heading remains visible
                # for target-species resolution. Hard-excluded Results/
                # Discussion sections are still skipped at their own boundary.
                walk(sec)
                continue
            walk(sec)

    body = root.find(".//body")
    walk(body if body is not None else root)
    return selected


def _target_species(
    *,
    article_title: str,
    sections: list[tuple[str, str]],
    lookup: dict[str, str],
    config: dict[str, Any],
) -> tuple[list[str], str, list[str]]:
    title_matches = v1._matched_species(article_title, lookup)
    target_patterns = _compile(
        [str(x) for x in config["section_policy"]["target_species_section_regex"]]
    )
    target_text = "\n".join(
        text
        for title, text in sections
        if any(pattern.search(title) for pattern in target_patterns)
    )
    section_matches = v1._matched_species(target_text, lookup) if target_text else []
    broad_matches = v1._matched_species(
        "\n".join(text for _, text in sections),
        lookup,
    )

    if title_matches:
        return title_matches, "article_title", broad_matches
    if section_matches:
        return section_matches, "target_species_section", broad_matches
    return [], "unresolved", broad_matches


def _normalise_frame(metadata: pd.DataFrame) -> pd.DataFrame:
    required = {
        "source_database",
        "source_record_id",
        "doi",
        "title",
        "publication_date",
        "publication_year",
        "query_family",
    }
    if missing := required - set(metadata.columns):
        raise typer.BadParameter(
            f"frozen metadata frame missing columns: {sorted(missing)}"
        )
    work = metadata.copy().fillna("")
    work["doi"] = work["doi"].map(v1._normalise_doi)
    work["title"] = work["title"].astype(str)
    work["publication_date"] = work["publication_date"].astype(str)
    if work["doi"].duplicated().any():
        duplicates = work.loc[work["doi"].ne("") & work["doi"].duplicated(), "doi"]
        if not duplicates.empty:
            raise typer.BadParameter(
                f"frozen metadata DOI must be unique after deduplication: {duplicates.iloc[0]}"
            )
    return work.reset_index(drop=True)


def _europe_pmc_map(
    metadata: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, dict[str, str]]:
    """Resolve OA PMC mirrors only for DOI records already in the frozen frame."""

    source = config["source"]["Europe_PMC"]
    dois = sorted({value for value in metadata["doi"].astype(str) if value})
    out: dict[str, dict[str, str]] = {}
    batch_size = 30
    for start in range(0, len(dois), batch_size):
        batch = dois[start : start + batch_size]
        doi_query = " OR ".join(f'DOI:"{doi}"' for doi in batch)
        query = (
            "(" + doi_query + ")"
            + f" AND FIRST_PDATE:[{source['date_start']} TO {source['date_end']}]"
            + " AND OPEN_ACCESS:Y"
        )
        params = urllib.parse.urlencode(
            {
                "query": query,
                "format": "json",
                "resultType": "lite",
                "pageSize": 1000,
            }
        )
        payload = v1._get_json(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + params
        )
        for item in payload.get("resultList", {}).get("result", []):
            doi = v1._normalise_doi(item.get("doi"))
            pmcid = str(item.get("pmcid", "") or "").strip()
            if doi in batch and pmcid and doi not in out:
                out[doi] = {
                    "doi": doi,
                    "pmcid": pmcid,
                    "publication_date": str(
                        item.get("firstPublicationDate")
                        or item.get("journalInfo", {}).get("printPublicationDate")
                        or ""
                    ),
                    "title": str(item.get("title", "") or ""),
                }
        time.sleep(0.10)
    return out


def _screen_one(
    frozen: dict[str, str],
    pmc: dict[str, str],
    *,
    lookup: dict[str, str],
    config: dict[str, Any],
) -> dict[str, Any]:
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/"
        + urllib.parse.quote(str(pmc["pmcid"]))
        + "/fullTextXML"
    )
    base = {
        "doi": frozen["doi"],
        "pmcid": str(pmc["pmcid"]),
        "publication_date": frozen["publication_date"],
        "title": frozen["title"],
        "source_database": frozen["source_database"],
        "source_record_id": frozen["source_record_id"],
    }
    try:
        raw = v1._get(url)
        root = ET.fromstring(raw)
    except Exception as exc:  # noqa: BLE001
        return {
            **base,
            "methods_section_sha256": "",
            "target_species": "",
            "target_species_source": "unresolved",
            "broad_methods_species_diagnostic": "",
            "coordinate_pairs_json": "[]",
            "has_supplementation_term": False,
            "has_natural_control_term": False,
            "has_field_term": False,
            "exclusion_flags": "",
            "field_status": "unresolved",
            "exclusion_status": "unresolved",
            "candidate_status": "xml_unavailable",
            "screen_error": type(exc).__name__,
        }

    sections = _selected_sections(root, config)
    if not sections:
        return {
            **base,
            "methods_section_sha256": "",
            "target_species": "",
            "target_species_source": "unresolved",
            "broad_methods_species_diagnostic": "",
            "coordinate_pairs_json": "[]",
            "has_supplementation_term": False,
            "has_natural_control_term": False,
            "has_field_term": False,
            "exclusion_flags": "",
            "field_status": "unresolved",
            "exclusion_status": "unresolved",
            "candidate_status": "methods_section_unresolved",
            "screen_error": "",
        }

    methods_text = "\n".join(text for _, text in sections)
    methods_sha = hashlib.sha256(methods_text.encode("utf-8")).hexdigest()
    target_species, target_source, broad_species = _target_species(
        article_title=frozen["title"],
        sections=sections,
        lookup=lookup,
        config=config,
    )
    coordinates = v1._coordinates(methods_text)
    flags = config["design_flags"]
    has_supplement = v1._contains_any(
        methods_text, [str(x) for x in flags["supplementation_terms"]]
    )
    has_control = v1._contains_any(
        methods_text, [str(x) for x in flags["natural_control_terms"]]
    )
    has_field = v1._contains_any(
        methods_text, [str(x) for x in flags["field_terms"]]
    )
    exclusion = [
        str(term)
        for term in flags["exclusion_terms"]
        if str(term).casefold() in methods_text.casefold()
    ]

    core_design = bool(target_species and has_supplement and has_control)
    if not core_design:
        candidate_status = "not_design_candidate"
    elif exclusion and not has_field:
        candidate_status = "excluded_by_conservative_environment_flag"
    elif exclusion:
        candidate_status = "design_candidate_conflict_unresolved"
    elif not has_field:
        candidate_status = "design_candidate_field_unresolved"
    else:
        candidate_status = "design_candidate_auto"

    return {
        **base,
        "methods_section_sha256": methods_sha,
        "target_species": "|".join(target_species),
        "target_species_source": target_source,
        "broad_methods_species_diagnostic": "|".join(broad_species),
        "coordinate_pairs_json": json.dumps(coordinates, separators=(",", ":")),
        "has_supplementation_term": has_supplement,
        "has_natural_control_term": has_control,
        "has_field_term": has_field,
        "exclusion_flags": "|".join(sorted(set(exclusion))),
        "field_status": "supported" if has_field else "unresolved",
        "exclusion_status": (
            "none"
            if not exclusion
            else "conflict_with_field"
            if has_field
            else "conservative_exclusion"
        ),
        "candidate_status": candidate_status,
        "screen_error": "",
    }


@app.command("run")
def run(
    metadata_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_csv: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    metadata = _normalise_frame(pd.read_csv(metadata_csv))
    traits = pd.read_csv(trait_states_csv)
    lookup = v1._species_lookup(traits)

    pmc_map = _europe_pmc_map(metadata, config)
    eligible = metadata.loc[
        metadata["doi"].ne("") & metadata["doi"].isin(pmc_map)
    ].copy()

    rows: list[dict[str, Any]] = []
    for index, frozen in enumerate(eligible.to_dict("records"), start=1):
        rows.append(
            _screen_one(
                frozen,
                pmc_map[frozen["doi"]],
                lookup=lookup,
                config=config,
            )
        )
        if index % 25 == 0:
            typer.echo(f"screened {index}/{len(eligible)} frozen-frame OA records")
        time.sleep(0.15)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "doi",
        "pmcid",
        "publication_date",
        "title",
        "source_database",
        "source_record_id",
        "methods_section_sha256",
        "target_species",
        "target_species_source",
        "broad_methods_species_diagnostic",
        "coordinate_pairs_json",
        "has_supplementation_term",
        "has_natural_control_term",
        "has_field_term",
        "exclusion_flags",
        "field_status",
        "exclusion_status",
        "candidate_status",
        "screen_error",
    ]
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    frame = pd.DataFrame(rows)
    counts = (
        frame["candidate_status"].value_counts(dropna=False).to_dict()
        if not frame.empty
        else {}
    )
    summary = {
        "n_frozen_metadata_records": int(len(metadata)),
        "n_frozen_records_with_oa_pmc_match": int(len(eligible)),
        "candidate_status_counts": {str(k): int(v) for k, v in counts.items()},
        "n_target_species_records": int(
            frame["target_species"].fillna("").ne("").sum()
        )
        if not frame.empty
        else 0,
        "n_records_with_explicit_coordinates": int(
            frame["coordinate_pairs_json"].fillna("[]").ne("[]").sum()
        )
        if not frame.empty
        else 0,
        "methods_text_materialized": False,
        "abstracts_materialized": False,
        "results_text_materialized": False,
        "effect_sizes_materialized": False,
    }
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
