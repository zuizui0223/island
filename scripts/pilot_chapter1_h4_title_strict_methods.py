"""Fast outcome-blind lower-bound pilot for prospective H4 Methods support.

This pilot only follows metadata records whose titles explicitly contain a frozen
pollen-supplementation phrase. It is diagnostic only and cannot exclude records from
the final Methods screen.
"""

from __future__ import annotations

import csv
import re
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd
import typer

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen_chapter1_h4_prospective_methods as methods_screen  # noqa: E402

app = typer.Typer(add_completion=False, no_args_is_help=True)

STRICT_TITLE = re.compile(
    r"pollen limit|pollen supplement|supplemental pollin|hand[- ]pollin|pollen addition",
    re.IGNORECASE,
)


def _pmcid_for_doi(doi: str) -> str:
    if not doi:
        return ""
    query = f'DOI:"{doi}" AND OPEN_ACCESS:Y'
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "resultType": "lite",
            "pageSize": 5,
        }
    )
    payload = methods_screen._get_json(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + params
    )
    for item in payload.get("resultList", {}).get("result", []):
        pmcid = str(item.get("pmcid", "") or "").strip()
        if pmcid:
            return pmcid
    return ""


@app.command("run")
def run(
    metadata_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_csv: Path = typer.Option(...),
) -> None:
    metadata = pd.read_csv(metadata_csv).fillna("")
    traits = pd.read_csv(trait_states_csv)
    config = methods_screen.load_config(config_path)
    lookup = methods_screen._species_lookup(traits)

    strict = metadata.loc[
        metadata["title"].astype(str).str.contains(STRICT_TITLE, na=False)
    ].copy()
    rows = []
    for record in strict.itertuples(index=False):
        doi = methods_screen._normalise_doi(record.doi)
        pmcid = _pmcid_for_doi(doi)
        if not pmcid:
            rows.append(
                {
                    "doi": doi,
                    "pmcid": "",
                    "publication_date": str(record.publication_date),
                    "title": str(record.title),
                    "methods_section_sha256": "",
                    "matched_species": "",
                    "coordinate_pairs_json": "[]",
                    "has_supplementation_term": False,
                    "has_natural_control_term": False,
                    "has_field_term": False,
                    "exclusion_flags": "",
                    "candidate_status": "no_open_access_pmc_match",
                    "screen_error": "",
                }
            )
            continue
        rows.append(
            methods_screen.screen_record(
                {
                    "doi": doi,
                    "pmcid": pmcid,
                    "publication_date": str(record.publication_date),
                    "title": str(record.title),
                },
                lookup=lookup,
                config=config,
            )
        )
        time.sleep(0.2)

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

    frame = pd.DataFrame(rows)
    typer.echo(f"strict_title_records={len(strict)}")
    if not frame.empty:
        typer.echo(frame["candidate_status"].value_counts(dropna=False).to_string())
        typer.echo(
            f"matched_species_records={frame['matched_species'].fillna('').ne('').sum()}"
        )


if __name__ == "__main__":
    app()
