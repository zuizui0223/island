"""Expand accepted plant names to auditable GBIF synonym search names.

The output is a discovery manifest, not a taxonomic decision table. Every alias
is tied back to the accepted species used by the analysis pipeline so evidence
found under an old name can be retained without changing the analysis taxon.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

GBIF_API = "https://api.gbif.org/v1"
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

ALIAS_COLUMNS = [
    "accepted_species",
    "search_name",
    "name_role",
    "gbif_usage_key",
    "gbif_accepted_key",
    "gbif_taxonomic_status",
    "gbif_canonical_name",
    "source",
]


def _text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def fetch_match(getter: JsonGetter, species: str) -> dict[str, Any]:
    return getter(f"{GBIF_API}/species/match", {"name": species, "strict": "true"})


def fetch_synonyms(getter: JsonGetter, key: str, limit: int) -> list[dict[str, Any]]:
    payload = getter(f"{GBIF_API}/species/{key}/synonyms", {"limit": limit})
    return list(payload.get("results") or [])[:limit]


def expand_species(
    getter: JsonGetter,
    species: str,
    *,
    max_synonyms: int,
) -> tuple[list[dict[str, str]], list[str]]:
    try:
        match = fetch_match(getter, species)
    except Exception as exc:  # noqa: BLE001
        return [], [f"match:{species}:{exc}"]

    usage_key = _text(match.get("usageKey") or match.get("speciesKey"))
    accepted_key = _text(match.get("acceptedUsageKey") or usage_key)
    rows = [
        {
            "accepted_species": species,
            "search_name": species,
            "name_role": "accepted_input",
            "gbif_usage_key": usage_key,
            "gbif_accepted_key": accepted_key,
            "gbif_taxonomic_status": _text(match.get("status")),
            "gbif_canonical_name": _text(match.get("canonicalName")),
            "source": "analysis_species_master",
        }
    ]
    if not accepted_key:
        return rows, []

    try:
        synonyms = fetch_synonyms(getter, accepted_key, max_synonyms)
    except Exception as exc:  # noqa: BLE001
        return rows, [f"synonyms:{species}:{exc}"]

    seen = {species.casefold()}
    for synonym in synonyms:
        name = _text(
            synonym.get("scientificName")
            or synonym.get("canonicalName")
            or synonym.get("name")
        )
        if not name or name.casefold() in seen:
            continue
        seen.add(name.casefold())
        rows.append(
            {
                "accepted_species": species,
                "search_name": name,
                "name_role": "gbif_synonym",
                "gbif_usage_key": _text(synonym.get("key") or synonym.get("usageKey")),
                "gbif_accepted_key": accepted_key,
                "gbif_taxonomic_status": _text(synonym.get("taxonomicStatus") or synonym.get("status")),
                "gbif_canonical_name": _text(synonym.get("canonicalName")),
                "source": "gbif_species_synonyms",
            }
        )
    return rows, []


def expand_names(
    species_df: pd.DataFrame,
    getter: JsonGetter,
    *,
    max_taxa: int,
    max_synonyms_per_species: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if "accepted_species" not in species_df.columns:
        raise typer.BadParameter("species table must have accepted_species")
    names = [
        name
        for name in dict.fromkeys(
            species_df["accepted_species"].fillna("").astype(str).str.strip()
        )
        if name
    ][:max_taxa]
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    species_with_synonyms: set[str] = set()
    for species in names:
        species_rows, species_errors = expand_species(
            getter,
            species,
            max_synonyms=max_synonyms_per_species,
        )
        rows.extend(species_rows)
        errors.extend(species_errors)
        if any(row["name_role"] == "gbif_synonym" for row in species_rows):
            species_with_synonyms.add(species)
    frame = pd.DataFrame(rows, columns=ALIAS_COLUMNS)
    report = {
        "source_lane": "gbif_taxon_name_expansion",
        "n_taxa_queried": len(names),
        "n_species_with_synonyms": len(species_with_synonyms),
        "n_search_names": len(frame),
        "n_synonym_names": int(frame["name_role"].eq("gbif_synonym").sum()) if len(frame) else 0,
        "n_lookup_errors": len(errors),
        "lookup_errors": errors[:50],
        "policy": (
            "Aliases expand discovery only. Evidence found under a synonym remains attached to the "
            "accepted_species from the analysis master; this output does not redefine taxonomy."
        ),
    }
    return frame, report


def _httpx_getter() -> JsonGetter:
    import httpx

    client = httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.1 (taxon-name-expansion)"},
    )

    def getter(url: str, params: dict[str, Any]) -> dict[str, Any]:
        response = client.get(url, params=params)
        response.raise_for_status()
        return response.json()

    return getter


@app.command("expand")
def expand_command(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(100, min=1, max=1000),
    max_synonyms_per_species: int = typer.Option(50, min=1, max=500),
) -> None:
    species_df = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = expand_names(
        species_df,
        _httpx_getter(),
        max_taxa=max_taxa,
        max_synonyms_per_species=max_synonyms_per_species,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "taxon_search_names.csv", index=False)
    (output_dir / "taxon_name_expansion_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
