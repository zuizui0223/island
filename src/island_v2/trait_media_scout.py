"""Collect taxon-matched GBIF occurrence images for visible floral-trait recovery.

This module only acquires media evidence and provenance. It does not infer a
trait value. Images are later scored by a human or vision model and compared
against other recovery lanes with ``trait_recovery_benchmark``.
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

MEDIA_COLUMNS = [
    "accepted_species",
    "gbif_taxon_key",
    "occurrence_key",
    "scientific_name",
    "basis_of_record",
    "country_code",
    "event_date",
    "media_identifier",
    "media_type",
    "media_format",
    "media_title",
    "media_description",
    "media_creator",
    "media_license",
    "media_reference",
]


def _text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    if isinstance(value, (list, tuple, set)):
        return " | ".join(_text(item) for item in value if _text(item))
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return " ".join(str(value).strip().split())


def fetch_taxon_key(getter: JsonGetter, species: str) -> str:
    payload = getter(f"{GBIF_API}/species/match", {"name": species, "strict": "true"})
    key = payload.get("usageKey") or payload.get("speciesKey") or ""
    return str(key) if key else ""


def fetch_occurrences(
    getter: JsonGetter,
    taxon_key: str,
    limit: int,
) -> list[dict[str, Any]]:
    payload = getter(
        f"{GBIF_API}/occurrence/search",
        {
            "taxon_key": taxon_key,
            "media_type": "StillImage",
            "occurrence_status": "PRESENT",
            "limit": limit,
        },
    )
    return list(payload.get("results") or [])


def _valid_media(occurrence: dict[str, Any]) -> list[dict[str, Any]]:
    valid: list[dict[str, Any]] = []
    for media in occurrence.get("media") or []:
        identifier = _text(media.get("identifier"))
        media_type = _text(media.get("type"))
        if not identifier:
            continue
        if media_type and media_type.lower() not in {"stillimage", "image"}:
            continue
        valid.append(media)
    return valid


def _media_row(
    species: str,
    taxon_key: str,
    occurrence: dict[str, Any],
    media: dict[str, Any],
) -> dict[str, str]:
    return {
        "accepted_species": species,
        "gbif_taxon_key": taxon_key,
        "occurrence_key": _text(occurrence.get("key")),
        "scientific_name": _text(occurrence.get("scientificName") or occurrence.get("species")),
        "basis_of_record": _text(occurrence.get("basisOfRecord")),
        "country_code": _text(occurrence.get("countryCode")),
        "event_date": _text(occurrence.get("eventDate")),
        "media_identifier": _text(media.get("identifier")),
        "media_type": _text(media.get("type")),
        "media_format": _text(media.get("format")),
        "media_title": _text(media.get("title")),
        "media_description": _text(media.get("description")),
        "media_creator": _text(media.get("creator")),
        "media_license": _text(media.get("license")),
        "media_reference": _text(media.get("references")),
    }


def media_rows_for_species(
    getter: JsonGetter,
    species: str,
    *,
    max_occurrences: int,
    max_images: int,
) -> tuple[list[dict[str, str]], list[str]]:
    errors: list[str] = []
    try:
        taxon_key = fetch_taxon_key(getter, species)
    except Exception as exc:  # noqa: BLE001
        return [], [f"match:{species}:{exc}"]
    if not taxon_key:
        return [], []
    try:
        occurrences = fetch_occurrences(getter, taxon_key, max_occurrences)
    except Exception as exc:  # noqa: BLE001
        return [], [f"occurrence:{species}:{exc}"]

    # Round-robin across occurrences: take one image from every observation before
    # taking a second image from any observation. GBIF observations often contain
    # several near-duplicate views (for example five leaf photos from one plant),
    # so this maximizes independent biological/photographic opportunities for a
    # visible flower within a fixed image budget.
    occurrence_media = [(occurrence, _valid_media(occurrence)) for occurrence in occurrences]
    occurrence_media = [(occurrence, media) for occurrence, media in occurrence_media if media]
    rows: list[dict[str, str]] = []
    seen_urls: set[str] = set()
    media_index = 0
    while len(rows) < max_images:
        added = False
        for occurrence, media_items in occurrence_media:
            if media_index >= len(media_items):
                continue
            media = media_items[media_index]
            identifier = _text(media.get("identifier"))
            if identifier in seen_urls:
                continue
            seen_urls.add(identifier)
            rows.append(_media_row(species, taxon_key, occurrence, media))
            added = True
            if len(rows) >= max_images:
                break
        if not added and all(media_index + 1 >= len(items) for _, items in occurrence_media):
            break
        media_index += 1
    return rows, errors


def collect_media(
    species_df: pd.DataFrame,
    getter: JsonGetter,
    *,
    max_taxa: int,
    max_occurrences: int,
    max_images_per_species: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if "accepted_species" not in species_df.columns:
        raise typer.BadParameter("species table must have accepted_species")
    names = [
        name
        for name in dict.fromkeys(species_df["accepted_species"].fillna("").astype(str).str.strip())
        if name
    ][:max_taxa]
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    species_with_media: set[str] = set()
    for species in names:
        species_rows, species_errors = media_rows_for_species(
            getter,
            species,
            max_occurrences=max_occurrences,
            max_images=max_images_per_species,
        )
        rows.extend(species_rows)
        errors.extend(species_errors)
        if species_rows:
            species_with_media.add(species)
    frame = pd.DataFrame(rows, columns=MEDIA_COLUMNS)
    report = {
        "source_lane": "gbif_occurrence_media",
        "selection_policy": "round_robin_distinct_occurrences_before_additional_images",
        "n_taxa_queried": len(names),
        "n_species_with_media": len(species_with_media),
        "n_species_without_media": len(names) - len(species_with_media),
        "n_media_rows": len(frame),
        "max_images_per_species": max_images_per_species,
        "n_lookup_errors": len(errors),
        "lookup_errors": errors[:50],
        "policy": (
            "Media rows are taxon-matched visual evidence only. No floral trait is inferred here; "
            "image availability is not trait availability and missing media is not biological absence."
        ),
    }
    return frame, report


def _httpx_getter() -> JsonGetter:
    import httpx

    client = httpx.Client(
        timeout=45.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.1 (trait-media-scout)"},
    )

    def getter(url: str, params: dict[str, Any]) -> dict[str, Any]:
        response = client.get(url, params=params)
        response.raise_for_status()
        return response.json()

    return getter


@app.command("collect")
def collect_command(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(50, min=1, max=500),
    max_occurrences: int = typer.Option(100, min=1, max=300),
    max_images_per_species: int = typer.Option(5, min=1, max=20),
) -> None:
    species_df = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = collect_media(
        species_df,
        _httpx_getter(),
        max_taxa=max_taxa,
        max_occurrences=max_occurrences,
        max_images_per_species=max_images_per_species,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "gbif_trait_media_manifest.csv", index=False)
    (output_dir / "gbif_trait_media_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
