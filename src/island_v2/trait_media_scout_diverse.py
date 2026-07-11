"""Collect diversified GBIF image evidence for visible floral-trait recovery.

Unlike the baseline media scout, this lane selects at most one image per
occurrence and prefers observations over preserved specimens. This prevents a
single observation with several near-duplicate images from consuming the whole
per-species budget.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.trait_media_scout import (
    MEDIA_COLUMNS,
    _httpx_getter,
    _text,
    fetch_occurrences,
    fetch_taxon_key,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
JsonGetter = Callable[[str, dict[str, Any]], dict[str, Any]]

FLOWER_TERMS = re.compile(
    r"\b(flower|flowers|flowering|floral|bloom|blossom|inflorescence|corolla)\b",
    re.I,
)
BASIS_PRIORITY = {
    "HUMAN_OBSERVATION": 0,
    "OBSERVATION": 1,
    "MACHINE_OBSERVATION": 2,
    "MATERIAL_SAMPLE": 3,
    "PRESERVED_SPECIMEN": 4,
}


def _occurrence_score(occurrence: dict[str, Any]) -> tuple[int, int, str]:
    basis = _text(occurrence.get("basisOfRecord"))
    media_text = " ".join(
        _text(value)
        for media in occurrence.get("media") or []
        for value in (media.get("title"), media.get("description"))
    )
    flower_rank = 0 if FLOWER_TERMS.search(media_text) else 1
    basis_rank = BASIS_PRIORITY.get(basis, 5)
    return flower_rank, basis_rank, _text(occurrence.get("key"))


def media_rows_for_species_diverse(
    getter: JsonGetter,
    species: str,
    *,
    max_occurrences: int,
    max_images: int,
) -> tuple[list[dict[str, str]], list[str]]:
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

    rows: list[dict[str, str]] = []
    seen_urls: set[str] = set()
    for occurrence in sorted(occurrences, key=_occurrence_score):
        chosen: dict[str, Any] | None = None
        for media in occurrence.get("media") or []:
            identifier = _text(media.get("identifier"))
            media_type = _text(media.get("type"))
            if not identifier or identifier in seen_urls:
                continue
            if media_type and media_type.lower() not in {"stillimage", "image"}:
                continue
            chosen = media
            break
        if chosen is None:
            continue
        identifier = _text(chosen.get("identifier"))
        seen_urls.add(identifier)
        rows.append(
            {
                "accepted_species": species,
                "gbif_taxon_key": taxon_key,
                "occurrence_key": _text(occurrence.get("key")),
                "scientific_name": _text(
                    occurrence.get("scientificName") or occurrence.get("species")
                ),
                "basis_of_record": _text(occurrence.get("basisOfRecord")),
                "country_code": _text(occurrence.get("countryCode")),
                "event_date": _text(occurrence.get("eventDate")),
                "media_identifier": identifier,
                "media_type": _text(chosen.get("type")),
                "media_format": _text(chosen.get("format")),
                "media_title": _text(chosen.get("title")),
                "media_description": _text(chosen.get("description")),
                "media_creator": _text(chosen.get("creator")),
                "media_license": _text(chosen.get("license")),
                "media_reference": _text(chosen.get("references")),
            }
        )
        if len(rows) >= max_images:
            break
    return rows, []


def collect_diverse_media(
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
        for name in dict.fromkeys(
            species_df["accepted_species"].fillna("").astype(str).str.strip()
        )
        if name
    ][:max_taxa]
    rows: list[dict[str, str]] = []
    errors: list[str] = []
    species_with_media: set[str] = set()
    for species in names:
        species_rows, species_errors = media_rows_for_species_diverse(
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
        "source_lane": "gbif_occurrence_media_diverse",
        "selection_strategy": "flower_metadata_then_basis_priority_one_image_per_occurrence",
        "n_taxa_queried": len(names),
        "n_species_with_media": len(species_with_media),
        "n_species_without_media": len(names) - len(species_with_media),
        "n_media_rows": len(frame),
        "max_images_per_species": max_images_per_species,
        "n_lookup_errors": len(errors),
        "lookup_errors": errors[:50],
        "policy": (
            "Selection diversifies visual evidence but does not infer traits. Missing media and "
            "image content remain acquisition outcomes, not biological absences."
        ),
    }
    return frame, report


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    max_taxa: int = typer.Option(50, min=1, max=500),
    max_occurrences: int = typer.Option(300, min=1, max=300),
    max_images_per_species: int = typer.Option(5, min=1, max=20),
) -> None:
    species_df = pd.read_csv(species_csv, dtype=str).fillna("")
    frame, report = collect_diverse_media(
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
