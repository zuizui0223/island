"""Collect auditable GBIF occurrence-media manifests for visible floral traits.

This module does not classify images. It measures and freezes the availability of
species-linked still images so a later vision/manual step can extract flower colour,
form and symmetry with image-level provenance.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

GBIF_OCCURRENCE_SEARCH = "https://api.gbif.org/v1/occurrence/search"
USER_AGENT = "island-floral-v2/0.1 GBIF visual evidence collector"
OUTPUT_COLUMNS = [
    "accepted_species",
    "gbif_occurrence_key",
    "scientific_name",
    "taxon_key",
    "basis_of_record",
    "country_code",
    "event_date",
    "media_identifier",
    "media_type",
    "media_format",
    "media_title",
    "media_description",
    "media_created",
    "media_creator",
    "media_license",
    "media_reference",
    "publisher",
    "dataset_key",
    "record_url",
    "evidence_scope",
    "media_id_sha256",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _media_rows(record: dict[str, object], species: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    occurrence_key = _text(record.get("key"))
    media_items = record.get("media") or []
    if not isinstance(media_items, list):
        return rows
    for media in media_items:
        if not isinstance(media, dict):
            continue
        identifier = _text(media.get("identifier"))
        media_type = _text(media.get("type"))
        if not identifier or media_type.lower() not in {"stillimage", "image"}:
            continue
        rows.append(
            {
                "accepted_species": species,
                "gbif_occurrence_key": occurrence_key,
                "scientific_name": _text(record.get("scientificName")),
                "taxon_key": _text(record.get("taxonKey")),
                "basis_of_record": _text(record.get("basisOfRecord")),
                "country_code": _text(record.get("countryCode")),
                "event_date": _text(record.get("eventDate")),
                "media_identifier": identifier,
                "media_type": media_type,
                "media_format": _text(media.get("format")),
                "media_title": _text(media.get("title")),
                "media_description": _text(media.get("description")),
                "media_created": _text(media.get("created")),
                "media_creator": _text(media.get("creator")),
                "media_license": _text(media.get("license")),
                "media_reference": _text(media.get("references")),
                "publisher": _text(record.get("publishingOrgKey")),
                "dataset_key": _text(record.get("datasetKey")),
                "record_url": f"https://www.gbif.org/occurrence/{occurrence_key}",
                "evidence_scope": "species_linked_occurrence_image",
                "media_id_sha256": _sha(identifier),
            }
        )
    return rows


def collect_species(
    client: httpx.Client,
    species: str,
    *,
    max_images: int = 10,
    page_limit: int = 100,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    offset = 0
    seen: set[str] = set()
    while len(rows) < max_images and offset < 500:
        response = client.get(
            GBIF_OCCURRENCE_SEARCH,
            params={
                "scientific_name": species,
                "media_type": "StillImage",
                "limit": page_limit,
                "offset": offset,
            },
            timeout=30,
        )
        response.raise_for_status()
        payload = response.json()
        results = payload.get("results") or []
        for record in results:
            for row in _media_rows(record, species):
                identifier = row["media_identifier"]
                if identifier in seen:
                    continue
                seen.add(identifier)
                rows.append(row)
                if len(rows) >= max_images:
                    break
            if len(rows) >= max_images:
                break
        if payload.get("endOfRecords") or not results:
            break
        offset += len(results)
    return rows


@app.command()
def collect(
    species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    species_column: str = typer.Option("accepted_species"),
    max_taxa: int = typer.Option(25, min=1, max=500),
    max_images_per_species: int = typer.Option(10, min=1, max=50),
) -> None:
    table = pd.read_csv(species_csv, dtype=str).fillna("")
    if species_column not in table.columns:
        raise typer.BadParameter(f"missing species column: {species_column}")
    species: list[str] = []
    seen_species: set[str] = set()
    for value in table[species_column]:
        name = _text(value)
        if name and name not in seen_species:
            seen_species.add(name)
            species.append(name)
        if len(species) >= max_taxa:
            break

    rows: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    with httpx.Client(headers={"User-Agent": USER_AGENT}) as client:
        for name in species:
            try:
                rows.extend(
                    collect_species(client, name, max_images=max_images_per_species)
                )
            except Exception as exc:
                errors.append({"accepted_species": name, "error": str(exc)})

    manifest = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    errors_frame = pd.DataFrame(errors, columns=["accepted_species", "error"])
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_dir / "gbif_visual_evidence_manifest.csv", index=False)
    errors_frame.to_csv(output_dir / "gbif_visual_evidence_errors.csv", index=False)
    counts = (
        manifest.groupby("accepted_species").size().sort_values(ascending=False).to_dict()
        if len(manifest)
        else {}
    )
    summary = {
        "n_species_requested": len(species),
        "n_species_with_images": int(manifest["accepted_species"].nunique()) if len(manifest) else 0,
        "n_images": len(manifest),
        "max_images_per_species": max_images_per_species,
        "n_species_errors": len(errors_frame),
        "images_per_species": counts,
        "policy": "Image availability only; no trait state is inferred by this collector.",
    }
    (output_dir / "gbif_visual_evidence_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    app()
