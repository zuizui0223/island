"""Fetch the latest open AusTraits plain-text release and standardize one dataset.

The official AusTraits R client resolves versions through the Zenodo concept-record
versions endpoint. This module follows the same public API without requiring R. It
prefers a ZIP archive containing CSV tables, discovers a long-format trait table by
column contract, and writes a compact standardized CSV for a requested dataset.

No trait values are inferred, recoded, or assigned units here.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Iterable

VERSIONS_API = "https://zenodo.org/api/records/3568418/versions?size=100"
USER_AGENT = "island-v2-free-bulk-trait-pilot/0.1"
TAXON_ALIASES = ("taxon_name", "accepted_name", "scientific_name", "species_name", "taxon")
TRAIT_ALIASES = ("trait_name", "trait")
VALUE_ALIASES = ("trait_value", "value")
UNIT_ALIASES = ("unit", "trait_unit", "unit_name")
DATASET_ALIASES = ("dataset_id", "dataset_key", "dataset", "source_id")


def _request_json(url: str) -> dict[str, object]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def _download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(request, timeout=300) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def _release_hits(payload: dict[str, object]) -> list[dict[str, object]]:
    hits = payload.get("hits", {})
    if not isinstance(hits, dict):
        return []
    raw = hits.get("hits", [])
    return [item for item in raw if isinstance(item, dict)] if isinstance(raw, list) else []


def select_latest_release(payload: dict[str, object]) -> dict[str, object]:
    hits = _release_hits(payload)
    if not hits:
        raise ValueError("Zenodo versions response contains no AusTraits releases")

    def key(item: dict[str, object]) -> tuple[str, str]:
        metadata = item.get("metadata", {})
        if not isinstance(metadata, dict):
            metadata = {}
        return (str(metadata.get("publication_date", "")), str(metadata.get("version", "")))

    return max(hits, key=key)


def select_plain_text_archive(release: dict[str, object]) -> dict[str, str]:
    files = release.get("files", [])
    if not isinstance(files, list):
        raise ValueError("AusTraits release has no file list")
    candidates: list[dict[str, str]] = []
    for item in files:
        if not isinstance(item, dict):
            continue
        key = str(item.get("key", ""))
        links = item.get("links", {})
        if not isinstance(links, dict):
            continue
        url = str(links.get("self", ""))
        if key.lower().endswith(".zip") and url:
            candidates.append({"key": key, "url": url})
    if not candidates:
        raise ValueError("latest AusTraits release contains no ZIP plain-text archive")

    def score(item: dict[str, str]) -> tuple[int, int, str]:
        name = item["key"].lower()
        return (
            int("csv" in name or "text" in name or "flat" in name),
            int("austraits" in name),
            name,
        )

    return max(candidates, key=score)


def _header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.reader(handle)
        return next(reader, [])


def _first_present(fields: Iterable[str], aliases: Iterable[str]) -> str | None:
    available = {field.strip() for field in fields}
    return next((alias for alias in aliases if alias in available), None)


def discover_long_trait_table(root: Path) -> tuple[Path, dict[str, str | None]]:
    candidates: list[tuple[int, Path, dict[str, str | None]]] = []
    for path in root.rglob("*.csv"):
        fields = _header(path)
        taxon = _first_present(fields, TAXON_ALIASES)
        trait = _first_present(fields, TRAIT_ALIASES)
        value = _first_present(fields, VALUE_ALIASES)
        if not taxon or not trait or not value:
            continue
        candidates.append((path.stat().st_size, path, {
            "taxon": taxon,
            "trait": trait,
            "value": value,
            "unit": _first_present(fields, UNIT_ALIASES),
            "dataset": _first_present(fields, DATASET_ALIASES),
        }))
    if not candidates:
        raise ValueError("no long-format trait CSV with taxon/trait/value columns found in AusTraits archive")
    _, path, columns = max(candidates, key=lambda item: item[0])
    return path, columns


def standardize_dataset(
    source_csv: Path,
    output_csv: Path,
    *,
    columns: dict[str, str | None],
    dataset_key: str | None,
) -> dict[str, object]:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    input_rows = 0
    written_rows = 0
    source_taxa: set[str] = set()
    source_traits: set[str] = set()
    source_units: set[str] = set()
    dataset_column = columns.get("dataset")
    unit_column = columns.get("unit")
    if dataset_key and not dataset_column:
        raise ValueError(
            f"dataset filter {dataset_key!r} requested but discovered trait table has no dataset column"
        )

    with source_csv.open(newline="", encoding="utf-8-sig", errors="replace") as source, output_csv.open(
        "w", newline="", encoding="utf-8"
    ) as target:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(
            target,
            fieldnames=[
                "taxon_name",
                "trait_name",
                "trait_value",
                "trait_unit",
                "dataset_key",
                "source_record_id",
            ],
        )
        writer.writeheader()
        for index, row in enumerate(reader, start=1):
            input_rows += 1
            if dataset_key and str(row.get(str(dataset_column), "")).strip() != dataset_key:
                continue
            taxon = " ".join(str(row.get(str(columns["taxon"]), "")).split())
            trait = str(row.get(str(columns["trait"]), "")).strip()
            value = str(row.get(str(columns["value"]), "")).strip()
            unit = str(row.get(str(unit_column), "")).strip() if unit_column else ""
            if not taxon or not trait or not value:
                continue
            source_taxa.add(taxon)
            source_traits.add(trait)
            if unit:
                source_units.add(unit)
            writer.writerow({
                "taxon_name": taxon,
                "trait_name": trait,
                "trait_value": value,
                "trait_unit": unit,
                "dataset_key": dataset_key or str(row.get(str(dataset_column), "")).strip(),
                "source_record_id": f"austraits:{dataset_key or 'all'}:{index}",
            })
            written_rows += 1

    return {
        "input_rows_in_discovered_table": input_rows,
        "standardized_rows": written_rows,
        "unique_taxa": len(source_taxa),
        "unique_traits": len(source_traits),
        "unique_nonempty_units": len(source_units),
    }


def run_fetch(
    output_csv: Path,
    manifest_json: Path,
    *,
    dataset_key: str | None,
    versions_api: str = VERSIONS_API,
) -> dict[str, object]:
    versions = _request_json(versions_api)
    release = select_latest_release(versions)
    archive = select_plain_text_archive(release)
    metadata = release.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}

    with tempfile.TemporaryDirectory(prefix="austraits_bulk_") as tmp:
        tmpdir = Path(tmp)
        archive_path = tmpdir / archive["key"]
        _download(archive["url"], archive_path)
        extract_dir = tmpdir / "extracted"
        extract_dir.mkdir()
        with zipfile.ZipFile(archive_path) as zipped:
            zipped.extractall(extract_dir)
        trait_table, columns = discover_long_trait_table(extract_dir)
        stats = standardize_dataset(
            trait_table,
            output_csv,
            columns=columns,
            dataset_key=dataset_key,
        )

    manifest = {
        "source": "AusTraits",
        "versions_api": versions_api,
        "release_version": metadata.get("version"),
        "publication_date": metadata.get("publication_date"),
        "doi": metadata.get("doi"),
        "archive_key": archive["key"],
        "dataset_key": dataset_key,
        "discovered_columns": columns,
        **stats,
    }
    manifest_json.parent.mkdir(parents=True, exist_ok=True)
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--dataset-key", default="eFLOWER_Dun_2022")
    parser.add_argument("--versions-api", default=VERSIONS_API)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_fetch(
        args.output_csv,
        args.manifest_json,
        dataset_key=args.dataset_key or None,
        versions_api=args.versions_api,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
