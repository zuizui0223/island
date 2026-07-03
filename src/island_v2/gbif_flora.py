"""Rebuild exact-polygon GBIF island flowering-plant candidate lists.

Stages are deliberately separated:
1. prepare-islands standardizes an island polygon source.
2. make-requests builds one exact-WKT GBIF download request per island.
3. submit/poll/collect interacts with GBIF asynchronously and writes raw names.

A failed request is never silently treated as an empty island flora.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import time
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer
from shapely import make_valid
from shapely.geometry import MultiPolygon, Polygon

app = typer.Typer(add_completion=False, no_args_is_help=True)
GBIF_API = "https://api.gbif.org/v1"


def now_utc() -> str:
    return datetime.now(UTC).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def geometry_id(geometry: Any) -> str:
    return hashlib.sha256(geometry.wkb).hexdigest()[:20]


def wkt_hash(wkt: str) -> str:
    return hashlib.sha256(wkt.encode("utf-8")).hexdigest()


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise typer.BadParameter(f"{label} missing columns: {', '.join(sorted(missing))}")


def as_polygon_parts(geometry: Any) -> list[Polygon]:
    """Repair geometry and return polygon analysis units only."""
    repaired = make_valid(geometry)
    if isinstance(repaired, Polygon):
        return [repaired]
    if isinstance(repaired, MultiPolygon):
        return list(repaired.geoms)
    return [part for part in getattr(repaired, "geoms", []) if isinstance(part, Polygon)]


def read_islands(path: Path, layer: str | None) -> gpd.GeoDataFrame:
    try:
        return gpd.read_file(path, layer=layer)
    except Exception as exc:  # driver-specific details vary
        raise typer.BadParameter(f"Could not read island geometry source {path}: {exc}") from exc


def stable_island_table(
    source: gpd.GeoDataFrame,
    source_label: str,
    min_area_km2: float,
    name_field: str | None,
    parent_id_field: str | None,
) -> gpd.GeoDataFrame:
    """Convert an island layer to one stable row per polygon part."""
    if source.crs is None:
        raise typer.BadParameter("Island source has no CRS; provide a georeferenced vector file.")
    source = source.to_crs(4326).copy()
    if name_field is not None and name_field not in source.columns:
        raise typer.BadParameter(f"name_field {name_field!r} not present in source.")
    if parent_id_field is not None and parent_id_field not in source.columns:
        raise typer.BadParameter(f"parent_id_field {parent_id_field!r} not present in source.")

    rows: list[dict[str, Any]] = []
    for row_index, record in source.iterrows():
        if record.geometry is None or record.geometry.is_empty:
            continue
        parent_feature_id = str(record[parent_id_field]) if parent_id_field else str(row_index)
        source_name = str(record[name_field]) if name_field else ""
        for part_index, part in enumerate(as_polygon_parts(record.geometry), start=1):
            area_km2 = float(gpd.GeoSeries([part], crs=4326).to_crs(6933).area.iloc[0] / 1_000_000)
            if area_km2 < min_area_km2:
                continue
            rows.append(
                {
                    "island_id": f"{source_label}_{geometry_id(part)}",
                    "source_label": source_label,
                    "parent_feature_id": parent_feature_id,
                    "part_index": part_index,
                    "island_name": source_name,
                    "area_km2": area_km2,
                    "geometry_sha256": geometry_id(part),
                    "geometry": part,
                }
            )
    result = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
    if result.empty:
        raise typer.BadParameter("No polygon units remained after repair and area filtering.")
    if result["island_id"].duplicated().any():
        raise RuntimeError("Geometry-derived island IDs are not unique.")
    return result.sort_values("island_id").reset_index(drop=True)


def and_predicate(predicates: list[dict[str, Any]]) -> dict[str, Any]:
    return {"type": "and", "predicates": predicates}


def gbif_payload(
    wkt: str,
    taxon_key: str,
    checklist_key: str | None,
    format_name: str,
    email: str | None,
    minimum_year: int | None,
) -> dict[str, Any]:
    predicates: list[dict[str, Any]] = [
        {"type": "within", "geometry": wkt},
        {"type": "equals", "key": "TAXON_KEY", "value": str(taxon_key)},
        {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
        {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
    ]
    if checklist_key:
        predicates[1]["checklistKey"] = checklist_key
    if minimum_year is not None:
        predicates.append({"type": "greaterThanOrEquals", "key": "YEAR", "value": str(minimum_year)})
    payload: dict[str, Any] = {
        "sendNotification": bool(email),
        "format": format_name,
        "predicate": and_predicate(predicates),
    }
    if email:
        payload["notificationAddresses"] = [email]
    if checklist_key:
        payload["checklistKey"] = checklist_key
    return payload


def resolve_taxon(
    client: httpx.Client,
    scientific_name: str,
    rank: str | None,
    checklist_key: str | None,
) -> dict[str, Any]:
    params: dict[str, str] = {"name": scientific_name}
    if rank:
        params["rank"] = rank
    if checklist_key:
        params["checklistKey"] = checklist_key
    response = client.get(f"{GBIF_API}/species/match", params=params)
    response.raise_for_status()
    matched = response.json()
    usage_key = matched.get("usageKey") or matched.get("acceptedUsageKey") or matched.get("key")
    if not usage_key:
        raise RuntimeError(f"GBIF could not resolve taxon {scientific_name!r}: {matched}")
    return {
        "query_name": scientific_name,
        "query_rank": rank,
        "checklist_key": checklist_key,
        "usage_key": str(usage_key),
        "matched_name": matched.get("scientificName"),
        "matched_rank": matched.get("rank"),
        "match_type": matched.get("matchType"),
        "confidence": matched.get("confidence"),
    }


def request_rows(
    islands: gpd.GeoDataFrame,
    taxon: dict[str, Any],
    checklist_key: str | None,
    format_name: str,
    email: str | None,
    minimum_year: int | None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in islands.itertuples(index=False):
        wkt = record.geometry.wkt
        payload = gbif_payload(
            wkt=wkt,
            taxon_key=taxon["usage_key"],
            checklist_key=checklist_key,
            format_name=format_name,
            email=email,
            minimum_year=minimum_year,
        )
        rows.append(
            {
                "island_id": record.island_id,
                "area_km2": record.area_km2,
                "geometry_wkt": wkt,
                "geometry_wkt_sha256": wkt_hash(wkt),
                "payload_json": json.dumps(payload, ensure_ascii=False),
                "request_status": "prepared",
                "download_key": "",
                "error": "",
            }
        )
    return rows


def read_request_manifest(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str).fillna("")
    require_columns(
        table,
        {"island_id", "payload_json", "request_status", "download_key", "error"},
        "request manifest",
    )
    return table


def gbif_credentials() -> tuple[str, str]:
    username = os.getenv("GBIF_USERNAME", "")
    password = os.getenv("GBIF_PASSWORD", "")
    if not username or not password:
        raise typer.BadParameter("GBIF_USERNAME and GBIF_PASSWORD are required to submit downloads.")
    return username, password


def write_csv_atomic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def normalise_species_name(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def locate_tabular_file(archive: zipfile.ZipFile) -> str:
    candidates = [name for name in archive.namelist() if name.lower().endswith((".csv", ".tsv", ".txt"))]
    if not candidates:
        raise RuntimeError("GBIF archive has no delimited species-list file.")
    return sorted(candidates, key=lambda name: ("species" not in name.lower(), len(name)))[0]


def species_from_archive(path: Path, island_id: str, download_key: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with zipfile.ZipFile(path) as archive:
        member = locate_tabular_file(archive)
        with archive.open(member) as raw:
            text = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
            sample = text.read(4096)
            text.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            reader = csv.DictReader(text, dialect=dialect)
            fields = {field.lower(): field for field in (reader.fieldnames or [])}
            candidates = [
                fields.get("species"),
                fields.get("acceptedscientificname"),
                fields.get("scientificname"),
                fields.get("canonicalname"),
            ]
            name_fields = [field for field in candidates if field]
            if not name_fields:
                raise RuntimeError(f"No scientific-name field in {member}: {reader.fieldnames}")
            for row in reader:
                raw_name = next(
                    (normalise_species_name(row.get(field)) for field in name_fields if normalise_species_name(row.get(field))),
                    "",
                )
                if raw_name:
                    rows.append(
                        {
                            "island_id": island_id,
                            "raw_scientific_name": raw_name,
                            "gbif_download_key": download_key,
                        }
                    )
    return rows


@app.command("prepare-islands")
def prepare_islands(
    input_path: Path = typer.Option(..., exists=True, help="Island vector dataset: FileGDB, GPKG, or shapefile."),
    output_dir: Path = typer.Option(..., help="Directory for standardized v2 island polygons."),
    source_label: str = typer.Option("global_islands", help="Stable label for this source release."),
    layer: str | None = typer.Option(None, help="Layer name for GDB/GPKG input."),
    min_area_km2: float = typer.Option(20.0, min=0.0),
    name_field: str | None = typer.Option(None),
    parent_id_field: str | None = typer.Option(None),
) -> None:
    """Standardize exact island polygons and write a source manifest."""
    source = read_islands(input_path, layer)
    islands = stable_island_table(source, source_label, min_area_km2, name_field, parent_id_field)
    output_dir.mkdir(parents=True, exist_ok=True)
    gpkg_path = output_dir / "islands_v2.gpkg"
    islands.to_file(gpkg_path, layer="islands", driver="GPKG")
    manifest = islands.drop(columns="geometry").copy()
    manifest["source_path"] = str(input_path)
    manifest["source_sha256"] = sha256_file(input_path) if input_path.is_file() else "directory_source"
    manifest["prepared_at_utc"] = now_utc()
    manifest.to_csv(output_dir / "island_manifest.csv", index=False)
    typer.echo(f"Prepared {len(islands)} island polygon units at {gpkg_path}")


@app.command("make-requests")
def make_requests(
    islands_gpkg: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
    taxon_name: str = typer.Option("Angiosperms"),
    taxon_rank: str | None = typer.Option(None),
    checklist_key: str | None = typer.Option(None),
    format_name: str = typer.Option("SPECIES_LIST"),
    minimum_year: int | None = typer.Option(None),
    email: str | None = typer.Option(None, help="GBIF notification address; omit for no email."),
) -> None:
    """Resolve the target clade and create exact-polygon GBIF download payloads."""
    islands = gpd.read_file(islands_gpkg, layer="islands")
    require_columns(islands, {"island_id", "area_km2"}, "island GPKG")
    with httpx.Client(timeout=60.0, headers={"User-Agent": "island-floral-v2/0.1"}) as client:
        taxon = resolve_taxon(client, taxon_name, taxon_rank, checklist_key)
    rows = request_rows(islands, taxon, checklist_key, format_name, email, minimum_year)
    frame = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_csv, index=False)
    resolution_path = output_csv.with_name(output_csv.stem + "_taxon_resolution.json")
    resolution_path.write_text(json.dumps(taxon, indent=2, ensure_ascii=False), encoding="utf-8")
    typer.echo(f"Prepared {len(frame)} exact-polygon GBIF requests at {output_csv}")


@app.command("submit")
def submit(
    request_manifest: Path = typer.Option(..., exists=True),
    max_requests: int = typer.Option(20, min=1, max=100),
    sleep_seconds: float = typer.Option(1.0, min=0.0),
) -> None:
    """Submit a bounded batch; request errors remain visible in the manifest."""
    username, password = gbif_credentials()
    table = read_request_manifest(request_manifest)
    todo = table.index[table["request_status"].eq("prepared")].tolist()[:max_requests]
    if not todo:
        typer.echo("No prepared requests remain.")
        return
    with httpx.Client(timeout=60.0, auth=(username, password), headers={"User-Agent": "island-floral-v2/0.1"}) as client:
        for index in todo:
            try:
                response = client.post(
                    f"{GBIF_API}/occurrence/download/request",
                    json=json.loads(table.at[index, "payload_json"]),
                )
                response.raise_for_status()
                table.at[index, "download_key"] = response.text.strip().strip('"')
                table.at[index, "request_status"] = "submitted"
                table.at[index, "error"] = ""
            except Exception as exc:
                table.at[index, "request_status"] = "submit_failed"
                table.at[index, "error"] = str(exc)
            write_csv_atomic(table, request_manifest)
            time.sleep(sleep_seconds)
    typer.echo(f"Submitted {len(todo)} GBIF requests; inspect {request_manifest}.")


@app.command("poll")
def poll(
    request_manifest: Path = typer.Option(..., exists=True),
) -> None:
    """Poll submitted keys and record status, DOI, and archive link."""
    table = read_request_manifest(request_manifest)
    if "doi" not in table.columns:
        table["doi"] = ""
    if "download_link" not in table.columns:
        table["download_link"] = ""
    pending = table.index[table["request_status"].isin(["submitted", "running"])].tolist()
    with httpx.Client(timeout=60.0, headers={"User-Agent": "island-floral-v2/0.1"}) as client:
        for index in pending:
            key = table.at[index, "download_key"]
            if not key:
                continue
            try:
                response = client.get(f"{GBIF_API}/occurrence/download/{key}")
                response.raise_for_status()
                status = response.json()
                table.at[index, "request_status"] = str(status.get("status", "unknown")).lower()
                table.at[index, "doi"] = str(status.get("doi", ""))
                table.at[index, "download_link"] = str(status.get("downloadLink", ""))
                table.at[index, "error"] = ""
            except Exception as exc:
                table.at[index, "error"] = str(exc)
    write_csv_atomic(table, request_manifest)
    typer.echo(f"Polled {len(pending)} GBIF downloads.")


@app.command("collect")
def collect(
    request_manifest: Path = typer.Option(..., exists=True),
    download_dir: Path = typer.Option(...),
    output_species_csv: Path = typer.Option(...),
) -> None:
    """Collect successful species-list archives into island-by-raw-name candidates."""
    table = read_request_manifest(request_manifest)
    require_columns(table, {"download_link"}, "request manifest")
    ready = table.loc[table["request_status"].eq("succeeded")].copy()
    if ready.empty:
        typer.echo("No succeeded GBIF downloads are ready to collect.")
        return
    download_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    with httpx.Client(timeout=300.0, follow_redirects=True, headers={"User-Agent": "island-floral-v2/0.1"}) as client:
        for record in ready.itertuples(index=False):
            if not record.download_link:
                continue
            archive_path = download_dir / f"{record.island_id}_{record.download_key}.zip"
            if not archive_path.exists():
                response = client.get(record.download_link)
                response.raise_for_status()
                archive_path.write_bytes(response.content)
            try:
                rows.extend(species_from_archive(archive_path, record.island_id, record.download_key))
            except Exception as exc:
                typer.echo(f"ERROR parsing {record.island_id} / {record.download_key}: {exc}")
    species = pd.DataFrame(rows).drop_duplicates()
    if not species.empty:
        species["acquisition_stage"] = "gbif_exact_polygon_species_list"
        species["review_status"] = "unresolved_taxonomy"
    output_species_csv.parent.mkdir(parents=True, exist_ok=True)
    species.to_csv(output_species_csv, index=False)
    typer.echo(f"Wrote {len(species)} raw island-species candidates to {output_species_csv}")


if __name__ == "__main__":
    app()
