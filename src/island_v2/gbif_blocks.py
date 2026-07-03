"""Prepare bounded regional GBIF occurrence downloads for the v2 island universe.

The island analysis geometry is never replaced. For GBIF acquisition only, a
small outward query catchment is built around each island to avoid missing
coastal records because of coordinate rounding or shoreline mismatch. Returned
records must subsequently be assigned back to the original island polygons.

This avoids both unsuitable extremes:
- one download request per island; and
- continental bounding-box downloads that make mainland records look like
  island data.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer
from shapely import make_valid, orient_polygons, to_wkt
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box
from shapely.ops import transform, unary_union

from island_v2.gbif_flora import GBIF_API, gbif_payload, require_columns, resolve_taxon

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Build spatially bounded, reproducible GBIF request blocks."""


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def geometry_hash(geometry: Any) -> str:
    return hashlib.sha256(geometry.wkb).hexdigest()


def polygonal_only(geometry: Any) -> Polygon | MultiPolygon:
    """Return a valid polygonal geometry or raise rather than silently changing scope."""
    repaired = make_valid(geometry)
    if isinstance(repaired, (Polygon, MultiPolygon)):
        return repaired
    if isinstance(repaired, GeometryCollection):
        polygons = [part for part in repaired.geoms if isinstance(part, Polygon)]
        if polygons:
            return MultiPolygon(polygons) if len(polygons) > 1 else polygons[0]
    raise ValueError(f"Expected a polygonal island geometry, got {repaired.geom_type}")


def canonical_gbif_polygon(geometry: Any) -> Polygon | MultiPolygon:
    """Return valid GBIF-compatible polygon rings.

    GBIF rejects WKT polygons whose exterior rings are clockwise.  All query
    geometries are therefore canonicalized to counter-clockwise exteriors and
    clockwise interior rings before their WKT, hash, or payload is generated.
    This does not alter the spatial footprint.
    """
    return orient_polygons(polygonal_only(geometry), exterior_cw=False)


def stable_wkt(geometry: Any) -> str:
    """Use fixed precision and GBIF-compatible ring order for request provenance."""
    return to_wkt(
        canonical_gbif_polygon(geometry),
        rounding_precision=5,
        trim=True,
        output_dimension=2,
    )


def _wrap_east(x: float, y: float, z: float | None = None) -> tuple[float, ...]:
    shifted = x + 360.0 if x < 0 else x
    return (shifted, y) if z is None else (shifted, y, z)


def _wrap_west(x: float, y: float, z: float | None = None) -> tuple[float, ...]:
    # >= (not >) so the seam vertex itself (exactly 180 after _wrap_east) moves
    # back to -180 along with the rest of this part; otherwise it is left
    # stranded at +180 and the "west-shifted" part still spans ~360 degrees.
    shifted = x - 360.0 if x >= 180.0 else x
    return (shifted, y) if z is None else (shifted, y, z)


def split_at_antimeridian(geometry: Polygon | MultiPolygon) -> list[Polygon | MultiPolygon]:
    """Split a query geometry that spans more than 180 degrees of longitude.

    GBIF's `within` WKT predicate has no antimeridian support: a ring with
    vertices on both sides of +-180 degrees is rejected outright, even when
    the underlying island is small and simply happens to straddle the seam
    (a known GSHHG artifact for some Pacific/Arctic landmasses, and also
    possible after buffering/unioning nearby islands into one query
    catchment). This recentres the geometry on the seam, cuts it there, and
    returns each ordinary-range remainder separately so it can be queried on
    its own. Non-crossing geometries are returned unchanged.
    """
    minx, _, maxx, _ = geometry.bounds
    if maxx - minx <= 180.0:
        return [geometry]

    shifted = transform(_wrap_east, geometry)
    sminx, _, smaxx, _ = shifted.bounds
    if smaxx - sminx > 180.0:
        raise ValueError(
            f"Query geometry spans {maxx - minx:.1f} degrees of longitude even after "
            "antimeridian recentring; refusing to build an invalid GBIF WKT from it."
        )

    huge = 1_000.0
    parts: list[Polygon | MultiPolygon] = []
    west_of_seam = shifted.intersection(box(-huge, -huge, 180.0, huge))
    if not west_of_seam.is_empty:
        parts.append(polygonal_only(make_valid(west_of_seam)))
    east_of_seam = shifted.intersection(box(180.0, -huge, huge, huge))
    if not east_of_seam.is_empty:
        parts.append(polygonal_only(make_valid(transform(_wrap_west, east_of_seam))))
    if not parts:
        raise ValueError("Antimeridian split of a query geometry produced no polygonal remainder.")
    return parts


def explode_antimeridian_catchments(query: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Replace any seam-crossing catchment with its per-side, GBIF-safe parts.

    `island_id` stays unique (required by `block_query_geometries`) by
    suffixing split parts; `analysis_island_id` keeps the link back to the
    single original island so coverage/QA can still be verified per island.
    """
    rows: list[dict[str, Any]] = []
    for record in query.to_dict("records"):
        parts = split_at_antimeridian(record["geometry"])
        if len(parts) == 1:
            row = dict(record)
            row["geometry"] = parts[0]
            row["analysis_island_id"] = record["island_id"]
            rows.append(row)
            continue
        for part_index, part in enumerate(parts, start=1):
            row = dict(record)
            row["geometry"] = canonical_gbif_polygon(part)
            row["analysis_island_id"] = record["island_id"]
            row["island_id"] = f"{record['island_id']}__antimeridian_part{part_index}"
            rows.append(row)
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=query.crs)


def grid_index(longitude: float, latitude: float, grid_degrees: float) -> tuple[int, int]:
    """Return stable 0-based 360 x 180-style geographic grid indices."""
    lon = min(179.999999, max(-180.0, longitude))
    lat = min(89.999999, max(-90.0, latitude))
    return int((lon + 180.0) // grid_degrees), int((lat + 90.0) // grid_degrees)


def grid_label(lon_index: int, lat_index: int, grid_degrees: float) -> str:
    west = -180 + lon_index * grid_degrees
    south = -90 + lat_index * grid_degrees
    return f"cell_w{west:+06.1f}_s{south:+05.1f}"


def resolve_explicit_taxon_key(
    client: httpx.Client,
    taxon_key: str,
    declared_name: str,
) -> dict[str, Any]:
    """Validate and record a caller-specified GBIF backbone taxon key.

    Explicit keys are used for broad, stable acquisition scopes when a common
    English group label does not resolve in a selected checklist. The response
    is saved verbatim enough to audit the exact GBIF concept used.
    """
    response = client.get(f"{GBIF_API}/species/{taxon_key}")
    response.raise_for_status()
    matched = response.json()
    usage_key = matched.get("key") or matched.get("nubKey") or taxon_key
    return {
        "query_name": declared_name,
        "query_rank": matched.get("rank"),
        "checklist_key": None,
        "usage_key": str(usage_key),
        "matched_name": matched.get("scientificName") or matched.get("canonicalName"),
        "matched_rank": matched.get("rank"),
        "match_type": "EXPLICIT_USAGE_KEY",
        "confidence": None,
        "resolution_method": "explicit_gbif_taxon_key",
    }


def resolve_query_taxon(
    client: httpx.Client,
    taxon_name: str,
    taxon_rank: str | None,
    checklist_key: str | None,
    taxon_key: str | None,
) -> dict[str, Any]:
    if taxon_key:
        return resolve_explicit_taxon_key(client, taxon_key, taxon_name)
    resolved = resolve_taxon(client, taxon_name, taxon_rank, checklist_key)
    resolved["resolution_method"] = "species_match"
    return resolved


def build_query_catchments(
    islands: gpd.GeoDataFrame,
    buffer_m: float,
    simplify_m: float,
) -> gpd.GeoDataFrame:
    """Build outward query envelopes while preserving original analysis geometry elsewhere.

    Buffering happens in EPSG:6933. The buffer is deliberately larger than the
    simplification tolerance so an original island edge is not intentionally
    excluded by simplification. Extra mainland records are removed later by
    joining returned points to the original exact island geometry.
    """
    if islands.crs is None:
        raise typer.BadParameter("Island GeoPackage has no CRS")
    projected = islands.to_crs(6933).copy()
    query_geometry = projected.geometry.buffer(buffer_m).simplify(
        simplify_m,
        preserve_topology=True,
    )
    query = gpd.GeoDataFrame(
        islands.drop(columns="geometry").copy(),
        geometry=gpd.GeoSeries(query_geometry, crs=6933).to_crs(4326),
        crs=4326,
    )
    query["geometry"] = query.geometry.map(canonical_gbif_polygon)
    query = explode_antimeridian_catchments(query)
    query["query_geometry_wkt"] = query.geometry.map(stable_wkt)
    query["query_geometry_sha256"] = query.geometry.map(geometry_hash)
    query["query_wkt_chars"] = query["query_geometry_wkt"].str.len()
    return query


def union_query_geometries(geometries: list[Any]) -> Polygon | MultiPolygon:
    return canonical_gbif_polygon(unary_union(geometries))


def block_query_geometries(
    catchments: gpd.GeoDataFrame,
    grid_degrees: float,
    max_wkt_chars: int,
    max_islands_per_block: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pack nearby island query envelopes into bounded WKT blocks.

    Greedy packing is deterministic: geographic grid cell, then descending WKT
    length and island ID. The actual multi-polygon union is computed at block
    flush time and its final WKT length is stored. A single complex island may
    occupy its own block; it is never silently dropped.
    """
    if max_wkt_chars < 10_000:
        raise typer.BadParameter("max_wkt_chars must be at least 10,000")
    if max_islands_per_block < 1:
        raise typer.BadParameter("max_islands_per_block must be positive")

    table = catchments.copy()
    centroids = gpd.GeoSeries(table.to_crs(6933).geometry.centroid, crs=6933).to_crs(4326)
    indices = [grid_index(point.x, point.y, grid_degrees) for point in centroids]
    table["grid_lon_index"] = [value[0] for value in indices]
    table["grid_lat_index"] = [value[1] for value in indices]
    table["regional_cell"] = [grid_label(lon, lat, grid_degrees) for lon, lat in indices]

    blocks: list[dict[str, Any]] = []
    members: list[dict[str, Any]] = []
    global_sequence = 0

    def flush(pending: list[dict[str, Any]], regional_cell: str, local_sequence: int) -> None:
        nonlocal global_sequence
        if not pending:
            return
        global_sequence += 1
        geometry = union_query_geometries([record["geometry"] for record in pending])
        query_wkt = stable_wkt(geometry)
        # block_id is derived only from (regional_cell, local_sequence), not the
        # shared global_sequence counter, so it stays stable across manifest
        # regeneration: a fix that adds/removes rows in one regional cell (e.g.
        # the antimeridian split above) must not renumber -- and thereby
        # invalidate the recorded provenance hash of -- every other cell's
        # already-submitted blocks. global_sequence/block_sequence is kept only
        # as a display-order hint.
        block_id = f"gbif_block_{regional_cell}_{local_sequence:03d}"
        blocks.append(
            {
                "island_id": block_id,
                "block_id": block_id,
                "regional_cell": regional_cell,
                "block_sequence": global_sequence,
                "n_islands": len(pending),
                "analysis_area_km2_sum": float(sum(record["area_km2"] for record in pending)),
                "query_geometry_wkt": query_wkt,
                "geometry_wkt": query_wkt,
                "geometry_wkt_sha256": sha256_text(query_wkt),
                "query_geometry_sha256": geometry_hash(geometry),
                "query_wkt_chars": len(query_wkt),
                "request_status": "prepared",
                "download_key": "",
                "error": "",
            }
        )
        for record in pending:
            members.append(
                {
                    "block_id": block_id,
                    "island_id": record["island_id"],
                    "analysis_island_id": record.get("analysis_island_id", record["island_id"]),
                    "area_km2": record["area_km2"],
                    "analysis_geometry_sha256": record["geometry_sha256"],
                    "query_geometry_sha256": record["query_geometry_sha256"],
                    "regional_cell": regional_cell,
                }
            )

    grouped = table.sort_values(
        ["regional_cell", "query_wkt_chars", "island_id"],
        ascending=[True, False, True],
    ).groupby("regional_cell", sort=True)
    for regional_cell, group in grouped:
        pending: list[dict[str, Any]] = []
        estimated_chars = 32
        local_sequence = 0
        for record in group.to_dict("records"):
            proposed = estimated_chars + int(record["query_wkt_chars"]) + 2
            if pending and (proposed > max_wkt_chars or len(pending) >= max_islands_per_block):
                local_sequence += 1
                flush(pending, str(regional_cell), local_sequence)
                pending = []
                estimated_chars = 32
            pending.append(record)
            estimated_chars += int(record["query_wkt_chars"]) + 2
        if pending:
            local_sequence += 1
            flush(pending, str(regional_cell), local_sequence)

    block_table = pd.DataFrame(blocks).sort_values("block_id").reset_index(drop=True)
    member_table = pd.DataFrame(members).sort_values(["block_id", "island_id"]).reset_index(drop=True)
    if block_table.empty or member_table.empty:
        raise RuntimeError("No GBIF request blocks were created")
    if member_table["island_id"].duplicated().any():
        raise RuntimeError("An island was assigned to more than one GBIF block")
    if len(member_table) != len(catchments):
        raise RuntimeError("Not every island was assigned to a GBIF block")
    if (block_table["query_wkt_chars"] > max_wkt_chars).any():
        raise RuntimeError(
            "At least one final GBIF query geometry exceeded max_wkt_chars; "
            "raise the configured limit or simplify the query envelopes further."
        )
    return block_table, member_table


def add_gbif_payloads(
    blocks: pd.DataFrame,
    taxon: dict[str, Any],
    checklist_key: str | None,
    email: str | None,
) -> pd.DataFrame:
    """Build coordinate-preserving SIMPLE_CSV payloads for each regional block."""
    result = blocks.copy()
    payloads = []
    for record in result.itertuples(index=False):
        payload = gbif_payload(
            wkt=record.query_geometry_wkt,
            taxon_key=taxon["usage_key"],
            checklist_key=checklist_key,
            format_name="SIMPLE_CSV",
            email=email,
            minimum_year=None,
        )
        payloads.append(json.dumps(payload, ensure_ascii=False, separators=(",", ":")))
    result["payload_json"] = payloads
    return result


@app.command("prepare")
def prepare(
    islands_gpkg: Path = typer.Option(..., exists=True, help="Prepared v2 island GeoPackage."),
    output_dir: Path = typer.Option(..., help="Directory for block manifest and members."),
    taxon_name: str = typer.Option("Tracheophyta"),
    taxon_rank: str | None = typer.Option(None),
    taxon_key: str | None = typer.Option(None, help="Explicit GBIF backbone usage key."),
    checklist_key: str | None = typer.Option(None),
    target_scope: str = typer.Option("angiosperms"),
    email: str | None = typer.Option(None, help="Optional GBIF notification email; not stored unless supplied."),
    grid_degrees: float = typer.Option(30.0, min=1.0, max=90.0),
    max_wkt_chars: int = typer.Option(120_000, min=10_000),
    max_islands_per_block: int = typer.Option(125, min=1),
    query_buffer_m: float = typer.Option(2_000.0, min=0.0),
    query_simplify_m: float = typer.Option(500.0, min=0.0),
) -> None:
    """Build bounded, coordinate-preserving GBIF regional request manifests."""
    islands = gpd.read_file(islands_gpkg, layer="islands")
    require_columns(islands, {"island_id", "area_km2", "geometry_sha256"}, "island GeoPackage")
    catchments = build_query_catchments(islands, buffer_m=query_buffer_m, simplify_m=query_simplify_m)
    blocks, members = block_query_geometries(
        catchments,
        grid_degrees=grid_degrees,
        max_wkt_chars=max_wkt_chars,
        max_islands_per_block=max_islands_per_block,
    )
    with httpx.Client(timeout=60.0, headers={"User-Agent": "island-floral-v2/0.2"}) as client:
        taxon = resolve_query_taxon(client, taxon_name, taxon_rank, checklist_key, taxon_key)
    manifest = add_gbif_payloads(blocks, taxon, checklist_key, email)

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_dir / "gbif_block_request_manifest.csv", index=False)
    members.to_csv(output_dir / "gbif_block_members.csv", index=False)
    (output_dir / "gbif_block_taxon_resolution.json").write_text(
        json.dumps(taxon, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    policy = {
        "created_at_utc": datetime.now(UTC).isoformat(),
        "format": "SIMPLE_CSV",
        "target_scope": target_scope,
        "query_taxon_name": taxon_name,
        "query_taxon_key": taxon.get("usage_key"),
        "query_taxon_rank": taxon.get("matched_rank"),
        "taxon_resolution_method": taxon.get("resolution_method"),
        "checklist_key": checklist_key,
        "n_islands": int(len(islands)),
        "n_blocks": int(len(manifest)),
        "grid_degrees": grid_degrees,
        "max_wkt_chars": max_wkt_chars,
        "max_islands_per_block": max_islands_per_block,
        "query_buffer_m": query_buffer_m,
        "query_simplify_m": query_simplify_m,
        "ring_orientation_rule": "All GBIF WKT exterior rings are counter-clockwise; holes are clockwise.",
        "assignment_rule": "Downloaded coordinates must be joined to original island geometry, not catchments.",
        "post_download_taxonomy_rule": "Retain angiosperms only after download; the broader query taxon is an acquisition proxy.",
        "email_in_payload": bool(email),
    }
    (output_dir / "gbif_block_policy.json").write_text(json.dumps(policy, indent=2), encoding="utf-8")
    typer.echo(
        f"Prepared {len(manifest)} regional GBIF occurrence-download blocks "
        f"for {len(islands)} islands at {output_dir}"
    )


if __name__ == "__main__":
    app()
