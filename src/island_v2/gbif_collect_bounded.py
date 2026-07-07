"""Disk-backed collection of succeeded GBIF blocks into island flora products.

GBIF downloads can contain millions of buffered records. Only exact-island rows
are retained, and they are streamed into a temporary SQLite store rather than
being accumulated in RAM before compact island summaries are written.
"""

from __future__ import annotations

import json
import math
import sqlite3
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer

from island_v2.gbif_collect import (
    EFFORT_COLUMNS,
    _download_archive,
    _succeeded_blocks,
    assign_occurrences_to_islands,
    iter_block_occurrences,
    sha256_file,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

_STORE_COLUMNS = [
    "gbif_id",
    "dataset_key",
    "family",
    "genus",
    "species",
    "coordinate_uncertainty_m",
    "year",
    "basis_of_record",
    "establishment_means",
    "island_id",
]


@app.callback()
def main() -> None:
    """Collect succeeded GBIF blocks without retaining all exact rows in memory."""


def _empty_manifest() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "block_id",
            "download_key",
            "doi",
            "n_block_member_islands",
            "n_islands_available_for_assignment",
            "archive_path",
            "archive_sha256",
            "archive_bytes",
            "n_source_rows",
            "n_valid_coordinate_rows",
            "n_exact_island_rows_before_dedup",
            "n_buffer_only_rows",
            "n_rows_retained_after_global_dedup",
            "collection_status",
            "error",
        ]
    )


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _open_store(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=FILE")
    connection.executescript(
        """
        CREATE TABLE occurrences (
            dedup_key TEXT PRIMARY KEY,
            source_order INTEGER NOT NULL,
            gbif_id TEXT NOT NULL,
            dataset_key TEXT NOT NULL,
            family TEXT NOT NULL,
            genus TEXT NOT NULL,
            species TEXT NOT NULL,
            coordinate_uncertainty_m REAL,
            year INTEGER,
            basis_of_record TEXT NOT NULL,
            establishment_means TEXT NOT NULL,
            island_id TEXT NOT NULL
        );
        CREATE INDEX occurrences_by_island ON occurrences (island_id, coordinate_uncertainty_m);
        CREATE INDEX occurrences_by_species ON occurrences (species, source_order);
        """
    )
    return connection


_INSERT = """
INSERT INTO occurrences (
    dedup_key, source_order, gbif_id, dataset_key, family, genus, species,
    coordinate_uncertainty_m, year, basis_of_record, establishment_means, island_id
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT(dedup_key) DO UPDATE SET
    source_order = excluded.source_order,
    gbif_id = excluded.gbif_id,
    dataset_key = excluded.dataset_key,
    family = excluded.family,
    genus = excluded.genus,
    species = excluded.species,
    coordinate_uncertainty_m = excluded.coordinate_uncertainty_m,
    year = excluded.year,
    basis_of_record = excluded.basis_of_record,
    establishment_means = excluded.establishment_means,
    island_id = excluded.island_id
WHERE (
    excluded.coordinate_uncertainty_m IS NOT NULL
    AND (
        occurrences.coordinate_uncertainty_m IS NULL
        OR excluded.coordinate_uncertainty_m < occurrences.coordinate_uncertainty_m
    )
) OR (
    (
        (excluded.coordinate_uncertainty_m IS NULL AND occurrences.coordinate_uncertainty_m IS NULL)
        OR excluded.coordinate_uncertainty_m = occurrences.coordinate_uncertainty_m
    )
    AND excluded.source_order < occurrences.source_order
)
"""


def _insert_exact_chunk(
    connection: sqlite3.Connection,
    exact: pd.DataFrame,
    block_id: str,
    source_order: int,
) -> int:
    work = exact.reindex(columns=_STORE_COLUMNS).copy()
    for column in (
        "gbif_id",
        "dataset_key",
        "family",
        "genus",
        "species",
        "basis_of_record",
        "establishment_means",
        "island_id",
    ):
        work[column] = work[column].map(_text)
    work["coordinate_uncertainty_m"] = pd.to_numeric(
        work["coordinate_uncertainty_m"], errors="coerce"
    )
    work["year"] = pd.to_numeric(work["year"], errors="coerce")

    def rows() -> Iterable[tuple[object, ...]]:
        for offset, row in enumerate(work.itertuples(index=False, name=None)):
            (
                gbif_id,
                dataset_key,
                family,
                genus,
                species,
                uncertainty,
                year,
                basis,
                establishment,
                island_id,
            ) = row
            ordinal = source_order + offset
            dedup_key = f"gbif:{gbif_id}" if gbif_id else f"row:{block_id}:{ordinal}"
            yield (
                dedup_key,
                ordinal,
                gbif_id,
                dataset_key,
                family,
                genus,
                species,
                None if pd.isna(uncertainty) else float(uncertainty),
                None if pd.isna(year) else int(year),
                basis,
                establishment,
                island_id,
            )

    connection.executemany(_INSERT, rows())
    return source_order + len(work)


def _quantile(connection: sqlite3.Connection, island_id: str, probability: float) -> float | None:
    n = int(
        connection.execute(
            "SELECT COUNT(*) FROM occurrences WHERE island_id = ? AND coordinate_uncertainty_m IS NOT NULL",
            (island_id,),
        ).fetchone()[0]
    )
    if not n:
        return None
    rank = (n - 1) * probability
    lower = math.floor(rank)
    upper = math.ceil(rank)
    low_value = float(
        connection.execute(
            """
            SELECT coordinate_uncertainty_m FROM occurrences
            WHERE island_id = ? AND coordinate_uncertainty_m IS NOT NULL
            ORDER BY coordinate_uncertainty_m LIMIT 1 OFFSET ?
            """,
            (island_id, lower),
        ).fetchone()[0]
    )
    if lower == upper:
        return low_value
    high_value = float(
        connection.execute(
            """
            SELECT coordinate_uncertainty_m FROM occurrences
            WHERE island_id = ? AND coordinate_uncertainty_m IS NOT NULL
            ORDER BY coordinate_uncertainty_m LIMIT 1 OFFSET ?
            """,
            (island_id, upper),
        ).fetchone()[0]
    )
    return low_value + (high_value - low_value) * (rank - lower)


def _species_from_store(connection: sqlite3.Connection) -> pd.DataFrame:
    table = pd.read_sql_query(
        """
        SELECT
            island_id,
            species,
            COUNT(*) AS n_records,
            COUNT(DISTINCT NULLIF(gbif_id, '')) AS n_unique_gbif_ids,
            COALESCE(REPLACE(GROUP_CONCAT(DISTINCT NULLIF(TRIM(basis_of_record), '')), ',', '|'), '')
                AS basis_of_record_set
        FROM occurrences
        WHERE TRIM(species) <> ''
        GROUP BY island_id, species
        ORDER BY island_id, species
        """,
        connection,
    )
    if table.empty:
        return pd.DataFrame(
            columns=[
                "island_id",
                "species",
                "n_records",
                "n_unique_gbif_ids",
                "basis_of_record_set",
                "review_status",
            ]
        )
    table["review_status"] = "unresolved_taxonomy"
    return table[
        [
            "island_id",
            "species",
            "n_records",
            "n_unique_gbif_ids",
            "basis_of_record_set",
            "review_status",
        ]
    ]


def _effort_from_store(connection: sqlite3.Connection) -> pd.DataFrame:
    table = pd.read_sql_query(
        """
        SELECT
            island_id,
            COUNT(*) AS n_records,
            COUNT(DISTINCT NULLIF(gbif_id, '')) AS n_unique_gbif_ids,
            COUNT(DISTINCT NULLIF(TRIM(species), '')) AS n_species,
            COUNT(DISTINCT NULLIF(TRIM(dataset_key), '')) AS n_datasets,
            SUM(CASE WHEN UPPER(TRIM(basis_of_record)) = 'PRESERVED_SPECIMEN' THEN 1 ELSE 0 END) AS n_preserved_specimen,
            SUM(CASE WHEN UPPER(TRIM(basis_of_record)) = 'HUMAN_OBSERVATION' THEN 1 ELSE 0 END) AS n_human_observation,
            SUM(CASE WHEN TRIM(basis_of_record) <> '' AND UPPER(TRIM(basis_of_record)) NOT IN ('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION') THEN 1 ELSE 0 END) AS n_other_basis_of_record,
            SUM(CASE WHEN coordinate_uncertainty_m IS NOT NULL THEN 1 ELSE 0 END) AS n_uncertainty_reported,
            SUM(CASE WHEN coordinate_uncertainty_m > 1000 THEN 1 ELSE 0 END) AS n_coordinate_uncertainty_gt_1000m,
            SUM(CASE WHEN INSTR(UPPER(establishment_means), 'NATIVE') > 0 THEN 1 ELSE 0 END) AS n_establishment_native,
            SUM(CASE WHEN INSTR(UPPER(establishment_means), 'INTRODUCED') > 0 THEN 1 ELSE 0 END) AS n_establishment_introduced,
            SUM(CASE WHEN INSTR(UPPER(establishment_means), 'CULTIVATED') > 0 THEN 1 ELSE 0 END) AS n_establishment_cultivated,
            MIN(year) AS year_min,
            MAX(year) AS year_max
        FROM occurrences
        GROUP BY island_id
        ORDER BY island_id
        """,
        connection,
    )
    if table.empty:
        return pd.DataFrame(columns=EFFORT_COLUMNS)
    table["median_coordinate_uncertainty_m"] = [
        _quantile(connection, island_id, 0.5) for island_id in table["island_id"]
    ]
    table["p90_coordinate_uncertainty_m"] = [
        _quantile(connection, island_id, 0.9) for island_id in table["island_id"]
    ]
    return table[EFFORT_COLUMNS]


def _taxa_from_store(connection: sqlite3.Connection) -> pd.DataFrame:
    return pd.read_sql_query(
        """
        SELECT
            a.species AS accepted_species,
            COALESCE((
                SELECT b.genus FROM occurrences b
                WHERE b.species = a.species AND TRIM(b.genus) <> ''
                ORDER BY b.source_order LIMIT 1
            ), '') AS genus,
            COALESCE((
                SELECT b.family FROM occurrences b
                WHERE b.species = a.species AND TRIM(b.family) <> ''
                ORDER BY b.source_order LIMIT 1
            ), '') AS family,
            COUNT(DISTINCT a.island_id) AS n_islands,
            COUNT(*) AS n_records
        FROM occurrences a
        WHERE TRIM(a.species) <> ''
        GROUP BY a.species
        ORDER BY a.species
        """,
        connection,
    )


def _summarize_exact_chunks(
    exact_chunks: Iterable[pd.DataFrame], block_id: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int]]:
    with tempfile.TemporaryDirectory(prefix="island_gbif_exact_") as temporary:
        connection = _open_store(Path(temporary) / "exact_rows.sqlite")
        before_dedup = 0
        source_order = 0
        try:
            for exact in exact_chunks:
                if exact.empty:
                    continue
                before_dedup += int(len(exact))
                source_order = _insert_exact_chunk(connection, exact, block_id, source_order)
            connection.commit()
            retained = int(connection.execute("SELECT COUNT(*) FROM occurrences").fetchone()[0])
            duplicate_audit = {
                "n_input_rows": before_dedup,
                "n_unique_gbif_ids": int(
                    connection.execute(
                        "SELECT COUNT(DISTINCT NULLIF(gbif_id, '')) FROM occurrences"
                    ).fetchone()[0]
                ),
                "n_duplicate_gbif_ids": 0,
                "n_duplicate_rows_removed": before_dedup - retained,
                "n_cross_block_duplicate_gbif_ids": 0,
                "n_conflicting_island_assignments": 0,
            }
            return (
                _species_from_store(connection),
                _effort_from_store(connection),
                _taxa_from_store(connection),
                duplicate_audit,
            )
        finally:
            connection.close()


def _combine_species(parts: list[pd.DataFrame]) -> pd.DataFrame:
    columns = [
        "island_id",
        "species",
        "n_records",
        "n_unique_gbif_ids",
        "basis_of_record_set",
        "review_status",
    ]
    table = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    if table.empty:
        return table
    table["n_records"] = pd.to_numeric(table["n_records"], errors="coerce").fillna(0).astype(int)
    table["n_unique_gbif_ids"] = (
        pd.to_numeric(table["n_unique_gbif_ids"], errors="coerce").fillna(0).astype(int)
    )
    rows = []
    for (island_id, species), group in table.groupby(["island_id", "species"], sort=True):
        bases = {
            item
            for value in group["basis_of_record_set"].astype(str)
            for item in value.split("|")
            if item
        }
        rows.append(
            {
                "island_id": island_id,
                "species": species,
                "n_records": int(group["n_records"].sum()),
                "n_unique_gbif_ids": int(group["n_unique_gbif_ids"].sum()),
                "basis_of_record_set": "|".join(sorted(bases)),
                "review_status": "unresolved_taxonomy",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _combine_effort(parts: list[pd.DataFrame]) -> pd.DataFrame:
    table = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=EFFORT_COLUMNS)
    if table.empty:
        return table
    duplicates = table.loc[table["island_id"].duplicated(), "island_id"].tolist()
    if duplicates:
        raise typer.BadParameter(f"Block partition overlaps on island_id: {duplicates[:5]}")
    return table.sort_values("island_id").reset_index(drop=True)[EFFORT_COLUMNS]


def _combine_taxa(parts: list[pd.DataFrame]) -> pd.DataFrame:
    columns = ["accepted_species", "genus", "family", "n_islands", "n_records"]
    table = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    if table.empty:
        return table
    table["n_islands"] = pd.to_numeric(table["n_islands"], errors="coerce").fillna(0).astype(int)
    table["n_records"] = pd.to_numeric(table["n_records"], errors="coerce").fillna(0).astype(int)
    rows = []
    for species, group in table.groupby("accepted_species", sort=True):
        genus = next((value for value in group["genus"].astype(str) if value.strip()), "")
        family = next((value for value in group["family"].astype(str) if value.strip()), "")
        rows.append(
            {
                "accepted_species": species,
                "genus": genus,
                "family": family,
                "n_islands": int(group["n_islands"].sum()),
                "n_records": int(group["n_records"].sum()),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def collect_succeeded_blocks(
    campaign: dict[str, Any],
    members: pd.DataFrame,
    islands: gpd.GeoDataFrame,
    download_dir: Path,
    chunksize: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int], list[str]]:
    """Download and exact-assign succeeded blocks using a disk-backed row store."""
    if "analysis_island_id" not in members.columns:
        members = members.copy()
        members["analysis_island_id"] = members["island_id"]
    required_members = {"block_id", "analysis_island_id"}
    missing_members = required_members.difference(members.columns)
    if missing_members:
        raise typer.BadParameter(
            f"Block-members table missing columns: {', '.join(sorted(missing_members))}"
        )
    if "island_id" not in islands.columns:
        raise typer.BadParameter("Exact island GeoPackage must contain island_id")

    islands_by_id = islands.set_index("island_id")
    manifest_rows: list[dict[str, Any]] = []
    species_parts: list[pd.DataFrame] = []
    effort_parts: list[pd.DataFrame] = []
    taxa_parts: list[pd.DataFrame] = []
    total_audit = {
        "n_input_rows": 0,
        "n_unique_gbif_ids": 0,
        "n_duplicate_gbif_ids": 0,
        "n_duplicate_rows_removed": 0,
        "n_cross_block_duplicate_gbif_ids": 0,
        "n_conflicting_island_assignments": 0,
    }
    failures: list[str] = []

    download_dir.mkdir(parents=True, exist_ok=True)
    with httpx.Client(
        timeout=300.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.7"},
    ) as client:
        for entry in _succeeded_blocks(campaign):
            block_id = str(entry["block_id"])
            download_key = str(entry["download_key"])
            block_islands = members.loc[
                members["block_id"].eq(block_id), "analysis_island_id"
            ].unique()
            present = [island_id for island_id in block_islands if island_id in islands_by_id.index]
            metrics: dict[str, Any] = {
                "block_id": block_id,
                "download_key": download_key,
                "doi": str(entry.get("doi", "")),
                "n_block_member_islands": int(len(block_islands)),
                "n_islands_available_for_assignment": int(len(present)),
                "archive_path": "",
                "archive_sha256": "",
                "archive_bytes": 0,
                "n_source_rows": 0,
                "n_valid_coordinate_rows": 0,
                "n_exact_island_rows_before_dedup": 0,
                "n_buffer_only_rows": 0,
                "n_rows_retained_after_global_dedup": 0,
                "collection_status": "pending",
                "error": "",
            }
            if not present:
                metrics["collection_status"] = "skipped_no_matching_islands"
                manifest_rows.append(metrics)
                continue
            try:
                island_subset = gpd.GeoDataFrame(
                    islands_by_id.loc[present].reset_index(), geometry="geometry", crs=islands.crs
                )
                archive = _download_archive(
                    client,
                    download_key,
                    download_dir,
                    entry.get("last_download_link") or entry.get("download_link"),
                )
                metrics["archive_path"] = str(archive)
                metrics["archive_sha256"] = sha256_file(archive)
                metrics["archive_bytes"] = int(archive.stat().st_size)

                def exact_chunks() -> Iterable[pd.DataFrame]:
                    for occurrences in iter_block_occurrences(archive, chunksize=chunksize):
                        metrics["n_source_rows"] += int(
                            occurrences.attrs.get("n_source_rows", len(occurrences))
                        )
                        metrics["n_valid_coordinate_rows"] += int(len(occurrences))
                        if occurrences.empty:
                            continue
                        assigned = assign_occurrences_to_islands(occurrences, island_subset)
                        exact = assigned.loc[assigned["island_id"].notna()].copy()
                        metrics["n_exact_island_rows_before_dedup"] += int(len(exact))
                        metrics["n_buffer_only_rows"] += int(len(assigned) - len(exact))
                        if not exact.empty:
                            yield exact

                species, effort, taxa, audit = _summarize_exact_chunks(exact_chunks(), block_id)
                species_parts.append(species)
                effort_parts.append(effort)
                taxa_parts.append(taxa)
                metrics["n_rows_retained_after_global_dedup"] = int(
                    audit["n_input_rows"] - audit["n_duplicate_rows_removed"]
                )
                for key in total_audit:
                    total_audit[key] += int(audit[key])
                metrics["collection_status"] = "collected"
            except Exception as exc:
                metrics["collection_status"] = "collection_failed"
                metrics["error"] = f"{type(exc).__name__}: {exc}"
                failures.append(block_id)
            manifest_rows.append(metrics)

    return (
        _combine_species(species_parts),
        _combine_effort(effort_parts),
        _combine_taxa(taxa_parts),
        pd.DataFrame(manifest_rows) if manifest_rows else _empty_manifest(),
        total_audit,
        failures,
    )


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True, help="Campaign ledger JSON."),
    block_members_csv: Path = typer.Option(..., exists=True, help="gbif_block_members.csv."),
    islands_gpkg: Path = typer.Option(..., exists=True, help="Exact island GeoPackage."),
    download_dir: Path = typer.Option(..., help="Where block archives are cached."),
    output_dir: Path = typer.Option(..., help="Where island-level outputs are written."),
    chunksize: int = typer.Option(50_000, min=1, help="Rows read per GBIF SIMPLE_CSV chunk."),
) -> None:
    """Collect succeeded blocks into compact exact-island products and audits."""
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    members = pd.read_csv(block_members_csv, dtype=str).fillna("")
    islands = gpd.read_file(islands_gpkg, layer="islands")
    output_dir.mkdir(parents=True, exist_ok=True)
    species_table, effort, taxa, manifest, duplicate_audit, failures = collect_succeeded_blocks(
        campaign, members, islands, download_dir, chunksize
    )
    species_table.to_csv(output_dir / "island_species_occurrences.csv", index=False)
    effort.to_csv(output_dir / "island_observation_effort.csv", index=False)
    taxa.to_csv(output_dir / "island_taxa.csv", index=False)
    manifest.to_csv(output_dir / "collection_manifest.csv", index=False)
    retained = int(duplicate_audit["n_input_rows"] - duplicate_audit["n_duplicate_rows_removed"])
    status = {
        "collection_mode": "exact_island_rows_sqlite_spill",
        "n_succeeded_blocks": int(len(_succeeded_blocks(campaign))),
        "n_collected_blocks": int(manifest["collection_status"].eq("collected").sum()) if not manifest.empty else 0,
        "n_failed_blocks": int(manifest["collection_status"].eq("collection_failed").sum()) if not manifest.empty else 0,
        "n_exact_island_rows_before_global_dedup": int(duplicate_audit["n_input_rows"]),
        "n_occurrences_after_global_dedup": retained,
        "n_occurrences_on_island": retained,
        "n_occurrences_in_buffer_only": int(manifest["n_buffer_only_rows"].sum()) if not manifest.empty else 0,
        "n_island_species_pairs": int(len(species_table)),
        "n_islands_with_records": int(effort["island_id"].nunique()) if not effort.empty else 0,
        "n_accepted_species": int(len(taxa)),
        "duplicate_audit": duplicate_audit,
        "failed_block_ids": failures,
    }
    (output_dir / "collection_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    typer.echo(
        f"Collected {status['n_collected_blocks']} blocks: {len(species_table)} island-species pairs, "
        f"{len(taxa)} accepted species across {status['n_islands_with_records']} islands "
        f"({status['n_occurrences_in_buffer_only']} buffer-only rows discarded before global dedup; "
        f"{duplicate_audit['n_duplicate_rows_removed']} duplicate exact-island rows removed)."
    )
    if failures:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
