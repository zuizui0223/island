"""Sharded full taxonomy normalization for Island Plant Database 2.0 alpha2."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.island_plant_database_taxonomy import (
    ALPHA1_TAXA_SUMMARY_SHA256,
    COL_XR_CHECKLIST_KEY,
    SOURCE_LICENSE,
    load_source_taxa,
    match_rows,
    sha256_file,
    summarize,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

SOURCE_TOTAL = 115_328
DEFAULT_SEGMENT_SIZE = 10_000
DEFAULT_SEGMENT_COUNT = 12


@app.callback()
def main() -> None:
    """Run and combine deterministic alpha2 taxonomy segments."""


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_json(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return _sha256_bytes(payload)


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            frame.to_csv(zipped, index=False)


def segment_bounds(
    segment_index: int,
    total: int = SOURCE_TOTAL,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
) -> tuple[int, int]:
    if segment_index < 0:
        raise ValueError("segment_index must be non-negative")
    if total < 1 or segment_size < 1:
        raise ValueError("total and segment_size must be positive")
    start = segment_index * segment_size
    end = min(total, start + segment_size)
    if start >= total:
        raise ValueError(f"segment_index {segment_index} starts beyond source total {total}")
    return start, end


def segment_filename(index: int, stem: str, suffix: str) -> str:
    return f"{stem}_segment_{index:02d}{suffix}"


def build_segment(
    source_taxa: Path,
    output_dir: Path,
    segment_index: int,
    *,
    segment_size: int = DEFAULT_SEGMENT_SIZE,
    batch_size: int = 1000,
    expected_source_sha256: str = ALPHA1_TAXA_SUMMARY_SHA256,
) -> dict[str, Any]:
    actual_sha = sha256_file(source_taxa)
    if actual_sha != expected_source_sha256:
        raise ValueError(
            f"taxonomy source SHA mismatch: expected {expected_source_sha256}, got {actual_sha}"
        )
    source = load_source_taxa(source_taxa)
    if len(source) != SOURCE_TOTAL:
        raise ValueError(f"taxonomy source row count changed: {len(source)} != {SOURCE_TOTAL}")
    start, end = segment_bounds(segment_index, len(source), segment_size)
    segment = source.iloc[start:end].reset_index(drop=True)
    crosswalk, raw_batches, metadata = match_rows(segment, batch_size=batch_size)
    if len(crosswalk) != len(segment):
        raise RuntimeError("segment crosswalk cardinality mismatch")
    if crosswalk["submitted_name"].duplicated().any():
        raise RuntimeError("segment crosswalk contains duplicate submitted names")

    output_dir.mkdir(parents=True, exist_ok=True)
    crosswalk_path = output_dir / segment_filename(segment_index, "taxonomy_crosswalk", ".csv.gz")
    metadata_path = output_dir / segment_filename(segment_index, "gbif_colxr_metadata", ".json")
    manifest_path = output_dir / segment_filename(segment_index, "segment_manifest", ".json")
    _write_gzip_csv(crosswalk, crosswalk_path)
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    raw_dir = output_dir / f"raw_segment_{segment_index:02d}"
    raw_dir.mkdir(exist_ok=True)
    for local_index, payload in enumerate(raw_batches):
        global_start = start + int(payload["start"])
        raw_payload = {
            "segment_index": segment_index,
            "global_start": global_start,
            "queries": payload["queries"],
            "results": payload["results"],
        }
        (raw_dir / f"batch_{global_start:06d}.json").write_text(
            json.dumps(raw_payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

    summary = summarize(crosswalk)
    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "target_version": "2.0.0-alpha2-taxonomy",
        "stage": "colxr_taxonomy_full_segment",
        "segment_index": segment_index,
        "segment_start": start,
        "segment_end_exclusive": end,
        "segment_rows": int(len(segment)),
        "segment_size_parameter": segment_size,
        "batch_size": batch_size,
        "source_database_version": "2.0.0-alpha1",
        "source_taxa_total": int(len(source)),
        "source_taxa_summary_sha256": actual_sha,
        "checklist_key": COL_XR_CHECKLIST_KEY,
        "source_license": SOURCE_LICENSE,
        "metadata_sha256": _sha256_json(metadata),
        "crosswalk_sha256": sha256_file(crosswalk_path),
        "summary": summary,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def _segment_manifests(segments_dir: Path) -> list[tuple[Path, dict[str, Any]]]:
    rows = []
    for path in sorted(segments_dir.glob("segment_manifest_segment_*.json")):
        rows.append((path, json.loads(path.read_text(encoding="utf-8"))))
    return rows


def combine_segments(
    segments_dir: Path,
    output_dir: Path,
    *,
    expected_segments: int = DEFAULT_SEGMENT_COUNT,
    expected_total: int = SOURCE_TOTAL,
) -> dict[str, Any]:
    manifests = _segment_manifests(segments_dir)
    if len(manifests) != expected_segments:
        raise ValueError(f"expected {expected_segments} segment manifests, found {len(manifests)}")

    indices = [int(item[1]["segment_index"]) for item in manifests]
    if indices != list(range(expected_segments)):
        raise ValueError(f"segment indices are not complete and ordered: {indices}")
    source_hashes = {item[1]["source_taxa_summary_sha256"] for item in manifests}
    metadata_hashes = {item[1]["metadata_sha256"] for item in manifests}
    checklist_keys = {item[1]["checklist_key"] for item in manifests}
    if source_hashes != {ALPHA1_TAXA_SUMMARY_SHA256}:
        raise ValueError(f"segment source hashes disagree: {source_hashes}")
    if len(metadata_hashes) != 1:
        raise ValueError(f"COL XR metadata changed across segments: {metadata_hashes}")
    if checklist_keys != {COL_XR_CHECKLIST_KEY}:
        raise ValueError(f"segment checklist keys disagree: {checklist_keys}")

    expected_start = 0
    parts: list[pd.DataFrame] = []
    for _, manifest in manifests:
        index = int(manifest["segment_index"])
        start = int(manifest["segment_start"])
        end = int(manifest["segment_end_exclusive"])
        if start != expected_start:
            raise ValueError(f"segment {index} starts at {start}, expected {expected_start}")
        if end <= start:
            raise ValueError(f"segment {index} has non-positive extent")
        expected_start = end
        crosswalk_path = segments_dir / segment_filename(index, "taxonomy_crosswalk", ".csv.gz")
        if not crosswalk_path.exists():
            raise FileNotFoundError(crosswalk_path)
        if sha256_file(crosswalk_path) != manifest["crosswalk_sha256"]:
            raise ValueError(f"segment {index} crosswalk SHA mismatch")
        frame = pd.read_csv(crosswalk_path, dtype=str).fillna("")
        if len(frame) != int(manifest["segment_rows"]):
            raise ValueError(f"segment {index} row count mismatch")
        parts.append(frame)
    if expected_start != expected_total:
        raise ValueError(f"segments cover {expected_start} rows, expected {expected_total}")

    full = pd.concat(parts, ignore_index=True)
    if len(full) != expected_total:
        raise ValueError(f"combined crosswalk has {len(full)} rows, expected {expected_total}")
    if full["submitted_name"].duplicated().any():
        dupes = full.loc[full["submitted_name"].duplicated(), "submitted_name"].head(10).tolist()
        raise ValueError(f"combined crosswalk has duplicate submitted names: {dupes}")
    if full["checklist_key"].ne(COL_XR_CHECKLIST_KEY).any():
        raise ValueError("combined crosswalk contains an unexpected checklist key")

    review = full.loc[~full["automatic_resolution_candidate"].astype(str).eq("true")].copy()
    output_dir.mkdir(parents=True, exist_ok=True)
    crosswalk_path = output_dir / "taxonomy_crosswalk.csv.gz"
    review_path = output_dir / "taxonomy_review_queue.csv.gz"
    _write_gzip_csv(full, crosswalk_path)
    _write_gzip_csv(review, review_path)

    first_metadata_path = segments_dir / segment_filename(0, "gbif_colxr_metadata", ".json")
    metadata = json.loads(first_metadata_path.read_text(encoding="utf-8"))
    metadata_path = output_dir / "gbif_colxr_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary = summarize(full)
    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha2-taxonomy",
        "stage": "colxr_taxonomy_full_crosswalk",
        "source_database_version": "2.0.0-alpha1",
        "source_taxa_total": expected_total,
        "source_taxa_summary_sha256": ALPHA1_TAXA_SUMMARY_SHA256,
        "checklist_key": COL_XR_CHECKLIST_KEY,
        "taxonomy": "Catalogue of Life Extended Release",
        "source_license": SOURCE_LICENSE,
        "segment_count": expected_segments,
        "metadata_sha256": next(iter(metadata_hashes)),
        "crosswalk_sha256": sha256_file(crosswalk_path),
        "review_queue_sha256": sha256_file(review_path),
        "summary": summary,
        "scientific_boundary": {
            "alpha1_mutated": False,
            "automatic_results_are_candidates": True,
            "review_required_results_promoted": False,
            "unmatched_results_promoted": False,
        },
    }
    (output_dir / "TAXONOMY_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


@app.command("segment")
def segment_command(
    source_taxa: Path = typer.Option(..., "--source-taxa", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
    segment_index: int = typer.Option(..., "--segment-index", min=0),
    segment_size: int = typer.Option(DEFAULT_SEGMENT_SIZE, "--segment-size", min=1),
    batch_size: int = typer.Option(1000, "--batch-size", min=1, max=1000),
    expected_source_sha256: str = typer.Option(ALPHA1_TAXA_SUMMARY_SHA256, "--expected-source-sha256"),
) -> None:
    manifest = build_segment(
        source_taxa,
        output_dir,
        segment_index,
        segment_size=segment_size,
        batch_size=batch_size,
        expected_source_sha256=expected_source_sha256,
    )
    typer.echo(json.dumps(manifest["summary"], sort_keys=True))


@app.command("combine")
def combine_command(
    segments_dir: Path = typer.Option(..., "--segments-dir", exists=True, file_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
    expected_segments: int = typer.Option(DEFAULT_SEGMENT_COUNT, "--expected-segments", min=1),
    expected_total: int = typer.Option(SOURCE_TOTAL, "--expected-total", min=1),
) -> None:
    manifest = combine_segments(
        segments_dir,
        output_dir,
        expected_segments=expected_segments,
        expected_total=expected_total,
    )
    typer.echo(json.dumps(manifest["summary"], sort_keys=True))


if __name__ == "__main__":
    app()
