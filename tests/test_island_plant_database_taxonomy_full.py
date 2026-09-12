from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.island_plant_database_taxonomy import (
    COL_XR_CHECKLIST_KEY,
    OUTPUT_COLUMNS,
    summarize,
)
from island_v2.island_plant_database_taxonomy_full import (
    _sha256_json,
    _write_gzip_csv,
    combine_segments,
    segment_bounds,
    segment_filename,
)


def _crosswalk(names: list[str], automatic: list[bool]) -> pd.DataFrame:
    rows = []
    for name, is_auto in zip(names, automatic, strict=True):
        row = {column: "" for column in OUTPUT_COLUMNS}
        row.update(
            {
                "submitted_name": name,
                "submitted_genus": name.split()[0],
                "submitted_family": "Testaceae",
                "matched_usage_key": f"usage-{name}",
                "matched_usage_name": name,
                "matched_canonical_name": name,
                "matched_status": "ACCEPTED",
                "matched_rank": "SPECIES",
                "synonym": "false",
                "candidate_accepted_key": f"accepted-{name}",
                "candidate_accepted_name": name,
                "candidate_accepted_rank": "SPECIES",
                "accepted_target_basis": "usage",
                "candidate_genus": name.split()[0],
                "candidate_family": "Testaceae",
                "candidate_kingdom": "Plantae",
                "match_type": "EXACT" if is_auto else "VARIANT",
                "confidence": "100" if is_auto else "90",
                "automatic_resolution_candidate": "true" if is_auto else "false",
                "resolution_status": "exact_accepted_candidate" if is_auto else "review_required",
                "review_status": "automatic_candidate_not_promoted" if is_auto else "needs_taxonomy_review",
                "checklist_key": COL_XR_CHECKLIST_KEY,
                "source_license": "CC-BY-4.0",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def _write_segment(
    root: Path,
    index: int,
    start: int,
    names: list[str],
    automatic: list[bool],
    *,
    source_hash: str,
    metadata: dict[str, object],
) -> None:
    frame = _crosswalk(names, automatic)
    crosswalk_path = root / segment_filename(index, "taxonomy_crosswalk", ".csv.gz")
    _write_gzip_csv(frame, crosswalk_path)
    metadata_path = root / segment_filename(index, "gbif_colxr_metadata", ".json")
    metadata_path.write_text(json.dumps(metadata, sort_keys=True), encoding="utf-8")
    from island_v2.island_plant_database_taxonomy import sha256_file

    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "target_version": "2.0.0-alpha2-taxonomy",
        "stage": "colxr_taxonomy_full_segment",
        "segment_index": index,
        "segment_start": start,
        "segment_end_exclusive": start + len(names),
        "segment_rows": len(names),
        "segment_size_parameter": len(names),
        "batch_size": 1000,
        "source_database_version": "2.0.0-alpha1",
        "source_taxa_total": 4,
        "source_taxa_summary_sha256": source_hash,
        "checklist_key": COL_XR_CHECKLIST_KEY,
        "source_license": "CC-BY-4.0",
        "metadata_sha256": _sha256_json(metadata),
        "crosswalk_sha256": sha256_file(crosswalk_path),
        "summary": summarize(frame),
    }
    (root / segment_filename(index, "segment_manifest", ".json")).write_text(
        json.dumps(manifest, sort_keys=True), encoding="utf-8"
    )


def test_segment_bounds_cover_last_partial_segment() -> None:
    assert segment_bounds(0, total=25, segment_size=10) == (0, 10)
    assert segment_bounds(1, total=25, segment_size=10) == (10, 20)
    assert segment_bounds(2, total=25, segment_size=10) == (20, 25)
    with pytest.raises(ValueError, match="starts beyond"):
        segment_bounds(3, total=25, segment_size=10)


def test_combine_segments_preserves_complete_unique_coverage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    segments = tmp_path / "segments"
    output = tmp_path / "combined"
    segments.mkdir()
    source_hash = "a" * 64
    metadata = {"checklistKey": COL_XR_CHECKLIST_KEY, "release": "test"}

    _write_segment(segments, 0, 0, ["Alpha one", "Beta two"], [True, False], source_hash=source_hash, metadata=metadata)
    _write_segment(segments, 1, 2, ["Gamma three", "Delta four"], [True, True], source_hash=source_hash, metadata=metadata)

    import island_v2.island_plant_database_taxonomy_full as full

    monkeypatch.setattr(full, "ALPHA1_TAXA_SUMMARY_SHA256", source_hash)
    manifest = combine_segments(segments, output, expected_segments=2, expected_total=4)

    crosswalk = pd.read_csv(output / "taxonomy_crosswalk.csv.gz", dtype=str).fillna("")
    review = pd.read_csv(output / "taxonomy_review_queue.csv.gz", dtype=str).fillna("")
    assert crosswalk["submitted_name"].tolist() == ["Alpha one", "Beta two", "Gamma three", "Delta four"]
    assert review["submitted_name"].tolist() == ["Beta two"]
    assert manifest["summary"]["n_input"] == 4
    assert manifest["summary"]["n_automatic_resolution_candidates"] == 3
    assert manifest["summary"]["n_review_required"] == 1
    assert manifest["scientific_boundary"]["alpha1_mutated"] is False


def test_combine_rejects_gap_or_overlap(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    segments = tmp_path / "segments"
    output = tmp_path / "combined"
    segments.mkdir()
    source_hash = "b" * 64
    metadata = {"checklistKey": COL_XR_CHECKLIST_KEY, "release": "test"}

    _write_segment(segments, 0, 0, ["Alpha one", "Beta two"], [True, True], source_hash=source_hash, metadata=metadata)
    _write_segment(segments, 1, 3, ["Gamma three"], [True], source_hash=source_hash, metadata=metadata)

    import island_v2.island_plant_database_taxonomy_full as full

    monkeypatch.setattr(full, "ALPHA1_TAXA_SUMMARY_SHA256", source_hash)
    with pytest.raises(ValueError, match="starts at 3, expected 2"):
        combine_segments(segments, output, expected_segments=2, expected_total=3)


def test_combine_rejects_taxonomy_metadata_drift(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    segments = tmp_path / "segments"
    output = tmp_path / "combined"
    segments.mkdir()
    source_hash = "c" * 64

    _write_segment(
        segments, 0, 0, ["Alpha one"], [True], source_hash=source_hash,
        metadata={"checklistKey": COL_XR_CHECKLIST_KEY, "release": "one"}
    )
    _write_segment(
        segments, 1, 1, ["Beta two"], [True], source_hash=source_hash,
        metadata={"checklistKey": COL_XR_CHECKLIST_KEY, "release": "two"}
    )

    import island_v2.island_plant_database_taxonomy_full as full

    monkeypatch.setattr(full, "ALPHA1_TAXA_SUMMARY_SHA256", source_hash)
    with pytest.raises(ValueError, match="metadata changed across segments"):
        combine_segments(segments, output, expected_segments=2, expected_total=2)
