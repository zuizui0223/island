"""Independent WFO 2026-06 audit for alpha3 taxonomy review candidates.

The audit is exact-name only and fail-closed. It never mutates alpha1/alpha2,
never applies fuzzy matching, and never promotes a taxon. It identifies exact
WFO species concepts that can be reviewed as independent rescue candidates.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

WFO_RELEASE = "2026-06"
WFO_ZENODO_RECORD = "20782718"
WFO_DOI = "10.5281/zenodo.20782718"
WFO_EXPECTED_MD5 = "0e4486945cd9f7af548ca87eb9a870ed"
WFO_MEMBER = "classification.csv"
WFO_LICENSE = "CC0-1.0"
EXPECTED_QUEUE_ROWS = 3_436
INFRA_RANKS = {"SUBSPECIES", "VARIETY", "FORM"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _key(value: object) -> str:
    return _text(value).casefold()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - upstream Zenodo publishes MD5 for this file
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            frame.to_csv(zipped, index=False)


def iter_classification(backbone_zip: Path) -> Iterator[dict[str, str]]:
    with zipfile.ZipFile(backbone_zip) as archive:
        bad = archive.testzip()
        if bad:
            raise ValueError(f"corrupt WFO archive member: {bad}")
        members = set(archive.namelist())
        if WFO_MEMBER not in members:
            raise ValueError(f"WFO archive missing {WFO_MEMBER}")
        with archive.open(WFO_MEMBER) as raw:
            handle = io.TextIOWrapper(raw, encoding="utf-8-sig", errors="replace", newline="")
            yield from csv.DictReader(handle, delimiter="\t")


def _row_summary(row: dict[str, str]) -> dict[str, str]:
    return {
        "taxon_id": _text(row.get("taxonID")),
        "scientific_name": _text(row.get("scientificName")),
        "family": _text(row.get("family")),
        "rank": _text(row.get("taxonRank")),
        "status": _text(row.get("taxonomicStatus")),
        "accepted_name_usage_id": _text(row.get("acceptedNameUsageID")),
        "kingdom": _text(row.get("kingdom")),
    }


def collect_exact_matches(
    backbone_zip: Path,
    queue: pd.DataFrame,
) -> tuple[dict[str, list[dict[str, str]]], dict[str, dict[str, str]], int]:
    query_keys = {_key(name) for name in queue["submitted_name"] if _key(name)}
    exact: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    accepted_ids: set[str] = set()
    rows_scanned = 0
    for raw in iter_classification(backbone_zip):
        rows_scanned += 1
        scientific_name = _text(raw.get("scientificName"))
        key = scientific_name.casefold()
        if key not in query_keys:
            continue
        row = _row_summary(raw)
        exact[key].append(row)
        if row["status"].casefold() == "synonym" and row["accepted_name_usage_id"]:
            accepted_ids.add(row["accepted_name_usage_id"].casefold())

    accepted_targets: dict[str, dict[str, str]] = {}
    if accepted_ids:
        for raw in iter_classification(backbone_zip):
            taxon_id = _text(raw.get("taxonID"))
            if taxon_id.casefold() in accepted_ids:
                accepted_targets[taxon_id.casefold()] = _row_summary(raw)
    return dict(exact), accepted_targets, rows_scanned


def _family_concordance(submitted_family: str, target_family: str) -> str:
    submitted = _text(submitted_family)
    target = _text(target_family)
    if not submitted:
        return "unknown_submitted_family"
    if not target:
        return "missing_wfo_family"
    return "match" if submitted.casefold() == target.casefold() else "conflict"


def _colxr_species_target(colxr_name: str, colxr_rank: str) -> str:
    """Return the species concept implied by a COL XR target, if identifiable.

    A species-rank target is already directly comparable.  An infraspecific
    canonical name still carries a parent species concept in its first two
    botanical name tokens. Higher-rank targets do not identify a species and
    therefore must not be treated as disagreement with an exact WFO species.
    """
    name = _text(colxr_name)
    rank = _text(colxr_rank).upper()
    if rank == "SPECIES":
        return name
    if rank in INFRA_RANKS and name:
        tokens = name.split()
        if (
            len(tokens) >= 2
            and tokens[0][:1].isupper()
            and tokens[1][:1].islower()
        ):
            return " ".join(tokens[:2])
    return ""


def _target_concordance(colxr_name: str, colxr_rank: str, wfo_name: str) -> tuple[str, str]:
    colxr_species = _colxr_species_target(colxr_name, colxr_rank)
    wfo = _text(wfo_name)
    if not colxr_species:
        return "no_colxr_species_target", ""
    if not wfo:
        return "missing_wfo_target", colxr_species
    state = "same_target" if colxr_species.casefold() == wfo.casefold() else "different_target"
    return state, colxr_species


def audit_queue(queue: pd.DataFrame, backbone_zip: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {
        "submitted_name",
        "submitted_family",
        "triage_class",
        "candidate_accepted_name",
        "candidate_accepted_rank",
        "second_backbone_eligible",
    }
    missing = sorted(required - set(queue.columns))
    if missing:
        raise ValueError(f"second-backbone queue missing columns: {missing}")
    if len(queue) != EXPECTED_QUEUE_ROWS:
        raise ValueError(f"second-backbone queue changed: {len(queue)} != {EXPECTED_QUEUE_ROWS}")
    if not queue["second_backbone_eligible"].astype(str).eq("true").all():
        raise ValueError("queue contains rows not qualified for second-backbone audit")
    if queue["submitted_name"].duplicated().any():
        raise ValueError("second-backbone queue contains duplicate submitted names")

    exact_by_name, accepted_targets, rows_scanned = collect_exact_matches(backbone_zip, queue)
    rows: list[dict[str, Any]] = []
    for source in queue.to_dict("records"):
        name = _text(source["submitted_name"])
        source_family = _text(source.get("submitted_family", ""))
        exact_rows = exact_by_name.get(name.casefold(), [])
        concepts: dict[str, dict[str, str]] = {}
        exact_species_rows = 0
        for match in exact_rows:
            if match["rank"].casefold() != "species":
                continue
            status = match["status"].casefold()
            if status not in {"accepted", "synonym"}:
                continue
            exact_species_rows += 1
            if status == "accepted":
                target = match
            else:
                target = accepted_targets.get(match["accepted_name_usage_id"].casefold(), {})
            if not target:
                continue
            if target.get("rank", "").casefold() != "species":
                continue
            if target.get("status", "").casefold() != "accepted":
                continue
            target_id = _key(target.get("taxon_id", ""))
            if target_id:
                concepts[target_id] = target

        concept_list = [concepts[key] for key in sorted(concepts)]
        target_state = ""
        colxr_species_target = ""
        if len(concept_list) == 0:
            wfo_status = "no_exact_accepted_species_concept"
            target = {}
        elif len(concept_list) > 1:
            wfo_status = "ambiguous_exact_name_multiple_species_targets"
            target = {}
        else:
            target = concept_list[0]
            family_state = _family_concordance(source_family, target.get("family", ""))
            target_state, colxr_species_target = _target_concordance(
                source.get("candidate_accepted_name", ""),
                source.get("candidate_accepted_rank", ""),
                target.get("scientific_name", ""),
            )
            if family_state == "conflict":
                wfo_status = "exact_species_family_conflict"
            elif family_state != "match":
                wfo_status = "exact_species_family_unverified"
            elif target_state == "different_target":
                wfo_status = "exact_species_colxr_target_conflict"
            else:
                wfo_status = "exact_species_rescue_candidate"

        target_name = _text(target.get("scientific_name", "")) if target else ""
        family_state = _family_concordance(source_family, target.get("family", "")) if target else ""
        rescue = wfo_status == "exact_species_rescue_candidate"
        rows.append(
            {
                **source,
                "wfo_exact_name_rows": len(exact_rows),
                "wfo_exact_species_rows": exact_species_rows,
                "wfo_unique_accepted_species_targets": len(concept_list),
                "wfo_candidate_target_ids": "|".join(
                    _text(item.get("taxon_id", "")) for item in concept_list
                ),
                "wfo_candidate_target_names": "|".join(
                    _text(item.get("scientific_name", "")) for item in concept_list
                ),
                "wfo_resolution_status": wfo_status,
                "wfo_accepted_taxon_id": _text(target.get("taxon_id", "")) if target else "",
                "wfo_accepted_name": target_name,
                "wfo_accepted_family": _text(target.get("family", "")) if target else "",
                "wfo_accepted_rank": _text(target.get("rank", "")) if target else "",
                "colxr_species_target_for_comparison": colxr_species_target,
                "wfo_family_concordance": family_state,
                "wfo_colxr_target_concordance": target_state,
                "wfo_rescue_candidate": "true" if rescue else "false",
                "wfo_review_status": "candidate_not_promoted" if rescue else "manual_review_required",
                "wfo_release": WFO_RELEASE,
                "wfo_doi": WFO_DOI,
                "wfo_license": WFO_LICENSE,
            }
        )

    audit = pd.DataFrame(rows)
    summary = {
        "n_input": int(len(audit)),
        "n_wfo_rescue_candidates": int(audit["wfo_rescue_candidate"].eq("true").sum()),
        "n_manual_review": int(audit["wfo_rescue_candidate"].ne("true").sum()),
        "wfo_resolution_status_counts": (
            audit["wfo_resolution_status"].value_counts().sort_index().to_dict()
        ),
        "wfo_rows_scanned": int(rows_scanned),
    }
    return audit, summary


def write_bundle(queue_path: Path, backbone_zip: Path, output_dir: Path) -> dict[str, Any]:
    actual_md5 = md5_file(backbone_zip)
    if actual_md5 != WFO_EXPECTED_MD5:
        raise ValueError(
            f"WFO 2026-06 archive MD5 mismatch: {actual_md5} != {WFO_EXPECTED_MD5}"
        )
    queue = pd.read_csv(queue_path, dtype=str).fillna("")
    audit, summary = audit_queue(queue, backbone_zip)
    rescue = audit.loc[audit["wfo_rescue_candidate"].eq("true")].copy()
    unresolved = audit.loc[audit["wfo_rescue_candidate"].ne("true")].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "wfo_second_backbone_audit.csv.gz": audit,
        "wfo_rescue_candidates.csv.gz": rescue,
        "wfo_unresolved_queue.csv.gz": unresolved,
    }
    files: dict[str, dict[str, Any]] = {}
    for name, frame in outputs.items():
        path = output_dir / name
        _write_gzip_csv(frame, path)
        files[name] = {"rows": int(len(frame)), "sha256": sha256_file(path)}

    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha3-wfo-second-backbone",
        "source_layer": "2.0.0-alpha3-taxonomy-triage",
        "source_queue_rows": int(len(queue)),
        "wfo_release": WFO_RELEASE,
        "wfo_zenodo_record": WFO_ZENODO_RECORD,
        "wfo_doi": WFO_DOI,
        "wfo_archive_md5": WFO_EXPECTED_MD5,
        "wfo_archive_sha256": sha256_file(backbone_zip),
        "wfo_license": WFO_LICENSE,
        "summary": summary,
        "files": files,
        "scientific_boundary": {
            "alpha1_mutated": False,
            "alpha2_taxonomy_mutated": False,
            "taxa_promoted": 0,
            "taxa_removed": 0,
            "fuzzy_matching_used": False,
            "wfo_results_are_candidates": True,
        },
    }
    (output_dir / "WFO_SECOND_BACKBONE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("audit")
def audit(
    queue: Path = typer.Option(..., "--queue", exists=True, dir_okay=False),
    wfo_zip: Path = typer.Option(..., "--wfo-zip", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
) -> None:
    manifest = write_bundle(queue, wfo_zip, output_dir)
    typer.echo(json.dumps(manifest["summary"], sort_keys=True))


if __name__ == "__main__":
    app()
