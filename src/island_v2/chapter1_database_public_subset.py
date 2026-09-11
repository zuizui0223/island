from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False)
KEYS = ["accepted_species", "axis"]
PUBLIC_STATUS = "redistributable"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_public_subset(
    species_axis_path: Path,
    cell_rights_path: Path,
) -> tuple[pd.DataFrame, dict[str, int]]:
    species_axis = pd.read_csv(species_axis_path, dtype=str).fillna("")
    rights = pd.read_csv(cell_rights_path, dtype=str).fillna("")

    for name, frame in (("species axis", species_axis), ("cell rights", rights)):
        missing = [column for column in KEYS if column not in frame.columns]
        if missing:
            raise ValueError(f"{name} missing key columns: {missing}")
        if frame.duplicated(KEYS).any():
            raise ValueError(f"{name} has duplicate accepted_species x axis keys")

    required_rights = {"release_status", "quality"}
    missing_rights = sorted(required_rights - set(rights.columns))
    if missing_rights:
        raise ValueError(f"cell rights missing columns: {missing_rights}")

    resolved = species_axis[
        species_axis["quality"].str.lower().isin({"high", "medium", "low"})
    ].copy()
    merged = resolved.merge(
        rights[KEYS + ["quality", "release_status"]],
        on=KEYS,
        how="left",
        validate="one_to_one",
        suffixes=("", "_rights"),
        indicator=True,
    )
    if (merged["_merge"] != "both").any():
        missing_count = int((merged["_merge"] != "both").sum())
        raise ValueError(f"cell rights missing {missing_count} resolved species-axis rows")
    if not merged["quality"].str.lower().eq(merged["quality_rights"].str.lower()).all():
        raise ValueError("quality mismatch between Database 1.0 and cell-rights ledger")

    public_keys = rights.loc[rights["release_status"].eq(PUBLIC_STATUS), KEYS]
    subset = species_axis.merge(public_keys, on=KEYS, how="inner", validate="one_to_one")
    if subset.empty:
        raise ValueError("public subset would be empty")
    if not subset["quality"].str.lower().isin({"high", "medium", "low"}).all():
        raise ValueError("public subset contains unresolved cells")

    counts = {
        "resolved_cells": int(len(resolved)),
        "redistributable_cells": int(len(subset)),
        "blocked_or_review_cells": int(len(resolved) - len(subset)),
    }
    return subset, counts


@app.command("build")
def build(
    source_dir: Path = typer.Option(..., "--source-dir", exists=True, file_okay=False),
    cell_rights: Path = typer.Option(..., "--cell-rights", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
    rights_run_id: int = typer.Option(..., "--rights-run-id"),
    rights_artifact_id: int = typer.Option(..., "--rights-artifact-id"),
    rights_artifact_digest: str = typer.Option(..., "--rights-artifact-digest"),
    source_database_sha256: str = typer.Option(..., "--source-database-sha256"),
    public_release_license: str = typer.Option("CC-BY-SA-4.0", "--public-release-license"),
) -> None:
    source_path = source_dir / "species_axis_coverage.csv.gz"
    if not source_path.exists():
        raise FileNotFoundError(f"missing immutable Database 1.0 species axis: {source_path}")

    subset, counts = build_public_subset(source_path, cell_rights)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / "species_axis_coverage.public.csv.gz"
    subset.to_csv(out, index=False, compression="gzip")

    manifest = {
        "schema_version": 1,
        "database_id": "chapter1_database_v1_public_subset",
        "derivation": "rights-filtered derivative of immutable Chapter 1 Database 1.0",
        "scientific_values_modified": False,
        "selection_rule": "release_status == redistributable in CELL_RELEASE_RIGHTS.csv.gz",
        "public_release_license": public_release_license,
        "source_database_species_axis_sha256": source_database_sha256,
        "rights_audit_run_id": rights_run_id,
        "rights_audit_artifact_id": rights_artifact_id,
        "rights_audit_artifact_digest": rights_artifact_digest,
        **counts,
        "files": [
            {
                "file": out.name,
                "sha256": sha256_file(out),
                "bytes": out.stat().st_size,
            }
        ],
    }
    (output_dir / "PUBLIC_SUBSET_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "README.md").write_text(
        "# Chapter 1 Database 1.0 public subset\n\n"
        "This is a rights-filtered derivative of the immutable scientific Database 1.0. "
        "No scientific trait value is changed. Only resolved species × axis cells whose "
        "post-release audit status is `redistributable` are included. Blocked and review-required "
        "cells are omitted rather than reinterpreted.\n\n"
        f"Included cells: **{counts['redistributable_cells']:,}** of "
        f"{counts['resolved_cells']:,}** resolved Database 1.0 cells.\n\n"
        f"Compilation licence: **{public_release_license}**. Preserve all row-level provenance "
        "and upstream attribution/licence obligations.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    app()
