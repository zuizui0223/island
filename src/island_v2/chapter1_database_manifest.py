from __future__ import annotations

import hashlib
import json
import shlex
import subprocess
from pathlib import Path
from typing import Literal

import pandas as pd
import typer
import yaml
from pydantic import BaseModel, Field, model_validator

app = typer.Typer(help="Validate and dispatch versioned Chapter 1 database snapshots.")

REQUIRED_COLUMNS = {
    "accepted_species",
    "axis",
    "trait_composition",
    "trait_names",
    "source_groups",
    "source_lineages",
    "quality",
}
RIGHTS_REGISTRY_REQUIRED_COLUMNS = {
    "source_lineage",
    "source_family",
    "provider",
    "source_url",
    "dataset_doi",
    "license",
    "redistribution_status",
    "rights_evidence",
}
ALLOWED_QUALITY = {"high", "medium", "low", "unresolved", ""}
ALLOWED_RIGHTS_STATUS = {"redistributable", "review_required", "not_redistributable"}
CANONICAL_WORKFLOW = "run-chapter1-progressive-trait-analysis.yml"


class GitHubArtifactSource(BaseModel):
    kind: Literal["github_artifact"] = "github_artifact"
    repository: str
    run_id: int = Field(gt=0)
    artifact_name: str = Field(min_length=1)
    relative_path: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class RightsRegistrySpec(BaseModel):
    """Optional rights-aware provenance registry shipped with a database snapshot.

    Database 1.0 predates this contract and is valid without a registry. New database
    releases can opt in immediately; Database 2.0+ is expected to do so before the
    active pointer is advanced.
    """

    relative_path: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    required_for_public_release: bool = True


class PreviousSnapshot(BaseModel):
    run_id: int = Field(gt=0)
    artifact_name: str = Field(min_length=1)
    relative_path: str = Field(min_length=1)


class PublicationRecord(BaseModel):
    zenodo_doi: str | None = None
    zenodo_record_url: str | None = None
    release_status: Literal["not_deposited", "draft", "published"] = "not_deposited"


class DatabaseSpec(BaseModel):
    database_id: str
    version: str
    analysis_contract: str
    denominator_species: int = Field(gt=0)
    axes: list[str]
    resolved_cells: int = Field(ge=0)
    axis_resolved_cells: dict[str, int]
    source: GitHubArtifactSource
    rights_registry: RightsRegistrySpec | None = None
    previous_snapshot: PreviousSnapshot | None = None
    publication: PublicationRecord = Field(default_factory=PublicationRecord)

    @model_validator(mode="after")
    def validate_dimensions(self) -> "DatabaseSpec":
        if len(self.axes) != len(set(self.axes)):
            raise ValueError("database.axes must be unique")
        denominator = self.denominator_species * len(self.axes)
        if self.resolved_cells > denominator:
            raise ValueError("resolved_cells exceeds species-axis denominator")
        if set(self.axis_resolved_cells) != set(self.axes):
            raise ValueError("axis_resolved_cells keys must match axes exactly")
        if sum(self.axis_resolved_cells.values()) != self.resolved_cells:
            raise ValueError("axis_resolved_cells must sum to resolved_cells")
        return self


class Chapter1DatabaseManifest(BaseModel):
    schema_version: Literal[1]
    database: DatabaseSpec


def load_manifest(path: Path) -> Chapter1DatabaseManifest:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    return Chapter1DatabaseManifest.model_validate(payload)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tokenize_lineages(series: pd.Series) -> set[str]:
    lineages: set[str] = set()
    for value in series.fillna("").astype(str):
        lineages.update(token.strip() for token in value.split("|") if token.strip())
    return lineages


def validate_species_axis(path: Path, manifest: Chapter1DatabaseManifest) -> dict[str, object]:
    db = manifest.database
    if sha256_file(path) != db.source.sha256:
        raise ValueError("species-axis SHA-256 does not match the database manifest")

    frame = pd.read_csv(path, dtype=str).fillna("")
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"species-axis file is missing required columns: {sorted(missing)}")

    expected_rows = db.denominator_species * len(db.axes)
    if len(frame) != expected_rows:
        raise ValueError(f"expected {expected_rows} species-axis rows, found {len(frame)}")

    if frame[["accepted_species", "axis"]].duplicated().any():
        raise ValueError("species-axis file contains duplicate accepted_species × axis rows")

    observed_axes = set(frame["axis"].unique())
    if observed_axes != set(db.axes):
        raise ValueError(f"axis mismatch: expected {sorted(db.axes)}, found {sorted(observed_axes)}")

    species_per_axis = frame.groupby("axis")["accepted_species"].nunique().to_dict()
    if any(count != db.denominator_species for count in species_per_axis.values()):
        raise ValueError("one or more axes do not contain the fixed species denominator")

    quality = frame["quality"].str.strip().str.lower()
    invalid_quality = sorted(set(quality).difference(ALLOWED_QUALITY))
    if invalid_quality:
        raise ValueError(f"unsupported quality labels: {invalid_quality}")

    resolved = quality.isin({"high", "medium", "low"})
    resolved_cells = int(resolved.sum())
    if resolved_cells != db.resolved_cells:
        raise ValueError(
            f"resolved-cell mismatch: manifest={db.resolved_cells}, file={resolved_cells}"
        )

    axis_counts = (
        frame.loc[resolved]
        .groupby("axis", observed=True)
        .size()
        .reindex(db.axes, fill_value=0)
        .astype(int)
        .to_dict()
    )
    if axis_counts != db.axis_resolved_cells:
        raise ValueError(
            f"axis resolved-cell mismatch: manifest={db.axis_resolved_cells}, file={axis_counts}"
        )

    return {
        "database_id": db.database_id,
        "version": db.version,
        "analysis_contract": db.analysis_contract,
        "rows": len(frame),
        "resolved_cells": resolved_cells,
        "axis_resolved_cells": axis_counts,
        "sha256": db.source.sha256,
    }


def validate_rights_registry(
    path: Path,
    manifest: Chapter1DatabaseManifest,
    species_axis_path: Path | None = None,
) -> dict[str, object]:
    spec = manifest.database.rights_registry
    if spec is None:
        raise ValueError("database manifest does not declare a rights_registry")
    if sha256_file(path) != spec.sha256:
        raise ValueError("rights-registry SHA-256 does not match the database manifest")

    registry = pd.read_csv(path, dtype=str).fillna("")
    missing = RIGHTS_REGISTRY_REQUIRED_COLUMNS.difference(registry.columns)
    if missing:
        raise ValueError(f"rights registry is missing required columns: {sorted(missing)}")
    if registry["source_lineage"].eq("").any():
        raise ValueError("rights registry contains blank source_lineage values")
    if registry["source_lineage"].duplicated().any():
        raise ValueError("rights registry contains duplicate source_lineage rows")

    statuses = set(registry["redistribution_status"].str.strip().str.lower())
    invalid = sorted(statuses.difference(ALLOWED_RIGHTS_STATUS))
    if invalid:
        raise ValueError(f"unsupported redistribution_status labels: {invalid}")

    if species_axis_path is not None:
        species_axis = pd.read_csv(species_axis_path, dtype=str).fillna("")
        resolved = species_axis[
            species_axis["quality"].str.strip().str.lower().isin({"high", "medium", "low"})
        ]
        required_lineages = tokenize_lineages(resolved["source_lineages"])
        registered_lineages = set(registry["source_lineage"])
        missing_lineages = sorted(required_lineages.difference(registered_lineages))
        if missing_lineages:
            preview = missing_lineages[:10]
            raise ValueError(
                f"rights registry does not cover {len(missing_lineages)} resolved source lineages; "
                f"examples={preview}"
            )

    redistributable = registry["redistribution_status"].str.lower().eq("redistributable")
    missing_license = redistributable & registry["license"].eq("")
    missing_evidence = redistributable & registry["rights_evidence"].eq("")
    if missing_license.any() or missing_evidence.any():
        raise ValueError(
            "redistributable rights-registry rows require both license and rights_evidence"
        )

    return {
        "rows": int(len(registry)),
        "sha256": spec.sha256,
        "redistributable_rows": int(redistributable.sum()),
        "review_required_rows": int(
            registry["redistribution_status"].str.lower().eq("review_required").sum()
        ),
        "not_redistributable_rows": int(
            registry["redistribution_status"].str.lower().eq("not_redistributable").sum()
        ),
    }


def dispatch_command(
    manifest: Chapter1DatabaseManifest,
    workflow_ref: str = "main",
) -> list[str]:
    db = manifest.database
    src = db.source
    cmd = [
        "gh",
        "workflow",
        "run",
        CANONICAL_WORKFLOW,
        "--repo",
        src.repository,
        "--ref",
        workflow_ref,
        "-f",
        f"trait_run_id={src.run_id}",
        "-f",
        f"trait_artifact_name={src.artifact_name}",
        "-f",
        f"species_axis_relative_path={src.relative_path}",
    ]
    if db.previous_snapshot is not None:
        prev = db.previous_snapshot
        cmd.extend(
            [
                "-f",
                f"previous_run_id={prev.run_id}",
                "-f",
                f"previous_artifact_name={prev.artifact_name}",
                "-f",
                f"previous_species_axis_relative_path={prev.relative_path}",
            ]
        )
    return cmd


@app.command("validate")
def validate_command(
    manifest_path: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    species_axis: Path | None = typer.Option(None, "--species-axis", dir_okay=False),
    rights_registry: Path | None = typer.Option(None, "--rights-registry", dir_okay=False),
) -> None:
    """Validate manifest metadata and optionally materialized database files."""
    manifest = load_manifest(manifest_path)
    report: dict[str, object] = {
        "manifest": str(manifest_path),
        "database_id": manifest.database.database_id,
        "version": manifest.database.version,
        "analysis_contract": manifest.database.analysis_contract,
        "manifest_valid": True,
    }
    if species_axis is not None:
        report["species_axis"] = validate_species_axis(species_axis, manifest)
    if rights_registry is not None:
        report["rights_registry"] = validate_rights_registry(
            rights_registry,
            manifest,
            species_axis_path=species_axis,
        )
    typer.echo(json.dumps(report, indent=2, sort_keys=True))


@app.command("dispatch")
def dispatch_command_cli(
    manifest_path: Path = typer.Option(..., "--manifest", exists=True, dir_okay=False),
    workflow_ref: str = typer.Option("main", "--ref"),
    execute: bool = typer.Option(False, "--execute", help="Actually dispatch via the GitHub CLI."),
) -> None:
    """Dispatch the frozen Chapter 1 workflow using only the selected database manifest."""
    manifest = load_manifest(manifest_path)
    if manifest.database.analysis_contract != "chapter1_progressive_analysis_v1":
        raise typer.BadParameter(
            "manifest analysis_contract is not chapter1_progressive_analysis_v1"
        )
    cmd = dispatch_command(manifest, workflow_ref=workflow_ref)
    typer.echo(shlex.join(cmd))
    if execute:
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    app()
