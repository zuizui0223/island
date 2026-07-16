"""Materialize source-specific GitHub artifacts into the all-master ledger layout."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

VALUE_COLUMNS = ("candidate_value", "standardized_value", "trait_value", "value")
EXCERPT_COLUMNS = ("evidence_excerpt", "source_excerpt")


@app.callback()
def main() -> None:
    """Materialize immutable source artifacts for exhaustive acquisition accounting."""


def _find_unique(root: Path, patterns: tuple[str, ...], label: str) -> Path:
    matches = {
        path.resolve() for pattern in patterns for path in root.glob(pattern) if path.is_file()
    }
    if len(matches) != 1:
        rendered = ", ".join(str(path) for path in sorted(matches)) or "none"
        raise ValueError(f"expected exactly one {label} candidate file; found {rendered}")
    return next(iter(matches))


def _first_column(frame: pd.DataFrame, choices: tuple[str, ...], label: str) -> str:
    for column in choices:
        if column in frame.columns:
            return column
    raise ValueError(f"{label} candidates are missing any of columns {choices}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _candidate_summary(path: Path, label: str) -> dict[str, Any]:
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"accepted_species", "trait_name", "source_url", "source_record_id"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} candidates are missing required columns: {missing}")
    value_column = _first_column(frame, VALUE_COLUMNS, label)
    excerpt_column = _first_column(frame, EXCERPT_COLUMNS, label)
    if frame.empty:
        raise ValueError(f"{label} candidate file is empty")
    evidence_columns = (
        "accepted_species",
        "trait_name",
        value_column,
        "source_url",
        "source_record_id",
        excerpt_column,
    )
    empty_counts = {
        column: int(frame[column].astype(str).str.strip().eq("").sum())
        for column in evidence_columns
    }
    if any(empty_counts.values()):
        raise ValueError(f"{label} candidates contain empty evidence fields: {empty_counts}")
    return {
        "n_candidate_rows": int(len(frame)),
        "n_species": int(frame["accepted_species"].nunique()),
        "n_species_trait_pairs": int(
            len(frame.drop_duplicates(["accepted_species", "trait_name"]))
        ),
        "species_by_trait": {
            str(trait): int(count)
            for trait, count in frame.groupby("trait_name")["accepted_species"]
            .nunique()
            .sort_values(ascending=False)
            .items()
        },
        "value_column": value_column,
        "excerpt_column": excerpt_column,
        "empty_evidence_fields": empty_counts,
    }


def _copy_source(source: Path, destination: Path, label: str) -> dict[str, Any]:
    summary = _candidate_summary(source, label)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    summary.update(
        {
            "artifact_candidate_path": str(source),
            "materialized_candidate_path": str(destination),
            "candidate_sha256": _sha256(destination),
        }
    )
    return summary


def materialize(
    *,
    eol_artifact_dir: Path,
    austraits_artifact_dir: Path,
    staging_root: Path,
    output_dir: Path,
    efloras_artifact_dir: Path | None = None,
    gift_artifact_dir: Path | None = None,
    usda_artifact_dir: Path | None = None,
) -> dict[str, Any]:
    """Copy exact source-backed candidates into paths read by acquisition-rate config."""
    eol_source = _find_unique(
        eol_artifact_dir,
        ("ingested/bulk_trait_candidates_*.csv", "**/ingested/bulk_trait_candidates_*.csv"),
        "EOL",
    )
    austraits_source = _find_unique(
        austraits_artifact_dir,
        (
            "**/data/v2/staging/traits/bulk/austraits_all/trait_candidates.csv.gz",
            "**/data/v2/staging/traits/bulk/austraits_all/trait_candidates.csv",
        ),
        "AusTraits",
    )
    sources = {
        "eol_traitbank": _copy_source(
            eol_source,
            staging_root / "eol_traitbank" / "trait_candidates.csv",
            "EOL",
        ),
        "austraits_all": _copy_source(
            austraits_source,
            staging_root / "bulk" / "austraits_all" / "trait_candidates.csv.gz",
            "AusTraits",
        ),
    }
    if efloras_artifact_dir is not None:
        efloras_source = _find_unique(
            efloras_artifact_dir,
            ("all_trait_candidates.csv", "**/all_trait_candidates.csv"),
            "eFloras",
        )
        sources["efloras"] = _copy_source(
            efloras_source,
            staging_root / "efloras" / "combined" / "trait_candidates.csv",
            "eFloras",
        )
    if gift_artifact_dir is not None:
        gift_source = _find_unique(
            gift_artifact_dir,
            (
                "**/data/v2/staging/traits/bulk/gift_v3_2/"
                "bulk_trait_candidates_gift_v3_2_direct.csv",
                "**/bulk_trait_candidates_gift_v3_2_direct.csv",
            ),
            "GIFT",
        )
        sources["gift_v3_2_direct"] = _copy_source(
            gift_source,
            staging_root / "bulk" / "gift_v3_2" / "trait_candidates.csv",
            "GIFT",
        )
    if usda_artifact_dir is not None:
        usda_source = _find_unique(
            usda_artifact_dir,
            (
                "**/data/v2/staging/traits/bulk/usda_plants_current/"
                "bulk_trait_candidates_usda_plants_current.csv",
                "**/bulk_trait_candidates_usda_plants_current.csv",
            ),
            "USDA PLANTS",
        )
        sources["usda_plants_current"] = _copy_source(
            usda_source,
            staging_root / "bulk" / "usda_plants_current" / "trait_candidates.csv",
            "USDA PLANTS",
        )

    manifest = {
        "version": "1.0",
        "source_count": len(sources),
        "sources": sources,
        "rules": {
            "population": "all master species remain in the downstream ledgers",
            "evidence": "candidate value, URL, excerpt, and source record ID must be nonempty",
            "inference": "materialization copies candidates and does not infer missing values",
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "materialization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


@app.command()
def run(
    eol_artifact_dir: Path = typer.Option(..., exists=True, file_okay=False),
    austraits_artifact_dir: Path = typer.Option(..., exists=True, file_okay=False),
    staging_root: Path = typer.Option(Path("data/v2/staging/traits")),
    output_dir: Path = typer.Option(...),
    efloras_artifact_dir: Path | None = typer.Option(None, file_okay=False),
    gift_artifact_dir: Path | None = typer.Option(None, file_okay=False),
    usda_artifact_dir: Path | None = typer.Option(None, file_okay=False),
) -> None:
    """Materialize candidate artifacts and emit a provenance/count manifest."""
    manifest = materialize(
        eol_artifact_dir=eol_artifact_dir,
        austraits_artifact_dir=austraits_artifact_dir,
        efloras_artifact_dir=efloras_artifact_dir,
        gift_artifact_dir=gift_artifact_dir,
        usda_artifact_dir=usda_artifact_dir,
        staging_root=staging_root,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
