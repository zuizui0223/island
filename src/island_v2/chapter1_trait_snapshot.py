"""Materialize a versioned Chapter 1 trait snapshot from species-by-axis coverage.

The scientific analysis contract is fixed elsewhere.  This module makes the trait
layer replaceable: every acquisition wave is converted to the same long-form ledgers,
coverage report, provenance manifest, and (optionally) a transition audit versus the
previous snapshot.

A global fill percentage is descriptive only.  Analysis eligibility is decided later
from response-specific island support (pilot >=30; confirmatory >=50), not from a
species-level completion target.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXES = (
    "flower_colour",
    "floral_structural_complexity",
    "reproductive_assurance",
)
TRAIT_TO_AXIS = {
    "flower_primary_color": "flower_colour",
    "floral_form": "floral_structural_complexity",
    "floral_symmetry": "floral_structural_complexity",
    "tube_depth_class": "floral_structural_complexity",
    "flower_size_class": "floral_structural_complexity",
    "inflorescence_display": "floral_structural_complexity",
    "self_incompatibility": "reproductive_assurance",
    "autonomous_selfing_capacity": "reproductive_assurance",
    "mating_system": "reproductive_assurance",
    "cleistogamy": "reproductive_assurance",
}
QUALITY_RANK = {"": 0, "unresolved": 0, "low": 1, "medium": 2, "high": 3}
DIRECT_QUALITY = {"high", "medium"}
ANALYSIS_QUALITY = {"high", "medium", "low"}
GZIP = {"method": "gzip", "mtime": 0}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _quality(value: object) -> str:
    value = str(value or "").strip().casefold()
    if value not in QUALITY_RANK:
        raise ValueError(f"unknown quality: {value!r}")
    return "" if value == "unresolved" else value


def _states(value: str, *, label: str) -> list[str]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label}: state set is not JSON: {value!r}") from exc
    if not isinstance(parsed, list) or not parsed:
        raise ValueError(f"{label}: state set must be a non-empty JSON list")
    states = [str(item).strip() for item in parsed if str(item).strip()]
    if not states or len(states) != len(parsed):
        raise ValueError(f"{label}: state set contains an empty value")
    return sorted(set(states))


def _parse_composition(value: object, *, axis: str, species: str) -> list[dict[str, str]]:
    text = str(value or "").strip()
    if not text:
        return []
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in text.split("|"):
        trait_name, separator, state_text = item.partition("=")
        trait_name = trait_name.strip()
        if not separator:
            raise ValueError(f"{species} {axis}: malformed composition item {item!r}")
        if TRAIT_TO_AXIS.get(trait_name) != axis:
            raise ValueError(
                f"{species} {axis}: trait {trait_name!r} is outside the declared axis"
            )
        if trait_name in seen:
            raise ValueError(f"{species} {axis}: duplicate trait {trait_name}")
        seen.add(trait_name)
        states = _states(state_text, label=f"{species} {trait_name}")
        rows.append(
            {
                "accepted_species": species,
                "axis": axis,
                "trait_name": trait_name,
                "normalized_value": "|".join(states),
            }
        )
    return rows


def validate_species_axis(frame: pd.DataFrame, *, expected_species: int) -> pd.DataFrame:
    required = {
        "accepted_species",
        "axis",
        "trait_composition",
        "trait_names",
        "source_groups",
        "source_lineages",
        "quality",
    }
    if missing := required - set(frame.columns):
        raise ValueError(f"species-axis snapshot lacks columns: {sorted(missing)}")

    work = frame.copy().fillna("")
    for column in required:
        work[column] = work[column].astype(str)
    work["quality"] = work["quality"].map(_quality)

    if work[["accepted_species", "axis"]].duplicated().any():
        raise ValueError("duplicate accepted_species x axis cell")
    if work["accepted_species"].eq("").any():
        raise ValueError("blank accepted_species")
    if work["accepted_species"].nunique() != expected_species:
        raise ValueError(
            f"expected {expected_species} species, got {work['accepted_species'].nunique()}"
        )
    if set(work["axis"]) != set(AXES):
        raise ValueError(f"axis mismatch: {sorted(set(work['axis']))}")
    counts = work.groupby("axis")["accepted_species"].nunique().to_dict()
    if any(counts.get(axis) != expected_species for axis in AXES):
        raise ValueError(f"non-exact axis denominator: {counts}")

    resolved = work["quality"].ne("")
    if work.loc[resolved, "trait_composition"].eq("").any():
        raise ValueError("resolved species-axis cell lacks trait_composition")
    if work.loc[~resolved, "trait_composition"].ne("").any():
        raise ValueError("unresolved species-axis cell contains trait_composition")
    if work.loc[resolved, "source_groups"].eq("").any():
        raise ValueError("resolved species-axis cell lacks source_groups")
    if work.loc[resolved, "source_lineages"].eq("").any():
        raise ValueError("resolved species-axis cell lacks source_lineages")

    # Parse every resolved composition now so downstream analyses never discover
    # an ontology/cross-axis error after the snapshot has been declared.
    for row in work.loc[resolved].itertuples(index=False):
        parsed = _parse_composition(
            row.trait_composition,
            axis=row.axis,
            species=row.accepted_species,
        )
        parsed_names = {item["trait_name"] for item in parsed}
        declared_names = {x for x in str(row.trait_names).split("|") if x}
        if parsed_names != declared_names:
            raise ValueError(
                f"{row.accepted_species} {row.axis}: trait_names do not match composition"
            )
    return work.sort_values(["accepted_species", "axis"]).reset_index(drop=True)


def materialize_long_ledger(frame: pd.DataFrame, *, direct_only: bool) -> pd.DataFrame:
    allowed = DIRECT_QUALITY if direct_only else ANALYSIS_QUALITY
    rows: list[dict[str, str]] = []
    for row in frame.loc[frame["quality"].isin(allowed)].itertuples(index=False):
        for parsed in _parse_composition(
            row.trait_composition,
            axis=row.axis,
            species=row.accepted_species,
        ):
            rows.append(
                {
                    **parsed,
                    "quality": row.quality,
                    "source_groups": row.source_groups,
                    "source_lineages": row.source_lineages,
                    "evidence_scope": "direct" if row.quality in DIRECT_QUALITY else "validated_low",
                }
            )
    columns = [
        "accepted_species",
        "axis",
        "trait_name",
        "normalized_value",
        "quality",
        "source_groups",
        "source_lineages",
        "evidence_scope",
    ]
    ledger = pd.DataFrame(rows, columns=columns)
    if not ledger.empty and ledger[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("long ledger contains duplicate accepted_species x trait_name")
    return ledger.sort_values(["accepted_species", "axis", "trait_name"]).reset_index(drop=True)


def coverage_summary(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for axis in AXES:
        part = frame.loc[frame["axis"].eq(axis)]
        resolved = part["quality"].ne("")
        row: dict[str, Any] = {
            "axis": axis,
            "denominator_species": int(part["accepted_species"].nunique()),
            "resolved_species": int(resolved.sum()),
            "fill_fraction": float(resolved.mean()),
            "high": int(part["quality"].eq("high").sum()),
            "medium": int(part["quality"].eq("medium").sum()),
            "low": int(part["quality"].eq("low").sum()),
            "unresolved": int(part["quality"].eq("").sum()),
        }
        rows.append(row)
    table = pd.DataFrame(rows)
    resolved_cells = int(frame["quality"].ne("").sum())
    summary = {
        "species": int(frame["accepted_species"].nunique()),
        "species_axis_cells": int(len(frame)),
        "resolved_species_axis_cells": resolved_cells,
        "global_species_axis_fill_fraction": resolved_cells / len(frame),
        "global_fill_fraction_is_descriptive_not_an_analysis_gate": True,
        "by_axis": {row["axis"]: row for row in rows},
    }
    return table, summary


def transition_audit(previous: pd.DataFrame, current: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    cols = ["accepted_species", "axis", "quality", "trait_composition"]
    joined = previous[cols].merge(
        current[cols],
        on=["accepted_species", "axis"],
        how="outer",
        suffixes=("_previous", "_current"),
        validate="one_to_one",
        indicator=True,
    )
    if not joined["_merge"].eq("both").all():
        raise ValueError("species-axis denominator changed between snapshots")

    previous_rank = joined["quality_previous"].map(QUALITY_RANK)
    current_rank = joined["quality_current"].map(QUALITY_RANK)
    joined["transition"] = "unchanged"
    joined.loc[(previous_rank == 0) & (current_rank > 0), "transition"] = "newly_resolved"
    joined.loc[(previous_rank > 0) & (current_rank == 0), "transition"] = "became_unresolved"
    joined.loc[current_rank > previous_rank, "transition"] = "quality_upgrade"
    joined.loc[(current_rank < previous_rank) & (current_rank > 0), "transition"] = "quality_downgrade"
    changed_value = (
        previous_rank.gt(0)
        & current_rank.gt(0)
        & joined["trait_composition_previous"].ne(joined["trait_composition_current"])
    )
    joined.loc[changed_value, "transition"] = "value_revision"
    counts = {
        str(key): int(value)
        for key, value in joined["transition"].value_counts().sort_index().items()
    }
    return joined.drop(columns=["_merge"]), counts


def build_snapshot(
    *,
    species_axis_csv: Path,
    output_dir: Path,
    expected_species: int,
    source_run_id: str,
    source_artifact_name: str,
    previous_species_axis_csv: Path | None = None,
    allow_evidence_revision: bool = False,
) -> dict[str, Any]:
    current = validate_species_axis(
        pd.read_csv(species_axis_csv, dtype=str).fillna(""),
        expected_species=expected_species,
    )
    coverage_table, coverage = coverage_summary(current)
    all_ledger = materialize_long_ledger(current, direct_only=False)
    direct_ledger = materialize_long_ledger(current, direct_only=True)

    transition_counts: dict[str, int] | None = None
    transition = pd.DataFrame()
    previous_sha256: str | None = None
    if previous_species_axis_csv is not None:
        previous = validate_species_axis(
            pd.read_csv(previous_species_axis_csv, dtype=str).fillna(""),
            expected_species=expected_species,
        )
        transition, transition_counts = transition_audit(previous, current)
        previous_sha256 = _sha256(previous_species_axis_csv)
        forbidden = sum(
            transition_counts.get(label, 0)
            for label in ("became_unresolved", "quality_downgrade", "value_revision")
        )
        if forbidden and not allow_evidence_revision:
            raise ValueError(
                "snapshot contains evidence loss/revision; rerun with explicit "
                "--allow-evidence-revision after auditing the transition table"
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    current.to_csv(
        output_dir / "chapter1_trait_snapshot_species_axis.csv.gz",
        index=False,
        compression=GZIP,
    )
    all_ledger.to_csv(
        output_dir / "chapter1_trait_ledger_all_analysis.csv.gz",
        index=False,
        compression=GZIP,
    )
    direct_ledger.to_csv(
        output_dir / "chapter1_trait_ledger_direct_only.csv.gz",
        index=False,
        compression=GZIP,
    )
    coverage_table.to_csv(output_dir / "chapter1_trait_snapshot_coverage.csv", index=False)
    if not transition.empty:
        transition.to_csv(
            output_dir / "chapter1_trait_snapshot_transition.csv.gz",
            index=False,
            compression=GZIP,
        )

    manifest = {
        "contract": "chapter1_trait_snapshot_v1",
        "source_workflow_run_id": str(source_run_id),
        "source_artifact_name": str(source_artifact_name),
        "source_species_axis_file": species_axis_csv.name,
        "source_species_axis_sha256": _sha256(species_axis_csv),
        "previous_species_axis_sha256": previous_sha256,
        "expected_species": expected_species,
        "axes": list(AXES),
        "coverage": coverage,
        "n_all_analysis_trait_rows": int(len(all_ledger)),
        "n_direct_only_trait_rows": int(len(direct_ledger)),
        "transition_counts": transition_counts,
        "allow_evidence_revision": bool(allow_evidence_revision),
        "analysis_eligibility_rule": "response-specific island support; global fill fraction is descriptive only",
    }
    (output_dir / "chapter1_trait_snapshot_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


@app.command()
def run(
    species_axis_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option()],
    expected_species: Annotated[int, typer.Option()] = 106_295,
    source_run_id: Annotated[str, typer.Option()] = "manual",
    source_artifact_name: Annotated[str, typer.Option()] = "manual",
    previous_species_axis_csv: Annotated[
        Path | None, typer.Option(exists=True, dir_okay=False)
    ] = None,
    allow_evidence_revision: Annotated[bool, typer.Option()] = False,
) -> None:
    manifest = build_snapshot(
        species_axis_csv=species_axis_csv,
        output_dir=output_dir,
        expected_species=expected_species,
        source_run_id=source_run_id,
        source_artifact_name=source_artifact_name,
        previous_species_axis_csv=previous_species_axis_csv,
        allow_evidence_revision=allow_evidence_revision,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
