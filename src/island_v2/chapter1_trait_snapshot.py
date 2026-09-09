"""Materialize a versioned Chapter 1 trait snapshot from species-by-axis coverage.

The hypothesis/model contract is fixed elsewhere.  This module makes the trait
layer replaceable: every acquisition wave becomes the same ledgers, coverage report,
provenance manifest, and optional transition audit.

Global fill percentage is descriptive only.  Analysis eligibility is response-specific
(pilot >=30 islands; confirmatory >=50 islands).
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
    text = str(value or "").strip().casefold()
    if text not in QUALITY_RANK:
        raise ValueError(f"unknown quality: {text!r}")
    return "" if text == "unresolved" else text


def _states(value: str, *, label: str) -> list[str]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label}: state set is not JSON: {value!r}") from exc
    if not isinstance(parsed, list) or not parsed:
        raise ValueError(f"{label}: state set must be a non-empty JSON list")
    states = [str(item).strip() for item in parsed]
    if any(not state for state in states):
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
        rows.append(
            {
                "accepted_species": species,
                "axis": axis,
                "trait_name": trait_name,
                "normalized_value": "|".join(
                    _states(state_text, label=f"{species} {trait_name}")
                ),
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

    materializable = work["trait_composition"].ne("")
    direct_claim_without_composition = work["quality"].isin(DIRECT_QUALITY) & ~materializable
    if direct_claim_without_composition.any():
        raise ValueError("direct-quality species-axis cell lacks trait_composition")
    if work.loc[materializable, "quality"].eq("").any():
        raise ValueError("trait composition lacks an analysis quality")
    if work.loc[materializable, "source_groups"].eq("").any():
        raise ValueError("resolved species-axis cell lacks source_groups")
    if work.loc[materializable, "source_lineages"].eq("").any():
        raise ValueError("resolved species-axis cell lacks source_lineages")

    for row in work.loc[materializable].itertuples(index=False):
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
    eligible = frame["quality"].isin(allowed) & frame["trait_composition"].ne("")
    for row in frame.loc[eligible].itertuples(index=False):
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
                    "evidence_scope": (
                        "direct" if row.quality in DIRECT_QUALITY else "validated_low"
                    ),
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
        source_reported = part["quality"].ne("")
        resolved = part["trait_composition"].ne("")
        nonmaterializable_low = part["quality"].eq("low") & ~resolved
        rows.append(
            {
                "axis": axis,
                "denominator_species": int(part["accepted_species"].nunique()),
                "resolved_species": int(resolved.sum()),
                "fill_fraction": float(resolved.mean()),
                "source_reported_filled_species": int(source_reported.sum()),
                "source_reported_fill_fraction": float(source_reported.mean()),
                "source_reported_low_without_trait_composition": int(
                    nonmaterializable_low.sum()
                ),
                "high": int(part["quality"].eq("high").sum()),
                "medium": int(part["quality"].eq("medium").sum()),
                "materializable_low": int((part["quality"].eq("low") & resolved).sum()),
                "unresolved_for_analysis": int((~resolved).sum()),
            }
        )
    table = pd.DataFrame(rows)
    resolved_cells = int(frame["trait_composition"].ne("").sum())
    source_reported_cells = int(frame["quality"].ne("").sum())
    return table, {
        "species": int(frame["accepted_species"].nunique()),
        "species_axis_cells": len(frame),
        "resolved_species_axis_cells": resolved_cells,
        "global_species_axis_fill_fraction": resolved_cells / len(frame),
        "source_reported_filled_species_axis_cells": source_reported_cells,
        "global_source_reported_fill_fraction": source_reported_cells / len(frame),
        "source_reported_low_without_trait_composition": int(
            (frame["quality"].eq("low") & frame["trait_composition"].eq("")).sum()
        ),
        "global_fill_fraction_is_descriptive_not_an_analysis_gate": True,
        "by_axis": {row["axis"]: row for row in rows},
    }


def transition_audit(
    previous: pd.DataFrame, current: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, int]]:
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
    previous_resolved = previous_rank.gt(0) & joined["trait_composition_previous"].ne("")
    current_resolved = current_rank.gt(0) & joined["trait_composition_current"].ne("")

    joined["transition"] = "unchanged"
    joined.loc[~previous_resolved & current_resolved, "transition"] = "newly_resolved"
    joined.loc[previous_resolved & ~current_resolved, "transition"] = "became_unresolved"
    joined.loc[
        previous_resolved & current_resolved & current_rank.gt(previous_rank),
        "transition",
    ] = "quality_upgrade"
    joined.loc[
        previous_resolved & current_resolved & current_rank.lt(previous_rank),
        "transition",
    ] = "quality_downgrade"
    joined.loc[
        previous_resolved
        & current_resolved
        & joined["trait_composition_previous"].ne(joined["trait_composition_current"]),
        "transition",
    ] = "value_revision"

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
        "n_all_analysis_trait_rows": len(all_ledger),
        "n_direct_only_trait_rows": len(direct_ledger),
        "transition_counts": transition_counts,
        "allow_evidence_revision": bool(allow_evidence_revision),
        "analysis_eligibility_rule": (
            "response-specific island support; global fill fraction is descriptive only"
        ),
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
