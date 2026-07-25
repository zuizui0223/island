"""Re-audit Validated Low genus consensus after new direct evidence arrives.

The baseline and augmented direct-evidence sets are evaluated against the same
currently unresolved species-axis tasks.  Only keys that become eligible after
adding the new evidence are emitted, preventing a re-run from claiming old
genus rules as a new contribution.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer

from island_v2.direct_web_reacquisition import (
    _coverage_snapshot,
    _file_sha256,
    _text,
)
from island_v2.validated_low_inference import Thresholds, infer_low_evidence

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _direct_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    evidence = frame.copy().fillna("")
    if "evidence_quality" not in evidence:
        evidence["evidence_quality"] = ""
    if "quality" in evidence:
        missing_quality = evidence["evidence_quality"].map(_text).eq("")
        evidence.loc[missing_quality, "evidence_quality"] = evidence.loc[
            missing_quality,
            "quality",
        ]
    for column in (
        "accepted_species",
        "trait_name",
        "normalized_value",
        "evidence_quality",
        "evidence_scope",
    ):
        if column not in evidence:
            evidence[column] = ""
    evidence["evidence_quality"] = evidence["evidence_quality"].map(
        lambda value: _text(value).casefold()
    )
    return evidence.loc[
        evidence["evidence_quality"].isin({"high", "medium"})
        & evidence["evidence_scope"].isin({"species_direct", "synonym_direct"})
    ].copy()


def incremental_genus_reaudit(
    *,
    unresolved: pd.DataFrame,
    baseline_evidence: pd.DataFrame,
    new_evidence: pd.DataFrame,
    thresholds: Thresholds | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return new Validated Low rows caused only by the added direct evidence."""
    thresholds = thresholds or Thresholds()
    tasks = unresolved[["accepted_species", "axis"]].copy().fillna("")
    tasks = tasks.drop_duplicates(["accepted_species", "axis"])
    baseline = _direct_evidence(baseline_evidence)
    added = _direct_evidence(new_evidence)
    augmented = pd.concat([baseline, added], ignore_index=True).drop_duplicates()

    base_low, _, base_rules = infer_low_evidence(tasks, baseline, thresholds)
    augmented_low, _, augmented_rules = infer_low_evidence(tasks, augmented, thresholds)
    base_keys = set(zip(base_low["accepted_species"], base_low["axis"]))
    augmented_keys = set(zip(augmented_low["accepted_species"], augmented_low["axis"]))
    incremental_keys = augmented_keys.difference(base_keys)
    incremental = augmented_low.loc[
        [
            (species, axis) in incremental_keys
            for species, axis in zip(
                augmented_low["accepted_species"],
                augmented_low["axis"],
            )
        ]
    ].copy()
    incremental = incremental.drop_duplicates(
        ["accepted_species", "axis", "trait_name", "normalized_value"]
    )

    report = {
        "contract_version": "incremental_validated_low_genus_reaudit_v1",
        "thresholds": {
            "min_species": thresholds.min_species,
            "min_dominance": thresholds.min_dominance,
            "min_masked_accuracy": thresholds.min_masked_accuracy,
            "colour_min_dominance": thresholds.colour_min_dominance,
            "reproductive_min_dominance": thresholds.reproductive_min_dominance,
        },
        "input_unresolved_species_axis": len(tasks),
        "baseline_direct_rows": len(baseline),
        "new_direct_rows": len(added),
        "baseline_hypothetical_low_species_axis": len(base_keys),
        "augmented_hypothetical_low_species_axis": len(augmented_keys),
        "incremental_validated_low_species_axis": len(incremental_keys),
        "incremental_by_axis": {
            axis: sum(candidate_axis == axis for _, candidate_axis in incremental_keys)
            for axis in sorted(tasks["axis"].unique())
        },
        "eligible_rules": {
            "baseline": (int(base_rules["eligible"].sum()) if "eligible" in base_rules else 0),
            "augmented": (
                int(augmented_rules["eligible"].sum()) if "eligible" in augmented_rules else 0
            ),
        },
        "policy": {
            "quality": "low",
            "evidence_scope": "genus_consensus",
            "species_direct_coverage_unchanged": True,
            "family_inference_used": False,
            "global_fallback_used": False,
            "baseline_rules_not_reclaimed_as_new": True,
        },
    }
    return incremental, augmented_rules, report


def run_reaudit(
    *,
    coverage_csv: Path,
    integrated_evidence_csv: Path,
    new_evidence_csvs: list[Path],
    output_dir: Path,
    thresholds: Thresholds | None = None,
) -> dict[str, Any]:
    if not new_evidence_csvs:
        raise ValueError("at least one new evidence CSV is required")
    thresholds = thresholds or Thresholds()
    coverage = pd.read_csv(coverage_csv, dtype=str).fillna("")
    quality_column = "combined_quality" if "combined_quality" in coverage else "after_quality"
    if quality_column not in coverage:
        raise ValueError("coverage lacks combined_quality and after_quality")
    unresolved = coverage.loc[
        coverage[quality_column].map(_text).eq(""),
        ["accepted_species", "axis"],
    ]
    baseline = pd.read_csv(integrated_evidence_csv, dtype=str).fillna("")
    new_parts = [pd.read_csv(path, dtype=str).fillna("") for path in new_evidence_csvs]
    incremental, rules, report = incremental_genus_reaudit(
        unresolved=unresolved,
        baseline_evidence=baseline,
        new_evidence=pd.concat(new_parts, ignore_index=True),
        thresholds=thresholds,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    accepted_path = output_dir / "incremental_validated_low_evidence.csv.gz"
    rules_path = output_dir / "augmented_validated_genus_rules.csv.gz"
    coverage_path = output_dir / "coverage_after_incremental_genus.csv.gz"
    incremental.to_csv(accepted_path, index=False, compression="gzip")
    rules.to_csv(rules_path, index=False, compression="gzip")
    incremental_keys = set(zip(incremental["accepted_species"], incremental["axis"]))
    coverage["quality_after_incremental_genus"] = coverage[quality_column].map(
        lambda value: _text(value).casefold()
    )
    promote = [
        not _text(quality) and (species, axis) in incremental_keys
        for species, axis, quality in zip(
            coverage["accepted_species"],
            coverage["axis"],
            coverage[quality_column],
        )
    ]
    coverage.loc[promote, "quality_after_incremental_genus"] = "low"
    coverage.to_csv(coverage_path, index=False, compression="gzip")
    before = _coverage_snapshot(coverage, quality_column)
    after = _coverage_snapshot(coverage, "quality_after_incremental_genus")
    report["coverage"] = {
        "before": before,
        "after": after,
        "net_filled_species_axis": (after["filled_species_axis"] - before["filled_species_axis"]),
        "net_unresolved_species_axis": (
            after["unresolved_species_axis"] - before["unresolved_species_axis"]
        ),
    }
    report["inputs"] = {
        "coverage_csv": {
            "path": str(coverage_csv),
            "sha256": _file_sha256(coverage_csv),
        },
        "integrated_evidence_csv": {
            "path": str(integrated_evidence_csv),
            "sha256": _file_sha256(integrated_evidence_csv),
        },
        "new_evidence_csvs": [
            {"path": str(path), "sha256": _file_sha256(path)} for path in new_evidence_csvs
        ],
    }
    report["outputs"] = {
        accepted_path.name: {
            "rows": len(incremental),
            "sha256": _file_sha256(accepted_path),
        },
        rules_path.name: {
            "rows": len(rules),
            "sha256": _file_sha256(rules_path),
        },
        coverage_path.name: {
            "rows": len(coverage),
            "sha256": _file_sha256(coverage_path),
        },
    }
    summary_path = output_dir / "incremental_genus_reaudit_summary.json"
    summary_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


@app.command("run")
def run_command(
    coverage_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    integrated_evidence_csv: Annotated[
        Path,
        typer.Option(exists=True, dir_okay=False),
    ],
    new_evidence_csv: Annotated[
        list[Path],
        typer.Option(
            "--new-evidence-csv",
            exists=True,
            dir_okay=False,
        ),
    ],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
    min_species: Annotated[int, typer.Option(min=2, max=100)] = 3,
    min_dominance: Annotated[float, typer.Option(min=0.5, max=1.0)] = 0.80,
    min_masked_accuracy: Annotated[
        float,
        typer.Option(min=0.5, max=1.0),
    ] = 0.75,
) -> None:
    report = run_reaudit(
        coverage_csv=coverage_csv,
        integrated_evidence_csv=integrated_evidence_csv,
        new_evidence_csvs=new_evidence_csv,
        output_dir=output_dir,
        thresholds=Thresholds(
            min_species=min_species,
            min_dominance=min_dominance,
            min_masked_accuracy=min_masked_accuracy,
        ),
    )
    typer.echo(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
