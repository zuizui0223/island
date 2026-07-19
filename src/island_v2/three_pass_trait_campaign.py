"""Plan a quality-ordered three-pass trait evidence campaign.

Pass 1 seeks high-quality species/synonym evidence from literature, floras, and
structured databases. Pass 2 seeks medium-quality species/synonym web
descriptions only for unresolved species-axis tasks. Pass 3 permits explicit
genus inference only for tasks still unresolved after passes 1 and 2. Family
inference is never created; remaining tasks are reported as unresolved.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXES = {
    "flower_colour": {"flower_primary_color", "flower_color"},
    "floral_structural_complexity": {
        "floral_form", "floral_symmetry", "tube_depth_class", "flower_size_class",
        "inflorescence_display", "flower_shape",
    },
    "reproductive_assurance": {
        "self_incompatibility", "autonomous_selfing_capacity", "mating_system", "cleistogamy",
    },
}
QUALITY_RANK = {"high": 3, "medium": 2, "low": 1}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _species_column(frame: pd.DataFrame) -> str:
    for candidate in ("accepted_species", "species", "accepted_name", "scientific_name"):
        if candidate in frame.columns:
            return candidate
    raise typer.BadParameter("master table has no species column")


def build_tasks(master: pd.DataFrame, evidence: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    species_col = _species_column(master)
    species = sorted({_text(value) for value in master[species_col] if _text(value)})

    ledger = evidence.copy()
    required = {"accepted_species", "trait_name", "evidence_quality", "evidence_scope"}
    for column in required:
        if column not in ledger.columns:
            ledger[column] = ""
    for column in ledger.columns:
        ledger[column] = ledger[column].map(_text)
    ledger = ledger.loc[
        ledger["accepted_species"].isin(species)
        & ledger["evidence_quality"].isin(QUALITY_RANK)
        & ~ledger["evidence_scope"].isin({"family_inference", "global_fallback"})
    ].copy()

    rows: list[dict[str, object]] = []
    counts = {"pass1_high": 0, "pass2_medium": 0, "pass3_low": 0, "unresolved": 0}
    for accepted_species in species:
        subset = ledger.loc[ledger["accepted_species"].eq(accepted_species)]
        genus = accepted_species.split()[0] if accepted_species else ""
        for axis, traits in AXES.items():
            axis_rows = subset.loc[subset["trait_name"].isin(traits)]
            high = axis_rows.loc[
                axis_rows["evidence_quality"].eq("high")
                & axis_rows["evidence_scope"].isin(DIRECT_SCOPES)
            ]
            medium = axis_rows.loc[
                axis_rows["evidence_quality"].eq("medium")
                & axis_rows["evidence_scope"].isin(DIRECT_SCOPES)
            ]
            low = axis_rows.loc[
                axis_rows["evidence_quality"].eq("low")
                & axis_rows["evidence_scope"].eq("genus_inference")
            ]

            if not high.empty:
                state, best_quality, next_pass = "resolved", "high", "none"
            elif not medium.empty:
                state, best_quality, next_pass = "resolved", "medium", "none"
            elif not low.empty:
                state, best_quality, next_pass = "resolved", "low", "none"
            else:
                state, best_quality = "unresolved", ""
                next_pass = "pass1_high"

            rows.append({
                "accepted_species": accepted_species,
                "genus": genus,
                "axis": axis,
                "state": state,
                "best_quality": best_quality,
                "next_pass": next_pass,
                "high_direct_present": not high.empty,
                "medium_direct_present": not medium.empty,
                "low_genus_present": not low.empty,
                "family_inference_allowed": False,
            })

    tasks = pd.DataFrame(rows)

    pass1 = tasks.loc[~tasks["high_direct_present"]].copy()
    pass1["target_quality"] = "high"
    pass1["pass"] = 1

    pass2 = tasks.loc[~tasks["high_direct_present"] & ~tasks["medium_direct_present"]].copy()
    pass2["target_quality"] = "medium"
    pass2["pass"] = 2

    pass3 = tasks.loc[
        ~tasks["high_direct_present"]
        & ~tasks["medium_direct_present"]
        & ~tasks["low_genus_present"]
    ].copy()
    pass3["target_quality"] = "low"
    pass3["pass"] = 3

    counts["pass1_high"] = int(len(pass1))
    counts["pass2_medium"] = int(len(pass2))
    counts["pass3_low"] = int(len(pass3))
    counts["unresolved"] = int((tasks["state"] == "unresolved").sum())
    summary = {
        "contract_version": "three_pass_species_axis_v1",
        "species": len(species),
        "species_axis_tasks": len(tasks),
        "passes": {
            "1": "high literature_flora_structured_database",
            "2": "medium species_or_synonym_web_description",
            "3": "low explicit_genus_inference",
        },
        "family_inference_allowed": False,
        "remaining_after_pass3": "unresolved",
        "task_counts": counts,
    }
    return pd.concat([pass1, pass2, pass3], ignore_index=True), summary


@app.command()
def build(
    master_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    evidence_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    evidence = pd.read_csv(evidence_csv, dtype=str).fillna("")
    tasks, summary = build_tasks(master, evidence)
    output_dir.mkdir(parents=True, exist_ok=True)
    for number, name in ((1, "high"), (2, "medium"), (3, "low")):
        tasks.loc[tasks["pass"].eq(number)].to_csv(
            output_dir / f"pass{number}_{name}_tasks.csv.gz", index=False
        )
    unresolved = tasks.loc[tasks["pass"].eq(3), ["accepted_species", "genus", "axis"]].drop_duplicates()
    unresolved.to_csv(output_dir / "unresolved_after_direct_evidence.csv.gz", index=False)
    (output_dir / "three_pass_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
