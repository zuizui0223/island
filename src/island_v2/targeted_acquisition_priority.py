"""Prioritize direct-trait acquisition for genera near validated-low eligibility."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXES = {
    "flower_colour": {"flower_primary_color", "flower_color"},
    "floral_structural_complexity": {
        "floral_form", "floral_symmetry", "tube_depth_class",
        "flower_size_class", "inflorescence_display", "flower_shape",
    },
    "reproductive_assurance": {
        "self_incompatibility", "autonomous_selfing_capacity",
        "mating_system", "cleistogamy",
    },
}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _genus(species: object) -> str:
    return _text(species).split(" ", 1)[0] if _text(species) else ""


def build_priority(
    evidence: pd.DataFrame,
    unresolved: pd.DataFrame,
    min_species: int = 3,
    max_target_species_per_genus_axis: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    evidence = evidence.copy()
    unresolved = unresolved.copy()

    for column in ("accepted_species", "trait_name", "evidence_quality", "evidence_scope"):
        if column not in evidence.columns:
            evidence[column] = ""
    if "accepted_species" not in unresolved.columns:
        unresolved["accepted_species"] = ""
    if "axis" not in unresolved.columns:
        unresolved["axis"] = ""

    evidence["accepted_species"] = evidence["accepted_species"].map(_text)
    evidence["trait_name"] = evidence["trait_name"].map(_text)
    evidence["genus"] = evidence["accepted_species"].map(_genus)
    unresolved["accepted_species"] = unresolved["accepted_species"].map(_text)
    unresolved["axis"] = unresolved["axis"].map(_text)
    unresolved["genus"] = unresolved["accepted_species"].map(_genus)

    direct = evidence.loc[
        evidence["evidence_quality"].isin({"high", "medium"})
        & evidence["evidence_scope"].isin(DIRECT_SCOPES)
        & evidence["genus"].ne("")
    ].copy()

    rows: list[dict[str, object]] = []
    target_frames: list[pd.DataFrame] = []

    for axis, traits in AXES.items():
        axis_unresolved = unresolved.loc[
            unresolved["axis"].eq(axis) & unresolved["genus"].ne("")
        ].drop_duplicates(["accepted_species", "axis"]).copy()
        if axis_unresolved.empty:
            continue

        axis_direct = direct.loc[direct["trait_name"].isin(traits)].copy()
        trait_counts = (
            axis_direct.groupby(["genus", "trait_name"])["accepted_species"]
            .nunique()
            .rename("n_direct_species")
            .reset_index()
        )
        if trait_counts.empty:
            continue

        best = (
            trait_counts.sort_values(
                ["genus", "n_direct_species", "trait_name"],
                ascending=[True, False, True],
            )
            .drop_duplicates("genus")
            .rename(columns={"trait_name": "best_supported_trait"})
        )
        unresolved_counts = (
            axis_unresolved.groupby("genus")["accepted_species"]
            .nunique()
            .rename("unresolved_species_in_genus_axis")
            .reset_index()
        )
        near = best.merge(unresolved_counts, on="genus", how="inner")
        near = near.loc[
            near["n_direct_species"].ge(1)
            & near["n_direct_species"].lt(min_species)
        ].copy()
        if near.empty:
            continue

        near["axis"] = axis
        near["additional_species_needed"] = min_species - near["n_direct_species"]
        near["leverage_score"] = (
            near["unresolved_species_in_genus_axis"]
            / near["additional_species_needed"]
        )
        near["target_traits"] = "|".join(sorted(traits))
        near = near.sort_values(
            ["leverage_score", "unresolved_species_in_genus_axis", "genus"],
            ascending=[False, False, True],
        )

        direct_species_by_genus = (
            axis_direct.groupby("genus")["accepted_species"].apply(set).to_dict()
        )
        for item in near.to_dict("records"):
            genus = str(item["genus"])
            needed = min(
                int(item["additional_species_needed"]),
                max_target_species_per_genus_axis,
            )
            candidates = axis_unresolved.loc[
                axis_unresolved["genus"].eq(genus)
                & ~axis_unresolved["accepted_species"].isin(
                    direct_species_by_genus.get(genus, set())
                )
            ].sort_values("accepted_species")
            chosen = candidates.head(needed).copy()
            if chosen.empty:
                continue
            for key, value in item.items():
                chosen[key] = value
            target_frames.append(chosen)
            rows.append(item)

    priority = pd.DataFrame(rows)
    if not priority.empty:
        priority = priority.drop_duplicates(["genus", "axis"])
        priority = priority.sort_values(
            ["leverage_score", "unresolved_species_in_genus_axis", "genus", "axis"],
            ascending=[False, False, True, True],
        ).reset_index(drop=True)
        priority.insert(0, "priority_rank", range(1, len(priority) + 1))

    targets = pd.concat(target_frames, ignore_index=True) if target_frames else unresolved.iloc[0:0].copy()
    if not targets.empty:
        rank_map = priority.set_index(["genus", "axis"])["priority_rank"]
        targets["priority_rank"] = [
            rank_map.loc[(g, a)] for g, a in zip(targets["genus"], targets["axis"])
        ]
        targets["acquisition_round"] = "near_threshold_genus_round_1"
        targets = targets.sort_values(
            ["priority_rank", "accepted_species"]
        ).reset_index(drop=True)

    summary = {
        "contract": "targeted_acquisition_priority_v1",
        "min_species": min_species,
        "near_threshold_genus_axis_groups": int(len(priority)),
        "target_species_axis_tasks": int(len(targets)),
        "potential_unresolved_species_axis_coverage": int(
            priority["unresolved_species_in_genus_axis"].sum()
        ) if not priority.empty else 0,
        "axes": {
            axis: {
                "genus_axis_groups": int((priority["axis"] == axis).sum()) if not priority.empty else 0,
                "target_tasks": int((targets["axis"] == axis).sum()) if not targets.empty else 0,
            }
            for axis in AXES
        },
    }
    return priority, targets, summary


@app.command()
def main(
    evidence_csv: Path = typer.Option(..., exists=True),
    unresolved_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    min_species: int = typer.Option(3),
    max_target_species_per_genus_axis: int = typer.Option(2),
) -> None:
    evidence = pd.read_csv(evidence_csv, dtype=str).fillna("")
    unresolved = pd.read_csv(unresolved_csv, dtype=str).fillna("")
    priority, targets, summary = build_priority(
        evidence,
        unresolved,
        min_species=min_species,
        max_target_species_per_genus_axis=max_target_species_per_genus_axis,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    priority.to_csv(output_dir / "near_threshold_genus_priority.csv.gz", index=False)
    targets.to_csv(output_dir / "targeted_acquisition_tasks.csv.gz", index=False)
    (output_dir / "targeted_acquisition_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
