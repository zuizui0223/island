"""Prioritize target-aligned direct-trait acquisition near validated-low eligibility."""
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
REQUIRED_DOMINANCE = {
    "flower_colour": 0.90,
    "floral_structural_complexity": 0.80,
    "reproductive_assurance": 0.95,
}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}
VALUE_COLUMNS = (
    "reported_value", "raw_candidate_value", "normalized_value",
    "provisional_candidate_value", "candidate_value", "value",
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _genus(species: object) -> str:
    text = _text(species)
    return text.split(" ", 1)[0] if text else ""


def _value(frame: pd.DataFrame) -> pd.Series:
    out = pd.Series("", index=frame.index, dtype="object")
    for column in VALUE_COLUMNS:
        if column in frame.columns:
            values = frame[column].map(_text)
            out = out.where(out.ne(""), values)
    return out


def _normalise_hints(candidate_hints: pd.DataFrame | None) -> pd.DataFrame:
    if candidate_hints is None or candidate_hints.empty:
        return pd.DataFrame(columns=["accepted_species", "trait_name", "candidate_value"])
    hints = candidate_hints.copy()
    aliases = {"species": "accepted_species", "trait": "trait_name", "field": "trait_name"}
    hints = hints.rename(
        columns={a: b for a, b in aliases.items() if a in hints.columns and b not in hints.columns}
    )
    for column in ("accepted_species", "trait_name"):
        if column not in hints.columns:
            hints[column] = ""
        hints[column] = hints[column].map(_text)
    hints["candidate_value"] = _value(hints)
    return hints.loc[
        hints["accepted_species"].ne("") & hints["trait_name"].ne(""),
        ["accepted_species", "trait_name", "candidate_value"],
    ].drop_duplicates()


def build_priority(
    evidence: pd.DataFrame,
    unresolved: pd.DataFrame,
    min_species: int = 3,
    max_target_species_per_genus_axis: int = 2,
    candidate_hints: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    evidence = evidence.copy()
    unresolved = unresolved.copy()
    for column in ("accepted_species", "trait_name", "evidence_quality", "evidence_scope"):
        if column not in evidence.columns:
            evidence[column] = ""
    for column in ("accepted_species", "axis"):
        if column not in unresolved.columns:
            unresolved[column] = ""

    evidence["accepted_species"] = evidence["accepted_species"].map(_text)
    evidence["trait_name"] = evidence["trait_name"].map(_text)
    evidence["value"] = _value(evidence)
    evidence["genus"] = evidence["accepted_species"].map(_genus)
    unresolved["accepted_species"] = unresolved["accepted_species"].map(_text)
    unresolved["axis"] = unresolved["axis"].map(_text)
    unresolved["genus"] = unresolved["accepted_species"].map(_genus)
    hints = _normalise_hints(candidate_hints)

    direct = evidence.loc[
        evidence["evidence_quality"].isin({"high", "medium"})
        & evidence["evidence_scope"].isin(DIRECT_SCOPES)
        & evidence["genus"].ne("")
        & evidence["value"].ne("")
    ].drop_duplicates(["accepted_species", "trait_name", "value"]).copy()

    priority_rows: list[dict[str, object]] = []
    target_frames: list[pd.DataFrame] = []
    rejected_no_hint = 0
    rejected_dominance = 0

    for axis, traits in AXES.items():
        axis_unresolved = unresolved.loc[
            unresolved["axis"].eq(axis) & unresolved["genus"].ne("")
        ].drop_duplicates(["accepted_species", "axis"]).copy()
        axis_direct = direct.loc[direct["trait_name"].isin(traits)].copy()
        if axis_unresolved.empty or axis_direct.empty:
            continue

        votes = (
            axis_direct.groupby(["genus", "trait_name", "value"])["accepted_species"]
            .nunique().rename("value_species").reset_index()
        )
        totals = (
            axis_direct.groupby(["genus", "trait_name"])["accepted_species"]
            .nunique().rename("n_direct_species").reset_index()
        )
        rule = votes.merge(totals, on=["genus", "trait_name"])
        rule["dominance"] = rule["value_species"] / rule["n_direct_species"]
        rule = (
            rule.sort_values(
                ["genus", "trait_name", "value_species", "value"],
                ascending=[True, True, False, True],
            ).drop_duplicates(["genus", "trait_name"])
        )
        best = (
            rule.sort_values(
                ["genus", "n_direct_species", "dominance", "trait_name"],
                ascending=[True, False, False, True],
            ).drop_duplicates("genus")
            .rename(columns={"trait_name": "best_supported_trait", "value": "target_value"})
        )
        unresolved_counts = (
            axis_unresolved.groupby("genus")["accepted_species"]
            .nunique().rename("unresolved_species_in_genus_axis").reset_index()
        )
        near = best.merge(unresolved_counts, on="genus", how="inner")
        near = near.loc[near["n_direct_species"].between(1, min_species - 1)].copy()
        if near.empty:
            continue
        near["axis"] = axis
        near["additional_species_needed"] = min_species - near["n_direct_species"]
        near["projected_dominance_if_supportive"] = (
            near["value_species"] + near["additional_species_needed"]
        ) / (near["n_direct_species"] + near["additional_species_needed"])
        near["required_dominance"] = REQUIRED_DOMINANCE[axis]
        bad = near["projected_dominance_if_supportive"].lt(near["required_dominance"])
        rejected_dominance += int(bad.sum())
        near = near.loc[~bad].copy()

        direct_species = (
            axis_direct.groupby(["genus", "trait_name"])["accepted_species"].apply(set).to_dict()
        )
        for item in near.to_dict("records"):
            genus = str(item["genus"])
            trait = str(item["best_supported_trait"])
            target_value = str(item["target_value"])
            needed = min(
                int(item["additional_species_needed"]), max_target_species_per_genus_axis
            )
            candidates = axis_unresolved.loc[
                axis_unresolved["genus"].eq(genus)
                & ~axis_unresolved["accepted_species"].isin(
                    direct_species.get((genus, trait), set())
                )
            ].copy()
            if not hints.empty:
                h = hints.loc[hints["trait_name"].eq(trait)].copy()
                h["hint_score"] = 1 + h["candidate_value"].eq(target_value).astype(int)
                hs = h.groupby("accepted_species")["hint_score"].max()
                candidates["hint_score"] = (
                    candidates["accepted_species"].map(hs).fillna(0).astype(int)
                )
                candidates = candidates.loc[candidates["hint_score"].gt(0)]
            else:
                candidates["hint_score"] = 0
            candidates = candidates.sort_values(
                ["hint_score", "accepted_species"], ascending=[False, True]
            )
            chosen = candidates.head(needed).copy()
            if len(chosen) < needed:
                rejected_no_hint += 1
                continue
            item["leverage_score"] = float(item["unresolved_species_in_genus_axis"]) / needed
            item["target_traits"] = trait
            item["target_value"] = target_value
            item["candidate_hint_required"] = not hints.empty
            for key, value in item.items():
                chosen[key] = value
            priority_rows.append(item)
            target_frames.append(chosen)

    priority = pd.DataFrame(priority_rows)
    if not priority.empty:
        priority = priority.drop_duplicates(["genus", "axis"])
        priority = priority.sort_values(
            ["leverage_score", "unresolved_species_in_genus_axis", "genus", "axis"],
            ascending=[False, False, True, True],
        ).reset_index(drop=True)
        priority.insert(0, "priority_rank", range(1, len(priority) + 1))
    targets = (
        pd.concat(target_frames, ignore_index=True)
        if target_frames else unresolved.iloc[0:0].copy()
    )
    if not targets.empty:
        rank_map = priority.set_index(["genus", "axis"])["priority_rank"]
        targets["priority_rank"] = [
            rank_map.loc[(g, a)] for g, a in zip(targets["genus"], targets["axis"])
        ]
        targets["acquisition_round"] = "target_aligned_genus_round_2"
        targets = targets.sort_values(
            ["priority_rank", "hint_score", "accepted_species"],
            ascending=[True, False, True],
        ).reset_index(drop=True)

    summary = {
        "contract": "targeted_acquisition_priority_v2",
        "min_species": min_species,
        "candidate_hints_used": not hints.empty,
        "near_threshold_genus_axis_groups": int(len(priority)),
        "target_species_axis_tasks": int(len(targets)),
        "potential_unresolved_species_axis_coverage": int(
            priority["unresolved_species_in_genus_axis"].sum()
        ) if not priority.empty else 0,
        "rejected_groups_without_target_aligned_candidates": rejected_no_hint,
        "rejected_groups_unable_to_meet_dominance": rejected_dominance,
        "axes": {
            axis: {
                "genus_axis_groups": int((priority["axis"] == axis).sum())
                if not priority.empty else 0,
                "target_tasks": int((targets["axis"] == axis).sum())
                if not targets.empty else 0,
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
    candidate_csv: Path | None = typer.Option(None),
) -> None:
    evidence = pd.read_csv(evidence_csv, dtype=str).fillna("")
    unresolved = pd.read_csv(unresolved_csv, dtype=str).fillna("")
    hints = pd.read_csv(candidate_csv, dtype=str).fillna("") if candidate_csv else None
    priority, targets, summary = build_priority(
        evidence,
        unresolved,
        min_species=min_species,
        max_target_species_per_genus_axis=max_target_species_per_genus_axis,
        candidate_hints=hints,
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
