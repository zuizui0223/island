"""Conservative genus-level trait inference with masked validation.

Creates traceable low-quality evidence only when direct high/medium species
records show strong within-genus consistency. Family inference and global
fallback are never used.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
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
RULE_COLUMNS = [
    "genus", "axis", "trait_name", "inferred_value", "n_direct_species",
    "dominant_species", "counterexample_species", "dominance",
    "required_dominance", "eligible_before_validation", "masked_n",
    "masked_correct", "masked_accuracy", "eligible",
]
ACCEPTED_COLUMNS = [
    "accepted_species", "genus", "axis", "trait_name", "normalized_value",
    "evidence_quality", "evidence_scope", "inference_method",
    "n_direct_species", "dominant_species", "counterexample_species",
    "dominance", "masked_n", "masked_accuracy", "family_inference_used",
    "global_fallback_used",
]


@dataclass(frozen=True)
class Thresholds:
    min_species: int = 3
    min_dominance: float = 0.80
    min_masked_accuracy: float = 0.75
    colour_min_dominance: float = 0.90
    reproductive_min_dominance: float = 0.95


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _empty_rules() -> pd.DataFrame:
    return pd.DataFrame(columns=RULE_COLUMNS)


def _empty_accepted() -> pd.DataFrame:
    return pd.DataFrame(columns=ACCEPTED_COLUMNS)


def _prepare(evidence: pd.DataFrame) -> pd.DataFrame:
    frame = evidence.copy()
    for column in ("accepted_species", "trait_name", "evidence_quality", "evidence_scope"):
        if column not in frame.columns:
            frame[column] = ""

    # Real acquisition artifacts use raw_candidate_value. Unit-level normalized
    # inputs use normalized_value. Prefer normalized values, then fall back to
    # the actual ledger columns without silently discarding all evidence.
    if "normalized_value" not in frame.columns:
        frame["normalized_value"] = ""
    for fallback in ("raw_candidate_value", "reported_value", "provisional_candidate_value"):
        if fallback in frame.columns:
            missing = frame["normalized_value"].map(_text).eq("")
            frame.loc[missing, "normalized_value"] = frame.loc[missing, fallback]

    for column in frame.columns:
        frame[column] = frame[column].map(_text)
    frame = frame.loc[
        frame["evidence_quality"].isin({"high", "medium"})
        & frame["evidence_scope"].isin(DIRECT_SCOPES)
        & frame["accepted_species"].str.contains(" ", regex=False)
        & frame["normalized_value"].ne("")
    ].copy()
    frame["genus"] = frame["accepted_species"].str.split().str[0]
    frame["axis"] = ""
    for axis, traits in AXES.items():
        frame.loc[frame["trait_name"].isin(traits), "axis"] = axis
    return frame.loc[frame["axis"].ne("")].copy()


def build_species_votes(evidence: pd.DataFrame) -> pd.DataFrame:
    frame = _prepare(evidence)
    if frame.empty:
        return frame
    frame["_quality_rank"] = frame["evidence_quality"].map({"medium": 1, "high": 2}).fillna(0)
    frame = frame.sort_values(
        ["accepted_species", "trait_name", "normalized_value", "_quality_rank"],
        ascending=[True, True, True, False],
    )
    return frame.drop_duplicates(
        ["accepted_species", "trait_name", "normalized_value"], keep="first"
    ).drop(columns="_quality_rank")


def _required_dominance(axis: str, thresholds: Thresholds) -> float:
    if axis == "flower_colour":
        return thresholds.colour_min_dominance
    if axis == "reproductive_assurance":
        return thresholds.reproductive_min_dominance
    return thresholds.min_dominance


def build_genus_consensus(votes: pd.DataFrame, thresholds: Thresholds = Thresholds()) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if votes.empty:
        return _empty_rules().drop(columns=["masked_n", "masked_correct", "masked_accuracy", "eligible"])
    species_values = votes.drop_duplicates(["accepted_species", "trait_name", "normalized_value"])
    for (genus, axis, trait), group in species_values.groupby(["genus", "axis", "trait_name"]):
        per_species = group.groupby("accepted_species")["normalized_value"].agg(lambda x: sorted(set(x)))
        unambiguous = per_species.loc[per_species.map(len).eq(1)].map(lambda x: x[0])
        n_species = int(len(unambiguous))
        if n_species == 0:
            continue
        counts = unambiguous.value_counts()
        dominant_value = str(counts.index[0])
        dominant_count = int(counts.iloc[0])
        dominance = dominant_count / n_species
        required = _required_dominance(axis, thresholds)
        rows.append({
            "genus": genus, "axis": axis, "trait_name": trait,
            "inferred_value": dominant_value, "n_direct_species": n_species,
            "dominant_species": dominant_count,
            "counterexample_species": n_species - dominant_count,
            "dominance": dominance, "required_dominance": required,
            "eligible_before_validation": n_species >= thresholds.min_species and dominance >= required,
        })
    columns = RULE_COLUMNS[:10]
    return pd.DataFrame(rows, columns=columns)


def masked_validate(votes: pd.DataFrame, consensus: pd.DataFrame) -> pd.DataFrame:
    columns = ["genus", "axis", "trait_name", "masked_n", "masked_correct", "masked_accuracy"]
    if votes.empty or consensus.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, object]] = []
    for _, rule in consensus.iterrows():
        subset = votes.loc[
            votes["genus"].eq(rule["genus"])
            & votes["trait_name"].eq(rule["trait_name"])
        ]
        per_species = subset.groupby("accepted_species")["normalized_value"].agg(lambda x: sorted(set(x)))
        per_species = per_species.loc[per_species.map(len).eq(1)].map(lambda x: x[0])
        outcomes: list[bool] = []
        for species, truth in per_species.items():
            training = per_species.drop(index=species)
            if training.empty:
                continue
            outcomes.append(str(training.value_counts().index[0]) == truth)
        rows.append({
            "genus": rule["genus"], "axis": rule["axis"], "trait_name": rule["trait_name"],
            "masked_n": len(outcomes), "masked_correct": int(sum(outcomes)),
            "masked_accuracy": (sum(outcomes) / len(outcomes)) if outcomes else 0.0,
        })
    return pd.DataFrame(rows, columns=columns)


def infer_low_evidence(
    unresolved: pd.DataFrame,
    evidence: pd.DataFrame,
    thresholds: Thresholds = Thresholds(),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tasks = unresolved.copy()
    for column in ("accepted_species", "axis"):
        if column not in tasks.columns:
            tasks[column] = ""
        tasks[column] = tasks[column].map(_text)
    tasks["genus"] = tasks["accepted_species"].str.split().str[0]

    votes = build_species_votes(evidence)
    consensus = build_genus_consensus(votes, thresholds)
    if consensus.empty:
        remaining = tasks.copy()
        remaining["state"] = "unresolved"
        remaining["reason"] = "no_usable_direct_evidence"
        return _empty_accepted(), remaining, _empty_rules()

    validation = masked_validate(votes, consensus)
    rules = consensus.merge(validation, on=["genus", "axis", "trait_name"], how="left")
    rules["masked_n"] = rules["masked_n"].fillna(0).astype(int)
    rules["masked_correct"] = rules["masked_correct"].fillna(0).astype(int)
    rules["masked_accuracy"] = rules["masked_accuracy"].fillna(0.0)
    rules["eligible"] = rules["eligible_before_validation"] & rules["masked_accuracy"].ge(thresholds.min_masked_accuracy)
    rules = rules.reindex(columns=RULE_COLUMNS)

    accepted_parts: list[pd.DataFrame] = []
    for axis in AXES:
        axis_tasks = tasks.loc[tasks["axis"].eq(axis)]
        axis_rules = rules.loc[rules["axis"].eq(axis) & rules["eligible"]]
        if axis_tasks.empty or axis_rules.empty:
            continue
        merged = axis_tasks.merge(axis_rules, on=["genus", "axis"], how="inner")
        merged["normalized_value"] = merged["inferred_value"]
        merged["evidence_quality"] = "low"
        merged["evidence_scope"] = "genus_consensus"
        merged["inference_method"] = "validated_genus_consensus"
        merged["family_inference_used"] = False
        merged["global_fallback_used"] = False
        accepted_parts.append(merged.reindex(columns=ACCEPTED_COLUMNS))

    accepted = pd.concat(accepted_parts, ignore_index=True) if accepted_parts else _empty_accepted()
    resolved_keys = set(zip(accepted["accepted_species"], accepted["axis"]))
    remaining = tasks.loc[
        [(species, axis) not in resolved_keys for species, axis in zip(tasks["accepted_species"], tasks["axis"])]
    ].copy()
    remaining["state"] = "unresolved"
    remaining["reason"] = "no_validated_genus_consensus"
    return accepted, remaining, rules


@app.command()
def build(
    unresolved_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    evidence_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
    min_species: int = 3,
    min_dominance: float = 0.80,
    min_masked_accuracy: float = 0.75,
) -> None:
    unresolved = pd.read_csv(unresolved_csv, dtype=str).fillna("")
    evidence = pd.read_csv(evidence_csv, dtype=str).fillna("")
    thresholds = Thresholds(
        min_species=min_species,
        min_dominance=min_dominance,
        min_masked_accuracy=min_masked_accuracy,
    )
    accepted, remaining, rules = infer_low_evidence(unresolved, evidence, thresholds)
    output_dir.mkdir(parents=True, exist_ok=True)
    accepted.to_csv(output_dir / "accepted_validated_low_evidence.csv.gz", index=False)
    remaining.to_csv(output_dir / "unresolved_after_validated_low.csv.gz", index=False)
    rules.to_csv(output_dir / "validated_genus_rules.csv.gz", index=False)
    summary = {
        "contract": "validated_low_inference_v1",
        "input_tasks": int(len(unresolved)),
        "accepted_low_records": int(len(accepted)),
        "resolved_species_axis_tasks": int(len(set(zip(accepted["accepted_species"], accepted["axis"])))),
        "remaining_tasks": int(len(remaining)),
        "eligible_rules": int(rules["eligible"].sum()) if "eligible" in rules else 0,
        "family_inference_used": False,
        "global_fallback_used": False,
    }
    (output_dir / "validated_low_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
