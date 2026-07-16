"""Leave-one-species-out validation of grounded taxonomic inference.

The completed trait matrix is rectangular, but most cells are inferred. This
module asks whether those inference tiers can predict species-direct values when
the target species is removed from all support pools. It mirrors the operational
cascade:

- genus/family count one species-direct winner per supporting species;
- a unique winner must meet both support and consensus thresholds;
- family is used only when genus evidence does not qualify;
- targets without a qualified taxonomic prediction remain unresolved.

Validation is descriptive and policy-auditing. It never promotes an inference
tier to direct evidence.
"""

from __future__ import annotations

import glob
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_yaml(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter(f"YAML config must be a mapping: {path}")
    return payload


def load_species_direct_cells(cascade_root: Path, traits: set[str]) -> pd.DataFrame:
    files = sorted(glob.glob(str(cascade_root / "shards" / "shard_*" / "fills.csv.gz")))
    if not files:
        raise typer.BadParameter(f"no fill shards found under {cascade_root}")
    parts: list[pd.DataFrame] = []
    usecols = [
        "accepted_species",
        "genus",
        "family",
        "trait_name",
        "filled_value",
        "fill_tier",
        "support_n",
        "value_distribution",
    ]
    for file_path in files:
        frame = pd.read_csv(file_path, usecols=usecols, dtype=str).fillna("")
        frame = frame.loc[
            frame["fill_tier"].eq("species_direct")
            & frame["trait_name"].isin(traits)
        ].copy()
        if not frame.empty:
            parts.append(frame)
    if not parts:
        raise typer.BadParameter("no species-direct cells found for requested traits")
    direct = pd.concat(parts, ignore_index=True)
    if direct.duplicated(["accepted_species", "trait_name"]).any():
        raise typer.BadParameter("species-direct matrix contains duplicate species x trait cells")
    return direct.sort_values(["trait_name", "accepted_species"]).reset_index(drop=True)


def _parse_distribution(value: str, fallback_value: str) -> dict[str, float]:
    try:
        payload = json.loads(str(value)) if str(value).strip() else {}
    except json.JSONDecodeError:
        payload = {}
    distribution: dict[str, float] = {}
    if isinstance(payload, dict):
        for key, raw in payload.items():
            try:
                numeric = float(raw)
            except (TypeError, ValueError):
                continue
            if np.isfinite(numeric) and numeric > 0:
                distribution[str(key)] = numeric
    if not distribution:
        distribution[str(fallback_value)] = 1.0
    return distribution


def prepare_direct_cells(direct: pd.DataFrame) -> pd.DataFrame:
    result = direct.copy()
    result["raw_distribution"] = [
        _parse_distribution(value, filled)
        for value, filled in zip(
            result["value_distribution"], result["filled_value"], strict=True
        )
    ]
    result["winner_weight"] = [
        float(distribution.get(str(filled), 0.0))
        for distribution, filled in zip(
            result["raw_distribution"], result["filled_value"], strict=True
        )
    ]
    if (result["winner_weight"] <= 0).any():
        raise typer.BadParameter("a species-direct winner has no positive reconstructed weight")
    return result


def _weighted_winner(distribution: dict[str, float]) -> tuple[str, float, float] | None:
    positive = {str(key): float(value) for key, value in distribution.items() if value > 0}
    if not positive:
        return None
    winner = sorted(positive, key=lambda value: (-positive[value], value))[0]
    total = float(sum(positive.values()))
    return winner, float(positive[winner]), float(positive[winner] / total)


def _winner_distribution(frame: pd.DataFrame) -> dict[str, float]:
    totals: dict[str, float] = {}
    for row in frame.itertuples(index=False):
        value = str(row.filled_value)
        totals[value] = totals.get(value, 0.0) + 1.0
    return totals


def _raw_distribution(frame: pd.DataFrame) -> dict[str, float]:
    totals: dict[str, float] = {}
    for distribution in frame["raw_distribution"]:
        for value, weight in distribution.items():
            totals[str(value)] = totals.get(str(value), 0.0) + float(weight)
    return totals


def _prediction_record(
    distribution: dict[str, float],
    prefix: str,
) -> dict[str, Any]:
    winner = _weighted_winner(distribution)
    if winner is None:
        return {
            f"{prefix}_prediction": "",
            f"{prefix}_support_species": 0,
            f"{prefix}_total_support_species": 0,
            f"{prefix}_winner_share": np.nan,
            f"{prefix}_top_tie_n": 0,
            f"{prefix}_distribution": "{}",
        }
    value, winner_weight, share = winner
    top_weight = max(distribution.values())
    top_tie_n = sum(float(weight) == float(top_weight) for weight in distribution.values())
    return {
        f"{prefix}_prediction": value,
        f"{prefix}_support_species": int(winner_weight),
        f"{prefix}_total_support_species": int(sum(distribution.values())),
        f"{prefix}_winner_share": float(share),
        f"{prefix}_top_tie_n": int(top_tie_n),
        f"{prefix}_distribution": json.dumps(
            {key: round(float(weight), 6) for key, weight in sorted(distribution.items())},
            sort_keys=True,
        ),
    }


def leave_one_species_out_predictions(
    direct: pd.DataFrame,
    genus_min_support: int,
    family_min_support: int,
    genus_min_consensus: float = 1.0,
    family_min_consensus: float = 0.8,
    family_min_supporting_genera: int = 2,
    inference_policies: dict[str, dict[str, Any]] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for trait_name, trait in direct.groupby("trait_name", sort=True):
        policy = (inference_policies or {}).get(str(trait_name)) or {}
        trait_genus_min = int(policy.get("min_genus_support", genus_min_support))
        trait_family_min = int(policy.get("min_family_support", family_min_support))
        trait_genus_consensus = float(
            policy.get("min_genus_consensus", genus_min_consensus)
        )
        trait_family_consensus = float(
            policy.get("min_family_consensus", family_min_consensus)
        )
        trait_family_genera_min = int(
            policy.get(
                "min_family_supporting_genera",
                family_min_supporting_genera,
            )
        )
        family_allowed = bool(policy.get("family_allowed", True))
        trait = trait.reset_index(drop=True)
        for target_index, target in trait.iterrows():
            remaining = trait.drop(index=target_index)
            genus_name = str(target["genus"])
            family_name = str(target["family"])
            genus = remaining.loc[
                remaining["genus"].astype(str).eq(genus_name) & pd.Series(genus_name != "", index=remaining.index)
            ]
            family = remaining.loc[
                remaining["family"].astype(str).eq(family_name) & pd.Series(family_name != "", index=remaining.index)
            ]
            genus_distribution = _winner_distribution(genus)
            family_distribution = _winner_distribution(family)
            record: dict[str, Any] = {
                "accepted_species": str(target["accepted_species"]),
                "genus": genus_name,
                "family": family_name,
                "trait_name": str(trait_name),
                "actual_value": str(target["filled_value"]),
                "direct_support_n": int(float(target["support_n"])) if str(target["support_n"]).strip() else 0,
                "operational_genus_min_support": trait_genus_min,
                "operational_family_min_support": trait_family_min,
                "operational_genus_min_consensus": trait_genus_consensus,
                "operational_family_min_consensus": trait_family_consensus,
                "operational_family_min_supporting_genera": trait_family_genera_min,
                "operational_family_allowed": family_allowed,
            }
            record.update(
                _prediction_record(genus_distribution, "genus")
            )
            record.update(
                _prediction_record(family_distribution, "family")
            )
            family_winner = str(record["family_prediction"])
            record["family_supporting_genera"] = int(
                family.loc[family["filled_value"].astype(str).eq(family_winner), "genus"]
                .astype(str)
                .loc[lambda values: values.ne("")]
                .nunique()
            )

            if (
                record["genus_prediction"]
                and int(record["genus_support_species"]) >= trait_genus_min
                and int(record["genus_top_tie_n"]) == 1
                and float(record["genus_winner_share"]) >= trait_genus_consensus
            ):
                cascade_tier = "genus_inference"
                prediction = record["genus_prediction"]
                support = record["genus_support_species"]
                winner_share = record["genus_winner_share"]
                distribution = record["genus_distribution"]
            elif (
                family_allowed
                and record["family_prediction"]
                and int(record["family_support_species"]) >= trait_family_min
                and int(record["family_top_tie_n"]) == 1
                and float(record["family_winner_share"]) >= trait_family_consensus
                and int(record["family_supporting_genera"])
                >= trait_family_genera_min
            ):
                cascade_tier = "family_inference"
                prediction = record["family_prediction"]
                support = record["family_support_species"]
                winner_share = record["family_winner_share"]
                distribution = record["family_distribution"]
            else:
                cascade_tier = "unresolved_no_evidence"
                prediction = ""
                support = 0
                winner_share = np.nan
                distribution = "{}"

            record.update(
                {
                    "cascade_prediction_tier": cascade_tier,
                    "cascade_prediction": prediction,
                    "cascade_support_species": int(support),
                    "cascade_winner_share": winner_share,
                    "cascade_distribution": distribution,
                    "cascade_correct": bool(prediction == record["actual_value"]),
                }
            )
            rows.append(record)
    return pd.DataFrame(rows)


def _classification_metrics(frame: pd.DataFrame, actual: str, predicted: str) -> dict[str, Any]:
    work = frame.loc[
        frame[actual].astype(str).ne("") & frame[predicted].astype(str).ne("")
    ].copy()
    if work.empty:
        return {
            "n": 0,
            "accuracy": np.nan,
            "balanced_accuracy": np.nan,
            "macro_f1": np.nan,
            "n_classes_actual": 0,
        }
    labels = sorted(set(work[actual].astype(str)) | set(work[predicted].astype(str)))
    recalls: list[float] = []
    f1_values: list[float] = []
    for label in labels:
        actual_positive = work[actual].astype(str).eq(label)
        predicted_positive = work[predicted].astype(str).eq(label)
        tp = int((actual_positive & predicted_positive).sum())
        fn = int((actual_positive & ~predicted_positive).sum())
        fp = int((~actual_positive & predicted_positive).sum())
        recall = float(tp / (tp + fn)) if tp + fn else np.nan
        precision = float(tp / (tp + fp)) if tp + fp else np.nan
        if np.isfinite(recall):
            recalls.append(recall)
        if np.isfinite(recall) and np.isfinite(precision) and recall + precision > 0:
            f1_values.append(float(2 * recall * precision / (recall + precision)))
        elif np.isfinite(recall) and np.isfinite(precision):
            f1_values.append(0.0)
    return {
        "n": int(len(work)),
        "accuracy": float((work[actual].astype(str) == work[predicted].astype(str)).mean()),
        "balanced_accuracy": float(np.mean(recalls)) if recalls else np.nan,
        "macro_f1": float(np.mean(f1_values)) if f1_values else np.nan,
        "n_classes_actual": int(work[actual].astype(str).nunique()),
    }


def summarize_prediction_layers(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for trait_name, trait in predictions.groupby("trait_name", sort=True):
        for layer, column in (
            ("genus_prediction_all_evaluable", "genus_prediction"),
            ("family_prediction_all_evaluable", "family_prediction"),
            ("cascade_policy_prediction", "cascade_prediction"),
        ):
            metrics = _classification_metrics(trait, "actual_value", column)
            rows.append(
                {
                    "trait_name": trait_name,
                    "validation_layer": layer,
                    **metrics,
                }
            )
        for cascade_tier, subset in trait.groupby("cascade_prediction_tier", sort=True):
            metrics = _classification_metrics(subset, "actual_value", "cascade_prediction")
            rows.append(
                {
                    "trait_name": trait_name,
                    "validation_layer": f"cascade_used_{cascade_tier}",
                    **metrics,
                }
            )
    return pd.DataFrame(rows)


def threshold_grid(
    direct: pd.DataFrame,
    genus_thresholds: list[int],
    family_thresholds: list[int],
    genus_min_consensus: float = 1.0,
    family_min_consensus: float = 0.8,
    family_min_supporting_genera: int = 2,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for genus_min in genus_thresholds:
        for family_min in family_thresholds:
            predictions = leave_one_species_out_predictions(
                direct,
                genus_min_support=genus_min,
                family_min_support=family_min,
                genus_min_consensus=genus_min_consensus,
                family_min_consensus=family_min_consensus,
                family_min_supporting_genera=family_min_supporting_genera,
            )
            for trait_name, trait in predictions.groupby("trait_name", sort=True):
                metrics = _classification_metrics(
                    trait, "actual_value", "cascade_prediction"
                )
                tier_counts = trait["cascade_prediction_tier"].value_counts().to_dict()
                rows.append(
                    {
                        "trait_name": trait_name,
                        "genus_min_support": int(genus_min),
                        "family_min_support": int(family_min),
                        **metrics,
                        "n_genus_predictions": int(
                            tier_counts.get("genus_inference", 0)
                        ),
                        "n_family_predictions": int(
                            tier_counts.get("family_inference", 0)
                        ),
                        "n_unresolved_predictions": int(
                            tier_counts.get("unresolved_no_evidence", 0)
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _support_bin(value: int, bins: list[list[int]]) -> str:
    for lower, upper in bins:
        if int(lower) <= value <= int(upper):
            return f"{int(lower)}-{int(upper)}" if lower != upper else str(int(lower))
    return "outside_grid"


def support_stratified_metrics(
    predictions: pd.DataFrame,
    bins: list[list[int]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for trait_name, trait in predictions.groupby("trait_name", sort=True):
        for layer in ("genus", "family"):
            pred_col = f"{layer}_prediction"
            support_col = f"{layer}_support_species"
            work = trait.loc[trait[pred_col].astype(str).ne("")].copy()
            if work.empty:
                continue
            work["support_bin"] = work[support_col].map(
                lambda value: _support_bin(int(value), bins)
            )
            for support_bin, group in work.groupby("support_bin", sort=False):
                metrics = _classification_metrics(group, "actual_value", pred_col)
                rows.append(
                    {
                        "trait_name": trait_name,
                        "inference_layer": layer,
                        "support_bin": support_bin,
                        **metrics,
                        "median_winner_share": float(
                            pd.to_numeric(
                                group[f"{layer}_winner_share"], errors="coerce"
                            ).median()
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _distribution_probability(distribution_json: str, positive_value: str) -> float:
    try:
        distribution = json.loads(str(distribution_json))
    except json.JSONDecodeError:
        return np.nan
    if not isinstance(distribution, dict):
        return np.nan
    weights = {
        str(key): float(value)
        for key, value in distribution.items()
        if float(value) > 0
    }
    total = float(sum(weights.values()))
    return float(weights.get(positive_value, 0.0) / total) if total > 0 else np.nan


def focused_binary_metrics(
    predictions: pd.DataFrame,
    trait_name: str,
    positive_value: str,
) -> pd.DataFrame:
    trait = predictions.loc[predictions["trait_name"].eq(trait_name)].copy()
    rows: list[dict[str, Any]] = []
    for layer, pred_col, dist_col in (
        ("genus_all_evaluable", "genus_prediction", "genus_distribution"),
        ("family_all_evaluable", "family_prediction", "family_distribution"),
        ("cascade_policy", "cascade_prediction", "cascade_distribution"),
    ):
        work = trait.loc[trait[pred_col].astype(str).ne("")].copy()
        actual_positive = work["actual_value"].astype(str).eq(positive_value)
        predicted_positive = work[pred_col].astype(str).eq(positive_value)
        tp = int((actual_positive & predicted_positive).sum())
        tn = int((~actual_positive & ~predicted_positive).sum())
        fp = int((~actual_positive & predicted_positive).sum())
        fn = int((actual_positive & ~predicted_positive).sum())
        probability = work[dist_col].map(
            lambda value: _distribution_probability(value, positive_value)
        )
        valid_probability = probability.notna()
        y = actual_positive.astype(float)
        rows.append(
            {
                "trait_name": trait_name,
                "positive_value": positive_value,
                "validation_layer": layer,
                "n": int(len(work)),
                "accuracy": float((tp + tn) / len(work)) if len(work) else np.nan,
                "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
                "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
                "positive_predictive_value": float(tp / (tp + fp)) if tp + fp else np.nan,
                "negative_predictive_value": float(tn / (tn + fn)) if tn + fn else np.nan,
                "actual_positive_share": float(actual_positive.mean()) if len(work) else np.nan,
                "predicted_positive_share": float(predicted_positive.mean()) if len(work) else np.nan,
                "probabilistic_brier": float(
                    np.mean((y.loc[valid_probability] - probability.loc[valid_probability]) ** 2)
                )
                if valid_probability.any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    cascade_root: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    validation_config: Path = typer.Option(
        Path("config/trait_inference_validation.yml"), exists=True
    ),
    cascade_config: Path = typer.Option(
        Path("config/trait_fill_cascade.yml"), exists=True
    ),
) -> None:
    validation = _load_yaml(validation_config)
    cascade = _load_yaml(cascade_config)
    traits = {str(value) for value in validation["validate_traits"]}
    direct = prepare_direct_cells(load_species_direct_cells(cascade_root, traits))
    genus_min = int(cascade.get("min_genus_support", 1))
    family_min = int(cascade.get("min_family_support", 3))
    family_genera_min = int(cascade.get("min_family_supporting_genera", 2))
    genus_consensus = float(cascade.get("min_genus_consensus", 1.0))
    family_consensus = float(cascade.get("min_family_consensus", 0.8))
    inference_policies = dict(cascade.get("inference_policies") or {})
    predictions = leave_one_species_out_predictions(
        direct,
        genus_min_support=genus_min,
        family_min_support=family_min,
        genus_min_consensus=genus_consensus,
        family_min_consensus=family_consensus,
        family_min_supporting_genera=family_genera_min,
        inference_policies=inference_policies,
    )
    layer_metrics = summarize_prediction_layers(predictions)
    grid = threshold_grid(
        direct,
        [int(value) for value in validation["cascade_policy"]["genus_support_grid"]],
        [int(value) for value in validation["cascade_policy"]["family_support_grid"]],
        genus_min_consensus=genus_consensus,
        family_min_consensus=family_consensus,
        family_min_supporting_genera=family_genera_min,
    )
    support_metrics = support_stratified_metrics(
        predictions,
        [[int(value) for value in pair] for pair in validation["reporting"]["support_bins"]],
    )
    focused = validation["focused_binary_diagnostic"]
    binary = focused_binary_metrics(
        predictions,
        str(focused["trait_name"]),
        str(focused["positive_value"]),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    direct.drop(columns=["raw_distribution"]).to_csv(
        output_dir / "species_direct_cells.csv", index=False
    )
    predictions.to_csv(output_dir / "trait_inference_loo_predictions.csv", index=False)
    layer_metrics.to_csv(output_dir / "trait_inference_loo_metrics.csv", index=False)
    grid.to_csv(output_dir / "trait_inference_threshold_grid.csv", index=False)
    support_metrics.to_csv(output_dir / "trait_inference_support_metrics.csv", index=False)
    binary.to_csv(output_dir / "focused_binary_validation.csv", index=False)

    manifest = {
        "contract": str(validation["contract"]),
        "n_species_direct_cells": int(len(direct)),
        "n_direct_species": int(direct["accepted_species"].nunique()),
        "trait_direct_counts": {
            str(key): int(value)
            for key, value in direct["trait_name"].value_counts().sort_index().items()
        },
        "operational_cascade_thresholds": {
            "min_genus_support": genus_min,
            "min_family_support": family_min,
            "min_family_supporting_genera": family_genera_min,
            "min_genus_consensus": genus_consensus,
            "min_family_consensus": family_consensus,
            "trait_specific_policy_count": len(inference_policies),
        },
        "operational_layer_metrics": layer_metrics.to_dict("records"),
        "focused_binary_metrics": binary.to_dict("records"),
        "threshold_grid_rows": int(len(grid)),
        "interpretation": (
            "Leave-one-species-out validation removes each target species from all taxonomic support. "
            "Inference tiers remain inferred regardless of validation accuracy."
        ),
    }
    (output_dir / "trait_inference_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
