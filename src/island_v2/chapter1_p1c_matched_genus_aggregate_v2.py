"""Aggregate the already completed P1c matched-genus permutation shards.

This module exists only because the original aggregator compared independently
recomputed observed floating-point statistics by exact equality.  The 2000
permutations are not regenerated.  Their frozen shard artifacts from workflow
run 34936193944 are reused unchanged.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


class P1CAggregateError(ValueError):
    """Raised when frozen shard artifacts violate an aggregation invariant."""


def _load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != "chapter1_p1_assembly_depth_defense_v1":
        raise P1CAggregateError("unexpected P1 contract")
    return value


def _strict_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    mapped = series.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )
    if mapped.isna().any():
        raise P1CAggregateError("invalid boolean values in shard valid column")
    return mapped.astype(bool)


def aggregate_existing_shards(
    *,
    shard_dir: Path,
    p1_config_path: Path,
    output_dir: Path,
    source_run_id: int = 34936193944,
    observed_tolerance: float = 1.0e-12,
) -> dict[str, Any]:
    p1 = _load_config(p1_config_path)
    spec = p1["p1c_matched_complexity_genus_null"]
    expected_total = int(spec["permutations"])
    expected_shards = int(spec["execution_shards"])
    minimum_valid = int(spec["minimum_valid_permutations"])
    master_seed = int(spec["random_seed"])
    seed_rule = str(spec["permutation_seed_rule"])
    threshold = float(spec["primary_pass_rule"]["max_randomization_p_value"])

    csv_files = sorted(shard_dir.rglob("p1c_null_shard_*.csv.gz"))
    manifest_files = sorted(shard_dir.rglob("p1c_null_shard_*_manifest.json"))
    if len(csv_files) != expected_shards:
        raise P1CAggregateError(
            f"expected {expected_shards} shard CSVs, found {len(csv_files)}"
        )
    if len(manifest_files) != expected_shards:
        raise P1CAggregateError(
            f"expected {expected_shards} shard manifests, found {len(manifest_files)}"
        )

    manifests = [json.loads(path.read_text(encoding="utf-8")) for path in manifest_files]
    for item in manifests:
        if item.get("contract") != "chapter1_p1c_matched_genus_null_shard_v1":
            raise P1CAggregateError("unexpected shard manifest contract")
        if int(item.get("n_rows", -1)) != int(spec["permutations_per_shard"]):
            raise P1CAggregateError("unexpected shard row count")
        if item.get("source_artifact_digest") != p1["pinned_primary_artifact"]["artifact_digest"]:
            raise P1CAggregateError("source artifact digest differs among shards")

    shard_ids = [int(item["shard_id"]) for item in manifests]
    if sorted(shard_ids) != list(range(expected_shards)):
        raise P1CAggregateError("shard IDs do not equal frozen schedule")

    observed_values = np.asarray(
        [float(item["observed_primary_statistic"]) for item in manifests], dtype=float
    )
    if not np.isfinite(observed_values).all():
        raise P1CAggregateError("non-finite observed statistic in shard manifest")
    observed_spread = float(observed_values.max() - observed_values.min())
    if observed_spread > float(observed_tolerance):
        raise P1CAggregateError(
            f"observed statistic spread {observed_spread} exceeds tolerance {observed_tolerance}"
        )
    observed = float(np.median(observed_values))

    frame = pd.concat([pd.read_csv(path) for path in csv_files], ignore_index=True)
    if len(frame) != expected_total:
        raise P1CAggregateError(
            f"expected {expected_total} permutation rows, found {len(frame)}"
        )
    ids = pd.to_numeric(frame["permutation_id"], errors="raise").astype(int)
    if ids.duplicated().any() or set(ids) != set(range(expected_total)):
        raise P1CAggregateError("permutation IDs do not equal frozen 0..1999 schedule")
    seeds = pd.to_numeric(frame["seed"], errors="raise").astype(int)
    expected_seeds = ids + master_seed
    if not np.array_equal(seeds.to_numpy(), expected_seeds.to_numpy()):
        raise P1CAggregateError("permutation seeds do not follow frozen schedule")

    valid_mask = _strict_bool(frame["valid"])
    valid = frame.loc[valid_mask].copy()
    n_valid = int(len(valid))
    evaluable = n_valid >= minimum_valid

    randomization_p: float | None = None
    exceedances: int | None = None
    quantiles: dict[str, float] = {}
    if n_valid:
        values = pd.to_numeric(
            valid["primary_median_conditional_attenuation"], errors="raise"
        )
        if not np.isfinite(values.to_numpy(float)).all():
            raise P1CAggregateError("non-finite primary statistic in valid permutation")
        for q in (0.025, 0.25, 0.5, 0.75, 0.975):
            quantiles[str(q)] = float(values.quantile(q))
        if evaluable:
            exceedances = int(values.ge(observed).sum())
            randomization_p = float((1 + exceedances) / (1 + n_valid))

    if not evaluable:
        status = "NOT_EVALUABLE_insufficient_valid_permutations"
    elif randomization_p is not None and randomization_p <= threshold:
        status = "true_genus_exceeds_matched_complexity_null"
    else:
        status = "true_genus_not_distinguishable_from_matched_complexity_null"

    output_dir.mkdir(parents=True, exist_ok=True)
    frame.sort_values("permutation_id").to_csv(
        output_dir / "p1c_matched_genus_null_permutations.csv.gz", index=False
    )
    result = {
        "contract": "chapter1_p1c_matched_genus_null_result_v2_aggregation_fix",
        "status": status,
        "source_permutation_run_id": int(source_run_id),
        "permutations_regenerated": False,
        "permutations_requested": expected_total,
        "valid_permutations": n_valid,
        "minimum_valid_permutations": minimum_valid,
        "observed_primary_median_conditional_attenuation": observed,
        "observed_statistic_spread_across_shards": observed_spread,
        "observed_statistic_equality_tolerance": float(observed_tolerance),
        "null_primary_statistic_quantiles": quantiles,
        "null_exceedances_greater_equal_observed": exceedances,
        "one_sided_randomization_p_value": randomization_p,
        "pass_threshold": threshold,
        "strong_assembly_depth_defense_pass": bool(
            evaluable and randomization_p is not None and randomization_p <= threshold
        ),
        "seed_rule": seed_rule,
        "master_seed": master_seed,
        "primary_scope": "direct_only_Palearctic_two_axis",
        "aggregation_fix": (
            "Replace exact float equality of independently recomputed observed statistics "
            "with a predeclared numerical-equivalence tolerance of 1e-12. No permutation "
            "was regenerated or changed."
        ),
        "claim_boundary": (
            "A positive result shows genus-specific taxonomic structure beyond matched "
            "arbitrary within-family partitions. It remains non-causal, does not make the "
            "family-to-genus incremental magnitude spatially precise, and does not exclude "
            "within-lineage evolution or identify a pollinator mechanism."
        ),
    }
    (output_dir / "chapter1_p1c_matched_genus_null_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


@app.command("aggregate")
def aggregate_command(
    shard_dir: Path = typer.Option(..., exists=True, file_okay=False),
    p1_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
    source_run_id: int = typer.Option(34936193944),
    observed_tolerance: float = typer.Option(1.0e-12),
) -> None:
    typer.echo(
        json.dumps(
            aggregate_existing_shards(
                shard_dir=shard_dir,
                p1_config_path=p1_config_path,
                output_dir=output_dir,
                source_run_id=source_run_id,
                observed_tolerance=observed_tolerance,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
