"""Final P1c matched-complexity genus randomization test.

The inferential contract is frozen in ``config/chapter1_p1_assembly_depth_defense.yml``.
This module runs deterministic permutation shards and aggregates them without any
adaptive retuning.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.chapter1_p1c_matched_genus_qualification import (
    fit_palearctic_profiles,
    load_yaml,
    make_matched_pseudo_taxonomy,
)
from island_v2.chapter1_pr138_lineage_representation_bridge import match_gift_species
from island_v2.chapter1_taxonomic_depth_decomposition import (
    build_source_group_contract,
    build_taxonomic_decomposition,
    file_sha256,
    validate_taxonomy,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


class P1CNullError(ValueError):
    """Raised when final P1c randomization invariants fail."""


def _load_inputs(
    artifact_root: Path,
    taxonomy_path: Path,
    p1_config_path: Path,
    pattern_config_path: Path,
    source_config_path: Path,
    explanation_config_path: Path,
) -> dict[str, Any]:
    p1 = load_yaml(p1_config_path)
    pattern = load_yaml(pattern_config_path)
    source_config = load_yaml(source_config_path)
    explanation = load_yaml(explanation_config_path)
    if p1.get("contract") != "chapter1_p1_assembly_depth_defense_v1":
        raise P1CNullError("unexpected P1 contract")
    if source_config.get("contract") != "chapter1_pr138_source_pool_sensitivity_v1":
        raise P1CNullError("unexpected source contract")
    spec = explanation["validations"]["V2_H3_taxonomic_depth"]["frozen_implementation"]
    taxonomy_hash = file_sha256(taxonomy_path)
    taxonomy = validate_taxonomy(
        pd.read_csv(taxonomy_path),
        str(spec["taxonomy_input"]["sha256_newline_canonicalized"]),
        taxonomy_hash,
    )
    return {
        "p1": p1,
        "pattern": pattern,
        "source_config": source_config,
        "spec": spec,
        "taxonomy": taxonomy,
        "taxonomy_hash": taxonomy_hash,
        "direct_scores": pd.read_csv(
            artifact_root / "syndrome/direct/species_syndrome_concordance.csv.gz"
        ),
        "status_flora": pd.read_csv(
            artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz"
        ),
        "gift_flora": pd.read_csv(
            artifact_root / "fixed/source/source/gift_native_mainland_flora.csv.gz"
        ),
        "assignments": pd.read_csv(
            artifact_root / "fixed/source/result/island_source_assignments.csv.gz"
        ),
        "covariates": pd.read_csv(
            artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv"
        ),
        "realm_assignment": pd.read_csv(
            artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv"
        ),
        "true_decomposition": pd.read_csv(
            artifact_root / "taxonomic-depth/direct_only/decomposition.csv"
        ),
    }


def _vector_norms(slopes: pd.DataFrame, axes: list[str]) -> pd.DataFrame:
    expected = set(axes)
    rows: list[dict[str, Any]] = []
    for (mode, stratum), group in slopes.groupby(["source_mode", "stratum"], sort=True):
        if set(group["syndrome"].astype(str)) != expected:
            continue
        values = group.set_index("syndrome").loc[axes, "distance_slope"].to_numpy(float)
        if not np.isfinite(values).all():
            continue
        rows.append(
            {
                "source_mode": str(mode),
                "stratum": str(stratum),
                "vector_norm": float(np.linalg.norm(values)),
            }
        )
    return pd.DataFrame(rows)


def _true_reference(data: dict[str, Any]) -> tuple[pd.DataFrame, float, dict[str, float]]:
    p1 = data["p1"]
    primary = p1["primary_context"]
    axes = [str(x) for x in primary["response_axes"]]
    modes = [str(x) for x in primary["source_modes"]]
    strata = [str(x) for x in primary["strata"]]
    family_slopes = fit_palearctic_profiles(
        data["true_decomposition"],
        data["covariates"],
        data["realm_assignment"],
        data["pattern"],
        stage="after_family_residual",
        source_modes=modes,
        strata=strata,
        axes=axes,
    )
    genus_slopes = fit_palearctic_profiles(
        data["true_decomposition"],
        data["covariates"],
        data["realm_assignment"],
        data["pattern"],
        stage="after_genus_residual",
        source_modes=modes,
        strata=strata,
        axes=axes,
    )
    family = _vector_norms(family_slopes, axes).rename(columns={"vector_norm": "family_norm"})
    genus = _vector_norms(genus_slopes, axes).rename(columns={"vector_norm": "genus_norm"})
    ref = family.merge(genus, on=["source_mode", "stratum"], validate="one_to_one")
    expected_profiles = len(modes) * len(strata)
    if len(ref) != expected_profiles:
        raise P1CNullError("true reference missing primary profiles")
    if (ref["family_norm"] <= 0).any():
        raise P1CNullError("non-positive family vector norm")
    ref["conditional_attenuation"] = 1.0 - ref["genus_norm"] / ref["family_norm"]
    observed = float(ref["conditional_attenuation"].median())
    by_stratum = {
        str(stratum): float(
            ref.loc[ref["stratum"].astype(str).eq(stratum), "conditional_attenuation"].median()
        )
        for stratum in strata
    }
    return ref, observed, by_stratum


def _run_one_permutation(
    data: dict[str, Any],
    family_reference: pd.DataFrame,
    *,
    permutation_id: int,
) -> dict[str, Any]:
    p1 = data["p1"]
    null = p1["p1c_matched_complexity_genus_null"]
    primary = p1["primary_context"]
    axes = [str(x) for x in primary["response_axes"]]
    modes = [str(x) for x in primary["source_modes"]]
    strata = [str(x) for x in primary["strata"]]
    seed = int(null["random_seed"]) + int(permutation_id)
    pseudo_taxonomy, partition_audit = make_matched_pseudo_taxonomy(
        data["taxonomy"], seed=seed
    )
    if not partition_audit["group_size_multiset_preserved_all_families"]:
        return {"permutation_id": permutation_id, "seed": seed, "valid": False}

    positions, availability, _, _ = build_source_group_contract(
        data["direct_scores"],
        pseudo_taxonomy,
        data["gift_flora"],
        axes=axes,
        minimum_source_scored_species=2,
        matched_taxa=data["matched_taxa"],
    )
    minimum = data["spec"]["minimum_support"]
    decomposition = build_taxonomic_decomposition(
        data["direct_scores"],
        pseudo_taxonomy,
        data["status_flora"],
        positions,
        availability,
        data["assignments"],
        data["covariates"],
        axes=axes,
        source_modes=modes,
        strata=strata,
        minimum_species=int(minimum["observed_species_per_island_response"]),
        minimum_families=int(minimum["represented_families"]),
        minimum_genera=int(minimum["represented_genera"]),
    )
    slopes = fit_palearctic_profiles(
        decomposition,
        data["covariates"],
        data["realm_assignment"],
        data["pattern"],
        stage="after_genus_residual",
        source_modes=modes,
        strata=strata,
        axes=axes,
    )
    pseudo = _vector_norms(slopes, axes).rename(columns={"vector_norm": "pseudo_genus_norm"})
    joined = family_reference[["source_mode", "stratum", "family_norm"]].merge(
        pseudo,
        on=["source_mode", "stratum"],
        how="inner",
        validate="one_to_one",
    )
    expected_profiles = len(modes) * len(strata)
    if len(joined) != expected_profiles:
        return {
            "permutation_id": permutation_id,
            "seed": seed,
            "valid": False,
            "n_profiles_fit": int(len(joined)),
        }
    joined["conditional_attenuation"] = (
        1.0 - joined["pseudo_genus_norm"] / joined["family_norm"]
    )
    row: dict[str, Any] = {
        "permutation_id": int(permutation_id),
        "seed": seed,
        "valid": True,
        "n_profiles_fit": int(len(joined)),
        "primary_median_conditional_attenuation": float(
            joined["conditional_attenuation"].median()
        ),
        "minimum_conditional_attenuation": float(joined["conditional_attenuation"].min()),
    }
    for stratum in strata:
        row[f"median_{stratum}"] = float(
            joined.loc[
                joined["stratum"].astype(str).eq(stratum), "conditional_attenuation"
            ].median()
        )
    for _, profile in joined.sort_values(["stratum", "source_mode"]).iterrows():
        key = f"attenuation__{profile['stratum']}__{profile['source_mode']}"
        row[key] = float(profile["conditional_attenuation"])
    return row


def run_shard(
    *,
    artifact_root: Path,
    taxonomy_path: Path,
    p1_config_path: Path,
    pattern_config_path: Path,
    source_config_path: Path,
    explanation_config_path: Path,
    shard_id: int,
    output_dir: Path,
) -> dict[str, Any]:
    data = _load_inputs(
        artifact_root,
        taxonomy_path,
        p1_config_path,
        pattern_config_path,
        source_config_path,
        explanation_config_path,
    )
    null = data["p1"]["p1c_matched_complexity_genus_null"]
    shards = int(null["execution_shards"])
    per_shard = int(null["permutations_per_shard"])
    total = int(null["permutations"])
    if shards * per_shard != total:
        raise P1CNullError("execution shards do not cover frozen permutation total")
    if shard_id < 0 or shard_id >= shards:
        raise P1CNullError("shard_id outside frozen range")

    family_reference, observed, observed_by_stratum = _true_reference(data)
    data["matched_taxa"] = match_gift_species(
        data["gift_flora"], data["taxonomy"][["accepted_species"]]
    )
    start = shard_id * per_shard
    stop = start + per_shard
    rows = [
        _run_one_permutation(data, family_reference, permutation_id=permutation_id)
        for permutation_id in range(start, stop)
    ]
    frame = pd.DataFrame(rows).sort_values("permutation_id")
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"p1c_null_shard_{shard_id:02d}.csv.gz"
    frame.to_csv(path, index=False)
    manifest = {
        "contract": "chapter1_p1c_matched_genus_null_shard_v1",
        "shard_id": int(shard_id),
        "permutation_start_inclusive": int(start),
        "permutation_stop_exclusive": int(stop),
        "n_rows": int(len(frame)),
        "n_valid": int(frame.get("valid", pd.Series(dtype=bool)).astype(bool).sum()),
        "observed_primary_statistic": observed,
        "observed_median_by_stratum": observed_by_stratum,
        "taxonomy_sha256_newline_canonicalized": data["taxonomy_hash"],
        "source_artifact_digest": str(
            data["p1"]["pinned_primary_artifact"]["artifact_digest"]
        ),
    }
    (output_dir / f"p1c_null_shard_{shard_id:02d}_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def aggregate_shards(
    *,
    shard_dir: Path,
    p1_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    p1 = load_yaml(p1_config_path)
    null = p1["p1c_matched_complexity_genus_null"]
    expected_total = int(null["permutations"])
    expected_shards = int(null["execution_shards"])
    files = sorted(shard_dir.glob("p1c_null_shard_*.csv.gz"))
    if len(files) != expected_shards:
        raise P1CNullError(f"expected {expected_shards} shard CSVs, found {len(files)}")
    frame = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    if len(frame) != expected_total:
        raise P1CNullError("aggregated permutation count differs from frozen total")
    ids = pd.to_numeric(frame["permutation_id"], errors="raise").astype(int)
    if set(ids) != set(range(expected_total)) or ids.duplicated().any():
        raise P1CNullError("permutation IDs do not equal frozen 0..1999 schedule")

    manifests = sorted(shard_dir.glob("p1c_null_shard_*_manifest.json"))
    if len(manifests) != expected_shards:
        raise P1CNullError("missing shard manifests")
    loaded_manifests = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    observed_values = {float(item["observed_primary_statistic"]) for item in loaded_manifests}
    if len(observed_values) != 1:
        raise P1CNullError("observed statistic differs among shards")
    observed = observed_values.pop()

    valid = frame.loc[frame["valid"].astype(bool)].copy()
    n_valid = int(len(valid))
    minimum_valid = int(null["minimum_valid_permutations"])
    evaluable = n_valid >= minimum_valid
    randomization_p = None
    exceedances = None
    if evaluable:
        values = pd.to_numeric(
            valid["primary_median_conditional_attenuation"], errors="raise"
        )
        exceedances = int(values.ge(observed).sum())
        randomization_p = float((1 + exceedances) / (1 + n_valid))
    threshold = float(null["primary_pass_rule"]["max_randomization_p_value"])
    if not evaluable:
        classification = "NOT_EVALUABLE_insufficient_valid_permutations"
    elif randomization_p is not None and randomization_p <= threshold:
        classification = "true_genus_exceeds_matched_complexity_null"
    else:
        classification = "true_genus_not_distinguishable_from_matched_complexity_null"

    output_dir.mkdir(parents=True, exist_ok=True)
    frame.sort_values("permutation_id").to_csv(
        output_dir / "p1c_matched_genus_null_permutations.csv.gz", index=False
    )
    quantiles: dict[str, float] = {}
    if n_valid:
        series = pd.to_numeric(
            valid["primary_median_conditional_attenuation"], errors="coerce"
        ).dropna()
        for q in (0.025, 0.25, 0.5, 0.75, 0.975):
            quantiles[str(q)] = float(series.quantile(q))
    result = {
        "contract": "chapter1_p1c_matched_genus_null_result_v1",
        "status": classification,
        "permutations_requested": expected_total,
        "valid_permutations": n_valid,
        "minimum_valid_permutations": minimum_valid,
        "observed_primary_median_conditional_attenuation": observed,
        "null_primary_statistic_quantiles": quantiles,
        "null_exceedances_greater_equal_observed": exceedances,
        "one_sided_randomization_p_value": randomization_p,
        "pass_threshold": threshold,
        "strong_assembly_depth_defense_pass": bool(
            evaluable and randomization_p is not None and randomization_p <= threshold
        ),
        "seed_rule": str(null["permutation_seed_rule"]),
        "master_seed": int(null["random_seed"]),
        "primary_scope": "direct_only_Palearctic_two_axis",
        "claim_boundary": (
            "A positive result shows taxonomic localization beyond matched arbitrary "
            "within-family partitions. It remains non-causal and does not exclude "
            "within-lineage evolution or identify a pollinator mechanism."
        ),
    }
    (output_dir / "chapter1_p1c_matched_genus_null_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


@app.command("run-shard")
def run_shard_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    taxonomy_path: Path = typer.Option(..., exists=True, dir_okay=False),
    p1_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    pattern_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    source_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    explanation_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    shard_id: int = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    result = run_shard(
        artifact_root=artifact_root,
        taxonomy_path=taxonomy_path,
        p1_config_path=p1_config_path,
        pattern_config_path=pattern_config_path,
        source_config_path=source_config_path,
        explanation_config_path=explanation_config_path,
        shard_id=shard_id,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(result, indent=2))


@app.command("aggregate")
def aggregate_command(
    shard_dir: Path = typer.Option(..., exists=True, file_okay=False),
    p1_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            aggregate_shards(
                shard_dir=shard_dir,
                p1_config_path=p1_config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
