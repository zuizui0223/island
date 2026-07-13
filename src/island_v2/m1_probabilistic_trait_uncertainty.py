"""Propagate genus/family trait uncertainty through the M1 Bombus association.

The hard broad tier treats inferred categorical values as known. This sensitivity
instead converts each species' stored value_distribution into a bumblebee-guild
probability, fixes species-direct cells at 0/1, excludes global fallback, and draws
one trait state per species per imputation. The sampled state is then shared across
all island occurrences of that species.

The resulting coefficient distribution is pooled with Rubin-style within- plus
between-imputation variance. Because the empirical taxonomic distributions are
not a formal missing-data posterior, the output remains a sensitivity analysis.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.m0_m4_full_analysis import fit_grouped_binomial_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)


@dataclass
class ProbabilityMapping:
    species_table: pd.DataFrame
    island_table: pd.DataFrame
    island_codes: np.ndarray
    species_codes: np.ndarray
    trials: np.ndarray


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("probabilistic trait-uncertainty config must be a mapping")
    return payload


def _distribution_probability(
    distribution_json: str,
    positive_value: str,
) -> float:
    try:
        payload = json.loads(str(distribution_json)) if str(distribution_json).strip() else {}
    except json.JSONDecodeError:
        payload = {}
    if not isinstance(payload, dict):
        return np.nan
    weights: dict[str, float] = {}
    for key, raw in payload.items():
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value) and value > 0:
            weights[str(key)] = value
    total = float(sum(weights.values()))
    return float(weights.get(positive_value, 0.0) / total) if total > 0 else np.nan


def build_probability_mapping(
    species_master: pd.DataFrame,
    reference_island_data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[ProbabilityMapping, dict[str, Any]]:
    outcome = config["outcome"]
    included = {str(value) for value in outcome["included_fill_tiers"]}
    required = {
        "island_id",
        "accepted_species",
        "trait_name",
        "filled_value",
        "fill_tier",
        "value_distribution",
    }
    missing = required.difference(species_master.columns)
    if missing:
        raise typer.BadParameter(f"species master missing columns: {sorted(missing)}")

    rows = species_master.loc[
        species_master["trait_name"].astype(str).eq(str(outcome["trait_name"]))
        & species_master["fill_tier"].astype(str).isin(included)
    ][list(required)].copy()
    rows = rows.dropna(subset=["island_id", "accepted_species"])
    rows = rows.drop_duplicates(["island_id", "accepted_species"])
    if rows.empty:
        raise typer.BadParameter("no direct/genus/family outcome rows remain")

    # The cascade is species-level. A species must have exactly one tier/value/distribution
    # across all islands before a species-level imputation is legitimate.
    consistency = rows.groupby("accepted_species").agg(
        n_fill_tiers=("fill_tier", "nunique"),
        n_values=("filled_value", "nunique"),
        n_distributions=("value_distribution", "nunique"),
    )
    bad = consistency.loc[
        (consistency["n_fill_tiers"] > 1)
        | (consistency["n_values"] > 1)
        | (consistency["n_distributions"] > 1)
    ]
    if len(bad):
        raise typer.BadParameter(
            "species trait inference is inconsistent across islands: "
            f"{bad.head(20).to_dict('index')}"
        )

    species = rows[
        ["accepted_species", "fill_tier", "filled_value", "value_distribution"]
    ].drop_duplicates("accepted_species")
    positive_value = str(outcome["positive_value"])
    probabilities: list[float] = []
    for row in species.itertuples(index=False):
        if str(row.fill_tier) == "species_direct":
            probability = 1.0 if str(row.filled_value) == positive_value else 0.0
        else:
            probability = _distribution_probability(
                str(row.value_distribution), positive_value
            )
        probabilities.append(probability)
    species["positive_probability"] = probabilities
    species = species.dropna(subset=["positive_probability"]).copy()
    if not species["positive_probability"].between(0, 1, inclusive="both").all():
        raise typer.BadParameter("trait probabilities fall outside [0, 1]")
    species = species.sort_values("accepted_species").reset_index(drop=True)
    species["species_code"] = np.arange(len(species), dtype=int)

    rows = rows.loc[rows["accepted_species"].isin(species["accepted_species"])].copy()
    rows = rows.merge(
        species[["accepted_species", "species_code"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )

    baseline_specs = {
        str(key): [str(value) for value in values]
        for key, values in config["baseline_specs"].items()
    }
    baseline_union = list(
        dict.fromkeys(value for values in baseline_specs.values() for value in values)
    )
    bombus = str(config["bombus_predictor"])
    spatial_cluster = str(config["spatial_cluster"])
    reference_columns = ["island_id", spatial_cluster, bombus, *baseline_union]
    missing_reference = set(reference_columns).difference(reference_island_data.columns)
    if missing_reference:
        raise typer.BadParameter(
            f"reference island data missing columns: {sorted(missing_reference)}"
        )
    reference = reference_island_data[reference_columns].drop_duplicates("island_id")
    if reference["island_id"].duplicated().any():
        raise typer.BadParameter("reference island data must contain one row per island")

    island_ids = sorted(set(rows["island_id"].astype(str)) & set(reference["island_id"].astype(str)))
    island_table = reference.loc[reference["island_id"].astype(str).isin(island_ids)].copy()
    island_table = island_table.sort_values("island_id").reset_index(drop=True)
    island_table["island_code"] = np.arange(len(island_table), dtype=int)
    island_lookup = island_table.set_index("island_id")["island_code"]
    rows["island_code"] = rows["island_id"].map(island_lookup)
    rows = rows.dropna(subset=["island_code"]).copy()
    rows["island_code"] = rows["island_code"].astype(int)

    island_codes = rows["island_code"].to_numpy(dtype=int)
    species_codes = rows["species_code"].to_numpy(dtype=int)
    trials = np.bincount(island_codes, minlength=len(island_table)).astype(float)
    island_table["trait_trials"] = trials

    metadata = {
        "n_species": int(len(species)),
        "n_island_species_rows": int(len(rows)),
        "n_islands": int(len(island_table)),
        "n_spatial_blocks": int(island_table[spatial_cluster].nunique()),
        "species_by_fill_tier": {
            str(key): int(value)
            for key, value in species["fill_tier"].value_counts().sort_index().items()
        },
        "island_species_rows_by_fill_tier": {
            str(key): int(value)
            for key, value in rows["fill_tier"].value_counts().sort_index().items()
        },
        "n_uncertain_species": int(
            species["positive_probability"].between(0, 1, inclusive="neither").sum()
        ),
        "mean_inferred_positive_probability": float(
            species.loc[
                ~species["fill_tier"].eq("species_direct"),
                "positive_probability",
            ].mean()
        ),
    }
    return (
        ProbabilityMapping(
            species_table=species,
            island_table=island_table,
            island_codes=island_codes,
            species_codes=species_codes,
            trials=trials,
        ),
        metadata,
    )


def aggregate_species_states(
    mapping: ProbabilityMapping,
    species_state: np.ndarray,
) -> np.ndarray:
    if len(species_state) != len(mapping.species_table):
        raise typer.BadParameter("species state vector length does not match species table")
    row_state = species_state[mapping.species_codes]
    return np.bincount(
        mapping.island_codes,
        weights=row_state,
        minlength=len(mapping.island_table),
    ).astype(float)


def _fit_snapshot(
    island_table: pd.DataFrame,
    successes: np.ndarray,
    trials: np.ndarray,
    baseline: list[str],
    bombus_predictor: str,
    spatial_cluster: str,
    z_threshold: float,
) -> dict[str, Any]:
    data = island_table.copy()
    data["trait_successes"] = successes
    data["trait_trials"] = trials
    data = data.loc[data["trait_trials"] > 0].copy()

    m0_coefficients, m0_fit = fit_grouped_binomial_clustered(
        data,
        "trait_successes",
        "trait_trials",
        baseline,
        spatial_cluster,
        z_threshold,
    )
    m1_coefficients, m1_fit = fit_grouped_binomial_clustered(
        data,
        "trait_successes",
        "trait_trials",
        [*baseline, bombus_predictor],
        spatial_cluster,
        z_threshold,
    )
    if m0_fit.get("status") not in {"fitted", "max_iter_reached"}:
        raise RuntimeError(f"M0 fit failed: {m0_fit}")
    if m1_fit.get("status") not in {"fitted", "max_iter_reached"}:
        raise RuntimeError(f"M1 fit failed: {m1_fit}")
    bombus_row = m1_coefficients.loc[
        m1_coefficients["predictor"].eq(bombus_predictor)
    ]
    if len(bombus_row) != 1:
        raise RuntimeError("M1 Bombus coefficient is missing or duplicated")
    row = bombus_row.iloc[0]
    return {
        "estimate": float(row["estimate_log_odds"]),
        "within_se": float(row["cluster_robust_se"]),
        "z_value": float(row["z_value"]),
        "n_islands": int(m1_fit["n_islands"]),
        "n_spatial_blocks": int(m1_fit["n_clusters"]),
        "m0_aic": float(m0_fit["aic"]),
        "m1_aic": float(m1_fit["aic"]),
        "delta_aic_m1_minus_m0": float(m1_fit["aic"] - m0_fit["aic"]),
        "m0_pseudo_r2": float(m0_fit["mcfadden_pseudo_r2"]),
        "m1_pseudo_r2": float(m1_fit["mcfadden_pseudo_r2"]),
    }


def rubin_style_pool(
    draws: pd.DataFrame,
    baseline_spec: str,
    z_threshold: float,
) -> dict[str, Any]:
    work = draws.loc[draws["baseline_spec"].eq(baseline_spec)].copy()
    work = work.dropna(subset=["estimate", "within_se"])
    m = int(len(work))
    if m < 2:
        raise typer.BadParameter(f"too few valid imputations for {baseline_spec}")
    estimates = work["estimate"].to_numpy(dtype=float)
    within_variances = np.square(work["within_se"].to_numpy(dtype=float))
    q_bar = float(np.mean(estimates))
    u_bar = float(np.mean(within_variances))
    between = float(np.var(estimates, ddof=1))
    total_variance = float(u_bar + (1.0 + 1.0 / m) * between)
    total_se = float(math.sqrt(total_variance))
    z_value = float(q_bar / total_se) if total_se > 0 else np.nan
    empirical_lower, empirical_upper = np.quantile(estimates, [0.025, 0.975])
    return {
        "baseline_spec": baseline_spec,
        "n_imputations": m,
        "pooled_estimate_log_odds": q_bar,
        "pooled_odds_ratio_per_predictor_sd": float(math.exp(q_bar)),
        "mean_within_imputation_variance": u_bar,
        "between_imputation_variance": between,
        "total_variance": total_variance,
        "pooled_total_se": total_se,
        "pooled_z_value": z_value,
        "pooled_two_sided_normal_p": float(
            math.erfc(abs(z_value) / math.sqrt(2.0))
        )
        if math.isfinite(z_value)
        else np.nan,
        "pooled_nominally_supported": bool(abs(z_value) >= z_threshold)
        if math.isfinite(z_value)
        else False,
        "empirical_estimate_2_5": float(empirical_lower),
        "empirical_estimate_97_5": float(empirical_upper),
        "proportion_positive_estimates": float(np.mean(estimates > 0)),
        "proportion_imputations_nominally_supported": float(
            np.mean(np.abs(work["z_value"].to_numpy(dtype=float)) >= z_threshold)
        ),
        "mean_delta_aic_m1_minus_m0": float(
            work["delta_aic_m1_minus_m0"].mean()
        ),
        "proportion_imputations_m1_lower_aic": float(
            np.mean(work["delta_aic_m1_minus_m0"].to_numpy(dtype=float) < 0)
        ),
        "monte_carlo_se_of_pooled_estimate": float(math.sqrt(between / m)),
    }


def run_imputations(
    mapping: ProbabilityMapping,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    imputation = config["imputation"]
    n_imputations = int(imputation["n_imputations"])
    seed = int(imputation["random_seed"])
    z_threshold = float(config["nominal_z_threshold"])
    bombus = str(config["bombus_predictor"])
    spatial_cluster = str(config["spatial_cluster"])
    baselines = {
        str(key): [str(value) for value in values]
        for key, values in config["baseline_specs"].items()
    }

    probabilities = mapping.species_table["positive_probability"].to_numpy(dtype=float)
    hard_state = mapping.species_table["filled_value"].astype(str).eq(
        str(config["outcome"]["positive_value"])
    ).astype(float).to_numpy()
    expected_successes = aggregate_species_states(mapping, probabilities)
    hard_successes = aggregate_species_states(mapping, hard_state)

    comparator_rows: list[dict[str, Any]] = []
    for comparator, successes in (
        ("expected_probability_counts", expected_successes),
        ("hard_modal_broad_without_global_fallback", hard_successes),
    ):
        for baseline_name, baseline in baselines.items():
            result = _fit_snapshot(
                mapping.island_table,
                successes,
                mapping.trials,
                baseline,
                bombus,
                spatial_cluster,
                z_threshold,
            )
            comparator_rows.append(
                {
                    "comparator": comparator,
                    "baseline_spec": baseline_name,
                    **result,
                }
            )

    rng = np.random.default_rng(seed)
    draw_rows: list[dict[str, Any]] = []
    for imputation_index in range(1, n_imputations + 1):
        species_state = (
            rng.random(len(probabilities)) < probabilities
        ).astype(float)
        successes = aggregate_species_states(mapping, species_state)
        for baseline_name, baseline in baselines.items():
            result = _fit_snapshot(
                mapping.island_table,
                successes,
                mapping.trials,
                baseline,
                bombus,
                spatial_cluster,
                z_threshold,
            )
            draw_rows.append(
                {
                    "imputation": imputation_index,
                    "baseline_spec": baseline_name,
                    **result,
                }
            )
    draws = pd.DataFrame(draw_rows)
    pooled = pd.DataFrame(
        [
            rubin_style_pool(draws, baseline_name, z_threshold)
            for baseline_name in baselines
        ]
    )
    return draws, pooled, pd.DataFrame(comparator_rows)


@app.command("run")
def run(
    species_master_csv: Path = typer.Option(..., exists=True),
    reference_island_data_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/m1_probabilistic_trait_uncertainty.yml"), exists=True
    ),
) -> None:
    config = _load_config(config_path)
    master = pd.read_csv(
        species_master_csv,
        usecols=[
            "island_id",
            "accepted_species",
            "trait_name",
            "filled_value",
            "fill_tier",
            "value_distribution",
        ],
    )
    reference = pd.read_csv(reference_island_data_csv)
    mapping, metadata = build_probability_mapping(master, reference, config)
    draws, pooled, comparators = run_imputations(mapping, config)

    minimum = int(config["reporting"].get("minimum_imputations", 100))
    if int(draws["imputation"].nunique()) < minimum:
        raise RuntimeError("fewer valid imputations than the declared minimum")

    output_dir.mkdir(parents=True, exist_ok=True)
    mapping.species_table.to_csv(
        output_dir / "probabilistic_species_trait_table.csv.gz",
        index=False,
        compression="gzip",
    )
    mapping.island_table.to_csv(
        output_dir / "probabilistic_analysis_islands.csv", index=False
    )
    draws.to_csv(output_dir / "imputation_m1_draws.csv", index=False)
    pooled.to_csv(output_dir / "imputation_m1_pooled.csv", index=False)
    comparators.to_csv(output_dir / "imputation_m1_comparators.csv", index=False)

    manifest = {
        "contract": str(config["analysis_contract"]),
        "analysis_role": str(config["analysis_role"]),
        **metadata,
        "n_imputations": int(draws["imputation"].nunique()),
        "pooled_results": pooled.to_dict("records"),
        "comparators": comparators.to_dict("records"),
        "interpretation": (
            "Species-direct cells are fixed; genus/family cells are sampled once per species per imputation from stored taxonomic value distributions. "
            "Global fallback is excluded. This propagates trait inference uncertainty but does not upgrade inferred cells to direct evidence."
        ),
    }
    (output_dir / "probabilistic_trait_uncertainty_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
