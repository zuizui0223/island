"""Build an outcome-blind source-exposure curation queue for PR136 H2a.

The queue prioritizes islands whose source-channel applicability is unresolved but
that already have enough *support* to contribute to a genus-fixed residual-vector
test once source exposure is reviewed. It deliberately ignores trait values,
residual signs, coefficients, p-values, and pollination-syndrome expectations.

Priority is balanced across biogeographic context and distance bins so curation
does not simply chase the best-observed northern islands. The queue decides what
to review next; it never assigns ``applicable`` or
``structurally_not_applicable`` itself.
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


def _distance_bin(values: pd.Series, n_bins: int) -> pd.Series:
    if n_bins < 1:
        raise typer.BadParameter("distance_bins must be positive")
    numeric = pd.to_numeric(values, errors="coerce")
    ranks = numeric.rank(method="first", pct=True)
    bins = np.minimum(np.floor(np.clip(ranks - 1e-12, 0, 1) * n_bins).astype(int), n_bins - 1)
    return pd.Series([f"D{x + 1}" for x in bins], index=values.index)


def build_source_exposure_queue(
    genus_null: pd.DataFrame,
    applicability: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    priority_strata = [str(x) for x in config["priority_strata"]]
    minimum_outcomes = int(config["minimum_outcomes_per_stratum"])
    context = str(config["context_column"])
    geography = str(config["geography_column"])
    n_bins = int(config["distance_bins"])
    per_cell = int(config["wave_size_per_context_distance_bin"])
    if per_cell < 1:
        raise typer.BadParameter("wave_size_per_context_distance_bin must be positive")

    required_null = {"island_id", "stratum", "outcome", "observed_n_species"}
    required_app = {"island_id", "applicability"}
    required_cov = {"island_id", context, geography}
    for table, required, label in (
        (genus_null, required_null, "genus-null"),
        (applicability, required_app, "applicability"),
        (covariates, required_cov, "covariates"),
    ):
        missing = required - set(table.columns)
        if missing:
            raise typer.BadParameter(f"{label} table missing columns: {sorted(missing)}")

    # Only support columns enter the queue. In particular,
    # deviation_observed_minus_null is deliberately never selected below.
    support = genus_null.loc[
        genus_null["stratum"].astype(str).isin(priority_strata),
        ["island_id", "stratum", "outcome", "observed_n_species"],
    ].copy()
    support["observed_n_species"] = pd.to_numeric(
        support["observed_n_species"], errors="coerce"
    )
    support = support.dropna(subset=["observed_n_species"])
    support = support.loc[support["observed_n_species"].gt(0)]

    by_stratum = (
        support.groupby(["island_id", "stratum"], as_index=False)
        .agg(
            n_supported_outcomes=("outcome", "nunique"),
            minimum_direct_trials=("observed_n_species", "min"),
            median_direct_trials=("observed_n_species", "median"),
        )
    )
    outcomes_wide = by_stratum.pivot(
        index="island_id", columns="stratum", values="n_supported_outcomes"
    )
    trials_wide = by_stratum.pivot(
        index="island_id", columns="stratum", values="minimum_direct_trials"
    )

    rows: list[dict[str, Any]] = []
    for island_id in sorted(set(by_stratum["island_id"].astype(str))):
        outcome_counts = [
            float(outcomes_wide.loc[island_id, stratum])
            if stratum in outcomes_wide.columns and pd.notna(outcomes_wide.loc[island_id, stratum])
            else 0.0
            for stratum in priority_strata
        ]
        trial_counts = [
            float(trials_wide.loc[island_id, stratum])
            if stratum in trials_wide.columns and pd.notna(trials_wide.loc[island_id, stratum])
            else 0.0
            for stratum in priority_strata
        ]
        row: dict[str, Any] = {
            "island_id": island_id,
            "minimum_supported_outcomes_across_priority_strata": int(min(outcome_counts)),
            "minimum_direct_trials_across_priority_strata": float(min(trial_counts)),
        }
        for stratum, count in zip(priority_strata, outcome_counts, strict=True):
            row[f"n_supported_outcomes__{stratum}"] = int(count)
        rows.append(row)
    queue = pd.DataFrame(rows)

    app_table = applicability[["island_id", "applicability"]].drop_duplicates("island_id")
    cov_table = covariates[["island_id", context, geography]].drop_duplicates("island_id")
    queue = queue.merge(app_table, on="island_id", how="left", validate="one_to_one")
    queue = queue.merge(cov_table, on="island_id", how="left", validate="one_to_one")
    queue["applicability"] = queue["applicability"].fillna("unresolved").astype(str)
    queue[context] = queue[context].fillna("").astype(str)
    queue[geography] = pd.to_numeric(queue[geography], errors="coerce")

    queue = queue.loc[
        queue["applicability"].eq("unresolved")
        & queue["minimum_supported_outcomes_across_priority_strata"].ge(minimum_outcomes)
        & queue[context].ne("")
        & queue[geography].notna()
    ].copy()
    if queue.empty:
        return queue, {
            "contract": "chapter1_pr136_source_exposure_queue_v1",
            "n_eligible_unresolved_islands": 0,
            "n_wave1_islands": 0,
        }

    queue["distance_support_bin"] = queue.groupby(context, group_keys=False)[geography].apply(
        lambda x: _distance_bin(x, n_bins)
    )
    queue = queue.sort_values(
        [
            context,
            "distance_support_bin",
            "minimum_supported_outcomes_across_priority_strata",
            "minimum_direct_trials_across_priority_strata",
            "island_id",
        ],
        ascending=[True, True, False, False, True],
    ).reset_index(drop=True)
    queue["rank_within_context_distance_bin"] = (
        queue.groupby([context, "distance_support_bin"]).cumcount() + 1
    )
    queue["curation_wave"] = np.ceil(
        queue["rank_within_context_distance_bin"] / per_cell
    ).astype(int)
    queue["wave1_priority"] = queue["curation_wave"].eq(1)
    queue["priority_basis"] = (
        "unresolved_source_exposure; genus_fixed_support_only; balanced_by_context_and_distance"
    )
    queue["trait_values_used_for_priority"] = False
    queue["residual_effect_direction_used_for_priority"] = False

    manifest = {
        "contract": "chapter1_pr136_source_exposure_queue_v1",
        "priority_strata": priority_strata,
        "minimum_outcomes_per_stratum": minimum_outcomes,
        "distance_bins": n_bins,
        "wave_size_per_context_distance_bin": per_cell,
        "n_eligible_unresolved_islands": int(len(queue)),
        "n_wave1_islands": int(queue["wave1_priority"].sum()),
        "eligible_by_context": {
            str(k): int(v) for k, v in queue[context].value_counts().sort_index().items()
        },
        "wave1_by_context": {
            str(k): int(v)
            for k, v in queue.loc[queue["wave1_priority"], context]
            .value_counts()
            .sort_index()
            .items()
        },
        "outcome_blind": True,
        "forbidden_priority_inputs": [
            "trait_state_value",
            "genus_fixed_residual_value",
            "effect_direction",
            "p_value",
            "pollination_syndrome_match",
        ],
        "interpretation": (
            "Curation order only. The queue does not infer source exposure or channel state."
        ),
    }
    return queue, manifest


@app.command("run")
def run(
    genus_null_csv: Path = typer.Option(..., exists=True),
    applicability_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_pr136_source_exposure_queue.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    queue, manifest = build_source_exposure_queue(
        pd.read_csv(genus_null_csv),
        pd.read_csv(applicability_csv),
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    queue.to_csv(output_dir / "source_exposure_curation_queue.csv", index=False)
    (output_dir / "source_exposure_curation_queue_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
