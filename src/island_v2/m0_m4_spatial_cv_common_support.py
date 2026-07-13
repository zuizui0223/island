"""Common-support repeated spatial CV across trait evidence tiers.

The ordinary tier-wise CV can differ because primary, broad, and sensitivity_all
contain different analyzable islands. This module intersects eligible island IDs
for each nested model comparison, assigns the same geographic blocks to the same
folds, and then evaluates every tier on exactly the same held-out islands.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.m0_m4_spatial_cv import (
    _comparison_specs,
    _complete_rows,
    assign_balanced_block_folds,
    fit_grouped_binomial,
    heldout_metrics,
    predict_probability,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

TIERS = ("primary", "broad", "sensitivity_all")


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("full-analysis config must be a mapping")
    return payload


def run_common_support_cv(
    tier_data: dict[str, pd.DataFrame],
    baseline: list[str],
    n_folds: int,
    repeats: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if set(tier_data) != set(TIERS):
        raise typer.BadParameter(f"tier_data must contain exactly {TIERS}")

    fold_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    for spec in _comparison_specs(baseline):
        union_predictors = list(dict.fromkeys([*spec["reduced"], *spec["full"]]))
        complete_by_tier = {
            tier: _complete_rows(
                table,
                spec["successes"],
                spec["trials"],
                union_predictors,
            )
            for tier, table in tier_data.items()
        }
        common_ids = set.intersection(
            *(set(table["island_id"].astype(str)) for table in complete_by_tier.values())
        )
        if len(common_ids) < 50:
            raise typer.BadParameter(
                f"{spec['comparison']}: fewer than 50 common-support islands ({len(common_ids)})"
            )
        common_by_tier = {
            tier: table.loc[table["island_id"].astype(str).isin(common_ids)].copy()
            for tier, table in complete_by_tier.items()
        }
        reference = common_by_tier["primary"][["island_id", "spatial_block"]].drop_duplicates()
        if reference["island_id"].duplicated().any():
            raise RuntimeError("primary common-support table has duplicate island IDs")
        for tier, table in common_by_tier.items():
            block_check = table[["island_id", "spatial_block"]].drop_duplicates().merge(
                reference,
                on="island_id",
                how="inner",
                suffixes=("_tier", "_reference"),
                validate="one_to_one",
            )
            if not block_check["spatial_block_tier"].astype(str).eq(
                block_check["spatial_block_reference"].astype(str)
            ).all():
                raise RuntimeError(f"{tier}: spatial block mismatch on common support")

        support_rows.append(
            {
                "comparison": spec["comparison"],
                "n_common_islands": int(len(common_ids)),
                "n_common_spatial_blocks": int(reference["spatial_block"].nunique()),
                **{
                    f"n_{tier}_eligible_before_intersection": int(len(table))
                    for tier, table in complete_by_tier.items()
                },
            }
        )

        for repeat in range(repeats):
            repeat_seed = seed + repeat * 1009
            mapping = assign_balanced_block_folds(reference, n_folds, repeat_seed)
            for tier, table in common_by_tier.items():
                fold_id = table["spatial_block"].astype(str).map(mapping)
                if fold_id.isna().any():
                    raise RuntimeError(f"{tier}: missing common fold assignment")
                for fold in range(n_folds):
                    test = table.loc[fold_id.eq(fold)].copy()
                    train = table.loc[~fold_id.eq(fold)].copy()
                    reduced_model = fit_grouped_binomial(
                        train,
                        spec["successes"],
                        spec["trials"],
                        spec["reduced"],
                    )
                    full_model = fit_grouped_binomial(
                        train,
                        spec["successes"],
                        spec["trials"],
                        spec["full"],
                    )
                    reduced = heldout_metrics(
                        test,
                        spec["successes"],
                        spec["trials"],
                        predict_probability(reduced_model, test),
                    )
                    full = heldout_metrics(
                        test,
                        spec["successes"],
                        spec["trials"],
                        predict_probability(full_model, test),
                    )
                    fold_rows.append(
                        {
                            "comparison": spec["comparison"],
                            "tier": tier,
                            "repeat": repeat + 1,
                            "fold": fold + 1,
                            "seed": repeat_seed,
                            "n_train_islands": int(len(train)),
                            "n_test_islands": int(len(test)),
                            "n_test_trials": full["n_trials"],
                            "delta_log_likelihood_full_minus_reduced": (
                                full["log_likelihood"] - reduced["log_likelihood"]
                            ),
                            "delta_brier_full_minus_reduced": (
                                full["weighted_brier"] - reduced["weighted_brier"]
                            ),
                        }
                    )

    folds = pd.DataFrame(fold_rows)
    repeat_rows: list[dict[str, Any]] = []
    for (comparison, tier, repeat), group in folds.groupby(
        ["comparison", "tier", "repeat"], sort=True
    ):
        repeat_rows.append(
            {
                "comparison": comparison,
                "tier": tier,
                "repeat": int(repeat),
                "n_folds": int(len(group)),
                "total_test_trials": float(group["n_test_trials"].sum()),
                "delta_log_likelihood_per_1000_trials": float(
                    1000.0
                    * group["delta_log_likelihood_full_minus_reduced"].sum()
                    / group["n_test_trials"].sum()
                ),
                "mean_delta_brier_full_minus_reduced": float(
                    np.average(
                        group["delta_brier_full_minus_reduced"],
                        weights=group["n_test_trials"],
                    )
                ),
            }
        )
    repeats_table = pd.DataFrame(repeat_rows)

    summary_rows: list[dict[str, Any]] = []
    for (comparison, tier), group in repeats_table.groupby(["comparison", "tier"], sort=True):
        ll = group["delta_log_likelihood_per_1000_trials"].to_numpy(dtype=float)
        brier = group["mean_delta_brier_full_minus_reduced"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "comparison": comparison,
                "tier": tier,
                "n_repeats": int(len(group)),
                "median_delta_log_likelihood_per_1000_trials": float(np.median(ll)),
                "min_delta_log_likelihood_per_1000_trials": float(np.min(ll)),
                "max_delta_log_likelihood_per_1000_trials": float(np.max(ll)),
                "proportion_repeats_full_model_better_loglik": float(np.mean(ll > 0)),
                "median_delta_brier_full_minus_reduced": float(np.median(brier)),
                "proportion_repeats_full_model_better_brier": float(np.mean(brier < 0)),
                "predictive_gain_consistent": bool(np.all(ll > 0) and np.all(brier < 0)),
            }
        )
    return folds, repeats_table, pd.DataFrame(summary_rows), pd.DataFrame(support_rows)


@app.command("run")
def run(
    primary_csv: Path = typer.Option(..., exists=True),
    broad_csv: Path = typer.Option(..., exists=True),
    sensitivity_all_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_full_analysis.yml"), exists=True),
    n_folds: int = typer.Option(10, min=5),
    repeats: int = typer.Option(5, min=1),
    seed: int = typer.Option(20260713),
) -> None:
    config = _load_config(config_path)
    tier_data = {
        "primary": pd.read_csv(primary_csv),
        "broad": pd.read_csv(broad_csv),
        "sensitivity_all": pd.read_csv(sensitivity_all_csv),
    }
    baseline = [str(value) for value in config["baseline_covariates"]]
    folds, repeat_table, summary, support = run_common_support_cv(
        tier_data,
        baseline,
        n_folds=n_folds,
        repeats=repeats,
        seed=seed,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    folds.to_csv(output_dir / "common_support_spatial_cv_fold_results.csv", index=False)
    repeat_table.to_csv(output_dir / "common_support_spatial_cv_repeat_results.csv", index=False)
    summary.to_csv(output_dir / "common_support_spatial_cv_summary.csv", index=False)
    support.to_csv(output_dir / "common_support_counts.csv", index=False)
    manifest = {
        "contract": "m0_m4_common_support_spatial_cv_v1",
        "n_folds": n_folds,
        "repeats": repeats,
        "seed": seed,
        "same_islands_across_tiers": True,
        "same_spatial_fold_assignment_across_tiers": True,
        "support": support.to_dict("records"),
        "results": summary.to_dict("records"),
    }
    (output_dir / "common_support_spatial_cv_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
