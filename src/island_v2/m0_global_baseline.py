"""Run the global M0 baseline without Bombus applicability gating.

M0 is the global descriptive/baseline model in the preregistered M0-M4 ladder.
It uses all islands with trait composition and available baseline covariates. Bombus
occurrence, applicability, and channel state are deliberately excluded from M0.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

CLIMATE_COLUMNS = [
    "bio1_median",
    "bio4_median",
    "bio5_median",
    "bio6_median",
    "bio12_median",
    "bio15_median",
    "bio18_median",
    "bio19_median",
    "elevation_median",
]


def _standardize(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.astype(float).copy()
    means = out.mean(axis=0)
    sds = out.std(axis=0, ddof=0).replace(0, np.nan)
    return (out - means) / sds


def build_climate_pcs(environment: pd.DataFrame) -> pd.DataFrame:
    missing = {"island_id", *CLIMATE_COLUMNS}.difference(environment.columns)
    if missing:
        raise typer.BadParameter(f"island environment missing columns: {sorted(missing)}")
    work = environment[["island_id", *CLIMATE_COLUMNS]].copy()
    for column in CLIMATE_COLUMNS:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    complete = work.dropna(subset=CLIMATE_COLUMNS).copy()
    if len(complete) < 3:
        raise typer.BadParameter("too few complete island climate rows for PCA")
    z = _standardize(complete[CLIMATE_COLUMNS])
    u, s, _ = np.linalg.svd(z.to_numpy(dtype=float), full_matrices=False)
    scores = u[:, :2] * s[:2]
    result = complete[["island_id"]].copy()
    result["climate_pc1"] = scores[:, 0]
    result["climate_pc2"] = scores[:, 1]
    return result


def _pivot_trait(composition: pd.DataFrame, trait_name: str) -> pd.DataFrame:
    work = composition.loc[composition["trait_name"].eq(trait_name)].copy()
    if work.empty:
        return pd.DataFrame(columns=["island_id"])
    work["species_share_within_trait"] = pd.to_numeric(
        work["species_share_within_trait"], errors="coerce"
    )
    return work.pivot_table(
        index="island_id",
        columns="filled_value",
        values="species_share_within_trait",
        aggfunc="sum",
        fill_value=0.0,
    ).reset_index()


def _fit_standardized_ols(data: pd.DataFrame, outcome: str, predictors: list[str]) -> dict[str, object]:
    work = data[[outcome, *predictors]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(work) <= len(predictors) + 1:
        return {"outcome": outcome, "status": "insufficient_complete_cases", "n": int(len(work))}
    y = work[outcome]
    y = (y - y.mean()) / y.std(ddof=0)
    xs = []
    for predictor in predictors:
        x = work[predictor]
        sd = x.std(ddof=0)
        xs.append(((x - x.mean()) / sd) if sd > 0 else pd.Series(0.0, index=x.index))
    X = np.column_stack([np.ones(len(work)), *[x.to_numpy(dtype=float) for x in xs]])
    beta, _, _, _ = np.linalg.lstsq(X, y.to_numpy(dtype=float), rcond=None)
    fitted = X @ beta
    resid = y.to_numpy(dtype=float) - fitted
    sse = float(np.sum(resid**2))
    sst = float(np.sum((y.to_numpy(dtype=float) - y.mean()) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    return {
        "outcome": outcome,
        "status": "fitted",
        "n": int(len(work)),
        "r2": float(r2),
        **{f"beta_{name}": float(value) for name, value in zip(predictors, beta[1:], strict=True)},
    }


@app.command("run")
def run(
    trait_composition_csv: Path = typer.Option(..., exists=True),
    island_environment_csv: Path = typer.Option(..., exists=True),
    island_metrics_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    analysis_tier: str = typer.Option("primary"),
) -> None:
    composition = pd.read_csv(trait_composition_csv)
    environment = pd.read_csv(island_environment_csv)
    metrics = pd.read_csv(island_metrics_csv)

    climate = build_climate_pcs(environment)
    base = metrics[["island_id", "n_trait_species"]].drop_duplicates("island_id")
    base = base.merge(climate, on="island_id", how="inner", validate="one_to_one")

    traits = {}
    for trait_name in (
        "pollination_functional_guild",
        "pollen_vector_mode",
        "autonomous_selfing_capacity",
        "self_incompatibility",
    ):
        traits[trait_name] = _pivot_trait(composition, trait_name)

    predictors = ["climate_pc1", "climate_pc2", "n_trait_species"]
    rows: list[dict[str, object]] = []
    outcome_manifest: list[str] = []
    for trait_name, table in traits.items():
        if table.empty:
            continue
        joined = base.merge(table, on="island_id", how="inner", validate="one_to_one")
        for outcome in [c for c in table.columns if c != "island_id"]:
            result = _fit_standardized_ols(joined, outcome, predictors)
            result["trait_name"] = trait_name
            result["analysis_tier"] = analysis_tier
            rows.append(result)
            outcome_manifest.append(f"{trait_name}:{outcome}")

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_dir / "m0_global_baseline_results.csv", index=False)
    climate.to_csv(output_dir / "m0_global_climate_pcs.csv", index=False)
    manifest = {
        "contract": "m0_global_baseline_v1",
        "analysis_tier": analysis_tier,
        "scope": "global_all_islands_with_trait_composition_and_complete_climate",
        "bombus_applicability_used": False,
        "bombus_channel_used": False,
        "n_islands_trait_metrics": int(base["island_id"].nunique()),
        "predictors": predictors,
        "missing_preregistered_baseline_blocks": [
            "geography/source_pool",
            "lineage_composition",
            "island_area",
            "isolation",
            "establishment",
            "observation_process",
        ],
        "outcomes": sorted(set(outcome_manifest)),
        "warning": "This is an executable global M0 baseline subset, not the final fully adjusted M0 until the missing preregistered baseline blocks are connected.",
    }
    (output_dir / "m0_global_baseline_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
