"""Observation-process audit for Chapter 1 WHEN/WHERE inference.

All islands in the frozen geographic universe are retained. Islands without flora or
trait records are explicit observation-process zeros; they are never interpreted as
biological zeros for floral or reproductive traits.

The module separates two stages:

1. flora observability: does an island have any realised flora record?
2. direct-trait observability: conditional on or irrespective of flora records, does
   the island have enough direct species-level evidence to enter the Chapter 1 trait
   composition surface?

Selection models use island area, climate, a mainland-distance/source-pool-accessibility
gradient, and biogeographic context. Their coefficients describe evidence availability,
not causal ecological filtering.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _build_support_table(
    all_islands: pd.DataFrame,
    island_species: pd.DataFrame,
    composition: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required_islands = {
        "island_id",
        config["context_column"],
        config["cluster_column"],
        config["distance_gradient"],
        *config["baseline_covariates"],
    }
    missing = required_islands - set(all_islands.columns)
    if missing:
        raise typer.BadParameter(f"all-island table missing columns: {sorted(missing)}")
    if all_islands["island_id"].duplicated().any():
        raise typer.BadParameter("all-island table contains duplicate island_id values")
    if "island_id" not in island_species.columns or "species" not in island_species.columns:
        raise typer.BadParameter("island-species table must contain island_id and species")
    required_comp = {"island_id", "trait_name", "trials", "stratum"}
    missing_comp = required_comp - set(composition.columns)
    if missing_comp:
        raise typer.BadParameter(f"composition table missing columns: {sorted(missing_comp)}")

    out = all_islands.copy()
    out["island_id"] = _text(out["island_id"])
    full_n = int(out["island_id"].nunique())

    species = island_species[["island_id", "species"]].copy()
    species["island_id"] = _text(species["island_id"])
    species["species"] = _text(species["species"])
    species = species.loc[species["island_id"].ne("") & species["species"].ne("")]
    flora_counts = (
        species.drop_duplicates(["island_id", "species"])
        .groupby("island_id")["species"]
        .nunique()
        .rename("n_flora_species_recorded")
    )
    out = out.merge(flora_counts, on="island_id", how="left", validate="one_to_one")
    out["n_flora_species_recorded"] = (
        pd.to_numeric(out["n_flora_species_recorded"], errors="coerce").fillna(0).astype(int)
    )
    out["flora_recorded"] = out["n_flora_species_recorded"].gt(0)

    comp = composition.loc[composition["stratum"].astype(str).eq("all_native")].copy()
    comp["trials"] = pd.to_numeric(comp["trials"], errors="coerce").fillna(0)
    comp = comp.loc[comp["trials"].gt(0)]
    comp["island_id"] = _text(comp["island_id"])
    comp["trait_name"] = _text(comp["trait_name"])

    domain_columns: list[str] = []
    for label, trait_name in config["trait_domains"].items():
        column = f"direct_trait_{label}"
        domain_columns.append(column)
        ids = set(comp.loc[comp["trait_name"].eq(str(trait_name)), "island_id"])
        out[column] = out["island_id"].isin(ids)

    any_ids = set(comp["island_id"])
    out["direct_trait_any"] = out["island_id"].isin(any_ids)
    out["n_direct_trait_domains"] = out[domain_columns].sum(axis=1).astype(int)
    out["direct_trait_all_three"] = out[domain_columns].all(axis=1)

    out["observation_stage"] = np.select(
        [
            ~out["flora_recorded"],
            out["flora_recorded"] & ~out["direct_trait_any"],
            out["direct_trait_any"] & ~out["direct_trait_all_three"],
            out["direct_trait_all_three"],
        ],
        [
            "no_flora_record",
            "flora_recorded_no_direct_trait",
            "direct_trait_partial",
            "direct_trait_all_three",
        ],
        default="unresolved",
    )
    if int(out["island_id"].nunique()) != full_n:
        raise RuntimeError("observation audit changed the full island universe")
    return out


def _design_matrix(
    data: pd.DataFrame,
    config: dict[str, Any],
    *,
    add_flora_support: bool,
) -> tuple[np.ndarray, list[str]]:
    baseline = [str(x) for x in config["baseline_covariates"]]
    distance = str(config["distance_gradient"])
    context_column = str(config["context_column"])
    reference = str(config["reference_context"])
    contexts = [str(x) for x in config["contexts"]]

    names = ["intercept"]
    columns = [np.ones(len(data), dtype=float)]

    scale_columns = [*baseline, distance]
    if add_flora_support:
        data = data.copy()
        data["log1p_flora_species_support"] = np.log1p(
            pd.to_numeric(data["n_flora_species_recorded"], errors="coerce").fillna(0)
        )
        scale_columns.append("log1p_flora_species_support")

    scaled: dict[str, np.ndarray] = {}
    for column in scale_columns:
        x = pd.to_numeric(data[column], errors="coerce").to_numpy(float)
        mean = float(np.nanmean(x))
        sd = float(np.nanstd(x))
        if not np.isfinite(sd) or sd <= 0:
            raise typer.BadParameter(f"constant or invalid observation predictor: {column}")
        z = (x - mean) / sd
        scaled[column] = z
        names.append(f"z_{column}")
        columns.append(z)

    context_values = data[context_column].astype(str)
    for context in contexts:
        if context == reference:
            continue
        indicator = context_values.eq(context).to_numpy(float)
        names.append(f"context[{context}]")
        columns.append(indicator)
        names.append(f"z_{distance}:context[{context}]")
        columns.append(scaled[distance] * indicator)

    return np.column_stack(columns), names


def _fit_selection_endpoint(
    support: pd.DataFrame,
    config: dict[str, Any],
    *,
    endpoint: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if endpoint == "flora_recorded":
        work = support.copy()
        outcome = "flora_recorded"
        add_flora_support = False
    elif endpoint == "direct_trait_any":
        work = support.copy()
        outcome = "direct_trait_any"
        add_flora_support = False
    elif endpoint == "direct_trait_any_given_flora":
        work = support.loc[support["flora_recorded"]].copy()
        outcome = "direct_trait_any"
        add_flora_support = True
    elif endpoint == "direct_trait_all_three_given_flora":
        work = support.loc[support["flora_recorded"]].copy()
        outcome = "direct_trait_all_three"
        add_flora_support = True
    else:
        raise typer.BadParameter(f"unknown observation endpoint: {endpoint}")

    needed = [
        *config["baseline_covariates"],
        config["distance_gradient"],
        config["context_column"],
        config["cluster_column"],
    ]
    work = work.dropna(subset=needed).copy()
    y = work[outcome].astype(bool).astype(float).to_numpy()
    if len(work) == 0 or np.all(y == y[0]):
        return pd.DataFrame(), {
            "endpoint": endpoint,
            "status": "not_fit_constant_or_empty",
            "n_islands": int(len(work)),
            "n_positive": int(y.sum()) if len(y) else 0,
        }

    X, names = _design_matrix(work, config, add_flora_support=add_flora_support)
    coef, fit, _ = _fit_grouped_binomial_design(
        y,
        np.ones(len(y), dtype=float),
        X,
        names,
        work[str(config["cluster_column"])].astype(str).to_numpy(),
    )
    coef.insert(0, "endpoint", endpoint)
    coef["odds_ratio"] = np.exp(coef["estimate_log_odds"])
    fit = {
        "endpoint": endpoint,
        "status": "fit",
        "n_islands": int(len(work)),
        "n_positive": int(y.sum()),
        "positive_fraction": float(y.mean()),
        **fit,
    }
    return coef, fit


def run_observation_bias_audit(
    all_islands: pd.DataFrame,
    island_species: pd.DataFrame,
    composition: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    support = _build_support_table(all_islands, island_species, composition, config)

    context = str(config["context_column"])
    context_summary = (
        support.groupby(context, dropna=False)
        .agg(
            n_all_islands=("island_id", "nunique"),
            n_flora_recorded=("flora_recorded", "sum"),
            n_direct_trait_any=("direct_trait_any", "sum"),
            n_direct_trait_all_three=("direct_trait_all_three", "sum"),
            median_flora_species_recorded=("n_flora_species_recorded", "median"),
        )
        .reset_index()
    )
    for numerator in ["n_flora_recorded", "n_direct_trait_any", "n_direct_trait_all_three"]:
        context_summary[f"fraction_{numerator.removeprefix('n_')}"] = (
            context_summary[numerator] / context_summary["n_all_islands"]
        )

    coefficients: list[pd.DataFrame] = []
    fits: list[dict[str, Any]] = []
    for endpoint in config["selection_models"]:
        coef, fit = _fit_selection_endpoint(support, config, endpoint=str(endpoint))
        if not coef.empty:
            coefficients.append(coef)
        fits.append(fit)
    coefficient_table = (
        pd.concat(coefficients, ignore_index=True) if coefficients else pd.DataFrame()
    )
    fit_table = pd.DataFrame(fits)
    return support, context_summary, coefficient_table, fit_table


@app.command("run")
def run(
    all_islands_csv: Path = typer.Option(..., exists=True),
    island_species_csv: Path = typer.Option(..., exists=True),
    composition_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_observation_bias.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    support, context_summary, coefficients, fits = run_observation_bias_audit(
        pd.read_csv(all_islands_csv),
        pd.read_csv(island_species_csv),
        pd.read_csv(composition_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    support.to_csv(output_dir / "all_island_observation_support.csv.gz", index=False)
    context_summary.to_csv(output_dir / "observation_support_by_context.csv", index=False)
    coefficients.to_csv(output_dir / "observation_selection_coefficients.csv", index=False)
    fits.to_csv(output_dir / "observation_selection_fit.csv", index=False)

    manifest = {
        "contract": "chapter1_observation_bias_v1",
        "n_full_islands": int(support["island_id"].nunique()),
        "n_no_flora_record": int((~support["flora_recorded"]).sum()),
        "n_flora_recorded": int(support["flora_recorded"].sum()),
        "n_direct_trait_any": int(support["direct_trait_any"].sum()),
        "n_direct_trait_all_three": int(support["direct_trait_all_three"].sum()),
        "no_record_is_biological_zero": False,
        "distance_gradient_interpretation": (
            "dispersal limitation plus mainland/source-pool accessibility"
        ),
        "pollinator_predictors": False,
    }
    (output_dir / "chapter1_observation_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
