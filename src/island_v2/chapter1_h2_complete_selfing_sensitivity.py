"""Strict complete-observation sensitivity for Chapter 1 H2.

Rebuild Direct-only selfing_core using only species with all three reproductive
components observed and informative:
- self_incompatibility
- mating_system
- autonomous_selfing_capacity

Then refit generalized_accessible ~ isolation + complete selfing_core + controls
on the corrected geography. This is a missing-component reliability sensitivity,
not a formal errors-in-variables correction.
"""
from __future__ import annotations

import json
from functools import reduce
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_v14_h2_decomposition import _bh, _clustered_ols

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})


def _tokens(value: object) -> set[str]:
    return {x.strip() for x in str(value or "").split("|") if x.strip()}


def _signed(value: object, preferred: set[str], opposed: set[str]) -> float:
    tokens = _tokens(value)
    if not tokens:
        return float("nan")
    if tokens <= preferred:
        return 1.0
    if tokens <= opposed:
        return -1.0
    return float("nan")


def build_species_score(
    state_audit: pd.DataFrame,
    syndrome_spec: dict[str, Any],
    *,
    minimum_informative_traits: int,
    require_all: bool,
) -> pd.DataFrame:
    audit = state_audit.loc[_truthy(state_audit["resolved_for_primary"])].copy()
    audit["accepted_species"] = audit["accepted_species"].astype(str)
    audit["trait_name"] = audit["trait_name"].astype(str)

    parts: list[pd.DataFrame] = []
    weights: dict[str, float] = {}
    for trait, spec in syndrome_spec["traits"].items():
        trait = str(trait)
        weight = float(spec["weight"])
        preferred = {str(x) for x in spec["preferred"]}
        opposed = {str(x) for x in spec["opposed"]}
        part = audit.loc[
            audit["trait_name"].eq(trait),
            ["accepted_species", "canonical_signature"],
        ].drop_duplicates("accepted_species")
        part = part.copy()
        part[trait] = [
            _signed(v, preferred, opposed) for v in part["canonical_signature"]
        ]
        parts.append(part[["accepted_species", trait]])
        weights[trait] = weight

    wide = reduce(
        lambda left, right: left.merge(right, on="accepted_species", how="outer"),
        parts,
    )
    traits = list(weights)
    values: list[float] = []
    informative: list[int] = []
    for row in wide.itertuples(index=False):
        record = row._asdict()
        num = 0.0
        den = 0.0
        n = 0
        for trait in traits:
            value = record.get(trait)
            if value is not None and pd.notna(value):
                weight = weights[trait]
                num += weight * float(value)
                den += weight
                n += 1
        informative.append(n)
        eligible = n == len(traits) if require_all else n >= minimum_informative_traits
        values.append((num / den + 1.0) / 2.0 if eligible and den > 0 else float("nan"))

    wide["n_informative_traits"] = informative
    wide["score"] = values
    return wide[["accepted_species", "score", "n_informative_traits"]]


def aggregate_island_score(
    flora: pd.DataFrame,
    species_scores: pd.DataFrame,
    *,
    score_name: str,
) -> pd.DataFrame:
    work = flora[["island_id", "accepted_species"]].copy()
    work["island_id"] = work["island_id"].astype(str)
    work["accepted_species"] = work["accepted_species"].astype(str)
    scores = species_scores.dropna(subset=["score"]).copy()
    scores["accepted_species"] = scores["accepted_species"].astype(str)
    joined = work.merge(scores[["accepted_species", "score"]], on="accepted_species", how="inner")
    joined = joined.drop_duplicates(["island_id", "accepted_species"])
    result = (
        joined.groupby("island_id", as_index=False)
        .agg(**{score_name: ("score", "mean"), f"{score_name}_n_species": ("accepted_species", "nunique")})
    )
    return result


def fit_models(
    flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    syndrome_config: dict[str, Any],
    h2_config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    minimum = int(syndrome_config["score_definition"]["minimum_informative_traits"])
    self_spec = syndrome_config["syndromes"]["selfing_core"]
    access_spec = syndrome_config["syndromes"]["generalized_accessible"]

    self_baseline_species = build_species_score(
        state_audit, self_spec,
        minimum_informative_traits=minimum,
        require_all=False,
    )
    self_complete_species = build_species_score(
        state_audit, self_spec,
        minimum_informative_traits=minimum,
        require_all=True,
    )
    access_species = build_species_score(
        state_audit, access_spec,
        minimum_informative_traits=minimum,
        require_all=False,
    )

    self_baseline = aggregate_island_score(
        flora, self_baseline_species, score_name="selfing_core_direct"
    )
    self_complete = aggregate_island_score(
        flora, self_complete_species, score_name="selfing_core_complete"
    )
    access = aggregate_island_score(
        flora, access_species, score_name="generalized_accessible_direct"
    )

    geography = str(h2_config["geography_column"])
    context = str(h2_config["context_column"])
    cluster = str(h2_config["cluster_column"])
    baseline_covariates = [str(x) for x in h2_config["baseline_covariates"]]
    needed = ["island_id", geography, context, cluster, *baseline_covariates]

    base = (
        access.merge(self_baseline, on="island_id", how="inner")
        .merge(self_complete, on="island_id", how="left")
        .merge(covariates[needed].drop_duplicates("island_id"), on="island_id", how="left")
    )

    rows: list[dict[str, Any]] = []
    for context_value in [str(x) for x in h2_config["contexts"]]:
        part = base.loc[base[context].astype(str).eq(context_value)].copy()
        for variant, mediator in (
            ("direct_min2_baseline_reconstruction", "selfing_core_direct"),
            ("direct_complete3_selfing", "selfing_core_complete"),
        ):
            result = _clustered_ols(
                part,
                response="generalized_accessible_direct",
                predictors=[geography, mediator, *baseline_covariates],
                cluster_column=cluster,
            )
            row: dict[str, Any] = {
                "variant": variant,
                "context": context_value,
                "response": "generalized_accessible_direct",
                "mediator": mediator,
                "status": result["status"],
                "n_islands": int(result.get("n_islands", 0)),
                "n_clusters": int(result.get("n_clusters", 0)),
            }
            if result["status"] == "fit":
                distance = result["coefficients"][f"z_{geography}"]
                mediation = result["coefficients"][f"z_{mediator}"]
                row.update(
                    distance_estimate=distance["estimate"],
                    distance_se=distance["se"],
                    distance_p=distance["p"],
                    selfing_estimate=mediation["estimate"],
                    selfing_se=mediation["se"],
                    selfing_p=mediation["p"],
                )
            rows.append(row)

    results = pd.DataFrame(rows)
    results["distance_q_within_variant"] = np.nan
    for variant in results["variant"].drop_duplicates():
        mask = results["variant"].eq(variant) & results["status"].eq("fit")
        results.loc[mask, "distance_q_within_variant"] = _bh(
            results.loc[mask, "distance_p"]
        )

    summary = {
        "contract": "chapter1_h2_complete_selfing_sensitivity_v1",
        "interpretation": (
            "Direct-only sensitivity requiring all three selfing_core components to be "
            "observed for every contributing mediator species. This addresses missing-"
            "component attenuation but is not a formal errors-in-variables correction."
        ),
        "n_direct_species_selfing_min2": int(self_baseline_species["score"].notna().sum()),
        "n_direct_species_selfing_complete3": int(self_complete_species["score"].notna().sum()),
        "n_direct_species_accessibility_min2": int(access_species["score"].notna().sum()),
        "n_islands_selfing_complete3": int(self_complete["island_id"].nunique()),
        "n_islands_accessibility": int(access["island_id"].nunique()),
    }
    return results, summary


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    h2_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    syndrome_config = yaml.safe_load(syndrome_config_path.read_text(encoding="utf-8"))
    h2_config = yaml.safe_load(h2_config_path.read_text(encoding="utf-8"))
    results, summary = fit_models(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        pd.read_csv(covariates_csv),
        syndrome_config,
        h2_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "h2_complete_selfing_sensitivity.csv", index=False)
    (output_dir / "h2_complete_selfing_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
