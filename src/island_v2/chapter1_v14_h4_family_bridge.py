"""Post-hoc H4 bridge aligned to the two Chapter 1 v14 H2 response families."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_global_distance import (
    MEASUREMENT_COLUMNS,
    _clustered_wls,
    _context_dummies,
    _lower,
    _measurement_dummies,
    _p2,
    _truthy,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_v14_h4_family_bridge_v1"


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected v14 H4 family-bridge contract")
    if value.get("inferential_role") != "posthoc_atomic_family_reconstruction_sensitivity":
        raise typer.BadParameter("v14 family bridge must remain post-hoc")
    return value


def build_family_scores(
    preflight_rows: pd.DataFrame,
    *,
    components: list[str],
    score_name: str,
    minimum_components: int,
    weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    required = {"species_key", "trait", "trait_state"}
    if missing := required - set(preflight_rows.columns):
        raise typer.BadParameter(f"preflight trait rows missing columns: {sorted(missing)}")

    work = preflight_rows.loc[
        preflight_rows["trait"].astype(str).isin(components),
        ["species_key", "trait", "trait_state"],
    ].copy()
    work["trait_state"] = pd.to_numeric(work["trait_state"], errors="coerce")
    work = work.dropna(subset=["species_key", "trait", "trait_state"])

    conflicts = (
        work.groupby(["species_key", "trait"])["trait_state"]
        .nunique()
        .gt(1)
    )
    if bool(conflicts.any()):
        raise typer.BadParameter("family score has conflicting frozen states")

    wide = (
        work.drop_duplicates(["species_key", "trait"])
        .pivot(index="species_key", columns="trait", values="trait_state")
        .reindex(columns=components)
    )
    wide["n_components"] = wide[components].notna().sum(axis=1)
    score_weights = {
        component: float((weights or {}).get(component, 1.0))
        for component in components
    }
    numerator = sum(
        wide[component].fillna(0.0) * score_weights[component]
        for component in components
    )
    denominator = sum(
        wide[component].notna().astype(float) * score_weights[component]
        for component in components
    )
    wide[score_name] = numerator / denominator
    wide = wide.loc[wide["n_components"].ge(int(minimum_components))].copy()
    return wide.reset_index()


def unique_effect_rows(rows: pd.DataFrame) -> pd.DataFrame:
    required = {
        "row_id",
        "species_key",
        "site_key",
        "study_key",
        "analysis_regime",
        "z_distance",
        "PL_Effect_Size",
        *MEASUREMENT_COLUMNS,
    }
    if missing := required - set(rows.columns):
        raise typer.BadParameter(f"matched effect rows missing columns: {sorted(missing)}")

    check_cols = [
        "species_key",
        "site_key",
        "study_key",
        "analysis_regime",
        "z_distance",
        "PL_Effect_Size",
        *MEASUREMENT_COLUMNS,
    ]
    duplicated = rows.loc[rows["row_id"].duplicated(keep=False)].copy()
    if not duplicated.empty:
        disagreement = (
            duplicated.groupby("row_id")[check_cols]
            .nunique(dropna=False)
            .gt(1)
            .any(axis=1)
        )
        if bool(disagreement.any()):
            raise typer.BadParameter("duplicated row_id has inconsistent outcome metadata")
    return rows[["row_id", *check_cols]].drop_duplicates("row_id").reset_index(drop=True)


def aggregate_family_cells(
    outcome_rows: pd.DataFrame,
    family_scores: pd.DataFrame,
    *,
    score_name: str,
    publication_total_weight: float,
) -> pd.DataFrame:
    joined = outcome_rows.merge(
        family_scores[["species_key", score_name]],
        on="species_key",
        how="inner",
        validate="many_to_one",
    )
    group_cols = [
        "study_key",
        "site_key",
        "species_key",
        "analysis_regime",
        "z_distance",
        score_name,
        *MEASUREMENT_COLUMNS,
    ]
    cells = (
        joined.groupby(group_cols, as_index=False, dropna=False)
        .agg(
            PL_Effect_Size=("PL_Effect_Size", "mean"),
            n_effect_rows=("PL_Effect_Size", "size"),
        )
        .reset_index(drop=True)
    )
    counts = cells.groupby("study_key")["site_key"].transform("size").astype(float)
    cells["analysis_weight"] = float(publication_total_weight) / counts
    return cells


def fit_global_family_score(cells: pd.DataFrame, *, score_name: str) -> dict[str, Any]:
    if cells.empty or cells[score_name].nunique() < 2:
        return {"evaluable": False, "reason": "family_score_has_no_variation"}

    distance = pd.to_numeric(cells["z_distance"], errors="coerce").to_numpy(float)
    score = pd.to_numeric(cells[score_name], errors="coerce").to_numpy(float)
    context_cols, context_names, _ = _context_dummies(cells)
    measure_cols, measure_names = _measurement_dummies(cells)
    columns = [
        np.ones(len(cells), dtype=float),
        *context_cols,
        distance,
        score,
        *measure_cols,
    ]
    names = [
        "intercept",
        *context_names,
        "z_distance",
        score_name,
        *measure_names,
    ]
    fit = _clustered_wls(cells, np.column_stack(columns), names)
    if not fit.get("evaluable"):
        return dict(fit)

    index = fit["names"].index(score_name)
    estimate = float(fit["beta"][index])
    se = float(fit["se"][index])
    z_value = estimate / se if se > 0 else float("nan")
    return {
        "evaluable": True,
        "estimate": estimate,
        "se": se,
        "z": z_value,
        "two_sided_p": _p2(z_value),
        "one_sided_negative_p": _lower(z_value),
        "n_cells": int(len(cells)),
        "n_publications": int(cells["study_key"].nunique()),
        "n_sites": int(cells["site_key"].nunique()),
        "n_species": int(cells["species_key"].nunique()),
    }


def analyse_family(
    preflight_rows: pd.DataFrame,
    matched_effect_rows: pd.DataFrame,
    *,
    family: str,
    family_config: dict[str, Any],
    analysis_config: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    components = [str(x) for x in family_config["components"]]
    score_name = str(family_config["score_name"])
    minimum_components = int(analysis_config["minimum_nonmissing_components_per_species"])
    scores = build_family_scores(
        preflight_rows,
        components=components,
        score_name=score_name,
        minimum_components=minimum_components,
        weights={
            str(k): float(v)
            for k, v in family_config.get("weights", {}).items()
        },
    )
    outcomes = unique_effect_rows(matched_effect_rows)
    cells = aggregate_family_cells(
        outcomes,
        scores,
        score_name=score_name,
        publication_total_weight=float(analysis_config["publication_total_weight"]),
    )

    result_rows: list[dict[str, Any]] = []
    specifications = {
        "primary": cells,
        "supplemental_only": cells.loc[
            cells["PL_Effect_Size_Type2"].astype(str).eq("Sup")
        ].copy(),
        "no_zero_constant": cells.loc[~_truthy(cells["Constant_added"])].copy(),
    }
    for analysis, part in specifications.items():
        result_rows.append(
            {
                "family": family,
                "score_name": score_name,
                "analysis": analysis,
                "inferential_role": "posthoc_atomic_family_reconstruction_sensitivity",
                **fit_global_family_score(part, score_name=score_name),
            }
        )

    support = {
        "family": family,
        "score_name": score_name,
        "components": components,
        "weights": {
            str(k): float(v)
            for k, v in family_config.get("weights", {}).items()
        },
        "role": str(family_config.get("role", "family_score")),
        "minimum_components": minimum_components,
        "n_score_species": int(len(scores)),
        "n_three_component_species": int(scores["n_components"].eq(len(components)).sum()),
        "score_mean": float(scores[score_name].mean()),
        "score_sd": float(scores[score_name].std(ddof=1)),
        "n_low_score_le_one_third": int(scores[score_name].le(1.0 / 3.0).sum()),
        "n_high_score_ge_two_thirds": int(scores[score_name].ge(2.0 / 3.0).sum()),
    }
    return result_rows, support


def run_family_bridge(
    reproductive_preflight: pd.DataFrame,
    reproductive_effect_rows: pd.DataFrame,
    architecture_preflight: pd.DataFrame,
    architecture_effect_rows: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    parents = {
        "reproductive_assurance": (reproductive_preflight, reproductive_effect_rows),
        "floral_architecture": (architecture_preflight, architecture_effect_rows),
    }
    results: list[dict[str, Any]] = []
    support: dict[str, Any] = {}

    for family, family_config in config["families"].items():
        parent_name = str(family_config["parent"])
        preflight_rows, effect_rows = parents[parent_name]
        rows, family_support = analyse_family(
            preflight_rows,
            effect_rows,
            family=str(family),
            family_config=family_config,
            analysis_config=config["analysis"],
        )
        results.extend(rows)
        support[str(family)] = family_support

    frame = pd.DataFrame(results)
    manifest = {
        "contract": CONTRACT,
        "inferential_role": "posthoc_atomic_family_reconstruction_sensitivity",
        "support": support,
        "all_primary_estimates_negative": bool(
            frame.loc[frame["analysis"].eq("primary"), "estimate"].lt(0).all()
        ),
        "historical_selection_identified": False,
        "mediation_identified": False,
        "prospective_replication": False,
    }
    return frame, manifest


@app.command("analyse")
def analyse(
    reproductive_preflight_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    reproductive_effect_rows_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    architecture_preflight_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    architecture_effect_rows_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    results, manifest = run_family_bridge(
        pd.read_csv(reproductive_preflight_csv),
        pd.read_csv(reproductive_effect_rows_csv),
        pd.read_csv(architecture_preflight_csv),
        pd.read_csv(architecture_effect_rows_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "h4_family_bridge_results.csv", index=False)
    (output_dir / "h4_family_bridge_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
