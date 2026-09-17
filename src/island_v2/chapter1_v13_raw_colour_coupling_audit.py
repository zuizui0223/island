"""Raw flower-colour conditioned floral-architecture coupling audit for Chapter 1 v13.

This analysis does not create pollination-syndrome scores.  It conditions directly on
reported flower colour and asks whether raw floral-form or tube-depth states associated
with predeclared functional pollination architectures are represented among species with
that colour, and whether that coupling changes with island isolation.

The estimand for each combination is therefore P(raw architecture state | raw colour,
architecture trait resolved) at the island-flora level.  It separates colour-frequency
change from colour-architecture coupling.  Labels such as ``bird`` or ``butterfly`` are
interpretation aids inherited from the frozen PR138 architecture definitions, not
observations of realized pollinators.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_v13_colour_architecture_audit import (
    FOCAL_COLOURS,
    parse_axis_traits,
)
from island_v2.chapter1_v13_raw_colour_audit import (
    _bh,
    _fit_one_colour,
    _selfing_table,
)
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)

COUPLING_SPECS: dict[str, dict[str, Any]] = {
    "red_pink__butterfly_form_given_colour": {
        "colour": "red_pink",
        "trait_name": "floral_form",
        "states": {"salverform", "tubular", "spurred", "funnel_trumpet"},
        "architecture_label": "butterfly_associated_form",
    },
    "yellow_orange__butterfly_form_given_colour": {
        "colour": "yellow_orange",
        "trait_name": "floral_form",
        "states": {"salverform", "tubular", "spurred", "funnel_trumpet"},
        "architecture_label": "butterfly_associated_form",
    },
    "blue_purple__butterfly_form_given_colour": {
        "colour": "blue_purple",
        "trait_name": "floral_form",
        "states": {"salverform", "tubular", "spurred", "funnel_trumpet"},
        "architecture_label": "butterfly_associated_form",
    },
    "red_pink__butterfly_deep_tube_given_colour": {
        "colour": "red_pink",
        "trait_name": "tube_depth_class",
        "states": {"deep"},
        "architecture_label": "butterfly_associated_deep_tube",
    },
    "yellow_orange__butterfly_deep_tube_given_colour": {
        "colour": "yellow_orange",
        "trait_name": "tube_depth_class",
        "states": {"deep"},
        "architecture_label": "butterfly_associated_deep_tube",
    },
    "blue_purple__butterfly_deep_tube_given_colour": {
        "colour": "blue_purple",
        "trait_name": "tube_depth_class",
        "states": {"deep"},
        "architecture_label": "butterfly_associated_deep_tube",
    },
    "red_pink__bird_form_given_colour": {
        "colour": "red_pink",
        "trait_name": "floral_form",
        "states": {"tubular", "funnel_trumpet", "bell_campanulate"},
        "architecture_label": "bird_associated_form",
    },
    "yellow_orange__bird_form_given_colour": {
        "colour": "yellow_orange",
        "trait_name": "floral_form",
        "states": {"tubular", "funnel_trumpet", "bell_campanulate"},
        "architecture_label": "bird_associated_form",
    },
    "red_pink__bird_deep_tube_given_colour": {
        "colour": "red_pink",
        "trait_name": "tube_depth_class",
        "states": {"intermediate", "deep"},
        "architecture_label": "bird_associated_intermediate_deep_tube",
    },
    "yellow_orange__bird_deep_tube_given_colour": {
        "colour": "yellow_orange",
        "trait_name": "tube_depth_class",
        "states": {"intermediate", "deep"},
        "architecture_label": "bird_associated_intermediate_deep_tube",
    },
    "blue_purple__large_bee_form_given_colour": {
        "colour": "blue_purple",
        "trait_name": "floral_form",
        "states": {"bell_campanulate", "bilabiate", "papilionaceous", "funnel_trumpet"},
        "architecture_label": "large_bee_associated_form",
    },
    "yellow_orange__large_bee_form_given_colour": {
        "colour": "yellow_orange",
        "trait_name": "floral_form",
        "states": {"bell_campanulate", "bilabiate", "papilionaceous", "funnel_trumpet"},
        "architecture_label": "large_bee_associated_form",
    },
    "blue_purple__large_bee_deep_tube_given_colour": {
        "colour": "blue_purple",
        "trait_name": "tube_depth_class",
        "states": {"intermediate", "deep"},
        "architecture_label": "large_bee_associated_intermediate_deep_tube",
    },
    "yellow_orange__large_bee_deep_tube_given_colour": {
        "colour": "yellow_orange",
        "trait_name": "tube_depth_class",
        "states": {"intermediate", "deep"},
        "architecture_label": "large_bee_associated_intermediate_deep_tube",
    },
}


def _species_raw_states(
    species_axis: pd.DataFrame,
    *,
    evidence_scope: str,
) -> pd.DataFrame:
    required = {"accepted_species", "axis", "trait_composition", "quality"}
    missing = required - set(species_axis.columns)
    if missing:
        raise typer.BadParameter(f"species-axis table missing columns: {sorted(missing)}")
    work = species_axis.loc[
        species_axis["axis"].astype(str).isin({"flower_colour", "floral_structural_complexity"})
    ].copy()
    if evidence_scope == "direct":
        work = work.loc[
            work["quality"].fillna("").astype(str).str.lower().isin({"high", "medium"})
        ].copy()
    elif evidence_scope != "all":
        raise typer.BadParameter("evidence_scope must be 'all' or 'direct'")

    rows: list[dict[str, Any]] = []
    for species, group in work.groupby("accepted_species", sort=False):
        colour_states: set[str] = set()
        form_states: set[str] = set()
        tube_states: set[str] = set()
        for row in group.itertuples(index=False):
            parsed = parse_axis_traits(row.trait_composition)
            if str(row.axis) == "flower_colour":
                colour_states |= parsed.get("flower_primary_color", set())
            elif str(row.axis) == "floral_structural_complexity":
                form_states |= parsed.get("floral_form", set())
                tube_states |= parsed.get("tube_depth_class", set())
        colour_states &= FOCAL_COLOURS
        if not colour_states:
            continue
        rows.append(
            {
                "accepted_species": str(species),
                "colour_states": colour_states,
                "floral_form": form_states,
                "tube_depth_class": tube_states,
            }
        )
    return pd.DataFrame(
        rows,
        columns=["accepted_species", "colour_states", "floral_form", "tube_depth_class"],
    )


def build_colour_conditioned_architecture_counts(
    species_axis: pd.DataFrame,
    status_flora: pd.DataFrame,
    *,
    evidence_scope: str,
    strata: list[str],
) -> pd.DataFrame:
    """Build P(architecture | raw colour, architecture trait resolved) island counts."""
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
        "floristic_status",
    }
    missing = required_flora - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")

    states = _species_raw_states(species_axis, evidence_scope=evidence_scope)
    if states.empty:
        return pd.DataFrame(
            columns=["island_id", "stratum", "combination", "successes", "trials", "share"]
        )
    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    joined = flora.merge(states, on="accepted_species", how="inner", validate="many_to_one")

    rows: list[dict[str, Any]] = []
    for stratum in strata:
        subset = joined.loc[stratum_mask(joined, stratum)].copy()
        subset = subset.drop_duplicates(["island_id", "accepted_species"])
        for island_id, part in subset.groupby("island_id", sort=False):
            for name, spec in COUPLING_SPECS.items():
                colour = str(spec["colour"])
                trait_name = str(spec["trait_name"])
                target_states = set(spec["states"])
                eligible = part.loc[
                    part["colour_states"].map(lambda values: colour in values)
                    & part[trait_name].map(bool)
                ].copy()
                trials = int(eligible["accepted_species"].nunique())
                if trials <= 0:
                    continue
                successes = int(
                    sum(bool(target_states & states) for states in eligible[trait_name])
                )
                rows.append(
                    {
                        "island_id": str(island_id),
                        "stratum": str(stratum),
                        "combination": name,
                        "successes": successes,
                        "trials": trials,
                        "share": successes / trials,
                    }
                )
    return pd.DataFrame(rows)


def fit_colour_conditioned_architecture_models(
    coupling_counts: pd.DataFrame,
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(config["geography_column"])
    context_column = str(config["context_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    contexts = [str(x) for x in config["contexts"]]
    strata = [str(x) for x in config["strata"]]
    tiers = {str(k): int(v) for k, v in config["support_tiers"].items()}

    required = {"island_id", "stratum", "combination", "successes", "trials"}
    missing = required - set(coupling_counts.columns)
    if missing:
        raise typer.BadParameter(f"coupling counts missing columns: {sorted(missing)}")
    needed_cov = ["island_id", geography, context_column, cluster, *baseline]
    missing = set(needed_cov) - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")

    data = coupling_counts.merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    ).merge(
        _selfing_table(island_scores),
        on=["island_id", "stratum"],
        how="left",
        validate="many_to_one",
    )
    for column in ["successes", "trials", "selfing_core", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context_column] = data[context_column].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)

    rows: list[dict[str, Any]] = []
    for stratum in strata:
        for context in contexts:
            for combination, spec in COUPLING_SPECS.items():
                base = data.loc[
                    data["stratum"].astype(str).eq(stratum)
                    & data[context_column].eq(context)
                    & data["combination"].astype(str).eq(combination)
                ].copy()
                for tier, threshold in tiers.items():
                    for model, conditional in (
                        ("unconditional", False),
                        ("conditional_selfing", True),
                    ):
                        required_columns = ["successes", "trials", geography, *baseline, cluster]
                        if conditional:
                            required_columns.append("selfing_core")
                        complete = base.dropna(subset=required_columns).copy()
                        n_islands = int(complete["island_id"].nunique())
                        if n_islands < threshold:
                            rows.append(
                                {
                                    "stratum": stratum,
                                    "context": context,
                                    "combination": combination,
                                    "architecture_label": spec["architecture_label"],
                                    "trait_name": spec["trait_name"],
                                    "support_tier": tier,
                                    "threshold": threshold,
                                    "model": model,
                                    "status": "not_testable",
                                    "n_unique_islands": n_islands,
                                }
                            )
                            continue
                        coef, fit = _fit_one_colour(
                            complete,
                            geography=geography,
                            baseline=baseline,
                            cluster=cluster,
                            conditional_selfing=conditional,
                        )
                        indexed = coef.set_index("predictor")
                        distance = indexed.loc[f"z_{geography}"]
                        row: dict[str, Any] = {
                            "stratum": stratum,
                            "context": context,
                            "combination": combination,
                            "architecture_label": spec["architecture_label"],
                            "trait_name": spec["trait_name"],
                            "support_tier": tier,
                            "threshold": threshold,
                            "model": model,
                            "status": "fit",
                            "n_unique_islands": int(fit["n_unique_islands"]),
                            "n_clusters": int(fit["n_clusters"]),
                            "median_trials": float(fit["median_trials"]),
                            "distance_estimate": float(distance["estimate_log_odds"]),
                            "distance_se": float(distance["cluster_robust_se"]),
                            "distance_p": float(distance["p_value"]),
                        }
                        if conditional:
                            selfing = indexed.loc["z_selfing_core"]
                            row.update(
                                {
                                    "selfing_core_estimate": float(selfing["estimate_log_odds"]),
                                    "selfing_core_se": float(selfing["cluster_robust_se"]),
                                    "selfing_core_p": float(selfing["p_value"]),
                                }
                            )
                        rows.append(row)

    result = pd.DataFrame(rows)
    if not result.empty and "distance_p" in result.columns:
        result["distance_q"] = np.nan
        fit_mask = result["status"].eq("fit")
        for _, index in result.loc[fit_mask].groupby(
            ["stratum", "context", "support_tier", "model"]
        ).groups.items():
            result.loc[index, "distance_q"] = _bh(result.loc[index, "distance_p"])
    return result


@app.command("run")
def run(
    species_axis_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    evidence_scope: str = typer.Option("all"),
    stratum: list[str] = typer.Option(["all_observed"], "--stratum"),
) -> None:
    config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    config = dict(config)
    config["strata"] = [str(x) for x in stratum]

    species_axis = pd.read_csv(species_axis_csv, dtype=str)
    status_flora = pd.read_csv(status_flora_csv)
    island_scores = pd.read_csv(island_scores_csv)
    covariates = pd.read_csv(covariates_csv)

    counts = build_colour_conditioned_architecture_counts(
        species_axis,
        status_flora,
        evidence_scope=evidence_scope,
        strata=config["strata"],
    )
    models = fit_colour_conditioned_architecture_models(
        counts,
        island_scores,
        covariates,
        config,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(
        output_dir / "raw_colour_conditioned_architecture_counts.csv.gz",
        index=False,
        compression="gzip",
    )
    models.to_csv(
        output_dir / "raw_colour_conditioned_architecture_model_results.csv",
        index=False,
    )
    manifest = {
        "contract": "chapter1_v13_raw_colour_coupling_audit_v1",
        "estimand": "P(raw architecture state | raw colour, architecture trait resolved)",
        "evidence_scope": evidence_scope,
        "strata": config["strata"],
        "combination_specs": {
            key: {
                "colour": value["colour"],
                "trait_name": value["trait_name"],
                "states": sorted(value["states"]),
                "architecture_label": value["architecture_label"],
            }
            for key, value in COUPLING_SPECS.items()
        },
        "weighted_score_used": False,
        "marginal_colour_frequency_removed_by_conditioning": True,
        "realized_pollinator_identity_claimed": False,
        "causal_mediation_claimed": False,
    }
    (output_dir / "raw_colour_coupling_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
