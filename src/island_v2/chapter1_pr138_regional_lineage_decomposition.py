"""Decompose PR138 regional plant-response branches into genus assembly and residual parts.

This module operates on the frozen source-adjusted PR138 branch inputs and frozen
observation-selection weights. It asks a stricter question than the observed branch test:
how much of the biogeographic response geometry is reproduced by genus composition, and
does a regional response-vector difference remain after subtracting that expectation?

For each of the two primary plant-side axes, every species is replaced by the mean score
of scored species in its genus. Island membership, floristic status, GIFT mainland
membership, source assignment, score missingness, geography, climate and spatial blocks
remain fixed. Mainland source expectations are recomputed from the genus-mean scores.

Observed source-adjusted scores are never reconstructed: the exact frozen upstream scores
are used, so the observed estimand can be checked bit-for-bit against the existing joint
source+selection workflow. The genus component and observed-minus-genus residual are then
fit with the same frozen IPW weights and the same context-interaction design.

This is an assembly decomposition. A residual is not proof of pollinator causation, and a
genus component does not identify why those lineages colonized or persisted.
"""

from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.chapter1_pr138_source_pool_sensitivity import (
    build_source_expectations,
    match_gift_species_to_frozen_scores,
)
from island_v2.chapter1_pr138_syndrome_analysis import _bh, _joint_wald
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXIS_COMPONENTS = {
    "accessibility_generalization": "generalized_accessible",
    "reproductive_assurance": "selfing_core",
}
PRIMARY_STRATA = ("all_native", "native_nonendemic")
FOCAL_CONTRASTS = {
    "analysis_regime": ("northern_midlatitude", "tropical"),
    "biogeographic_realm": ("Palearctic", "Neotropical"),
}
OUTCOMES = ("observed", "genus_expected", "residual")


def _standardize(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (values - mean) / sd


def _genus_species_scores(species_scores: pd.DataFrame, component: str) -> pd.DataFrame:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    missing = required - set(species_scores.columns)
    if missing:
        raise ValueError(f"species scores missing columns: {sorted(missing)}")
    work = species_scores.loc[
        species_scores["syndrome"].astype(str).eq(component),
        ["accepted_species", "syndrome_concordance"],
    ].drop_duplicates("accepted_species").copy()
    work["score"] = pd.to_numeric(work["syndrome_concordance"], errors="coerce")
    work = work.dropna(subset=["score"])
    work["accepted_species"] = work["accepted_species"].astype(str)
    work["genus"] = work["accepted_species"].str.split().str[0]
    work = work.loc[work["genus"].ne("")].copy()
    work["genus_expected_score"] = work.groupby("genus")["score"].transform("mean")
    return work[["accepted_species", "genus", "score", "genus_expected_score"]]


def _island_genus_expectation(
    status_flora: pd.DataFrame,
    genus_scores: pd.DataFrame,
    *,
    stratum: str,
) -> pd.DataFrame:
    subset = status_flora.loc[
        stratum_mask(status_flora, stratum), ["island_id", "accepted_species"]
    ].drop_duplicates()
    merged = subset.merge(
        genus_scores[["accepted_species", "genus_expected_score"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    return (
        merged.groupby("island_id", as_index=False)
        .agg(
            island_genus_expected=("genus_expected_score", "mean"),
            island_genus_expected_n_species=("accepted_species", "nunique"),
        )
    )


def _mainland_genus_scores(
    matched_full_scores: pd.DataFrame,
    genus_scores: pd.DataFrame,
    *,
    axis: str,
    minimum_species: int,
) -> pd.DataFrame:
    matched = matched_full_scores.loc[
        matched_full_scores["syndrome"].astype(str).eq(AXIS_COMPONENTS[axis]),
        ["entity_ID", "accepted_species"],
    ].drop_duplicates()
    merged = matched.merge(
        genus_scores[["accepted_species", "genus_expected_score"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    out = (
        merged.groupby("entity_ID", as_index=False)
        .agg(
            mainland_syndrome_score=("genus_expected_score", "mean"),
            n_trait_scored_species=("accepted_species", "nunique"),
        )
    )
    out["syndrome"] = axis
    out["source_score_eligible"] = out["n_trait_scored_species"].ge(int(minimum_species))
    return out


def build_regional_lineage_responses(
    source_adjusted_scores: pd.DataFrame,
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    source_assignments: pd.DataFrame,
    source_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return exact observed branch rows plus deterministic genus expectation/residual."""

    minimum_species = int(
        source_config["source_region_scores"][
            "minimum_trait_scored_species_per_region_syndrome"
        ]
    )
    minimum_sources = int(
        source_config["response"]["source_expectation_requires_minimum_source_regions"]
    )
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]

    # Important: matching ambiguity is resolved against the full syndrome-score species
    # universe, exactly as in the frozen source-pool workflow, before filtering to an axis.
    matched_full, match_audit = match_gift_species_to_frozen_scores(gift_flora, species_scores)

    parts: list[pd.DataFrame] = []
    for axis, component in AXIS_COMPONENTS.items():
        genus_scores = _genus_species_scores(species_scores, component)
        mainland = _mainland_genus_scores(
            matched_full,
            genus_scores,
            axis=axis,
            minimum_species=minimum_species,
        )
        source_expectation = build_source_expectations(
            source_assignments,
            mainland,
            min_source_regions=minimum_sources,
        )
        source_expectation = source_expectation.rename(
            columns={"source_expectation": "source_genus_expected"}
        )

        for stratum in PRIMARY_STRATA:
            island_expectation = _island_genus_expectation(
                status_flora,
                genus_scores,
                stratum=stratum,
            )
            observed = source_adjusted_scores.loc[
                source_adjusted_scores["syndrome"].astype(str).eq(component)
                & source_adjusted_scores["stratum"].astype(str).eq(stratum),
                ["source_mode", "island_id", "syndrome_score", "n_species"],
            ].copy()
            observed = observed.rename(columns={"syndrome_score": "observed"})
            observed["axis"] = axis
            observed["stratum"] = stratum
            observed = observed.merge(
                island_expectation,
                on="island_id",
                how="inner",
                validate="many_to_one",
            )
            observed = observed.merge(
                source_expectation.loc[
                    source_expectation["syndrome"].astype(str).eq(axis),
                    ["source_mode", "island_id", "source_genus_expected", "n_source_regions"],
                ],
                on=["source_mode", "island_id"],
                how="inner",
                validate="many_to_one",
            )
            observed = observed.loc[observed["source_mode"].astype(str).isin(source_modes)].copy()
            observed["genus_expected"] = (
                pd.to_numeric(observed["island_genus_expected"], errors="coerce")
                - pd.to_numeric(observed["source_genus_expected"], errors="coerce")
            )
            observed["residual"] = (
                pd.to_numeric(observed["observed"], errors="coerce")
                - observed["genus_expected"]
            )
            parts.append(observed)

    responses = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return responses, match_audit


def _prepare_weighted_axis(
    responses: pd.DataFrame,
    weights: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    axis: str,
    stratum: str,
    source_mode: str,
    selection_mode: str,
    context_layer: str,
    contexts: list[str],
    outcome: str,
    pattern_config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    score = responses.loc[
        responses["axis"].eq(axis)
        & responses["stratum"].eq(stratum)
        & responses["source_mode"].eq(source_mode),
        ["island_id", outcome],
    ].copy()
    weight = weights.loc[
        weights["source_mode"].eq(source_mode)
        & weights["context_layer"].eq(context_layer)
        & weights["stratum"].eq(stratum)
        & weights["syndrome"].eq(axis)
        & weights["selection_mode"].eq(selection_mode),
        ["island_id", "analysis_weight"],
    ].copy()
    needed = ["island_id", context_layer, cluster, geography, *baseline]
    work = score.merge(weight, on="island_id", how="inner", validate="one_to_one")
    work = work.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    work = work.loc[work[context_layer].astype(str).isin(contexts)].copy()
    for column in [outcome, "analysis_weight", geography, *baseline]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cluster] = work[cluster].fillna("").astype(str)
    work = work.dropna(subset=[outcome, "analysis_weight", geography, *baseline])
    return work.loc[work[cluster].ne("")].reset_index(drop=True)


def _fit_context_slope(
    work: pd.DataFrame,
    *,
    outcome: str,
    pattern_config: dict[str, Any],
) -> dict[str, Any]:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    design = np.column_stack(
        [
            np.ones(len(work), dtype=float),
            _standardize(work[geography]),
            *[_standardize(work[x]) for x in baseline],
        ]
    )
    coefficients, _, fit = _fit_weighted_clustered_design(
        work[outcome].to_numpy(float),
        work["analysis_weight"].to_numpy(float),
        design,
        names,
        work[cluster].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return {"status": str(fit.get("status", "fit_failed"))}
    row = coefficients.set_index("predictor").loc[f"z_{geography}"]
    return {
        "status": "fit",
        "distance_slope": float(row["estimate"]),
        "cluster_robust_se": float(row["cluster_robust_se"]),
        "p_value": float(row["p_value"]),
        "n_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
    }


def _fit_between_vector(
    responses: pd.DataFrame,
    weights: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    stratum: str,
    source_mode: str,
    selection_mode: str,
    context_layer: str,
    context_a: str,
    context_b: str,
    outcome: str,
    pattern_config: dict[str, Any],
) -> dict[str, Any]:
    """Exact PR138 `_between_contexts` design for the two primary axes."""

    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    parts: list[pd.DataFrame] = []
    for axis in sorted(AXIS_COMPONENTS):
        work = _prepare_weighted_axis(
            responses,
            weights,
            covariates,
            axis=axis,
            stratum=stratum,
            source_mode=source_mode,
            selection_mode=selection_mode,
            context_layer=context_layer,
            contexts=[context_a, context_b],
            outcome=outcome,
            pattern_config=pattern_config,
        )
        work["axis"] = axis
        parts.append(work)
    stacked = pd.concat(parts, ignore_index=True)
    b_indicator = stacked[context_layer].astype(str).eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    interaction_names: list[str] = []
    for axis in sorted(AXIS_COMPONENTS):
        mask = stacked["axis"].eq(axis).to_numpy()
        indicator = mask.astype(float)
        names.append(f"syndrome[{axis}]")
        columns.append(indicator)
        for predictor in baseline:
            values = np.zeros(len(stacked), dtype=float)
            values[mask] = _standardize(stacked.loc[mask, predictor])
            names.append(f"syndrome[{axis}]:z_{predictor}")
            columns.append(values)
        names.append(f"syndrome[{axis}]:context[{context_b}]")
        columns.append(indicator * b_indicator)
        geography_z = np.zeros(len(stacked), dtype=float)
        geography_z[mask] = _standardize(stacked.loc[mask, geography])
        names.append(f"syndrome[{axis}]:z_{geography}")
        columns.append(geography_z)
        interaction = f"syndrome[{axis}]:z_{geography}:context[{context_b}]"
        names.append(interaction)
        columns.append(geography_z * b_indicator)
        interaction_names.append(interaction)

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        stacked[outcome].to_numpy(float),
        stacked["analysis_weight"].to_numpy(float),
        np.column_stack(columns),
        names,
        stacked[cluster].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return {"status": str(fit.get("status", "fit_failed"))}
    index = coefficients.set_index("predictor")
    positions = [names.index(name) for name in interaction_names]
    vector = np.array([float(index.loc[name, "estimate"]) for name in interaction_names])
    vector_covariance = covariance[np.ix_(positions, positions)]
    statistic, df, p_value = _joint_wald(vector, vector_covariance)
    result: dict[str, Any] = {
        "status": "fit",
        "joint_context_difference_chisq": statistic,
        "joint_context_difference_df": df,
        "p_value": p_value,
        "n_unique_islands": int(stacked["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
    }
    for axis, name in zip(sorted(AXIS_COMPONENTS), interaction_names, strict=True):
        row = index.loc[name]
        result[f"{axis}_interaction"] = float(row["estimate"])
        result[f"{axis}_interaction_se"] = float(row["cluster_robust_se"])
        result[f"{axis}_interaction_p"] = float(row["p_value"])
    return result


def run_regional_lineage_decomposition(
    source_adjusted_scores: pd.DataFrame,
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    source_assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    frozen_weights: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    responses, match_audit = build_regional_lineage_responses(
        source_adjusted_scores,
        species_scores,
        status_flora,
        gift_flora,
        source_assignments,
        source_config,
    )
    source_modes = sorted(responses["source_mode"].dropna().astype(str).unique())
    selection_modes = sorted(frozen_weights["selection_mode"].dropna().astype(str).unique())

    layer_covariates: dict[str, pd.DataFrame] = {"analysis_regime": covariates.copy()}
    layer_covariates["biogeographic_realm"] = covariates.merge(
        realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )

    slope_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for layer, (context_a, context_b) in FOCAL_CONTRASTS.items():
        cov = layer_covariates[layer]
        for stratum in PRIMARY_STRATA:
            for source_mode in source_modes:
                for selection_mode in selection_modes:
                    for axis in sorted(AXIS_COMPONENTS):
                        for context in (context_a, context_b):
                            for outcome in OUTCOMES:
                                work = _prepare_weighted_axis(
                                    responses,
                                    frozen_weights,
                                    cov,
                                    axis=axis,
                                    stratum=stratum,
                                    source_mode=source_mode,
                                    selection_mode=selection_mode,
                                    context_layer=layer,
                                    contexts=[context],
                                    outcome=outcome,
                                    pattern_config=pattern_config,
                                )
                                result = _fit_context_slope(
                                    work,
                                    outcome=outcome,
                                    pattern_config=pattern_config,
                                )
                                slope_rows.append(
                                    {
                                        "context_layer": layer,
                                        "stratum": stratum,
                                        "source_mode": source_mode,
                                        "selection_mode": selection_mode,
                                        "context": context,
                                        "axis": axis,
                                        "outcome": outcome,
                                        **result,
                                    }
                                )
                    row: dict[str, Any] = {
                        "context_layer": layer,
                        "stratum": stratum,
                        "source_mode": source_mode,
                        "selection_mode": selection_mode,
                        "context_a": context_a,
                        "context_b": context_b,
                    }
                    for outcome in OUTCOMES:
                        result = _fit_between_vector(
                            responses,
                            frozen_weights,
                            cov,
                            stratum=stratum,
                            source_mode=source_mode,
                            selection_mode=selection_mode,
                            context_layer=layer,
                            context_a=context_a,
                            context_b=context_b,
                            outcome=outcome,
                            pattern_config=pattern_config,
                        )
                        for key, value in result.items():
                            row[f"{outcome}_{key}"] = value
                    between_rows.append(row)

    slopes = pd.DataFrame(slope_rows)
    between = pd.DataFrame(between_rows)
    between["residual_q_across_source_modes"] = np.nan
    fit = between["residual_status"].eq("fit")
    if fit.any():
        between.loc[fit, "residual_q_across_source_modes"] = (
            between.loc[fit]
            .groupby(["context_layer", "stratum", "selection_mode"], group_keys=False)[
                "residual_p_value"
            ]
            .apply(_bh)
        )
    between["residual_vector_supported"] = (
        pd.to_numeric(between["residual_q_across_source_modes"], errors="coerce") <= 0.05
    ).fillna(False)

    manifest = {
        "contract": "chapter1_pr138_regional_lineage_decomposition_v1",
        "axes": AXIS_COMPONENTS,
        "primary_strata": list(PRIMARY_STRATA),
        "focal_contrasts": FOCAL_CONTRASTS,
        "observed_scores_reconstructed": False,
        "observed_scores_source": "frozen source_adjusted_pathway_island_scores",
        "genus_expectation": "replace each scored species by mean score of scored species in same genus",
        "gift_matching_ambiguity_universe": "full frozen species-syndrome score universe before axis filtering",
        "source_assignments_recomputed": False,
        "source_assignments_outcome_blind": True,
        "source_scores_recomputed_for_genus_expectation": True,
        "selection_weights_recomputed": False,
        "selection_weights_source": "frozen joint source+selection workflow",
        "causal_pollination_claimed": False,
        "claim_ceiling": (
            "Residual regional response-vector heterogeneity beyond genus composition can be "
            "reported. It does not identify pollinator loss, selection, or the process causing "
            "lineage turnover."
        ),
    }
    return {
        "responses": responses,
        "match_audit": match_audit,
        "slopes": slopes,
        "between": between,
        "manifest": manifest,
    }


@app.command("run")
def run(
    source_adjusted_scores_csv: Annotated[Path, typer.Option(exists=True)],
    species_scores_csv: Annotated[Path, typer.Option(exists=True)],
    status_flora_csv: Annotated[Path, typer.Option(exists=True)],
    gift_flora_csv: Annotated[Path, typer.Option(exists=True)],
    source_assignments_csv: Annotated[Path, typer.Option(exists=True)],
    covariates_csv: Annotated[Path, typer.Option(exists=True)],
    realm_assignment_csv: Annotated[Path, typer.Option(exists=True)],
    frozen_weights_csv: Annotated[Path, typer.Option(exists=True)],
    pattern_config_path: Annotated[Path, typer.Option(exists=True)],
    source_config_path: Annotated[Path, typer.Option(exists=True)],
    output_dir: Annotated[Path, typer.Option()],
) -> None:
    pattern = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    source = yaml.safe_load(source_config_path.read_text(encoding="utf-8"))
    outputs = run_regional_lineage_decomposition(
        pd.read_csv(source_adjusted_scores_csv),
        pd.read_csv(species_scores_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(source_assignments_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pd.read_csv(frozen_weights_csv),
        pattern,
        source,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs["responses"].to_csv(output_dir / "regional_lineage_responses.csv.gz", index=False)
    outputs["match_audit"].to_csv(output_dir / "regional_lineage_gift_match_audit.csv", index=False)
    outputs["slopes"].to_csv(output_dir / "regional_lineage_context_slopes.csv", index=False)
    outputs["between"].to_csv(output_dir / "regional_lineage_between_context.csv", index=False)
    manifest = outputs["manifest"]
    (output_dir / "regional_lineage_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
