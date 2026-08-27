"""Area moderation of Chapter 1 source-pool and plant-response branches.

The module tests a predeclared continuous distance-by-area interaction.  It does
not threshold islands into small/large classes and it does not observe genetic
founder effects, pollinator mobility, or effective pollination service.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_global_branching import build_branch_scores
from island_v2.chapter1_pr136_biogeographic_residual import (
    _fit_weighted_clustered_design,
)
from island_v2.chapter1_pr138_syndrome_analysis import _bh, _joint_wald


def _standardize(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (values - mean) / sd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_text(path: Path) -> str:
    """Hash tracked config text after canonical newline normalization."""
    canonical = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _branch_family(
    island_scores: pd.DataFrame,
    branching_config: dict[str, Any],
    family_spec: dict[str, Any],
) -> pd.DataFrame:
    scores = build_branch_scores(island_scores, branching_config)
    responses = [str(x) for x in family_spec["responses"]]
    out = scores.loc[
        scores["syndrome"].isin(responses),
        ["island_id", "stratum", "syndrome", "syndrome_score"],
    ].rename(columns={"syndrome": "response", "syndrome_score": "response_score"})
    out["source_mode"] = "not_applicable"
    return out


def _lineage_family(
    lineage_scores: pd.DataFrame,
    family_spec: dict[str, Any],
) -> pd.DataFrame:
    required = {
        "island_id",
        "stratum",
        "source_mode",
        "evidence_scope",
        "source_matching",
        "minimum_represented_genera",
        *[str(x) for x in family_spec["responses"]],
    }
    missing = required - set(lineage_scores.columns)
    if missing:
        raise ValueError(f"lineage scores missing columns: {sorted(missing)}")
    selected = lineage_scores.loc[
        lineage_scores["evidence_scope"].eq(str(family_spec["evidence_scope"]))
        & lineage_scores["source_matching"].eq(str(family_spec["source_matching"]))
        & lineage_scores["minimum_represented_genera"].eq(
            int(family_spec["minimum_represented_genera"])
        )
    ].copy()
    return selected.melt(
        id_vars=["island_id", "stratum", "source_mode"],
        value_vars=[str(x) for x in family_spec["responses"]],
        var_name="response",
        value_name="response_score",
    )


def build_families(
    island_scores: pd.DataFrame,
    lineage_scores: pd.DataFrame,
    branching_config: dict[str, Any],
    area_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    families: dict[str, pd.DataFrame] = {}
    for family, spec in area_config["families"].items():
        kind = str(spec["kind"])
        if kind == "branch_scores":
            table = _branch_family(island_scores, branching_config, spec)
        elif kind == "lineage_scores":
            table = _lineage_family(lineage_scores, spec)
        else:
            raise ValueError(f"unknown family kind: {kind}")
        table = table.copy()
        table["family"] = str(family)
        families[str(family)] = table
    return families


def _prepare(
    family_data: pd.DataFrame,
    covariates: pd.DataFrame,
    area_config: dict[str, Any],
    *,
    context_column: str,
) -> pd.DataFrame:
    distance = str(area_config["distance_column"])
    area = str(area_config["area_column"])
    cluster = str(area_config["cluster_column"])
    controls = [str(x) for x in area_config["control_columns"]]
    required = {"island_id", distance, area, cluster, context_column, *controls}
    missing = required - set(covariates.columns)
    if missing:
        raise ValueError(f"covariates missing columns: {sorted(missing)}")
    data = family_data.merge(
        covariates[list(required)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["response_score", distance, area, *controls]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context_column] = data[context_column].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    return data


def _fit_within(
    data: pd.DataFrame,
    area_config: dict[str, Any],
    *,
    family: str,
    source_mode: str,
    stratum: str,
    context: str,
    context_column: str,
    context_layer: str,
    support_tier: str,
    threshold: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    distance = str(area_config["distance_column"])
    area = str(area_config["area_column"])
    cluster = str(area_config["cluster_column"])
    controls = [str(x) for x in area_config["control_columns"]]
    needed = ["response_score", distance, area, cluster, *controls]
    work = data.loc[
        data["stratum"].eq(stratum) & data[context_column].eq(context)
    ].dropna(subset=needed).copy()
    counts = work.groupby("response")["island_id"].nunique()
    responses = sorted(counts.loc[counts.ge(threshold)].index.astype(str))
    base = {
        "family": family,
        "source_mode": source_mode,
        "context_layer": context_layer,
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context": context,
        "n_retained_responses": len(responses),
        "retained_responses": "|".join(responses),
    }
    if len(responses) < 2:
        return pd.DataFrame(), {**base, "status": "not_testable"}
    work = work.loc[work["response"].isin(responses)].copy()

    names: list[str] = []
    columns: list[np.ndarray] = []
    term_names: dict[str, dict[str, str]] = {}
    for response in responses:
        mask = work["response"].eq(response).to_numpy()
        indicator = mask.astype(float)
        names.append(f"response[{response}]")
        columns.append(indicator)
        for control in controls:
            values = np.zeros(len(work), dtype=float)
            values[mask] = _standardize(work.loc[mask, control])
            names.append(f"response[{response}]:z_{control}")
            columns.append(values)
        z_distance = np.zeros(len(work), dtype=float)
        z_area = np.zeros(len(work), dtype=float)
        z_distance[mask] = _standardize(work.loc[mask, distance])
        z_area[mask] = _standardize(work.loc[mask, area])
        distance_name = f"response[{response}]:z_distance"
        area_name = f"response[{response}]:z_area"
        interaction_name = f"response[{response}]:z_distance:z_area"
        names.extend([distance_name, area_name, interaction_name])
        columns.extend([z_distance, z_area, z_distance * z_area])
        term_names[response] = {
            "distance": distance_name,
            "area": area_name,
            "distance_x_area": interaction_name,
        }

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["response_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return pd.DataFrame(), {**base, "status": str(fit.get("status", "fit_failed"))}

    indexed = coefficients.set_index("predictor")
    interaction_names = [term_names[x]["distance_x_area"] for x in responses]
    indices = [names.index(x) for x in interaction_names]
    vector = np.array([float(indexed.loc[x, "estimate"]) for x in interaction_names])
    stat, df, p_value = _joint_wald(vector, covariance[np.ix_(indices, indices)])

    rows: list[dict[str, Any]] = []
    for response in responses:
        term = term_names[response]
        d = indexed.loc[term["distance"]]
        a = indexed.loc[term["area"]]
        interaction = indexed.loc[term["distance_x_area"]]
        distance_estimate = float(d["estimate"])
        interaction_estimate = float(interaction["estimate"])
        rows.append(
            {
                **base,
                "response": response,
                "n_islands": int(counts.loc[response]),
                "distance_estimate_at_mean_area": distance_estimate,
                "distance_se": float(d["cluster_robust_se"]),
                "distance_p": float(d["p_value"]),
                "area_estimate_at_mean_distance": float(a["estimate"]),
                "area_se": float(a["cluster_robust_se"]),
                "area_p": float(a["p_value"]),
                "distance_x_area_estimate": interaction_estimate,
                "distance_x_area_se": float(interaction["cluster_robust_se"]),
                "distance_x_area_p": float(interaction["p_value"]),
                "distance_slope_at_small_area_z_minus1": (
                    distance_estimate - interaction_estimate
                ),
                "distance_slope_at_large_area_z_plus1": (
                    distance_estimate + interaction_estimate
                ),
            }
        )
    return pd.DataFrame(rows), {
        **base,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_distance_x_area_chisq": stat,
        "joint_distance_x_area_df": df,
        "p_value": p_value,
    }


def _fit_between(
    data: pd.DataFrame,
    area_config: dict[str, Any],
    *,
    family: str,
    source_mode: str,
    stratum: str,
    context_a: str,
    context_b: str,
    context_column: str,
    context_layer: str,
    support_tier: str,
    threshold: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    distance = str(area_config["distance_column"])
    area = str(area_config["area_column"])
    cluster = str(area_config["cluster_column"])
    controls = [str(x) for x in area_config["control_columns"]]
    needed = ["response_score", distance, area, cluster, *controls]
    work = data.loc[
        data["stratum"].eq(stratum)
        & data[context_column].isin([context_a, context_b])
    ].dropna(subset=needed).copy()
    counts = work.groupby(["response", context_column])["island_id"].nunique().unstack(
        fill_value=0
    )
    for context in [context_a, context_b]:
        if context not in counts.columns:
            counts[context] = 0
    responses = sorted(
        counts.index[
            (counts[context_a] >= threshold) & (counts[context_b] >= threshold)
        ].astype(str)
    )
    base = {
        "family": family,
        "source_mode": source_mode,
        "context_layer": context_layer,
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context_a": context_a,
        "context_b": context_b,
        "n_retained_responses": len(responses),
        "retained_responses": "|".join(responses),
    }
    if len(responses) < 2:
        return pd.DataFrame(), {**base, "status": "not_testable"}
    work = work.loc[work["response"].isin(responses)].copy()
    b_indicator = work[context_column].eq(context_b).to_numpy(float)

    names: list[str] = []
    columns: list[np.ndarray] = []
    difference_names: dict[str, str] = {}
    for response in responses:
        mask = work["response"].eq(response).to_numpy()
        indicator = mask.astype(float)
        names.append(f"response[{response}]")
        columns.append(indicator)
        for control in controls:
            values = np.zeros(len(work), dtype=float)
            values[mask] = _standardize(work.loc[mask, control])
            names.append(f"response[{response}]:z_{control}")
            columns.append(values)
        names.append(f"response[{response}]:context[{context_b}]")
        columns.append(indicator * b_indicator)
        z_distance = np.zeros(len(work), dtype=float)
        z_area = np.zeros(len(work), dtype=float)
        z_distance[mask] = _standardize(work.loc[mask, distance])
        z_area[mask] = _standardize(work.loc[mask, area])
        interaction = z_distance * z_area
        names.extend(
            [
                f"response[{response}]:z_distance",
                f"response[{response}]:z_area",
                f"response[{response}]:z_distance:z_area",
                f"response[{response}]:z_distance:context[{context_b}]",
                f"response[{response}]:z_area:context[{context_b}]",
            ]
        )
        columns.extend(
            [
                z_distance,
                z_area,
                interaction,
                z_distance * b_indicator,
                z_area * b_indicator,
            ]
        )
        difference_name = (
            f"response[{response}]:z_distance:z_area:context[{context_b}]"
        )
        names.append(difference_name)
        columns.append(interaction * b_indicator)
        difference_names[response] = difference_name

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["response_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return pd.DataFrame(), {**base, "status": str(fit.get("status", "fit_failed"))}
    indexed = coefficients.set_index("predictor")
    names_to_test = [difference_names[x] for x in responses]
    indices = [names.index(x) for x in names_to_test]
    vector = np.array([float(indexed.loc[x, "estimate"]) for x in names_to_test])
    stat, df, p_value = _joint_wald(vector, covariance[np.ix_(indices, indices)])
    rows = []
    for response in responses:
        row = indexed.loc[difference_names[response]]
        rows.append(
            {
                **base,
                "response": response,
                "distance_x_area_difference_b_minus_a": float(row["estimate"]),
                "difference_se": float(row["cluster_robust_se"]),
                "difference_p": float(row["p_value"]),
                "n_islands_a": int(counts.loc[response, context_a]),
                "n_islands_b": int(counts.loc[response, context_b]),
            }
        )
    return pd.DataFrame(rows), {
        **base,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_context_difference_chisq": stat,
        "joint_context_difference_df": df,
        "p_value": p_value,
    }


def _apply_inference(
    coefficients: pd.DataFrame,
    within: pd.DataFrame,
    between_coefficients: pd.DataFrame,
    between: pd.DataFrame,
    *,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    family = [
        "family",
        "context_layer",
        "stratum",
        "support_tier",
        "source_mode",
    ]
    if not within.empty:
        fit = within["status"].eq("fit")
        within["q_vector_family"] = np.nan
        within.loc[fit, "q_vector_family"] = (
            within.loc[fit].groupby(family, group_keys=False)["p_value"].transform(_bh)
        )
        within["area_moderation_vector_supported"] = (
            within["q_vector_family"].le(alpha).fillna(False)
        )
    if not coefficients.empty:
        fit_keys = [*family, "context"]
        coefficients["distance_q"] = coefficients.groupby(
            family, group_keys=False
        )["distance_p"].transform(_bh)
        coefficients["distance_x_area_q"] = coefficients.groupby(
            family, group_keys=False
        )["distance_x_area_p"].transform(_bh)
        gate = within[fit_keys + ["area_moderation_vector_supported"]]
        coefficients = coefficients.merge(gate, on=fit_keys, how="left", validate="many_to_one")
        coefficients["distance_axis_supported"] = coefficients["distance_q"].le(alpha)
        coefficients["interaction_axis_supported"] = coefficients[
            "distance_x_area_q"
        ].le(alpha)
        coefficients["area_moderation_state"] = "no_supported_axis_moderation"
        supported = (
            coefficients["area_moderation_vector_supported"].fillna(False)
            & coefficients["distance_axis_supported"].fillna(False)
            & coefficients["interaction_axis_supported"].fillna(False)
        )
        opposite = (
            coefficients["distance_estimate_at_mean_area"]
            * coefficients["distance_x_area_estimate"]
        ).lt(0)
        coefficients.loc[
            supported & opposite, "area_moderation_state"
        ] = "distance_effect_stronger_on_smaller_islands"
        coefficients.loc[
            supported & ~opposite, "area_moderation_state"
        ] = "distance_effect_stronger_on_larger_islands"
        vector_only = (
            coefficients["area_moderation_vector_supported"].fillna(False) & ~supported
        )
        coefficients.loc[
            vector_only, "area_moderation_state"
        ] = "supported_vector_without_axiswise_amplification"

    if not between.empty:
        fit = between["status"].eq("fit")
        between["q_between_family"] = np.nan
        between.loc[fit, "q_between_family"] = (
            between.loc[fit].groupby(family, group_keys=False)["p_value"].transform(_bh)
        )
        between["area_moderation_difference_supported"] = between[
            "q_between_family"
        ].le(alpha).fillna(False)
    if not between_coefficients.empty:
        between_coefficients["difference_q"] = between_coefficients.groupby(
            family, group_keys=False
        )["difference_p"].transform(_bh)
        keys = [*family, "context_a", "context_b"]
        gate = between[keys + ["area_moderation_difference_supported"]]
        between_coefficients = between_coefficients.merge(
            gate, on=keys, how="left", validate="many_to_one"
        )
        between_coefficients["difference_axis_supported"] = (
            between_coefficients["area_moderation_difference_supported"].fillna(False)
            & between_coefficients["difference_q"].le(alpha)
        )
    return coefficients, within, between_coefficients, between


def run_area_capacity_analysis(
    island_scores: pd.DataFrame,
    lineage_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    branching_config: dict[str, Any],
    area_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    families = build_families(
        island_scores, lineage_scores, branching_config, area_config
    )
    primary_strata = [str(x) for x in area_config["primary_strata"]]
    tiers = {str(k): int(v) for k, v in area_config["support_tiers"].items()}
    coefficient_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_coefficient_parts: list[pd.DataFrame] = []
    between_rows: list[dict[str, Any]] = []

    for layer, layer_spec in area_config["context_layers"].items():
        context_column = str(layer_spec["column"])
        layer_covariates = covariates.copy()
        if context_column not in layer_covariates.columns:
            assignment = realm_assignment[["island_id", context_column]].drop_duplicates(
                "island_id"
            )
            layer_covariates = layer_covariates.merge(
                assignment,
                on="island_id",
                how="left",
                validate="one_to_one",
            )
        contexts = [str(x) for x in layer_spec["contexts"]]
        contrasts = [[str(x) for x in pair] for pair in layer_spec["direct_contrasts"]]
        for family, family_data in families.items():
            prepared = _prepare(
                family_data,
                layer_covariates,
                area_config,
                context_column=context_column,
            )
            for source_mode in prepared["source_mode"].dropna().astype(str).unique():
                source_data = prepared.loc[prepared["source_mode"].eq(source_mode)].copy()
                for stratum, (support_tier, threshold) in itertools.product(
                    primary_strata, tiers.items()
                ):
                    for context in contexts:
                        coefficients, result = _fit_within(
                            source_data,
                            area_config,
                            family=family,
                            source_mode=source_mode,
                            stratum=stratum,
                            context=context,
                            context_column=context_column,
                            context_layer=str(layer),
                            support_tier=support_tier,
                            threshold=threshold,
                        )
                        if not coefficients.empty:
                            coefficient_parts.append(coefficients)
                        within_rows.append(result)
                    for context_a, context_b in contrasts:
                        differences, result = _fit_between(
                            source_data,
                            area_config,
                            family=family,
                            source_mode=source_mode,
                            stratum=stratum,
                            context_a=context_a,
                            context_b=context_b,
                            context_column=context_column,
                            context_layer=str(layer),
                            support_tier=support_tier,
                            threshold=threshold,
                        )
                        if not differences.empty:
                            between_coefficient_parts.append(differences)
                        between_rows.append(result)

    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True)
        if coefficient_parts
        else pd.DataFrame()
    )
    within = pd.DataFrame(within_rows)
    between_coefficients = (
        pd.concat(between_coefficient_parts, ignore_index=True)
        if between_coefficient_parts
        else pd.DataFrame()
    )
    between = pd.DataFrame(between_rows)
    return _apply_inference(
        coefficients,
        within,
        between_coefficients,
        between,
        alpha=float(area_config["multiple_testing"]["alpha"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--island-scores-csv", type=Path, required=True)
    parser.add_argument("--lineage-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--branching-config-path", type=Path, required=True)
    parser.add_argument("--area-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    branching_config = yaml.safe_load(
        args.branching_config_path.read_text(encoding="utf-8")
    )
    area_config = yaml.safe_load(args.area_config_path.read_text(encoding="utf-8"))
    coefficients, within, between_coefficients, between = run_area_capacity_analysis(
        pd.read_csv(args.island_scores_csv),
        pd.read_csv(args.lineage_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        branching_config,
        area_config,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(args.output_dir / "area_moderation_coefficients.csv", index=False)
    within.to_csv(args.output_dir / "area_moderation_within_context.csv", index=False)
    between_coefficients.to_csv(
        args.output_dir / "area_moderation_between_coefficients.csv", index=False
    )
    between.to_csv(args.output_dir / "area_moderation_between_context.csv", index=False)
    manifest = {
        "contract": str(area_config["contract"]),
        "model": str(area_config["model"]),
        "input_sha256": {
            "island_scores": _sha256(args.island_scores_csv),
            "lineage_scores": _sha256(args.lineage_scores_csv),
            "covariates": _sha256(args.covariates_csv),
            "realm_assignment": _sha256(args.realm_assignment_csv),
            "branching_config": _sha256_text(args.branching_config_path),
            "area_config": _sha256_text(args.area_config_path),
        },
        "n_coefficient_rows": len(coefficients),
        "n_within_rows": len(within),
        "n_between_coefficient_rows": len(between_coefficients),
        "n_between_rows": len(between),
        "claim_ceiling": deepcopy(area_config["claim_ceiling"]),
        "continuous_area_no_threshold": True,
        "joint_vector_gate_before_axis_classification": True,
    }
    (args.output_dir / "area_moderation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
