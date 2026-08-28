"""V3 falsification of area moderation by response-measurement support.

This module treats species/lineage support only as a precision sensitivity.  It
does not condition the primary biological model on richness and does not turn
area into a direct measurement of founder effects or pollination service.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_area_capacity_moderation import _fit_within, _sha256_text
from island_v2.chapter1_global_branching import build_branch_scores
from island_v2.chapter1_pr138_syndrome_analysis import _bh


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _plant_family(
    island_scores: pd.DataFrame,
    branching_config: dict[str, Any],
) -> pd.DataFrame:
    required = {
        "island_id",
        "stratum",
        "syndrome",
        "n_species",
        "mean_informative_traits",
    }
    missing = required - set(island_scores.columns)
    if missing:
        raise ValueError(f"plant scores missing support columns: {sorted(missing)}")
    branches = build_branch_scores(island_scores, branching_config)
    rows: list[pd.DataFrame] = []
    for axis in ["accessibility_generalization", "reproductive_assurance"]:
        spec = branching_config["branch_axes"][axis]
        components = [str(x) for x in spec["components"]]
        component_support = island_scores.loc[
            island_scores["syndrome"].isin(components),
            [
                "island_id",
                "stratum",
                "syndrome",
                "n_species",
                "mean_informative_traits",
            ],
        ].copy()
        support = component_support.groupby(["island_id", "stratum"], as_index=False).agg(
            n_response_species=("n_species", "min"),
            mean_informative_traits=("mean_informative_traits", "min"),
            n_components=("syndrome", "nunique"),
        )
        support = support.loc[support["n_components"].eq(len(components))]
        selected = branches.loc[
            branches["syndrome"].eq(axis),
            ["island_id", "stratum", "syndrome", "syndrome_score"],
        ].rename(columns={"syndrome": "response", "syndrome_score": "response_score"})
        selected = selected.merge(
            support.drop(columns="n_components"),
            on=["island_id", "stratum"],
            how="inner",
            validate="one_to_one",
        )
        rows.append(selected)
    out = pd.concat(rows, ignore_index=True)
    out["family"] = "plant_side_branches"
    out["source_mode"] = "not_applicable"
    out["n_response_groups"] = np.nan
    return out


def _lineage_family(lineage_scores: pd.DataFrame, evidence_scope: str) -> pd.DataFrame:
    responses = ["entry_enrichment", "loading_increment"]
    required = {
        "island_id",
        "stratum",
        "source_mode",
        "evidence_scope",
        "source_matching",
        "minimum_represented_genera",
        "n_represented_species",
        "n_represented_genera",
        *responses,
    }
    missing = required - set(lineage_scores.columns)
    if missing:
        raise ValueError(f"lineage scores missing support columns: {sorted(missing)}")
    selected = lineage_scores.loc[
        lineage_scores["evidence_scope"].eq(evidence_scope)
        & lineage_scores["source_matching"].eq("prevalence_richness")
        & lineage_scores["minimum_represented_genera"].eq(5)
    ].copy()
    out = selected.melt(
        id_vars=[
            "island_id",
            "stratum",
            "source_mode",
            "n_represented_species",
            "n_represented_genera",
        ],
        value_vars=responses,
        var_name="response",
        value_name="response_score",
    ).rename(
        columns={
            "n_represented_species": "n_response_species",
            "n_represented_genera": "n_response_groups",
        }
    )
    out["family"] = "source_lineage_broad"
    out["mean_informative_traits"] = np.nan
    return out


def build_support_families(
    all_island_scores: pd.DataFrame,
    direct_island_scores: pd.DataFrame,
    lineage_scores: pd.DataFrame,
    branching_config: dict[str, Any],
) -> pd.DataFrame:
    parts = []
    for scope, scores, lineage_scope in [
        ("all_analysis_eligible", all_island_scores, "broad"),
        ("direct_only", direct_island_scores, "broad_direct"),
    ]:
        plant = _plant_family(scores, branching_config)
        plant["evidence_scope"] = scope
        lineage = _lineage_family(lineage_scores, lineage_scope)
        lineage["evidence_scope"] = scope
        parts.extend([plant, lineage])
    out = pd.concat(parts, ignore_index=True)
    out["n_response_species"] = pd.to_numeric(out["n_response_species"], errors="coerce")
    return out


def add_capped_information_weights(
    data: pd.DataFrame,
    *,
    cap_species: int,
) -> pd.DataFrame:
    out = data.copy()
    raw = out["n_response_species"].clip(lower=1, upper=cap_species)
    means = raw.groupby(out["response"]).transform("mean")
    out["analysis_weight"] = raw / means
    return out


def common_support_area_overlap(
    data: pd.DataFrame,
    *,
    area_column: str,
    minimum_per_half: int,
    lower_quantile: float = 0.05,
    upper_quantile: float = 0.95,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    islands = (
        data[["island_id", "n_response_species", area_column]].drop_duplicates("island_id").dropna()
    )
    if islands.empty:
        return data.iloc[0:0].copy(), {"status": "not_testable_no_support"}
    median = float(islands["n_response_species"].median())
    low = islands.loc[islands["n_response_species"].le(median)]
    high = islands.loc[islands["n_response_species"].gt(median)]
    base = {
        "median_support": median,
        "n_low_before": len(low),
        "n_high_before": len(high),
    }
    if len(low) < minimum_per_half or len(high) < minimum_per_half:
        return data.iloc[0:0].copy(), {**base, "status": "not_testable_support_halves"}
    low_range = (
        float(low[area_column].quantile(lower_quantile)),
        float(low[area_column].quantile(upper_quantile)),
    )
    high_range = (
        float(high[area_column].quantile(lower_quantile)),
        float(high[area_column].quantile(upper_quantile)),
    )
    lower = max(low_range[0], high_range[0])
    upper = min(low_range[1], high_range[1])
    if not lower < upper:
        return data.iloc[0:0].copy(), {
            **base,
            "status": "not_testable_no_area_overlap",
            "retained_area_lower": lower,
            "retained_area_upper": upper,
        }
    retained_ids = set(
        islands.loc[islands[area_column].between(lower, upper), "island_id"].astype(str)
    )
    retained = data.loc[data["island_id"].astype(str).isin(retained_ids)].copy()
    retained_islands = retained[["island_id", "n_response_species"]].drop_duplicates("island_id")
    return retained, {
        **base,
        "status": "retained",
        "retained_area_lower": lower,
        "retained_area_upper": upper,
        "n_retained": int(retained_islands["island_id"].nunique()),
        "n_low_retained": int(retained_islands["n_response_species"].le(median).sum()),
        "n_high_retained": int(retained_islands["n_response_species"].gt(median).sum()),
    }


def _context_covariates(
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    context_column: str,
) -> pd.DataFrame:
    if context_column in covariates.columns:
        return covariates.copy()
    assignment = realm_assignment[["island_id", context_column]].drop_duplicates("island_id")
    return covariates.merge(assignment, on="island_id", how="left", validate="one_to_one")


def _prepare_family(
    family: pd.DataFrame,
    covariates: pd.DataFrame,
    area_config: dict[str, Any],
    context_column: str,
) -> pd.DataFrame:
    columns = [
        "island_id",
        area_config["distance_column"],
        area_config["area_column"],
        area_config["cluster_column"],
        context_column,
        *area_config["control_columns"],
    ]
    data = family.merge(
        covariates[columns].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = [
        "response_score",
        "n_response_species",
        area_config["distance_column"],
        area_config["area_column"],
        *area_config["control_columns"],
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data


def _apply_sensitivity_inference(
    coefficients: pd.DataFrame,
    omnibus: pd.DataFrame,
    *,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if coefficients.empty or omnibus.empty:
        return coefficients, omnibus
    family = [
        "family",
        "evidence_scope",
        "context_layer",
        "stratum",
        "support_tier",
        "source_mode",
    ]
    fit = omnibus["status"].eq("fit")
    omnibus["q_vector_sensitivity_family"] = np.nan
    omnibus.loc[fit, "q_vector_sensitivity_family"] = (
        omnibus.loc[fit].groupby(family, group_keys=False)["p_value"].transform(_bh)
    )
    omnibus["area_moderation_vector_supported"] = (
        omnibus["q_vector_sensitivity_family"].le(alpha).fillna(False)
    )
    coefficients["distance_q_sensitivity_family"] = coefficients.groupby(family, group_keys=False)[
        "distance_p"
    ].transform(_bh)
    coefficients["interaction_q_sensitivity_family"] = coefficients.groupby(
        family, group_keys=False
    )["distance_x_area_p"].transform(_bh)
    gate_keys = [*family, "sensitivity_mode", "context"]
    coefficients = coefficients.merge(
        omnibus[gate_keys + ["area_moderation_vector_supported"]],
        on=gate_keys,
        how="left",
        validate="many_to_one",
    )
    coefficients["distance_axis_supported"] = coefficients["distance_q_sensitivity_family"].le(
        alpha
    )
    coefficients["interaction_axis_supported"] = coefficients[
        "interaction_q_sensitivity_family"
    ].le(alpha)
    coefficients["small_island_amplification_supported"] = (
        coefficients["area_moderation_vector_supported"].fillna(False)
        & coefficients["distance_axis_supported"].fillna(False)
        & coefficients["interaction_axis_supported"].fillna(False)
        & (
            coefficients["distance_estimate_at_mean_area"]
            * coefficients["distance_x_area_estimate"]
        ).lt(0)
    )
    return coefficients, omnibus


def run_support_sensitivities(
    families: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    area_config: dict[str, Any],
    v3_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    coefficient_parts: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    overlap_rows: list[dict[str, Any]] = []
    cap = int(v3_config["weighting_modes"]["capped_information"]["cap_species"])
    minimum_half = int(
        v3_config["common_support"]["minimum_islands_per_support_half_before_trimming"]
    )
    tiers = {str(k): int(v) for k, v in area_config["support_tiers"].items()}
    for layer, layer_spec in area_config["context_layers"].items():
        context_column = str(layer_spec["column"])
        layer_covariates = _context_covariates(covariates, realm_assignment, context_column)
        for (family_name, evidence_scope), family in families.groupby(
            ["family", "evidence_scope"], sort=True
        ):
            prepared = _prepare_family(family, layer_covariates, area_config, context_column)
            for source_mode in sorted(prepared["source_mode"].astype(str).unique()):
                source_data = prepared.loc[
                    prepared["source_mode"].astype(str).eq(source_mode)
                ].copy()
                for stratum in area_config["primary_strata"]:
                    for context in layer_spec["contexts"]:
                        cell = source_data.loc[
                            source_data["stratum"].eq(stratum)
                            & source_data[context_column].eq(context)
                        ].copy()
                        for sensitivity_mode in [
                            "equal_island",
                            "capped_information",
                            "common_support",
                        ]:
                            fit_data = cell.copy()
                            weight_column = None
                            overlap = {"status": "not_applicable"}
                            if sensitivity_mode == "capped_information":
                                fit_data = add_capped_information_weights(fit_data, cap_species=cap)
                                weight_column = "analysis_weight"
                            elif sensitivity_mode == "common_support":
                                fit_data, overlap = common_support_area_overlap(
                                    fit_data,
                                    area_column=str(area_config["area_column"]),
                                    minimum_per_half=minimum_half,
                                )
                            overlap_rows.append(
                                {
                                    "family": family_name,
                                    "evidence_scope": evidence_scope,
                                    "source_mode": source_mode,
                                    "context_layer": layer,
                                    "stratum": stratum,
                                    "context": context,
                                    "sensitivity_mode": sensitivity_mode,
                                    **overlap,
                                }
                            )
                            for support_tier, threshold in tiers.items():
                                coefficients, result = _fit_within(
                                    fit_data,
                                    area_config,
                                    family=str(family_name),
                                    source_mode=source_mode,
                                    stratum=str(stratum),
                                    context=str(context),
                                    context_column=context_column,
                                    context_layer=str(layer),
                                    support_tier=support_tier,
                                    threshold=threshold,
                                    weight_column=weight_column,
                                )
                                result.update(
                                    evidence_scope=evidence_scope,
                                    sensitivity_mode=sensitivity_mode,
                                )
                                omnibus_rows.append(result)
                                if not coefficients.empty:
                                    coefficients["evidence_scope"] = evidence_scope
                                    coefficients["sensitivity_mode"] = sensitivity_mode
                                    coefficient_parts.append(coefficients)
    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True) if coefficient_parts else pd.DataFrame()
    )
    omnibus = pd.DataFrame(omnibus_rows)
    coefficients, omnibus = _apply_sensitivity_inference(
        coefficients,
        omnibus,
        alpha=float(area_config["multiple_testing"]["alpha"]),
    )
    return coefficients, omnibus, pd.DataFrame(overlap_rows)


def support_overlap_diagnostics(
    families: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    area_config: dict[str, Any],
) -> pd.DataFrame:
    rows = []
    area = str(area_config["area_column"])
    for layer, layer_spec in area_config["context_layers"].items():
        context_column = str(layer_spec["column"])
        layer_covariates = _context_covariates(covariates, realm_assignment, context_column)
        data = families.merge(
            layer_covariates[["island_id", context_column, area]].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="many_to_one",
        )
        group_columns = [
            "family",
            "evidence_scope",
            "source_mode",
            "stratum",
            "response",
            context_column,
        ]
        for keys, group in data.groupby(group_columns, dropna=False):
            work = (
                group[["island_id", area, "n_response_species"]]
                .dropna()
                .drop_duplicates("island_id")
            )
            rho = np.nan
            if (
                len(work) >= 3
                and work[area].nunique() > 1
                and work["n_response_species"].nunique() > 1
            ):
                rho = float(work[area].rank().corr(work["n_response_species"].rank()))
            rows.append(
                {
                    "context_layer": layer,
                    **dict(zip(group_columns, keys, strict=True)),
                    "n_islands": len(work),
                    "spearman_area_support": rho,
                    "median_response_species": float(work["n_response_species"].median())
                    if len(work)
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def heteroskedastic_null_simulation(
    families: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    area_config: dict[str, Any],
    v3_config: dict[str, Any],
) -> pd.DataFrame:
    draws = int(v3_config["heteroskedastic_null_simulation"]["draws"])
    seed = int(v3_config["heteroskedastic_null_simulation"]["seed"])
    rng = np.random.default_rng(seed)
    context_column = str(area_config["context_layers"]["biogeographic_realm"]["column"])
    cov = _context_covariates(covariates, realm_assignment, context_column)
    lineage = families.loc[
        families["family"].eq("source_lineage_broad") & families["response"].eq("entry_enrichment")
    ].copy()
    data = _prepare_family(lineage, cov, area_config, context_column)
    distance = str(area_config["distance_column"])
    area = str(area_config["area_column"])
    controls = [str(x) for x in area_config["control_columns"]]
    rows = []
    for keys, group in data.loc[data[context_column].eq("Palearctic")].groupby(
        ["evidence_scope", "source_mode", "stratum"], sort=True
    ):
        needed = ["response_score", "n_response_species", distance, area, *controls]
        work = group.dropna(subset=needed).copy()
        base = dict(zip(["evidence_scope", "source_mode", "stratum"], keys, strict=True))
        if len(work) < 30 or work["n_response_species"].le(0).any():
            rows.append({**base, "status": "not_testable", "n_islands": len(work)})
            continue
        standardized = []
        invalid = False
        for column in [distance, area, *controls]:
            values = work[column].to_numpy(float)
            sd = float(np.std(values, ddof=0))
            if not np.isfinite(sd) or sd <= 0:
                invalid = True
                break
            standardized.append((values - float(np.mean(values))) / sd)
        if invalid:
            rows.append(
                {**base, "status": "not_testable_constant_predictor", "n_islands": len(work)}
            )
            continue
        z_distance, z_area, *z_controls = standardized
        null_design = np.column_stack([np.ones(len(work)), z_distance, z_area, *z_controls])
        full_design = np.column_stack([null_design, z_distance * z_area])
        y = work["response_score"].to_numpy(float)
        null_beta = np.linalg.lstsq(null_design, y, rcond=None)[0]
        fitted = null_design @ null_beta
        rmse = float(np.sqrt(np.mean((y - fitted) ** 2)))
        observed = float(np.linalg.lstsq(full_design, y, rcond=None)[0][-1])
        support = work["n_response_species"].to_numpy(float)
        sigma = rmse * np.sqrt(float(np.mean(support)) / support)
        errors = rng.normal(size=(draws, len(work))) * sigma
        simulated_y = fitted[None, :] + errors
        projection = np.linalg.pinv(full_design)
        simulated = (projection @ simulated_y.T)[-1]
        exceedance = float((1 + np.sum(np.abs(simulated) >= abs(observed))) / (draws + 1))
        rows.append(
            {
                **base,
                "status": "simulated",
                "n_islands": len(work),
                "draws": draws,
                "observed_distance_x_area": observed,
                "null_rmse": rmse,
                "null_exceedance_probability": exceedance,
                "passes_false_positive_threshold": exceedance
                <= float(v3_config["inference"]["simulation_false_positive_threshold"]),
            }
        )
    return pd.DataFrame(rows)


def classify_primary(
    coefficients: pd.DataFrame,
    simulations: pd.DataFrame,
    area_config: dict[str, Any],
) -> pd.DataFrame:
    target = coefficients.loc[
        coefficients["family"].eq("source_lineage_broad")
        & coefficients["context_layer"].eq("biogeographic_realm")
        & coefficients["context"].eq("Palearctic")
        & coefficients["response"].eq("entry_enrichment")
    ].copy()
    rows = []
    keys = ["stratum", "support_tier", "source_mode"]
    all_keys = target.loc[
        target["evidence_scope"].eq("all_analysis_eligible"), keys
    ].drop_duplicates()
    for key_values in all_keys.itertuples(index=False, name=None):
        selector = pd.Series(True, index=target.index)
        for column, value in zip(keys, key_values, strict=True):
            selector &= target[column].eq(value)
        all_scope = target.loc[selector & target["evidence_scope"].eq("all_analysis_eligible")]
        direct = target.loc[selector & target["evidence_scope"].eq("direct_only")]
        expected_modes = {"equal_island", "capped_information", "common_support"}
        modes_complete = set(all_scope["sensitivity_mode"]) == expected_modes
        all_supported = bool(
            modes_complete and all_scope["small_island_amplification_supported"].fillna(False).all()
        )
        direct_testable = not direct.empty
        direct_contradiction = bool(
            direct_testable
            and (direct["distance_estimate_at_mean_area"] * direct["distance_x_area_estimate"])
            .ge(0)
            .any()
        )
        simulation = simulations.loc[
            simulations["evidence_scope"].eq("all_analysis_eligible")
            & simulations["source_mode"].eq(key_values[2])
            & simulations["stratum"].eq(key_values[0])
        ]
        simulation_passes = bool(
            len(simulation) == 1
            and simulation["passes_false_positive_threshold"].fillna(False).iloc[0]
        )
        robust = all_supported and not direct_contradiction and simulation_passes
        rows.append(
            {
                **dict(zip(keys, key_values, strict=True)),
                "all_three_sensitivity_modes_testable": modes_complete,
                "all_three_modes_formally_support_small_island_amplification": all_supported,
                "direct_only_testable": direct_testable,
                "direct_only_contradicts": direct_contradiction,
                "heteroskedastic_null_passes": simulation_passes,
                "classification": (
                    "area_capacity_interpretation_strengthened"
                    if robust
                    else "retain_area_as_measurement_sensitive_modifier_only"
                ),
            }
        )
    return pd.DataFrame(rows)


def run(
    all_island_scores: pd.DataFrame,
    direct_island_scores: pd.DataFrame,
    lineage_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    branching_config: dict[str, Any],
    area_config: dict[str, Any],
    explanation_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    v3 = explanation_config["validations"]["V3_H4_area_support_artifact"]["frozen_implementation"]
    families = build_support_families(
        all_island_scores, direct_island_scores, lineage_scores, branching_config
    )
    coefficients, omnibus, overlap = run_support_sensitivities(
        families, covariates, realm_assignment, area_config, v3
    )
    diagnostics = support_overlap_diagnostics(families, covariates, realm_assignment, area_config)
    simulations = heteroskedastic_null_simulation(
        families, covariates, realm_assignment, area_config, v3
    )
    classification = classify_primary(coefficients, simulations, area_config)
    return {
        "support_families": families,
        "coefficients": coefficients,
        "omnibus": omnibus,
        "overlap": overlap,
        "support_diagnostics": diagnostics,
        "simulations": simulations,
        "classification": classification,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-island-scores-csv", type=Path, required=True)
    parser.add_argument("--direct-island-scores-csv", type=Path, required=True)
    parser.add_argument("--lineage-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--branching-config-path", type=Path, required=True)
    parser.add_argument("--area-config-path", type=Path, required=True)
    parser.add_argument("--explanation-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    branching = yaml.safe_load(args.branching_config_path.read_text(encoding="utf-8"))
    area = yaml.safe_load(args.area_config_path.read_text(encoding="utf-8"))
    explanation = yaml.safe_load(args.explanation_config_path.read_text(encoding="utf-8"))
    outputs = run(
        pd.read_csv(args.all_island_scores_csv),
        pd.read_csv(args.direct_island_scores_csv),
        pd.read_csv(args.lineage_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        branching,
        area,
        explanation,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    names = {
        "support_families": "V3_support_families.csv.gz",
        "coefficients": "V3_area_support_coefficients.csv",
        "omnibus": "V3_area_support_omnibus.csv",
        "overlap": "V3_common_support_overlap.csv",
        "support_diagnostics": "V3_support_overlap_diagnostics.csv",
        "simulations": "V3_heteroskedastic_null_simulation.csv",
        "classification": "V3_primary_classification.csv",
    }
    for key, filename in names.items():
        outputs[key].to_csv(
            args.output_dir / filename,
            index=False,
            compression="gzip" if filename.endswith(".gz") else None,
        )
    manifest = {
        "contract": "chapter1_V3_area_support_artifact_v1",
        "input_sha256": {
            "all_island_scores": _sha256(args.all_island_scores_csv),
            "direct_island_scores": _sha256(args.direct_island_scores_csv),
            "lineage_scores": _sha256(args.lineage_scores_csv),
            "covariates": _sha256(args.covariates_csv),
            "realm_assignment": _sha256(args.realm_assignment_csv),
            "branching_config": _sha256_text(args.branching_config_path),
            "area_config": _sha256_text(args.area_config_path),
            "explanation_config": _sha256_text(args.explanation_config_path),
        },
        "rows": {key: len(value) for key, value in outputs.items()},
        "classification_counts": outputs["classification"]["classification"]
        .value_counts()
        .to_dict()
        if not outputs["classification"].empty
        else {},
        "support_is_not_a_causal_control": True,
        "no_outcome_defined_overlap": True,
        "claim_ceiling": deepcopy(area["claim_ceiling"]),
    }
    (args.output_dir / "chapter1_v3_area_support_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
