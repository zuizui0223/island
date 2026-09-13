"""Matched context x taxonomic-depth audit for Chapter 1.

This module is a secondary publication-promotion audit. It re-expresses the frozen
V4 floral-architecture decomposition using the same four architecture components at
two nested analysis stages:

1. source-adjusted pre-genus architecture response; and
2. source-matched beyond-genus residual.

For every source mode, evidence scope and floristic stratum, the analysis forms the
paired island-level stage difference and directly tests whether its isolation slope
differs between northern-midlatitude and tropical contexts. The output remains a
plant-architecture result. It does not identify pollinator identity, historical loss,
effective service, in-situ evolution, or a causal mechanism.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_context_analysis import _chi_square_sf_integer_df
from island_v2.chapter1_pr136_biogeographic_residual import (
    _fit_weighted_clustered_design,
)
from island_v2.chapter1_pr138_syndrome_analysis import _bh

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _two_sided_normal_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return float("nan")
    return math.erfc(abs(float(z_value)) / math.sqrt(2.0))


def _joint_wald(
    beta: np.ndarray,
    covariance: np.ndarray,
) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    statistic = float(beta @ np.linalg.pinv(covariance) @ beta)
    return statistic, rank, _chi_square_sf_integer_df(statistic, rank)


def _validate_unique(
    frame: pd.DataFrame,
    keys: list[str],
    *,
    label: str,
) -> None:
    duplicated = frame.duplicated(keys, keep=False)
    if duplicated.any():
        examples = frame.loc[duplicated, keys].head(5).to_dict("records")
        raise ValueError(f"{label} has duplicate keys: {examples}")


def build_paired_stage_delta(
    pre_genus: pd.DataFrame,
    genus_long: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    source_mode: str,
    stratum: str,
    components: list[str],
    contexts: list[str],
    context_column: str,
    cluster_column: str,
    geography_column: str,
    baseline_covariates: list[str],
) -> pd.DataFrame:
    """Build common-support island x component post-minus-pre stage differences."""

    required_pre = {
        "island_id",
        "stratum",
        "source_mode",
        "syndrome",
        "syndrome_score",
    }
    required_genus = {
        "island_id",
        "stratum",
        "source_mode",
        "architecture_component",
        "lineage_outcome",
        "lineage_value",
    }
    required_covariates = {
        "island_id",
        context_column,
        cluster_column,
        geography_column,
        *baseline_covariates,
    }
    if missing := required_pre - set(pre_genus.columns):
        raise ValueError(f"pre-genus table missing columns: {sorted(missing)}")
    if missing := required_genus - set(genus_long.columns):
        raise ValueError(f"genus table missing columns: {sorted(missing)}")
    if missing := required_covariates - set(covariates.columns):
        raise ValueError(f"covariates missing columns: {sorted(missing)}")

    pre = pre_genus.loc[
        pre_genus["source_mode"].astype(str).eq(source_mode)
        & pre_genus["stratum"].astype(str).eq(stratum)
        & pre_genus["syndrome"].astype(str).isin(components),
        ["island_id", "syndrome", "syndrome_score"],
    ].copy()
    pre = pre.rename(
        columns={
            "syndrome": "architecture_component",
            "syndrome_score": "pre_genus_response",
        }
    )
    pre["island_id"] = pre["island_id"].astype(str)
    pre["architecture_component"] = pre["architecture_component"].astype(str)
    pre["pre_genus_response"] = pd.to_numeric(
        pre["pre_genus_response"], errors="coerce"
    )
    pre = pre.dropna(subset=["pre_genus_response"])
    _validate_unique(
        pre,
        ["island_id", "architecture_component"],
        label="pre-genus table",
    )

    post = genus_long.loc[
        genus_long["source_mode"].astype(str).eq(source_mode)
        & genus_long["stratum"].astype(str).eq(stratum)
        & genus_long["architecture_component"].astype(str).isin(components)
        & genus_long["lineage_outcome"].astype(str).eq("beyond_genus_residual"),
        [
            "island_id",
            "architecture_component",
            "lineage_value",
            "n_represented_species",
            "n_represented_genera",
        ],
    ].copy()
    post["island_id"] = post["island_id"].astype(str)
    post["architecture_component"] = post["architecture_component"].astype(str)
    post["beyond_genus_response"] = pd.to_numeric(
        post["lineage_value"], errors="coerce"
    )
    post = post.drop(columns="lineage_value").dropna(
        subset=["beyond_genus_response"]
    )
    _validate_unique(
        post,
        ["island_id", "architecture_component"],
        label="beyond-genus table",
    )

    paired = pre.merge(
        post,
        on=["island_id", "architecture_component"],
        how="inner",
        validate="one_to_one",
    )

    covariate_columns = [
        "island_id",
        context_column,
        cluster_column,
        geography_column,
        *baseline_covariates,
    ]
    cov = covariates[covariate_columns].copy()
    cov["island_id"] = cov["island_id"].astype(str)
    _validate_unique(cov, ["island_id"], label="covariate table")
    paired = paired.merge(
        cov,
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric_columns = [geography_column, *baseline_covariates]
    for column in numeric_columns:
        paired[column] = pd.to_numeric(paired[column], errors="coerce")
    paired[context_column] = paired[context_column].fillna("").astype(str)
    paired[cluster_column] = paired[cluster_column].fillna("").astype(str)
    paired = paired.loc[paired[context_column].isin(contexts)].copy()
    paired = paired.dropna(subset=numeric_columns)
    paired = paired.loc[paired[cluster_column].ne("")].copy()
    paired["stage_delta"] = (
        paired["beyond_genus_response"] - paired["pre_genus_response"]
    )
    paired["source_mode"] = source_mode
    paired["stratum"] = stratum
    return paired


def _standardize_island_predictors(
    paired: pd.DataFrame,
    *,
    geography_column: str,
    baseline_covariates: list[str],
) -> pd.DataFrame:
    """Standardize predictors once on the unique paired-island pool."""

    columns = [geography_column, *baseline_covariates]
    island_values = paired[["island_id", *columns]].drop_duplicates("island_id")
    if island_values["island_id"].duplicated().any():
        raise ValueError("island predictors are not unique after deduplication")
    standardized = island_values[["island_id"]].copy()
    for column in columns:
        values = pd.to_numeric(island_values[column], errors="coerce")
        mean = float(values.mean())
        sd = float(values.std(ddof=0))
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant or invalid predictor: {column}")
        standardized[column] = (values - mean) / sd
    return standardized


def fit_matched_depth_interaction(
    paired: pd.DataFrame,
    *,
    components: list[str],
    contexts: list[str],
    context_column: str,
    cluster_column: str,
    geography_column: str,
    baseline_covariates: list[str],
    minimum_paired_islands: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fit the four-component context x taxonomic-stage interaction.

    Stage is encoded through ``stage_delta = beyond_genus - pre_genus``. The
    context contrast of the distance slope of this delta is algebraically the
    context x taxonomic-stage x distance interaction for each architecture
    component.
    """

    counts = (
        paired.groupby(["architecture_component", context_column])["island_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    base: dict[str, Any] = {
        "status": "not_testable",
        "minimum_paired_islands": int(minimum_paired_islands),
    }
    for component in components:
        for context in contexts:
            value = (
                int(counts.loc[component, context])
                if component in counts.index and context in counts.columns
                else 0
            )
            base[f"n_{context}_{component}"] = value
            if value < minimum_paired_islands:
                base["failure_reason"] = (
                    f"{component} x {context} paired support {value} "
                    f"< {minimum_paired_islands}"
                )
                return pd.DataFrame(), base

    standardized = _standardize_island_predictors(
        paired,
        geography_column=geography_column,
        baseline_covariates=baseline_covariates,
    )
    work = paired.drop(
        columns=[geography_column, *baseline_covariates]
    ).merge(
        standardized,
        on="island_id",
        how="left",
        validate="many_to_one",
    )

    names: list[str] = []
    design_columns: list[np.ndarray] = []
    slope_names: dict[tuple[str, str], str] = {}

    for component in components:
        for context in contexts:
            mask = (
                work["architecture_component"].astype(str).eq(component)
                & work[context_column].astype(str).eq(context)
            ).to_numpy()
            indicator = mask.astype(float)

            intercept_name = f"{component}|{context}|intercept"
            names.append(intercept_name)
            design_columns.append(indicator)

            for predictor in baseline_covariates:
                name = f"{component}|{context}|z_{predictor}"
                names.append(name)
                design_columns.append(
                    indicator * work[predictor].to_numpy(float)
                )

            slope_name = f"{component}|{context}|z_{geography_column}"
            names.append(slope_name)
            design_columns.append(
                indicator * work[geography_column].to_numpy(float)
            )
            slope_names[(component, context)] = slope_name

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["stage_delta"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(design_columns),
        names,
        work[cluster_column].to_numpy(str),
    )
    if coefficients.empty:
        return pd.DataFrame(), {
            **base,
            "status": str(fit.get("status", "fit_failed")),
            "failure_reason": str(fit.get("status", "fit_failed")),
        }

    coefficient_lookup = coefficients.set_index("predictor")
    name_index = {name: index for index, name in enumerate(names)}
    component_rows: list[dict[str, Any]] = []
    contrast_vectors: list[np.ndarray] = []
    contrast_estimates: list[float] = []

    context_a, context_b = contexts
    for component in components:
        a_name = slope_names[(component, context_a)]
        b_name = slope_names[(component, context_b)]
        a_index = name_index[a_name]
        b_index = name_index[b_name]
        a_row = coefficient_lookup.loc[a_name]
        b_row = coefficient_lookup.loc[b_name]

        contrast = np.zeros(len(names), dtype=float)
        contrast[b_index] = 1.0
        contrast[a_index] = -1.0
        estimate = float(b_row["estimate"] - a_row["estimate"])
        variance = float(contrast @ covariance @ contrast)
        stderr = math.sqrt(max(variance, 0.0))
        z_value = estimate / stderr if stderr > 0 else float("nan")
        p_value = _two_sided_normal_p(z_value)

        contrast_vectors.append(contrast)
        contrast_estimates.append(estimate)
        component_rows.append(
            {
                "architecture_component": component,
                f"{context_a}_stage_delta_slope": float(a_row["estimate"]),
                f"{context_a}_stage_delta_se": float(
                    a_row["cluster_robust_se"]
                ),
                f"{context_b}_stage_delta_slope": float(b_row["estimate"]),
                f"{context_b}_stage_delta_se": float(
                    b_row["cluster_robust_se"]
                ),
                "context_difference_in_stage_delta_slope": estimate,
                "cluster_robust_se": stderr,
                "ci95_low": estimate - 1.96 * stderr,
                "ci95_high": estimate + 1.96 * stderr,
                "z_value": z_value,
                "p_value": p_value,
                f"n_islands_{context_a}": int(counts.loc[component, context_a]),
                f"n_islands_{context_b}": int(counts.loc[component, context_b]),
            }
        )

    contrast_matrix = np.vstack(contrast_vectors)
    contrast_covariance = contrast_matrix @ covariance @ contrast_matrix.T
    statistic, df, p_value = _joint_wald(
        np.asarray(contrast_estimates, dtype=float),
        contrast_covariance,
    )

    summary = {
        **base,
        "status": "fit",
        "failure_reason": "",
        "n_rows": int(len(work)),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_context_by_stage_chisq": statistic,
        "joint_context_by_stage_df": df,
        "p_value": p_value,
        "context_a": context_a,
        "context_b": context_b,
    }
    return pd.DataFrame(component_rows), summary


def run_depth_audit(
    *,
    all_pre_genus: pd.DataFrame,
    all_genus_long: pd.DataFrame,
    direct_pre_genus: pd.DataFrame,
    direct_genus_long: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    spec = config["matched_depth_test"]
    components = [str(x) for x in spec["components"]]
    contexts = [str(x) for x in spec["contexts"]]
    strata = [str(x) for x in spec["floristic_strata"]]
    source_modes = [str(x) for x in spec["source_modes"]]
    evidence_scopes = [str(x) for x in spec["evidence_scopes"]]
    context_column = str(spec["context_column"])
    cluster_column = str(spec["cluster_column"])
    geography_column = str(spec["geography_column"])
    baseline_covariates = [str(x) for x in spec["baseline_covariates"]]
    minimum = int(spec["minimum_paired_islands_per_component_context"])
    alpha = float(spec["alpha"])

    tables = {
        "all_analysis_eligible": (all_pre_genus, all_genus_long),
        "direct_only": (direct_pre_genus, direct_genus_long),
    }
    if set(evidence_scopes) != set(tables):
        raise ValueError(
            "evidence scopes in config must be exactly "
            "all_analysis_eligible and direct_only"
        )

    component_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    for evidence_scope in evidence_scopes:
        pre_genus, genus_long = tables[evidence_scope]
        for stratum in strata:
            for source_mode in source_modes:
                paired = build_paired_stage_delta(
                    pre_genus,
                    genus_long,
                    covariates,
                    source_mode=source_mode,
                    stratum=stratum,
                    components=components,
                    contexts=contexts,
                    context_column=context_column,
                    cluster_column=cluster_column,
                    geography_column=geography_column,
                    baseline_covariates=baseline_covariates,
                )
                component, summary = fit_matched_depth_interaction(
                    paired,
                    components=components,
                    contexts=contexts,
                    context_column=context_column,
                    cluster_column=cluster_column,
                    geography_column=geography_column,
                    baseline_covariates=baseline_covariates,
                    minimum_paired_islands=minimum,
                )
                summary.update(
                    {
                        "evidence_scope": evidence_scope,
                        "stratum": stratum,
                        "source_mode": source_mode,
                    }
                )
                summary_rows.append(summary)
                if not component.empty:
                    component.insert(0, "source_mode", source_mode)
                    component.insert(0, "stratum", stratum)
                    component.insert(0, "evidence_scope", evidence_scope)
                    component_parts.append(component)

    components_frame = (
        pd.concat(component_parts, ignore_index=True)
        if component_parts
        else pd.DataFrame()
    )
    summary_frame = pd.DataFrame(summary_rows)
    summary_frame["q_across_source_modes"] = np.nan
    fit_mask = summary_frame["status"].eq("fit")
    if fit_mask.any():
        summary_frame.loc[
            fit_mask, "q_across_source_modes"
        ] = summary_frame.loc[fit_mask].groupby(
            ["evidence_scope", "stratum"], group_keys=False
        )["p_value"].transform(_bh)
    summary_frame["source_mode_supported"] = (
        summary_frame["q_across_source_modes"].le(alpha).fillna(False)
    )

    support = (
        summary_frame.loc[summary_frame["status"].eq("fit")]
        .groupby(["evidence_scope", "stratum"], as_index=False)
        .agg(
            n_source_modes=("source_mode", "nunique"),
            n_supported_source_modes=("source_mode_supported", "sum"),
            max_q=("q_across_source_modes", "max"),
        )
    )
    support_lookup = {
        (str(row.evidence_scope), str(row.stratum)): {
            "n_source_modes": int(row.n_source_modes),
            "n_supported_source_modes": int(row.n_supported_source_modes),
            "max_q": float(row.max_q),
        }
        for row in support.itertuples(index=False)
    }
    expected_modes = len(source_modes)

    def _full(scope: str, current_stratum: str) -> bool:
        row = support_lookup.get((scope, current_stratum), {})
        return (
            row.get("n_source_modes") == expected_modes
            and row.get("n_supported_source_modes") == expected_modes
        )

    nne_full = all(
        _full(scope, "native_nonendemic") for scope in evidence_scopes
    )
    all_native_full = all(
        _full(scope, "all_native") for scope in evidence_scopes
    )
    if nne_full and all_native_full:
        classification = "full_hierarchical_depth_promotion"
    elif nne_full:
        classification = "bounded_hierarchical_depth_signal"
    else:
        classification = "no_promotion"

    support_summary = []
    for scope in evidence_scopes:
        for current_stratum in strata:
            row = support_lookup.get((scope, current_stratum))
            support_summary.append(
                {
                    "evidence_scope": scope,
                    "stratum": current_stratum,
                    "n_source_modes": int(row["n_source_modes"]) if row else 0,
                    "n_supported_source_modes": (
                        int(row["n_supported_source_modes"]) if row else 0
                    ),
                    "max_q": float(row["max_q"]) if row else None,
                }
            )

    manifest = {
        "contract": str(config["contract"]),
        "status": "matched_depth_audit_complete",
        "components": components,
        "contexts": contexts,
        "source_modes": source_modes,
        "evidence_scopes": evidence_scopes,
        "floristic_strata": strata,
        "minimum_paired_islands_per_component_context": minimum,
        "alpha": alpha,
        "support_summary": support_summary,
        "publication_gate_classification": classification,
        "claim_boundary": (
            "This result compares matched plant-architecture responses across "
            "taxonomic stages and biogeographic contexts. It does not identify "
            "pollinator identity, historical loss, effective service, in-situ "
            "evolution, or causal mechanism."
        ),
    }
    return components_frame, summary_frame, manifest


@app.command("run")
def run_command(
    all_pre_genus_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    all_genus_long_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    direct_pre_genus_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    direct_genus_long_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    covariates_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_el_hierarchical_depth_audit.yml"),
        exists=True,
        dir_okay=False,
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    expected_contract = "chapter1_el_hierarchical_depth_audit_v1"
    if str(config.get("contract", "")) != expected_contract:
        raise typer.BadParameter(
            f"expected contract {expected_contract}, got {config.get('contract')}"
        )

    components, summary, manifest = run_depth_audit(
        all_pre_genus=pd.read_csv(all_pre_genus_csv),
        all_genus_long=pd.read_csv(all_genus_long_csv),
        direct_pre_genus=pd.read_csv(direct_pre_genus_csv),
        direct_genus_long=pd.read_csv(direct_genus_long_csv),
        covariates=pd.read_csv(covariates_csv),
        config=config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    components.to_csv(
        output_dir / "matched_depth_component_effects.csv",
        index=False,
    )
    summary.to_csv(
        output_dir / "matched_depth_joint_tests.csv",
        index=False,
    )
    (output_dir / "matched_depth_manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2, allow_nan=False))


if __name__ == "__main__":
    app()
