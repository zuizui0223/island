"""P2 post-baseline audit of component non-concordance in the floral-island response.

The prospective contract is in ``config/chapter1_p2_component_nonconcordance.yml``.
This module never changes the frozen H2 trait definitions.  It asks whether the same-layer
North--Tropical response-vector difference survives common island support and, as a stronger
sensitivity analysis, a common co-observed species denominator across the two primary axes.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.chapter1_pr138_syndrome_analysis import (
    _between_contexts,
    _joint_wald,
    _prepare,
    _standardize,
    _within_context,
)
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)


class P2Error(ValueError):
    """Raised when a frozen P2 invariant fails."""


def _load_yaml(path: Path) -> dict[str, Any]:
    obj = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise P2Error(f"expected YAML mapping: {path}")
    return obj


def _scope_dir(scope: str) -> str:
    mapping = {"all_analysis_eligible": "all", "direct_only": "direct"}
    if scope not in mapping:
        raise P2Error(f"unknown evidence scope: {scope}")
    return mapping[scope]


def _angle_degrees(a: np.ndarray, b: np.ndarray) -> float:
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na <= 0 or nb <= 0 or not np.isfinite([na, nb]).all():
        return float("nan")
    cosine = float(np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0))
    return float(np.degrees(np.arccos(cosine)))


def _determinant(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def _common_island_scores(
    branch_scores: pd.DataFrame,
    *,
    axes: list[str],
    stratum: str,
) -> pd.DataFrame:
    work = branch_scores.loc[
        branch_scores["stratum"].astype(str).eq(stratum)
        & branch_scores["syndrome"].astype(str).isin(axes)
    ].copy()
    scores = work.pivot(index="island_id", columns="syndrome", values="syndrome_score")
    support = work.pivot(index="island_id", columns="syndrome", values="n_species")
    missing_axes = set(axes) - set(scores.columns)
    if missing_axes:
        return pd.DataFrame(columns=branch_scores.columns)
    common = scores[axes].notna().all(axis=1) & support[axes].notna().all(axis=1)
    island_ids = scores.index[common]
    out = work.loc[work["island_id"].astype(str).isin(island_ids.astype(str))].copy()
    counts = out.groupby("island_id")["syndrome"].nunique()
    valid = counts.index[counts.eq(len(axes))]
    return out.loc[out["island_id"].isin(valid)].copy()


def _build_coobserved_species_scores(
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    *,
    source_axes: dict[str, str],
    strata: list[str],
    minimum_species: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    source_names = list(source_axes.values())
    subset = species_scores.loc[
        species_scores["syndrome"].astype(str).isin(source_names),
        ["accepted_species", "syndrome", "syndrome_concordance"],
    ].copy()
    wide = subset.pivot_table(
        index="accepted_species", columns="syndrome", values="syndrome_concordance", aggfunc="first"
    )
    if not set(source_names).issubset(wide.columns):
        raise P2Error("species-score table lacks one or more primary source axes")
    common = wide[source_names].dropna().copy()
    reverse = {source: target for target, source in source_axes.items()}
    common = common.rename(columns=reverse).reset_index()

    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    merged = flora.merge(common, on="accepted_species", how="inner", validate="many_to_one")
    output: list[pd.DataFrame] = []
    diagnostics: dict[str, int] = {"n_coobserved_species": int(len(common))}
    target_axes = list(source_axes)
    for stratum in strata:
        part = merged.loc[stratum_mask(merged, stratum)].copy()
        part = part.drop_duplicates(["island_id", "accepted_species"])
        if part.empty:
            diagnostics[f"{stratum}_islands"] = 0
            continue
        grouped = part.groupby("island_id", as_index=False).agg(
            n_species=("accepted_species", "nunique"),
            **{axis: (axis, "mean") for axis in target_axes},
        )
        grouped = grouped.loc[grouped["n_species"].ge(int(minimum_species))].copy()
        diagnostics[f"{stratum}_islands"] = int(len(grouped))
        for axis in target_axes:
            out = grouped[["island_id", "n_species", axis]].rename(columns={axis: "syndrome_score"})
            out["syndrome"] = axis
            out["stratum"] = stratum
            output.append(out)
    if not output:
        return pd.DataFrame(columns=["island_id", "n_species", "syndrome_score", "syndrome", "stratum"]), diagnostics
    return pd.concat(output, ignore_index=True), diagnostics


def _between_design(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    axes: list[str],
    pattern_config: dict[str, Any],
) -> tuple[pd.DataFrame, np.ndarray, list[str], list[str], list[str]]:
    geography = str(pattern_config["geography_column"])
    context = str(pattern_config["context_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    work = data.loc[
        data["stratum"].astype(str).eq(stratum)
        & data[context].astype(str).isin([context_a, context_b])
        & data["syndrome"].astype(str).isin(axes)
    ].copy()
    if work.empty:
        raise P2Error("empty between-context design")
    b_indicator = work[context].eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    main_slope_names: list[str] = []
    interaction_names: list[str] = []
    for axis in axes:
        mask = work["syndrome"].astype(str).eq(axis).to_numpy()
        if not mask.any():
            raise P2Error(f"axis absent from between design: {axis}")
        indicator = mask.astype(float)
        names.append(f"syndrome[{axis}]")
        columns.append(indicator)
        for predictor in baseline:
            z = np.zeros(len(work))
            z[mask] = _standardize(work.loc[mask, predictor])
            names.append(f"syndrome[{axis}]:z_{predictor}")
            columns.append(z)
        names.append(f"syndrome[{axis}]:context[{context_b}]")
        columns.append(indicator * b_indicator)
        z_geo = np.zeros(len(work))
        z_geo[mask] = _standardize(work.loc[mask, geography])
        main = f"syndrome[{axis}]:z_{geography}"
        names.append(main)
        columns.append(z_geo)
        main_slope_names.append(main)
        interaction = f"syndrome[{axis}]:z_{geography}:context[{context_b}]"
        names.append(interaction)
        columns.append(z_geo * b_indicator)
        interaction_names.append(interaction)
    return work, np.column_stack(columns), names, main_slope_names, interaction_names


def _fit_between_details(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    axes: list[str],
    threshold: int,
    pattern_config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame]:
    context = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    counts = (
        data.loc[
            data["stratum"].astype(str).eq(stratum)
            & data[context].astype(str).isin([context_a, context_b])
            & data["syndrome"].astype(str).isin(axes)
        ]
        .groupby(["syndrome", context])["island_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    for ctx in [context_a, context_b]:
        if ctx not in counts.columns:
            counts[ctx] = 0
    if any(int(counts.loc[axis, ctx]) < int(threshold) for axis in axes for ctx in [context_a, context_b]):
        return {"status": "not_testable"}, pd.DataFrame()

    work, X, names, main_names, interaction_names = _between_design(
        data,
        stratum=stratum,
        context_a=context_a,
        context_b=context_b,
        axes=axes,
        pattern_config=pattern_config,
    )
    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["syndrome_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        X,
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return {"status": str(fit.get("status", "fit_failed"))}, pd.DataFrame()
    beta = coefficients.set_index("predictor")["estimate"]
    main_idx = [names.index(name) for name in main_names]
    int_idx = [names.index(name) for name in interaction_names]
    main = np.asarray([float(beta[name]) for name in main_names])
    delta = np.asarray([float(beta[name]) for name in interaction_names])
    vector_b = main + delta
    delta_cov = covariance[np.ix_(int_idx, int_idx)]
    stat, df, p_value = _joint_wald(delta, delta_cov)

    T = np.zeros((4, len(names)), dtype=float)
    for j, (mi, di) in enumerate(zip(main_idx, int_idx, strict=True)):
        T[j, mi] = 1.0
        T[2 + j, mi] = 1.0
        T[2 + j, di] = 1.0
    slope_cov = T @ covariance @ T.T
    labels = [
        f"{context_a}:{axes[0]}",
        f"{context_a}:{axes[1]}",
        f"{context_b}:{axes[0]}",
        f"{context_b}:{axes[1]}",
    ]
    cov_rows = []
    for i, row_label in enumerate(labels):
        for j, col_label in enumerate(labels):
            cov_rows.append({"row_parameter": row_label, "column_parameter": col_label, "covariance": float(slope_cov[i, j])})
    return {
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "context_difference_chisq": float(stat),
        "context_difference_df": int(df),
        "p_value": float(p_value),
        f"{context_a}_{axes[0]}_slope_joint": float(main[0]),
        f"{context_a}_{axes[1]}_slope_joint": float(main[1]),
        f"{context_b}_{axes[0]}_slope_joint": float(vector_b[0]),
        f"{context_b}_{axes[1]}_slope_joint": float(vector_b[1]),
        f"delta_{axes[0]}": float(delta[0]),
        f"delta_{axes[1]}": float(delta[1]),
    }, pd.DataFrame(cov_rows)


def _light_within_vector(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_value: str,
    axes: list[str],
    pattern_config: dict[str, Any],
) -> np.ndarray:
    geography = str(pattern_config["geography_column"])
    context = str(pattern_config["context_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    slopes: list[float] = []
    for axis in axes:
        part = data.loc[
            data["stratum"].astype(str).eq(stratum)
            & data[context].astype(str).eq(context_value)
            & data["syndrome"].astype(str).eq(axis)
        ].copy()
        if len(part) < len(baseline) + 3:
            raise P2Error("insufficient rows for within-context bootstrap fit")
        columns = [np.ones(len(part), dtype=float)]
        for predictor in baseline:
            columns.append(_standardize(part[predictor]))
        columns.append(_standardize(part[geography]))
        X = np.column_stack(columns)
        y = part["syndrome_score"].to_numpy(float)
        estimate, *_ = np.linalg.lstsq(X, y, rcond=None)
        slopes.append(float(estimate[-1]))
    return np.asarray(slopes, dtype=float)


def _light_between_delta(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    axes: list[str],
    pattern_config: dict[str, Any],
) -> np.ndarray:
    work, X, names, _, interaction_names = _between_design(
        data,
        stratum=stratum,
        context_a=context_a,
        context_b=context_b,
        axes=axes,
        pattern_config=pattern_config,
    )
    estimate, *_ = np.linalg.lstsq(X, work["syndrome_score"].to_numpy(float), rcond=None)
    lookup = {name: float(value) for name, value in zip(names, estimate, strict=True)}
    return np.asarray([lookup[name] for name in interaction_names], dtype=float)


def _bootstrap_profile(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    axes: list[str],
    pattern_config: dict[str, Any],
    draws: int,
    seed: int,
) -> pd.DataFrame:
    context = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    base = data.loc[
        data["stratum"].astype(str).eq(stratum)
        & data[context].astype(str).isin([context_a, context_b])
        & data["syndrome"].astype(str).isin(axes)
    ].copy()
    island_blocks = base[["island_id", cluster]].drop_duplicates()
    if island_blocks["island_id"].duplicated().any():
        raise P2Error("an island maps to multiple spatial blocks")
    blocks = np.asarray(sorted(island_blocks[cluster].astype(str).unique()))
    if len(blocks) < 2:
        raise P2Error("too few spatial blocks for bootstrap")
    rows_by_block = {block: base.loc[base[cluster].astype(str).eq(block)].copy() for block in blocks}
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, Any]] = []
    for draw in range(int(draws)):
        sampled = rng.choice(blocks, size=len(blocks), replace=True)
        pieces: list[pd.DataFrame] = []
        for instance, block in enumerate(sampled):
            piece = rows_by_block[str(block)].copy()
            piece[cluster] = f"bootstrap_{draw:04d}_{instance:04d}_{block}"
            piece["island_id"] = piece["island_id"].astype(str) + f"__b{instance:04d}"
            pieces.append(piece)
        boot = pd.concat(pieces, ignore_index=True)
        try:
            a = _light_within_vector(
                boot, stratum=stratum, context_value=context_a, axes=axes, pattern_config=pattern_config
            )
            b = _light_within_vector(
                boot, stratum=stratum, context_value=context_b, axes=axes, pattern_config=pattern_config
            )
            delta = _light_between_delta(
                boot,
                stratum=stratum,
                context_a=context_a,
                context_b=context_b,
                axes=axes,
                pattern_config=pattern_config,
            )
            rows.append(
                {
                    "draw": draw,
                    "valid": True,
                    f"{context_a}_{axes[0]}": float(a[0]),
                    f"{context_a}_{axes[1]}": float(a[1]),
                    f"{context_b}_{axes[0]}": float(b[0]),
                    f"{context_b}_{axes[1]}": float(b[1]),
                    f"delta_{axes[0]}": float(delta[0]),
                    f"delta_{axes[1]}": float(delta[1]),
                    "determinant": _determinant(a, b),
                    "dot_product": float(np.dot(a, b)),
                    "angle_degrees": _angle_degrees(a, b),
                }
            )
        except (ValueError, np.linalg.LinAlgError, P2Error):
            rows.append({"draw": draw, "valid": False})
    return pd.DataFrame(rows)


def _percentile_interval(series: pd.Series, confidence: float) -> tuple[float, float]:
    x = pd.to_numeric(series, errors="coerce").dropna().to_numpy(float)
    alpha = (1.0 - float(confidence)) / 2.0
    if len(x) == 0:
        return float("nan"), float("nan")
    lo, hi = np.quantile(x, [alpha, 1.0 - alpha])
    return float(lo), float(hi)


def run_p2(
    *,
    artifact_root: Path,
    p2_config_path: Path,
    numerical_spec_path: Path,
    pattern_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    p2 = _load_yaml(p2_config_path)
    numerical = _load_yaml(numerical_spec_path)
    pattern = _load_yaml(pattern_config_path)
    if p2.get("contract") != "chapter1_p2_component_nonconcordance_v1":
        raise P2Error("unexpected P2 contract")
    if numerical.get("contract") != "chapter1_p2_component_nonconcordance_numerical_spec_v1":
        raise P2Error("unexpected P2 numerical contract")
    if str(pattern.get("context_column")) != "analysis_regime":
        raise P2Error("P2 requires the frozen analysis_regime pattern config")

    axes = [str(x) for x in p2["primary_axes"]["display_names"]]
    source_axes = {str(k): str(v) for k, v in p2["primary_axes"]["source_species_scores"].items()}
    scopes = [str(x) for x in p2["scopes"]]
    strata = [str(x) for x in p2["strata"]]
    context_a = str(p2["primary_context"]["context_a"])
    context_b = str(p2["primary_context"]["context_b"])
    threshold = int(p2["minimum_islands_per_axis_per_context"])
    draws = int(numerical["bootstrap"]["draws"])
    seed = int(numerical["bootstrap"]["seed"])
    confidence = float(numerical["bootstrap"]["confidence_level"])
    minimum_valid = int(numerical["bootstrap"]["minimum_valid_draws"])

    covariates = pd.read_csv(artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv")
    status_flora = pd.read_csv(artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz")

    output_dir.mkdir(parents=True, exist_ok=True)
    support_rows: list[dict[str, Any]] = []
    p2a_rows: list[dict[str, Any]] = []
    p2a_cov: list[pd.DataFrame] = []
    p2a_bootstrap: list[pd.DataFrame] = []
    p2a_summary: list[dict[str, Any]] = []
    p2b_rows: list[dict[str, Any]] = []
    p2b_support: list[dict[str, Any]] = []
    realm_rows: list[pd.DataFrame] = []

    for scope_index, scope in enumerate(scopes):
        directory = _scope_dir(scope)
        branch_scores = pd.read_csv(artifact_root / f"branching/{directory}/island_plant_side_branch_scores.csv.gz")
        species_scores = pd.read_csv(artifact_root / f"syndrome/{directory}/species_syndrome_concordance.csv.gz")
        frozen_slopes = pd.read_csv(artifact_root / f"branching/{directory}/global_branch_distance_slopes.csv")
        frozen_between = pd.read_csv(artifact_root / f"branching/{directory}/global_branch_between_context_omnibus.csv")

        realm = frozen_between.loc[
            frozen_between["context_layer"].astype(str).eq("biogeographic_realm")
            & frozen_between["axis_set"].astype(str).eq("universal_plant_response")
            & frozen_between["support_tier"].astype(str).eq("confirmatory")
            & frozen_between["context_a"].astype(str).eq("Palearctic")
            & frozen_between["context_b"].astype(str).eq("Neotropical")
            & frozen_between["stratum"].astype(str).isin(strata)
        ].copy()
        if not realm.empty:
            realm.insert(0, "evidence_scope", scope)
            realm_rows.append(realm)

        coobserved, coobs_diag = _build_coobserved_species_scores(
            species_scores,
            status_flora,
            source_axes=source_axes,
            strata=strata,
            minimum_species=int(p2["p2b_coobserved_species_sensitivity"]["minimum_species_per_island_per_axis"]),
        )

        for stratum_index, stratum in enumerate(strata):
            original = branch_scores.loc[
                branch_scores["stratum"].astype(str).eq(stratum)
                & branch_scores["syndrome"].astype(str).isin(axes)
            ].copy()
            original_context = original.merge(
                covariates[["island_id", "analysis_regime"]].drop_duplicates("island_id"),
                on="island_id",
                how="left",
                validate="many_to_one",
            )
            common_scores = _common_island_scores(branch_scores, axes=axes, stratum=stratum)
            prepared = _prepare(common_scores, covariates, pattern)
            for ctx in [context_a, context_b]:
                for axis in axes:
                    before = original_context.loc[
                        original_context["analysis_regime"].astype(str).eq(ctx)
                        & original_context["syndrome"].astype(str).eq(axis),
                        "island_id",
                    ].nunique()
                    after = prepared.loc[
                        prepared["analysis_regime"].astype(str).eq(ctx)
                        & prepared["syndrome"].astype(str).eq(axis),
                        "island_id",
                    ].nunique()
                    support_rows.append(
                        {
                            "evidence_scope": scope,
                            "stratum": stratum,
                            "context": ctx,
                            "axis": axis,
                            "original_islands": int(before),
                            "common_island_support": int(after),
                        }
                    )

            within_vectors: dict[str, np.ndarray] = {}
            for ctx in [context_a, context_b]:
                slopes, within = _within_context(
                    prepared,
                    stratum=stratum,
                    context_value=ctx,
                    support_tier="confirmatory",
                    threshold=threshold,
                    pattern_config=pattern,
                    syndrome_config={},
                )
                if str(within.get("status", "")) == "fit" and set(slopes["syndrome"].astype(str)) == set(axes):
                    indexed = slopes.set_index("syndrome")
                    within_vectors[ctx] = indexed.loc[axes, "distance_slope"].to_numpy(float)
            between = _between_contexts(
                prepared,
                stratum=stratum,
                context_a=context_a,
                context_b=context_b,
                support_tier="confirmatory",
                threshold=threshold,
                pattern_config=pattern,
            )
            details, cov = _fit_between_details(
                prepared,
                stratum=stratum,
                context_a=context_a,
                context_b=context_b,
                axes=axes,
                threshold=threshold,
                pattern_config=pattern,
            )
            row = {
                "evidence_scope": scope,
                "stratum": stratum,
                **{f"direct_{k}": v for k, v in between.items()},
                **details,
            }
            if context_a in within_vectors and context_b in within_vectors:
                a = within_vectors[context_a]
                b = within_vectors[context_b]
                row.update(
                    {
                        f"{context_a}_{axes[0]}_slope_within": float(a[0]),
                        f"{context_a}_{axes[1]}_slope_within": float(a[1]),
                        f"{context_b}_{axes[0]}_slope_within": float(b[0]),
                        f"{context_b}_{axes[1]}_slope_within": float(b[1]),
                        "determinant": _determinant(a, b),
                        "dot_product": float(np.dot(a, b)),
                        "angle_degrees": _angle_degrees(a, b),
                    }
                )
                frozen = frozen_slopes.loc[
                    frozen_slopes["context_layer"].astype(str).eq("analysis_regime")
                    & frozen_slopes["axis_set"].astype(str).eq("universal_plant_response")
                    & frozen_slopes["support_tier"].astype(str).eq("confirmatory")
                    & frozen_slopes["stratum"].astype(str).eq(stratum)
                    & frozen_slopes["context"].astype(str).isin([context_a, context_b])
                    & frozen_slopes["syndrome"].astype(str).isin(axes)
                ].copy()
                sign_match = True
                for ctx, vec in [(context_a, a), (context_b, b)]:
                    indexed = frozen.loc[frozen["context"].astype(str).eq(ctx)].set_index("syndrome")
                    if set(axes) - set(indexed.index):
                        sign_match = False
                        continue
                    frozen_vec = indexed.loc[axes, "distance_slope"].to_numpy(float)
                    sign_match = sign_match and bool(np.all(np.sign(vec) == np.sign(frozen_vec)))
                row["within_direction_matches_frozen_h2"] = bool(sign_match)
            else:
                row["within_direction_matches_frozen_h2"] = False
            p2a_rows.append(row)
            if not cov.empty:
                cov.insert(0, "stratum", stratum)
                cov.insert(0, "evidence_scope", scope)
                p2a_cov.append(cov)

            if str(details.get("status")) == "fit" and context_a in within_vectors and context_b in within_vectors:
                boot = _bootstrap_profile(
                    prepared,
                    stratum=stratum,
                    context_a=context_a,
                    context_b=context_b,
                    axes=axes,
                    pattern_config=pattern,
                    draws=draws,
                    seed=seed + scope_index * 1000 + stratum_index * 100,
                )
                boot.insert(0, "stratum", stratum)
                boot.insert(0, "evidence_scope", scope)
                p2a_bootstrap.append(boot)
                valid = boot.loc[boot["valid"].eq(True)].copy()
                det_lo, det_hi = _percentile_interval(valid["determinant"], confidence)
                dot_lo, dot_hi = _percentile_interval(valid["dot_product"], confidence)
                angle_lo, angle_hi = _percentile_interval(valid["angle_degrees"], confidence)
                d0_lo, d0_hi = _percentile_interval(valid[f"delta_{axes[0]}"], confidence)
                d1_lo, d1_hi = _percentile_interval(valid[f"delta_{axes[1]}"], confidence)
                p2a_summary.append(
                    {
                        "evidence_scope": scope,
                        "stratum": stratum,
                        "valid_draws": int(len(valid)),
                        "minimum_valid_draws": minimum_valid,
                        "determinant_ci_low": det_lo,
                        "determinant_ci_high": det_hi,
                        "determinant_ci_excludes_zero": bool(det_lo > 0 or det_hi < 0),
                        "dot_product_ci_low": dot_lo,
                        "dot_product_ci_high": dot_hi,
                        "angle_ci_low": angle_lo,
                        "angle_ci_high": angle_hi,
                        f"delta_{axes[0]}_ci_low": d0_lo,
                        f"delta_{axes[0]}_ci_high": d0_hi,
                        f"delta_{axes[1]}_ci_low": d1_lo,
                        f"delta_{axes[1]}_ci_high": d1_hi,
                    }
                )

            co_prepared = _prepare(
                _common_island_scores(coobserved, axes=axes, stratum=stratum), covariates, pattern
            )
            p2b_between = _between_contexts(
                co_prepared,
                stratum=stratum,
                context_a=context_a,
                context_b=context_b,
                support_tier="confirmatory",
                threshold=int(p2["p2b_coobserved_species_sensitivity"]["minimum_islands_per_context"]),
                pattern_config=pattern,
            )
            p2b_details, _ = _fit_between_details(
                co_prepared,
                stratum=stratum,
                context_a=context_a,
                context_b=context_b,
                axes=axes,
                threshold=int(p2["p2b_coobserved_species_sensitivity"]["minimum_islands_per_context"]),
                pattern_config=pattern,
            )
            p2b_rows.append(
                {
                    "evidence_scope": scope,
                    "stratum": stratum,
                    "n_coobserved_species": int(coobs_diag["n_coobserved_species"]),
                    **{f"direct_{k}": v for k, v in p2b_between.items()},
                    **p2b_details,
                }
            )
            for ctx in [context_a, context_b]:
                counts = co_prepared.loc[
                    co_prepared["analysis_regime"].astype(str).eq(ctx)
                ].groupby("syndrome")["island_id"].nunique()
                p2b_support.append(
                    {
                        "evidence_scope": scope,
                        "stratum": stratum,
                        "context": ctx,
                        "n_coobserved_species": int(coobs_diag["n_coobserved_species"]),
                        "accessibility_islands": int(counts.get(axes[0], 0)),
                        "reproductive_islands": int(counts.get(axes[1], 0)),
                        "same_island_denominator": bool(counts.get(axes[0], 0) == counts.get(axes[1], 0)),
                    }
                )

    p2a_df = pd.DataFrame(p2a_rows)
    p2a_summary_df = pd.DataFrame(p2a_summary)
    p2b_df = pd.DataFrame(p2b_rows)
    realm_df = pd.concat(realm_rows, ignore_index=True) if realm_rows else pd.DataFrame()

    primary = p2a_df.loc[
        p2a_df["evidence_scope"].eq("direct_only")
        & p2a_df["stratum"].eq("native_nonendemic")
    ]
    primary_boot = p2a_summary_df.loc[
        p2a_summary_df["evidence_scope"].eq("direct_only")
        & p2a_summary_df["stratum"].eq("native_nonendemic")
    ]
    if len(primary) != 1 or len(primary_boot) != 1:
        raise P2Error("primary P2a profile missing")
    prow = primary.iloc[0]
    brow = primary_boot.iloc[0]
    p2a_supported = bool(
        str(prow.get("status", "")) == "fit"
        and float(prow.get("p_value", np.nan)) <= 0.05
        and bool(prow.get("within_direction_matches_frozen_h2", False))
    )
    strong_nonconcordance = bool(
        p2a_supported
        and int(brow["valid_draws"]) >= minimum_valid
        and bool(brow["determinant_ci_excludes_zero"])
    )

    primary_b = p2b_df.loc[
        p2b_df["evidence_scope"].eq("direct_only")
        & p2b_df["stratum"].eq("native_nonendemic")
    ]
    p2b_executable = bool(
        len(primary_b) == 1 and str(primary_b.iloc[0].get("status", "")) == "fit"
    )
    if strong_nonconcordance and p2b_executable:
        status = "strong_component_nonconcordance_with_common_species_sensitivity_executable"
    elif strong_nonconcordance:
        status = "strong_component_nonconcordance_common_species_sensitivity_underpowered"
    elif p2a_supported:
        status = "vector_branching_retained_but_noncollinearity_not_established"
    else:
        status = "p2_common_support_does_not_defend_h2_branching"

    pd.DataFrame(support_rows).to_csv(output_dir / "p2a_common_island_support.csv", index=False)
    p2a_df.to_csv(output_dir / "p2a_common_island_vector_results.csv", index=False)
    (pd.concat(p2a_cov, ignore_index=True) if p2a_cov else pd.DataFrame()).to_csv(
        output_dir / "p2a_joint_slope_covariance.csv", index=False
    )
    (pd.concat(p2a_bootstrap, ignore_index=True) if p2a_bootstrap else pd.DataFrame()).to_csv(
        output_dir / "p2a_paired_block_bootstrap.csv.gz", index=False
    )
    p2a_summary_df.to_csv(output_dir / "p2a_paired_block_bootstrap_summary.csv", index=False)
    p2b_df.to_csv(output_dir / "p2b_coobserved_species_vector_results.csv", index=False)
    pd.DataFrame(p2b_support).to_csv(output_dir / "p2b_coobserved_species_support.csv", index=False)
    realm_df.to_csv(output_dir / "p2c_realm_boundary_audit.csv", index=False)

    result = {
        "contract": "chapter1_p2_component_nonconcordance_result_v1",
        "status": status,
        "source_primary_analysis": p2["source_primary_analysis"],
        "primary_profile": "direct_only__native_nonendemic",
        "p2a_common_island_vector_supported": p2a_supported,
        "p2a_strong_nonconcordance_determinant_ci_excludes_zero": strong_nonconcordance,
        "p2a_primary_p_value": float(prow.get("p_value", np.nan)),
        "p2a_primary_determinant": float(prow.get("determinant", np.nan)),
        "p2a_primary_determinant_ci": [float(brow["determinant_ci_low"]), float(brow["determinant_ci_high"])],
        "p2a_primary_angle_degrees": float(prow.get("angle_degrees", np.nan)),
        "p2a_primary_angle_ci": [float(brow["angle_ci_low"]), float(brow["angle_ci_high"])],
        "p2a_primary_valid_bootstrap_draws": int(brow["valid_draws"]),
        "p2b_common_species_sensitivity_executable": p2b_executable,
        "p2b_primary_p_value": None if not p2b_executable else float(primary_b.iloc[0].get("p_value", np.nan)),
        "p2c_realm_boundary_role": "claim_boundary_only",
        "claim_boundary": (
            "P2 tests same-layer North-Tropical component geometry. It does not convert the separate "
            "Palearctic within-context result into a supported Palearctic-versus-tropical direct contrast, "
            "does not establish temporal asynchrony, and does not identify a pollinator mechanism."
        ),
    }
    (output_dir / "chapter1_p2_component_nonconcordance_result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    summary = [
        "# Chapter 1 P2 component non-concordance",
        "",
        f"Status: `{status}`",
        "",
        f"Primary common-island direct-only native-nonendemic vector p = {result['p2a_primary_p_value']:.6g}.",
        f"Determinant = {result['p2a_primary_determinant']:.6g}; 95% block-bootstrap interval "
        f"[{result['p2a_primary_determinant_ci'][0]:.6g}, {result['p2a_primary_determinant_ci'][1]:.6g}].",
        f"Angle = {result['p2a_primary_angle_degrees']:.2f} degrees; 95% interval "
        f"[{result['p2a_primary_angle_ci'][0]:.2f}, {result['p2a_primary_angle_ci'][1]:.2f}].",
        f"Common-species sensitivity executable: {p2b_executable}.",
        "",
        "Palearctic-versus-Neotropical realm tests are reported separately as a claim boundary and are not "
        "substituted for the same-layer North-Tropical H2 contrast.",
    ]
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    return result


@app.command()
def run(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    p2_config_path: Path = typer.Option(
        Path("config/chapter1_p2_component_nonconcordance.yml"), exists=True, dir_okay=False
    ),
    numerical_spec_path: Path = typer.Option(
        Path("config/chapter1_p2_component_nonconcordance_numerical_spec.yml"), exists=True, dir_okay=False
    ),
    pattern_config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True, dir_okay=False
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            run_p2(
                artifact_root=artifact_root,
                p2_config_path=p2_config_path,
                numerical_spec_path=numerical_spec_path,
                pattern_config_path=pattern_config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
