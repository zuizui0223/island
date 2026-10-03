"""Post-hoc biological-domain H1-H4 stress test.

This module deliberately does not replace the final confirmatory H1.  It asks whether
the frozen component-trait database tells the same ecological story after traits are
grouped into biologically interpretable domains:

A. reproductive assurance;
B. floral accessibility / specialization;
C. signal / display evidence.

Domain C is only partly directionally identifiable from the current ontology. Flower
size is ordered and can support a strict size-reduction contrast. Colour is hue rather
than measured contrast, and inflorescence_display is architecture rather than total
display area, so those two are retained as explicitly labelled proxies/sensitivities.

Reviewer guards are built into the analysis:
- finite-cluster t and cluster sign-flip inference for island models;
- trait-specific and common-species-support denominators;
- leave-one-species-out genus residual sensitivity;
- strict native / introduced support diagnostics;
- exact-species GloPL matching only;
- finite-publication inference for H4;
- BH correction within declared post-hoc families;
- no causal mediation, pollinator identity, or within-lineage evolutionary claim.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy.stats import t as student_t

from island_v2.chapter1_final_inference_audit import (
    finite_cluster_two_sided_p,
)
from island_v2.chapter1_h1_final_directional import synthesize_regions
from island_v2.chapter1_h1_three_axis_raw_state import build_valid_state_ledger
from island_v2.chapter1_h5_glopl_floral_architecture_moderation import (
    normalize_species_name,
)
from island_v2.chapter1_h5_glopl_global_distance import (
    MEASUREMENT_COLUMNS,
    _clustered_wls,
    _context_dummies,
    _measurement_dummies,
    _site_key,
    assign_context,
    build_study_key,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_biological_domains_reanalysis_v1"


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not bool(ok.any()):
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected biological-domain reanalysis contract")
    return config


def response_specs(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    specs: dict[str, dict[str, Any]] = {}
    for domain in ("reproductive_assurance", "accessibility_specialization"):
        item = config["domains"][domain]
        specs[domain] = {
            "domain": domain,
            "role": str(item["role"]),
            "minimum_components": int(item["minimum_components_per_species"]),
            "components": {
                str(trait): {str(k): float(v) for k, v in mapping.items()}
                for trait, mapping in item["components"].items()
            },
        }
    signal = config["domains"]["signal_display"]
    for response, item in signal["components"].items():
        specs[str(response)] = {
            "domain": "signal_display",
            "role": str(item["role"]),
            "minimum_components": 1,
            "components": {
                str(item["source_trait"]): {
                    str(k): float(v)
                    for k, v in item["values"].items()
                }
            },
        }
    return specs


def build_species_response_scores(
    species_axis: pd.DataFrame,
    ontology: dict[str, Any],
    raw_axis_config: dict[str, Any],
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ledger, coverage = build_valid_state_ledger(
        species_axis,
        ontology,
        raw_axis_config,
        evidence_scope=evidence_scope,
    )
    specs = response_specs(config)
    rows: list[dict[str, Any]] = []
    for response, spec in specs.items():
        for trait, mapping in spec["components"].items():
            part = ledger.loc[
                ledger["trait_name"].astype(str).eq(trait),
                ["accepted_species", "trait_name", "state"],
            ].copy()
            if part.empty:
                continue
            part["state_score"] = part["state"].map(mapping)
            part = part.dropna(subset=["state_score"])
            if part.empty:
                continue
            trait_score = (
                part.groupby(
                    ["accepted_species", "trait_name"],
                    as_index=False,
                )
                .agg(trait_score=("state_score", "mean"))
            )
            trait_score["response"] = response
            trait_score["domain"] = spec["domain"]
            trait_score["role"] = spec["role"]
            rows.extend(trait_score.to_dict(orient="records"))

    trait_scores = pd.DataFrame(rows)
    if trait_scores.empty:
        raise typer.BadParameter("no species response scores could be built")

    grouped = (
        trait_scores.groupby(
            ["accepted_species", "response", "domain", "role"],
            as_index=False,
        )
        .agg(
            score=("trait_score", "mean"),
            n_components=("trait_name", "nunique"),
        )
    )
    minima = {
        response: int(spec["minimum_components"])
        for response, spec in specs.items()
    }
    grouped["minimum_components"] = grouped["response"].map(minima)
    grouped = grouped.loc[
        grouped["n_components"].ge(grouped["minimum_components"])
    ].copy()
    grouped["evidence_scope"] = evidence_scope

    support = (
        grouped.groupby(["response", "domain", "role"], as_index=False)
        .agg(
            n_species=("accepted_species", "nunique"),
            mean_score=("score", "mean"),
            sd_score=("score", "std"),
            median_components=("n_components", "median"),
        )
    )
    support.insert(0, "evidence_scope", evidence_scope)
    return grouped, support


def _flora_scope_mask(status_flora: pd.DataFrame, scope: str) -> pd.Series:
    if scope == "all_observed":
        return pd.Series(True, index=status_flora.index)
    floristic = status_flora.get(
        "floristic_status",
        pd.Series("", index=status_flora.index),
    ).fillna("").astype(str).str.casefold()
    origin = status_flora.get(
        "origin_status",
        pd.Series("", index=status_flora.index),
    ).fillna("").astype(str).str.casefold()
    if scope == "source_native":
        return floristic.isin({"native_nonendemic", "endemic"}) | origin.eq(
            "native"
        )
    if scope == "source_introduced":
        return origin.str.contains("introduced", regex=False) | floristic.str.contains(
            "introduced", regex=False
        )
    raise ValueError(f"unknown flora scope {scope}")


def build_island_scores(
    status_flora: pd.DataFrame,
    species_scores: pd.DataFrame,
    config: dict[str, Any],
    *,
    flora_scope: str,
    support_mode: str,
) -> pd.DataFrame:
    mask = _flora_scope_mask(status_flora, flora_scope)
    flora = status_flora.loc[
        mask,
        ["island_id", "accepted_species"],
    ].copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.drop_duplicates(["island_id", "accepted_species"])

    scores = species_scores.copy()
    if support_mode == "common_species":
        required = [str(x) for x in config["common_support_responses"]]
        wide = scores.pivot_table(
            index="accepted_species",
            columns="response",
            values="score",
            aggfunc="first",
        )
        common_species = wide.dropna(subset=required).index.astype(str)
        scores = scores.loc[
            scores["accepted_species"].astype(str).isin(common_species)
            & scores["response"].astype(str).isin(required)
        ].copy()
    elif support_mode != "trait_specific":
        raise ValueError(f"unknown support mode {support_mode}")

    joined = flora.merge(
        scores[
            [
                "accepted_species",
                "response",
                "domain",
                "role",
                "score",
            ]
        ],
        on="accepted_species",
        how="inner",
        validate="many_to_many",
    )
    out = (
        joined.groupby(
            ["island_id", "response", "domain", "role"],
            as_index=False,
        )
        .agg(
            island_score=("score", "mean"),
            n_scored_species=("accepted_species", "nunique"),
        )
    )
    out = out.loc[
        out["n_scored_species"].ge(
            int(config["minimum_species_per_island_score"])
        )
    ].copy()
    out["flora_scope"] = flora_scope
    out["support_mode"] = support_mode
    return out


def build_genus_residual_species_scores(
    species_scores: pd.DataFrame,
) -> pd.DataFrame:
    work = species_scores.copy()
    work["genus"] = (
        work["accepted_species"].astype(str).str.split().str[0].fillna("")
    )
    stats = (
        work.groupby(["response", "genus"], as_index=False)
        .agg(
            genus_sum=("score", "sum"),
            genus_n=("score", "size"),
        )
    )
    work = work.merge(
        stats,
        on=["response", "genus"],
        how="left",
        validate="many_to_one",
    )
    valid = work["genus_n"].gt(1) & work["genus"].ne("")
    work = work.loc[valid].copy()
    work["loo_genus_mean"] = (
        work["genus_sum"] - work["score"]
    ) / (work["genus_n"] - 1.0)
    work["score"] = work["score"] - work["loo_genus_mean"]
    work["role"] = work["role"].astype(str) + "_genus_LOO_residual"
    return work.drop(columns=["genus_sum", "genus_n", "loo_genus_mean"])


def _wild_signflip_p(
    estimate: float,
    se: float,
    influences: np.ndarray,
    *,
    replications: int,
    seed: int,
) -> float:
    if (
        not math.isfinite(estimate)
        or not math.isfinite(se)
        or se <= 0
        or len(influences) < 2
    ):
        return float("nan")
    rng = np.random.default_rng(seed)
    observed = estimate / se
    exceed = 0
    done = 0
    while done < replications:
        n = min(2000, replications - done)
        signs = rng.choice(
            np.array([-1.0, 1.0]),
            size=(n, len(influences)),
        )
        t_star = (signs @ influences) / se
        exceed += int(np.count_nonzero(t_star >= observed))
        done += n
    return float((exceed + 1.0) / (replications + 1.0))


def fit_clustered_score(
    frame: pd.DataFrame,
    *,
    response_column: str,
    predictors: list[str],
    target_predictor: str,
    cluster_column: str,
    replications: int,
    seed: int,
) -> dict[str, Any]:
    required = [response_column, *predictors, cluster_column]
    work = frame[required].copy()
    for column in [response_column, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cluster_column] = work[cluster_column].fillna("").astype(str)
    work = work.dropna(subset=[response_column, *predictors])
    work = work.loc[work[cluster_column].ne("")].copy()
    if len(work) < 20:
        return {"status": "not_testable", "n_islands": int(len(work))}

    names = ["intercept"]
    columns = [np.ones(len(work), dtype=float)]
    for predictor in predictors:
        try:
            columns.append(_z(work[predictor]))
        except ValueError:
            return {
                "status": "not_testable",
                "reason": f"constant_predictor:{predictor}",
                "n_islands": int(len(work)),
            }
        names.append(f"z_{predictor}")
    X = np.column_stack(columns)
    y = pd.to_numeric(work[response_column], errors="coerce").to_numpy(float)
    bread = np.linalg.pinv(X.T @ X)
    beta = bread @ X.T @ y
    residual = y - X @ beta

    labels = work[cluster_column].astype(str).to_numpy()
    unique = np.unique(labels)
    target_name = f"z_{target_predictor}"
    target_index = names.index(target_name)
    influences: list[float] = []
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for label in unique:
        mask = labels == label
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
        influences.append(float(bread[target_index] @ score))
    n = len(work)
    k = X.shape[1]
    g = len(unique)
    correction = 1.0
    if g > 1 and n > k:
        correction = (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    covariance = correction * bread @ meat @ bread
    estimate = float(beta[target_index])
    se = float(
        math.sqrt(max(float(covariance[target_index, target_index]), 0.0))
    )
    t_value = estimate / se if se > 0 else float("nan")
    df = g - 1
    p_one = (
        float(student_t.sf(t_value, df=df))
        if df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    p_two = (
        float(2.0 * student_t.sf(abs(t_value), df=df))
        if df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    inf = np.asarray(influences, dtype=float) * math.sqrt(correction)
    p_wild = _wild_signflip_p(
        estimate,
        se,
        inf,
        replications=replications,
        seed=seed,
    )
    squared = np.square(inf)
    total = float(np.sum(squared))
    effective_clusters = (
        float(total * total / np.sum(np.square(squared)))
        if total > 0 and np.sum(np.square(squared)) > 0
        else float("nan")
    )
    max_share = (
        float(np.max(squared) / total)
        if total > 0
        else float("nan")
    )
    return {
        "status": "fit",
        "estimate": estimate,
        "se": se,
        "t_value": t_value,
        "df": int(df),
        "p_one_sided_positive": p_one,
        "p_two_sided": p_two,
        "p_wild_positive": p_wild,
        "n_islands": int(n),
        "n_clusters": int(g),
        "effective_clusters": effective_clusters,
        "max_cluster_variance_share": max_share,
    }


def run_h1(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    flora_scope: str,
    support_mode: str,
    analysis_layer: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    context_col = str(config["context_column"])
    cluster_col = str(config["cluster_column"])
    geography = str(config["geography_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    cov_cols = [
        "island_id",
        geography,
        context_col,
        cluster_col,
        *baseline,
    ]
    work = island_scores.merge(
        covariates[cov_cols].drop_duplicates("island_id"),
        on="island_id",
        how="left",
    )
    rows: list[dict[str, Any]] = []
    counter = 0
    for response in sorted(work["response"].astype(str).unique()):
        for context in config["contexts"]:
            part = work.loc[
                work["response"].astype(str).eq(response)
                & work[context_col].astype(str).eq(str(context))
            ].copy()
            if part["island_id"].nunique() < int(config["minimum_islands"]):
                rows.append(
                    {
                        "evidence_scope": evidence_scope,
                        "flora_scope": flora_scope,
                        "support_mode": support_mode,
                        "analysis_layer": analysis_layer,
                        "response": response,
                        "context": context,
                        "status": "not_testable",
                        "n_islands": int(part["island_id"].nunique()),
                    }
                )
                continue
            result = fit_clustered_score(
                part,
                response_column="island_score",
                predictors=[geography, *baseline],
                target_predictor=geography,
                cluster_column=cluster_col,
                replications=int(config["wild_cluster_replications"]),
                seed=int(config["wild_cluster_seed"]) + counter,
            )
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "flora_scope": flora_scope,
                    "support_mode": support_mode,
                    "analysis_layer": analysis_layer,
                    "response": response,
                    "context": context,
                    **result,
                }
            )
            counter += 1
    regional = pd.DataFrame(rows)

    synth_rows: list[dict[str, Any]] = []
    alpha = float(config["alpha"])
    for response, part in regional.loc[
        regional["status"].eq("fit")
    ].groupby("response"):
        meta_input = part.rename(
            columns={"se": "cluster_robust_se"}
        )
        synthesis = synthesize_regions(
            meta_input,
            contexts=[str(x) for x in config["contexts"]],
            alpha=alpha,
        )
        synth_rows.append(
            {
                "evidence_scope": evidence_scope,
                "flora_scope": flora_scope,
                "support_mode": support_mode,
                "analysis_layer": analysis_layer,
                "response": str(response),
                **synthesis,
            }
        )
    synthesis = pd.DataFrame(synth_rows)
    if not synthesis.empty and "H1a_one_sided_p" in synthesis:
        synthesis["H1a_q_all_responses"] = _bh(
            synthesis["H1a_one_sided_p"]
        )
        core = synthesis["response"].isin(
            [str(x) for x in config["primary_responses"]]
        )
        synthesis["H1a_q_core_two"] = np.nan
        synthesis.loc[core, "H1a_q_core_two"] = _bh(
            synthesis.loc[core, "H1a_one_sided_p"]
        )
    return regional, synthesis


def run_h2(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    support_mode: str,
) -> pd.DataFrame:
    wide = island_scores.pivot_table(
        index="island_id",
        columns="response",
        values="island_score",
        aggfunc="first",
    ).reset_index()
    geography = str(config["geography_column"])
    context_col = str(config["context_column"])
    cluster_col = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    cov_cols = [
        "island_id",
        geography,
        context_col,
        cluster_col,
        *baseline,
    ]
    work = wide.merge(
        covariates[cov_cols].drop_duplicates("island_id"),
        on="island_id",
        how="left",
    )
    rows: list[dict[str, Any]] = []
    counter = 0
    for model in config["H2_models"]:
        response = str(model["response"])
        conditions = [str(x) for x in model["condition_on"]]
        if response not in work.columns or any(
            x not in work.columns for x in conditions
        ):
            continue
        for context in config["contexts"]:
            part = work.loc[
                work[context_col].astype(str).eq(str(context))
            ].copy()
            result = fit_clustered_score(
                part,
                response_column=response,
                predictors=[geography, *conditions, *baseline],
                target_predictor=geography,
                cluster_column=cluster_col,
                replications=int(config["wild_cluster_replications"]),
                seed=int(config["wild_cluster_seed"]) + 5000 + counter,
            )
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "support_mode": support_mode,
                    "flora_scope": "all_observed",
                    "response": response,
                    "condition_on": "|".join(conditions),
                    "analysis_role": str(model["role"]),
                    "context": context,
                    **result,
                }
            )
            counter += 1
    out = pd.DataFrame(rows)
    if not out.empty:
        fit = out["status"].eq("fit")
        out["q_two_sided_posthoc_family"] = np.nan
        out.loc[fit, "q_two_sided_posthoc_family"] = _bh(
            out.loc[fit, "p_two_sided"]
        )
    return out


def prepare_glopl_rows(
    glopl_csv: Path,
    corrected_site_distances: Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    if _sha256(glopl_csv) != str(config["GloPL"]["sha256"]):
        raise typer.BadParameter("GloPL SHA-256 mismatch")
    columns = [
        "Species_accepted_names",
        "Latitude",
        "Longitude",
        "DOI",
        "Author",
        "Year",
        "PL_Effect_Size",
        *MEASUREMENT_COLUMNS,
    ]
    raw = pd.read_csv(
        glopl_csv,
        usecols=columns,
        dtype=str,
        encoding="latin-1",
    ).fillna("")
    raw["site_key"] = _site_key(raw["Latitude"], raw["Longitude"])
    raw["study_key"] = build_study_key(raw)
    raw["analysis_regime"] = assign_context(
        raw["Latitude"],
        {
            "contexts": {
                "tropical_abs_lat_lt": 23.4366,
                "northern_high_lat_ge": 66.5634,
            }
        },
    )
    raw["species_key"] = raw["Species_accepted_names"].map(
        normalize_species_name
    )
    raw["PL_Effect_Size"] = pd.to_numeric(
        raw["PL_Effect_Size"], errors="coerce"
    )
    distances = pd.read_csv(corrected_site_distances)
    distances["site_key"] = distances["site_key"].astype(str)
    work = raw.merge(
        distances[["site_key", "spherical_distance_km"]],
        on="site_key",
        how="left",
        validate="many_to_one",
    )
    work["log_distance"] = np.log1p(
        pd.to_numeric(
            work["spherical_distance_km"], errors="coerce"
        )
    )
    mean = float(config["GloPL"]["distance_mean"])
    sd = float(config["GloPL"]["distance_sd"])
    work["z_distance"] = (work["log_distance"] - mean) / sd
    valid = (
        work["PL_Effect_Size"].notna()
        & work["z_distance"].notna()
        & work["study_key"].astype(str).ne("")
        & work["site_key"].astype(str).ne("")
        & work["species_key"].astype(str).ne("")
        & work["analysis_regime"].isin(config["contexts"])
    )
    return work.loc[valid].copy()


def _h4_cells(
    glopl_rows: pd.DataFrame,
    species_scores: pd.DataFrame,
    *,
    response: str,
    publication_total_weight: float,
) -> pd.DataFrame:
    scores = species_scores.loc[
        species_scores["response"].astype(str).eq(response),
        ["accepted_species", "score"],
    ].copy()
    scores["species_key"] = scores["accepted_species"].map(
        normalize_species_name
    )
    scores = scores.dropna(subset=["species_key", "score"])
    conflicts = scores.groupby("species_key")["score"].nunique().gt(1)
    if bool(conflicts.any()):
        raise typer.BadParameter(
            f"conflicting exact species scores for {response}"
        )
    scores = scores[["species_key", "score"]].drop_duplicates("species_key")
    joined = glopl_rows.merge(
        scores,
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
        "score",
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
    cells["analysis_weight"] = publication_total_weight / counts
    return cells


def _fit_h4_score(cells: pd.DataFrame) -> dict[str, Any]:
    if cells.empty or cells["score"].nunique() < 2:
        return {"status": "not_testable"}
    distance = pd.to_numeric(cells["z_distance"], errors="coerce").to_numpy(float)
    score = pd.to_numeric(cells["score"], errors="coerce").to_numpy(float)
    context_cols, context_names, _ = _context_dummies(cells)
    measure_cols, measure_names = _measurement_dummies(cells)
    X = np.column_stack(
        [
            np.ones(len(cells), dtype=float),
            *context_cols,
            distance,
            score,
            *measure_cols,
        ]
    )
    names = [
        "intercept",
        *context_names,
        "z_distance",
        "score",
        *measure_names,
    ]
    fit = _clustered_wls(cells, X, names)
    if not fit.get("evaluable"):
        return {"status": "not_testable", "reason": fit.get("reason", "")}
    idx = fit["names"].index("score")
    estimate = float(fit["beta"][idx])
    se = float(fit["se"][idx])
    df = int(fit["n_publications"]) - 1
    t_value = estimate / se if se > 0 else float("nan")
    p_two = (
        float(2.0 * student_t.sf(abs(t_value), df=df))
        if df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    p_negative = (
        float(student_t.cdf(t_value, df=df))
        if df > 0 and math.isfinite(t_value)
        else float("nan")
    )
    return {
        "status": "fit",
        "estimate": estimate,
        "se": se,
        "df": df,
        "p_two_sided_finite_publication": p_two,
        "p_one_sided_negative_finite_publication": p_negative,
        "n_cells": int(fit["n_cells"]),
        "n_publications": int(fit["n_publications"]),
        "n_sites": int(fit["n_sites"]),
        "n_species": int(cells["species_key"].nunique()),
    }


def run_h4(
    glopl_rows: pd.DataFrame,
    direct_species_scores: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    responses = [
        *[str(x) for x in config["GloPL"]["h4_primary_responses"]],
        *[str(x) for x in config["GloPL"]["h4_exploratory_responses"]],
    ]
    rows: list[dict[str, Any]] = []
    for response in responses:
        cells = _h4_cells(
            glopl_rows,
            direct_species_scores,
            response=response,
            publication_total_weight=float(
                config["GloPL"]["publication_total_weight"]
            ),
        )
        for analysis, part in {
            "primary": cells,
            "supplemental_only": cells.loc[
                cells["PL_Effect_Size_Type2"].astype(str).eq("Sup")
            ].copy(),
            "no_zero_constant": cells.loc[
                ~cells["Constant_added"].astype(str).str.strip().str.casefold().isin(
                    {"true", "1", "yes", "y", "t"}
                )
            ].copy(),
        }.items():
            rows.append(
                {
                    "response": response,
                    "analysis": analysis,
                    "role": (
                        "posthoc_exact_species_functional_triangulation"
                        if response in config["GloPL"]["h4_primary_responses"]
                        else "exploratory_signal_display_functional_association"
                    ),
                    **_fit_h4_score(part),
                }
            )
    out = pd.DataFrame(rows)
    primary = out["analysis"].eq("primary") & out["status"].eq("fit")
    out["q_two_sided_primary_family"] = np.nan
    out.loc[primary, "q_two_sided_primary_family"] = _bh(
        out.loc[primary, "p_two_sided_finite_publication"]
    )
    return out


def h3_summary(h3_json: Path) -> dict[str, Any]:
    payload = json.loads(h3_json.read_text(encoding="utf-8"))
    item = payload["corrected"]["global_gradient"]
    estimate = float(item["distance_slope"])
    se = float(item["distance_slope_se"])
    n_publications = int(item["n_publications"])
    return {
        "estimate": estimate,
        "se": se,
        "n_publications": n_publications,
        "finite_publication_two_sided_p": finite_cluster_two_sided_p(
            estimate,
            se,
            n_publications,
        ),
        "classification_dependency": False,
        "interpretation": (
            "H3 is unchanged by trait-domain reclassification and remains an "
            "independent isolation-associated pollen-limitation gradient."
        ),
    }


def summarize(
    h1_synthesis: pd.DataFrame,
    h2: pd.DataFrame,
    h3: dict[str, Any],
    h4: pd.DataFrame,
) -> dict[str, Any]:
    h1_primary = h1_synthesis.loc[
        h1_synthesis["flora_scope"].eq("all_observed")
        & h1_synthesis["support_mode"].eq("trait_specific")
        & h1_synthesis["analysis_layer"].eq("observed")
    ].copy()
    h2_fit = h2.loc[
        h2["status"].eq("fit")
        & h2["support_mode"].eq("trait_specific")
    ].copy()
    h4_primary = h4.loc[h4["analysis"].eq("primary")].copy()
    return {
        "contract": CONTRACT,
        "inferential_role": "posthoc_biological_reclassification_stress_test",
        "H1_global": json.loads(
            h1_primary.to_json(orient="records", double_precision=15)
        ),
        "H2_supported_cells_q_le_0_05": json.loads(
            h2_fit.loc[
                h2_fit["q_two_sided_posthoc_family"].le(0.05)
            ].to_json(orient="records", double_precision=15)
        ),
        "H3": h3,
        "H4_primary": json.loads(
            h4_primary.to_json(orient="records", double_precision=15)
        ),
        "claim_boundary": {
            "replaces_confirmatory_H1": False,
            "causal_mediation_identified": False,
            "named_pollinator_mechanism_identified": False,
            "within_lineage_evolution_identified": False,
            "signal_display_has_single_confirmatory_direction": False,
        },
    }


@app.command("run")
def run(
    species_axis_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    glopl_csv: Path = typer.Option(..., exists=True),
    ontology_path: Path = typer.Option(..., exists=True),
    raw_axis_config_path: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    h3_json: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    raw_config = yaml.safe_load(
        raw_axis_config_path.read_text(encoding="utf-8")
    )
    ontology = yaml.safe_load(ontology_path.read_text(encoding="utf-8"))
    species_axis = pd.read_csv(species_axis_csv, dtype=str).fillna("")
    status_flora = pd.read_csv(status_flora_csv, dtype=str).fillna("")
    covariates = pd.read_csv(covariates_csv)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_regional: list[pd.DataFrame] = []
    all_synthesis: list[pd.DataFrame] = []
    all_h2: list[pd.DataFrame] = []
    support_parts: list[pd.DataFrame] = []
    species_scores_by_scope: dict[str, pd.DataFrame] = {}

    for evidence_scope in config["evidence_scopes"]:
        species_scores, support = build_species_response_scores(
            species_axis,
            ontology,
            raw_config,
            config,
            evidence_scope=str(evidence_scope),
        )
        species_scores_by_scope[str(evidence_scope)] = species_scores
        support_parts.append(support)

        for flora_scope in (
            "all_observed",
            "source_native",
            "source_introduced",
        ):
            trait_scores = build_island_scores(
                status_flora,
                species_scores,
                config,
                flora_scope=flora_scope,
                support_mode="trait_specific",
            )
            regional, synthesis = run_h1(
                trait_scores,
                covariates,
                config,
                evidence_scope=str(evidence_scope),
                flora_scope=flora_scope,
                support_mode="trait_specific",
                analysis_layer="observed",
            )
            all_regional.append(regional)
            all_synthesis.append(synthesis)

        common = build_island_scores(
            status_flora,
            species_scores,
            config,
            flora_scope="all_observed",
            support_mode="common_species",
        )
        regional, synthesis = run_h1(
            common,
            covariates,
            config,
            evidence_scope=str(evidence_scope),
            flora_scope="all_observed",
            support_mode="common_species",
            analysis_layer="observed",
        )
        all_regional.append(regional)
        all_synthesis.append(synthesis)

        residual_species = build_genus_residual_species_scores(
            species_scores
        )
        residual_island = build_island_scores(
            status_flora,
            residual_species,
            config,
            flora_scope="all_observed",
            support_mode="trait_specific",
        )
        regional, synthesis = run_h1(
            residual_island,
            covariates,
            config,
            evidence_scope=str(evidence_scope),
            flora_scope="all_observed",
            support_mode="trait_specific",
            analysis_layer="genus_LOO_residual",
        )
        all_regional.append(regional)
        all_synthesis.append(synthesis)

        trait_scores = build_island_scores(
            status_flora,
            species_scores,
            config,
            flora_scope="all_observed",
            support_mode="trait_specific",
        )
        all_h2.append(
            run_h2(
                trait_scores,
                covariates,
                config,
                evidence_scope=str(evidence_scope),
                support_mode="trait_specific",
            )
        )
        all_h2.append(
            run_h2(
                common,
                covariates,
                config,
                evidence_scope=str(evidence_scope),
                support_mode="common_species",
            )
        )

    regional = pd.concat(all_regional, ignore_index=True)
    synthesis = pd.concat(all_synthesis, ignore_index=True)
    h2 = pd.concat(all_h2, ignore_index=True)
    support = pd.concat(support_parts, ignore_index=True)

    glopl_rows = prepare_glopl_rows(
        glopl_csv,
        Path(config["GloPL"]["corrected_site_distances"]),
        config,
    )
    h4 = run_h4(
        glopl_rows,
        species_scores_by_scope["direct_only"],
        config,
    )
    h3 = h3_summary(h3_json)
    summary = summarize(synthesis, h2, h3, h4)

    support.to_csv(output_dir / "species_score_support.csv", index=False)
    regional.to_csv(output_dir / "H1_regional.csv", index=False)
    synthesis.to_csv(output_dir / "H1_global_synthesis.csv", index=False)
    h2.to_csv(output_dir / "H2_conditional.csv", index=False)
    h4.to_csv(output_dir / "H4_exact_species.csv", index=False)
    (output_dir / "H3_summary.json").write_text(
        json.dumps(h3, indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
