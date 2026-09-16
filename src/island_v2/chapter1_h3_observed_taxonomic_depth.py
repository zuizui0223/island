"""Source-free taxonomic-depth decomposition for the broad Chapter 1 H2 response.

This analysis fills the missing all-observed H2 x taxonomic-depth cell. It keeps the
same observed island-species records across observed/family/genus stages and uses
leave-one-species-out (LOO) family and genus means estimated from the fixed scored
species pool. No GIFT source assignment or floristic-status restriction enters the
primary decomposition.

The result is a representation-depth diagnostic. It cannot identify native assembly,
introduction, historical colonisation, in-situ evolution, or a pollinator mechanism.
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

from island_v2.chapter1_all_data_probability import (
    _bh,
    _chi_square_sf_integer_df,
    _classify_signature,
    _truthy,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

STAGES = ("observed_score", "after_family_residual", "after_genus_residual")


def build_atomic_taxonomic_residuals(
    state_audit: pd.DataFrame,
    taxonomy: pd.DataFrame,
    probability_config: dict[str, Any],
    depth_config: dict[str, Any],
) -> pd.DataFrame:
    required_audit = {
        "accepted_species",
        "trait_name",
        "resolved_for_primary",
        "canonical_signature",
    }
    required_taxonomy = {"accepted_species", "family", "genus"}
    if missing := required_audit - set(state_audit.columns):
        raise typer.BadParameter(f"state audit missing columns: {sorted(missing)}")
    if missing := required_taxonomy - set(taxonomy.columns):
        raise typer.BadParameter(f"taxonomy missing columns: {sorted(missing)}")

    tax = taxonomy[["accepted_species", "family", "genus"]].copy().fillna("")
    for column in tax.columns:
        tax[column] = tax[column].astype(str).str.strip()
    tax = tax.loc[tax["accepted_species"].ne("")].drop_duplicates("accepted_species")
    if tax["accepted_species"].duplicated().any():
        raise typer.BadParameter("taxonomy contains duplicate accepted_species")

    audit = state_audit.loc[_truthy(state_audit["resolved_for_primary"])].copy()
    audit["accepted_species"] = audit["accepted_species"].astype(str)
    audit["trait_name"] = audit["trait_name"].astype(str)

    outcomes = [str(x) for x in depth_config["atomic_outcomes"]]
    minimum = int(
        depth_config["species_taxonomic_expectation"][
            "minimum_scored_species_per_group_including_focal"
        ]
    )
    pieces: list[pd.DataFrame] = []
    for outcome in outcomes:
        spec = probability_config["broad_outcomes"][outcome]
        part = audit.loc[
            audit["trait_name"].eq(str(spec["trait_name"])),
            ["accepted_species", "canonical_signature"],
        ].drop_duplicates("accepted_species")
        positive = {str(x) for x in spec["positive_states"]}
        negative = {str(x) for x in spec["negative_states"]}
        part = part.copy()
        part["state"] = [
            _classify_signature(value, positive, negative)
            for value in part["canonical_signature"]
        ]
        part = part.dropna(subset=["state"])[["accepted_species", "state"]]
        part = part.merge(tax, on="accepted_species", how="left", validate="many_to_one")
        part = part.loc[part["family"].ne("") & part["genus"].ne("")].copy()
        if part.empty:
            continue

        for level in ("family", "genus"):
            stats = part.groupby(level)["state"].agg(["sum", "count"])
            part[f"{level}_sum"] = part[level].map(stats["sum"])
            part[f"{level}_n"] = part[level].map(stats["count"]).astype(int)
            denom = part[f"{level}_n"] - 1
            part[f"{level}_loo_mean"] = np.where(
                part[f"{level}_n"].ge(minimum) & denom.gt(0),
                (part[f"{level}_sum"] - part["state"]) / denom,
                np.nan,
            )

        part = part.dropna(subset=["family_loo_mean", "genus_loo_mean"]).copy()
        if part.empty:
            continue
        part["after_family_residual"] = part["state"] - part["family_loo_mean"]
        part["after_genus_residual"] = part["state"] - part["genus_loo_mean"]
        part["observed_score"] = part["state"]
        part["outcome"] = outcome
        pieces.append(
            part[
                [
                    "accepted_species",
                    "outcome",
                    "family",
                    "genus",
                    "observed_score",
                    "family_loo_mean",
                    "genus_loo_mean",
                    "after_family_residual",
                    "after_genus_residual",
                    "family_n",
                    "genus_n",
                ]
            ]
        )
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def build_island_stage_scores(
    status_flora: pd.DataFrame,
    species_residuals: pd.DataFrame,
) -> pd.DataFrame:
    required = {"island_id", "accepted_species"}
    if missing := required - set(status_flora.columns):
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    flora = status_flora[["island_id", "accepted_species"]].drop_duplicates().copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    pieces: list[pd.DataFrame] = []
    for outcome, species in species_residuals.groupby("outcome", sort=True):
        joined = flora.merge(species, on="accepted_species", how="inner", validate="many_to_one")
        if joined.empty:
            continue
        summary = (
            joined.groupby("island_id", as_index=False)
            .agg(
                observed_score=("observed_score", "mean"),
                after_family_residual=("after_family_residual", "mean"),
                after_genus_residual=("after_genus_residual", "mean"),
                n_species=("accepted_species", "nunique"),
                n_families=("family", "nunique"),
                n_genera=("genus", "nunique"),
            )
        )
        summary["outcome"] = str(outcome)
        pieces.append(summary)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _fit_ols_component(
    y: np.ndarray,
    design: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
) -> dict[str, Any]:
    X = np.asarray(design, dtype=float)
    y = np.asarray(y, dtype=float)
    beta = np.linalg.pinv(X.T @ X) @ (X.T @ y)
    residual = y - X @ beta
    bread = np.linalg.pinv(X.T @ X)
    score = X * residual[:, None]
    return {
        "theta": beta,
        "bread": bread,
        "score": score,
        "names": names,
        "clusters": np.asarray(clusters).astype(str),
    }


def _assemble_cluster_covariance(
    fits: list[dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    offsets: list[int] = []
    names: list[str] = []
    size = 0
    for fit in fits:
        offsets.append(size)
        names.extend(fit["names"])
        size += len(fit["names"])
    bread = np.zeros((size, size), dtype=float)
    scores: dict[str, np.ndarray] = {}
    total_rows = 0
    for fit, offset in zip(fits, offsets, strict=True):
        dim = len(fit["names"])
        bread[offset : offset + dim, offset : offset + dim] = fit["bread"]
        labels = fit["clusters"]
        total_rows += len(labels)
        for cluster in np.unique(labels):
            if cluster not in scores:
                scores[cluster] = np.zeros(size, dtype=float)
            scores[cluster][offset : offset + dim] += fit["score"][labels == cluster].sum(axis=0)
    meat = np.zeros((size, size), dtype=float)
    for score in scores.values():
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    g = len(scores)
    if g > 1 and total_rows > size:
        covariance *= (g / (g - 1.0)) * ((total_rows - 1.0) / (total_rows - size))
    theta = np.concatenate([fit["theta"] for fit in fits])
    return theta, covariance, names


def _prepare_model_data(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    depth_config: dict[str, Any],
) -> pd.DataFrame:
    context = str(depth_config["context_column"])
    cluster = str(depth_config["cluster_column"])
    distance = str(depth_config["distance_column"])
    controls = [str(x) for x in depth_config["controls"]]
    required = {"island_id", context, cluster, distance, *controls}
    if missing := required - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    data = island_scores.merge(
        covariates[list(required)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = [*STAGES, "n_species", distance, *controls]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context] = data[context].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=[*STAGES, distance, *controls])
    contexts = [str(x) for x in depth_config["contexts"]]
    return data.loc[data[context].isin(contexts) & data[cluster].ne("")].copy()


def fit_pairwise_stage(
    data: pd.DataFrame,
    depth_config: dict[str, Any],
    stage: str,
) -> tuple[pd.DataFrame, dict[str, Any], np.ndarray]:
    context = str(depth_config["context_column"])
    cluster = str(depth_config["cluster_column"])
    distance = str(depth_config["distance_column"])
    controls = [str(x) for x in depth_config["controls"]]
    context_a, context_b = [str(x) for x in depth_config["contexts"]]
    threshold = int(depth_config["island_model"]["minimum_islands_per_context_outcome"])
    outcomes = [str(x) for x in depth_config["atomic_outcomes"]]

    support = data.groupby(["outcome", context])["island_id"].nunique().unstack(fill_value=0)
    for value in (context_a, context_b):
        if value not in support.columns:
            support[value] = 0
    retained = [
        outcome
        for outcome in outcomes
        if outcome in support.index
        and int(support.loc[outcome, context_a]) >= threshold
        and int(support.loc[outcome, context_b]) >= threshold
    ]
    if len(retained) < int(depth_config["island_model"]["minimum_outcomes_per_vector"]):
        return pd.DataFrame(), {
            "stage": stage,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }, np.array([], dtype=float)

    fits: list[dict[str, Any]] = []
    interaction_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = data.loc[data["outcome"].eq(outcome)].copy()
        indicator = part[context].eq(context_b).to_numpy(float)
        columns = [np.ones(len(part), dtype=float), indicator]
        names = [f"{outcome}:intercept", f"{outcome}:context[{context_b}]"]
        for predictor in controls:
            z = _z(part[predictor])
            columns.extend([z, z * indicator])
            names.extend(
                [
                    f"{outcome}:z_{predictor}",
                    f"{outcome}:z_{predictor}:context[{context_b}]",
                ]
            )
        z_distance = _z(part[distance])
        columns.extend([z_distance, z_distance * indicator])
        interaction_name = f"{outcome}:z_{distance}:context[{context_b}]"
        names.extend([f"{outcome}:z_{distance}", interaction_name])
        fit = _fit_ols_component(
            part[stage].to_numpy(float),
            np.column_stack(columns),
            names,
            part[cluster].to_numpy(str),
        )
        fits.append(fit)
        interaction_indices.append(offset + names.index(interaction_name))
        rows.append(
            {
                "stage": stage,
                "outcome": outcome,
                "context_a": context_a,
                "context_b": context_b,
                "n_islands_context_a": int(support.loc[outcome, context_a]),
                "n_islands_context_b": int(support.loc[outcome, context_b]),
                "median_species_support": float(part["n_species"].median()),
            }
        )
        offset += len(names)

    theta, covariance, _ = _assemble_cluster_covariance(fits)
    vector = theta[interaction_indices]
    vector_cov = covariance[np.ix_(interaction_indices, interaction_indices)]
    stderr = np.sqrt(np.clip(np.diag(vector_cov), 0.0, None))
    for row, estimate, se in zip(rows, vector, stderr, strict=True):
        z_value = float(estimate / se) if se > 0 else float("nan")
        row.update(
            {
                "slope_difference_b_minus_a": float(estimate),
                "cluster_robust_se": float(se),
                "p_value": math.erfc(abs(z_value) / math.sqrt(2.0))
                if math.isfinite(z_value)
                else float("nan"),
            }
        )
    rank = int(np.linalg.matrix_rank(vector_cov))
    statistic = float(vector @ np.linalg.pinv(vector_cov) @ vector) if rank > 0 else float("nan")
    omnibus = {
        "stage": stage,
        "status": "fit",
        "context_a": context_a,
        "context_b": context_b,
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(data.loc[data["outcome"].isin(retained), "island_id"].nunique()),
        "n_clusters": int(data.loc[data["outcome"].isin(retained), cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(statistic, rank),
        "vector_norm": float(np.linalg.norm(vector)),
    }
    return pd.DataFrame(rows), omnibus, vector


def _bootstrap_attenuation(
    data: pd.DataFrame,
    depth_config: dict[str, Any],
) -> pd.DataFrame:
    spec = depth_config["attenuation"]["bootstrap"]
    if not bool(spec.get("enabled", False)):
        return pd.DataFrame()
    cluster = str(depth_config["cluster_column"])
    draws = int(spec["draws"])
    rng = np.random.default_rng(int(spec["seed"]))
    labels = sorted(data[cluster].astype(str).unique())
    by_cluster = {label: data.loc[data[cluster].astype(str).eq(label)].copy() for label in labels}
    rows: list[dict[str, Any]] = []
    for draw in range(draws):
        sampled = rng.choice(labels, size=len(labels), replace=True)
        pieces: list[pd.DataFrame] = []
        for replicate, label in enumerate(sampled):
            part = by_cluster[str(label)].copy()
            part[cluster] = f"{label}__boot{replicate}"
            pieces.append(part)
        boot = pd.concat(pieces, ignore_index=True)
        vectors: dict[str, np.ndarray] = {}
        ok = True
        for stage in STAGES:
            _, result, vector = fit_pairwise_stage(boot, depth_config, stage)
            if result.get("status") != "fit" or len(vector) == 0:
                ok = False
                break
            vectors[stage] = vector
        if not ok:
            continue
        observed = float(np.linalg.norm(vectors["observed_score"]))
        family = float(np.linalg.norm(vectors["after_family_residual"]))
        genus = float(np.linalg.norm(vectors["after_genus_residual"]))
        if observed <= 0:
            continue
        rows.append(
            {
                "draw": draw,
                "observed_norm": observed,
                "family_norm": family,
                "genus_norm": genus,
                "family_attenuation": 1.0 - family / observed,
                "genus_attenuation": 1.0 - genus / observed,
            }
        )
    return pd.DataFrame(rows)


def classify_depth(omnibus: pd.DataFrame, alpha: float) -> str:
    if omnibus.empty:
        return "not_testable"
    indexed = omnibus.set_index("stage")
    if any(stage not in indexed.index for stage in STAGES):
        return "not_testable"
    observed = float(indexed.loc["observed_score", "p_value"])
    family = float(indexed.loc["after_family_residual", "p_value"])
    genus = float(indexed.loc["after_genus_residual", "p_value"])
    if not math.isfinite(observed) or observed >= alpha:
        return "observed_common_support_not_reproduced"
    if not math.isfinite(family) or family >= alpha:
        return "compatible_with_family_or_deeper_structure"
    if not math.isfinite(genus) or genus >= alpha:
        return "compatible_with_genus_structuring"
    return "below_genus_or_non_taxonomic_residual_retained"


def run_analysis(
    state_audit: pd.DataFrame,
    taxonomy: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    depth_config: dict[str, Any],
    evidence_scope: str,
) -> dict[str, Any]:
    residuals = build_atomic_taxonomic_residuals(
        state_audit, taxonomy, probability_config, depth_config
    )
    island_scores = build_island_stage_scores(status_flora, residuals)
    data = _prepare_model_data(island_scores, covariates, depth_config)
    slope_parts: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    vectors: dict[str, np.ndarray] = {}
    for stage in STAGES:
        slopes, result, vector = fit_pairwise_stage(data, depth_config, stage)
        if not slopes.empty:
            slope_parts.append(slopes)
        omnibus_rows.append(result)
        vectors[stage] = vector
    omnibus = pd.DataFrame(omnibus_rows)
    slopes = pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame()

    norms = {
        stage: float(np.linalg.norm(vector)) if len(vector) else float("nan")
        for stage, vector in vectors.items()
    }
    observed_norm = norms["observed_score"]
    attenuation = pd.DataFrame(
        [
            {
                "evidence_scope": evidence_scope,
                "observed_norm": observed_norm,
                "family_norm": norms["after_family_residual"],
                "genus_norm": norms["after_genus_residual"],
                "family_attenuation": (
                    1.0 - norms["after_family_residual"] / observed_norm
                    if math.isfinite(observed_norm) and observed_norm > 0
                    else float("nan")
                ),
                "genus_attenuation": (
                    1.0 - norms["after_genus_residual"] / observed_norm
                    if math.isfinite(observed_norm) and observed_norm > 0
                    else float("nan")
                ),
            }
        ]
    )

    bootstrap = _bootstrap_attenuation(data, depth_config)
    if not bootstrap.empty:
        for column in ("family_attenuation", "genus_attenuation"):
            attenuation[f"{column}_bootstrap_median"] = float(bootstrap[column].median())
            attenuation[f"{column}_bootstrap_ci_low"] = float(bootstrap[column].quantile(0.025))
            attenuation[f"{column}_bootstrap_ci_high"] = float(bootstrap[column].quantile(0.975))
        attenuation["bootstrap_valid_draws"] = int(len(bootstrap))

    classification = classify_depth(
        omnibus, float(depth_config["interpretation_gate"]["alpha"])
    )
    support = (
        residuals.groupby("outcome", as_index=False)
        .agg(
            n_species_common_support=("accepted_species", "nunique"),
            n_families=("family", "nunique"),
            n_genera=("genus", "nunique"),
            median_family_size=("family_n", "median"),
            median_genus_size=("genus_n", "median"),
        )
    )
    manifest = {
        "contract": depth_config["contract"],
        "evidence_scope": evidence_scope,
        "flora_scope": "all_observed",
        "source_pool_used": False,
        "taxonomic_expectation": "leave_one_species_out",
        "common_species_support_across_stages": True,
        "classification": classification,
        "n_species_outcome_rows": int(len(residuals)),
        "n_island_outcome_rows": int(len(island_scores)),
        "n_unique_islands": int(island_scores["island_id"].nunique()) if not island_scores.empty else 0,
        "claim_boundary": (
            "Taxonomic representation depth only; no native assembly, introduction, "
            "historical source filtering, in-situ evolution or pollinator cause is identified."
        ),
    }
    return {
        "species_residuals": residuals,
        "island_scores": island_scores,
        "support": support,
        "slopes": slopes,
        "omnibus": omnibus,
        "attenuation": attenuation,
        "bootstrap": bootstrap,
        "manifest": manifest,
    }


@app.command("run")
def run(
    state_audit_csv: Path = typer.Option(..., exists=True),
    taxonomy_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    depth_config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    probability_config = yaml.safe_load(probability_config_path.read_text(encoding="utf-8"))
    depth_config = yaml.safe_load(depth_config_path.read_text(encoding="utf-8"))
    if depth_config.get("contract") != "chapter1_h3_observed_taxonomic_depth_v1":
        raise typer.BadParameter("unexpected H3 observed-taxonomic-depth contract")
    outputs = run_analysis(
        pd.read_csv(state_audit_csv),
        pd.read_csv(taxonomy_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        probability_config,
        depth_config,
        evidence_scope,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs["species_residuals"].to_csv(
        output_dir / "h3_species_taxonomic_residuals.csv.gz", index=False, compression="gzip"
    )
    outputs["island_scores"].to_csv(
        output_dir / "h3_island_stage_scores.csv.gz", index=False, compression="gzip"
    )
    outputs["support"].to_csv(output_dir / "h3_common_support.csv", index=False)
    outputs["slopes"].to_csv(output_dir / "h3_stage_pairwise_outcome_slopes.csv", index=False)
    outputs["omnibus"].to_csv(output_dir / "h3_stage_pairwise_omnibus.csv", index=False)
    outputs["attenuation"].to_csv(output_dir / "h3_attenuation_summary.csv", index=False)
    outputs["bootstrap"].to_csv(output_dir / "h3_attenuation_bootstrap.csv.gz", index=False, compression="gzip")
    (output_dir / "h3_manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n", encoding="utf-8"
    )
    omnibus = outputs["omnibus"].set_index("stage") if not outputs["omnibus"].empty else pd.DataFrame()
    attenuation = outputs["attenuation"].iloc[0].to_dict() if not outputs["attenuation"].empty else {}
    summary = [
        "# H3 broad observed-flora taxonomic-depth result",
        "",
        f"- evidence scope: `{evidence_scope}`",
        f"- classification: **{outputs['manifest']['classification']}**",
        f"- common-support islands: **{outputs['manifest']['n_unique_islands']}**",
    ]
    for stage in STAGES:
        if not omnibus.empty and stage in omnibus.index:
            row = omnibus.loc[stage]
            summary.append(
                f"- {stage}: p={float(row['p_value']):.6g}, vector norm={float(row['vector_norm']):.6g}"
            )
    if attenuation:
        summary.extend(
            [
                f"- family attenuation: {float(attenuation.get('family_attenuation', float('nan'))):.4f}",
                f"- genus attenuation: {float(attenuation.get('genus_attenuation', float('nan'))):.4f}",
            ]
        )
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    typer.echo(json.dumps(outputs["manifest"], indent=2))


if __name__ == "__main__":
    app()
