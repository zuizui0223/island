"""Primary PR136 H2a test: source-channel exposure x geography on genus-fixed residual traits.

This module tests the necessary-condition prediction in PR #136:

    Island isolation does not create a universal floral syndrome; it creates a
    syndrome only when it removes a pollination channel that the regional flora
    could otherwise use.

The response is the island-level genus-fixed residual trait share from
``genus_fixed_trait_null.py``. Source-channel exposure is outcome-blind and comes
from the reviewed applicability registry. ``applicable`` means the focal channel
is regionally/source available; ``structurally_not_applicable`` means the source
region lacks that channel and is therefore a falsification regime, not a deficit.

The primary H2a test is the joint vector of geography x source-exposure
interactions across outcomes that meet the same island-support threshold in both
source-exposed and structurally-absent groups. A significant slope in one group
and a nonsignificant slope in the other is never counted as evidence of
contingency.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _joint_wald(beta: np.ndarray, covariance: np.ndarray) -> tuple[float, int, float]:
    rank = int(np.linalg.matrix_rank(covariance))
    if rank <= 0:
        return float("nan"), 0, float("nan")
    stat = float(beta @ np.linalg.pinv(covariance) @ beta)
    return stat, rank, _chi_square_sf_integer_df(stat, rank)


def _fit_weighted_clustered_design(
    y: np.ndarray,
    weights: np.ndarray,
    design: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    n, p = design.shape
    unique_clusters = np.unique(clusters)
    n_clusters = int(len(unique_clusters))
    if n < max(10, p + 3) or n_clusters < 2:
        return (
            pd.DataFrame(),
            np.empty((0, 0)),
            {"status": "insufficient_complete_rows", "n_rows": n, "n_clusters": n_clusters},
        )

    xtwx = design.T @ (weights[:, None] * design)
    bread = np.linalg.pinv(xtwx)
    beta = bread @ (design.T @ (weights * y))
    residual = y - design @ beta

    meat = np.zeros((p, p), dtype=float)
    for cluster in unique_clusters:
        mask = clusters == cluster
        score = design[mask].T @ (weights[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if n_clusters > 1 and n > p:
        covariance *= (n_clusters / (n_clusters - 1.0)) * ((n - 1.0) / (n - p))

    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    rows: list[dict[str, Any]] = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        rows.append(
            {
                "predictor": name,
                "estimate": float(estimate),
                "cluster_robust_se": float(stderr),
                "z_value": z,
                "p_value": _normal_two_sided_p(z) if math.isfinite(z) else float("nan"),
            }
        )
    return (
        pd.DataFrame(rows),
        covariance,
        {"status": "fit", "n_rows": n, "n_clusters": n_clusters},
    )


def _prepare(
    genus_null: pd.DataFrame,
    applicability: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required_null = {
        "island_id",
        "outcome",
        "stratum",
        "observed_n_species",
        "deviation_observed_minus_null",
    }
    required_app = {"island_id", "applicability"}
    geography = str(config["geography_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    required_cov = {"island_id", geography, cluster, *baseline}

    missing = required_null - set(genus_null.columns)
    if missing:
        raise typer.BadParameter(f"genus-null table missing columns: {sorted(missing)}")
    missing = required_app - set(applicability.columns)
    if missing:
        raise typer.BadParameter(f"applicability table missing columns: {sorted(missing)}")
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariate table missing columns: {sorted(missing)}")

    app_table = applicability[["island_id", "applicability"]].drop_duplicates("island_id").copy()
    if app_table["island_id"].duplicated().any():
        raise typer.BadParameter("applicability table must have one row per island")
    app_table["applicability"] = app_table["applicability"].fillna("").astype(str)
    app_table["source_exposed"] = np.where(
        app_table["applicability"].eq("applicable"),
        1.0,
        np.where(app_table["applicability"].eq("structurally_not_applicable"), 0.0, np.nan),
    )

    data = genus_null.merge(app_table, on="island_id", how="left", validate="many_to_one")
    data = data.merge(
        covariates[["island_id", geography, cluster, *baseline]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    numeric = [
        "observed_n_species",
        "deviation_observed_minus_null",
        "source_exposed",
        geography,
        *baseline,
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[cluster] = data[cluster].fillna("").astype(str)
    data = data.dropna(subset=numeric)
    data = data.loc[data["observed_n_species"].gt(0) & data[cluster].ne("")].copy()
    return data


def _fit_stratum_tier(
    data: pd.DataFrame,
    *,
    stratum: str,
    support_tier: str,
    threshold: int,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    geography = str(config["geography_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    min_outcomes = int(config.get("minimum_outcomes_per_vector", 2))

    work = data.loc[data["stratum"].astype(str).eq(stratum)].copy()
    counts = (
        work.groupby(["outcome", "source_exposed"])["island_id"]
        .nunique()
        .unstack(fill_value=0)
    )
    for exposure in (0.0, 1.0):
        if exposure not in counts.columns:
            counts[exposure] = 0
    retained = sorted(
        counts.index[counts[0.0].ge(threshold) & counts[1.0].ge(threshold)].astype(str)
    )

    support_rows: list[dict[str, Any]] = []
    for outcome in sorted(work["outcome"].dropna().astype(str).unique()):
        n_structural = int(counts.loc[outcome, 0.0]) if outcome in counts.index else 0
        n_exposed = int(counts.loc[outcome, 1.0]) if outcome in counts.index else 0
        support_rows.append(
            {
                "stratum": stratum,
                "support_tier": support_tier,
                "threshold": threshold,
                "outcome": outcome,
                "n_structurally_absent_islands": n_structural,
                "n_source_exposed_islands": n_exposed,
                "retained_in_joint_vector": outcome in retained,
            }
        )
    support = pd.DataFrame(support_rows)

    if len(retained) < min_outcomes:
        return (
            pd.DataFrame(),
            support,
            {
                "stratum": stratum,
                "support_tier": support_tier,
                "threshold": threshold,
                "status": "not_testable",
                "n_retained_outcomes": len(retained),
                "retained_outcomes": "|".join(retained),
            },
        )

    work = work.loc[work["outcome"].astype(str).isin(retained)].copy()
    names: list[str] = []
    columns: list[np.ndarray] = []
    interaction_names: list[str] = []
    geography_names: list[str] = []
    exposed_slope_vectors: list[tuple[str, str]] = []

    for outcome in retained:
        mask = work["outcome"].astype(str).eq(outcome).to_numpy()
        indicator = mask.astype(float)
        names.append(f"outcome[{outcome}]")
        columns.append(indicator)

        for predictor in baseline:
            vector = np.zeros(len(work), dtype=float)
            vector[mask] = _standardize(work.loc[mask, predictor])
            names.append(f"outcome[{outcome}]:z_{predictor}")
            columns.append(vector)

        z_geography = np.zeros(len(work), dtype=float)
        z_geography[mask] = _standardize(work.loc[mask, geography])
        geo_name = f"outcome[{outcome}]:z_{geography}"
        names.append(geo_name)
        columns.append(z_geography)
        geography_names.append(geo_name)

        source = work["source_exposed"].to_numpy(float)
        source_name = f"outcome[{outcome}]:source_exposed"
        names.append(source_name)
        columns.append(indicator * source)

        interaction_name = f"outcome[{outcome}]:z_{geography}:source_exposed"
        names.append(interaction_name)
        columns.append(z_geography * source)
        interaction_names.append(interaction_name)
        exposed_slope_vectors.append((geo_name, interaction_name))

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["deviation_observed_minus_null"].to_numpy(float),
        work["observed_n_species"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return (
            coefficients,
            support,
            {
                "stratum": stratum,
                "support_tier": support_tier,
                "threshold": threshold,
                "status": fit["status"],
                "n_retained_outcomes": len(retained),
                "retained_outcomes": "|".join(retained),
                "n_clusters": int(fit.get("n_clusters", 0)),
            },
        )

    coefficients.insert(0, "support_tier", support_tier)
    coefficients.insert(0, "stratum", stratum)

    beta = coefficients.set_index("predictor")["estimate"]
    interaction_indices = [names.index(x) for x in interaction_names]
    interaction_vector = np.array([float(beta[x]) for x in interaction_names])
    interaction_cov = covariance[np.ix_(interaction_indices, interaction_indices)]
    joint_stat, joint_df, joint_p = _joint_wald(interaction_vector, interaction_cov)

    slope_rows: list[dict[str, Any]] = []
    for outcome, (geo_name, interaction_name) in zip(
        retained, exposed_slope_vectors, strict=True
    ):
        geo_idx = names.index(geo_name)
        int_idx = names.index(interaction_name)
        structural = float(beta[geo_name])
        structural_se = float(math.sqrt(max(covariance[geo_idx, geo_idx], 0.0)))
        exposed = structural + float(beta[interaction_name])
        exposed_var = float(
            covariance[geo_idx, geo_idx]
            + covariance[int_idx, int_idx]
            + 2.0 * covariance[geo_idx, int_idx]
        )
        exposed_se = math.sqrt(max(exposed_var, 0.0))
        slope_rows.extend(
            [
                {
                    "stratum": stratum,
                    "support_tier": support_tier,
                    "outcome": outcome,
                    "channel_regime": "structurally_absent",
                    "geography_slope": structural,
                    "cluster_robust_se": structural_se,
                },
                {
                    "stratum": stratum,
                    "support_tier": support_tier,
                    "outcome": outcome,
                    "channel_regime": "source_exposed",
                    "geography_slope": exposed,
                    "cluster_robust_se": exposed_se,
                },
            ]
        )

    result = {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_structurally_absent_islands": int(
            work.loc[work["source_exposed"].eq(0.0), "island_id"].nunique()
        ),
        "n_source_exposed_islands": int(
            work.loc[work["source_exposed"].eq(1.0), "island_id"].nunique()
        ),
        "n_clusters": int(fit["n_clusters"]),
        "joint_exposure_interaction_chisq": joint_stat,
        "joint_exposure_interaction_df": joint_df,
        "joint_exposure_interaction_p": joint_p,
        "source_exposure_contingency_supported": bool(
            math.isfinite(joint_p) and joint_p <= float(config.get("alpha", 0.05))
        ),
    }
    return coefficients, support, {**result, "slope_rows": slope_rows}


def run_source_exposure_test(
    genus_null: pd.DataFrame,
    applicability: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = _prepare(genus_null, applicability, covariates, config)
    strata = [str(x) for x in config["strata"]]
    support_tiers = {str(k): int(v) for k, v in config["support_tiers"].items()}

    coefficient_parts: list[pd.DataFrame] = []
    support_parts: list[pd.DataFrame] = []
    joint_rows: list[dict[str, Any]] = []
    slope_rows: list[dict[str, Any]] = []

    for stratum in strata:
        for tier, threshold in support_tiers.items():
            coef, support, result = _fit_stratum_tier(
                data,
                stratum=stratum,
                support_tier=tier,
                threshold=threshold,
                config=config,
            )
            if not coef.empty:
                coefficient_parts.append(coef)
            support_parts.append(support)
            slope_rows.extend(result.pop("slope_rows", []))
            joint_rows.append(result)

    return (
        pd.concat(coefficient_parts, ignore_index=True)
        if coefficient_parts
        else pd.DataFrame(),
        pd.concat(support_parts, ignore_index=True) if support_parts else pd.DataFrame(),
        pd.DataFrame(joint_rows),
        pd.DataFrame(slope_rows),
    )


@app.command("run")
def run(
    genus_null_csv: Path = typer.Option(..., exists=True),
    applicability_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_pr136_source_exposure.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    coefficients, support, joint, slopes = run_source_exposure_test(
        pd.read_csv(genus_null_csv),
        pd.read_csv(applicability_csv),
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output_dir / "source_exposure_coefficients.csv", index=False)
    support.to_csv(output_dir / "source_exposure_support.csv", index=False)
    joint.to_csv(output_dir / "source_exposure_joint_test.csv", index=False)
    slopes.to_csv(output_dir / "source_exposure_regime_slopes.csv", index=False)

    confirmatory = joint.loc[joint["support_tier"].astype(str).eq("confirmatory")]
    manifest = {
        "contract": "chapter1_pr136_source_exposure_v1",
        "response": "observed_direct_trait_share_minus_same_genus_null_expectation",
        "source_exposed": "applicability == applicable",
        "structurally_absent": "applicability == structurally_not_applicable",
        "unresolved_applicability_excluded": True,
        "primary_test": "joint geography x source-channel-exposure interaction vector",
        "invalid_shortcut": "significant in one regime and nonsignificant in another",
        "n_confirmatory_rows": int(len(confirmatory)),
        "n_confirmatory_testable": int(confirmatory["status"].eq("fit").sum()),
        "n_confirmatory_supported": int(
            confirmatory.get(
                "source_exposure_contingency_supported", pd.Series(dtype=bool)
            )
            .fillna(False)
            .sum()
        ),
        "claim_ceiling": (
            "H2a necessary-condition test only; source exposure is not realized deficit "
            "and does not establish historical pollinator loss"
        ),
    }
    (output_dir / "chapter1_pr136_source_exposure_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
