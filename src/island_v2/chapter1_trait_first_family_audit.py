"""Compare trait-first family summaries with the historical species-first concordance route.

The audit keeps the predeclared biological trait directions and historical family
weights, but changes the order of aggregation. Atomic trait shares are first estimated
within each island and then combined into accessibility/generalisation and reproductive
assurance family scores. This preserves trait-specific denominators before family
compression and avoids species-level renormalisation over heterogeneous missing traits.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import build_broad_counts
from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.chapter1_context_analysis import _chi_square_sf_integer_df

app = typer.Typer(add_completion=False, no_args_is_help=True)

FAMILY_WEIGHTS = {
    "accessibility_generalization": {
        "generalized_form": 1.0,
        "actinomorphic_symmetry": 0.75,
        "shallow_open_tube": 1.0,
    },
    "reproductive_assurance": {
        "self_compatibility": 1.0,
        "selfing_mating_system": 1.0,
        "autonomous_selfing": 1.0,
    },
}


def build_trait_first_scores(
    counts: pd.DataFrame,
    *,
    require_all_components: bool,
) -> pd.DataFrame:
    work = counts.loc[counts["stratum"].eq("all_observed")].copy()
    work["share"] = pd.to_numeric(work["successes"], errors="coerce") / pd.to_numeric(
        work["trials"], errors="coerce"
    )
    rows: list[dict[str, Any]] = []
    for island_id, island in work.groupby("island_id", sort=False):
        indexed = island.drop_duplicates("outcome").set_index("outcome")
        for family, weights in FAMILY_WEIGHTS.items():
            available = [name for name in weights if name in indexed.index and np.isfinite(indexed.loc[name, "share"])]
            if require_all_components and len(available) != len(weights):
                continue
            if not available:
                continue
            denom = float(sum(weights[name] for name in available))
            score = float(sum(weights[name] * float(indexed.loc[name, "share"]) for name in available) / denom)
            trials = int(sum(int(indexed.loc[name, "trials"]) for name in available))
            rows.append(
                {
                    "island_id": str(island_id),
                    "family": family,
                    "family_score": score,
                    "n_component_traits": len(available),
                    "sum_trait_trials": trials,
                }
            )
    return pd.DataFrame(rows)


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant predictor")
    return (x - mean) / sd


def fit_north_tropical_family_difference(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    context = str(config["context_column"])
    geography = str(config["geography_column"])
    cluster = str(config["cluster_column"])
    baseline = [str(x) for x in config["baseline_covariates"]]
    required = ["island_id", context, geography, cluster, *baseline]
    data = scores.merge(
        covariates[required].drop_duplicates("island_id"),
        on="island_id", how="left", validate="many_to_one",
    )
    numeric = ["family_score", geography, *baseline]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna(subset=numeric)
    data = data.loc[
        data[context].isin(["northern_midlatitude", "tropical"])
        & data[cluster].fillna("").astype(str).ne("")
    ].copy()

    families = list(FAMILY_WEIGHTS)
    support = data.groupby(["family", context])["island_id"].nunique().unstack(fill_value=0)
    threshold = int(config["minimum_islands_per_outcome"])
    retained = [
        family for family in families
        if family in support.index
        and int(support.loc[family, "northern_midlatitude"]) >= threshold
        and int(support.loc[family, "tropical"]) >= threshold
    ]
    if len(retained) != 2:
        return pd.DataFrame(), {"status": "not_testable", "retained": "|".join(retained)}
    work = data.loc[data["family"].isin(retained)].copy()
    tropical = work[context].eq("tropical").to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    interactions: list[str] = []
    for family in retained:
        mask = work["family"].eq(family).to_numpy()
        indicator = mask.astype(float)
        names.extend([f"{family}:intercept", f"{family}:context[tropical]"])
        columns.extend([indicator, indicator * tropical])
        for predictor in baseline:
            z = np.zeros(len(work), dtype=float)
            z[mask] = _standardize(work.loc[mask, predictor])
            names.extend([f"{family}:z_{predictor}", f"{family}:z_{predictor}:context[tropical]"])
            columns.extend([z, z * tropical])
        z_geo = np.zeros(len(work), dtype=float)
        z_geo[mask] = _standardize(work.loc[mask, geography])
        interaction = f"{family}:z_{geography}:context[tropical]"
        names.extend([f"{family}:z_{geography}", interaction])
        columns.extend([z_geo, z_geo * tropical])
        interactions.append(interaction)

    coef, covariance, fit = _fit_weighted_clustered_design(
        work["family_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster].astype(str).to_numpy(),
    )
    if coef.empty:
        return coef, {"status": fit.get("status", "fit_failed")}
    indexed = coef.set_index("predictor")
    indices = [names.index(name) for name in interactions]
    vector = np.array([float(indexed.loc[name, "estimate"]) for name in interactions])
    vcov = covariance[np.ix_(indices, indices)]
    rank = int(np.linalg.matrix_rank(vcov))
    stat = float(vector @ np.linalg.pinv(vcov) @ vector) if rank > 0 else float("nan")
    rows = []
    for family, name in zip(retained, interactions, strict=True):
        r = indexed.loc[name]
        rows.append(
            {
                "family": family,
                "slope_difference_tropical_minus_north": float(r["estimate"]),
                "cluster_robust_se": float(r["cluster_robust_se"]),
                "p_value": float(r["p_value"]),
                "n_islands_north": int(support.loc[family, "northern_midlatitude"]),
                "n_islands_tropical": int(support.loc[family, "tropical"]),
            }
        )
    return pd.DataFrame(rows), {
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster].nunique()),
        "joint_wald_chisq": stat,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(stat, rank),
    }


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    counts = build_broad_counts(
        pd.read_csv(status_flora_csv), pd.read_csv(state_audit_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries=[]
    for label, require_all in (("available_components", False), ("complete_three_components", True)):
        scores = build_trait_first_scores(counts, require_all_components=require_all)
        component, summary = fit_north_tropical_family_difference(
            scores, pd.read_csv(covariates_csv), config
        )
        scores.to_csv(output_dir / f"trait_first_scores_{label}.csv.gz", index=False, compression="gzip")
        component.to_csv(output_dir / f"trait_first_components_{label}.csv", index=False)
        summaries.append({"evidence_scope": evidence_scope, "support_rule": label, **summary})
    summary_frame=pd.DataFrame(summaries)
    summary_frame.to_csv(output_dir / "trait_first_family_audit.csv", index=False)
    typer.echo(summary_frame.to_csv(index=False))


if __name__ == "__main__":
    app()
