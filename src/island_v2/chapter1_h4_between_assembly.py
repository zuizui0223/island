"""Taxonomic depth of the North--Tropical H3 moderation difference.

This companion to ``chapter1_h4_atomic_assembly`` asks whether the direct
North--Tropical distance x area x context vector survives family and genus adjustment
on exactly the same source-matched observed species.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import _bh
from island_v2.chapter1_h4_atomic_assembly import (
    _fit_clustered_ols,
    _joint,
    _z_masked,
    run_h4,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def fit_between_stage(
    decomposition: pd.DataFrame,
    covariates: pd.DataFrame,
    ladder_config: dict[str, Any],
    h4_config: dict[str, Any],
    *,
    source_mode: str,
    stratum: str,
    stage: str,
    context_a: str,
    context_b: str,
) -> dict[str, Any]:
    outcomes = [str(x) for x in h4_config["atomic_outcomes"]]
    context = str(ladder_config["context_column"])
    cluster = str(ladder_config["cluster_column"])
    distance = str(ladder_config["geography_column"])
    area = str(ladder_config["area_column"])
    controls = [str(x) for x in ladder_config["control_columns"]]
    threshold = int(h4_config["minimum_support"]["islands_per_context_outcome"])

    work = decomposition.loc[
        decomposition["source_mode"].eq(source_mode)
        & decomposition["stratum"].eq(stratum)
        & decomposition["syndrome"].isin(outcomes)
    ].copy()
    required_cov = ["island_id", context, cluster, distance, area, *controls]
    work = work.merge(
        covariates[required_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    work = work.loc[work[context].astype(str).isin([context_a, context_b])].copy()
    for column in [stage, distance, area, *controls]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna(subset=[stage, distance, area, *controls])
    support = work.groupby(["syndrome", context])["island_id"].nunique().unstack(fill_value=0)
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
    if len(retained) < 2:
        return {
            "source_mode": source_mode,
            "stratum": stratum,
            "stage": stage,
            "context_a": context_a,
            "context_b": context_b,
            "status": "not_testable",
            "n_retained_outcomes": len(retained),
            "retained_outcomes": "|".join(retained),
        }

    work = work.loc[work["syndrome"].isin(retained)].reset_index(drop=True)
    b = work[context].astype(str).eq(context_b).to_numpy(float)
    columns: list[np.ndarray] = []
    names: list[str] = []
    triple_indices: list[int] = []
    for outcome in retained:
        mask = work["syndrome"].eq(outcome).to_numpy()
        indicator = mask.astype(float)
        columns.extend([indicator, indicator * b])
        names.extend([f"{outcome}:intercept", f"{outcome}:context[{context_b}]"])
        for predictor in controls:
            z = _z_masked(work, mask, predictor)
            columns.extend([z, z * b])
            names.extend(
                [f"{outcome}:z_{predictor}", f"{outcome}:z_{predictor}:context[{context_b}]"]
            )
        za = _z_masked(work, mask, area)
        zd = _z_masked(work, mask, distance)
        columns.extend([za, za * b, zd, zd * b, za * zd, za * zd * b])
        names.extend(
            [
                f"{outcome}:z_{area}",
                f"{outcome}:z_{area}:context[{context_b}]",
                f"{outcome}:z_{distance}",
                f"{outcome}:z_{distance}:context[{context_b}]",
                f"{outcome}:z_{distance}:z_{area}",
                f"{outcome}:z_{distance}:z_{area}:context[{context_b}]",
            ]
        )
        triple_indices.append(len(names) - 1)

    beta, covariance, n_clusters = _fit_clustered_ols(
        work[stage].to_numpy(float),
        np.column_stack(columns),
        work[cluster].to_numpy(str),
    )
    vector = beta[triple_indices]
    cov = covariance[np.ix_(triple_indices, triple_indices)]
    statistic, df, p_value = _joint(vector, cov)
    return {
        "source_mode": source_mode,
        "stratum": stratum,
        "stage": stage,
        "context_a": context_a,
        "context_b": context_b,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(n_clusters),
        "triple_joint_wald": statistic,
        "triple_joint_df": df,
        "triple_p_value": p_value,
    }


def classify_between_depth(stage_results: pd.DataFrame, alpha: float) -> pd.DataFrame:
    work = stage_results.loc[stage_results["status"].eq("fit")].copy()
    if work.empty:
        return pd.DataFrame()
    work["triple_q_across_source_modes"] = work.groupby(
        ["stratum", "stage"], group_keys=False
    )["triple_p_value"].transform(_bh)
    work["triple_supported"] = work["triple_q_across_source_modes"].le(alpha)
    rows = []
    for (stratum, source_mode), group in work.groupby(["stratum", "source_mode"], sort=True):
        status = group.set_index("stage")["triple_supported"].to_dict()
        observed = bool(status.get("observed_score", False))
        family = bool(status.get("after_family_residual", False))
        genus = bool(status.get("after_genus_residual", False))
        if not observed:
            classification = "observed_between_context_not_supported"
        elif not family:
            classification = "between_context_difference_compatible_with_family_or_deeper_sorting"
        elif not genus:
            classification = "between_context_difference_compatible_with_genus_level_assembly"
        else:
            classification = "between_context_residual_below_genus_retained"
        rows.append(
            {
                "stratum": stratum,
                "source_mode": source_mode,
                "classification": classification,
            }
        )
    per_mode = pd.DataFrame(rows)
    aggregate = []
    for stratum, group in per_mode.groupby("stratum", sort=True):
        values = sorted(set(group["classification"]))
        aggregate.append(
            {
                "stratum": stratum,
                "source_mode_classifications": "|".join(
                    f"{row.source_mode}:{row.classification}" for row in group.itertuples()
                ),
                "robust_classification": values[0] if len(values) == 1 else "source_definition_sensitive",
            }
        )
    return per_mode.merge(pd.DataFrame(aggregate), on="stratum", how="left")


@app.command("run")
def run(
    state_audit_csv: Path = typer.Option(..., exists=True),
    taxonomy_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    source_assignments_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    ladder_config_path: Path = typer.Option(..., exists=True),
    h4_config_path: Path = typer.Option(..., exists=True),
    source_config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    probability_config = yaml.safe_load(probability_config_path.read_text(encoding="utf-8"))
    ladder_config = yaml.safe_load(ladder_config_path.read_text(encoding="utf-8"))
    h4_config = yaml.safe_load(h4_config_path.read_text(encoding="utf-8"))
    source_config = yaml.safe_load(source_config_path.read_text(encoding="utf-8"))
    _, decomposition, _, _, _ = run_h4(
        pd.read_csv(state_audit_csv),
        pd.read_csv(taxonomy_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(gift_flora_csv),
        pd.read_csv(source_assignments_csv),
        pd.read_csv(covariates_csv),
        probability_config,
        ladder_config,
        h4_config,
        source_config,
    )
    covariates = pd.read_csv(covariates_csv)
    contexts = [str(x) for x in h4_config["contexts"]]
    context_a, context_b = contexts[0], contexts[1]
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]
    strata = [str(h4_config["primary_stratum"]), *[str(x) for x in h4_config["sensitivity_strata"]]]
    rows = []
    for source_mode in source_modes:
        for stratum in strata:
            for stage in [str(x) for x in h4_config["stages"]]:
                rows.append(
                    fit_between_stage(
                        decomposition,
                        covariates,
                        ladder_config,
                        h4_config,
                        source_mode=source_mode,
                        stratum=stratum,
                        stage=stage,
                        context_a=context_a,
                        context_b=context_b,
                    )
                )
    results = pd.DataFrame(rows)
    classification = classify_between_depth(results, float(ladder_config["alpha"]))
    for frame in (results, classification):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "h4_between_stage_results.csv", index=False)
    classification.to_csv(output_dir / "h4_between_depth_classification.csv", index=False)


if __name__ == "__main__":
    app()
