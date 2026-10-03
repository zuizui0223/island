# ruff: noqa: B008
"""Taxonomic representation depth of the current seven-response Chapter 1 H1.

This analysis is a localization diagnostic, not a causal decomposition. It asks whether
the current within-region isolation vectors remain after subtracting species-level
leave-one-species-out family and genus expectations on the exact same taxonomically
eligible species.

Two flora scopes are kept separate:
- all_observed;
- WCVP regional_native_compatible.

The current beta-binomial H1 is first re-fitted on the exact common taxonomic support.
Only then are observed/family-residual/genus-residual island means analysed with the
continuous clustered decomposition required for signed taxonomic residuals.
"""
from __future__ import annotations

import argparse
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2 import chapter1_all_data_probability as h1
from island_v2.chapter1_h3_observed_taxonomic_depth import (
    STAGES,
    _assemble_cluster_covariance,
    _fit_ols_component,
    _z,
    build_atomic_taxonomic_residuals,
    build_island_stage_scores,
)
from island_v2.chapter1_wcvp_native_compatibility import (
    build_regional_native_compatible_flora,
)


def _read_yaml(path: Path) -> dict[str, Any]:
    obj = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"invalid YAML object: {path}")
    return obj


def _common_counts(
    flora: pd.DataFrame,
    residuals: pd.DataFrame,
    *,
    stratum: str,
) -> pd.DataFrame:
    occurrence = flora[["island_id", "accepted_species"]].drop_duplicates().copy()
    occurrence["island_id"] = occurrence["island_id"].astype(str)
    occurrence["accepted_species"] = occurrence["accepted_species"].astype(str)
    rows: list[pd.DataFrame] = []
    for outcome, species in residuals.groupby("outcome", sort=False):
        states = species[["accepted_species", "observed_score"]].drop_duplicates(
            "accepted_species"
        )
        joined = occurrence.merge(
            states, on="accepted_species", how="inner", validate="many_to_one"
        )
        if joined.empty:
            continue
        summary = (
            joined.groupby("island_id", as_index=False)
            .agg(successes=("observed_score", "sum"), trials=("observed_score", "size"))
        )
        summary["outcome"] = str(outcome)
        summary["stratum"] = str(stratum)
        rows.append(summary)
    return (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(
            columns=["island_id", "successes", "trials", "outcome", "stratum"]
        )
    )


def _prepare_stage_data(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    distance = str(config["distance_column"])
    controls = [str(x) for x in config["controls"]]
    required = ["island_id", context, cluster, distance, *controls]
    if missing := set(required) - set(covariates.columns):
        raise ValueError(f"covariates missing columns: {sorted(missing)}")
    data = island_scores.merge(
        covariates[required].drop_duplicates("island_id"),
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
    return data.loc[
        data[context].isin([str(x) for x in config["contexts"]])
        & data[cluster].ne("")
    ].copy()


def _fit_context_stage(
    data: pd.DataFrame,
    config: dict[str, Any],
    *,
    context_value: str,
    stage: str,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, float]]:
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    distance = str(config["distance_column"])
    controls = [str(x) for x in config["controls"]]
    outcomes = [str(x) for x in config["atomic_outcomes"]]
    threshold = int(config["island_model"]["minimum_islands_per_context_outcome"])
    minimum_outcomes = int(config["island_model"]["minimum_outcomes_per_vector"])

    work = data.loc[data[context].eq(context_value)].copy()
    counts = work.groupby("outcome")["island_id"].nunique()
    retained = [x for x in outcomes if int(counts.get(x, 0)) >= threshold]
    if len(retained) < minimum_outcomes:
        return (
            pd.DataFrame(),
            {
                "context": context_value,
                "stage": stage,
                "status": "not_testable",
                "n_retained_outcomes": len(retained),
                "retained_outcomes": "|".join(retained),
            },
            {},
        )

    fits: list[dict[str, Any]] = []
    slope_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in retained:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        for predictor in controls:
            columns.append(_z(part[predictor]))
            names.append(f"{outcome}:z_{predictor}")
        columns.append(_z(part[distance]))
        slope_name = f"{outcome}:z_{distance}"
        names.append(slope_name)
        fit = _fit_ols_component(
            part[stage].to_numpy(float),
            np.column_stack(columns),
            names,
            part[cluster].astype(str).to_numpy(),
        )
        fits.append(fit)
        slope_indices.append(offset + names.index(slope_name))
        rows.append(
            {
                "context": context_value,
                "stage": stage,
                "outcome": outcome,
                "n_islands": int(part["island_id"].nunique()),
                "n_species_median": float(part["n_species"].median()),
            }
        )
        offset += len(names)

    theta, covariance, _ = _assemble_cluster_covariance(fits)
    slopes = theta[slope_indices]
    slope_cov = covariance[np.ix_(slope_indices, slope_indices)]
    ses = np.sqrt(np.clip(np.diag(slope_cov), 0.0, None))
    for row, estimate, se in zip(rows, slopes, ses, strict=True):
        z = float(estimate / se) if se > 0 else float("nan")
        row.update(
            distance_estimate=float(estimate),
            cluster_robust_se=float(se),
            p_value=(math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")),
        )

    rank = int(np.linalg.matrix_rank(slope_cov))
    statistic = (
        float(slopes @ np.linalg.pinv(slope_cov) @ slopes)
        if rank > 0
        else float("nan")
    )
    omnibus = {
        "context": context_value,
        "stage": stage,
        "status": "fit",
        "n_retained_outcomes": len(retained),
        "retained_outcomes": "|".join(retained),
        "n_unique_islands": int(
            work.loc[work["outcome"].isin(retained), "island_id"].nunique()
        ),
        "n_clusters": int(
            work.loc[work["outcome"].isin(retained), cluster].nunique()
        ),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": h1._chi_square_sf_integer_df(statistic, rank),
        "vector_norm": float(np.linalg.norm(slopes)),
    }
    vector = {outcome: float(value) for outcome, value in zip(retained, slopes, strict=True)}
    return pd.DataFrame(rows), omnibus, vector


def _run_stage_models(
    data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[tuple[str, str], dict[str, float]]]:
    slopes: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    vectors: dict[tuple[str, str], dict[str, float]] = {}
    for stage in STAGES:
        stage_rows: list[dict[str, Any]] = []
        for context_value in [str(x) for x in config["contexts"]]:
            s, o, vector = _fit_context_stage(
                data, config, context_value=context_value, stage=stage
            )
            if not s.empty:
                slopes.append(s)
            stage_rows.append(o)
            vectors[(context_value, stage)] = vector
        stage_frame = pd.DataFrame(stage_rows)
        if "p_value" in stage_frame:
            stage_frame["q_value"] = h1._bh(stage_frame["p_value"])
            stage_frame["vector_supported"] = (
                stage_frame["q_value"]
                .le(float(config["multiplicity"]["alpha"]))
                .fillna(False)
            )
        omnibus_rows.extend(stage_frame.to_dict(orient="records"))
    return (
        pd.concat(slopes, ignore_index=True) if slopes else pd.DataFrame(),
        pd.DataFrame(omnibus_rows),
        vectors,
    )


def _geometry(
    vectors: dict[tuple[str, str], dict[str, float]],
    omnibus: pd.DataFrame,
    gate: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    gate_idx = gate.set_index("context") if not gate.empty else pd.DataFrame()
    omni_idx = (
        omnibus.set_index(["context", "stage"]) if not omnibus.empty else pd.DataFrame()
    )
    for context in [str(x) for x in config["contexts"]]:
        obs = vectors.get((context, "observed_score"), {})
        fam = vectors.get((context, "after_family_residual"), {})
        gen = vectors.get((context, "after_genus_residual"), {})
        common = [
            x for x in config["atomic_outcomes"]
            if x in obs and x in fam and x in gen
        ]
        if not common:
            rows.append({"context": context, "status": "not_testable"})
            continue
        a = np.array([obs[x] for x in common], dtype=float)
        b = np.array([fam[x] for x in common], dtype=float)
        c = np.array([gen[x] for x in common], dtype=float)
        an = float(np.linalg.norm(a))
        bn = float(np.linalg.norm(b))
        cn = float(np.linalg.norm(c))

        def cosine(x: np.ndarray, y: np.ndarray) -> float:
            den = float(np.linalg.norm(x) * np.linalg.norm(y))
            return float(x @ y / den) if den > 0 else float("nan")

        def projection(x: np.ndarray, reference: np.ndarray) -> float:
            den = float(reference @ reference)
            return float(x @ reference / den) if den > 0 else float("nan")

        gate_q = (
            float(gate_idx.loc[context, "q_value"])
            if not gate_idx.empty and context in gate_idx.index
            else float("nan")
        )
        genus_q = (
            float(omni_idx.loc[(context, "after_genus_residual"), "q_value"])
            if not omni_idx.empty and (context, "after_genus_residual") in omni_idx.index
            else float("nan")
        )
        family_q = (
            float(omni_idx.loc[(context, "after_family_residual"), "q_value"])
            if not omni_idx.empty and (context, "after_family_residual") in omni_idx.index
            else float("nan")
        )
        genus_cos = cosine(c, a)
        alpha = float(config["multiplicity"]["alpha"])
        if not math.isfinite(gate_q) or gate_q > alpha:
            classification = "current_H1_common_support_not_reproduced"
        elif math.isfinite(genus_q) and genus_q <= alpha:
            classification = (
                "below_genus_or_non_taxonomic_component_retained"
                if math.isfinite(genus_cos) and genus_cos > 0
                else "post_genus_association_retained_but_reoriented"
            )
        elif math.isfinite(family_q) and family_q <= alpha:
            classification = "compatible_with_genus_structuring"
        else:
            classification = "compatible_with_family_or_deeper_structuring"
        rows.append(
            {
                "context": context,
                "status": "fit",
                "n_common_outcomes": len(common),
                "common_outcomes": "|".join(common),
                "common_support_H1_q": gate_q,
                "observed_norm": an,
                "family_norm": bn,
                "genus_norm": cn,
                "family_attenuation": (1.0 - bn / an if an > 0 else float("nan")),
                "genus_attenuation": (1.0 - cn / an if an > 0 else float("nan")),
                "family_cosine_to_observed": cosine(b, a),
                "genus_cosine_to_observed": genus_cos,
                "family_projection_fraction": projection(b, a),
                "genus_projection_fraction": projection(c, a),
                "family_stage_q": family_q,
                "genus_stage_q": genus_q,
                "classification": classification,
            }
        )
    return pd.DataFrame(rows)


def _analyse_scope(
    flora: pd.DataFrame,
    residuals: pd.DataFrame,
    covariates: pd.DataFrame,
    h1_config: dict[str, Any],
    depth_config: dict[str, Any],
    *,
    flora_scope: str,
    evidence_scope: str,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    counts = _common_counts(flora, residuals, stratum=flora_scope)
    gate_config = deepcopy(h1_config)
    gate_config["strata"] = [flora_scope]
    gate_config["between_contexts"] = []
    gate_config["model_outcomes"] = [str(x) for x in depth_config["atomic_outcomes"]]
    gate_slopes, _, gate, _ = h1.run_probability_analysis(
        counts, covariates, gate_config
    )
    island_scores = build_island_stage_scores(flora, residuals)
    data = _prepare_stage_data(island_scores, covariates, depth_config)
    slopes, omnibus, vectors = _run_stage_models(data, depth_config)
    geometry = _geometry(vectors, omnibus, gate, depth_config)

    counts.to_csv(output_dir / "common_support_counts.csv.gz", index=False, compression="gzip")
    gate_slopes.to_csv(output_dir / "common_support_h1_within_slopes.csv", index=False)
    gate.to_csv(output_dir / "common_support_h1_within_omnibus.csv", index=False)
    island_scores.to_csv(
        output_dir / "taxonomic_island_stage_scores.csv.gz",
        index=False,
        compression="gzip",
    )
    slopes.to_csv(output_dir / "taxonomic_stage_slopes.csv", index=False)
    omnibus.to_csv(output_dir / "taxonomic_stage_omnibus.csv", index=False)
    geometry.to_csv(output_dir / "taxonomic_depth_geometry.csv", index=False)

    return {
        "flora_scope": flora_scope,
        "evidence_scope": evidence_scope,
        "n_flora_rows": int(len(flora)),
        "n_flora_islands": int(flora["island_id"].nunique()),
        "n_residual_species": int(residuals["accepted_species"].nunique()),
        "n_residual_species_outcome_rows": int(len(residuals)),
        "common_support_H1_supported_contexts": int(
            gate.get("vector_supported", pd.Series(dtype=bool)).fillna(False).sum()
        ),
        "post_genus_positive_orientation_supported_contexts": int(
            geometry["classification"]
            .eq("below_genus_or_non_taxonomic_component_retained")
            .sum()
        ),
        "geometry": geometry.to_dict(orient="records"),
    }


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--status-flora-csv", type=Path, required=True)
    p.add_argument("--state-audit-csv", type=Path, required=True)
    p.add_argument("--taxonomy-csv", type=Path, required=True)
    p.add_argument("--covariates-csv", type=Path, required=True)
    p.add_argument("--wcvp-ranges-csv", type=Path, required=True)
    p.add_argument("--island-tdwg-csv", type=Path, required=True)
    p.add_argument("--h1-config-path", type=Path, required=True)
    p.add_argument("--depth-config-path", type=Path, required=True)
    p.add_argument("--evidence-scope", required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = _parse_args()
    h1_config = _read_yaml(args.h1_config_path)
    depth_config = _read_yaml(args.depth_config_path)
    if depth_config.get("contract") != "chapter1_current_h1_taxonomic_depth_v1":
        raise ValueError("unexpected taxonomic-depth contract")

    status = pd.read_csv(args.status_flora_csv)
    state_audit = pd.read_csv(args.state_audit_csv)
    taxonomy = pd.read_csv(args.taxonomy_csv)
    covariates = pd.read_csv(args.covariates_csv)
    residuals = build_atomic_taxonomic_residuals(
        state_audit, taxonomy, h1_config, depth_config
    )
    if residuals.empty:
        raise ValueError("no taxonomically eligible species residuals")

    upgraded, wcvp_audit = build_regional_native_compatible_flora(
        status,
        pd.read_csv(args.wcvp_ranges_csv),
        pd.read_csv(args.island_tdwg_csv),
    )
    regional_native = upgraded.loc[
        upgraded["origin_status"].astype(str).eq("native")
    ].copy()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    residuals.to_csv(
        args.output_dir / "species_taxonomic_residuals.csv.gz",
        index=False,
        compression="gzip",
    )
    support = (
        residuals.groupby("outcome", as_index=False)
        .agg(
            n_species=("accepted_species", "nunique"),
            n_families=("family", "nunique"),
            n_genera=("genus", "nunique"),
            median_family_size=("family_n", "median"),
            median_genus_size=("genus_n", "median"),
        )
    )
    support.to_csv(args.output_dir / "taxonomic_species_support.csv", index=False)

    summaries = []
    summaries.append(
        _analyse_scope(
            status,
            residuals,
            covariates,
            h1_config,
            depth_config,
            flora_scope="all_observed",
            evidence_scope=args.evidence_scope,
            output_dir=args.output_dir / "all_observed",
        )
    )
    summaries.append(
        _analyse_scope(
            regional_native,
            residuals,
            covariates,
            h1_config,
            depth_config,
            flora_scope="regional_native_compatible",
            evidence_scope=args.evidence_scope,
            output_dir=args.output_dir / "regional_native_compatible",
        )
    )

    manifest = {
        "contract": depth_config["contract"],
        "evidence_scope": args.evidence_scope,
        "inferential_role": "taxonomic_representation_depth_of_current_seven_response_H1",
        "stages": list(STAGES),
        "wcvp_coverage_audit": wcvp_audit,
        "summaries": summaries,
        "claim_boundary": (
            "Residual support below genus is compatible with within-genus species sorting, "
            "within-lineage change, or non-taxonomic covariation; disappearance after genus "
            "residualization is compatible with genus structuring but does not identify a "
            "historical assembly mechanism."
        ),
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
