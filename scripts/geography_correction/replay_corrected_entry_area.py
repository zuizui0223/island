"""Replay frozen source-lineage distance × area moderation on corrected geography.

This is a measurement-repair replay of the historical Chapter 1 area-capacity model.
It keeps the frozen tropical/northern source-lineage design and replaces only the
geographic distance exposure with the corrected 24 September 2026 coastline distance.

The source-lineage family uses the historical floral functional-position bridge and is
not a direct decomposition of the current seven-response H1.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import chi2

PREDICTORS = [
    "log_distance_to_continent_km",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
]
CONTROLS = ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
RESPONSES = ["entry_enrichment", "loading_increment"]
CONTEXTS = ["northern_midlatitude", "tropical"]
SOURCE_MODES = ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]


def _p2(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.dropna().astype(float)
    if finite.empty:
        return result
    ordered = finite.sort_values().index.tolist()
    m = len(ordered)
    raw = [float(finite.loc[idx]) * m / rank for rank, idx in enumerate(ordered, start=1)]
    adjusted = [0.0] * m
    running = 1.0
    for i in range(m - 1, -1, -1):
        running = min(running, raw[i])
        adjusted[i] = min(1.0, running)
    for idx, value in zip(ordered, adjusted, strict=True):
        result.loc[idx] = value
    return result


def _standardize(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (values - mean) / sd


def _fit_clustered(
    y: np.ndarray,
    design: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    n, p = design.shape
    unique_clusters = np.unique(clusters)
    if n < max(10, p + 3) or len(unique_clusters) < 2:
        return pd.DataFrame(), np.empty((0, 0)), {
            "status": "insufficient_complete_rows",
            "n_rows": int(n),
            "n_clusters": int(len(unique_clusters)),
        }
    bread = np.linalg.pinv(design.T @ design)
    beta = bread @ (design.T @ y)
    residual = y - design @ beta
    meat = np.zeros((p, p), dtype=float)
    for cluster in unique_clusters:
        mask = clusters == cluster
        score = design[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    g = len(unique_clusters)
    if g > 1 and n > p:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - p))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    rows = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        rows.append(
            {
                "predictor": name,
                "estimate": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": _p2(z) if math.isfinite(z) else float("nan"),
            }
        )
    return pd.DataFrame(rows), covariance, {
        "status": "fit",
        "n_rows": int(n),
        "n_clusters": int(g),
    }


def _lineage_family(scores: pd.DataFrame) -> pd.DataFrame:
    required = {
        "island_id",
        "stratum",
        "source_mode",
        "evidence_scope",
        "source_matching",
        "minimum_represented_genera",
        *RESPONSES,
    }
    missing = required - set(scores.columns)
    if missing:
        raise ValueError(f"lineage scores missing columns: {sorted(missing)}")
    selected = scores.loc[
        scores["evidence_scope"].astype(str).eq("broad")
        & scores["source_matching"].astype(str).eq("prevalence_richness")
        & pd.to_numeric(scores["minimum_represented_genera"], errors="coerce").eq(5)
    ].copy()
    return selected.melt(
        id_vars=["island_id", "stratum", "source_mode"],
        value_vars=RESPONSES,
        var_name="response",
        value_name="response_score",
    )


def _prepare(family: pd.DataFrame, covariates: pd.DataFrame) -> pd.DataFrame:
    needed = ["island_id", "analysis_regime", "spatial_block", *PREDICTORS]
    missing = set(needed) - set(covariates.columns)
    if missing:
        raise ValueError(f"covariates missing columns: {sorted(missing)}")
    data = family.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["response_score", *PREDICTORS]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data["analysis_regime"] = data["analysis_regime"].fillna("").astype(str)
    data["spatial_block"] = data["spatial_block"].fillna("").astype(str)
    return data


def _fit_context(
    data: pd.DataFrame,
    *,
    source_mode: str,
    context: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    needed = ["response_score", "spatial_block", *PREDICTORS]
    work = data.loc[
        data["source_mode"].astype(str).eq(source_mode)
        & data["stratum"].astype(str).eq("native_nonendemic")
        & data["analysis_regime"].astype(str).eq(context)
        & data["response"].astype(str).isin(RESPONSES)
    ].dropna(subset=needed).copy()

    counts = work.groupby("response")["island_id"].nunique()
    if any(int(counts.get(response, 0)) < 50 for response in RESPONSES):
        return pd.DataFrame(), {
            "source_mode": source_mode,
            "context": context,
            "status": "not_testable",
        }

    names: list[str] = []
    columns: list[np.ndarray] = []
    terms: dict[str, dict[str, str]] = {}
    for response in RESPONSES:
        mask = work["response"].astype(str).eq(response).to_numpy()
        indicator = mask.astype(float)
        names.append(f"response[{response}]")
        columns.append(indicator)
        for control in CONTROLS:
            values = np.zeros(len(work), dtype=float)
            values[mask] = _standardize(work.loc[mask, control])
            names.append(f"response[{response}]:z_{control}")
            columns.append(values)
        z_distance = np.zeros(len(work), dtype=float)
        z_area = np.zeros(len(work), dtype=float)
        z_distance[mask] = _standardize(work.loc[mask, "log_distance_to_continent_km"])
        z_area[mask] = _standardize(work.loc[mask, "log_island_area_km2"])
        d_name = f"response[{response}]:z_distance"
        a_name = f"response[{response}]:z_area"
        i_name = f"response[{response}]:z_distance:z_area"
        names.extend([d_name, a_name, i_name])
        columns.extend([z_distance, z_area, z_distance * z_area])
        terms[response] = {"distance": d_name, "area": a_name, "interaction": i_name}

    coefficients, covariance, fit = _fit_clustered(
        work["response_score"].to_numpy(float),
        np.column_stack(columns),
        names,
        work["spatial_block"].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {
            "source_mode": source_mode,
            "context": context,
            "status": fit["status"],
        }

    indexed = coefficients.set_index("predictor")
    interaction_names = [terms[x]["interaction"] for x in RESPONSES]
    indices = [names.index(x) for x in interaction_names]
    vector = np.array([float(indexed.loc[x, "estimate"]) for x in interaction_names])
    cov = covariance[np.ix_(indices, indices)]
    rank = int(np.linalg.matrix_rank(cov))
    statistic = float(vector @ np.linalg.pinv(cov) @ vector)
    p_value = float(chi2.sf(statistic, rank))

    rows: list[dict[str, Any]] = []
    for response in RESPONSES:
        d = indexed.loc[terms[response]["distance"]]
        i = indexed.loc[terms[response]["interaction"]]
        a = indexed.loc[terms[response]["area"]]
        d_est = float(d["estimate"])
        i_est = float(i["estimate"])
        rows.append(
            {
                "source_mode": source_mode,
                "context": context,
                "response": response,
                "n_islands": int(counts.loc[response]),
                "distance_estimate_at_mean_area": d_est,
                "distance_se": float(d["cluster_robust_se"]),
                "distance_p": float(d["p_value"]),
                "area_estimate_at_mean_distance": float(a["estimate"]),
                "area_se": float(a["cluster_robust_se"]),
                "area_p": float(a["p_value"]),
                "distance_x_area_estimate": i_est,
                "distance_x_area_se": float(i["cluster_robust_se"]),
                "distance_x_area_p": float(i["p_value"]),
                "distance_slope_at_small_area_z_minus1": d_est - i_est,
                "distance_slope_at_large_area_z_plus1": d_est + i_est,
            }
        )
    return pd.DataFrame(rows), {
        "source_mode": source_mode,
        "context": context,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_distance_x_area_chisq": statistic,
        "joint_distance_x_area_df": rank,
        "p_value": p_value,
    }


def _run(scores: pd.DataFrame, covariates: pd.DataFrame, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    family = _lineage_family(scores)
    data = _prepare(family, covariates)
    coeffs: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    for source_mode in SOURCE_MODES:
        for context in CONTEXTS:
            frame, result = _fit_context(data, source_mode=source_mode, context=context)
            if not frame.empty:
                coeffs.append(frame)
            omnibus_rows.append(result)
    coefficients = pd.concat(coeffs, ignore_index=True)
    omnibus = pd.DataFrame(omnibus_rows)

    fit = omnibus["status"].eq("fit")
    omnibus["q_vector_family"] = np.nan
    for source_mode, idx in omnibus.loc[fit].groupby("source_mode").groups.items():
        omnibus.loc[idx, "q_vector_family"] = _bh(omnibus.loc[idx, "p_value"])
    omnibus["area_moderation_vector_supported"] = omnibus["q_vector_family"].le(0.05).fillna(False)

    coefficients["distance_q"] = np.nan
    coefficients["distance_x_area_q"] = np.nan
    for source_mode, idx in coefficients.groupby("source_mode").groups.items():
        coefficients.loc[idx, "distance_q"] = _bh(coefficients.loc[idx, "distance_p"])
        coefficients.loc[idx, "distance_x_area_q"] = _bh(
            coefficients.loc[idx, "distance_x_area_p"]
        )
    gate = omnibus[
        ["source_mode", "context", "area_moderation_vector_supported"]
    ]
    coefficients = coefficients.merge(
        gate,
        on=["source_mode", "context"],
        how="left",
        validate="many_to_one",
    )
    coefficients["distance_axis_supported"] = coefficients["distance_q"].le(0.05)
    coefficients["interaction_axis_supported"] = coefficients["distance_x_area_q"].le(0.05)
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

    coefficients.insert(0, "geography", label)
    omnibus.insert(0, "geography", label)
    return coefficients, omnibus


def _validate_original(
    replay_coeff: pd.DataFrame,
    replay_omnibus: pd.DataFrame,
    frozen_coeff: pd.DataFrame,
    frozen_omnibus: pd.DataFrame,
) -> dict[str, Any]:
    fc = frozen_coeff.loc[
        frozen_coeff["family"].astype(str).eq("source_lineage_broad")
        & frozen_coeff["context_layer"].astype(str).eq("analysis_regime")
        & frozen_coeff["stratum"].astype(str).eq("native_nonendemic")
        & frozen_coeff["support_tier"].astype(str).eq("confirmatory")
        & frozen_coeff["source_mode"].astype(str).isin(SOURCE_MODES)
        & frozen_coeff["context"].astype(str).isin(CONTEXTS)
        & frozen_coeff["response"].astype(str).isin(RESPONSES)
    ].copy()
    fo = frozen_omnibus.loc[
        frozen_omnibus["family"].astype(str).eq("source_lineage_broad")
        & frozen_omnibus["context_layer"].astype(str).eq("analysis_regime")
        & frozen_omnibus["stratum"].astype(str).eq("native_nonendemic")
        & frozen_omnibus["support_tier"].astype(str).eq("confirmatory")
        & frozen_omnibus["source_mode"].astype(str).isin(SOURCE_MODES)
        & frozen_omnibus["context"].astype(str).isin(CONTEXTS)
    ].copy()

    cm = replay_coeff.merge(
        fc,
        on=["source_mode", "context", "response"],
        validate="one_to_one",
        suffixes=("_replay", "_frozen"),
    )
    om = replay_omnibus.merge(
        fo,
        on=["source_mode", "context"],
        validate="one_to_one",
        suffixes=("_replay", "_frozen"),
    )
    if len(cm) != 16 or len(om) != 8:
        raise AssertionError(f"unexpected replay support: coeff={len(cm)}, omnibus={len(om)}")

    checks: dict[str, float] = {}
    pairs = [
        ("distance", "distance_estimate_at_mean_area_replay", "distance_estimate_at_mean_area_frozen"),
        ("distance_se", "distance_se_replay", "distance_se_frozen"),
        ("interaction", "distance_x_area_estimate_replay", "distance_x_area_estimate_frozen"),
        ("interaction_se", "distance_x_area_se_replay", "distance_x_area_se_frozen"),
        ("distance_q", "distance_q_replay", "distance_q_frozen"),
        ("interaction_q", "distance_x_area_q_replay", "distance_x_area_q_frozen"),
    ]
    for name, left, right in pairs:
        a = pd.to_numeric(cm[left], errors="raise").to_numpy(float)
        b = pd.to_numeric(cm[right], errors="raise").to_numpy(float)
        checks[f"max_abs_{name}_difference"] = float(np.max(np.abs(a - b)))
        if not np.allclose(a, b, atol=1e-10, rtol=1e-8):
            raise AssertionError(f"coefficient replay mismatch: {name}")
    a = pd.to_numeric(om["q_vector_family_replay"], errors="raise").to_numpy(float)
    b = pd.to_numeric(om["q_vector_family_frozen"], errors="raise").to_numpy(float)
    checks["max_abs_vector_q_difference"] = float(np.max(np.abs(a - b)))
    if not np.allclose(a, b, atol=1e-10, rtol=1e-8):
        raise AssertionError("vector q replay mismatch")
    return {"status": "pass", "n_coefficient_rows": 16, "n_omnibus_rows": 8, **checks}


def _tropical_summary(coefficients: pd.DataFrame, omnibus: pd.DataFrame) -> pd.DataFrame:
    tropical = coefficients.loc[coefficients["context"].eq("tropical")].copy()
    vector = omnibus.loc[omnibus["context"].eq("tropical")].set_index("source_mode")
    rows: list[dict[str, Any]] = []
    for source_mode in SOURCE_MODES:
        entry = tropical.loc[
            tropical["source_mode"].eq(source_mode)
            & tropical["response"].eq("entry_enrichment")
        ].iloc[0]
        rows.append(
            {
                "source_mode": source_mode,
                "vector_q": float(vector.loc[source_mode, "q_vector_family"]),
                "vector_supported": bool(
                    vector.loc[source_mode, "area_moderation_vector_supported"]
                ),
                "entry_distance": float(entry["distance_estimate_at_mean_area"]),
                "entry_distance_q": float(entry["distance_q"]),
                "entry_distance_x_area": float(entry["distance_x_area_estimate"]),
                "entry_interaction_q": float(entry["distance_x_area_q"]),
                "entry_small_island_slope": float(
                    entry["distance_slope_at_small_area_z_minus1"]
                ),
                "entry_large_island_slope": float(
                    entry["distance_slope_at_large_area_z_plus1"]
                ),
                "formal_state": str(entry["area_moderation_state"]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage-scores-csv", type=Path, required=True)
    parser.add_argument("--frozen-coefficients-csv", type=Path, required=True)
    parser.add_argument("--frozen-omnibus-csv", type=Path, required=True)
    parser.add_argument("--original-covariates-csv", type=Path, required=True)
    parser.add_argument("--corrected-covariates-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    scores = pd.read_csv(args.lineage_scores_csv)
    frozen_coeff = pd.read_csv(args.frozen_coefficients_csv)
    frozen_omnibus = pd.read_csv(args.frozen_omnibus_csv)
    original_cov = pd.read_csv(args.original_covariates_csv)
    corrected_cov = pd.read_csv(args.corrected_covariates_csv)

    original_coeff, original_omnibus = _run(scores, original_cov, "original")
    gate = _validate_original(
        original_coeff, original_omnibus, frozen_coeff, frozen_omnibus
    )
    corrected_coeff, corrected_omnibus = _run(scores, corrected_cov, "corrected")
    summary = _tropical_summary(corrected_coeff, corrected_omnibus)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    original_coeff.to_csv(args.output_dir / "original_coefficients.csv", index=False)
    original_omnibus.to_csv(args.output_dir / "original_omnibus.csv", index=False)
    corrected_coeff.to_csv(args.output_dir / "corrected_coefficients.csv", index=False)
    corrected_omnibus.to_csv(args.output_dir / "corrected_omnibus.csv", index=False)
    summary.to_csv(args.output_dir / "tropical_entry_area_summary.csv", index=False)
    result = {
        "contract": "chapter1_corrected_entry_area_replay_v1",
        "inferential_role": "measurement_repair_replay_of_frozen_area_capacity_model",
        "original_reproduction_gate": gate,
        "scope": {
            "family": "source_lineage_broad",
            "stratum": "native_nonendemic",
            "evidence_scope": "broad",
            "source_matching": "prevalence_richness",
            "minimum_represented_genera": 5,
            "contexts": CONTEXTS,
            "source_modes": SOURCE_MODES,
            "responses": RESPONSES,
        },
        "claim_boundary": (
            "Area is a composite capacity proxy. This replay does not identify target "
            "size, habitat diversity, persistence, colonisation probability, extinction, "
            "pollinator service, or within-lineage evolution."
        ),
    }
    (args.output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_csv(index=False))


if __name__ == "__main__":
    main()
