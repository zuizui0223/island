"""Replay the frozen tropical native-nonendemic lineage entry/loading bridge on corrected geography.

This is a geography measurement-repair replay, not a new model search. It reuses the
frozen island-level lineage representation scores, source modes, source matching,
floristic stratum, evidence scopes, outcomes, equal-island OLS model, covariates and
spatial-block clustered covariance. Only the geographic distance exposure is replaced
by the corrected 24 September 2026 coastline distance.

The historical functional position is (-large_bee_like + generalized_accessible) / 2.
This diagnostic is therefore not a direct decomposition of the current seven-response H1.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

PREDICTORS = [
    "log_distance_to_continent_km",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
]
OUTCOMES = ["entry_enrichment", "species_enrichment", "loading_increment"]
EVIDENCE_SCOPES = ["broad", "broad_direct"]
SOURCE_MODES = ["geo_k5", "geo_k10", "geo_k20", "geo50_climate10"]


def _p2(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.dropna().astype(float)
    if finite.empty:
        return out
    order = finite.sort_values().index.tolist()
    n = len(order)
    ranked = [float(finite.loc[idx]) * n / rank for rank, idx in enumerate(order, start=1)]
    adjusted = [0.0] * n
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i])
        adjusted[i] = min(1.0, running)
    for idx, value in zip(order, adjusted, strict=True):
        out.loc[idx] = value
    return out


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def fit_equal_island_clustered(
    data: pd.DataFrame,
    *,
    response: str,
) -> dict[str, Any]:
    needed = {response, "spatial_block", *PREDICTORS}
    missing = needed - set(data.columns)
    if missing:
        raise typer.BadParameter(f"model data missing columns: {sorted(missing)}")
    work = data[[response, "spatial_block", *PREDICTORS]].copy()
    for column in [response, *PREDICTORS]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work["spatial_block"] = work["spatial_block"].fillna("").astype(str)
    work = work.dropna().loc[lambda x: x["spatial_block"].ne("")].copy()
    columns = [np.ones(len(work), dtype=float)]
    names = ["intercept"]
    for predictor in PREDICTORS:
        columns.append(_z(work[predictor]))
        names.append(predictor)
    X = np.column_stack(columns)
    y = work[response].to_numpy(float)
    bread = np.linalg.pinv(X.T @ X)
    beta = bread @ (X.T @ y)
    residual = y - X @ beta
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    labels = work["spatial_block"].to_numpy(str)
    for cluster in np.unique(labels):
        mask = labels == cluster
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(work)
    p = X.shape[1]
    g = len(np.unique(labels))
    if g > 1 and n > p:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - p))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    idx = names.index("log_distance_to_continent_km")
    z = float(beta[idx] / se[idx]) if se[idx] > 0 else float("nan")
    return {
        "estimate": float(beta[idx]),
        "se": float(se[idx]),
        "p_value": _p2(z) if math.isfinite(z) else float("nan"),
        "n_islands": int(n),
        "n_clusters": int(g),
    }


def _merge_covariates(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
) -> pd.DataFrame:
    keep = ["island_id", "analysis_regime", "spatial_block", *PREDICTORS]
    missing = set(keep) - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    return scores.merge(
        covariates[keep].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )


def _fit_surface(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    geography_label: str,
) -> pd.DataFrame:
    data = _merge_covariates(scores, covariates)
    rows: list[dict[str, Any]] = []
    for evidence_scope in EVIDENCE_SCOPES:
        for source_mode in SOURCE_MODES:
            subset = data.loc[
                data["evidence_scope"].astype(str).eq(evidence_scope)
                & data["stratum"].astype(str).eq("native_nonendemic")
                & data["source_mode"].astype(str).eq(source_mode)
                & data["source_matching"].astype(str).eq("prevalence_richness")
                & pd.to_numeric(data["minimum_represented_genera"], errors="coerce").eq(5)
                & data["analysis_regime"].astype(str).eq("tropical")
            ].copy()
            for outcome in OUTCOMES:
                fit = fit_equal_island_clustered(subset, response=outcome)
                rows.append(
                    {
                        "geography": geography_label,
                        "evidence_scope": evidence_scope,
                        "source_mode": source_mode,
                        "outcome": outcome,
                        **fit,
                    }
                )
    out = pd.DataFrame(rows)
    out["q_across_source_modes"] = np.nan
    for (evidence_scope, outcome), index in out.groupby(
        ["evidence_scope", "outcome"]
    ).groups.items():
        out.loc[index, "q_across_source_modes"] = _bh(out.loc[index, "p_value"])
    return out


def _validate_original(
    replay: pd.DataFrame,
    frozen: pd.DataFrame,
) -> dict[str, Any]:
    old = frozen.loc[
        frozen["stratum"].astype(str).eq("native_nonendemic")
        & frozen["context_layer"].astype(str).eq("analysis_regime")
        & frozen["context"].astype(str).eq("tropical")
        & frozen["source_matching"].astype(str).eq("prevalence_richness")
        & pd.to_numeric(frozen["minimum_represented_genera"], errors="coerce").eq(5)
        & frozen["evidence_scope"].astype(str).isin(EVIDENCE_SCOPES)
        & frozen["source_mode"].astype(str).isin(SOURCE_MODES)
        & frozen["outcome"].astype(str).isin(OUTCOMES)
    ].copy()
    merged = replay.merge(
        old[
            [
                "evidence_scope",
                "source_mode",
                "outcome",
                "distance_slope",
                "cluster_robust_se",
                "p_value",
                "n_islands",
            ]
        ],
        on=["evidence_scope", "source_mode", "outcome"],
        how="inner",
        validate="one_to_one",
        suffixes=("_replay", "_frozen"),
    )
    if len(merged) != len(EVIDENCE_SCOPES) * len(SOURCE_MODES) * len(OUTCOMES):
        raise AssertionError(f"unexpected frozen replay row count: {len(merged)}")
    checks = {}
    for name, left, right in [
        ("estimate", "estimate", "distance_slope"),
        ("se", "se", "cluster_robust_se"),
        ("p", "p_value_replay", "p_value_frozen"),
    ]:
        diff = np.abs(
            pd.to_numeric(merged[left], errors="raise").to_numpy(float)
            - pd.to_numeric(merged[right], errors="raise").to_numpy(float)
        )
        checks[f"max_abs_{name}_difference"] = float(np.max(diff))
        if not np.allclose(
            pd.to_numeric(merged[left], errors="raise").to_numpy(float),
            pd.to_numeric(merged[right], errors="raise").to_numpy(float),
            atol=1e-10,
            rtol=1e-8,
        ):
            raise AssertionError(f"frozen replay mismatch for {name}: max={np.max(diff)}")
    if not np.array_equal(
        pd.to_numeric(merged["n_islands_replay"], errors="raise").to_numpy(int),
        pd.to_numeric(merged["n_islands_frozen"], errors="raise").to_numpy(int),
    ):
        raise AssertionError("frozen replay island counts differ")
    return {"status": "pass", "n_rows": int(len(merged)), **checks}


def _classify(corrected: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (evidence_scope, outcome), part in corrected.groupby(
        ["evidence_scope", "outcome"], sort=False
    ):
        signs = np.sign(part["estimate"].to_numpy(float))
        same_nonzero_sign = bool(np.all(signs == signs[0]) and signs[0] != 0)
        all_fdr = bool(part["q_across_source_modes"].le(0.05).all())
        rows.append(
            {
                "evidence_scope": evidence_scope,
                "outcome": outcome,
                "all_four_source_modes_same_sign": same_nonzero_sign,
                "direction": (
                    "positive"
                    if same_nonzero_sign and signs[0] > 0
                    else "negative"
                    if same_nonzero_sign and signs[0] < 0
                    else "mixed"
                ),
                "all_four_source_modes_fdr_supported": all_fdr,
                "source_mode_robust": bool(same_nonzero_sign and all_fdr),
                "min_estimate": float(part["estimate"].min()),
                "max_estimate": float(part["estimate"].max()),
                "max_q": float(part["q_across_source_modes"].max()),
            }
        )
    return pd.DataFrame(rows)


@app.command()
def main(
    scores_csv: Path = typer.Option(..., exists=True),
    frozen_slopes_csv: Path = typer.Option(..., exists=True),
    original_covariates_csv: Path = typer.Option(..., exists=True),
    corrected_covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    scores = pd.read_csv(scores_csv)
    frozen = pd.read_csv(frozen_slopes_csv)
    original = pd.read_csv(original_covariates_csv)
    corrected = pd.read_csv(corrected_covariates_csv)

    original_replay = _fit_surface(scores, original, geography_label="original")
    gate = _validate_original(original_replay, frozen)
    corrected_replay = _fit_surface(scores, corrected, geography_label="corrected")
    classification = _classify(corrected_replay)

    output_dir.mkdir(parents=True, exist_ok=True)
    original_replay.to_csv(output_dir / "original_replay.csv", index=False)
    corrected_replay.to_csv(output_dir / "corrected_replay.csv", index=False)
    classification.to_csv(output_dir / "classification.csv", index=False)
    result = {
        "contract": "chapter1_corrected_lineage_entry_loading_replay_v1",
        "inferential_role": "measurement_repair_replay_of_frozen_lineage_representation_bridge",
        "functional_position": "(-large_bee_like + generalized_accessible) / 2",
        "scope": {
            "context": "tropical",
            "stratum": "native_nonendemic",
            "evidence_scopes": EVIDENCE_SCOPES,
            "source_modes": SOURCE_MODES,
            "source_matching": "prevalence_richness",
            "minimum_represented_genera": 5,
        },
        "original_reproduction_gate": gate,
        "claim_boundary": (
            "This is an assembly diagnostic for a historical frozen floral functional "
            "position, not a direct decomposition of the current seven-response H1. "
            "Entry/loading associations do not identify colonisation, extinction, "
            "pollinator loss or within-lineage evolution."
        ),
    }
    (output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(classification.to_csv(index=False))


if __name__ == "__main__":
    app()
