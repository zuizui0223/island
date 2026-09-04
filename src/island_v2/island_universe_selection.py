"""Audit and model selection from the physical-island universe into Chapter 1.

This module is outcome-blind. It treats absence from the observed flora/trait
table as missing observation support, never as biological absence or a zero
trait value.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REGION_MAPPING = {
    "northern_midlatitude": "Northern",
    "northern_high_latitude": "Northern high latitude",
    "tropical": "Tropical",
    "southern_extratropical": "Southern",
}

STAGES = (
    "candidate_gshhg",
    "raw_exact_gbif_record",
    "accepted_species_record",
    "strict_species",
    "any_resolved_configured_trait",
    "analysis_after_covariates",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def _island_ids(table: pd.DataFrame) -> set[str]:
    return set(table["island_id"].dropna().astype(str).str.strip()).difference({""})


def _assert_subset(child: set[str], parent: set[str], label: str) -> None:
    outside = sorted(child.difference(parent))
    if outside:
        raise ValueError(f"{label} is not nested in its parent set: {outside[:5]}")


def _fit_ridge_logistic(
    design: np.ndarray,
    outcome: np.ndarray,
    ridge: float,
    max_iterations: int = 100,
    tolerance: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray, int, bool]:
    """Fit deterministic ridge logistic regression with an unpenalized intercept."""
    if ridge < 0:
        raise ValueError("ridge must be non-negative")
    beta = np.zeros(design.shape[1], dtype=float)
    penalty = np.eye(design.shape[1], dtype=float) * ridge
    penalty[0, 0] = 0.0
    converged = False
    for iteration in range(1, max_iterations + 1):
        eta = np.clip(design @ beta, -30.0, 30.0)
        probability = 1.0 / (1.0 + np.exp(-eta))
        weight = np.clip(probability * (1.0 - probability), 1e-9, None)
        gradient = design.T @ (outcome - probability) - penalty @ beta
        information = design.T @ (design * weight[:, None]) + penalty
        step = np.linalg.solve(information, gradient)
        beta += step
        if float(np.max(np.abs(step))) < tolerance:
            converged = True
            break
    probability = 1.0 / (1.0 + np.exp(-np.clip(design @ beta, -30.0, 30.0)))
    return beta, probability, iteration, converged


def _selection_design(
    candidates: pd.DataFrame,
    climate_columns: list[str],
) -> tuple[np.ndarray, list[str]]:
    area = np.log1p(pd.to_numeric(candidates["area_km2"], errors="coerce").clip(lower=0))
    distance = np.log1p(
        pd.to_numeric(candidates["distance_to_continent_km"], errors="coerce").clip(lower=0)
    )
    numeric = pd.DataFrame({"log_area": area, "log_distance": distance})
    for column in climate_columns:
        values = pd.to_numeric(candidates[column], errors="coerce")
        numeric[column] = values.fillna(values.median())
        if values.isna().any():
            numeric[f"{column}_missing"] = values.isna().astype(float)
    for column in ["log_area", "log_distance", *climate_columns]:
        sd = float(numeric[column].std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            numeric[column] = 0.0
        else:
            numeric[column] = (numeric[column] - numeric[column].mean()) / sd

    region = pd.get_dummies(
        candidates["analysis_region_band"], prefix="region", drop_first=True, dtype=float
    )
    frame = pd.concat([numeric.reset_index(drop=True), region.reset_index(drop=True)], axis=1)
    design = np.column_stack([np.ones(len(frame)), frame.to_numpy(dtype=float)])
    return design, ["intercept", *frame.columns.astype(str).tolist()]


def audit_island_universe(
    candidates: pd.DataFrame,
    effort: pd.DataFrame,
    occurrences: pd.DataFrame,
    strict_species: pd.DataFrame,
    composition: pd.DataFrame,
    *,
    ridge: float = 1.0,
    ipw_cap: float = 10.0,
    overlap_lower: float = 0.1,
    overlap_upper: float = 0.9,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return island selection rows, stage/region/quartile summaries, and model audit."""
    _require_columns(
        candidates,
        {
            "island_id",
            "analysis_regime",
            "area_km2",
            "distance_to_continent_km",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        },
        "candidate covariates",
    )
    _require_columns(effort, {"island_id"}, "observation effort")
    _require_columns(occurrences, {"island_id", "species"}, "species occurrences")
    _require_columns(strict_species, {"accepted_species"}, "strict species")
    _require_columns(composition, {"island_id"}, "trait composition")
    if candidates["island_id"].astype(str).duplicated().any():
        raise ValueError("candidate covariates must contain one row per island_id")
    if not 0 < overlap_lower < overlap_upper < 1:
        raise ValueError("overlap bounds must satisfy 0 < lower < upper < 1")
    if ipw_cap <= 0:
        raise ValueError("ipw_cap must be positive")

    candidate_ids = _island_ids(candidates)
    effort_ids = _island_ids(effort)
    occurrence_ids = _island_ids(occurrences)
    accepted = set(
        strict_species["accepted_species"].dropna().astype(str).str.strip()
    ).difference({""})
    occurrence_work = occurrences[["island_id", "species"]].copy()
    occurrence_work["island_id"] = occurrence_work["island_id"].astype(str).str.strip()
    occurrence_work["species"] = occurrence_work["species"].astype(str).str.strip()
    strict_ids = set(
        occurrence_work.loc[occurrence_work["species"].isin(accepted), "island_id"]
    ).difference({""})
    composition_ids = _island_ids(composition)

    _assert_subset(effort_ids, candidate_ids, "effort islands")
    _assert_subset(occurrence_ids, effort_ids, "accepted-species islands")
    _assert_subset(strict_ids, occurrence_ids, "strict-species islands")
    _assert_subset(composition_ids, strict_ids, "resolved-trait islands")

    work = candidates.copy()
    work["island_id"] = work["island_id"].astype(str)
    work["analysis_region_band"] = work["analysis_regime"].map(REGION_MAPPING)
    if work["analysis_region_band"].isna().any():
        bad = sorted(work.loc[work["analysis_region_band"].isna(), "analysis_regime"].unique())
        raise ValueError(f"unmapped analysis regimes: {bad}")

    work["selection_stage"] = "no_exact_gbif_record"
    work.loc[work["island_id"].isin(effort_ids), "selection_stage"] = (
        "raw_no_accepted_species"
    )
    work.loc[work["island_id"].isin(occurrence_ids), "selection_stage"] = (
        "accepted_no_strict_species"
    )
    work.loc[work["island_id"].isin(strict_ids), "selection_stage"] = (
        "strict_no_resolved_trait"
    )
    work.loc[work["island_id"].isin(composition_ids), "selection_stage"] = "analysis"
    work["analysis_included"] = work["island_id"].isin(composition_ids)

    stage_sets = [
        candidate_ids,
        effort_ids,
        occurrence_ids,
        strict_ids,
        composition_ids,
        composition_ids.intersection(candidate_ids),
    ]
    stage_rows = []
    for index, (stage, ids) in enumerate(zip(STAGES, stage_sets, strict=True)):
        previous = stage_sets[index - 1] if index else None
        stage_rows.append(
            {
                "stage": stage,
                "n_islands": len(ids),
                "n_lost_from_previous": 0 if previous is None else len(previous.difference(ids)),
                "retained_fraction_of_candidates": len(ids) / len(candidate_ids),
            }
        )
    stage_summary = pd.DataFrame(stage_rows)

    region_summary = (
        work.groupby("analysis_region_band", as_index=False)
        .agg(
            candidate_islands=("island_id", "nunique"),
            analysis_islands=("analysis_included", "sum"),
        )
        .sort_values("analysis_region_band")
    )
    region_summary["retained_fraction"] = (
        region_summary["analysis_islands"] / region_summary["candidate_islands"]
    )

    quartile_tables = []
    for column, label in (
        ("area_km2", "area"),
        ("distance_to_continent_km", "distance"),
    ):
        ranked = pd.to_numeric(work[column], errors="coerce").rank(method="first")
        quartile = pd.qcut(ranked, 4, labels=["Q1", "Q2", "Q3", "Q4"])
        table = (
            work.assign(quartile=quartile)
            .groupby("quartile", observed=True, as_index=False)
            .agg(
                candidate_islands=("island_id", "nunique"),
                analysis_islands=("analysis_included", "sum"),
            )
        )
        table.insert(0, "axis", label)
        table["retained_fraction"] = table["analysis_islands"] / table["candidate_islands"]
        quartile_tables.append(table)
    quartile_summary = pd.concat(quartile_tables, ignore_index=True)

    climate_columns = ["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    design, feature_names = _selection_design(work, climate_columns)
    outcome = work["analysis_included"].astype(float).to_numpy()
    beta, propensity, iterations, converged = _fit_ridge_logistic(design, outcome, ridge)
    work["selection_propensity"] = propensity
    prevalence = float(outcome.mean())
    work["stabilized_ipw"] = np.where(
        work["analysis_included"],
        np.minimum(prevalence / np.clip(propensity, 1e-9, None), ipw_cap),
        np.nan,
    )
    work["selection_overlap"] = (
        work["analysis_included"]
        & work["selection_propensity"].between(overlap_lower, overlap_upper, inclusive="both")
    )

    included = work.loc[work["analysis_included"]]
    manifest = {
        "contract": "chapter1_island_universe_selection_v1",
        "stage_counts": dict(zip(stage_summary["stage"], stage_summary["n_islands"], strict=True)),
        "selection_model": {
            "outcome": "entry_into_any_resolved_configured_trait_analysis",
            "outcome_blind_to_trait_values": True,
            "features": feature_names,
            "coefficients": {name: float(value) for name, value in zip(feature_names, beta, strict=True)},
            "ridge": ridge,
            "iterations": iterations,
            "converged": converged,
            "candidate_prevalence": prevalence,
            "ipw_cap": ipw_cap,
            "overlap_bounds": [overlap_lower, overlap_upper],
            "n_included_in_overlap": int(included["selection_overlap"].sum()),
            "propensity_quantiles_included": {
                str(key): float(value)
                for key, value in included["selection_propensity"]
                .quantile([0.0, 0.05, 0.5, 0.95, 1.0])
                .items()
            },
            "stabilized_ipw_quantiles": {
                str(key): float(value)
                for key, value in included["stabilized_ipw"]
                .quantile([0.0, 0.5, 0.95, 0.99, 1.0])
                .items()
            },
        },
        "claim_ceiling": (
            "observation-selection audit only; no missing island is a biological absence "
            "or a zero trait value"
        ),
    }
    return work, stage_summary, region_summary, quartile_summary, manifest


@app.command("run")
def run(
    candidates_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    effort_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    occurrences_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    strict_species_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    composition_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option()],
    ridge: Annotated[float, typer.Option()] = 1.0,
    ipw_cap: Annotated[float, typer.Option()] = 10.0,
    overlap_lower: Annotated[float, typer.Option()] = 0.1,
    overlap_upper: Annotated[float, typer.Option()] = 0.9,
) -> None:
    """Write the frozen Chapter 1 island-universe selection audit."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise typer.BadParameter(f"output directory must be absent or empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "candidates": candidates_csv,
        "effort": effort_csv,
        "occurrences": occurrences_csv,
        "strict_species": strict_species_csv,
        "composition": composition_csv,
    }
    tables = {name: pd.read_csv(path, low_memory=False) for name, path in paths.items()}
    island, stages, regions, quartiles, manifest = audit_island_universe(
        tables["candidates"],
        tables["effort"],
        tables["occurrences"],
        tables["strict_species"],
        tables["composition"],
        ridge=ridge,
        ipw_cap=ipw_cap,
        overlap_lower=overlap_lower,
        overlap_upper=overlap_upper,
    )
    manifest["inputs"] = {
        name: {"path": str(path.resolve()), "sha256": _sha256(path)}
        for name, path in paths.items()
    }
    island.to_csv(output_dir / "island_universe_selection.csv.gz", index=False)
    stages.to_csv(output_dir / "island_universe_stage_summary.csv", index=False)
    regions.to_csv(output_dir / "island_universe_region_summary.csv", index=False)
    quartiles.to_csv(output_dir / "island_universe_quartile_summary.csv", index=False)
    (output_dir / "island_universe_selection_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest["stage_counts"], indent=2))


if __name__ == "__main__":
    app()
