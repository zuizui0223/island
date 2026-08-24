"""Rank targeted trait acquisition for corroborated endemic flora.

This module is intentionally island-first. It only emits targets when the
corrected endemic status stratum exists on at least the confirmatory number of
islands but direct trait evidence does not yet support that many islands.

Ranking is lexicographic rather than a tuned composite score:
1. how many currently unsupported endemic islands the species can newly support;
2. whether those islands extend the current isolation support envelope;
3. island-occurrence record count as a coarse evidence-recoverability proxy;
4. stable accepted-species name tie break.

A target is a search priority only. It never creates or infers a trait value.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.native_status_trait_analysis import CONTRASTS, direct_binary_species

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def corroborated_endemic_flora(status: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required - set(status.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")
    work = status.copy()
    for column in required:
        work[column] = _text(work[column])
    work = work.loc[
        work["origin_status"].str.lower().eq("native")
        & work["endemic_status"].str.lower().eq("endemic")
    ]
    return work[["island_id", "accepted_species"]].drop_duplicates()


def _outcome_direct_species(evidence: pd.DataFrame, outcome: str) -> set[str]:
    if outcome not in CONTRASTS:
        raise typer.BadParameter(f"unknown outcome: {outcome}")
    mapped = direct_binary_species(evidence, **CONTRASTS[outcome])
    return set(mapped["accepted_species"].astype(str))


def rank_endemic_targets(
    status: pd.DataFrame,
    evidence: pd.DataFrame,
    island_taxa: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    confirmatory_min_islands: int = 50,
    regimes: tuple[str, ...] = ("northern_midlatitude", "tropical", "southern_extratropical"),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    endemic = corroborated_endemic_flora(status)
    required_cov = {"island_id", "analysis_regime", "log_distance_to_continent_km"}
    missing_cov = required_cov - set(covariates.columns)
    if missing_cov:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing_cov)}")
    endemic = endemic.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    endemic = endemic.loc[endemic["analysis_regime"].isin(regimes)].copy()

    required_taxa = {"accepted_species", "n_records", "n_islands"}
    missing_taxa = required_taxa - set(island_taxa.columns)
    if missing_taxa:
        raise typer.BadParameter(f"island taxa missing columns: {sorted(missing_taxa)}")
    taxa = island_taxa[list(required_taxa)].drop_duplicates("accepted_species").copy()
    taxa["accepted_species"] = _text(taxa["accepted_species"])
    taxa["n_records"] = pd.to_numeric(taxa["n_records"], errors="coerce").fillna(0).astype(int)
    taxa["n_islands"] = pd.to_numeric(taxa["n_islands"], errors="coerce").fillna(0).astype(int)

    summary_rows: list[dict[str, Any]] = []
    island_rows: list[dict[str, Any]] = []
    species_rows: list[dict[str, Any]] = []

    for regime in regimes:
        region = endemic.loc[endemic["analysis_regime"].eq(regime)].copy()
        ceiling_islands = set(region["island_id"])
        for outcome in CONTRASTS:
            covered_species = _outcome_direct_species(evidence, outcome)
            supported_islands = set(
                region.loc[region["accepted_species"].isin(covered_species), "island_id"]
            )
            unsupported = ceiling_islands - supported_islands
            ceiling = len(ceiling_islands)
            supported = len(supported_islands)
            gap = max(0, confirmatory_min_islands - supported)
            if ceiling < confirmatory_min_islands:
                decision = "status_limited"
            elif supported >= confirmatory_min_islands:
                decision = "model_ready"
            elif supported >= 30:
                decision = "targeted_trait_acquisition"
            else:
                decision = "recoverability_test_before_acquisition"

            supported_distance = region.loc[
                region["island_id"].isin(supported_islands),
                "log_distance_to_continent_km",
            ].dropna()
            q05 = float(supported_distance.quantile(0.05)) if len(supported_distance) else np.nan
            q95 = float(supported_distance.quantile(0.95)) if len(supported_distance) else np.nan

            summary_rows.append(
                {
                    "regime": regime,
                    "outcome": outcome,
                    "n_endemic_status_islands": ceiling,
                    "n_direct_supported_islands": supported,
                    "gap_to_50": gap,
                    "n_currently_unsupported_islands": len(unsupported),
                    "decision": decision,
                    "supported_isolation_q05": q05,
                    "supported_isolation_q95": q95,
                }
            )
            if decision not in {"targeted_trait_acquisition", "recoverability_test_before_acquisition"}:
                continue

            unsupported_rows = region.loc[region["island_id"].isin(unsupported)].copy()
            for island_id, group in unsupported_rows.groupby("island_id", sort=True):
                distance = pd.to_numeric(group["log_distance_to_continent_km"], errors="coerce").iloc[0]
                edge = bool(
                    np.isfinite(distance)
                    and (
                        (np.isfinite(q05) and distance < q05)
                        or (np.isfinite(q95) and distance > q95)
                    )
                )
                island_rows.append(
                    {
                        "regime": regime,
                        "outcome": outcome,
                        "island_id": island_id,
                        "log_distance_to_continent_km": distance,
                        "extends_supported_isolation_5_95": edge,
                        "n_candidate_endemic_species": int(
                            (~group["accepted_species"].isin(covered_species)).sum()
                        ),
                    }
                )

            candidates = unsupported_rows.loc[
                ~unsupported_rows["accepted_species"].isin(covered_species)
            ].copy()
            if candidates.empty:
                continue
            candidate_summary = (
                candidates.groupby("accepted_species", as_index=False)
                .agg(
                    unsupported_islands_reached=("island_id", "nunique"),
                    isolation_min=("log_distance_to_continent_km", "min"),
                    isolation_max=("log_distance_to_continent_km", "max"),
                )
            )
            edge_species = (
                candidates.assign(
                    edge=lambda x: (
                        (np.isfinite(q05) & x["log_distance_to_continent_km"].lt(q05))
                        | (np.isfinite(q95) & x["log_distance_to_continent_km"].gt(q95))
                    )
                )
                .groupby("accepted_species")["edge"]
                .any()
                .rename("extends_supported_isolation_5_95")
            )
            candidate_summary = candidate_summary.merge(
                edge_species, on="accepted_species", how="left", validate="one_to_one"
            ).merge(taxa, on="accepted_species", how="left", validate="one_to_one")
            candidate_summary["n_records"] = candidate_summary["n_records"].fillna(0).astype(int)
            candidate_summary["n_islands"] = candidate_summary["n_islands"].fillna(0).astype(int)
            candidate_summary = candidate_summary.sort_values(
                [
                    "unsupported_islands_reached",
                    "extends_supported_isolation_5_95",
                    "n_records",
                    "accepted_species",
                ],
                ascending=[False, False, False, True],
            ).reset_index(drop=True)
            candidate_summary["priority_rank"] = np.arange(1, len(candidate_summary) + 1)
            candidate_summary.insert(0, "outcome", outcome)
            candidate_summary.insert(0, "regime", regime)
            candidate_summary["gap_to_50"] = gap
            species_rows.extend(candidate_summary.to_dict("records"))

    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(island_rows),
        pd.DataFrame(species_rows),
    )


@app.command("run")
def run(
    status_ledger_csv: Path = typer.Option(..., exists=True),
    direct_evidence_csv: Path = typer.Option(..., exists=True),
    island_taxa_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    confirmatory_min_islands: int = typer.Option(50),
) -> None:
    status = pd.read_csv(status_ledger_csv)
    evidence = pd.read_csv(direct_evidence_csv)
    taxa = pd.read_csv(island_taxa_csv)
    covariates = pd.read_csv(covariates_csv)
    summary, islands, species = rank_endemic_targets(
        status,
        evidence,
        taxa,
        covariates,
        confirmatory_min_islands=confirmatory_min_islands,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "endemic_trait_support_decisions.csv", index=False)
    islands.to_csv(output_dir / "endemic_unsupported_islands.csv", index=False)
    species.to_csv(output_dir / "endemic_trait_acquisition_targets.csv", index=False)
    manifest = {
        "contract": "corroborated_endemic_trait_targeting_v1",
        "confirmatory_min_islands": confirmatory_min_islands,
        "n_decision_rows": int(len(summary)),
        "n_target_species_rows": int(len(species)),
        "policy": (
            "Targets are emitted only for trait-limited endemic strata. Ranking is "
            "island gain, isolation-edge gain, then GBIF island-record count."
        ),
        "decisions": summary.to_dict("records"),
    }
    (output_dir / "endemic_trait_target_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
