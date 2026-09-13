"""Effort-matched source-side GloBI breadth analysis.

V2 leaves the outcome-blind genus predictor unchanged and improves only the expected
source comparator by matching genera on source prevalence, source richness class and
GloBI independent-reference effort class before island outcomes are compared.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_globi_source_breadth import (
    _all_effort_thresholds,
    fit_breadth_models,
    sha256_file,
)
from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _availability_matrices,
    _representation_count_matrix,
    _richness_bins,
    _source_assignment_matrix,
    broad_source_availability,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != "chapter1_globi_source_breadth_v2":
        raise typer.BadParameter("unexpected GloBI source-breadth v2 contract")
    return config


def reference_effort_bins(n_references: np.ndarray) -> np.ndarray:
    values = np.asarray(n_references, dtype=float)
    bins = np.full(len(values), -1, dtype=np.int16)
    ok = np.isfinite(values) & (values >= 1)
    bins[ok] = np.floor(np.log2(values[ok])).astype(np.int16)
    return bins


def compute_effort_matched_enrichment(
    prevalence: np.ndarray,
    source_richness: np.ndarray,
    island_species_counts: np.ndarray,
    positions: np.ndarray,
    effort_bins: np.ndarray,
    *,
    matching: str,
    minimum_represented_genera: int,
) -> dict[str, float | int] | None:
    candidate = (prevalence > 0) & np.isfinite(positions) & (effort_bins >= 0)
    counts = island_species_counts.copy()
    counts[~candidate] = 0
    represented = counts > 0
    n_represented = int(represented.sum())
    n_species = int(counts.sum())
    if n_represented < int(minimum_represented_genera) or n_species <= 0:
        return None

    observed_entry = float(np.mean(positions[represented]))
    observed_species = float(np.sum(positions * counts) / n_species)
    expected_entry_sum = 0.0
    expected_species_sum = 0.0
    richness_bins = _richness_bins(source_richness)

    if matching == "prevalence_richness":
        class_matrix = np.column_stack([prevalence[candidate], richness_bins[candidate]])
    elif matching == "prevalence_richness_effort":
        class_matrix = np.column_stack(
            [prevalence[candidate], richness_bins[candidate], effort_bins[candidate]]
        )
    else:
        raise ValueError(f"unknown v2 matching scheme: {matching}")

    for key in np.unique(class_matrix, axis=0):
        class_mask = candidate & (prevalence == int(key[0])) & (richness_bins == int(key[1]))
        if matching == "prevalence_richness_effort":
            class_mask &= effort_bins == int(key[2])
        represented_genera = int(np.sum(represented & class_mask))
        represented_species = int(np.sum(counts[class_mask]))
        if represented_genera == 0 and represented_species == 0:
            continue
        class_mean = float(np.mean(positions[class_mask]))
        expected_entry_sum += represented_genera * class_mean
        expected_species_sum += represented_species * class_mean

    expected_entry = expected_entry_sum / n_represented
    expected_species = expected_species_sum / n_species
    entry = observed_entry - expected_entry
    species = observed_species - expected_species
    return {
        "n_represented_genera": n_represented,
        "n_represented_species": n_species,
        "observed_entry": observed_entry,
        "expected_entry": expected_entry,
        "entry_enrichment": entry,
        "observed_species": observed_species,
        "expected_species": expected_species,
        "species_enrichment": species,
        "loading_increment": species - entry,
    }


def build_island_breadth_enrichment_v2(
    *,
    genus_breadth: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required = {
        "genus",
        "n_independent_references",
        "sampled_channel_count",
        "effective_channel_number",
    }
    if missing := required - set(genus_breadth.columns):
        raise ValueError(f"genus breadth lacks v2 columns: {sorted(missing)}")
    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: index for index, island in enumerate(islands)}
    source_modes = [str(x) for x in config["source_matching"]["source_modes"]]
    strata = [str(x) for x in config["source_matching"]["strata"]]
    matchings = [
        str(config["source_matching"]["primary_matching"]),
        str(config["source_matching"]["sensitivity_matching"]),
    ]
    minimum = int(config["source_matching"]["minimum_represented_genera"])
    metrics = [str(x) for x in config["model"]["response_metrics"]]
    parts: list[pd.DataFrame] = []

    for threshold in _all_effort_thresholds(config):
        eligible = genus_breadth.copy()
        eligible["n_independent_references"] = pd.to_numeric(
            eligible["n_independent_references"], errors="coerce"
        )
        eligible = eligible.loc[eligible["n_independent_references"].ge(threshold)].copy()
        eligible = eligible.dropna(subset=metrics)
        genera = sorted(set(eligible["genus"].astype(str)) - {""})
        if not genera:
            continue
        genus_index = {genus: index for index, genus in enumerate(genera)}
        ordered = eligible.drop_duplicates("genus").set_index("genus").loc[genera]
        positions = {metric: ordered[metric].to_numpy(float) for metric in metrics}
        effort_bins = reference_effort_bins(ordered["n_independent_references"].to_numpy(float))

        availability = broad_source_availability(gift_flora, set(genera))
        entities = sorted(
            set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
            | set(pd.to_numeric(availability["entity_ID"], errors="coerce").dropna().astype(int))
        )
        entity_index = {entity: index for index, entity in enumerate(entities)}
        presence, richness = _availability_matrices(availability, entity_index, genus_index)
        counts_by_stratum = {
            stratum: _representation_count_matrix(
                status_flora,
                island_index,
                genus_index,
                stratum=stratum,
                species_filter=None,
            ).toarray()
            for stratum in strata
        }
        for source_mode in source_modes:
            assignment = _source_assignment_matrix(
                assignments, island_index, entity_index, source_mode=source_mode
            )
            prevalence = (assignment @ presence).toarray().astype(np.int16)
            source_richness = (assignment @ richness).toarray().astype(np.float32)
            for stratum in strata:
                counts = counts_by_stratum[stratum]
                for matching in matchings:
                    for metric in metrics:
                        rows: list[dict[str, Any]] = []
                        for island_position, island_id in enumerate(islands):
                            result = compute_effort_matched_enrichment(
                                prevalence[island_position],
                                source_richness[island_position],
                                counts[island_position],
                                positions[metric],
                                effort_bins,
                                matching=matching,
                                minimum_represented_genera=minimum,
                            )
                            if result is None:
                                continue
                            rows.append(
                                {
                                    "island_id": island_id,
                                    "metric": metric,
                                    "min_independent_references": threshold,
                                    "stratum": stratum,
                                    "source_mode": source_mode,
                                    "source_matching": matching,
                                    **result,
                                }
                            )
                        if rows:
                            parts.append(pd.DataFrame(rows))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def classify_primary_v2(slopes: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    primary = config["primary_prediction"]
    modes = [str(x) for x in config["source_matching"]["source_modes"]]
    threshold = int(primary["effort_threshold"])
    sensitivities = [int(x) for x in config["sampling_effort_gate"]["sensitivity_min_independent_references"]]
    primary_matching = str(primary["source_matching"])
    matching_sensitivity = str(config["source_matching"]["sensitivity_matching"])
    rows: list[dict[str, Any]] = []
    for context in [str(x) for x in config["model"]["contexts"]]:
        for stratum in [str(x) for x in config["source_matching"]["strata"]]:
            target = slopes.loc[
                slopes["metric"].astype(str).eq(str(primary["metric"]))
                & slopes["outcome"].astype(str).eq(str(primary["outcome"]))
                & slopes["context"].astype(str).eq(context)
                & slopes["stratum"].astype(str).eq(stratum)
            ].copy()
            main = target.loc[
                target["source_matching"].astype(str).eq(primary_matching)
                & target["min_independent_references"].eq(threshold)
            ]
            main_ok = (
                set(main["source_mode"].astype(str)) == set(modes)
                and main["support_class"].eq("confirmatory").all()
                and main["distance_slope"].gt(0).all()
                and main["q_value"].le(0.05).all()
            )
            effort_sign_ok = True
            for sensitivity in sensitivities:
                part = target.loc[
                    target["source_matching"].astype(str).eq(primary_matching)
                    & target["min_independent_references"].eq(sensitivity)
                ]
                if set(part["source_mode"].astype(str)) != set(modes) or not part["distance_slope"].gt(0).all():
                    effort_sign_ok = False
            matching_sensitivity_rows = target.loc[
                target["source_matching"].astype(str).eq(matching_sensitivity)
                & target["min_independent_references"].eq(threshold)
            ]
            matching_sign_ok = (
                set(matching_sensitivity_rows["source_mode"].astype(str)) == set(modes)
                and matching_sensitivity_rows["distance_slope"].gt(0).all()
            )
            promoted = bool(main_ok and effort_sign_ok and matching_sign_ok)
            rows.append(
                {
                    "context": context,
                    "stratum": stratum,
                    "n_primary_source_modes": int(main["source_mode"].astype(str).nunique()),
                    "primary_all_positive": bool(len(main) and main["distance_slope"].gt(0).all()),
                    "primary_all_q_le_005": bool(len(main) and main["q_value"].le(0.05).all()),
                    "effort_sensitivity_no_sign_reversal": bool(effort_sign_ok),
                    "matching_sensitivity_no_sign_reversal": bool(matching_sign_ok),
                    "classification": (
                        "robust_effort_matched_source_breadth_enrichment"
                        if promoted
                        else "not_promoted"
                    ),
                }
            )
    return pd.DataFrame(rows)


def analyse_source_breadth_v2(
    *,
    genus_breadth_csv: Path,
    gift_flora_csv: Path,
    assignments_csv: Path,
    status_flora_csv: Path,
    covariates_csv: Path,
    config_path: Path,
    output_dir: Path,
    predictor_sha256: str = "",
) -> dict[str, Any]:
    config = load_config(config_path)
    observed_sha = sha256_file(genus_breadth_csv)
    if predictor_sha256 and observed_sha != predictor_sha256:
        raise ValueError("outcome-blind genus predictor SHA changed before v2 island join")
    breadth = pd.read_csv(genus_breadth_csv)
    gift = pd.read_csv(gift_flora_csv)
    assignments = pd.read_csv(assignments_csv)
    status = pd.read_csv(status_flora_csv)
    covariates = pd.read_csv(covariates_csv)
    enrichment = build_island_breadth_enrichment_v2(
        genus_breadth=breadth,
        gift_flora=gift,
        assignments=assignments,
        status_flora=status,
        covariates=covariates,
        config=config,
    )
    slopes = fit_breadth_models(enrichment, covariates, config)
    classes = classify_primary_v2(slopes, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    enrichment.to_csv(output_dir / "globi_source_breadth_v2_island_enrichment.csv.gz", index=False, compression="gzip")
    slopes.to_csv(output_dir / "globi_source_breadth_v2_context_slopes.csv", index=False)
    classes.to_csv(output_dir / "globi_source_breadth_v2_primary_classification.csv", index=False)
    manifest = {
        "contract": config["contract"],
        "status": "effort_matched_secondary_H3_analysis_complete",
        "predictor_sha256": observed_sha,
        "n_genus_breadth_rows": int(len(breadth)),
        "n_island_enrichment_rows": int(len(enrichment)),
        "n_slope_rows": int(len(slopes)),
        "primary_promoted_cells": int(
            classes["classification"].eq("robust_effort_matched_source_breadth_enrichment").sum()
        ),
        "N1_rescued": False,
        "N2_opened": False,
        "claim_boundary": str(config["claim_ceiling"]),
    }
    (output_dir / "globi_source_breadth_v2_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command()
def run(
    genus_breadth_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    gift_flora_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    assignments_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    status_flora_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    covariates_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(Path("config/chapter1_globi_source_breadth_v2.yml"), exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
    predictor_sha256: str = typer.Option(""),
) -> None:
    typer.echo(
        json.dumps(
            analyse_source_breadth_v2(
                genus_breadth_csv=genus_breadth_csv,
                gift_flora_csv=gift_flora_csv,
                assignments_csv=assignments_csv,
                status_flora_csv=status_flora_csv,
                covariates_csv=covariates_csv,
                config_path=config_path,
                output_dir=output_dir,
                predictor_sha256=predictor_sha256,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
