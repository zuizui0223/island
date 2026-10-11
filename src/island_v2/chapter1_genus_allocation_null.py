"""Mechanical null for genus allocation in the source-matched decomposition.

This analysis is deliberately outside the active Chapter 1 submission surface.

The source-matched decomposition writes the trait-specific isolation slope of total
species sorting as the exact sum of a genus-structure slope and a within-genus species
sorting slope. Because the within-genus expectation conditions on genus x source-
prevalence cells, sparse cells and phylogenetic trait conservation can mechanically push
variation into the genus term.

This module quantifies that mechanical allocation using the *real* source pools. For each
context, evidence scope, atomic trait and frozen source mode, binary trait states are
shuffled among source-evaluable species within genus while preserving island membership,
source assignment, source prevalence, genus membership and the number of positive states
per genus. The island assemblage itself is held fixed.

Inference is trait-specific. No direction-free joint Wald statistic is used.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import _bh
from island_v2.chapter1_current_h1_source_lineage_entry import (
    SOURCE_MODES,
    build_outcome_states,
    match_gift_species,
)
from island_v2.chapter1_species_sorting_identifiability import _outcome_matrices

app = typer.Typer(add_completion=False, no_args_is_help=True)


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if (
        not isinstance(cfg, dict)
        or cfg.get("contract") != "chapter1_genus_allocation_null_v1"
    ):
        raise typer.BadParameter("unexpected genus-allocation null contract")
    return cfg


def _z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("predictor has zero or non-finite standard deviation")
    return (values - mean) / sd


def _model_slope_weights(
    *,
    covariates: pd.DataFrame,
    islands: list[str],
    n_observed: np.ndarray,
    minimum_islands: int,
) -> tuple[np.ndarray, int] | None:
    """Return FWL weights for the standardized isolation coefficient.

    For a response vector y on the retained islands, beta_distance = weights @ y.
    This is algebraically identical to equal-island OLS with an intercept, standardized
    island area, climate PC1-PC4 and standardized log distance.
    """

    columns = [
        "island_id",
        "log_island_area_km2",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
        "log_distance_to_continent_km",
    ]
    cov = (
        covariates[columns]
        .drop_duplicates("island_id")
        .assign(island_id=lambda x: x["island_id"].astype(str))
        .set_index("island_id")
        .reindex(islands)
    )
    numeric = cov.drop(columns=[]).copy()
    for column in columns[1:]:
        numeric[column] = pd.to_numeric(numeric[column], errors="coerce")

    matrix = numeric[columns[1:]].to_numpy(float)
    complete = np.isfinite(matrix).all(axis=1)
    mask = (np.asarray(n_observed, dtype=int) > 0) & complete
    if int(mask.sum()) < int(minimum_islands):
        return None

    work = numeric.loc[mask, columns[1:]]
    controls = np.column_stack(
        [
            np.ones(len(work), dtype=float),
            _z(work["log_island_area_km2"].to_numpy(float)),
            _z(work["climate_pc1"].to_numpy(float)),
            _z(work["climate_pc2"].to_numpy(float)),
            _z(work["climate_pc3"].to_numpy(float)),
            _z(work["climate_pc4"].to_numpy(float)),
        ]
    )
    distance = _z(work["log_distance_to_continent_km"].to_numpy(float))
    gamma = np.linalg.lstsq(controls, distance, rcond=None)[0]
    residual = distance - controls @ gamma
    denominator = float(residual @ residual)
    if not math.isfinite(denominator) or denominator <= 1e-12:
        return None

    weights = np.zeros(len(islands), dtype=float)
    weights[mask] = residual / denominator
    return weights, int(mask.sum())


def _csr_slots(matrix) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return row, column and data arrays for nonzero CSR entries."""
    matrix = matrix.tocsr(copy=True)
    matrix.sort_indices()
    rows = np.repeat(np.arange(matrix.shape[0], dtype=np.int32), np.diff(matrix.indptr))
    return rows, matrix.indices.astype(np.int32), matrix.data


def _match_observed_to_candidate(
    *,
    candidate_island: np.ndarray,
    candidate_species: np.ndarray,
    candidate_prevalence: np.ndarray,
    observed_island: np.ndarray,
    observed_species: np.ndarray,
    n_species: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    candidate_key = (
        candidate_island.astype(np.int64) * int(n_species)
        + candidate_species.astype(np.int64)
    )
    observed_key = (
        observed_island.astype(np.int64) * int(n_species)
        + observed_species.astype(np.int64)
    )

    candidate_order = np.argsort(candidate_key, kind="mergesort")
    sorted_key = candidate_key[candidate_order]
    positions = np.searchsorted(sorted_key, observed_key)
    valid = positions < len(sorted_key)
    matched = np.zeros(len(observed_key), dtype=bool)
    matched[valid] = sorted_key[positions[valid]] == observed_key[valid]
    if not np.any(matched):
        return (
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int16),
        )

    candidate_position = candidate_order[positions[matched]]
    return (
        observed_island[matched].astype(np.int32),
        observed_species[matched].astype(np.int32),
        candidate_prevalence[candidate_position].astype(np.int16),
    )


def _group_indices(raw_keys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unique, inverse = np.unique(raw_keys.astype(np.int64), return_inverse=True)
    return unique, inverse.astype(np.int32)


def _map_groups(unique_keys: np.ndarray, raw_keys: np.ndarray) -> np.ndarray:
    positions = np.searchsorted(unique_keys, raw_keys.astype(np.int64))
    if np.any(positions >= len(unique_keys)):
        raise RuntimeError("observed slot maps outside candidate group universe")
    if np.any(unique_keys[positions] != raw_keys.astype(np.int64)):
        raise RuntimeError("observed slot has no matching candidate group")
    return positions.astype(np.int32)


def _linear_slope_coefficients(
    *,
    n_species: int,
    n_islands: int,
    candidate_island: np.ndarray,
    candidate_species: np.ndarray,
    candidate_prevalence: np.ndarray,
    observed_island: np.ndarray,
    observed_species: np.ndarray,
    observed_prevalence: np.ndarray,
    genus_codes: np.ndarray,
    slope_weights: np.ndarray,
) -> dict[str, Any]:
    """Build species-level linear coefficients for each decomposition slope."""

    n_observed = np.bincount(observed_island, minlength=n_islands).astype(float)
    if np.any(n_observed[observed_island] <= 0):
        raise RuntimeError("observed island denominator is zero")

    max_prevalence = int(candidate_prevalence.max())
    width = max_prevalence + 1

    source_raw = (
        candidate_island.astype(np.int64) * width
        + candidate_prevalence.astype(np.int64)
    )
    source_keys, source_inverse = _group_indices(source_raw)
    source_count = np.bincount(source_inverse, minlength=len(source_keys)).astype(float)
    observed_source_raw = (
        observed_island.astype(np.int64) * width
        + observed_prevalence.astype(np.int64)
    )
    observed_source_group = _map_groups(source_keys, observed_source_raw)
    observed_source_count = np.bincount(
        observed_source_group, minlength=len(source_keys)
    ).astype(float)

    n_genera = int(genus_codes.max()) + 1 if len(genus_codes) else 0
    candidate_genus = genus_codes[candidate_species].astype(np.int64)
    observed_genus = genus_codes[observed_species].astype(np.int64)
    genus_raw = (
        (
            candidate_island.astype(np.int64) * n_genera
            + candidate_genus
        )
        * width
        + candidate_prevalence.astype(np.int64)
    )
    genus_keys, genus_inverse = _group_indices(genus_raw)
    genus_count = np.bincount(genus_inverse, minlength=len(genus_keys)).astype(float)
    observed_genus_raw = (
        (
            observed_island.astype(np.int64) * n_genera
            + observed_genus
        )
        * width
        + observed_prevalence.astype(np.int64)
    )
    observed_genus_group = _map_groups(genus_keys, observed_genus_raw)
    observed_genus_count = np.bincount(
        observed_genus_group, minlength=len(genus_keys)
    ).astype(float)

    raw_weight = slope_weights[observed_island] / n_observed[observed_island]
    raw_coef = np.bincount(
        observed_species,
        weights=raw_weight,
        minlength=n_species,
    ).astype(float)

    candidate_observed_denominator = n_observed[candidate_island]
    source_numerator = (
        slope_weights[candidate_island]
        * observed_source_count[source_inverse]
    )
    source_denominator = (
        candidate_observed_denominator
        * source_count[source_inverse]
    )
    source_candidate_weight = np.divide(
        source_numerator,
        source_denominator,
        out=np.zeros_like(source_numerator, dtype=float),
        where=source_denominator > 0,
    )
    source_coef = np.bincount(
        candidate_species,
        weights=source_candidate_weight,
        minlength=n_species,
    ).astype(float)

    genus_numerator = (
        slope_weights[candidate_island]
        * observed_genus_count[genus_inverse]
    )
    genus_denominator = (
        candidate_observed_denominator
        * genus_count[genus_inverse]
    )
    genus_candidate_weight = np.divide(
        genus_numerator,
        genus_denominator,
        out=np.zeros_like(genus_numerator, dtype=float),
        where=genus_denominator > 0,
    )
    genus_coef = np.bincount(
        candidate_species,
        weights=genus_candidate_weight,
        minlength=n_species,
    ).astype(float)

    return {
        "raw": raw_coef,
        "source": source_coef,
        "genus_expectation": genus_coef,
        "genus_structure": genus_coef - source_coef,
        "within_genus": raw_coef - genus_coef,
        "total_sorting": raw_coef - source_coef,
        "n_observed": n_observed,
        "genus_inverse": genus_inverse,
        "genus_count": genus_count,
        "observed_genus_group": observed_genus_group,
    }


def _allocation_share(beta_genus: float, beta_within: float) -> float:
    denominator = abs(float(beta_genus)) + abs(float(beta_within))
    if denominator <= 1e-12:
        return float("nan")
    return abs(float(beta_genus)) / denominator


def _signed_fraction(beta_genus: float, beta_total: float) -> float:
    if abs(float(beta_total)) <= 1e-12:
        return float("nan")
    return float(beta_genus) / float(beta_total)


def _shuffle_groups(
    *,
    candidate_species: np.ndarray,
    genus_codes: np.ndarray,
    states: np.ndarray,
) -> list[np.ndarray]:
    """Return variable within-genus source-union species groups."""
    union_species = np.unique(candidate_species)
    groups: list[np.ndarray] = []
    for genus in np.unique(genus_codes[union_species]):
        indices = union_species[genus_codes[union_species] == genus]
        if len(indices) < 2:
            continue
        values = states[indices]
        if float(np.min(values)) == float(np.max(values)):
            continue
        groups.append(indices.astype(np.int32))
    return groups


def _cell_seed(base_seed: int, key: str) -> int:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return (int(base_seed) + int(digest[:8], 16)) % (2**32 - 1)


def _estimability(
    *,
    states: np.ndarray,
    candidate_species: np.ndarray,
    genus_inverse: np.ndarray,
    genus_count: np.ndarray,
    observed_genus_group: np.ndarray,
) -> dict[str, float | int]:
    state_sum = np.bincount(
        genus_inverse,
        weights=states[candidate_species],
        minlength=len(genus_count),
    )
    cell_n_ge_2 = genus_count >= 2
    cell_variable = (state_sum > 0) & (state_sum < genus_count)
    estimable_cell = cell_n_ge_2 & cell_variable

    observed_cell = observed_genus_group
    slot_n_ge_2 = cell_n_ge_2[observed_cell]
    slot_variable = cell_variable[observed_cell]
    slot_estimable = estimable_cell[observed_cell]

    unique_observed_cells = np.unique(observed_cell)
    n_slots = int(len(observed_cell))
    n_cells = int(len(unique_observed_cells))
    return {
        "observed_source_candidate_slots": n_slots,
        "slots_with_cell_n_ge_2": int(slot_n_ge_2.sum()),
        "slots_with_state_variation": int(slot_variable.sum()),
        "estimable_slots": int(slot_estimable.sum()),
        "estimable_slot_fraction": (
            float(slot_estimable.mean()) if n_slots else float("nan")
        ),
        "observed_genus_prevalence_cells": n_cells,
        "cells_with_n_ge_2": int(cell_n_ge_2[unique_observed_cells].sum()),
        "cells_with_state_variation": int(cell_variable[unique_observed_cells].sum()),
        "estimable_cells": int(estimable_cell[unique_observed_cells].sum()),
        "estimable_cell_fraction": (
            float(estimable_cell[unique_observed_cells].mean())
            if n_cells
            else float("nan")
        ),
    }


def run_cell(
    *,
    states: np.ndarray,
    genus_codes: np.ndarray,
    source_presence,
    observed,
    assignment,
    islands: list[str],
    covariates: pd.DataFrame,
    n_permutations: int,
    base_seed: int,
    cell_key: str,
    minimum_islands: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    prevalence = (assignment @ source_presence).tocsr()
    prevalence.sort_indices()
    observed = observed.tocsr(copy=True)
    observed.sort_indices()

    candidate_island, candidate_species, candidate_prevalence_raw = _csr_slots(
        prevalence
    )
    candidate_prevalence = candidate_prevalence_raw.astype(np.int16)
    observed_island_all, observed_species_all, _ = _csr_slots(observed)

    if len(candidate_species) == 0 or len(observed_species_all) == 0:
        return (
            {"status": "not_testable", "reason": "empty_source_or_observed_support"},
            {},
        )

    observed_island, observed_species, observed_prevalence = (
        _match_observed_to_candidate(
            candidate_island=candidate_island,
            candidate_species=candidate_species,
            candidate_prevalence=candidate_prevalence,
            observed_island=observed_island_all,
            observed_species=observed_species_all,
            n_species=len(states),
        )
    )
    if len(observed_species) == 0:
        return (
            {"status": "not_testable", "reason": "no_source_evaluable_observed_slots"},
            {},
        )

    n_observed = np.bincount(
        observed_island, minlength=len(islands)
    ).astype(int)
    model = _model_slope_weights(
        covariates=covariates,
        islands=islands,
        n_observed=n_observed,
        minimum_islands=minimum_islands,
    )
    if model is None:
        return (
            {
                "status": "not_testable",
                "reason": "minimum_island_or_model_support_gate",
                "n_islands_with_source_slots": int(np.sum(n_observed > 0)),
            },
            {},
        )
    slope_weights, n_model_islands = model

    coeff = _linear_slope_coefficients(
        n_species=len(states),
        n_islands=len(islands),
        candidate_island=candidate_island,
        candidate_species=candidate_species,
        candidate_prevalence=candidate_prevalence,
        observed_island=observed_island,
        observed_species=observed_species,
        observed_prevalence=observed_prevalence,
        genus_codes=genus_codes,
        slope_weights=slope_weights,
    )

    beta_genus = float(coeff["genus_structure"] @ states)
    beta_within = float(coeff["within_genus"] @ states)
    beta_total = float(coeff["total_sorting"] @ states)
    closure_error = float(beta_total - beta_genus - beta_within)
    if abs(closure_error) > 1e-10:
        raise RuntimeError(f"slope closure failed: {closure_error}")

    observed_share = _allocation_share(beta_genus, beta_within)
    observed_signed_fraction = _signed_fraction(beta_genus, beta_total)

    groups = _shuffle_groups(
        candidate_species=candidate_species,
        genus_codes=genus_codes,
        states=states,
    )
    rng = np.random.default_rng(_cell_seed(base_seed, cell_key))
    null_share = np.full(n_permutations, np.nan, dtype=float)
    null_genus = np.full(n_permutations, np.nan, dtype=float)
    null_total = np.full(n_permutations, np.nan, dtype=float)
    null_within = np.full(n_permutations, np.nan, dtype=float)

    for index in range(n_permutations):
        shuffled = states.copy()
        for group in groups:
            shuffled[group] = rng.permutation(states[group])
        bg = float(coeff["genus_structure"] @ shuffled)
        bw = float(coeff["within_genus"] @ shuffled)
        bt = float(coeff["total_sorting"] @ shuffled)
        null_genus[index] = bg
        null_within[index] = bw
        null_total[index] = bt
        null_share[index] = _allocation_share(bg, bw)

    valid_share = null_share[np.isfinite(null_share)]
    if math.isfinite(observed_share) and len(valid_share):
        p_share = float(
            (1 + np.sum(valid_share >= observed_share - 1e-15))
            / (1 + len(valid_share))
        )
    else:
        p_share = float("nan")

    p_abs_genus = float(
        (1 + np.sum(np.abs(null_genus) >= abs(beta_genus) - 1e-15))
        / (1 + n_permutations)
    )

    result = {
        "status": "fit",
        "n_model_islands": n_model_islands,
        "n_source_candidate_species_union": int(
            len(np.unique(candidate_species))
        ),
        "n_variable_shuffle_genera": int(len(groups)),
        "observed_beta_total_sorting": beta_total,
        "observed_beta_genus_structure": beta_genus,
        "observed_beta_within_genus": beta_within,
        "slope_closure_error": closure_error,
        "observed_absolute_genus_allocation_share": observed_share,
        "observed_signed_genus_fraction": observed_signed_fraction,
        "n_permutations": int(n_permutations),
        "null_share_n_valid": int(len(valid_share)),
        "null_share_mean": (
            float(np.mean(valid_share)) if len(valid_share) else float("nan")
        ),
        "null_share_median": (
            float(np.median(valid_share)) if len(valid_share) else float("nan")
        ),
        "null_share_q025": (
            float(np.quantile(valid_share, 0.025))
            if len(valid_share)
            else float("nan")
        ),
        "null_share_q975": (
            float(np.quantile(valid_share, 0.975))
            if len(valid_share)
            else float("nan")
        ),
        "p_observed_share_exceeds_null": p_share,
        "null_abs_genus_mean": float(np.mean(np.abs(null_genus))),
        "p_observed_abs_genus_exceeds_null": p_abs_genus,
        "null_total_mean": float(np.mean(null_total)),
        "null_within_mean": float(np.mean(null_within)),
    }
    estimability = _estimability(
        states=states,
        candidate_species=candidate_species,
        genus_inverse=coeff["genus_inverse"],
        genus_count=coeff["genus_count"],
        observed_genus_group=coeff["observed_genus_group"],
    )
    return result, estimability


def _add_bh(table: pd.DataFrame) -> pd.DataFrame:
    result = table.copy()
    result["q_share_within_trait_family"] = np.nan
    result["q_abs_genus_within_trait_family"] = np.nan
    fit = result["status"].astype(str).eq("fit")
    for _, index in result.loc[fit].groupby(
        ["evidence_scope", "context", "source_mode"]
    ).groups.items():
        result.loc[index, "q_share_within_trait_family"] = _bh(
            result.loc[index, "p_observed_share_exceeds_null"]
        )
        result.loc[index, "q_abs_genus_within_trait_family"] = _bh(
            result.loc[index, "p_observed_abs_genus_exceeds_null"]
        )
    return result


@app.command()
def main(
    state_audit_all_csv: Path = typer.Option(..., exists=True),
    state_audit_direct_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    source_assignments_csv: Path = typer.Option(..., exists=True),
    corrected_covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    null_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(null_config_path)
    probability_config = yaml.safe_load(
        probability_config_path.read_text(encoding="utf-8")
    )
    all_audit = pd.read_csv(state_audit_all_csv)
    direct_audit = pd.read_csv(state_audit_direct_csv)
    status_flora = pd.read_csv(status_flora_csv)
    gift_flora = pd.read_csv(gift_flora_csv)
    assignments = pd.read_csv(source_assignments_csv)
    covariates = pd.read_csv(corrected_covariates_csv)
    covariates["island_id"] = covariates["island_id"].astype(str)

    all_states = build_outcome_states(all_audit, probability_config)
    matched_gift = match_gift_species(
        gift_flora,
        all_states["accepted_species"],
    )

    contexts = [str(x) for x in config["scope"]["contexts"]]
    outcomes = [str(x) for x in config["atomic_outcomes"]]
    floristic_status = str(config["scope"]["floristic_status"])
    source_modes = [str(x) for x in config["scope"]["source_modes"]]
    if tuple(source_modes) != tuple(SOURCE_MODES):
        raise typer.BadParameter("null contract source modes differ from frozen modes")

    n_permutations = int(config["permutation_null"]["n_permutations"])
    base_seed = int(config["permutation_null"]["seed"])
    minimum_islands = int(config["model"]["minimum_islands_per_outcome"])

    result_rows: list[dict[str, Any]] = []
    estimability_rows: list[dict[str, Any]] = []

    for evidence_scope, audit in (
        ("all_analysis_eligible", all_audit),
        ("direct_only", direct_audit),
    ):
        states = build_outcome_states(audit, probability_config)
        scope_species = set(states["accepted_species"].astype(str))
        scope_matched = matched_gift.loc[
            matched_gift["accepted_species"].astype(str).isin(scope_species)
        ].copy()

        for context in contexts:
            context_covariates = covariates.loc[
                covariates["analysis_regime"].astype(str).eq(context)
            ].copy()
            islands = sorted(
                context_covariates["island_id"].astype(str).unique()
            )
            if not islands:
                continue

            for outcome in outcomes:
                (
                    state_array,
                    genus_codes,
                    source_presence,
                    observed,
                    assignment_matrices,
                    species,
                ) = _outcome_matrices(
                    states=states,
                    outcome=outcome,
                    matched_gift=scope_matched,
                    status_flora=status_flora,
                    islands=islands,
                    floristic_status=floristic_status,
                    assignments=assignments,
                )
                if not species:
                    continue

                for source_mode in source_modes:
                    key = "|".join(
                        [evidence_scope, context, outcome, source_mode]
                    )
                    result, estimability = run_cell(
                        states=state_array,
                        genus_codes=genus_codes,
                        source_presence=source_presence,
                        observed=observed,
                        assignment=assignment_matrices[source_mode],
                        islands=islands,
                        covariates=context_covariates,
                        n_permutations=n_permutations,
                        base_seed=base_seed,
                        cell_key=key,
                        minimum_islands=minimum_islands,
                    )
                    result_rows.append(
                        {
                            "evidence_scope": evidence_scope,
                            "context": context,
                            "outcome": outcome,
                            "source_mode": source_mode,
                            **result,
                        }
                    )
                    if estimability:
                        estimability_rows.append(
                            {
                                "evidence_scope": evidence_scope,
                                "context": context,
                                "outcome": outcome,
                                "source_mode": source_mode,
                                **estimability,
                            }
                        )

    results = _add_bh(pd.DataFrame(result_rows))
    estimability = pd.DataFrame(estimability_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "genus_allocation_null_results.csv", index=False)
    estimability.to_csv(output_dir / "within_genus_estimability.csv", index=False)

    fit = results.loc[results["status"].astype(str).eq("fit")].copy()
    summary_rows: list[dict[str, Any]] = []
    if not fit.empty:
        for (scope, context, outcome), part in fit.groupby(
            ["evidence_scope", "context", "outcome"]
        ):
            summary_rows.append(
                {
                    "evidence_scope": scope,
                    "context": context,
                    "outcome": outcome,
                    "n_source_modes_fit": int(len(part)),
                    "observed_share_min": float(
                        part["observed_absolute_genus_allocation_share"].min()
                    ),
                    "observed_share_max": float(
                        part["observed_absolute_genus_allocation_share"].max()
                    ),
                    "null_share_mean_min": float(part["null_share_mean"].min()),
                    "null_share_mean_max": float(part["null_share_mean"].max()),
                    "p_share_max": float(
                        part["p_observed_share_exceeds_null"].max()
                    ),
                    "q_share_max": float(
                        part["q_share_within_trait_family"].max()
                    ),
                    "all_four_modes_observed_share_above_null_mean": bool(
                        len(part) == 4
                        and (
                            part["observed_absolute_genus_allocation_share"]
                            > part["null_share_mean"]
                        ).all()
                    ),
                    "all_four_modes_share_BH_supported": bool(
                        len(part) == 4
                        and part["q_share_within_trait_family"].le(0.05).all()
                    ),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dir / "genus_allocation_null_summary.csv", index=False)

    result_manifest = {
        "contract": config["contract"],
        "status": "complete",
        "n_result_rows": int(len(results)),
        "n_fit_rows": int(len(fit)),
        "n_estimability_rows": int(len(estimability)),
        "n_permutations": n_permutations,
        "submission_inference": False,
        "interpretation_rule": (
            "Observed genus allocation may be interpreted as exceeding the mechanical "
            "null only trait-by-trait and source-mode-by-source-mode. Four source modes "
            "are correlated sensitivity definitions, not independent replications."
        ),
        "claim_boundary": (
            "This null can diagnose mechanical allocation caused by real source-pool "
            "architecture, trait conservation and sparse genus-by-prevalence cells. "
            "It does not identify a causal genus-level assembly process."
        ),
    }
    (output_dir / "RESULT.json").write_text(
        json.dumps(result_manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(summary.to_csv(index=False))


if __name__ == "__main__":
    app()
