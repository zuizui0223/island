"""Three-axis raw-state H1 reanalysis using the original species×axis cells.

This analysis keeps the Chapter 1 database at its original three measurement domains:
flower colour, floral structural complexity, and reproductive assurance.  Component
traits are measurement items, not seven complete trait axes.

Every resolved species×axis cell contributes through any ontology-valid component state
that is actually observed.  Missing components are not coded as zero and complete cases
are never required.  Raw multistate memberships are retained.  Formal tests are
multivariate beta-binomial Wald tests of all estimable raw state-prevalence responses
within each axis and geographic region.

Floristic-status partitions are fitted with the identical axis model so that introduced
species cannot be treated as a separate methodological problem.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import (
    _assemble_cluster_covariance,
    _bh,
    _chi_square_sf_integer_df,
    _fit_single_beta_binomial,
    _fit_within,
    _normal_two_sided_p,
    _prepare,
    _standardize,
)
from island_v2.chapter1_wcvp_partition_diagnostic import classify_wcvp_partitions

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _parse_composition(
    text: object,
    *,
    trait_to_axis: dict[str, str],
    allowed: dict[str, set[str]],
) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    for item in str(text or "").split("|"):
        trait, sep, raw_states = item.partition("=")
        trait = trait.strip()
        if not sep or trait not in trait_to_axis:
            continue
        try:
            states = json.loads(raw_states)
        except json.JSONDecodeError:
            continue
        if not isinstance(states, list):
            continue
        for state in states:
            state = str(state).strip()
            if state and state != "unresolved" and state in allowed.get(trait, set()):
                rows.append((trait_to_axis[trait], trait, state))
    return sorted(set(rows))


def build_valid_state_ledger(
    species_axis: pd.DataFrame,
    ontology: dict[str, Any],
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    qualities = set(config["evidence_scopes"][evidence_scope])
    trait_to_axis = {
        str(trait): str(axis)
        for axis, spec in config["axes"].items()
        for trait in spec["traits"]
    }
    ontology_traits = ontology.get("traits", {})
    allowed: dict[str, set[str]] = {}
    for trait in trait_to_axis:
        spec = ontology_traits.get(trait, {})
        allowed[trait] = {
            str(x)
            for x in spec.get("allowed_values", [])
            if str(x) != "unresolved"
        }
        if not allowed[trait]:
            raise ValueError(f"no ontology states for {trait}")

    work = species_axis.copy().fillna("")
    required = {"accepted_species", "axis", "trait_composition", "quality"}
    missing = required - set(work.columns)
    if missing:
        raise ValueError(f"species-axis snapshot missing {sorted(missing)}")
    work = work.loc[
        work["quality"].astype(str).isin(qualities)
        & work["trait_composition"].astype(str).ne("")
        & work["axis"].astype(str).isin(config["axes"])
    ].copy()

    state_rows: list[dict[str, str]] = []
    coverage_rows: list[dict[str, Any]] = []
    for row in work.itertuples(index=False):
        valid = _parse_composition(
            row.trait_composition,
            trait_to_axis=trait_to_axis,
            allowed=allowed,
        )
        for axis, trait, state in valid:
            state_rows.append(
                {
                    "accepted_species": str(row.accepted_species),
                    "axis": axis,
                    "trait_name": trait,
                    "state": state,
                }
            )
        coverage_rows.append(
            {
                "accepted_species": str(row.accepted_species),
                "axis": str(row.axis),
                "has_valid_component_state": bool(valid),
                "n_valid_state_memberships": len(valid),
                "n_valid_component_traits": len({trait for _, trait, _ in valid}),
            }
        )

    ledger = pd.DataFrame(
        state_rows,
        columns=["accepted_species", "axis", "trait_name", "state"],
    ).drop_duplicates()
    cell = pd.DataFrame(coverage_rows)
    audit = (
        cell.groupby("axis", as_index=False)
        .agg(
            resolved_axis_cells=("accepted_species", "size"),
            ontology_valid_axis_cells=("has_valid_component_state", "sum"),
            ontology_invalid_only_axis_cells=(
                "has_valid_component_state",
                lambda x: int((~x.astype(bool)).sum()),
            ),
            n_valid_state_memberships=("n_valid_state_memberships", "sum"),
            n_valid_component_trait_records=("n_valid_component_traits", "sum"),
        )
    )
    audit.insert(0, "evidence_scope", evidence_scope)
    return ledger, audit


def flora_scopes(
    status_flora: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    classified = classify_wcvp_partitions(status_flora, wcvp_ranges, island_tdwg)
    masks = {
        "all_observed": pd.Series(True, index=classified.index),
        "source_native": classified["wcvp_partition"].eq("source_native"),
        "source_introduced": classified["wcvp_partition"].eq("source_introduced"),
        "regional_native_compatible": classified["wcvp_partition"].isin(
            {"source_native", "wcvp_compatible_unresolved"}
        ),
        "regionally_incompatible_or_introduced": classified["wcvp_partition"].isin(
            {"source_introduced", "wcvp_incompatible_unresolved"}
        ),
    }
    return {
        label: classified.loc[
            mask,
            ["island_id", "accepted_species", "origin_status", "floristic_status"],
        ].drop_duplicates(["island_id", "accepted_species"])
        for label, mask in masks.items()
    }


def build_axis_state_counts(
    flora: pd.DataFrame,
    state_ledger: pd.DataFrame,
    *,
    stratum: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = flora[["island_id", "accepted_species"]].copy()
    base["island_id"] = base["island_id"].astype(str)
    base["accepted_species"] = base["accepted_species"].astype(str)
    base = base.drop_duplicates(["island_id", "accepted_species"])

    trait_species = state_ledger[
        ["accepted_species", "axis", "trait_name"]
    ].drop_duplicates()
    state_species = state_ledger[
        ["accepted_species", "axis", "trait_name", "state"]
    ].drop_duplicates()

    trait_join = base.merge(
        trait_species,
        on="accepted_species",
        how="inner",
        validate="many_to_many",
    )
    denominators = (
        trait_join.groupby(["island_id", "axis", "trait_name"], as_index=False)
        .agg(trials=("accepted_species", "nunique"))
    )

    state_join = base.merge(
        state_species,
        on="accepted_species",
        how="inner",
        validate="many_to_many",
    )
    successes = (
        state_join.groupby(
            ["island_id", "axis", "trait_name", "state"],
            as_index=False,
        )
        .agg(successes=("accepted_species", "nunique"))
    )

    observed_states = (
        state_ledger[["axis", "trait_name", "state"]]
        .drop_duplicates()
        .sort_values(["axis", "trait_name", "state"])
    )
    global_state_species = (
        state_ledger.groupby(["axis", "trait_name", "state"], as_index=False)
        .agg(state_unique_species_global=("accepted_species", "nunique"))
    )

    parts: list[pd.DataFrame] = []
    for (axis, trait), denom in denominators.groupby(["axis", "trait_name"], sort=False):
        states = observed_states.loc[
            observed_states["axis"].eq(axis)
            & observed_states["trait_name"].eq(trait),
            "state",
        ].tolist()
        if not states:
            continue
        expanded = denom.assign(_key=1).merge(
            pd.DataFrame({"state": states, "_key": 1}),
            on="_key",
            how="inner",
        ).drop(columns="_key")
        succ = successes.loc[
            successes["axis"].eq(axis)
            & successes["trait_name"].eq(trait)
        ]
        expanded = expanded.merge(
            succ[["island_id", "state", "successes"]],
            on=["island_id", "state"],
            how="left",
            validate="one_to_one",
        )
        expanded = expanded.merge(
            global_state_species.loc[
                global_state_species["axis"].eq(axis)
                & global_state_species["trait_name"].eq(trait),
                ["state", "state_unique_species_global"],
            ],
            on="state",
            how="left",
            validate="many_to_one",
        )
        expanded["successes"] = expanded["successes"].fillna(0).astype(int)
        expanded["outcome"] = (
            expanded["trait_name"].astype(str)
            + "::"
            + expanded["state"].astype(str)
        )
        expanded["stratum"] = stratum
        parts.append(
            expanded[
                [
                    "island_id",
                    "successes",
                    "trials",
                    "outcome",
                    "stratum",
                    "axis",
                    "trait_name",
                    "state",
                    "state_unique_species_global",
                ]
            ]
        )

    counts = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    coverage = pd.DataFrame(
        [
            {
                "flora_scope": stratum,
                "n_flora_rows": int(len(base)),
                "n_flora_species": int(base["accepted_species"].nunique()),
                "n_islands": int(base["island_id"].nunique()),
                "n_trait_resolved_island_species_rows": int(len(trait_join)),
                "n_state_membership_island_species_rows": int(len(state_join)),
            }
        ]
    )
    return counts, coverage


def _axis_fit(
    prepared: pd.DataFrame,
    counts: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    flora_scope: str,
    axis: str,
    context_value: str,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    model = config["model"]
    context = str(model["context_column"])
    part = prepared.loc[
        prepared["stratum"].eq(flora_scope)
        & prepared[context].eq(context_value)
    ].copy()
    axis_outcomes = counts.loc[counts["axis"].eq(axis), "outcome"].drop_duplicates()
    part = part.loc[part["outcome"].isin(axis_outcomes)].copy()

    support_rows: list[dict[str, Any]] = []
    eligible: list[str] = []
    for outcome, out in part.groupby("outcome", sort=False):
        n_islands = int(out["island_id"].nunique())
        positive = int(out.loc[out["successes"].gt(0), "island_id"].nunique())
        negative = int(out.loc[out["successes"].lt(out["trials"]), "island_id"].nunique())
        state_species = int(
            pd.to_numeric(
                out["state_unique_species_global"], errors="coerce"
            ).max()
        )
        ok = (
            n_islands >= int(model["minimum_islands_per_state"])
            and positive >= int(model["minimum_positive_islands_per_state"])
            and negative >= int(model["minimum_negative_islands_per_state"])
            and state_species >= int(model["minimum_unique_species_per_state"])
        )
        support_rows.append(
            {
                "evidence_scope": evidence_scope,
                "flora_scope": flora_scope,
                "axis": axis,
                "context": context_value,
                "outcome": str(outcome),
                "n_islands": n_islands,
                "positive_islands": positive,
                "negative_islands": negative,
                "state_unique_species_global": state_species,
                "eligible": ok,
            }
        )
        if ok:
            eligible.append(str(outcome))

    support = pd.DataFrame(support_rows)
    base_cfg = {
        "geography_column": str(model["geography_column"]),
        "context_column": str(model["context_column"]),
        "cluster_column": str(model["cluster_column"]),
        "baseline_covariates": [str(x) for x in model["baseline_covariates"]],
        "max_iter": int(model["max_iter"]),
        "minimum_outcomes_per_vector": int(model["minimum_states_per_axis_test"]),
        "model_outcomes": eligible,
    }
    if len(eligible) < int(model["minimum_states_per_axis_test"]):
        return pd.DataFrame(), {
            "evidence_scope": evidence_scope,
            "flora_scope": flora_scope,
            "axis": axis,
            "context": context_value,
            "status": "not_testable",
            "n_retained_outcomes": len(eligible),
            "retained_outcomes": "|".join(eligible),
        }, support

    slopes, omnibus = _fit_within(
        prepared,
        stratum=flora_scope,
        context_value=context_value,
        threshold=int(model["minimum_islands_per_state"]),
        config=base_cfg,
    )
    retry_used = False
    max_abs_retry_delta = 0.0
    initial_converged = bool(omnibus.get("all_optimizers_converged", False))
    if (
        not initial_converged
        and int(model.get("retry_max_iter", model["max_iter"]))
        > int(model["max_iter"])
    ):
        retry_cfg = dict(base_cfg)
        retry_cfg["max_iter"] = int(model["retry_max_iter"])
        retry_slopes, retry_omnibus = _fit_within(
            prepared,
            stratum=flora_scope,
            context_value=context_value,
            threshold=int(model["minimum_islands_per_state"]),
            config=retry_cfg,
        )
        retry_used = True
        if not slopes.empty and not retry_slopes.empty:
            paired = slopes[["outcome", "geography_slope_log_odds"]].merge(
                retry_slopes[["outcome", "geography_slope_log_odds"]],
                on="outcome",
                suffixes=("_initial", "_retry"),
                validate="one_to_one",
            )
            max_abs_retry_delta = float(
                (
                    paired["geography_slope_log_odds_initial"]
                    - paired["geography_slope_log_odds_retry"]
                )
                .abs()
                .max()
            )
        slopes, omnibus = retry_slopes, retry_omnibus
    if not slopes.empty:
        slopes.insert(0, "axis", axis)
        slopes.insert(0, "flora_scope", flora_scope)
        slopes.insert(0, "evidence_scope", evidence_scope)
    omnibus = {
        "evidence_scope": evidence_scope,
        "flora_scope": flora_scope,
        "axis": axis,
        "initial_all_optimizers_converged": initial_converged,
        "retry_used": retry_used,
        "max_abs_retry_slope_delta": max_abs_retry_delta,
        **omnibus,
    }
    return slopes, omnibus, support


def _fit_status_contrast(
    prepared: pd.DataFrame,
    counts: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
    axis: str,
    context_value: str,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    model = config["model"]
    status_cfg = config["status_contrast"]
    reference = str(status_cfg["reference"])
    comparison = str(status_cfg["comparison"])
    context = str(model["context_column"])
    geography = str(model["geography_column"])
    cluster = str(model["cluster_column"])
    baseline = [str(x) for x in model["baseline_covariates"]]

    axis_outcomes = counts.loc[counts["axis"].eq(axis), "outcome"].drop_duplicates()
    work = prepared.loc[
        prepared["stratum"].isin([reference, comparison])
        & prepared[context].eq(context_value)
        & prepared["outcome"].isin(axis_outcomes)
    ].copy()

    support_rows: list[dict[str, Any]] = []
    eligible: list[str] = []
    for outcome, out in work.groupby("outcome", sort=False):
        ok = True
        row: dict[str, Any] = {
            "evidence_scope": evidence_scope,
            "axis": axis,
            "context": context_value,
            "outcome": str(outcome),
        }
        for scope in (reference, comparison):
            piece = out.loc[out["stratum"].eq(scope)]
            n_islands = int(piece["island_id"].nunique())
            positive = int(
                piece.loc[piece["successes"].gt(0), "island_id"].nunique()
            )
            negative = int(
                piece.loc[
                    piece["successes"].lt(piece["trials"]), "island_id"
                ].nunique()
            )
            state_species = (
                int(
                    pd.to_numeric(
                        piece["state_unique_species_global"], errors="coerce"
                    ).max()
                )
                if not piece.empty
                else 0
            )
            row[f"{scope}_n_islands"] = n_islands
            row[f"{scope}_positive_islands"] = positive
            row[f"{scope}_negative_islands"] = negative
            row[f"{scope}_state_unique_species_global"] = state_species
            ok = ok and (
                n_islands >= int(model["minimum_islands_per_state"])
                and positive >= int(model["minimum_positive_islands_per_state"])
                and negative >= int(model["minimum_negative_islands_per_state"])
                and state_species >= int(model["minimum_unique_species_per_state"])
            )
        row["eligible"] = ok
        support_rows.append(row)
        if ok:
            eligible.append(str(outcome))

    support = pd.DataFrame(support_rows)
    if len(eligible) < int(model["minimum_states_per_axis_test"]):
        return pd.DataFrame(), {
            "evidence_scope": evidence_scope,
            "axis": axis,
            "context": context_value,
            "reference": reference,
            "comparison": comparison,
            "status": "not_testable",
            "n_retained_outcomes": len(eligible),
            "retained_outcomes": "|".join(eligible),
        }, support

    fits: list[dict[str, Any]] = []
    cluster_parts: list[np.ndarray] = []
    interaction_indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in eligible:
        part = work.loc[work["outcome"].eq(outcome)].copy()
        comparison_indicator = part["stratum"].eq(comparison).to_numpy(float)
        columns = [np.ones(len(part), dtype=float), comparison_indicator]
        names = [f"{outcome}:intercept", f"{outcome}:status[{comparison}]"]
        for predictor in baseline:
            z = _standardize(part[predictor])
            columns.extend([z, z * comparison_indicator])
            names.extend(
                [
                    f"{outcome}:z_{predictor}",
                    f"{outcome}:z_{predictor}:status[{comparison}]",
                ]
            )
        z_geo = _standardize(part[geography])
        interaction_name = (
            f"{outcome}:z_{geography}:status[{comparison}]"
        )
        columns.extend([z_geo, z_geo * comparison_indicator])
        names.extend([f"{outcome}:z_{geography}", interaction_name])
        fit = _fit_single_beta_binomial(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            np.column_stack(columns),
            names,
            max_iter=int(model.get("retry_max_iter", model["max_iter"])),
        )
        fits.append(fit)
        cluster_parts.append(part[cluster].astype(str).to_numpy())
        interaction_indices.append(offset + names.index(interaction_name))
        rows.append(
            {
                "evidence_scope": evidence_scope,
                "axis": axis,
                "context": context_value,
                "reference": reference,
                "comparison": comparison,
                "outcome": outcome,
                "n_islands_reference": int(
                    part.loc[
                        part["stratum"].eq(reference), "island_id"
                    ].nunique()
                ),
                "n_islands_comparison": int(
                    part.loc[
                        part["stratum"].eq(comparison), "island_id"
                    ].nunique()
                ),
                "kappa": float(fit["kappa"]),
                "optimizer_success": bool(fit["success"]),
            }
        )
        offset += len(fit["names"])

    covariance, _, theta = _assemble_cluster_covariance(fits, cluster_parts)
    vector = theta[interaction_indices]
    vector_cov = covariance[np.ix_(interaction_indices, interaction_indices)]
    vector_se = np.sqrt(np.clip(np.diag(vector_cov), 0.0, None))
    for row, estimate, stderr in zip(rows, vector, vector_se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        row.update(
            {
                "slope_difference_comparison_minus_reference": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": _normal_two_sided_p(z),
            }
        )

    rank = int(np.linalg.matrix_rank(vector_cov))
    statistic = (
        float(vector @ np.linalg.pinv(vector_cov) @ vector)
        if rank > 0
        else float("nan")
    )
    omnibus = {
        "evidence_scope": evidence_scope,
        "axis": axis,
        "context": context_value,
        "reference": reference,
        "comparison": comparison,
        "status": "fit",
        "n_retained_outcomes": len(eligible),
        "retained_outcomes": "|".join(eligible),
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster].nunique()),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(statistic, rank),
        "all_optimizers_converged": all(bool(f["success"]) for f in fits),
    }
    return pd.DataFrame(rows), omnibus, support


def _status_vector_similarity(slopes: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    status_cfg = config["status_contrast"]
    reference = str(status_cfg["reference"])
    comparison = str(status_cfg["comparison"])
    rows: list[dict[str, Any]] = []
    for evidence_scope in config["evidence_scopes"]:
        for axis in config["axes"]:
            for context_value in config["model"]["contexts"]:
                ref = slopes.loc[
                    slopes["evidence_scope"].eq(evidence_scope)
                    & slopes["flora_scope"].eq(reference)
                    & slopes["axis"].eq(axis)
                    & slopes["context"].eq(context_value),
                    ["outcome", "geography_slope_log_odds"],
                ]
                comp = slopes.loc[
                    slopes["evidence_scope"].eq(evidence_scope)
                    & slopes["flora_scope"].eq(comparison)
                    & slopes["axis"].eq(axis)
                    & slopes["context"].eq(context_value),
                    ["outcome", "geography_slope_log_odds"],
                ]
                paired = ref.merge(
                    comp,
                    on="outcome",
                    suffixes=("_reference", "_comparison"),
                    validate="one_to_one",
                )
                if len(paired) < 2:
                    continue
                x = paired["geography_slope_log_odds_reference"].to_numpy(float)
                y = paired["geography_slope_log_odds_comparison"].to_numpy(float)
                norm = float(np.linalg.norm(x) * np.linalg.norm(y))
                cosine = float(np.dot(x, y) / norm) if norm > 0 else float("nan")
                correlation = (
                    float(np.corrcoef(x, y)[0, 1])
                    if len(paired) >= 3
                    else float("nan")
                )
                sign_concordance = float(np.mean(np.sign(x) == np.sign(y)))
                rows.append(
                    {
                        "evidence_scope": evidence_scope,
                        "axis": axis,
                        "context": context_value,
                        "reference": reference,
                        "comparison": comparison,
                        "n_common_states": int(len(paired)),
                        "cosine_similarity": cosine,
                        "pearson_correlation": correlation,
                        "sign_concordance": sign_concordance,
                    }
                )
    return pd.DataFrame(rows)


def run_three_axis_analysis(
    species_axis: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    ontology: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    scopes = flora_scopes(status_flora, wcvp_ranges, island_tdwg)
    all_slopes: list[pd.DataFrame] = []
    all_omnibus: list[dict[str, Any]] = []
    all_support: list[pd.DataFrame] = []
    all_cell_audits: list[pd.DataFrame] = []
    all_flora_coverage: list[pd.DataFrame] = []
    count_parts: list[pd.DataFrame] = []
    status_slopes: list[pd.DataFrame] = []
    status_omnibus: list[dict[str, Any]] = []
    status_support: list[pd.DataFrame] = []

    model = config["model"]
    fit_cfg = {
        "geography_column": str(model["geography_column"]),
        "context_column": str(model["context_column"]),
        "cluster_column": str(model["cluster_column"]),
        "baseline_covariates": [str(x) for x in model["baseline_covariates"]],
    }

    for evidence_scope in config["evidence_scopes"]:
        ledger, cell_audit = build_valid_state_ledger(
            species_axis,
            ontology,
            config,
            evidence_scope=evidence_scope,
        )
        all_cell_audits.append(cell_audit)
        scope_counts: list[pd.DataFrame] = []

        for flora_scope in config["flora_scopes"]:
            flora = scopes[flora_scope]
            counts, coverage = build_axis_state_counts(
                flora,
                ledger,
                stratum=flora_scope,
            )
            coverage.insert(0, "evidence_scope", evidence_scope)
            all_flora_coverage.append(coverage)
            if counts.empty:
                continue
            counts.insert(0, "evidence_scope", evidence_scope)
            count_parts.append(counts)
            scope_counts.append(counts)

            prepared = _prepare(
                counts.drop(columns="evidence_scope"),
                covariates,
                fit_cfg,
            )
            for axis in config["axes"]:
                for context_value in model["contexts"]:
                    slopes, omnibus, support = _axis_fit(
                        prepared,
                        counts,
                        config,
                        evidence_scope=evidence_scope,
                        flora_scope=flora_scope,
                        axis=str(axis),
                        context_value=str(context_value),
                    )
                    if not slopes.empty:
                        all_slopes.append(slopes)
                    all_omnibus.append(omnibus)
                    if not support.empty:
                        all_support.append(support)

        if scope_counts:
            combined_counts = pd.concat(scope_counts, ignore_index=True)
            contrast_names = {
                str(config["status_contrast"]["reference"]),
                str(config["status_contrast"]["comparison"]),
            }
            contrast_counts = combined_counts.loc[
                combined_counts["stratum"].isin(contrast_names)
            ].copy()
            prepared_contrast = _prepare(
                contrast_counts.drop(columns="evidence_scope"),
                covariates,
                fit_cfg,
            )
            for axis in config["axes"]:
                for context_value in model["contexts"]:
                    slopes, omnibus, support = _fit_status_contrast(
                        prepared_contrast,
                        contrast_counts,
                        config,
                        evidence_scope=evidence_scope,
                        axis=str(axis),
                        context_value=str(context_value),
                    )
                    if not slopes.empty:
                        status_slopes.append(slopes)
                    status_omnibus.append(omnibus)
                    if not support.empty:
                        status_support.append(support)

    omnibus = pd.DataFrame(all_omnibus)
    if not omnibus.empty and "p_value" in omnibus.columns:
        omnibus["q_value"] = (
            omnibus.groupby(
                ["evidence_scope", "flora_scope"],
                group_keys=False,
            )["p_value"]
            .transform(_bh)
        )
        omnibus["axis_supported"] = (
            omnibus["q_value"].le(float(model["alpha"])).fillna(False)
        )

    status_omnibus_frame = pd.DataFrame(status_omnibus)
    if not status_omnibus_frame.empty and "p_value" in status_omnibus_frame.columns:
        status_omnibus_frame["q_value"] = (
            status_omnibus_frame.groupby(
                ["evidence_scope"],
                group_keys=False,
            )["p_value"]
            .transform(_bh)
        )
        status_omnibus_frame["status_vectors_differ"] = (
            status_omnibus_frame["q_value"]
            .le(float(model["alpha"]))
            .fillna(False)
        )

    axis_slopes = (
        pd.concat(all_slopes, ignore_index=True)
        if all_slopes
        else pd.DataFrame()
    )
    similarity = (
        _status_vector_similarity(axis_slopes, config)
        if not axis_slopes.empty
        else pd.DataFrame()
    )

    return {
        "axis_slopes": axis_slopes,
        "axis_omnibus": omnibus,
        "state_support": (
            pd.concat(all_support, ignore_index=True)
            if all_support
            else pd.DataFrame()
        ),
        "status_contrast_slopes": (
            pd.concat(status_slopes, ignore_index=True)
            if status_slopes
            else pd.DataFrame()
        ),
        "status_contrast_omnibus": status_omnibus_frame,
        "status_contrast_support": (
            pd.concat(status_support, ignore_index=True)
            if status_support
            else pd.DataFrame()
        ),
        "status_vector_similarity": similarity,
        "cell_coverage_audit": pd.concat(all_cell_audits, ignore_index=True),
        "flora_coverage": pd.concat(all_flora_coverage, ignore_index=True),
        "counts": (
            pd.concat(count_parts, ignore_index=True)
            if count_parts
            else pd.DataFrame()
        ),
    }


@app.command("run")
def run(
    species_axis_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    wcvp_ranges_csv: Path = typer.Option(..., exists=True),
    island_tdwg_csv: Path = typer.Option(..., exists=True),
    ontology_path: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    ontology = yaml.safe_load(ontology_path.read_text(encoding="utf-8"))
    out = run_three_axis_analysis(
        pd.read_csv(species_axis_csv, dtype=str).fillna(""),
        pd.read_csv(status_flora_csv, dtype=str).fillna(""),
        pd.read_csv(covariates_csv),
        pd.read_csv(wcvp_ranges_csv, dtype=str).fillna(""),
        pd.read_csv(island_tdwg_csv, dtype=str).fillna(""),
        ontology,
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    out["axis_slopes"].to_csv(output_dir / "axis_state_slopes.csv", index=False)
    out["axis_omnibus"].to_csv(output_dir / "axis_omnibus.csv", index=False)
    out["state_support"].to_csv(output_dir / "state_support.csv", index=False)
    out["cell_coverage_audit"].to_csv(
        output_dir / "cell_coverage_audit.csv", index=False
    )
    out["flora_coverage"].to_csv(output_dir / "flora_coverage.csv", index=False)
    out["counts"].to_csv(
        output_dir / "axis_state_counts.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    manifest = {
        "contract": config["contract"],
        "analysis": "three original species-by-axis domains using ontology-valid raw states",
        "resolved_species_axis_cells": int(
            out["cell_coverage_audit"]
            .loc[
                out["cell_coverage_audit"]["evidence_scope"].eq(
                    "all_analysis_eligible"
                ),
                "resolved_axis_cells",
            ]
            .sum()
        ),
        "ontology_valid_species_axis_cells": int(
            out["cell_coverage_audit"]
            .loc[
                out["cell_coverage_audit"]["evidence_scope"].eq(
                    "all_analysis_eligible"
                ),
                "ontology_valid_axis_cells",
            ]
            .sum()
        ),
        "claim_ceiling": config["claim_ceiling"],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
