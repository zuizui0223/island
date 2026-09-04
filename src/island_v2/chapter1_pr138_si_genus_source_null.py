"""Genus-fixed source-pool null for the PR138 Palearctic SI-only floral response.

The PR138 exact-SI analysis asks whether the source-adjusted attraction/accessibility
response persists when measured selfing capacity is unavailable. This module adds the
next decomposition layer: can that SI-restricted distance response be reproduced by
between-genus turnover alone?

For each of the two frozen attraction axes (``large_bee_like`` and
``generalized_accessible``), species scores are permuted only among scored SI species
within the same genus. Island occurrence, floristic status, score availability, GIFT
mainland occurrence, source-region eligibility, source assignments, and all geographic /
climatic covariates remain fixed. Island and mainland source scores are recomputed for
every draw before source-pool subtraction and fitting the Palearctic distance response.

The null is a lineage-sensitivity decomposition, not a causal evolutionary null. Failure
to exceed it means that between-genus composition is sufficient to reproduce the
observed assemblage slope; it does not show that ecological filtering of whole lineages
is absent.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import unicodedata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

AXES = ("large_bee_like", "generalized_accessible")
STRATA = ("all_native", "native_nonendemic")
EVIDENCE_MODES = ("all_analysis_eligible", "direct_only")


def _normalise_name(value: object) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).replace("×", " ")
    return " ".join(text.split()).casefold()


def _binomial_key(value: object) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).replace("×", " ")
    tokens = re.findall(r"[A-Za-z][A-Za-z.-]*", text)
    if len(tokens) < 2:
        return ""
    return f"{tokens[0].casefold()} {tokens[1].casefold()}"


def _genus(value: object) -> str:
    tokens = re.findall(r"[A-Za-z][A-Za-z.-]*", str(value or ""))
    return tokens[0] if len(tokens) >= 2 else ""


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def exact_si_species(trait_ledger: pd.DataFrame) -> set[str]:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(trait_ledger.columns)
    if missing:
        raise ValueError(f"trait ledger missing columns: {sorted(missing)}")
    frame = trait_ledger.copy()
    frame["accepted_species"] = frame["accepted_species"].fillna("").astype(str)
    frame["trait_name"] = frame["trait_name"].fillna("").astype(str)
    frame["normalized_value"] = frame["normalized_value"].fillna("").astype(str)
    species = set(
        frame.loc[
            frame["trait_name"].eq("self_incompatibility")
            & frame["normalized_value"].eq("SI"),
            "accepted_species",
        ]
    )
    species.discard("")
    return species


def axis_species_scores(
    species_scores: pd.DataFrame,
    si_species: set[str],
    axis: str,
) -> pd.DataFrame:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    missing = required - set(species_scores.columns)
    if missing:
        raise ValueError(f"species syndrome scores missing columns: {sorted(missing)}")
    frame = species_scores.loc[
        species_scores["syndrome"].astype(str).eq(axis)
        & species_scores["accepted_species"].astype(str).isin(si_species),
        ["accepted_species", "syndrome_concordance"],
    ].copy()
    frame["accepted_species"] = frame["accepted_species"].astype(str)
    frame["score"] = pd.to_numeric(frame["syndrome_concordance"], errors="coerce")
    frame = frame.dropna(subset=["score"])
    # A species is one frozen syndrome score. Contradictory duplicates fail closed.
    collapsed = frame.groupby("accepted_species", as_index=False).agg(
        n_scores=("score", "nunique"), score=("score", "first")
    )
    collapsed = collapsed.loc[collapsed["n_scores"].eq(1), ["accepted_species", "score"]]
    collapsed["genus"] = collapsed["accepted_species"].map(_genus)
    return (
        collapsed.loc[collapsed["genus"].ne("")]
        .sort_values(["genus", "accepted_species"])
        .reset_index(drop=True)
    )


def draw_within_genus_scores(
    scored_species: pd.DataFrame,
    *,
    n_draws: int,
    rng: np.random.Generator,
) -> np.ndarray:
    observed = scored_species["score"].to_numpy(float)
    draws = np.repeat(observed[:, None], n_draws, axis=1)
    for _, group in scored_species.groupby("genus", sort=True):
        idx = group.index.to_numpy(int)
        if len(idx) < 2:
            continue
        order = np.argsort(rng.random((len(idx), n_draws)), axis=0)
        draws[idx, :] = observed[idx][order]
    return draws


def match_gift_species(
    gift_flora: pd.DataFrame,
    scored_species: pd.DataFrame,
) -> pd.DataFrame:
    required = {"entity_ID", "work_species"}
    missing = required - set(gift_flora.columns)
    if missing:
        raise ValueError(f"GIFT flora missing columns: {sorted(missing)}")
    species = scored_species[["accepted_species"]].drop_duplicates().copy()
    species["norm"] = species["accepted_species"].map(_normalise_name)
    species["binomial"] = species["accepted_species"].map(_binomial_key)
    exact_groups = species.groupby("norm")["accepted_species"].agg(lambda x: sorted(set(x)))
    exact = {key: vals[0] for key, vals in exact_groups.items() if key and len(vals) == 1}
    bin_groups = species.groupby("binomial")["accepted_species"].agg(lambda x: sorted(set(x)))
    binomial = {key: vals[0] for key, vals in bin_groups.items() if key and len(vals) == 1}

    # Match each distinct GIFT name once, then join back to entity memberships.
    # The GIFT table contains ~1.7M entity-species rows but far fewer unique names.
    names = gift_flora[["work_species"]].drop_duplicates().copy()
    names["norm"] = names["work_species"].map(_normalise_name)
    names["binomial"] = names["work_species"].map(_binomial_key)
    names["accepted_species"] = names["norm"].map(exact).fillna("")
    missing_match = names["accepted_species"].eq("") & names["binomial"].ne("")
    names.loc[missing_match, "accepted_species"] = (
        names.loc[missing_match, "binomial"].map(binomial).fillna("")
    )
    matched_names = names.loc[
        names["accepted_species"].ne(""), ["work_species", "accepted_species"]
    ]
    flora = gift_flora[["entity_ID", "work_species"]].merge(
        matched_names, on="work_species", how="inner", validate="many_to_one"
    )
    return flora.drop_duplicates(["entity_ID", "accepted_species"])


def _stratum_mask(frame: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all_native":
        return frame["origin_status"].astype(str).eq("native")
    if stratum == "native_nonendemic":
        return frame["floristic_status"].astype(str).eq("native_nonendemic")
    raise ValueError(f"unsupported stratum: {stratum}")


def _aggregate_memberships(
    memberships: pd.DataFrame,
    scored_species: pd.DataFrame,
    draws: np.ndarray,
    *,
    group_column: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Aggregate group x species memberships with a sparse incidence matrix."""
    species_to_index = pd.Series(
        np.arange(len(scored_species), dtype=int),
        index=scored_species["accepted_species"].to_numpy(),
    )
    work = memberships[[group_column, "accepted_species"]].drop_duplicates().copy()
    work["species_index"] = work["accepted_species"].map(species_to_index)
    work = work.dropna(subset=["species_index"]).copy()
    if work.empty:
        return (
            pd.DataFrame(columns=[group_column, "observed_score", "n_species"]),
            np.empty((0, draws.shape[1])),
        )
    work["species_index"] = work["species_index"].astype(int)
    group_codes, groups = pd.factorize(work[group_column], sort=True)
    incidence = sparse.csr_matrix(
        (
            np.ones(len(work), dtype=np.float64),
            (group_codes, work["species_index"].to_numpy(int)),
        ),
        shape=(len(groups), len(scored_species)),
    )
    counts = np.asarray(incidence.sum(axis=1)).ravel()
    observed = np.asarray(incidence @ scored_species["score"].to_numpy(float)).ravel() / counts
    null = np.asarray(incidence @ draws) / counts[:, None]
    summary = pd.DataFrame(
        {
            group_column: groups.astype(object),
            "observed_score": observed,
            "n_species": counts.astype(int),
        }
    )
    return summary, null


def _source_expectation_matrices(
    assignments: pd.DataFrame,
    region_scores: pd.DataFrame,
    region_null: np.ndarray,
    *,
    min_region_species: int,
    min_source_regions: int,
) -> dict[str, tuple[pd.DataFrame, np.ndarray]]:
    eligible_mask = region_scores["n_species"].ge(min_region_species).to_numpy(bool)
    if not eligible_mask.any():
        return {}
    eligible_ids = set(region_scores.loc[eligible_mask, "entity_ID"].astype(int))
    entity_to_row = pd.Series(
        np.arange(len(region_scores), dtype=int),
        index=region_scores["entity_ID"].astype(int).to_numpy(),
    )
    out: dict[str, tuple[pd.DataFrame, np.ndarray]] = {}
    for mode, mode_assignments in assignments.groupby("source_mode", sort=True):
        work = mode_assignments.loc[
            mode_assignments["entity_ID"].astype(int).isin(eligible_ids),
            ["island_id", "entity_ID"],
        ].drop_duplicates().copy()
        if work.empty:
            continue
        island_codes, island_ids = pd.factorize(work["island_id"].astype(str), sort=True)
        region_rows = work["entity_ID"].astype(int).map(entity_to_row).to_numpy(int)
        incidence = sparse.csr_matrix(
            (np.ones(len(work), dtype=np.float64), (island_codes, region_rows)),
            shape=(len(island_ids), len(region_scores)),
        )
        counts = np.asarray(incidence.sum(axis=1)).ravel()
        eligible_islands = counts >= int(min_source_regions)
        if not eligible_islands.any():
            continue
        obs = np.asarray(incidence @ region_scores["observed_score"].to_numpy(float)).ravel()
        obs = obs / counts
        null = np.asarray(incidence @ region_null) / counts[:, None]
        rows = pd.DataFrame(
            {
                "island_id": island_ids.astype(object),
                "source_expectation": obs,
                "n_source_regions": counts.astype(int),
            }
        ).loc[eligible_islands].reset_index(drop=True)
        out[str(mode)] = (rows, null[eligible_islands, :])
    return out


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _design(
    work: pd.DataFrame, geography: str, baseline: list[str]
) -> tuple[np.ndarray, list[str]]:
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    X = np.column_stack(
        [
            np.ones(len(work), dtype=float),
            _standardize(work[geography]),
            *[_standardize(work[x]) for x in baseline],
        ]
    )
    return X, names


def _fit_clustered(
    y: np.ndarray,
    X: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    n, p = X.shape
    unique = np.unique(clusters.astype(str))
    if n < max(10, p + 3) or len(unique) < 2:
        return pd.DataFrame(), {
            "status": "insufficient_complete_rows",
            "n_rows": n,
            "n_clusters": len(unique),
        }
    bread = np.linalg.pinv(X.T @ X)
    beta = bread @ (X.T @ y)
    residual = y - X @ beta
    meat = np.zeros((p, p), dtype=float)
    cluster_labels = clusters.astype(str)
    for cluster in unique:
        mask = cluster_labels == cluster
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if len(unique) > 1 and n > p:
        covariance *= (len(unique) / (len(unique) - 1.0)) * ((n - 1.0) / (n - p))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    rows = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        p_value = math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")
        rows.append(
            {
                "predictor": name,
                "estimate": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": p_value,
            }
        )
    return pd.DataFrame(rows), {
        "status": "fit",
        "n_rows": n,
        "n_clusters": int(len(unique)),
    }


def _prepare_axis(
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    scored_species: pd.DataFrame,
    draws: np.ndarray,
    palearctic_ids: set[str],
    *,
    min_region_species: int,
    min_source_regions: int,
) -> dict[str, Any]:
    gift_matches = match_gift_species(gift_flora, scored_species)
    region_scores, region_null = _aggregate_memberships(
        gift_matches.rename(columns={"entity_ID": "group"}),
        scored_species,
        draws,
        group_column="group",
    )
    region_scores = region_scores.rename(columns={"group": "entity_ID"})
    region_scores["entity_ID"] = region_scores["entity_ID"].astype(int)
    source = _source_expectation_matrices(
        assignments,
        region_scores,
        region_null,
        min_region_species=min_region_species,
        min_source_regions=min_source_regions,
    )

    island: dict[str, tuple[pd.DataFrame, np.ndarray]] = {}
    flora = status_flora.loc[status_flora["island_id"].astype(str).isin(palearctic_ids)].copy()
    for stratum in STRATA:
        memberships = flora.loc[_stratum_mask(flora, stratum), ["island_id", "accepted_species"]]
        scores, null = _aggregate_memberships(
            memberships,
            scored_species,
            draws,
            group_column="island_id",
        )
        island[stratum] = (scores, null)

    genus_sizes = scored_species.groupby("genus")["accepted_species"].nunique()
    permutable_genera = set(genus_sizes.loc[genus_sizes.ge(2)].index)
    audit = {
        "n_scored_si_species": int(len(scored_species)),
        "n_genera": int(scored_species["genus"].nunique()),
        "n_permutable_genera": int(len(permutable_genera)),
        "n_permutable_species": int(scored_species["genus"].isin(permutable_genera).sum()),
        "n_eligible_source_regions": int(region_scores["n_species"].ge(min_region_species).sum()),
    }
    return {"island": island, "source": source, "audit": audit}


def _combine_axis_scenario(
    axis_payloads: dict[str, dict[str, Any]],
    *,
    stratum: str,
    source_mode: str,
    covariates: pd.DataFrame,
    geography: str,
    baseline: list[str],
    cluster: str,
) -> dict[str, Any] | None:
    assembled: dict[str, pd.DataFrame] = {}
    null_by_axis: dict[str, np.ndarray] = {}
    source_null_by_axis: dict[str, np.ndarray] = {}
    for axis in AXES:
        payload = axis_payloads[axis]
        if source_mode not in payload["source"]:
            return None
        island_scores, island_null = payload["island"][stratum]
        source_scores, source_null = payload["source"][source_mode]
        island = island_scores.reset_index(drop=True).copy()
        source = source_scores.reset_index(drop=True).copy()
        island["island_row"] = np.arange(len(island), dtype=int)
        source["source_row"] = np.arange(len(source), dtype=int)
        merged = island.merge(source, on="island_id", how="inner", validate="one_to_one")
        if merged.empty:
            return None
        merged[f"observed_adjusted_{axis}"] = (
            merged["observed_score"] - merged["source_expectation"]
        )
        assembled[axis] = merged[
            ["island_id", "island_row", "source_row", f"observed_adjusted_{axis}"]
        ]
        null_by_axis[axis] = island_null
        source_null_by_axis[axis] = source_null

    joint = assembled[AXES[0]].merge(
        assembled[AXES[1]], on="island_id", how="inner", validate="one_to_one"
    )
    cov_cols = ["island_id", geography, cluster, *baseline]
    joint = joint.merge(
        covariates[cov_cols].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    needed = [geography, cluster, *baseline]
    for column in [geography, *baseline]:
        joint[column] = pd.to_numeric(joint[column], errors="coerce")
    joint[cluster] = joint[cluster].fillna("").astype(str)
    joint = joint.dropna(subset=needed).loc[lambda x: x[cluster].ne("")].reset_index(drop=True)
    if joint.empty:
        return None

    observed = (
        -joint[f"observed_adjusted_{AXES[0]}"].to_numpy(float)
        + joint[f"observed_adjusted_{AXES[1]}"].to_numpy(float)
    ) / 2.0

    null_components: dict[str, np.ndarray] = {}
    for axis in AXES:
        frame = assembled[axis]
        index = joint.set_index("island_id")
        lookup = frame.set_index("island_id")
        island_rows = lookup.loc[index.index, "island_row"].to_numpy(int)
        source_rows = lookup.loc[index.index, "source_row"].to_numpy(int)
        source_null = source_null_by_axis[axis]
        null_components[axis] = (
            null_by_axis[axis][island_rows, :] - source_null[source_rows, :]
        )
    null_response = (-null_components[AXES[0]] + null_components[AXES[1]]) / 2.0
    X, names = _design(joint, geography, baseline)
    coefficients, fit = _fit_clustered(
        observed,
        X,
        names,
        joint[cluster].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return {"status": fit["status"], "n_unique_islands": int(len(joint))}
    distance_name = f"z_{geography}"
    distance = coefficients.set_index("predictor").loc[distance_name]
    beta_null = np.linalg.pinv(X.T @ X) @ (X.T @ null_response)
    distance_index = names.index(distance_name)
    null_slopes = beta_null[distance_index, :]
    null_mean_response = null_response.mean(axis=1)
    residual = observed - null_mean_response
    residual_coefficients, residual_fit = _fit_clustered(
        residual,
        X,
        names,
        joint[cluster].astype(str).to_numpy(),
    )
    residual_distance = residual_coefficients.set_index("predictor").loc[distance_name]
    observed_slope = float(distance["estimate"])
    empirical_p = float((1 + np.sum(null_slopes >= observed_slope)) / (1 + len(null_slopes)))
    return {
        "status": "fit",
        "n_unique_islands": int(joint["island_id"].nunique()),
        "n_clusters": int(residual_fit["n_clusters"]),
        "observed_distance_estimate": observed_slope,
        "observed_distance_se": float(distance["cluster_robust_se"]),
        "observed_distance_p": float(distance["p_value"]),
        "null_distance_mean": float(np.mean(null_slopes)),
        "null_distance_q025": float(np.quantile(null_slopes, 0.025)),
        "null_distance_q975": float(np.quantile(null_slopes, 0.975)),
        "residual_distance_estimate": float(residual_distance["estimate"]),
        "residual_distance_se": float(residual_distance["cluster_robust_se"]),
        "residual_distance_p": float(residual_distance["p_value"]),
        "empirical_one_sided_p_observed_gt_null": empirical_p,
    }


def run_evidence_mode(
    trait_ledger: pd.DataFrame,
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
    *,
    evidence_mode: str,
    n_draws: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    si = exact_si_species(trait_ledger)
    realm = realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id")
    palearctic_ids = set(
        realm.loc[
            realm["biogeographic_realm"].astype(str).eq("Palearctic"), "island_id"
        ].astype(str)
    )
    geography = str(pattern_config["geography_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    cluster = str(pattern_config["cluster_column"])
    min_region_species = int(
        source_config["source_region_scores"]["minimum_trait_scored_species_per_region_syndrome"]
    )
    min_source_regions = int(
        source_config["response"]["source_expectation_requires_minimum_source_regions"]
    )
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]

    rng = np.random.default_rng(seed)
    axis_payloads: dict[str, dict[str, Any]] = {}
    audit_rows: list[dict[str, Any]] = []
    for axis in AXES:
        scored = axis_species_scores(species_scores, si, axis)
        draws = draw_within_genus_scores(scored, n_draws=n_draws, rng=rng)
        payload = _prepare_axis(
            status_flora,
            gift_flora,
            assignments,
            scored,
            draws,
            palearctic_ids,
            min_region_species=min_region_species,
            min_source_regions=min_source_regions,
        )
        axis_payloads[axis] = payload
        audit_rows.append(
            {
                "evidence_mode": evidence_mode,
                "axis": axis,
                "n_exact_si_species": int(len(si)),
                **payload["audit"],
            }
        )

    rows: list[dict[str, Any]] = []
    for stratum in STRATA:
        for source_mode in source_modes:
            result = _combine_axis_scenario(
                axis_payloads,
                stratum=stratum,
                source_mode=source_mode,
                covariates=covariates,
                geography=geography,
                baseline=baseline,
                cluster=cluster,
            )
            if result is None:
                rows.append(
                    {
                        "evidence_mode": evidence_mode,
                        "stratum": stratum,
                        "source_mode": source_mode,
                        "status": "not_testable",
                    }
                )
            else:
                rows.append(
                    {
                        "evidence_mode": evidence_mode,
                        "stratum": stratum,
                        "source_mode": source_mode,
                        **result,
                    }
                )
    results = pd.DataFrame(rows)
    fit = results["status"].eq("fit")
    results["residual_distance_q"] = np.nan
    results["empirical_one_sided_q"] = np.nan
    if fit.any():
        results.loc[fit, "residual_distance_q"] = (
            results.loc[fit]
            .groupby(["evidence_mode", "stratum"], group_keys=False)["residual_distance_p"]
            .apply(_bh)
        )
        results.loc[fit, "empirical_one_sided_q"] = (
            results.loc[fit]
            .groupby(["evidence_mode", "stratum"], group_keys=False)[
                "empirical_one_sided_p_observed_gt_null"
            ]
            .apply(_bh)
        )
    return results, pd.DataFrame(audit_rows)


def run_si_genus_source_null(
    all_ledger: pd.DataFrame,
    all_scores: pd.DataFrame,
    direct_ledger: pd.DataFrame,
    direct_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
    *,
    n_draws: int = 1000,
    seed: int = 20260827,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    if n_draws < 100:
        raise ValueError("n_draws must be at least 100")
    pairs = {
        "all_analysis_eligible": (all_ledger, all_scores),
        "direct_only": (direct_ledger, direct_scores),
    }
    result_parts: list[pd.DataFrame] = []
    audit_parts: list[pd.DataFrame] = []
    for offset, mode in enumerate(EVIDENCE_MODES):
        ledger, scores = pairs[mode]
        results, audit = run_evidence_mode(
            ledger,
            scores,
            status_flora,
            gift_flora,
            assignments,
            covariates,
            realm_assignment,
            pattern_config,
            source_config,
            evidence_mode=mode,
            n_draws=n_draws,
            seed=seed + offset * 100_003,
        )
        result_parts.append(results)
        audit_parts.append(audit)
    manifest = {
        "contract": "chapter1_pr138_si_genus_source_null_v1",
        "target_population": "Palearctic exact-SI island floras",
        "response": (
            "source-adjusted attraction_shift = "
            "(-large_bee_like + generalized_accessible) / 2"
        ),
        "evidence_modes": list(EVIDENCE_MODES),
        "n_draws": int(n_draws),
        "seed": int(seed),
        "null": "axis-specific syndrome scores permuted among scored exact-SI species within genus",
        "island_species_membership_fixed": True,
        "gift_mainland_species_membership_fixed": True,
        "score_availability_fixed": True,
        "source_region_eligibility_fixed": True,
        "source_assignments_fixed_and_outcome_blind": True,
        "source_pool_recomputed_each_draw": True,
        "causal_pollinator_claimed": False,
        "claim_ceiling": (
            "A residual slope above this null would indicate response structure beyond measured "
            "between-genus composition. Failure to exceed the null means genus composition is "
            "sufficient to reproduce the slope; it does not rule out ecological filtering of "
            "whole genera or identify a pollination mechanism."
        ),
    }
    return {
        "results": pd.concat(result_parts, ignore_index=True),
        "audit": pd.concat(audit_parts, ignore_index=True),
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-ledger-csv", type=Path, required=True)
    parser.add_argument("--all-scores-csv", type=Path, required=True)
    parser.add_argument("--direct-ledger-csv", type=Path, required=True)
    parser.add_argument("--direct-scores-csv", type=Path, required=True)
    parser.add_argument("--status-flora-csv", type=Path, required=True)
    parser.add_argument("--gift-flora-csv", type=Path, required=True)
    parser.add_argument("--source-assignments-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--source-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()

    pattern = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    source = yaml.safe_load(args.source_config_path.read_text(encoding="utf-8"))
    outputs = run_si_genus_source_null(
        pd.read_csv(args.all_ledger_csv),
        pd.read_csv(args.all_scores_csv),
        pd.read_csv(args.direct_ledger_csv),
        pd.read_csv(args.direct_scores_csv),
        pd.read_csv(args.status_flora_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        pattern,
        source,
        n_draws=args.n_draws,
        seed=args.seed,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs["results"].to_csv(args.output_dir / "si_genus_source_null_results.csv", index=False)
    outputs["audit"].to_csv(args.output_dir / "si_genus_source_null_audit.csv", index=False)
    (args.output_dir / "si_genus_source_null_manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outputs["manifest"], indent=2))


if __name__ == "__main__":
    main()
