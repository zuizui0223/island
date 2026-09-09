"""Leave-one-genus influence audit for the PR138 Palearctic exact-SI response.

The genus-fixed source-pool null shows that the exact-SI Palearctic attraction/accessibility
slope does not exceed between-genus expectations. This module asks whether that assembly
signal is itself driven by one or a few genera.

Each scored SI genus is omitted in turn. Island floras, GIFT mainland memberships, source
assignments, floristic status, geography, climate and spatial blocks remain fixed. Island
and mainland syndrome scores are recomputed after each omission before source-pool
subtraction and refitting the same Palearctic distance model.

This is an influence diagnostic, not a causal genus-selection model.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from island_v2.chapter1_pr138_si_genus_source_null import (
    AXES,
    STRATA,
    _bh,
    _fit_clustered,
    axis_species_scores,
    exact_si_species,
    match_gift_species,
)

FULL_SENTINEL = "__full__"
EVIDENCE_MODES = ("all_analysis_eligible", "direct_only")


def _standardize(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _group_genus_arrays(
    memberships: pd.DataFrame,
    scored_species: pd.DataFrame,
    *,
    group_column: str,
    group_index: dict[object, int],
    genus_index: dict[str, int],
) -> tuple[np.ndarray, np.ndarray]:
    work = memberships[[group_column, "accepted_species"]].drop_duplicates().merge(
        scored_species[["accepted_species", "genus", "score"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    work = work.loc[work[group_column].isin(group_index)].copy()
    count = np.zeros((len(group_index), len(genus_index)), dtype=np.float32)
    total = np.zeros((len(group_index), len(genus_index)), dtype=np.float64)
    if work.empty:
        return count, total
    rows = work[group_column].map(group_index).to_numpy(int)
    cols = work["genus"].map(genus_index).to_numpy(int)
    scores = work["score"].to_numpy(float)
    np.add.at(count, (rows, cols), 1.0)
    np.add.at(total, (rows, cols), scores)
    return count, total


def _scores_after_omission(
    count: np.ndarray,
    total: np.ndarray,
    genus_position: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    n = count.sum(axis=1, dtype=np.float64)
    s = total.sum(axis=1)
    if genus_position is not None:
        n = n - count[:, genus_position]
        s = s - total[:, genus_position]
    score = np.full(len(n), np.nan, dtype=float)
    valid = n > 0
    score[valid] = s[valid] / n[valid]
    return n, score


def _prepare_axis(
    scored_species: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    *,
    island_index: dict[str, int],
    region_index: dict[int, int],
    genus_index: dict[str, int],
) -> dict[str, Any]:
    gift_matches = match_gift_species(gift_flora, scored_species).copy()
    gift_matches["entity_ID"] = gift_matches["entity_ID"].astype(int)
    region_count, region_sum = _group_genus_arrays(
        gift_matches,
        scored_species,
        group_column="entity_ID",
        group_index=region_index,
        genus_index=genus_index,
    )
    island: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for stratum in STRATA:
        if stratum == "all_native":
            mask = status_flora["origin_status"].astype(str).eq("native")
        elif stratum == "native_nonendemic":
            mask = status_flora["floristic_status"].astype(str).eq("native_nonendemic")
        else:
            raise ValueError(f"unsupported stratum: {stratum}")
        island[stratum] = _group_genus_arrays(
            status_flora.loc[mask, ["island_id", "accepted_species"]],
            scored_species,
            group_column="island_id",
            group_index=island_index,
            genus_index=genus_index,
        )
    return {"region": (region_count, region_sum), "island": island}


def _source_matrices(
    source_assignments: pd.DataFrame,
    *,
    island_index: dict[str, int],
    region_index: dict[int, int],
    source_modes: list[str],
) -> dict[str, sparse.csr_matrix]:
    matrices: dict[str, sparse.csr_matrix] = {}
    for mode in source_modes:
        work = source_assignments.loc[
            source_assignments["source_mode"].astype(str).eq(mode),
            ["island_id", "entity_ID"],
        ].drop_duplicates()
        work = work.loc[
            work["island_id"].astype(str).isin(island_index)
            & work["entity_ID"].astype(int).isin(region_index)
        ].copy()
        rows = work["island_id"].astype(str).map(island_index).to_numpy(int)
        cols = work["entity_ID"].astype(int).map(region_index).to_numpy(int)
        matrices[mode] = sparse.csr_matrix(
            (np.ones(len(work), dtype=float), (rows, cols)),
            shape=(len(island_index), len(region_index)),
        )
    return matrices


def _source_expectation(
    matrix: sparse.csr_matrix,
    region_count: np.ndarray,
    region_score: np.ndarray,
    *,
    min_region_species: int,
    min_source_regions: int,
) -> np.ndarray:
    eligible = region_count >= int(min_region_species)
    n_sources = np.asarray(matrix @ eligible.astype(float)).ravel()
    values = np.nan_to_num(region_score, nan=0.0) * eligible
    total = np.asarray(matrix @ values).ravel()
    out = np.full(matrix.shape[0], np.nan, dtype=float)
    valid = n_sources >= int(min_source_regions)
    out[valid] = total[valid] / n_sources[valid]
    return out


def _fit_distance(
    response: np.ndarray,
    covariates: pd.DataFrame,
    *,
    geography: str,
    baseline: list[str],
    cluster: str,
) -> dict[str, Any]:
    arrays = [covariates[geography].to_numpy(float)] + [
        covariates[column].to_numpy(float) for column in baseline
    ]
    valid = np.isfinite(response)
    for values in arrays:
        valid &= np.isfinite(values)
    n = int(valid.sum())
    if n < 50:
        return {"status": "not_testable", "n_unique_islands": n}
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    design = np.column_stack(
        [np.ones(n, dtype=float), *[_standardize(values[valid]) for values in arrays]]
    )
    coefficients, fit = _fit_clustered(
        response[valid],
        design,
        names,
        covariates.loc[valid, cluster].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return {"status": str(fit.get("status", "fit_failed")), "n_unique_islands": n}
    distance = coefficients.set_index("predictor").loc[f"z_{geography}"]
    return {
        "status": "fit",
        "n_unique_islands": n,
        "n_clusters": int(fit["n_clusters"]),
        "distance_estimate": float(distance["estimate"]),
        "distance_se": float(distance["cluster_robust_se"]),
        "distance_p": float(distance["p_value"]),
    }


def run_evidence_mode(
    trait_ledger: pd.DataFrame,
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    source_assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
    *,
    evidence_mode: str,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    if evidence_mode not in EVIDENCE_MODES:
        raise ValueError(f"unsupported evidence mode: {evidence_mode}")
    si_species = exact_si_species(trait_ledger)
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

    palearctic_ids = set(
        realm_assignment.loc[
            realm_assignment["biogeographic_realm"].astype(str).eq("Palearctic"),
            "island_id",
        ].astype(str)
    )
    cov = covariates.loc[
        covariates["island_id"].astype(str).isin(palearctic_ids),
        ["island_id", geography, cluster, *baseline],
    ].copy()
    for column in [geography, *baseline]:
        cov[column] = pd.to_numeric(cov[column], errors="coerce")
    cov[cluster] = cov[cluster].fillna("").astype(str)
    cov = (
        cov.dropna(subset=[geography, *baseline])
        .loc[lambda x: x[cluster].ne("")]
        .drop_duplicates("island_id")
        .sort_values("island_id")
        .reset_index(drop=True)
    )
    island_ids = cov["island_id"].astype(str).tolist()
    island_index = {value: idx for idx, value in enumerate(island_ids)}
    flora = status_flora.loc[
        status_flora["island_id"].astype(str).isin(island_index)
    ].copy()
    assignments = source_assignments.loc[
        source_assignments["island_id"].astype(str).isin(island_index)
    ].copy()
    assignments["entity_ID"] = assignments["entity_ID"].astype(int)
    region_ids = sorted(assignments["entity_ID"].unique())
    region_index = {value: idx for idx, value in enumerate(region_ids)}

    scored_by_axis = {
        axis: axis_species_scores(species_scores, si_species, axis) for axis in AXES
    }
    genera = sorted(set().union(*[set(frame["genus"]) for frame in scored_by_axis.values()]))
    genus_index = {value: idx for idx, value in enumerate(genera)}
    payload = {
        axis: _prepare_axis(
            scored_by_axis[axis],
            flora,
            gift_flora,
            island_index=island_index,
            region_index=region_index,
            genus_index=genus_index,
        )
        for axis in AXES
    }
    source_matrix = _source_matrices(
        assignments,
        island_index=island_index,
        region_index=region_index,
        source_modes=source_modes,
    )

    rows: list[dict[str, Any]] = []
    omission_specs = [(FULL_SENTINEL, None)] + [
        (genus, genus_index[genus]) for genus in genera
    ]
    for omitted_genus, genus_position in omission_specs:
        adjusted: dict[str, dict[tuple[str, str], np.ndarray]] = {
            axis: {} for axis in AXES
        }
        for axis in AXES:
            region_count, region_score = _scores_after_omission(
                *payload[axis]["region"], genus_position
            )
            source_expectations = {
                source_mode: _source_expectation(
                    source_matrix[source_mode],
                    region_count,
                    region_score,
                    min_region_species=min_region_species,
                    min_source_regions=min_source_regions,
                )
                for source_mode in source_modes
            }
            for stratum in STRATA:
                _, island_score = _scores_after_omission(
                    *payload[axis]["island"][stratum], genus_position
                )
                for source_mode in source_modes:
                    adjusted[axis][(stratum, source_mode)] = (
                        island_score - source_expectations[source_mode]
                    )
        for stratum in STRATA:
            for source_mode in source_modes:
                response = (
                    -adjusted[AXES[0]][(stratum, source_mode)]
                    + adjusted[AXES[1]][(stratum, source_mode)]
                ) / 2.0
                rows.append(
                    {
                        "evidence_mode": evidence_mode,
                        "omitted_genus": omitted_genus,
                        "stratum": stratum,
                        "source_mode": source_mode,
                        **_fit_distance(
                            response,
                            cov,
                            geography=geography,
                            baseline=baseline,
                            cluster=cluster,
                        ),
                    }
                )

    results = pd.DataFrame(rows)
    results["distance_q_across_source_modes"] = np.nan
    fit = results["status"].eq("fit")
    if fit.any():
        results.loc[fit, "distance_q_across_source_modes"] = (
            results.loc[fit]
            .groupby(["omitted_genus", "stratum"], group_keys=False)["distance_p"]
            .apply(_bh)
        )

    full = results.loc[
        results["omitted_genus"].eq(FULL_SENTINEL),
        [
            "stratum",
            "source_mode",
            "distance_estimate",
            "distance_q_across_source_modes",
            "n_unique_islands",
        ],
    ].rename(
        columns={
            "distance_estimate": "full_distance_estimate",
            "distance_q_across_source_modes": "full_distance_q",
            "n_unique_islands": "full_n_unique_islands",
        }
    )
    influence = results.loc[~results["omitted_genus"].eq(FULL_SENTINEL)].merge(
        full,
        on=["stratum", "source_mode"],
        how="left",
        validate="many_to_one",
    )
    influence["deletion_delta_estimate"] = (
        influence["distance_estimate"] - influence["full_distance_estimate"]
    )
    influence["absolute_deletion_delta"] = influence["deletion_delta_estimate"].abs()
    influence["positive_direction"] = influence["distance_estimate"] > 0
    influence["fdr_supported"] = influence["distance_q_across_source_modes"] < 0.05

    ranking = (
        influence.groupby("omitted_genus", as_index=False)
        .agg(
            mean_absolute_delta=("absolute_deletion_delta", "mean"),
            max_absolute_delta=("absolute_deletion_delta", "max"),
            minimum_estimate=("distance_estimate", "min"),
            maximum_estimate=("distance_estimate", "max"),
            n_positive=("positive_direction", "sum"),
            n_fdr_supported=("fdr_supported", "sum"),
            minimum_islands=("n_unique_islands", "min"),
        )
        .sort_values(["mean_absolute_delta", "max_absolute_delta"], ascending=False)
        .reset_index(drop=True)
    )
    ranking["influence_rank"] = np.arange(1, len(ranking) + 1)
    for axis in AXES:
        counts = (
            scored_by_axis[axis]
            .groupby("genus")["accepted_species"]
            .nunique()
            .rename(f"n_species_{axis}")
        )
        ranking = ranking.merge(
            counts,
            left_on="omitted_genus",
            right_index=True,
            how="left",
        )
    ranking = ranking.fillna({f"n_species_{axis}": 0 for axis in AXES})

    manifest = {
        "contract": "chapter1_pr138_si_leave_one_genus_influence_v1",
        "evidence_mode": evidence_mode,
        "target_population": "Palearctic exact-SI island floras",
        "response": "source-adjusted attraction_shift",
        "n_scored_si_genera_tested": int(len(genera)),
        "all_scored_si_genera_omitted_once": True,
        "island_membership_fixed": True,
        "gift_mainland_membership_fixed": True,
        "source_assignments_fixed_and_outcome_blind": True,
        "source_scores_recomputed_after_each_omission": True,
        "causal_genus_effect_claimed": False,
        "claim_ceiling": (
            "Persistence after each omission shows that no single scored SI genus is "
            "necessary for the observed direction. It does not identify which assembly "
            "process generated the between-genus pattern."
        ),
    }
    return {
        "results": results,
        "influence": influence,
        "ranking": ranking,
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
    args = parser.parse_args()

    status_flora = pd.read_csv(args.status_flora_csv)
    gift_flora = pd.read_csv(args.gift_flora_csv)
    assignments = pd.read_csv(args.source_assignments_csv)
    covariates = pd.read_csv(args.covariates_csv)
    realm = pd.read_csv(args.realm_assignment_csv)
    pattern = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    source = yaml.safe_load(args.source_config_path.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    evidence_inputs = {
        "all_analysis_eligible": (
            pd.read_csv(args.all_ledger_csv),
            pd.read_csv(args.all_scores_csv),
        ),
        "direct_only": (
            pd.read_csv(args.direct_ledger_csv),
            pd.read_csv(args.direct_scores_csv),
        ),
    }
    manifests: dict[str, Any] = {}
    for evidence_mode, (ledger, scores) in evidence_inputs.items():
        outputs = run_evidence_mode(
            ledger,
            scores,
            status_flora,
            gift_flora,
            assignments,
            covariates,
            realm,
            pattern,
            source,
            evidence_mode=evidence_mode,
        )
        outputs["results"].to_csv(
            args.output_dir / f"{evidence_mode}_si_leave_one_genus_results.csv", index=False
        )
        outputs["influence"].to_csv(
            args.output_dir / f"{evidence_mode}_si_leave_one_genus_influence.csv", index=False
        )
        outputs["ranking"].to_csv(
            args.output_dir / f"{evidence_mode}_si_leave_one_genus_ranking.csv", index=False
        )
        manifests[evidence_mode] = outputs["manifest"]
    (args.output_dir / "si_leave_one_genus_manifest.json").write_text(
        json.dumps(manifests, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifests, indent=2))


if __name__ == "__main__":
    main()
