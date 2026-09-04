"""Restrict the attraction/accessibility response to non-selfing reproductive subsets.

This is the strongest available decomposition of the two candidate explanations for the
northern floral response using the frozen PR138 data. The attraction/accessibility score is
recomputed only from species classified by the all-evidence reproductive ledger as:

- self-incompatible (SI),
- predominantly outcrossing, or
- both SI and predominantly outcrossing.

The same external GIFT source-pool subtraction and the same geography/area/climate model
are then applied. Ambiguous reproductive states are excluded rather than coerced.
Persistence of the attraction response in SI species rules out measured selfing syndrome
as a necessary explanation, but does not by itself identify pollinator-selection causally.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_pr136_biogeographic_residual import (
    _fit_weighted_clustered_design,
)
from island_v2.chapter1_pr138_source_pool_sensitivity import (
    build_adjusted_island_scores,
    build_mainland_syndrome_scores,
    build_source_expectations,
    match_gift_species_to_frozen_scores,
)
from island_v2.chapter1_pr138_syndrome_analysis import (
    build_island_syndrome_scores,
    score_species_syndromes,
)


def _standardize(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (values - mean) / sd


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
    adjusted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def reproductive_restriction_sets(
    trait_ledger: pd.DataFrame,
) -> dict[str, set[str]]:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(trait_ledger.columns)
    if missing:
        raise ValueError(f"trait ledger missing columns: {sorted(missing)}")

    ledger = trait_ledger.copy()
    ledger["accepted_species"] = ledger["accepted_species"].fillna("").astype(str)
    ledger["trait_name"] = ledger["trait_name"].fillna("").astype(str)
    ledger["normalized_value"] = ledger["normalized_value"].fillna("").astype(str)

    si = set(
        ledger.loc[
            ledger["trait_name"].eq("self_incompatibility")
            & ledger["normalized_value"].eq("SI"),
            "accepted_species",
        ]
    )
    outcrossing = set(
        ledger.loc[
            ledger["trait_name"].eq("mating_system")
            & ledger["normalized_value"].eq("predominantly_outcrossing"),
            "accepted_species",
        ]
    )
    si.discard("")
    outcrossing.discard("")
    return {
        "si_only": si,
        "predominantly_outcrossing": outcrossing,
        "si_and_predominantly_outcrossing": si & outcrossing,
    }


def _attraction_table(adjusted_scores: pd.DataFrame) -> pd.DataFrame:
    required = {
        "source_mode",
        "island_id",
        "stratum",
        "syndrome",
        "syndrome_score",
    }
    missing = required - set(adjusted_scores.columns)
    if missing:
        raise ValueError(f"adjusted scores missing columns: {sorted(missing)}")
    keep = adjusted_scores.loc[
        adjusted_scores["syndrome"].isin(
            ["large_bee_like", "generalized_accessible"]
        ),
        ["source_mode", "island_id", "stratum", "syndrome", "syndrome_score"],
    ].copy()
    pivot = (
        keep.drop_duplicates(
            ["source_mode", "island_id", "stratum", "syndrome"]
        )
        .pivot(
            index=["source_mode", "island_id", "stratum"],
            columns="syndrome",
            values="syndrome_score",
        )
        .reset_index()
    )
    for column in ["large_bee_like", "generalized_accessible"]:
        if column not in pivot.columns:
            pivot[column] = np.nan
    pivot["attraction_shift"] = (
        -pd.to_numeric(pivot["large_bee_like"], errors="coerce")
        + pd.to_numeric(pivot["generalized_accessible"], errors="coerce")
    ) / 2.0
    return pivot


def _fit_distance(
    work: pd.DataFrame,
    *,
    geography: str,
    baseline: list[str],
    cluster: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    needed = ["attraction_shift", geography, *baseline, cluster]
    complete = work.dropna(subset=needed).copy()
    if complete.empty:
        return pd.DataFrame(), {"status": "no_complete_rows"}
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    columns = [
        np.ones(len(complete), dtype=float),
        _standardize(complete[geography]),
        *[_standardize(complete[x]) for x in baseline],
    ]
    coef, _, fit = _fit_weighted_clustered_design(
        pd.to_numeric(complete["attraction_shift"], errors="coerce").to_numpy(float),
        np.ones(len(complete), dtype=float),
        np.column_stack(columns),
        names,
        complete[cluster].astype(str).to_numpy(),
    )
    return coef, {
        **fit,
        "n_unique_islands": int(complete["island_id"].nunique()),
    }


def analyse_restricted_attraction(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    *,
    context_column: str,
    contexts: list[str],
    context_layer: str,
    restriction: str,
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    strata = [str(x) for x in pattern_config["strata"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}

    table = _attraction_table(adjusted_scores)
    needed_cov = ["island_id", geography, context_column, cluster, *baseline]
    missing = set(needed_cov) - set(covariates.columns)
    if missing:
        raise ValueError(f"covariates missing columns: {sorted(missing)}")
    data = table.merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["attraction_shift", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context_column] = data[context_column].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)

    rows: list[dict[str, Any]] = []
    for source_mode in data["source_mode"].dropna().astype(str).drop_duplicates():
        mode = data.loc[data["source_mode"].eq(source_mode)].copy()
        for stratum in strata:
            for context in contexts:
                base = mode.loc[
                    mode["stratum"].eq(stratum)
                    & mode[context_column].eq(context)
                ].copy()
                n_complete = int(
                    base.dropna(
                        subset=["attraction_shift", geography, *baseline, cluster]
                    )["island_id"].nunique()
                )
                for tier, threshold in tiers.items():
                    row: dict[str, Any] = {
                        "restriction": restriction,
                        "context_layer": context_layer,
                        "source_mode": source_mode,
                        "stratum": stratum,
                        "context": context,
                        "support_tier": tier,
                        "threshold": threshold,
                        "n_unique_islands": n_complete,
                    }
                    if n_complete < threshold:
                        rows.append({**row, "status": "not_testable"})
                        continue
                    coef, fit = _fit_distance(
                        base,
                        geography=geography,
                        baseline=baseline,
                        cluster=cluster,
                    )
                    if coef.empty:
                        rows.append(
                            {
                                **row,
                                "status": str(fit.get("status", "fit_failed")),
                            }
                        )
                        continue
                    distance = coef.set_index("predictor").loc[f"z_{geography}"]
                    rows.append(
                        {
                            **row,
                            "status": "fit",
                            "n_clusters": int(fit["n_clusters"]),
                            "distance_estimate": float(distance["estimate"]),
                            "distance_se": float(distance["cluster_robust_se"]),
                            "distance_p": float(distance["p_value"]),
                        }
                    )

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["distance_q"] = np.nan
    fit_mask = result["status"].eq("fit")
    if fit_mask.any():
        family = ["restriction", "source_mode", "stratum", "support_tier"]
        result.loc[fit_mask, "distance_q"] = (
            result.loc[fit_mask]
            .groupby(family, group_keys=False)["distance_p"]
            .apply(_bh)
        )
    return result


def run_restrictions(
    trait_ledger: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    source_assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
    source_config: dict[str, Any],
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    restrictions = reproductive_restriction_sets(trait_ledger)
    species_scores = score_species_syndromes(trait_ledger, syndrome_config)
    strata = [str(x) for x in source_config["strata"]]
    min_species = int(
        source_config["source_region_scores"][
            "minimum_trait_scored_species_per_region_syndrome"
        ]
    )
    min_sources = int(
        source_config["response"][
            "source_expectation_requires_minimum_source_regions"
        ]
    )

    result_parts: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    adjusted_parts: list[pd.DataFrame] = []
    for restriction, species_set in restrictions.items():
        restricted_species_scores = species_scores.loc[
            species_scores["accepted_species"].isin(species_set)
        ].copy()
        island_scores = build_island_syndrome_scores(
            status_flora,
            restricted_species_scores,
            strata,
        )
        matched, _ = match_gift_species_to_frozen_scores(
            gift_flora,
            restricted_species_scores,
        )
        mainland_scores = build_mainland_syndrome_scores(
            matched,
            min_species=min_species,
        )
        expectations = build_source_expectations(
            source_assignments,
            mainland_scores,
            min_source_regions=min_sources,
        )
        adjusted = build_adjusted_island_scores(
            island_scores,
            expectations,
            strata,
        )
        adjusted.insert(0, "restriction", restriction)
        adjusted_parts.append(adjusted)

        scored_attraction_species = restricted_species_scores.loc[
            restricted_species_scores["syndrome"].isin(
                ["large_bee_like", "generalized_accessible"]
            ),
            "accepted_species",
        ].nunique()
        eligible_mainland = mainland_scores.loc[
            mainland_scores["syndrome"].isin(
                ["large_bee_like", "generalized_accessible"]
            )
            & mainland_scores["source_score_eligible"].astype(bool)
        ]
        audit_rows.append(
            {
                "restriction": restriction,
                "n_reproductively_classified_species": len(species_set),
                "n_species_with_attraction_syndrome_score": int(
                    scored_attraction_species
                ),
                "n_eligible_mainland_region_syndrome_rows": int(
                    len(eligible_mainland)
                ),
                "n_islands_with_adjusted_scores": int(
                    adjusted["island_id"].nunique()
                ),
            }
        )

        regime = analyse_restricted_attraction(
            adjusted,
            covariates,
            pattern_config,
            context_column=str(pattern_config["context_column"]),
            contexts=[str(x) for x in pattern_config["contexts"]],
            context_layer="analysis_regime",
            restriction=restriction,
        )
        result_parts.append(regime)

        realm_cov = covariates.merge(
            realm_assignment[
                ["island_id", "biogeographic_realm"]
            ].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="one_to_one",
        )
        realm_contexts = sorted(
            realm_cov["biogeographic_realm"]
            .dropna()
            .astype(str)
            .loc[lambda x: x.ne("")]
            .unique()
        )
        realm = analyse_restricted_attraction(
            adjusted,
            realm_cov,
            pattern_config,
            context_column="biogeographic_realm",
            contexts=realm_contexts,
            context_layer="biogeographic_realm",
            restriction=restriction,
        )
        result_parts.append(realm)

    results = (
        pd.concat(result_parts, ignore_index=True)
        if result_parts
        else pd.DataFrame()
    )
    adjusted_all = (
        pd.concat(adjusted_parts, ignore_index=True)
        if adjusted_parts
        else pd.DataFrame()
    )
    audit = pd.DataFrame(audit_rows)
    manifest = {
        "contract": "chapter1_pr138_outcrossing_restriction_v1",
        "restriction_definitions": {
            "si_only": "self_incompatibility exactly SI",
            "predominantly_outcrossing": (
                "mating_system exactly predominantly_outcrossing"
            ),
            "si_and_predominantly_outcrossing": (
                "both exact SI and exact predominantly_outcrossing"
            ),
        },
        "ambiguous_reproductive_states_coerced": False,
        "trait_evidence_layer": "all_analysis_eligible",
        "attraction_shift": "(-large_bee_like + generalized_accessible) / 2",
        "source_pool_adjustment_recomputed_within_reproductive_subset": True,
        "causal_pollinator_claimed": False,
        "interpretation_ceiling": (
            "Persistence in SI taxa shows measured selfing is not necessary for the "
            "attraction/accessibility response; it does not identify relaxed attraction "
            "selection or any pollinator guild causally."
        ),
    }
    return {
        "results": results,
        "audit": audit,
        "adjusted_scores": adjusted_all,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-ledger-csv", type=Path, required=True)
    parser.add_argument("--status-flora-csv", type=Path, required=True)
    parser.add_argument("--gift-flora-csv", type=Path, required=True)
    parser.add_argument("--source-assignments-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--syndrome-config-path", type=Path, required=True)
    parser.add_argument("--source-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    outputs = run_restrictions(
        pd.read_csv(args.trait_ledger_csv),
        pd.read_csv(args.status_flora_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.syndrome_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.source_config_path.read_text(encoding="utf-8")),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in {
        "results": "outcrossing_restriction_results.csv",
        "audit": "outcrossing_restriction_audit.csv",
        "adjusted_scores": "outcrossing_restriction_adjusted_scores.csv.gz",
    }.items():
        frame = outputs[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(args.output_dir / filename, index=False)
    manifest = outputs["manifest"]
    assert isinstance(manifest, dict)
    (args.output_dir / "outcrossing_restriction_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
