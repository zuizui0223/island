"""Source-pool-standardized sensitivity for PR138.

The primary PR138 estimand remains the observed island assemblage. This module is a
secondary attack on the alternative explanation that regional distance responses merely
reflect different mainland source-flora composition.

Source assignment is outcome-blind. The same frozen PR138 species syndrome scores are
used on islands and on external GIFT mainland native floras.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import unicodedata
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_pr138_syndrome_analysis import _between_contexts, _bh, _prepare, _within_context


def _normalise_name(value: object) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).replace("×", " ")
    return " ".join(text.split()).casefold()


def _binomial_key(value: object) -> str:
    tokens = re.findall(r"[A-Za-z][A-Za-z.-]*", unicodedata.normalize("NFKC", str(value or "")).replace("×", " "))
    if len(tokens) < 2:
        return ""
    return f"{tokens[0].casefold()} {tokens[1].casefold()}"


def match_gift_species_to_frozen_scores(
    gift_flora: pd.DataFrame, species_scores: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_flora = {"entity_ID", "work_species"}
    required_scores = {"accepted_species", "syndrome", "syndrome_concordance"}
    if missing := required_flora - set(gift_flora.columns):
        raise ValueError(f"GIFT flora missing columns: {sorted(missing)}")
    if missing := required_scores - set(species_scores.columns):
        raise ValueError(f"species scores missing columns: {sorted(missing)}")

    unique_scores = species_scores[["accepted_species"]].drop_duplicates().copy()
    unique_scores["norm"] = unique_scores["accepted_species"].map(_normalise_name)
    unique_scores["binomial"] = unique_scores["accepted_species"].map(_binomial_key)

    exact_groups = unique_scores.groupby("norm")["accepted_species"].agg(lambda x: sorted(set(x)))
    exact = {key: values[0] for key, values in exact_groups.items() if key and len(values) == 1}
    binomial_groups = unique_scores.groupby("binomial")["accepted_species"].agg(lambda x: sorted(set(x)))
    binomial = {key: values[0] for key, values in binomial_groups.items() if key and len(values) == 1}

    flora = gift_flora[["entity_ID", "work_species"]].drop_duplicates().copy()
    match_rows: list[dict[str, Any]] = []
    for row in flora.itertuples(index=False):
        norm = _normalise_name(row.work_species)
        key = _binomial_key(row.work_species)
        accepted = exact.get(norm)
        method = "exact_normalized" if accepted else ""
        if accepted is None and key:
            accepted = binomial.get(key)
            if accepted is not None:
                method = "unique_binomial_fallback"
        match_rows.append(
            {
                "entity_ID": int(row.entity_ID),
                "work_species": str(row.work_species),
                "accepted_species": accepted or "",
                "match_method": method or "unmatched",
            }
        )
    matches = pd.DataFrame(match_rows)
    matched = matches.loc[matches["accepted_species"].ne("")].copy()
    expanded = matched.merge(species_scores, on="accepted_species", how="inner", validate="many_to_many")
    audit = (
        matches.groupby("match_method", as_index=False)
        .agg(n_gift_entity_species_rows=("work_species", "size"), n_entities=("entity_ID", "nunique"))
        .sort_values("match_method")
    )
    return expanded, audit


def build_mainland_syndrome_scores(
    matched_scores: pd.DataFrame, *, min_species: int = 10
) -> pd.DataFrame:
    work = matched_scores.drop_duplicates(["entity_ID", "accepted_species", "syndrome"]).copy()
    out = (
        work.groupby(["entity_ID", "syndrome"], as_index=False)
        .agg(
            mainland_syndrome_score=("syndrome_concordance", "mean"),
            n_trait_scored_species=("accepted_species", "nunique"),
        )
    )
    out["source_score_eligible"] = out["n_trait_scored_species"].ge(int(min_species))
    return out


def _haversine_vector(lat: float, lon: float, source_lat: np.ndarray, source_lon: np.ndarray) -> np.ndarray:
    radius = 6371.0088
    lat1 = math.radians(float(lat))
    lon1 = math.radians(float(lon))
    lat2 = np.radians(source_lat.astype(float))
    lon2 = np.radians(source_lon.astype(float))
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + math.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 2.0 * radius * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def build_source_assignments(
    covariates: pd.DataFrame,
    source_regions: pd.DataFrame,
    source_config: dict[str, Any],
) -> pd.DataFrame:
    required_cov = {"island_id", "island_latitude", "island_longitude", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"}
    required_src = {"entity_ID", "source_latitude", "source_longitude", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"}
    if missing := required_cov - set(covariates.columns):
        raise ValueError(f"island covariates missing source-assignment columns: {sorted(missing)}")
    if missing := required_src - set(source_regions.columns):
        raise ValueError(f"source regions missing assignment columns: {sorted(missing)}")

    modes = source_config["source_assignment"]["primary_modes"]
    src = source_regions.copy()
    for column in ["source_latitude", "source_longitude", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]:
        src[column] = pd.to_numeric(src[column], errors="coerce")
    src_geo = src.dropna(subset=["source_latitude", "source_longitude"]).drop_duplicates("entity_ID").reset_index(drop=True)
    src_climate = src_geo.dropna(subset=["climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]).copy()
    if len(src_geo) < max(int(spec.get("k", 0)) for spec in modes.values()):
        raise ValueError("Too few GIFT source regions for the prespecified geographic source modes")

    rows: list[dict[str, Any]] = []
    geo_lats = src_geo["source_latitude"].to_numpy(float)
    geo_lons = src_geo["source_longitude"].to_numpy(float)
    for island in covariates.itertuples(index=False):
        lat = pd.to_numeric(pd.Series([getattr(island, "island_latitude")]), errors="coerce").iloc[0]
        lon = pd.to_numeric(pd.Series([getattr(island, "island_longitude")]), errors="coerce").iloc[0]
        if not (math.isfinite(lat) and math.isfinite(lon)):
            continue
        distances = _haversine_vector(float(lat), float(lon), geo_lats, geo_lons)
        geo_order = np.argsort(distances)
        for mode, spec in modes.items():
            rule = str(spec["rule"])
            if rule.endswith("nearest_source_region_representative_points"):
                k = int(spec["k"])
                selected = geo_order[:k]
                for rank, idx in enumerate(selected, start=1):
                    rows.append(
                        {
                            "island_id": str(island.island_id),
                            "source_mode": str(mode),
                            "source_rank": rank,
                            "entity_ID": int(src_geo.iloc[idx]["entity_ID"]),
                            "source_distance_km": float(distances[idx]),
                            "climate_distance": np.nan,
                        }
                    )
            elif rule == "ten_climate_nearest_among_fifty_geographically_nearest_source_regions":
                island_pc = np.array([getattr(island, f"climate_pc{i}") for i in range(1, 5)], dtype=float)
                if not np.isfinite(island_pc).all() or src_climate.empty:
                    continue
                candidate_k = min(int(spec["geographic_candidate_k"]), len(src_climate))
                climate_dist_geo = _haversine_vector(
                    float(lat),
                    float(lon),
                    src_climate["source_latitude"].to_numpy(float),
                    src_climate["source_longitude"].to_numpy(float),
                )
                candidate_idx = np.argsort(climate_dist_geo)[:candidate_k]
                candidate = src_climate.iloc[candidate_idx].copy()
                matrix = candidate[[f"climate_pc{i}" for i in range(1, 5)]].to_numpy(float)
                climate_distance = np.sqrt(((matrix - island_pc[None, :]) ** 2).sum(axis=1))
                order = np.argsort(climate_distance)[: int(spec["k"])]
                for rank, local_idx in enumerate(order, start=1):
                    source_row = candidate.iloc[local_idx]
                    original_position = candidate_idx[local_idx]
                    rows.append(
                        {
                            "island_id": str(island.island_id),
                            "source_mode": str(mode),
                            "source_rank": rank,
                            "entity_ID": int(source_row["entity_ID"]),
                            "source_distance_km": float(climate_dist_geo[original_position]),
                            "climate_distance": float(climate_distance[local_idx]),
                        }
                    )
            else:
                raise ValueError(f"Unknown frozen source assignment rule: {rule}")
    return pd.DataFrame(rows)


def build_source_expectations(
    assignments: pd.DataFrame,
    mainland_scores: pd.DataFrame,
    *,
    min_source_regions: int = 3,
) -> pd.DataFrame:
    scores = mainland_scores.loc[mainland_scores["source_score_eligible"].astype(bool)].copy()
    merged = assignments.merge(scores, on="entity_ID", how="left", validate="many_to_many")
    merged = merged.dropna(subset=["mainland_syndrome_score"])
    out = (
        merged.groupby(["island_id", "source_mode", "syndrome"], as_index=False)
        .agg(
            source_expectation=("mainland_syndrome_score", "mean"),
            n_source_regions=("entity_ID", "nunique"),
            median_source_distance_km=("source_distance_km", "median"),
            median_source_trait_species=("n_trait_scored_species", "median"),
        )
    )
    out["source_expectation_eligible"] = out["n_source_regions"].ge(int(min_source_regions))
    return out


def build_adjusted_island_scores(
    island_scores: pd.DataFrame,
    expectations: pd.DataFrame,
    strata: list[str],
) -> pd.DataFrame:
    islands = island_scores.loc[island_scores["stratum"].isin(strata)].copy()
    merged = islands.merge(expectations, on=["island_id", "syndrome"], how="inner", validate="many_to_many")
    merged = merged.loc[merged["source_expectation_eligible"].astype(bool)].copy()
    merged["observed_island_syndrome_score"] = pd.to_numeric(merged["syndrome_score"], errors="coerce")
    merged["syndrome_score"] = merged["observed_island_syndrome_score"] - pd.to_numeric(
        merged["source_expectation"], errors="coerce"
    )
    merged["source_adjusted_response"] = merged["syndrome_score"]
    return merged


def _analyse_context_layer(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    *,
    context_layer: str,
    source_modes: list[str],
    strata: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contexts = [str(x) for x in pattern_config["contexts"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}
    slope_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for source_mode in source_modes:
        mode_scores = adjusted_scores.loc[adjusted_scores["source_mode"].eq(source_mode)].copy()
        data = _prepare(mode_scores, covariates, pattern_config)
        for stratum in strata:
            for tier, threshold in tiers.items():
                for context in contexts:
                    slopes, result = _within_context(
                        data,
                        stratum=stratum,
                        context_value=context,
                        support_tier=tier,
                        threshold=threshold,
                        pattern_config=pattern_config,
                        syndrome_config={},
                    )
                    if not slopes.empty:
                        slopes.insert(0, "source_mode", source_mode)
                        slopes.insert(0, "context_layer", context_layer)
                        slope_parts.append(slopes)
                    result.update({"source_mode": source_mode, "context_layer": context_layer})
                    within_rows.append(result)
                for context_a, context_b in itertools.combinations(contexts, 2):
                    result = _between_contexts(
                        data,
                        stratum=stratum,
                        context_a=context_a,
                        context_b=context_b,
                        support_tier=tier,
                        threshold=threshold,
                        pattern_config=pattern_config,
                    )
                    result.update({"source_mode": source_mode, "context_layer": context_layer})
                    between_rows.append(result)
    slopes = pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame()
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    family = ["source_mode", "stratum", "support_tier"]
    if "p_value" in within.columns:
        within["q_within_stratum_tier"] = within.groupby(family, group_keys=False)["p_value"].transform(_bh)
        within["syndrome_vector_supported"] = within["q_within_stratum_tier"].le(0.05).fillna(False)
    if "northern_classic_projection_one_sided_p" in within.columns:
        within["q_classic_projection"] = within.groupby(family, group_keys=False)[
            "northern_classic_projection_one_sided_p"
        ].transform(_bh)
        within["classic_projection_supported"] = (
            within["q_classic_projection"].le(0.05)
            & pd.to_numeric(within["northern_classic_projection"], errors="coerce").gt(0)
        ).fillna(False)
    if "p_value" in between.columns:
        between["q_within_stratum_tier"] = between.groupby(family, group_keys=False)["p_value"].transform(_bh)
        between["regional_syndrome_vector_difference_supported"] = between[
            "q_within_stratum_tier"
        ].le(0.05).fillna(False)
    return slopes, within, between


def run(
    species_scores: pd.DataFrame,
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    source_regions: pd.DataFrame,
    gift_flora: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
    realm_assignment: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    min_species = int(source_config["source_region_scores"]["minimum_trait_scored_species_per_region_syndrome"])
    min_sources = int(source_config["response"]["source_expectation_requires_minimum_source_regions"])
    strata = [str(x) for x in source_config["strata"]]
    source_modes = list(source_config["source_assignment"]["primary_modes"])

    matched, match_audit = match_gift_species_to_frozen_scores(gift_flora, species_scores)
    mainland_scores = build_mainland_syndrome_scores(matched, min_species=min_species)
    assignments = build_source_assignments(covariates, source_regions, source_config)
    expectations = build_source_expectations(assignments, mainland_scores, min_source_regions=min_sources)
    adjusted = build_adjusted_island_scores(island_scores, expectations, strata)

    broad_config = deepcopy(pattern_config)
    broad_config["strata"] = strata
    broad_slopes, broad_within, broad_between = _analyse_context_layer(
        adjusted,
        covariates,
        broad_config,
        context_layer="analysis_regime",
        source_modes=source_modes,
        strata=strata,
    )

    realm_slopes = pd.DataFrame()
    realm_within = pd.DataFrame()
    realm_between = pd.DataFrame()
    if realm_assignment is not None and not realm_assignment.empty:
        required_realm = {"island_id", "biogeographic_realm"}
        if missing := required_realm - set(realm_assignment.columns):
            raise ValueError(f"realm assignment missing columns: {sorted(missing)}")
        realm_cov = covariates.merge(
            realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="one_to_one",
        )
        realm_config = deepcopy(pattern_config)
        realm_config["context_column"] = "biogeographic_realm"
        realm_config["contexts"] = sorted(
            realm_cov["biogeographic_realm"].dropna().astype(str).loc[lambda s: s.ne("")].unique()
        )
        realm_config["strata"] = strata
        realm_slopes, realm_within, realm_between = _analyse_context_layer(
            adjusted,
            realm_cov,
            realm_config,
            context_layer="biogeographic_realm",
            source_modes=source_modes,
            strata=strata,
        )

    match_counts = match_audit.set_index("match_method")["n_gift_entity_species_rows"].to_dict() if not match_audit.empty else {}
    manifest = {
        "contract": "chapter1_pr138_source_pool_sensitivity_v1",
        "source_modes": source_modes,
        "strata": strata,
        "n_gift_entity_species_rows": int(len(gift_flora)),
        "n_exact_species_matches": int(match_counts.get("exact_normalized", 0)),
        "n_unique_binomial_fallback_matches": int(match_counts.get("unique_binomial_fallback", 0)),
        "n_unmatched_gift_species_rows": int(match_counts.get("unmatched", 0)),
        "n_mainland_entity_syndrome_scores": int(len(mainland_scores)),
        "n_source_assignments": int(len(assignments)),
        "n_source_expectations": int(len(expectations)),
        "n_source_adjusted_island_scores": int(len(adjusted)),
        "source_assignment_uses_island_traits": False,
        "same_frozen_species_scores_used_for_island_and_mainland": True,
        "primary_observed_assemblage_estimand_replaced": False,
        "claim_boundary": "Source-pool sensitivity only; persistence is not proof of historical ancestry, in-situ evolution, or a pollinator mechanism.",
    }
    return {
        "match_audit": match_audit,
        "mainland_scores": mainland_scores,
        "assignments": assignments,
        "expectations": expectations,
        "adjusted_scores": adjusted,
        "broad_slopes": broad_slopes,
        "broad_within": broad_within,
        "broad_between": broad_between,
        "realm_slopes": realm_slopes,
        "realm_within": realm_within,
        "realm_between": realm_between,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-scores-csv", type=Path, required=True)
    parser.add_argument("--island-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--source-regions-csv", type=Path, required=True)
    parser.add_argument("--gift-flora-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--source-config-path", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pattern_config = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    source_config = yaml.safe_load(args.source_config_path.read_text(encoding="utf-8"))
    realm = pd.read_csv(args.realm_assignment_csv) if args.realm_assignment_csv else None
    outputs = run(
        pd.read_csv(args.species_scores_csv),
        pd.read_csv(args.island_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.source_regions_csv),
        pd.read_csv(args.gift_flora_csv),
        pattern_config,
        source_config,
        realm,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in {
        "match_audit": "gift_trait_species_match_audit.csv",
        "mainland_scores": "mainland_syndrome_scores.csv.gz",
        "assignments": "island_source_assignments.csv.gz",
        "expectations": "island_source_expectations.csv.gz",
        "adjusted_scores": "source_adjusted_island_syndrome_scores.csv.gz",
        "broad_slopes": "source_adjusted_regime_slopes.csv",
        "broad_within": "source_adjusted_regime_within.csv",
        "broad_between": "source_adjusted_regime_between.csv",
        "realm_slopes": "source_adjusted_realm_slopes.csv",
        "realm_within": "source_adjusted_realm_within.csv",
        "realm_between": "source_adjusted_realm_between.csv",
    }.items():
        frame = outputs[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(args.output_dir / filename, index=False)
    manifest = outputs["manifest"]
    assert isinstance(manifest, dict)
    (args.output_dir / "source_pool_sensitivity_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
