"""Source-pool lineage representation bridge for PR138 and izu-core.

This secondary analysis asks why genus composition can reproduce most of the PR138
floral-accessibility geography. It tests whether genera with different mainland-derived
functional positions are represented non-randomly on islands after conditioning on
source availability, then separates binary genus entry from species loading within
represented genera.

The source functional coordinate is fixed from GIFT mainland species only::

    (-large_bee_like + generalized_accessible) / 2

For each island and source definition, candidate genera are grouped by the number of
assigned source regions in which they occur and, in the primary matching scheme, by a
coarse source-species-richness bin. The observed number of represented genera and island
species in every availability class is held fixed. Expected functional position is the
class mean, so the enrichment statistic cannot be generated merely by having fewer
observed genera or by losing low-prevalence source genera.

This is an assembly diagnostic. It does not identify historical source ancestry,
pollinator loss, colonisation, extinction, diversification, or in-situ evolution.
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
SOURCE_MATCHING = ("prevalence_richness", "prevalence_only")
OUTCOMES = ("entry_enrichment", "species_enrichment", "loading_increment")
FOCAL_CONTEXTS = {
    "analysis_regime": ("northern_midlatitude", "tropical"),
    "biogeographic_realm": ("Palearctic", "Neotropical"),
}


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
    text = unicodedata.normalize("NFKC", str(value or "")).replace("×", " ")
    tokens = re.findall(r"[A-Za-z][A-Za-z.-]*", text)
    return tokens[0] if tokens else ""


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


def exact_si_species(trait_ledger: pd.DataFrame) -> set[str]:
    required = {"accepted_species", "trait_name", "normalized_value"}
    if missing := required - set(trait_ledger.columns):
        raise ValueError(f"trait ledger missing columns: {sorted(missing)}")
    frame = trait_ledger.copy()
    mask = (
        frame["trait_name"].fillna("").astype(str).eq("self_incompatibility")
        & frame["normalized_value"].fillna("").astype(str).eq("SI")
    )
    return set(frame.loc[mask, "accepted_species"].dropna().astype(str))


def functional_species_positions(
    species_scores: pd.DataFrame,
    species_filter: set[str] | None = None,
) -> pd.DataFrame:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    if missing := required - set(species_scores.columns):
        raise ValueError(f"species scores missing columns: {sorted(missing)}")
    work = species_scores.loc[
        species_scores["syndrome"].astype(str).isin(AXES),
        ["accepted_species", "syndrome", "syndrome_concordance"],
    ].copy()
    if species_filter is not None:
        work = work.loc[
            work["accepted_species"].astype(str).isin(species_filter)
        ].copy()
    work["score"] = pd.to_numeric(work["syndrome_concordance"], errors="coerce")
    work = work.dropna(subset=["score"])
    wide = work.pivot_table(
        index="accepted_species", columns="syndrome", values="score", aggfunc="first"
    )
    wide = wide.dropna(subset=list(AXES)).copy()
    wide["functional_position"] = (
        -wide["large_bee_like"] + wide["generalized_accessible"]
    ) / 2.0
    out = wide[["functional_position"]].reset_index()
    out["accepted_species"] = out["accepted_species"].astype(str)
    out["genus"] = out["accepted_species"].map(_genus)
    return out.loc[out["genus"].ne("")].reset_index(drop=True)


def match_gift_species(
    gift_flora: pd.DataFrame, species_universe: pd.DataFrame
) -> pd.DataFrame:
    required_gift = {"entity_ID", "work_species"}
    required_species = {"accepted_species"}
    if missing := required_gift - set(gift_flora.columns):
        raise ValueError(f"GIFT flora missing columns: {sorted(missing)}")
    if missing := required_species - set(species_universe.columns):
        raise ValueError(f"species universe missing columns: {sorted(missing)}")

    species = species_universe[["accepted_species"]].drop_duplicates().copy()
    species["accepted_species"] = species["accepted_species"].astype(str)
    species["norm"] = species["accepted_species"].map(_normalise_name)
    species["binomial"] = species["accepted_species"].map(_binomial_key)
    exact_groups = species.groupby("norm")["accepted_species"].agg(
        lambda x: sorted(set(x))
    )
    exact = {
        key: vals[0] for key, vals in exact_groups.items() if key and len(vals) == 1
    }
    binomial_groups = species.groupby("binomial")["accepted_species"].agg(
        lambda x: sorted(set(x))
    )
    binomial = {
        key: vals[0]
        for key, vals in binomial_groups.items()
        if key and len(vals) == 1
    }

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
    return (
        gift_flora[["entity_ID", "work_species"]]
        .merge(matched_names, on="work_species", how="inner", validate="many_to_one")
        .drop_duplicates(["entity_ID", "accepted_species"])
    )


def source_genus_positions(
    gift_flora: pd.DataFrame,
    species_scores: pd.DataFrame,
    species_filter: set[str] | None = None,
    *,
    minimum_source_scored_species: int = 1,
    matched_gift: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    functional = functional_species_positions(species_scores, species_filter)
    matched = (
        matched_gift.copy()
        if matched_gift is not None
        else match_gift_species(gift_flora, species_scores[["accepted_species"]])
    )
    source_species = (
        matched.merge(
            functional, on="accepted_species", how="inner", validate="many_to_one"
        )
        .drop_duplicates("accepted_species")
        .copy()
    )
    positions = source_species.groupby("genus", as_index=False).agg(
        source_functional_position=("functional_position", "mean"),
        n_source_scored_species=("accepted_species", "nunique"),
    )
    positions = positions.loc[
        positions["n_source_scored_species"].ge(int(minimum_source_scored_species))
    ].reset_index(drop=True)
    return positions, matched, functional


def gift_source_functional_counts(
    matched_gift: pd.DataFrame,
    functional_species: pd.DataFrame,
) -> dict[str, int]:
    """Separate GIFT entity-species memberships from unique matched species."""
    required_matched = {"entity_ID", "accepted_species"}
    required_functional = {"accepted_species"}
    if missing := required_matched - set(matched_gift.columns):
        raise ValueError(f"matched GIFT flora missing columns: {sorted(missing)}")
    if missing := required_functional - set(functional_species.columns):
        raise ValueError(f"functional species missing columns: {sorted(missing)}")

    eligible = set(functional_species["accepted_species"].dropna().astype(str))
    memberships = matched_gift.loc[
        matched_gift["accepted_species"].astype(str).isin(eligible),
        ["entity_ID", "accepted_species"],
    ].drop_duplicates(["entity_ID", "accepted_species"])
    return {
        "n_gift_source_functional_memberships": len(memberships),
        "n_unique_gift_source_functional_species": int(
            memberships["accepted_species"].astype(str).nunique()
        ),
    }


def broad_source_availability(
    gift_flora: pd.DataFrame, genera: set[str]
) -> pd.DataFrame:
    work = gift_flora[["entity_ID", "work_species"]].drop_duplicates().copy()
    work["genus"] = work["work_species"].map(_genus)
    work = work.loc[work["genus"].isin(genera) & work["genus"].ne("")].copy()
    return work.groupby(["entity_ID", "genus"], as_index=False).agg(
        source_species_richness=("work_species", "nunique")
    )


def restricted_source_availability(
    matched_gift: pd.DataFrame,
    functional_species: pd.DataFrame,
    genera: set[str],
) -> pd.DataFrame:
    work = matched_gift.merge(
        functional_species[["accepted_species", "genus"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    work = work.loc[work["genus"].isin(genera)].drop_duplicates(
        ["entity_ID", "accepted_species"]
    )
    return work.groupby(["entity_ID", "genus"], as_index=False).agg(
        source_species_richness=("accepted_species", "nunique")
    )


def _stratum_mask(frame: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all_native":
        return frame["origin_status"].astype(str).eq("native")
    if stratum == "native_nonendemic":
        return frame["floristic_status"].astype(str).eq("native_nonendemic")
    raise ValueError(f"unsupported stratum: {stratum}")


def _representation_count_matrix(
    status_flora: pd.DataFrame,
    island_index: dict[str, int],
    genus_index: dict[str, int],
    *,
    stratum: str,
    species_filter: set[str] | None,
) -> sparse.csr_matrix:
    work = status_flora.loc[
        _stratum_mask(status_flora, stratum), ["island_id", "accepted_species"]
    ].drop_duplicates().copy()
    if species_filter is not None:
        work = work.loc[
            work["accepted_species"].astype(str).isin(species_filter)
        ].copy()
    work["island_id"] = work["island_id"].astype(str)
    work["genus"] = work["accepted_species"].map(_genus)
    work = work.loc[
        work["island_id"].isin(island_index) & work["genus"].isin(genus_index)
    ].copy()
    grouped = work.groupby(["island_id", "genus"], as_index=False).agg(
        n_island_species=("accepted_species", "nunique")
    )
    rows = grouped["island_id"].map(island_index).to_numpy(int)
    cols = grouped["genus"].map(genus_index).to_numpy(int)
    return sparse.csr_matrix(
        (grouped["n_island_species"].to_numpy(float), (rows, cols)),
        shape=(len(island_index), len(genus_index)),
    )


def _availability_matrices(
    availability: pd.DataFrame,
    entity_index: dict[int, int],
    genus_index: dict[str, int],
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    work = availability.loc[
        availability["entity_ID"].astype(int).isin(entity_index)
        & availability["genus"].astype(str).isin(genus_index)
    ].copy()
    rows = work["entity_ID"].astype(int).map(entity_index).to_numpy(int)
    cols = work["genus"].astype(str).map(genus_index).to_numpy(int)
    presence = sparse.csr_matrix(
        (np.ones(len(work), dtype=float), (rows, cols)),
        shape=(len(entity_index), len(genus_index)),
    )
    richness = sparse.csr_matrix(
        (
            pd.to_numeric(work["source_species_richness"], errors="coerce")
            .fillna(0)
            .to_numpy(float),
            (rows, cols),
        ),
        shape=presence.shape,
    )
    return presence, richness


def _richness_bins(values: np.ndarray) -> np.ndarray:
    bins = np.zeros_like(values, dtype=np.int8)
    bins[(values > 1) & (values <= 2)] = 1
    bins[(values > 2) & (values <= 4)] = 2
    bins[(values > 4) & (values <= 8)] = 3
    bins[values > 8] = 4
    return bins


def compute_island_enrichment(
    prevalence: np.ndarray,
    source_richness: np.ndarray,
    island_species_counts: np.ndarray,
    positions: np.ndarray,
    *,
    matching: str,
    minimum_represented_genera: int,
) -> dict[str, float | int] | None:
    candidate = prevalence > 0
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

    if matching == "prevalence_only":
        class_pairs = [(int(x), -1) for x in np.unique(prevalence[candidate])]
        richness_bins = np.zeros_like(prevalence, dtype=np.int8)
    elif matching == "prevalence_richness":
        richness_bins = _richness_bins(source_richness)
        class_pairs = [
            (int(x[0]), int(x[1]))
            for x in np.unique(
                np.column_stack([prevalence[candidate], richness_bins[candidate]]),
                axis=0,
            )
        ]
    else:
        raise ValueError(f"unknown source matching scheme: {matching}")

    for source_prevalence, richness_bin in class_pairs:
        class_mask = candidate & (prevalence == source_prevalence)
        if matching == "prevalence_richness":
            class_mask &= richness_bins == richness_bin
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
        "n_candidate_genera": int(candidate.sum()),
        "observed_entry_mean": observed_entry,
        "expected_entry_mean": expected_entry,
        "entry_enrichment": entry,
        "observed_species_mean": observed_species,
        "expected_species_mean": expected_species,
        "species_enrichment": species,
        "loading_increment": species - entry,
    }


def _source_assignment_matrix(
    assignments: pd.DataFrame,
    island_index: dict[str, int],
    entity_index: dict[int, int],
    *,
    source_mode: str,
) -> sparse.csr_matrix:
    work = assignments.loc[
        assignments["source_mode"].astype(str).eq(source_mode),
        ["island_id", "entity_ID"],
    ].drop_duplicates().copy()
    work["island_id"] = work["island_id"].astype(str)
    work["entity_ID"] = pd.to_numeric(work["entity_ID"], errors="coerce")
    work = work.dropna(subset=["entity_ID"])
    work["entity_ID"] = work["entity_ID"].astype(int)
    work = work.loc[
        work["island_id"].isin(island_index) & work["entity_ID"].isin(entity_index)
    ]
    rows = work["island_id"].map(island_index).to_numpy(int)
    cols = work["entity_ID"].map(entity_index).to_numpy(int)
    return sparse.csr_matrix(
        (np.ones(len(work), dtype=float), (rows, cols)),
        shape=(len(island_index), len(entity_index)),
    )


def build_evidence_scores(
    *,
    evidence_scope: str,
    species_scores: pd.DataFrame,
    trait_ledger: pd.DataFrame | None,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    source_modes: list[str],
    minimum_source_scored_species: int,
    minimum_represented_genera: list[int],
    restricted_to_exact_si: bool,
    raw_gift_availability: bool,
    matched_gift: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    species_filter = None
    if restricted_to_exact_si:
        if trait_ledger is None:
            raise ValueError("exact-SI evidence requires a trait ledger")
        species_filter = exact_si_species(trait_ledger)

    positions, matched, functional = source_genus_positions(
        gift_flora,
        species_scores,
        species_filter,
        minimum_source_scored_species=minimum_source_scored_species,
        matched_gift=matched_gift,
    )
    genera = sorted(positions["genus"].astype(str).unique())
    genus_index = {genus: idx for idx, genus in enumerate(genera)}
    position_array = (
        positions.set_index("genus")
        .loc[genera, "source_functional_position"]
        .to_numpy(float)
    )
    if raw_gift_availability:
        availability = broad_source_availability(gift_flora, set(genera))
    else:
        availability = restricted_source_availability(matched, functional, set(genera))

    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: idx for idx, island in enumerate(islands)}
    entities = sorted(
        set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
        | set(pd.to_numeric(availability["entity_ID"], errors="coerce").dropna().astype(int))
    )
    entity_index = {entity: idx for idx, entity in enumerate(entities)}
    presence, richness = _availability_matrices(availability, entity_index, genus_index)

    parts: list[pd.DataFrame] = []
    for stratum in STRATA:
        counts = _representation_count_matrix(
            status_flora,
            island_index,
            genus_index,
            stratum=stratum,
            species_filter=species_filter,
        ).toarray()
        for source_mode in source_modes:
            assignment_matrix = _source_assignment_matrix(
                assignments,
                island_index,
                entity_index,
                source_mode=source_mode,
            )
            prevalence = (assignment_matrix @ presence).toarray().astype(np.int16)
            source_richness = (assignment_matrix @ richness).toarray().astype(np.float32)
            for matching in SOURCE_MATCHING:
                for minimum_repr in minimum_represented_genera:
                    rows: list[dict[str, Any]] = []
                    for island_idx, island_id in enumerate(islands):
                        result = compute_island_enrichment(
                            prevalence[island_idx],
                            source_richness[island_idx],
                            counts[island_idx],
                            position_array,
                            matching=matching,
                            minimum_represented_genera=minimum_repr,
                        )
                        if result is None:
                            continue
                        result.update(
                            {
                                "island_id": island_id,
                                "evidence_scope": evidence_scope,
                                "stratum": stratum,
                                "source_mode": source_mode,
                                "source_matching": matching,
                                "minimum_represented_genera": minimum_repr,
                            }
                        )
                        rows.append(result)
                    if rows:
                        parts.append(pd.DataFrame(rows))

    scores = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    source_counts = gift_source_functional_counts(matched, functional)
    manifest = {
        "evidence_scope": evidence_scope,
        "restricted_to_exact_si": restricted_to_exact_si,
        "raw_gift_availability": raw_gift_availability,
        "minimum_source_scored_species_per_genus": int(minimum_source_scored_species),
        "minimum_represented_genera": [int(x) for x in minimum_represented_genera],
        "n_functional_species": len(functional),
        **source_counts,
        "n_source_position_genera": len(positions),
        "n_source_availability_rows": len(availability),
    }
    return scores, manifest


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _fit_clustered(
    y: np.ndarray,
    X: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    n, p = X.shape
    cluster_labels = clusters.astype(str)
    unique_clusters = np.unique(cluster_labels)
    if n < max(20, p + 3) or len(unique_clusters) < 5:
        return pd.DataFrame(), np.empty((0, 0)), {
            "status": "insufficient_complete_rows",
            "n_rows": int(n),
            "n_clusters": len(unique_clusters),
        }
    bread = np.linalg.pinv(X.T @ X)
    beta = bread @ (X.T @ y)
    residual = y - X @ beta
    meat = np.zeros((p, p), dtype=float)
    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        score = X[mask].T @ residual[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if len(unique_clusters) > 1 and n > p:
        covariance *= (len(unique_clusters) / (len(unique_clusters) - 1.0)) * (
            (n - 1.0) / (n - p)
        )
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    rows = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        p_value = (
            math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")
        )
        rows.append(
            {
                "predictor": name,
                "estimate": float(estimate),
                "cluster_robust_se": float(stderr),
                "p_value": p_value,
            }
        )
    return pd.DataFrame(rows), covariance, {
        "status": "fit",
        "n_rows": int(n),
        "n_clusters": len(unique_clusters),
    }


def _covariates_with_realm(
    covariates: pd.DataFrame, realm_assignment: pd.DataFrame
) -> pd.DataFrame:
    return covariates.merge(
        realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates(
            "island_id"
        ),
        on="island_id",
        how="left",
        validate="one_to_one",
    )


def fit_context_slopes(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    cov = _covariates_with_realm(covariates, realm_assignment)
    family = [
        "evidence_scope",
        "stratum",
        "source_mode",
        "source_matching",
        "minimum_represented_genera",
    ]
    rows: list[dict[str, Any]] = []
    for key, subset in island_scores.groupby(family, sort=True):
        base_meta = dict(zip(family, key, strict=True))
        for context_layer, contexts in FOCAL_CONTEXTS.items():
            for context in contexts:
                work = subset.merge(
                    cov[["island_id", context_layer, cluster, geography, *baseline]],
                    on="island_id",
                    how="left",
                    validate="one_to_one",
                )
                work = work.loc[work[context_layer].astype(str).eq(context)].copy()
                for outcome in OUTCOMES:
                    complete = work.dropna(
                        subset=[outcome, cluster, geography, *baseline]
                    ).copy()
                    result = {
                        **base_meta,
                        "context_layer": context_layer,
                        "context": context,
                        "outcome": outcome,
                        "status": "insufficient_complete_rows",
                        "n_islands": len(complete),
                        "n_clusters": int(
                            complete[cluster].fillna("").astype(str).nunique()
                        ),
                        "distance_slope": np.nan,
                        "cluster_robust_se": np.nan,
                        "p_value": np.nan,
                    }
                    if len(complete) >= 20:
                        names = [
                            "intercept",
                            f"z_{geography}",
                            *[f"z_{x}" for x in baseline],
                        ]
                        X = np.column_stack(
                            [
                                np.ones(len(complete), dtype=float),
                                _standardize(complete[geography]),
                                *[_standardize(complete[x]) for x in baseline],
                            ]
                        )
                        coefficients, _, fit = _fit_clustered(
                            complete[outcome].to_numpy(float),
                            X,
                            names,
                            complete[cluster].astype(str).to_numpy(),
                        )
                        result["status"] = str(fit["status"])
                        result["n_clusters"] = int(fit["n_clusters"])
                        if not coefficients.empty:
                            row = coefficients.set_index("predictor").loc[
                                f"z_{geography}"
                            ]
                            result.update(
                                {
                                    "distance_slope": float(row["estimate"]),
                                    "cluster_robust_se": float(
                                        row["cluster_robust_se"]
                                    ),
                                    "p_value": float(row["p_value"]),
                                }
                            )
                    rows.append(result)
    slopes = pd.DataFrame(rows)
    q_family = [
        "evidence_scope",
        "stratum",
        "source_matching",
        "minimum_represented_genera",
        "context_layer",
        "context",
        "outcome",
    ]
    slopes["q_across_source_modes"] = slopes.groupby(
        q_family, group_keys=False
    )["p_value"].transform(_bh)
    slopes["positive_supported"] = (
        slopes["distance_slope"].gt(0)
        & slopes["q_across_source_modes"].le(0.05)
    ).fillna(False)
    slopes["negative_supported"] = (
        slopes["distance_slope"].lt(0)
        & slopes["q_across_source_modes"].le(0.05)
    ).fillna(False)
    return slopes


def fit_between_context_differences(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> pd.DataFrame:
    """Directly test primary broad distance-slope differences between contexts."""

    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    cov = _covariates_with_realm(covariates, realm_assignment)
    primary = island_scores.loc[
        island_scores["evidence_scope"].astype(str).eq("broad")
        & island_scores["source_matching"].astype(str).eq("prevalence_richness")
        & pd.to_numeric(
            island_scores["minimum_represented_genera"], errors="coerce"
        ).eq(5)
    ].copy()
    family = [
        "evidence_scope",
        "stratum",
        "source_mode",
        "source_matching",
        "minimum_represented_genera",
    ]
    rows: list[dict[str, Any]] = []
    for key, subset in primary.groupby(family, sort=True):
        base_meta = dict(zip(family, key, strict=True))
        for context_layer, (context_a, context_b) in FOCAL_CONTEXTS.items():
            work = subset.merge(
                cov[["island_id", context_layer, cluster, geography, *baseline]],
                on="island_id",
                how="left",
                validate="one_to_one",
            )
            work = work.loc[
                work[context_layer].astype(str).isin([context_a, context_b])
            ].copy()
            for outcome in OUTCOMES:
                complete = work.dropna(
                    subset=[outcome, cluster, geography, *baseline]
                ).copy()
                n_a = int(complete[context_layer].astype(str).eq(context_a).sum())
                n_b = int(complete[context_layer].astype(str).eq(context_b).sum())
                result = {
                    **base_meta,
                    "context_layer": context_layer,
                    "context_a": context_a,
                    "context_b": context_b,
                    "outcome": outcome,
                    "status": "insufficient_complete_rows",
                    "n_context_a": n_a,
                    "n_context_b": n_b,
                    "n_clusters": int(
                        complete[cluster].fillna("").astype(str).nunique()
                    ),
                    "context_a_distance_slope": np.nan,
                    "context_b_distance_slope": np.nan,
                    "slope_difference_b_minus_a": np.nan,
                    "cluster_robust_se": np.nan,
                    "p_value": np.nan,
                }
                if n_a < 20 or n_b < 20:
                    rows.append(result)
                    continue

                names: list[str] = []
                columns: list[np.ndarray] = []
                for context in (context_a, context_b):
                    mask = complete[context_layer].astype(str).eq(context).to_numpy()
                    indicator = mask.astype(float)
                    names.append(f"intercept[{context}]")
                    columns.append(indicator)
                    for predictor in baseline:
                        values = np.zeros(len(complete), dtype=float)
                        values[mask] = _standardize(complete.loc[mask, predictor])
                        names.append(f"z_{predictor}[{context}]")
                        columns.append(values)
                    values = np.zeros(len(complete), dtype=float)
                    values[mask] = _standardize(complete.loc[mask, geography])
                    names.append(f"z_{geography}[{context}]")
                    columns.append(values)

                coefficients, covariance, fit = _fit_clustered(
                    complete[outcome].to_numpy(float),
                    np.column_stack(columns),
                    names,
                    complete[cluster].astype(str).to_numpy(),
                )
                result["status"] = str(fit["status"])
                result["n_clusters"] = int(fit["n_clusters"])
                if coefficients.empty:
                    rows.append(result)
                    continue
                index = coefficients.set_index("predictor")
                name_a = f"z_{geography}[{context_a}]"
                name_b = f"z_{geography}[{context_b}]"
                pos_a = names.index(name_a)
                pos_b = names.index(name_b)
                beta_a = float(index.loc[name_a, "estimate"])
                beta_b = float(index.loc[name_b, "estimate"])
                difference = beta_b - beta_a
                variance = float(
                    covariance[pos_a, pos_a]
                    + covariance[pos_b, pos_b]
                    - 2.0 * covariance[pos_a, pos_b]
                )
                stderr = math.sqrt(max(0.0, variance))
                z = difference / stderr if stderr > 0 else float("nan")
                p_value = (
                    math.erfc(abs(z) / math.sqrt(2.0))
                    if math.isfinite(z)
                    else float("nan")
                )
                result.update(
                    {
                        "context_a_distance_slope": beta_a,
                        "context_b_distance_slope": beta_b,
                        "slope_difference_b_minus_a": difference,
                        "cluster_robust_se": stderr,
                        "p_value": p_value,
                    }
                )
                rows.append(result)

    between = pd.DataFrame(rows)
    q_family = [
        "evidence_scope",
        "stratum",
        "source_matching",
        "minimum_represented_genera",
        "context_layer",
        "outcome",
    ]
    between["q_across_source_modes"] = between.groupby(
        q_family, group_keys=False
    )["p_value"].transform(_bh)
    between["difference_supported"] = between[
        "q_across_source_modes"
    ].le(0.05).fillna(False)
    return between


def run(
    species_scores: pd.DataFrame,
    direct_species_scores: pd.DataFrame,
    trait_ledger: pd.DataFrame,
    direct_trait_ledger: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
) -> dict[str, Any]:
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]
    definitions = [
        {
            "evidence_scope": "broad",
            "species_scores": species_scores,
            "trait_ledger": None,
            "minimum_source_scored_species": 1,
            "minimum_represented_genera": [5, 10],
            "restricted_to_exact_si": False,
            "raw_gift_availability": True,
        },
        {
            "evidence_scope": "broad_source_min2",
            "species_scores": species_scores,
            "trait_ledger": None,
            "minimum_source_scored_species": 2,
            "minimum_represented_genera": [5],
            "restricted_to_exact_si": False,
            "raw_gift_availability": True,
        },
        {
            "evidence_scope": "exact_si_all",
            "species_scores": species_scores,
            "trait_ledger": trait_ledger,
            "minimum_source_scored_species": 1,
            "minimum_represented_genera": [3, 5],
            "restricted_to_exact_si": True,
            "raw_gift_availability": False,
        },
        {
            "evidence_scope": "exact_si_direct",
            "species_scores": direct_species_scores,
            "trait_ledger": direct_trait_ledger,
            "minimum_source_scored_species": 1,
            "minimum_represented_genera": [3],
            "restricted_to_exact_si": True,
            "raw_gift_availability": False,
        },
    ]
    direct_species = set(
        direct_species_scores["accepted_species"].dropna().astype(str)
    )
    all_species = set(species_scores["accepted_species"].dropna().astype(str))
    if not direct_species.issubset(all_species):
        raise ValueError(
            "direct-only syndrome species must be a subset of the all-analysis universe"
        )
    shared_matched = match_gift_species(
        gift_flora, species_scores[["accepted_species"]]
    )

    score_parts: list[pd.DataFrame] = []
    evidence_manifests: list[dict[str, Any]] = []
    for definition in definitions:
        scores_out, manifest = build_evidence_scores(
            evidence_scope=str(definition["evidence_scope"]),
            species_scores=definition["species_scores"],
            trait_ledger=definition["trait_ledger"],
            status_flora=status_flora,
            gift_flora=gift_flora,
            assignments=assignments,
            covariates=covariates,
            source_modes=source_modes,
            minimum_source_scored_species=int(
                definition["minimum_source_scored_species"]
            ),
            minimum_represented_genera=list(
                definition["minimum_represented_genera"]
            ),
            restricted_to_exact_si=bool(definition["restricted_to_exact_si"]),
            raw_gift_availability=bool(definition["raw_gift_availability"]),
            matched_gift=shared_matched,
        )
        score_parts.append(scores_out)
        evidence_manifests.append(manifest)
    island_scores = pd.concat(score_parts, ignore_index=True)
    slopes = fit_context_slopes(
        island_scores, covariates, realm_assignment, pattern_config
    )
    between = fit_between_context_differences(
        island_scores, covariates, realm_assignment, pattern_config
    )
    manifest = {
        "contract": "chapter1_pr138_lineage_representation_bridge_v1",
        "functional_position": "(-large_bee_like + generalized_accessible) / 2",
        "source_position_uses_gift_mainland_species_only": True,
        "source_assignment_is_frozen_and_outcome_blind": True,
        "primary_source_matching": "prevalence_richness",
        "primary_minimum_represented_genera": 5,
        "source_richness_bins": ["1", "2", "3-4", "5-8", "9+"],
        "outcomes": list(OUTCOMES),
        "entry_definition": (
            "one vote per represented genus after source-availability matching"
        ),
        "species_definition": (
            "one vote per represented island species using its genus source position"
        ),
        "loading_definition": "species_enrichment - entry_enrichment",
        "direct_between_context_tests": (
            "primary broad prevalence_richness minimum-five-genus design only"
        ),
        "evidence_manifests": evidence_manifests,
        "claim_boundary": (
            "Assembly diagnostic only. Enrichment does not identify pollinator loss, "
            "historical source ancestry, colonisation, extinction, diversification, "
            "or in-situ evolution."
        ),
    }
    return {
        "island_scores": island_scores,
        "slopes": slopes,
        "between": between,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-scores-csv", type=Path, required=True)
    parser.add_argument("--direct-species-scores-csv", type=Path, required=True)
    parser.add_argument("--trait-ledger-csv", type=Path, required=True)
    parser.add_argument("--direct-trait-ledger-csv", type=Path, required=True)
    parser.add_argument("--status-flora-csv", type=Path, required=True)
    parser.add_argument("--gift-flora-csv", type=Path, required=True)
    parser.add_argument("--source-assignments-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--source-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pattern_config = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    source_config = yaml.safe_load(args.source_config_path.read_text(encoding="utf-8"))
    outputs = run(
        pd.read_csv(args.species_scores_csv),
        pd.read_csv(args.direct_species_scores_csv),
        pd.read_csv(args.trait_ledger_csv),
        pd.read_csv(args.direct_trait_ledger_csv),
        pd.read_csv(args.status_flora_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        pattern_config,
        source_config,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs["island_scores"].to_csv(
        args.output_dir / "lineage_representation_island_scores.csv.gz", index=False
    )
    outputs["slopes"].to_csv(
        args.output_dir / "lineage_representation_context_slopes.csv", index=False
    )
    outputs["between"].to_csv(
        args.output_dir / "lineage_representation_between_context.csv", index=False
    )
    (args.output_dir / "lineage_representation_manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outputs["manifest"], indent=2))


if __name__ == "__main__":
    main()
