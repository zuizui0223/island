"""Direct source-matched lineage representation decomposition of current Chapter 1 H1.

This analysis reuses the frozen GIFT mainland source flora and island source assignments
from the historical lineage-representation bridge, but replaces the historical scalar
floral functional position with the seven current H1 binary responses.

For each H1 outcome, GIFT mainland-native scored species define a source genus position
(mean oriented H1 state). Island source-available genera are matched by frozen source
prevalence and source-richness classes. The analysis separates:

1. genus entry enrichment;
2. species-weighted enrichment; and
3. additional within-represented-genus loading.

The same exact source-backed native-nonendemic tropical flora and corrected 24 September
2026 geography are used throughout. This is an assembly representation diagnostic, not
a causal colonisation, extinction, pollinator-loss, or evolutionary analysis.
"""
from __future__ import annotations

import json
import math
import re
import unicodedata
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy import sparse

from island_v2.chapter1_all_data_probability import (
    _bh,
    _classify_signature,
    _truthy,
    run_probability_analysis,
)
from island_v2.chapter1_h3_observed_taxonomic_depth import (
    _assemble_cluster_covariance,
    _fit_ols_component,
    _z,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

SOURCE_MODES = ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10")
DECOMPOSITION_OUTCOMES = (
    "entry_enrichment",
    "species_enrichment",
    "loading_increment",
)


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


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if (
        not isinstance(cfg, dict)
        or cfg.get("contract") != "chapter1_current_h1_source_lineage_entry_v1"
    ):
        raise typer.BadParameter("unexpected current H1 source-lineage contract")
    return cfg


def build_outcome_states(
    state_audit: pd.DataFrame,
    probability_config: dict[str, Any],
) -> pd.DataFrame:
    required = {
        "accepted_species",
        "trait_name",
        "resolved_for_primary",
        "canonical_signature",
    }
    if missing := required - set(state_audit.columns):
        raise typer.BadParameter(f"state audit missing columns: {sorted(missing)}")
    audit = state_audit.loc[_truthy(state_audit["resolved_for_primary"])].copy()
    audit["accepted_species"] = audit["accepted_species"].astype(str)
    audit["trait_name"] = audit["trait_name"].astype(str)
    parts: list[pd.DataFrame] = []
    for outcome in probability_config["model_outcomes"]:
        spec = probability_config["broad_outcomes"][outcome]
        part = audit.loc[
            audit["trait_name"].eq(str(spec["trait_name"])),
            ["accepted_species", "canonical_signature"],
        ].drop_duplicates("accepted_species")
        positive = {str(x) for x in spec["positive_states"]}
        negative = {str(x) for x in spec["negative_states"]}
        part = part.copy()
        part["state"] = [
            _classify_signature(value, positive, negative)
            for value in part["canonical_signature"]
        ]
        part = part.dropna(subset=["state"])[["accepted_species", "state"]]
        part["outcome"] = str(outcome)
        part["genus"] = part["accepted_species"].map(_genus)
        part = part.loc[part["genus"].ne("")]
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def match_gift_species(
    gift_flora: pd.DataFrame,
    species_universe: pd.Series,
) -> pd.DataFrame:
    required = {"entity_ID", "work_species"}
    if missing := required - set(gift_flora.columns):
        raise typer.BadParameter(f"GIFT flora missing columns: {sorted(missing)}")

    species = pd.DataFrame(
        {"accepted_species": species_universe.dropna().astype(str).drop_duplicates()}
    )
    species["norm"] = species["accepted_species"].map(_normalise_name)
    species["binomial"] = species["accepted_species"].map(_binomial_key)

    exact_groups = species.groupby("norm")["accepted_species"].agg(
        lambda values: sorted(set(values))
    )
    exact = {
        key: values[0]
        for key, values in exact_groups.items()
        if key and len(values) == 1
    }
    binomial_groups = species.groupby("binomial")["accepted_species"].agg(
        lambda values: sorted(set(values))
    )
    binomial = {
        key: values[0]
        for key, values in binomial_groups.items()
        if key and len(values) == 1
    }

    names = gift_flora[["work_species"]].drop_duplicates().copy()
    names["norm"] = names["work_species"].map(_normalise_name)
    names["binomial"] = names["work_species"].map(_binomial_key)
    names["accepted_species"] = names["norm"].map(exact).fillna("")
    need_binomial = names["accepted_species"].eq("") & names["binomial"].ne("")
    names.loc[need_binomial, "accepted_species"] = (
        names.loc[need_binomial, "binomial"].map(binomial).fillna("")
    )
    matched_names = names.loc[
        names["accepted_species"].ne(""), ["work_species", "accepted_species"]
    ]
    return (
        gift_flora[["entity_ID", "work_species"]]
        .merge(matched_names, on="work_species", how="inner", validate="many_to_one")
        .drop_duplicates(["entity_ID", "accepted_species"])
    )


def build_source_genus_positions(
    matched_gift: pd.DataFrame,
    states: pd.DataFrame,
    *,
    minimum_scored_species: int,
) -> pd.DataFrame:
    source = (
        matched_gift.merge(
            states[["accepted_species", "outcome", "state", "genus"]],
            on="accepted_species",
            how="inner",
            validate="many_to_many",
        )
        .drop_duplicates(["accepted_species", "outcome"])
        .copy()
    )
    positions = source.groupby(["outcome", "genus"], as_index=False).agg(
        source_h1_position=("state", "mean"),
        n_source_scored_species=("accepted_species", "nunique"),
    )
    return positions.loc[
        positions["n_source_scored_species"].ge(int(minimum_scored_species))
    ].reset_index(drop=True)


def build_raw_source_availability(gift_flora: pd.DataFrame) -> pd.DataFrame:
    work = gift_flora[["entity_ID", "work_species"]].drop_duplicates().copy()
    work["genus"] = work["work_species"].map(_genus)
    work = work.loc[work["genus"].ne("")]
    return work.groupby(["entity_ID", "genus"], as_index=False).agg(
        source_species_richness=("work_species", "nunique")
    )


def build_island_genus_counts(
    status_flora: pd.DataFrame,
    tropical_islands: set[str],
) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "floristic_status"}
    if missing := required - set(status_flora.columns):
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    work = status_flora.loc[
        status_flora["floristic_status"].astype(str).eq("native_nonendemic"),
        ["island_id", "accepted_species"],
    ].drop_duplicates()
    work["island_id"] = work["island_id"].astype(str)
    work = work.loc[work["island_id"].isin(tropical_islands)]
    work["genus"] = work["accepted_species"].map(_genus)
    work = work.loc[work["genus"].ne("")]
    return work.groupby(["island_id", "genus"], as_index=False).agg(
        n_island_species=("accepted_species", "nunique")
    )


def _richness_bins_array(values: np.ndarray) -> np.ndarray:
    bins = np.zeros_like(values, dtype=np.int8)
    bins[(values > 1) & (values <= 2)] = 1
    bins[(values > 2) & (values <= 4)] = 2
    bins[(values > 4) & (values <= 8)] = 3
    bins[values > 8] = 4
    return bins


def compute_island_enrichment_arrays(
    prevalence: np.ndarray,
    source_richness: np.ndarray,
    island_species_counts: np.ndarray,
    positions: np.ndarray,
    *,
    minimum_represented_genera: int,
) -> dict[str, float | int] | None:
    candidate = prevalence > 0
    counts = island_species_counts.astype(float).copy()
    counts[~candidate] = 0.0
    represented = counts > 0
    n_represented = int(represented.sum())
    n_species = int(counts.sum())
    if n_represented < int(minimum_represented_genera) or n_species <= 0:
        return None

    observed_entry = float(np.mean(positions[represented]))
    observed_species = float(np.sum(positions * counts) / n_species)
    richness_bins = _richness_bins_array(source_richness)

    expected_entry_sum = 0.0
    expected_species_sum = 0.0
    classes = np.unique(
        np.column_stack([prevalence[candidate], richness_bins[candidate]]),
        axis=0,
    )
    for source_prevalence, richness_bin in classes:
        class_mask = (
            candidate
            & (prevalence == source_prevalence)
            & (richness_bins == richness_bin)
        )
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


def _matrix_inputs(
    *,
    positions: pd.DataFrame,
    source_availability: pd.DataFrame,
    island_counts: pd.DataFrame,
    assignments: pd.DataFrame,
    tropical_islands: set[str],
) -> tuple[
    list[str],
    dict[str, int],
    np.ndarray,
    dict[str, tuple[np.ndarray, np.ndarray]],
]:
    genera = sorted(positions["genus"].dropna().astype(str).unique())
    genus_index = {genus: idx for idx, genus in enumerate(genera)}
    islands = sorted(str(x) for x in tropical_islands)
    island_index = {island: idx for idx, island in enumerate(islands)}

    avail = source_availability.loc[
        source_availability["genus"].astype(str).isin(genus_index)
    ].copy()
    avail["entity_ID"] = pd.to_numeric(avail["entity_ID"], errors="coerce")
    avail = avail.dropna(subset=["entity_ID"])
    avail["entity_ID"] = avail["entity_ID"].astype(int)

    assign = assignments.loc[
        assignments["source_mode"].astype(str).isin(SOURCE_MODES)
        & assignments["island_id"].astype(str).isin(island_index),
        ["island_id", "source_mode", "entity_ID"],
    ].drop_duplicates().copy()
    assign["island_id"] = assign["island_id"].astype(str)
    assign["entity_ID"] = pd.to_numeric(assign["entity_ID"], errors="coerce")
    assign = assign.dropna(subset=["entity_ID"])
    assign["entity_ID"] = assign["entity_ID"].astype(int)

    entities = sorted(set(avail["entity_ID"]) | set(assign["entity_ID"]))
    entity_index = {entity: idx for idx, entity in enumerate(entities)}

    arows = avail["entity_ID"].map(entity_index).to_numpy(int)
    acols = avail["genus"].astype(str).map(genus_index).to_numpy(int)
    presence = sparse.csr_matrix(
        (np.ones(len(avail), dtype=float), (arows, acols)),
        shape=(len(entity_index), len(genus_index)),
    )
    richness = sparse.csr_matrix(
        (
            pd.to_numeric(avail["source_species_richness"], errors="coerce")
            .fillna(0)
            .to_numpy(float),
            (arows, acols),
        ),
        shape=presence.shape,
    )

    counts = island_counts.loc[
        island_counts["island_id"].astype(str).isin(island_index)
        & island_counts["genus"].astype(str).isin(genus_index)
    ].copy()
    crows = counts["island_id"].astype(str).map(island_index).to_numpy(int)
    ccols = counts["genus"].astype(str).map(genus_index).to_numpy(int)
    count_matrix = sparse.csr_matrix(
        (
            pd.to_numeric(counts["n_island_species"], errors="coerce")
            .fillna(0)
            .to_numpy(float),
            (crows, ccols),
        ),
        shape=(len(island_index), len(genus_index)),
    ).toarray()

    mode_matrices: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for source_mode in SOURCE_MODES:
        mode = assign.loc[assign["source_mode"].astype(str).eq(source_mode)].copy()
        rows = mode["island_id"].map(island_index).to_numpy(int)
        cols = mode["entity_ID"].map(entity_index).to_numpy(int)
        matrix = sparse.csr_matrix(
            (np.ones(len(mode), dtype=float), (rows, cols)),
            shape=(len(island_index), len(entity_index)),
        )
        prevalence = (matrix @ presence).toarray().astype(np.int16)
        source_richness = (matrix @ richness).toarray().astype(np.float32)
        mode_matrices[source_mode] = (prevalence, source_richness)
    return islands, genus_index, count_matrix, mode_matrices


def build_island_scores(
    *,
    states: pd.DataFrame,
    positions: pd.DataFrame,
    source_availability: pd.DataFrame,
    island_counts: pd.DataFrame,
    assignments: pd.DataFrame,
    tropical_islands: set[str],
    evidence_scope: str,
    minimum_represented_genera: int,
) -> pd.DataFrame:
    islands, genus_index, count_matrix, mode_matrices = _matrix_inputs(
        positions=positions,
        source_availability=source_availability,
        island_counts=island_counts,
        assignments=assignments,
        tropical_islands=tropical_islands,
    )

    parts: list[pd.DataFrame] = []
    for outcome in sorted(states["outcome"].astype(str).unique()):
        pos = positions.loc[positions["outcome"].astype(str).eq(outcome)].copy()
        pos = pos.loc[pos["genus"].astype(str).isin(genus_index)]
        if pos.empty:
            continue
        pos = pos.drop_duplicates("genus").copy()
        pos["genus_idx"] = pos["genus"].astype(str).map(genus_index)
        pos = pos.sort_values("genus_idx")
        indices = pos["genus_idx"].to_numpy(int)
        position_array = pos["source_h1_position"].to_numpy(float)

        for source_mode in SOURCE_MODES:
            prevalence_full, richness_full = mode_matrices[source_mode]
            rows: list[dict[str, Any]] = []
            for island_idx, island_id in enumerate(islands):
                result = compute_island_enrichment_arrays(
                    prevalence_full[island_idx, indices],
                    richness_full[island_idx, indices],
                    count_matrix[island_idx, indices],
                    position_array,
                    minimum_represented_genera=minimum_represented_genera,
                )
                if result is None:
                    continue
                rows.append(
                    {
                        "island_id": str(island_id),
                        "evidence_scope": evidence_scope,
                        "outcome": outcome,
                        "source_mode": source_mode,
                        **result,
                    }
                )
            if rows:
                parts.append(pd.DataFrame(rows))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

def build_parent_counts(
    *,
    states: pd.DataFrame,
    positions: pd.DataFrame,
    status_flora: pd.DataFrame,
    score_support: pd.DataFrame,
    source_mode: str,
) -> pd.DataFrame:
    flora = status_flora.loc[
        status_flora["floristic_status"].astype(str).eq("native_nonendemic"),
        ["island_id", "accepted_species"],
    ].drop_duplicates().copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    support = score_support.loc[
        score_support["source_mode"].astype(str).eq(source_mode),
        ["island_id", "outcome"],
    ].drop_duplicates()
    pieces: list[pd.DataFrame] = []
    for outcome, outcome_states in states.groupby("outcome", sort=False):
        genera = set(
            positions.loc[
                positions["outcome"].astype(str).eq(str(outcome)), "genus"
            ].astype(str)
        )
        if not genera:
            continue
        st = outcome_states[["accepted_species", "state", "genus"]].copy()
        st = st.loc[st["genus"].astype(str).isin(genera)]
        joined = flora.merge(
            st[["accepted_species", "state"]],
            on="accepted_species",
            how="inner",
            validate="many_to_one",
        )
        allowed = support.loc[support["outcome"].astype(str).eq(str(outcome)), "island_id"]
        joined = joined.loc[joined["island_id"].isin(set(allowed.astype(str)))]
        if joined.empty:
            continue
        summary = joined.groupby("island_id", as_index=False).agg(
            successes=("state", "sum"),
            trials=("state", "size"),
        )
        summary["outcome"] = str(outcome)
        summary["stratum"] = str(source_mode)
        pieces.append(summary)
    return pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()


def fit_parent_h1_gate(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    *,
    source_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if counts.empty:
        return pd.DataFrame(), pd.DataFrame(
            [
                {
                    "source_mode": source_mode,
                    "status": "not_testable",
                    "p_value": np.nan,
                    "q_value": np.nan,
                    "vector_supported": False,
                }
            ]
        )
    cfg = deepcopy(probability_config)
    cfg["contexts"] = ["tropical"]
    cfg["between_contexts"] = []
    cfg["strata"] = [source_mode]
    slopes, _, omnibus, _ = run_probability_analysis(counts, covariates, cfg)
    if not slopes.empty:
        slopes.insert(0, "source_mode", source_mode)
    if not omnibus.empty:
        omnibus.insert(0, "source_mode", source_mode)
    return slopes, omnibus


def fit_joint_source_vector(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    source_mode: str,
    response: str,
    minimum_islands_per_outcome: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    needed_cov = [
        "island_id",
        "spatial_block",
        "log_distance_to_continent_km",
        "log_island_area_km2",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
    ]
    data = scores.loc[
        scores["source_mode"].astype(str).eq(source_mode)
    ].merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    predictors = [
        "log_island_area_km2",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
        "log_distance_to_continent_km",
    ]
    fits: list[dict[str, Any]] = []
    indices: list[int] = []
    rows: list[dict[str, Any]] = []
    offset = 0
    for outcome in sorted(data["outcome"].astype(str).unique()):
        part = data.loc[data["outcome"].astype(str).eq(outcome)].copy()
        part = part.dropna(subset=[response, "spatial_block", *predictors])
        if len(part) < int(minimum_islands_per_outcome):
            continue
        columns = [np.ones(len(part), dtype=float)]
        names = [f"{outcome}:intercept"]
        try:
            for predictor in predictors:
                columns.append(_z(part[predictor]))
                names.append(f"{outcome}:z_{predictor}")
        except ValueError:
            continue
        distance_name = f"{outcome}:z_log_distance_to_continent_km"
        fit = _fit_ols_component(
            part[response].to_numpy(float),
            np.column_stack(columns),
            names,
            part["spatial_block"].astype(str).to_numpy(),
        )
        fits.append(fit)
        indices.append(offset + names.index(distance_name))
        rows.append(
            {
                "source_mode": source_mode,
                "response": response,
                "outcome": outcome,
                "n_islands": int(len(part)),
                "n_clusters": int(part["spatial_block"].nunique()),
            }
        )
        offset += len(names)

    if len(fits) < 2:
        return pd.DataFrame(rows), {
            "source_mode": source_mode,
            "response": response,
            "status": "not_testable",
            "n_retained_outcomes": len(fits),
            "retained_outcomes": "|".join(row["outcome"] for row in rows),
            "p_value": np.nan,
            "vector_norm": np.nan,
            "mean_slope": np.nan,
            "n_positive_slopes": 0,
        }

    theta, covariance, _ = _assemble_cluster_covariance(fits)
    vector = theta[indices]
    vector_cov = covariance[np.ix_(indices, indices)]
    stderr = np.sqrt(np.clip(np.diag(vector_cov), 0.0, None))
    for row, estimate, se in zip(rows, vector, stderr, strict=True):
        z = float(estimate / se) if se > 0 else float("nan")
        row.update(
            {
                "estimate": float(estimate),
                "cluster_robust_se": float(se),
                "p_value": (
                    math.erfc(abs(z) / math.sqrt(2.0))
                    if math.isfinite(z)
                    else np.nan
                ),
            }
        )
    rank = int(np.linalg.matrix_rank(vector_cov))
    statistic = (
        float(vector @ np.linalg.pinv(vector_cov) @ vector)
        if rank > 0
        else float("nan")
    )
    from island_v2.chapter1_all_data_probability import _chi_square_sf_integer_df

    return pd.DataFrame(rows), {
        "source_mode": source_mode,
        "response": response,
        "status": "fit",
        "n_retained_outcomes": int(len(vector)),
        "retained_outcomes": "|".join(row["outcome"] for row in rows),
        "joint_wald_chisq": statistic,
        "joint_df": rank,
        "p_value": _chi_square_sf_integer_df(statistic, rank),
        "vector_norm": float(np.linalg.norm(vector)),
        "mean_slope": float(np.mean(vector)),
        "n_positive_slopes": int(np.sum(vector > 0)),
        "n_negative_slopes": int(np.sum(vector < 0)),
    }


def _summarize_route(
    omnibus: pd.DataFrame,
    parent_gate: pd.DataFrame,
    *,
    alpha: float,
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    parent = parent_gate.copy()
    if not parent.empty:
        parent["q_across_source_modes"] = _bh(parent["p_value"])
    for response in DECOMPOSITION_OUTCOMES:
        part = omnibus.loc[omnibus["response"].eq(response)].copy()
        part["q_across_source_modes"] = _bh(part["p_value"])
        result[response] = {
            "all_four_modes_testable": bool(len(part) == 4 and part["status"].eq("fit").all()),
            "all_four_modes_fdr_supported": bool(
                len(part) == 4 and part["q_across_source_modes"].le(alpha).all()
            ),
            "all_four_modes_positive_mean": bool(
                len(part) == 4 and part["mean_slope"].gt(0).all()
            ),
            "min_positive_atomic_slopes": (
                int(part["n_positive_slopes"].min()) if len(part) else 0
            ),
            "max_q": (
                float(part["q_across_source_modes"].max()) if len(part) else np.nan
            ),
            "min_mean_slope": (
                float(part["mean_slope"].min()) if len(part) else np.nan
            ),
            "max_mean_slope": (
                float(part["mean_slope"].max()) if len(part) else np.nan
            ),
        }
    parent_supported = bool(
        len(parent) == 4
        and parent.get("q_across_source_modes", pd.Series(dtype=float)).le(alpha).all()
    )
    entry = result["entry_enrichment"]
    loading = result["loading_increment"]
    result["decision"] = {
        "all_four_parent_H1_gates_supported": parent_supported,
        "genus_entry_route_supported": bool(
            parent_supported
            and entry["all_four_modes_fdr_supported"]
            and entry["all_four_modes_positive_mean"]
        ),
        "additional_within_genus_loading_supported": bool(
            loading["all_four_modes_fdr_supported"]
            and loading["all_four_modes_positive_mean"]
        ),
    }
    return result


def run_scope(
    *,
    evidence_scope: str,
    state_audit: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    matched_gift: pd.DataFrame,
    source_availability: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    probability_config: dict[str, Any],
    analysis_config: dict[str, Any],
) -> dict[str, Any]:
    states = build_outcome_states(state_audit, probability_config)
    positions = build_source_genus_positions(
        matched_gift,
        states,
        minimum_scored_species=int(
            analysis_config["source_pool"]["minimum_scored_source_species_per_genus"]
        ),
    )
    tropical_islands = set(
        covariates.loc[
            covariates["analysis_regime"].astype(str).eq(
                str(analysis_config["scope"]["context"])
            ),
            "island_id",
        ].astype(str)
    )
    island_counts = build_island_genus_counts(status_flora, tropical_islands)
    scores = build_island_scores(
        states=states,
        positions=positions,
        source_availability=source_availability,
        island_counts=island_counts,
        assignments=assignments,
        tropical_islands=tropical_islands,
        evidence_scope=evidence_scope,
        minimum_represented_genera=int(
            analysis_config["source_pool"]["minimum_represented_genera"]
        ),
    )

    gate_slope_parts: list[pd.DataFrame] = []
    gate_rows: list[pd.DataFrame] = []
    atomic_parts: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    for source_mode in SOURCE_MODES:
        parent_counts = build_parent_counts(
            states=states,
            positions=positions,
            status_flora=status_flora,
            score_support=scores,
            source_mode=source_mode,
        )
        gate_slopes, gate = fit_parent_h1_gate(
            parent_counts,
            covariates,
            probability_config,
            source_mode=source_mode,
        )
        if not gate_slopes.empty:
            gate_slope_parts.append(gate_slopes)
        if not gate.empty:
            gate_rows.append(gate)

        for response in DECOMPOSITION_OUTCOMES:
            atomic, omnibus = fit_joint_source_vector(
                scores,
                covariates,
                source_mode=source_mode,
                response=response,
                minimum_islands_per_outcome=int(
                    probability_config["minimum_islands_per_outcome"]
                ),
            )
            if not atomic.empty:
                atomic.insert(0, "evidence_scope", evidence_scope)
                atomic_parts.append(atomic)
            omnibus["evidence_scope"] = evidence_scope
            omnibus_rows.append(omnibus)

    gate_table = pd.concat(gate_rows, ignore_index=True) if gate_rows else pd.DataFrame()
    if not gate_table.empty:
        gate_table["q_across_source_modes"] = _bh(gate_table["p_value"])
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["q_across_source_modes"] = np.nan
    for response, index in omnibus.groupby("response").groups.items():
        omnibus.loc[index, "q_across_source_modes"] = _bh(omnibus.loc[index, "p_value"])

    route = _summarize_route(
        omnibus,
        gate_table,
        alpha=float(analysis_config["multiplicity"]["alpha"]),
    )
    manifest = {
        "evidence_scope": evidence_scope,
        "n_state_species": int(states["accepted_species"].nunique()),
        "n_source_genus_positions": int(len(positions)),
        "n_source_scored_genera": int(positions["genus"].nunique()),
        "n_island_score_rows": int(len(scores)),
        "n_islands_any_score": int(scores["island_id"].nunique()) if not scores.empty else 0,
        "route_summary": route,
    }
    return {
        "states": states,
        "positions": positions,
        "scores": scores,
        "parent_gate_slopes": (
            pd.concat(gate_slope_parts, ignore_index=True)
            if gate_slope_parts
            else pd.DataFrame()
        ),
        "parent_gate": gate_table,
        "atomic_slopes": (
            pd.concat(atomic_parts, ignore_index=True)
            if atomic_parts
            else pd.DataFrame()
        ),
        "omnibus": omnibus,
        "manifest": manifest,
    }


@app.command()
def main(
    state_audit_all_csv: Path = typer.Option(..., exists=True),
    state_audit_direct_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    gift_flora_csv: Path = typer.Option(..., exists=True),
    source_assignments_csv: Path = typer.Option(..., exists=True),
    corrected_covariates_csv: Path = typer.Option(..., exists=True),
    probability_config_path: Path = typer.Option(..., exists=True),
    analysis_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    probability_config = yaml.safe_load(
        probability_config_path.read_text(encoding="utf-8")
    )
    analysis_config = load_config(analysis_config_path)
    status_flora = pd.read_csv(status_flora_csv)
    gift_flora = pd.read_csv(gift_flora_csv)
    assignments = pd.read_csv(source_assignments_csv)
    covariates = pd.read_csv(corrected_covariates_csv)
    all_audit = pd.read_csv(state_audit_all_csv)
    direct_audit = pd.read_csv(state_audit_direct_csv)

    all_states = build_outcome_states(all_audit, probability_config)
    matched_gift = match_gift_species(
        gift_flora,
        all_states["accepted_species"],
    )
    source_availability = build_raw_source_availability(gift_flora)

    outputs: dict[str, Any] = {}
    for label, audit in [
        ("all_analysis_eligible", all_audit),
        ("direct_only", direct_audit),
    ]:
        outputs[label] = run_scope(
            evidence_scope=label,
            state_audit=audit,
            status_flora=status_flora,
            gift_flora=gift_flora,
            matched_gift=matched_gift,
            source_availability=source_availability,
            assignments=assignments,
            covariates=covariates,
            probability_config=probability_config,
            analysis_config=analysis_config,
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    compact: list[dict[str, Any]] = []
    for label, result in outputs.items():
        scope_dir = output_dir / label
        scope_dir.mkdir(parents=True, exist_ok=True)
        for key in (
            "positions",
            "scores",
            "parent_gate_slopes",
            "parent_gate",
            "atomic_slopes",
            "omnibus",
        ):
            frame = result[key]
            frame.to_csv(scope_dir / f"{key}.csv", index=False)
        (scope_dir / "manifest.json").write_text(
            json.dumps(result["manifest"], indent=2) + "\n",
            encoding="utf-8",
        )
        route = result["manifest"]["route_summary"]
        for response in DECOMPOSITION_OUTCOMES:
            compact.append(
                {
                    "evidence_scope": label,
                    "response": response,
                    **route[response],
                }
            )
        compact.append(
            {
                "evidence_scope": label,
                "response": "decision",
                **route["decision"],
            }
        )
    pd.DataFrame(compact).to_csv(output_dir / "route_summary.csv", index=False)

    result_manifest = {
        "contract": analysis_config["contract"],
        "source_matching_reused": True,
        "current_H1_trait_recoding_reused": True,
        "context": analysis_config["scope"]["context"],
        "floristic_status": analysis_config["scope"]["floristic_status"],
        "source_modes": list(SOURCE_MODES),
        "claim_boundary": (
            "Direct representation decomposition of current H1. Genus entry can reflect "
            "arrival, establishment, persistence, extinction, habitat filtering or biotic "
            "interactions; loading does not identify within-lineage evolution."
        ),
        "scope_manifests": {
            label: outputs[label]["manifest"] for label in outputs
        },
    }
    (output_dir / "RESULT.json").write_text(
        json.dumps(result_manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(pd.DataFrame(compact).to_csv(index=False))


if __name__ == "__main__":
    app()
