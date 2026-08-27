"""Primary source-matched lineage-representation bridge for Chapter 1 -> izu-core."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.chapter1_pr138_si_genus_source_null import exact_si_species
from island_v2.chapter1_pr138_source_pool_sensitivity import match_gift_species_to_frozen_scores
from island_v2.chapter1_pr138_syndrome_analysis import _bh

AXES = ("large_bee_like", "generalized_accessible")
STRATA = ("all_native", "native_nonendemic")
CONTEXTS = {
    "analysis_regime": ("northern_midlatitude", "tropical"),
    "biogeographic_realm": ("Palearctic", "Neotropical"),
}
BASELINE = (
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
)
GEOGRAPHY = "log_distance_to_continent_km"
CLUSTER = "spatial_block"


@dataclass(frozen=True)
class SourceMatrices:
    island_ids: np.ndarray
    genera: np.ndarray
    n_regions: sparse.csr_matrix
    supply: sparse.csr_matrix
    position_num: sparse.csr_matrix


def build_species_positions(
    scores: pd.DataFrame, eligible: set[str] | None = None
) -> pd.DataFrame:
    work = scores.loc[
        scores["syndrome"].astype(str).isin(AXES),
        ["accepted_species", "syndrome", "syndrome_concordance"],
    ].copy()
    work["syndrome_concordance"] = pd.to_numeric(
        work["syndrome_concordance"], errors="coerce"
    )
    work = work.dropna(subset=["syndrome_concordance"])
    if eligible is not None:
        work = work.loc[work["accepted_species"].astype(str).isin(eligible)].copy()
    pivot = work.pivot_table(
        index="accepted_species",
        columns="syndrome",
        values="syndrome_concordance",
        aggfunc="first",
    )
    if any(axis not in pivot for axis in AXES):
        return pd.DataFrame(
            columns=["accepted_species", "genus", "functional_position"]
        )
    pivot = pivot.dropna(subset=list(AXES)).copy()
    pivot["functional_position"] = (
        -pivot["large_bee_like"] + pivot["generalized_accessible"]
    ) / 2.0
    out = pivot.reset_index()[["accepted_species", "functional_position"]]
    out["genus"] = out["accepted_species"].astype(str).str.split().str[0]
    return out.loc[out["genus"].ne("")].reset_index(drop=True)


def match_gift(
    gift: pd.DataFrame, positions: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    fake = positions[["accepted_species", "functional_position"]].rename(
        columns={"functional_position": "syndrome_concordance"}
    ).copy()
    fake["syndrome"] = "functional_position"
    matched, audit = match_gift_species_to_frozen_scores(gift, fake)
    matched = matched[
        ["entity_ID", "accepted_species", "syndrome_concordance"]
    ].drop_duplicates(["entity_ID", "accepted_species"])
    matched = matched.rename(
        columns={"syndrome_concordance": "functional_position"}
    )
    matched["genus"] = matched["accepted_species"].astype(str).str.split().str[0]
    return matched, audit


def build_source_matrices(
    matched: pd.DataFrame,
    assignments: pd.DataFrame,
    island_ids: np.ndarray,
) -> dict[str, SourceMatrices]:
    region_genus = matched.groupby(["entity_ID", "genus"], as_index=False).agg(
        n_region_species=("accepted_species", "nunique"),
        position_num=("functional_position", "sum"),
    )
    genera = np.array(
        sorted(region_genus["genus"].astype(str).unique()), dtype=object
    )
    entities = np.array(
        sorted(region_genus["entity_ID"].astype(int).unique()), dtype=int
    )
    gmap = pd.Series(np.arange(len(genera)), index=genera)
    emap = pd.Series(np.arange(len(entities)), index=entities)
    r = region_genus["entity_ID"].astype(int).map(emap).to_numpy(int)
    c = region_genus["genus"].astype(str).map(gmap).to_numpy(int)
    shape = (len(entities), len(genera))
    incidence = sparse.csr_matrix((np.ones(len(region_genus)), (r, c)), shape=shape)
    richness = sparse.csr_matrix(
        (region_genus["n_region_species"].to_numpy(float), (r, c)), shape=shape
    )
    position = sparse.csr_matrix(
        (region_genus["position_num"].to_numpy(float), (r, c)), shape=shape
    )
    imap = pd.Series(np.arange(len(island_ids)), index=island_ids.astype(str))
    out: dict[str, SourceMatrices] = {}
    for mode, frame in assignments.groupby("source_mode", sort=True):
        frame = frame[["island_id", "entity_ID"]].drop_duplicates()
        frame = frame.loc[
            frame["island_id"].astype(str).isin(imap.index)
            & frame["entity_ID"].astype(int).isin(emap.index)
        ].copy()
        ar = frame["island_id"].astype(str).map(imap).to_numpy(int)
        ac = frame["entity_ID"].astype(int).map(emap).to_numpy(int)
        assignment = sparse.csr_matrix(
            (np.ones(len(frame)), (ar, ac)), shape=(len(island_ids), len(entities))
        )
        out[str(mode)] = SourceMatrices(
            island_ids,
            genera,
            (assignment @ incidence).tocsr(),
            (assignment @ richness).tocsr(),
            (assignment @ position).tocsr(),
        )
    return out


def island_count_matrix(
    status: pd.DataFrame,
    island_ids: np.ndarray,
    genera: np.ndarray,
    stratum: str,
    eligible: set[str] | None,
) -> sparse.csr_matrix:
    if stratum == "all_native":
        mask = status["origin_status"].astype(str).eq("native")
    else:
        mask = status["floristic_status"].astype(str).eq("native_nonendemic")
    work = status.loc[mask, ["island_id", "accepted_species"]].drop_duplicates().copy()
    if eligible is not None:
        work = work.loc[work["accepted_species"].astype(str).isin(eligible)].copy()
    work["genus"] = work["accepted_species"].astype(str).str.split().str[0]
    imap = pd.Series(np.arange(len(island_ids)), index=island_ids.astype(str))
    gmap = pd.Series(np.arange(len(genera)), index=genera.astype(str))
    work = work.loc[
        work["island_id"].astype(str).isin(imap.index)
        & work["genus"].isin(gmap.index)
    ]
    grouped = work.groupby(["island_id", "genus"], as_index=False).agg(
        n=("accepted_species", "nunique")
    )
    r = grouped["island_id"].astype(str).map(imap).to_numpy(int)
    c = grouped["genus"].map(gmap).to_numpy(int)
    return sparse.csr_matrix(
        (grouped["n"].to_numpy(float), (r, c)),
        shape=(len(island_ids), len(genera)),
    )


def cell_residual(
    position: np.ndarray, n_regions: np.ndarray, supply: np.ndarray
) -> np.ndarray:
    richness_class = np.floor(
        np.log2(np.clip(supply, 1.0, None))
    ).astype(np.int64)
    keys = n_regions.astype(np.int64) * 100 + richness_class
    _, inverse = np.unique(keys, return_inverse=True)
    means = np.bincount(inverse, weights=position) / np.bincount(inverse)
    return position - means[inverse]


def build_responses(
    source: SourceMatrices,
    counts: sparse.csr_matrix,
    scope: str,
    stratum: str,
    mode: str,
    min_genera: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for i, island_id in enumerate(source.island_ids):
        region_row = source.n_regions.getrow(i)
        if region_row.nnz == 0:
            continue
        idx = region_row.indices
        n_regions = region_row.data
        supply = source.supply[i, idx].toarray().ravel()
        numerator = source.position_num[i, idx].toarray().ravel()
        ok = np.isfinite(supply) & np.isfinite(numerator) & (supply > 0)
        idx = idx[ok]
        n_regions = n_regions[ok]
        supply = supply[ok]
        numerator = numerator[ok]
        if len(idx) == 0:
            continue
        loading = counts[i, idx].toarray().ravel()
        represented = loading > 0
        n_represented = int(represented.sum())
        if n_represented < min_genera:
            continue
        residual = cell_residual(numerator / supply, n_regions, supply)
        entry = float(np.mean(residual[represented]))
        species = float(
            np.average(residual[represented], weights=loading[represented])
        )
        rows.append(
            {
                "evidence_scope": scope,
                "stratum": stratum,
                "source_mode": mode,
                "source_matching": "prevalence_richness",
                "minimum_represented_genera": min_genera,
                "island_id": str(island_id),
                "entry_enrichment": entry,
                "species_enrichment": species,
                "loading_increment": species - entry,
                "n_represented_genera": n_represented,
                "n_source_candidate_genera": int(len(idx)),
            }
        )
    return pd.DataFrame(rows)


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant predictor")
    return (x - float(np.mean(x))) / sd


def fit_slope(
    response: pd.DataFrame,
    covariates: pd.DataFrame,
    layer: str,
    context: str,
    outcome: str,
) -> dict[str, Any]:
    columns = ["island_id", layer, GEOGRAPHY, *BASELINE, CLUSTER]
    work = response[["island_id", outcome]].merge(
        covariates[columns].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    work = work.loc[work[layer].astype(str).eq(context)].dropna(
        subset=[outcome, GEOGRAPHY, *BASELINE, CLUSTER]
    ).copy()
    names = ["intercept", f"z_{GEOGRAPHY}", *[f"z_{x}" for x in BASELINE]]
    design = np.column_stack(
        [np.ones(len(work)), _z(work[GEOGRAPHY]), *[_z(work[x]) for x in BASELINE]]
    )
    coefficients, _, fit = _fit_weighted_clustered_design(
        work[outcome].to_numpy(float),
        np.ones(len(work)),
        design,
        names,
        work[CLUSTER].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return {
            "status": str(fit.get("status", "fit_failed")),
            "n_islands": int(work["island_id"].nunique()),
        }
    row = coefficients.set_index("predictor").loc[f"z_{GEOGRAPHY}"]
    return {
        "status": "fit",
        "n_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "distance_slope": float(row["estimate"]),
        "cluster_robust_se": float(row["cluster_robust_se"]),
        "p_value": float(row["p_value"]),
    }


def fit_between(
    response: pd.DataFrame,
    covariates: pd.DataFrame,
    layer: str,
    context_a: str,
    context_b: str,
    outcome: str,
) -> dict[str, Any]:
    columns = ["island_id", layer, GEOGRAPHY, *BASELINE, CLUSTER]
    work = response[["island_id", outcome]].merge(
        covariates[columns].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    work = work.loc[
        work[layer].astype(str).isin([context_a, context_b])
    ].dropna(subset=[outcome, GEOGRAPHY, *BASELINE, CLUSTER]).copy()
    parts: list[pd.DataFrame] = []
    for context in (context_a, context_b):
        part = work.loc[work[layer].astype(str).eq(context)].copy()
        part["z_distance"] = _z(part[GEOGRAPHY])
        for baseline in BASELINE:
            part[f"z_{baseline}"] = _z(part[baseline])
        parts.append(part)
    work = pd.concat(parts, ignore_index=True)
    design_columns: list[np.ndarray] = []
    names: list[str] = []
    for context in (context_a, context_b):
        indicator = work[layer].astype(str).eq(context).to_numpy(float)
        names += [f"intercept[{context}]", f"distance[{context}]"]
        design_columns += [indicator, indicator * work["z_distance"].to_numpy(float)]
        for baseline in BASELINE:
            names.append(f"{baseline}[{context}]")
            design_columns.append(
                indicator * work[f"z_{baseline}"].to_numpy(float)
            )
    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work[outcome].to_numpy(float),
        np.ones(len(work)),
        np.column_stack(design_columns),
        names,
        work[CLUSTER].astype(str).to_numpy(),
    )
    if coefficients.empty:
        return {"status": str(fit.get("status", "fit_failed"))}
    index = coefficients.set_index("predictor")
    ia = names.index(f"distance[{context_a}]")
    ib = names.index(f"distance[{context_b}]")
    slope_a = float(index.loc[f"distance[{context_a}]", "estimate"])
    slope_b = float(index.loc[f"distance[{context_b}]", "estimate"])
    difference = slope_b - slope_a
    se = math.sqrt(
        max(
            0.0,
            float(
                covariance[ib, ib]
                + covariance[ia, ia]
                - 2.0 * covariance[ib, ia]
            ),
        )
    )
    p_value = (
        math.erfc(abs(difference / se) / math.sqrt(2.0))
        if se > 0
        else float("nan")
    )
    return {
        "status": "fit",
        "n_context_a": int(parts[0]["island_id"].nunique()),
        "n_context_b": int(parts[1]["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "context_a_distance_slope": slope_a,
        "context_b_distance_slope": slope_b,
        "slope_difference_b_minus_a": difference,
        "cluster_robust_se": se,
        "p_value": p_value,
    }


def run_scope(
    scope: str,
    all_scores: pd.DataFrame,
    direct_scores: pd.DataFrame,
    all_ledger: pd.DataFrame,
    direct_ledger: pd.DataFrame,
    status: pd.DataFrame,
    gift: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm: pd.DataFrame,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    all_si = exact_si_species(all_ledger)
    direct_si = exact_si_species(direct_ledger)
    if scope == "broad":
        positions = build_species_positions(all_scores)
        eligible = None
        min_genera = 5
    elif scope == "exact_si_all":
        positions = build_species_positions(all_scores, all_si)
        eligible = all_si
        min_genera = 3
    elif scope == "exact_si_direct":
        positions = build_species_positions(direct_scores, direct_si)
        eligible = direct_si
        min_genera = 3
    else:
        raise ValueError(f"unsupported scope: {scope}")

    covariates = covariates.merge(
        realm[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    focal = covariates["analysis_regime"].astype(str).isin(
        CONTEXTS["analysis_regime"]
    ) | covariates["biogeographic_realm"].astype(str).isin(
        CONTEXTS["biogeographic_realm"]
    )
    island_ids = (
        covariates.loc[focal, "island_id"]
        .astype(str)
        .drop_duplicates()
        .to_numpy(object)
    )
    matched, match_audit = match_gift(gift, positions)
    source = build_source_matrices(matched, assignments, island_ids)
    parts: list[pd.DataFrame] = []
    for stratum in STRATA:
        exemplar = next(iter(source.values()))
        counts = island_count_matrix(
            status, island_ids, exemplar.genera, stratum, eligible
        )
        for mode, matrices in source.items():
            parts.append(
                build_responses(
                    matrices, counts, scope, stratum, mode, min_genera
                )
            )
    responses = pd.concat(parts, ignore_index=True)
    outcomes = ("entry_enrichment", "species_enrichment", "loading_increment")
    slope_rows: list[dict[str, Any]] = []
    for (stratum, mode), frame in responses.groupby(
        ["stratum", "source_mode"], sort=True
    ):
        for layer, contexts in CONTEXTS.items():
            for context in contexts:
                for outcome in outcomes:
                    slope_rows.append(
                        {
                            "evidence_scope": scope,
                            "stratum": stratum,
                            "source_mode": mode,
                            "source_matching": "prevalence_richness",
                            "minimum_represented_genera": min_genera,
                            "context_layer": layer,
                            "context": context,
                            "outcome": outcome,
                            **fit_slope(frame, covariates, layer, context, outcome),
                        }
                    )
    slopes = pd.DataFrame(slope_rows)
    family = ["stratum", "context_layer", "context", "outcome"]
    slopes["q_across_source_modes"] = slopes.groupby(
        family, group_keys=False
    )["p_value"].transform(_bh)
    slopes["positive_supported"] = (
        slopes["q_across_source_modes"].le(0.05)
        & slopes["distance_slope"].gt(0)
    )

    between_rows: list[dict[str, Any]] = []
    if scope == "broad":
        for (stratum, mode), frame in responses.groupby(
            ["stratum", "source_mode"], sort=True
        ):
            for layer, contexts in CONTEXTS.items():
                for outcome in outcomes:
                    between_rows.append(
                        {
                            "evidence_scope": scope,
                            "stratum": stratum,
                            "source_mode": mode,
                            "source_matching": "prevalence_richness",
                            "minimum_represented_genera": min_genera,
                            "context_layer": layer,
                            "context_a": contexts[0],
                            "context_b": contexts[1],
                            "outcome": outcome,
                            **fit_between(
                                frame,
                                covariates,
                                layer,
                                contexts[0],
                                contexts[1],
                                outcome,
                            ),
                        }
                    )
    between = pd.DataFrame(between_rows)
    if not between.empty:
        between["q_across_source_modes"] = between.groupby(
            ["stratum", "context_layer", "outcome"], group_keys=False
        )["p_value"].transform(_bh)
        between["difference_supported"] = between[
            "q_across_source_modes"
        ].le(0.05)

    manifest = {
        "contract": "chapter1_lineage_representation_bridge_primary_v1",
        "scope": scope,
        "functional_position": "(-large_bee_like + generalized_accessible) / 2",
        "source_position_uses_gift_mainland_memberships_only": True,
        "source_assignment_is_frozen_and_outcome_blind": True,
        "availability_matching": (
            "exact source-region count + floor(log2 summed source region-species supply)"
        ),
        "claim_boundary": (
            "Assembly discriminator only; no historical pollinator, colonisation, "
            "establishment, persistence, in-situ evolution, or causal geographic identification."
        ),
    }
    return {
        "responses": responses,
        "slopes": slopes,
        "between": between,
        "match_audit": match_audit,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scope",
        choices=["broad", "exact_si_all", "exact_si_direct"],
        required=True,
    )
    for name in [
        "all-scores",
        "direct-scores",
        "all-trait-ledger",
        "direct-trait-ledger",
        "status-flora",
        "gift-flora",
        "source-assignments",
        "covariates",
        "realm-assignment",
    ]:
        parser.add_argument(f"--{name}-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outputs = run_scope(
        args.scope,
        pd.read_csv(args.all_scores_csv),
        pd.read_csv(args.direct_scores_csv),
        pd.read_csv(args.all_trait_ledger_csv),
        pd.read_csv(args.direct_trait_ledger_csv),
        pd.read_csv(args.status_flora_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in {
        "responses": "lineage_representation_island_responses.csv.gz",
        "slopes": "lineage_representation_context_slopes.csv",
        "between": "lineage_representation_between_context.csv",
        "match_audit": "lineage_representation_gift_match_audit.csv",
    }.items():
        frame = outputs[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(args.output_dir / filename, index=False)
    manifest = outputs["manifest"]
    assert isinstance(manifest, dict)
    (args.output_dir / "lineage_representation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
