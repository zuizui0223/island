"""Literal additive source-lineage decomposition of the current seven-response H1.

V2 repairs one alignment issue discovered in the first current-H1 source-lineage audit:
the historical lineage bridge weighted genera by all island species, whereas current H1
is defined on trait-resolved species separately for each atomic response.

Here each island/outcome/source-mode row uses only trait-resolved, source-candidate,
source-genus-scored native non-endemic species. The decomposition is exact:

raw_h1_mean
= source_expectation
+ genus_entry_enrichment
+ within_genus_loading
+ within_genus_trait_residual

Both row-level and regression-slope closure are mandatory gates.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy import sparse

from island_v2.chapter1_all_data_probability import _bh
from island_v2.chapter1_current_h1_source_lineage_entry import (
    SOURCE_MODES,
    build_outcome_states,
    build_raw_source_availability,
    build_source_genus_positions,
    fit_joint_source_vector,
    match_gift_species,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

COMPONENTS = (
    "raw_h1_mean",
    "source_expectation",
    "genus_entry_enrichment",
    "within_genus_loading",
    "within_genus_trait_residual",
    "species_enrichment",
)


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if (
        not isinstance(cfg, dict)
        or cfg.get("contract")
        != "chapter1_current_h1_source_lineage_decomposition_v2"
    ):
        raise typer.BadParameter("unexpected current H1 source-lineage v2 contract")
    return cfg


def _richness_bins(values: np.ndarray) -> np.ndarray:
    bins = np.zeros_like(values, dtype=np.int8)
    bins[(values > 1) & (values <= 2)] = 1
    bins[(values > 2) & (values <= 4)] = 2
    bins[(values > 4) & (values <= 8)] = 3
    bins[values > 8] = 4
    return bins


def _build_assignment_and_availability(
    *,
    source_availability: pd.DataFrame,
    assignments: pd.DataFrame,
    islands: list[str],
    genera: list[str],
) -> tuple[
    dict[str, tuple[np.ndarray, np.ndarray]],
    dict[str, int],
]:
    island_index = {island: i for i, island in enumerate(islands)}
    genus_index = {genus: i for i, genus in enumerate(genera)}

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
    entity_index = {entity: i for i, entity in enumerate(entities)}

    arow = avail["entity_ID"].map(entity_index).to_numpy(int)
    acol = avail["genus"].astype(str).map(genus_index).to_numpy(int)
    presence = sparse.csr_matrix(
        (np.ones(len(avail), dtype=float), (arow, acol)),
        shape=(len(entity_index), len(genus_index)),
    )
    richness = sparse.csr_matrix(
        (
            pd.to_numeric(avail["source_species_richness"], errors="coerce")
            .fillna(0)
            .to_numpy(float),
            (arow, acol),
        ),
        shape=presence.shape,
    )

    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for source_mode in SOURCE_MODES:
        part = assign.loc[assign["source_mode"].astype(str).eq(source_mode)]
        rows = part["island_id"].map(island_index).to_numpy(int)
        cols = part["entity_ID"].map(entity_index).to_numpy(int)
        assignment_matrix = sparse.csr_matrix(
            (np.ones(len(part), dtype=float), (rows, cols)),
            shape=(len(island_index), len(entity_index)),
        )
        out[source_mode] = (
            (assignment_matrix @ presence).toarray().astype(np.int16),
            (assignment_matrix @ richness).toarray().astype(np.float32),
        )
    return out, genus_index


def _build_island_state_matrices(
    *,
    states: pd.DataFrame,
    status_flora: pd.DataFrame,
    outcome: str,
    islands: list[str],
    genera: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    island_index = {island: i for i, island in enumerate(islands)}
    genus_index = {genus: i for i, genus in enumerate(genera)}
    focal_states = states.loc[
        states["outcome"].astype(str).eq(outcome),
        ["accepted_species", "state", "genus"],
    ].copy()
    flora = status_flora.loc[
        status_flora["floristic_status"].astype(str).eq("native_nonendemic"),
        ["island_id", "accepted_species"],
    ].drop_duplicates().copy()
    flora["island_id"] = flora["island_id"].astype(str)
    joined = flora.merge(
        focal_states,
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    joined = joined.loc[
        joined["island_id"].isin(island_index)
        & joined["genus"].astype(str).isin(genus_index)
    ].copy()
    grouped = joined.groupby(["island_id", "genus"], as_index=False).agg(
        n_scored_species=("accepted_species", "nunique"),
        state_sum=("state", "sum"),
    )
    rows = grouped["island_id"].map(island_index).to_numpy(int)
    cols = grouped["genus"].astype(str).map(genus_index).to_numpy(int)
    counts = sparse.csr_matrix(
        (
            grouped["n_scored_species"].to_numpy(float),
            (rows, cols),
        ),
        shape=(len(island_index), len(genus_index)),
    ).toarray()
    state_sum = sparse.csr_matrix(
        (
            grouped["state_sum"].to_numpy(float),
            (rows, cols),
        ),
        shape=(len(island_index), len(genus_index)),
    ).toarray()
    return counts, state_sum


def _decompose_one(
    *,
    prevalence: np.ndarray,
    source_richness: np.ndarray,
    counts: np.ndarray,
    state_sum: np.ndarray,
    positions: np.ndarray,
    minimum_represented_genera: int,
) -> dict[str, float | int] | None:
    candidate = prevalence > 0
    c = counts.astype(float).copy()
    s = state_sum.astype(float).copy()
    c[~candidate] = 0.0
    s[~candidate] = 0.0
    represented = c > 0
    n_represented = int(represented.sum())
    n_species = int(c.sum())
    if n_represented < int(minimum_represented_genera) or n_species <= 0:
        return None

    raw_h1_mean = float(s.sum() / n_species)
    observed_entry = float(np.mean(positions[represented]))
    observed_species_genus_position = float(np.sum(positions * c) / n_species)

    richness_bins = _richness_bins(source_richness)
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
        represented_species = int(np.sum(c[class_mask]))
        if represented_genera == 0 and represented_species == 0:
            continue
        class_mean = float(np.mean(positions[class_mask]))
        expected_entry_sum += represented_genera * class_mean
        expected_species_sum += represented_species * class_mean

    expected_entry = expected_entry_sum / n_represented
    source_expectation = expected_species_sum / n_species
    genus_entry_enrichment = observed_entry - expected_entry
    species_enrichment = observed_species_genus_position - source_expectation
    within_genus_loading = species_enrichment - genus_entry_enrichment
    within_genus_trait_residual = raw_h1_mean - observed_species_genus_position

    reconstructed = (
        source_expectation
        + genus_entry_enrichment
        + within_genus_loading
        + within_genus_trait_residual
    )
    identity_error = raw_h1_mean - reconstructed
    return {
        "n_represented_genera": n_represented,
        "n_scored_species": n_species,
        "n_candidate_genera": int(candidate.sum()),
        "raw_h1_mean": raw_h1_mean,
        "source_expectation": source_expectation,
        "observed_entry_mean": observed_entry,
        "expected_entry_mean": expected_entry,
        "observed_species_genus_position": observed_species_genus_position,
        "genus_entry_enrichment": genus_entry_enrichment,
        "species_enrichment": species_enrichment,
        "within_genus_loading": within_genus_loading,
        "within_genus_trait_residual": within_genus_trait_residual,
        "identity_error": identity_error,
    }


def build_decomposition_scores(
    *,
    states: pd.DataFrame,
    positions: pd.DataFrame,
    source_availability: pd.DataFrame,
    status_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    evidence_scope: str,
) -> pd.DataFrame:
    tropical_islands = sorted(
        covariates.loc[
            covariates["analysis_regime"].astype(str).eq(
                str(config["scope"]["context"])
            ),
            "island_id",
        ]
        .astype(str)
        .unique()
    )
    minimum = int(config["source_pool"]["minimum_represented_genera"])
    rows: list[dict[str, Any]] = []

    for outcome in [str(x) for x in config["atomic_outcomes"]]:
        pos = positions.loc[
            positions["outcome"].astype(str).eq(outcome)
        ].drop_duplicates("genus").copy()
        if pos.empty:
            continue
        genera = sorted(pos["genus"].astype(str).unique())
        pos_map = pos.set_index("genus")["source_h1_position"]
        position_array = np.array([float(pos_map.loc[g]) for g in genera], dtype=float)

        mode_matrices, _ = _build_assignment_and_availability(
            source_availability=source_availability,
            assignments=assignments,
            islands=tropical_islands,
            genera=genera,
        )
        counts, state_sum = _build_island_state_matrices(
            states=states,
            status_flora=status_flora,
            outcome=outcome,
            islands=tropical_islands,
            genera=genera,
        )

        for source_mode in SOURCE_MODES:
            prevalence, richness = mode_matrices[source_mode]
            for i, island_id in enumerate(tropical_islands):
                result = _decompose_one(
                    prevalence=prevalence[i],
                    source_richness=richness[i],
                    counts=counts[i],
                    state_sum=state_sum[i],
                    positions=position_array,
                    minimum_represented_genera=minimum,
                )
                if result is None:
                    continue
                rows.append(
                    {
                        "island_id": island_id,
                        "evidence_scope": evidence_scope,
                        "outcome": outcome,
                        "source_mode": source_mode,
                        **result,
                    }
                )
    return pd.DataFrame(rows)


def _component_models(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    atomic_parts: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    min_islands = int(config["model"]["minimum_islands_per_outcome"])
    for source_mode in SOURCE_MODES:
        for component in COMPONENTS:
            atomic, omnibus = fit_joint_source_vector(
                scores,
                covariates,
                source_mode=source_mode,
                response=component,
                minimum_islands_per_outcome=min_islands,
            )
            if not atomic.empty:
                atomic.insert(0, "evidence_scope", evidence_scope)
                atomic_parts.append(atomic)
            omnibus["evidence_scope"] = evidence_scope
            omnibus["component"] = component
            omnibus_rows.append(omnibus)
    atomic_all = (
        pd.concat(atomic_parts, ignore_index=True)
        if atomic_parts
        else pd.DataFrame()
    )
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["q_across_source_modes"] = np.nan
    for component, index in omnibus.groupby("component").groups.items():
        omnibus.loc[index, "q_across_source_modes"] = _bh(
            omnibus.loc[index, "p_value"]
        )
    return atomic_all, omnibus


def _slope_identity_gate(
    atomic: pd.DataFrame,
    *,
    tolerance: float,
) -> dict[str, Any]:
    if atomic.empty:
        return {"status": "fail", "reason": "no_atomic_results"}
    index = atomic.pivot_table(
        index=["evidence_scope", "source_mode", "outcome"],
        columns="response",
        values="estimate",
        aggfunc="first",
    )
    required = {
        "raw_h1_mean",
        "source_expectation",
        "genus_entry_enrichment",
        "within_genus_loading",
        "within_genus_trait_residual",
    }
    if missing := required - set(index.columns):
        return {"status": "fail", "missing_components": sorted(missing)}
    reconstructed = (
        index["source_expectation"]
        + index["genus_entry_enrichment"]
        + index["within_genus_loading"]
        + index["within_genus_trait_residual"]
    )
    error = (index["raw_h1_mean"] - reconstructed).abs()
    maximum = float(error.max()) if len(error) else float("nan")
    return {
        "status": "pass" if math.isfinite(maximum) and maximum <= tolerance else "fail",
        "max_absolute_slope_identity_error": maximum,
        "tolerance": tolerance,
        "n_atomic_identities": int(len(error)),
    }


def _summarize(
    scores: pd.DataFrame,
    omnibus: pd.DataFrame,
    slope_gate: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    tolerance = float(config["additive_decomposition"]["identity_tolerance"])
    row_error = pd.to_numeric(scores["identity_error"], errors="coerce").abs()
    row_max = float(row_error.max()) if len(row_error) else float("nan")
    alpha = float(config["multiplicity"]["alpha"])
    components: dict[str, Any] = {}
    for component in COMPONENTS:
        part = omnibus.loc[omnibus["component"].eq(component)].copy()
        components[component] = {
            "all_four_modes_testable": bool(
                len(part) == 4 and part["status"].eq("fit").all()
            ),
            "all_four_modes_fdr_supported": bool(
                len(part) == 4 and part["q_across_source_modes"].le(alpha).all()
            ),
            "all_four_modes_positive_mean": bool(
                len(part) == 4 and part["mean_slope"].gt(0).all()
            ),
            "all_four_modes_negative_mean": bool(
                len(part) == 4 and part["mean_slope"].lt(0).all()
            ),
            "min_mean_slope": (
                float(part["mean_slope"].min()) if len(part) else np.nan
            ),
            "max_mean_slope": (
                float(part["mean_slope"].max()) if len(part) else np.nan
            ),
            "max_q": (
                float(part["q_across_source_modes"].max()) if len(part) else np.nan
            ),
            "min_positive_atomic_slopes": (
                int(part["n_positive_slopes"].min()) if len(part) else 0
            ),
            "max_positive_atomic_slopes": (
                int(part["n_positive_slopes"].max()) if len(part) else 0
            ),
        }
    gates_pass = bool(
        math.isfinite(row_max)
        and row_max <= tolerance
        and slope_gate.get("status") == "pass"
    )
    return {
        "row_identity_gate": {
            "status": "pass" if row_max <= tolerance else "fail",
            "max_absolute_identity_error": row_max,
            "tolerance": tolerance,
            "n_rows": int(len(scores)),
        },
        "slope_identity_gate": slope_gate,
        "all_identity_gates_pass": gates_pass,
        "components": components,
    }


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
    config: dict[str, Any],
) -> dict[str, Any]:
    states = build_outcome_states(state_audit, probability_config)
    positions = build_source_genus_positions(
        matched_gift,
        states,
        minimum_scored_species=int(
            config["source_pool"]["minimum_scored_source_species_per_genus"]
        ),
    )
    scores = build_decomposition_scores(
        states=states,
        positions=positions,
        source_availability=source_availability,
        status_flora=status_flora,
        assignments=assignments,
        covariates=covariates,
        config=config,
        evidence_scope=evidence_scope,
    )
    atomic, omnibus = _component_models(
        scores,
        covariates,
        config,
        evidence_scope=evidence_scope,
    )
    slope_gate = _slope_identity_gate(
        atomic,
        tolerance=float(config["additive_decomposition"]["identity_tolerance"]),
    )
    summary = _summarize(scores, omnibus, slope_gate, config)
    return {
        "states": states,
        "positions": positions,
        "scores": scores,
        "atomic": atomic,
        "omnibus": omnibus,
        "summary": summary,
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
    config = load_config(analysis_config_path)
    all_audit = pd.read_csv(state_audit_all_csv)
    direct_audit = pd.read_csv(state_audit_direct_csv)
    status_flora = pd.read_csv(status_flora_csv)
    gift_flora = pd.read_csv(gift_flora_csv)
    assignments = pd.read_csv(source_assignments_csv)
    covariates = pd.read_csv(corrected_covariates_csv)

    all_states = build_outcome_states(all_audit, probability_config)
    matched_gift = match_gift_species(gift_flora, all_states["accepted_species"])
    source_availability = build_raw_source_availability(gift_flora)

    scopes: dict[str, Any] = {}
    for label, audit in [
        ("all_analysis_eligible", all_audit),
        ("direct_only", direct_audit),
    ]:
        scopes[label] = run_scope(
            evidence_scope=label,
            state_audit=audit,
            status_flora=status_flora,
            gift_flora=gift_flora,
            matched_gift=matched_gift,
            source_availability=source_availability,
            assignments=assignments,
            covariates=covariates,
            probability_config=probability_config,
            config=config,
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = []
    for label, result in scopes.items():
        scope_dir = output_dir / label
        scope_dir.mkdir(parents=True, exist_ok=True)
        result["positions"].to_csv(scope_dir / "positions.csv", index=False)
        result["scores"].to_csv(scope_dir / "scores.csv", index=False)
        result["atomic"].to_csv(scope_dir / "atomic_slopes.csv", index=False)
        result["omnibus"].to_csv(scope_dir / "omnibus.csv", index=False)
        (scope_dir / "summary.json").write_text(
            json.dumps(result["summary"], indent=2) + "\n",
            encoding="utf-8",
        )
        for component, item in result["summary"]["components"].items():
            summary_rows.append(
                {
                    "evidence_scope": label,
                    "component": component,
                    **item,
                }
            )
    pd.DataFrame(summary_rows).to_csv(
        output_dir / "component_summary.csv", index=False
    )

    manifest = {
        "contract": config["contract"],
        "v1_alignment_issue_repaired": True,
        "source_assignments_reused": True,
        "current_H1_trait_recoding_reused": True,
        "scope_summaries": {
            label: scopes[label]["summary"] for label in scopes
        },
        "claim_ceiling": config["claim_ceiling"],
    }
    (output_dir / "RESULT.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(pd.DataFrame(summary_rows).to_csv(index=False))


if __name__ == "__main__":
    app()
