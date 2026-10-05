"""Species-level assembly sorting and within-lineage identifiability audit.

The current Chapter 1 trait state is defined once per accepted species. This route
therefore asks the strongest question the existing data can answer without pretending
that species-level states are population measurements:

1. among source-available, trait-resolved species, are particular trait states
   differentially represented on islands as isolation increases?
2. does the active input contain the population-specific trait dimension required to
   estimate within-species change?

For each island/outcome/source-mode row the exact identity is

raw_h1_mean
= source_species_expectation
+ species_sorting_enrichment

where the expectation preserves the observed number of island species in each
source-prevalence class. A within-species term is deliberately NOT inserted unless
population/locality-specific trait values exist.
"""
from __future__ import annotations

import json
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
    fit_joint_source_vector,
    match_gift_species,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

COMPONENTS = (
    "raw_h1_mean",
    "source_species_expectation",
    "species_sorting_enrichment",
)


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if (
        not isinstance(cfg, dict)
        or cfg.get("contract") != "chapter1_species_sorting_identifiability_v1"
    ):
        raise typer.BadParameter("unexpected species-sorting identifiability contract")
    return cfg


def _truthy_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def audit_identifiability(
    *,
    state_audit_all: pd.DataFrame,
    state_audit_direct: pd.DataFrame,
    trait_ledger_all: pd.DataFrame | None,
    trait_ledger_direct: pd.DataFrame | None,
    status_flora: pd.DataFrame,
    matched_gift: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame]:
    population_columns = {
        "population_id",
        "study_location",
        "locality",
        "latitude",
        "longitude",
        "island_id",
        "source_entity_id",
        "population_role",
        "geographic_scope",
    }

    def one(label: str, audit: pd.DataFrame, ledger: pd.DataFrame | None) -> dict[str, Any]:
        active = audit.loc[_truthy_series(audit["resolved_for_primary"])].copy()
        columns = set(active.columns)
        ledger_columns = set() if ledger is None else set(ledger.columns)
        present = sorted((columns | ledger_columns) & population_columns)
        duplicate_states = int(
            active.duplicated(["accepted_species", "trait_name"], keep=False).sum()
        )
        return {
            "evidence_scope": label,
            "n_state_rows": int(len(active)),
            "n_state_species": int(active["accepted_species"].nunique()),
            "duplicate_species_trait_rows": duplicate_states,
            "population_or_location_columns_present": present,
            "one_state_per_species_trait": bool(duplicate_states == 0),
            "within_species_change_estimable": bool(
                present and not duplicate_states
            ),
        }

    scopes = [
        one("all_analysis_eligible", state_audit_all, trait_ledger_all),
        one("direct_only", state_audit_direct, trait_ledger_direct),
    ]
    # The active H1 state tables have no population axis. Even if some upstream
    # source metadata happen to mention a locality, the analysis state is still
    # collapsed to accepted_species × trait_name.
    for row in scopes:
        row["within_species_change_estimable"] = False
        row["reason"] = (
            "active H1 trait states are accepted_species-level and contain no "
            "population-specific island/source trait state"
        )

    flora = status_flora.loc[
        status_flora["floristic_status"].astype(str).eq(
            str(config["scope"]["floristic_status"])
        ),
        ["island_id", "accepted_species"],
    ].drop_duplicates()
    island_counts = (
        flora.groupby("accepted_species", as_index=False)
        .agg(n_native_nonendemic_islands=("island_id", "nunique"))
    )
    gift_counts = (
        matched_gift.groupby("accepted_species", as_index=False)
        .agg(n_mainland_source_entities=("entity_ID", "nunique"))
    )
    all_active = state_audit_all.loc[
        _truthy_series(state_audit_all["resolved_for_primary"])
    ]
    direct_active = state_audit_direct.loc[
        _truthy_series(state_audit_direct["resolved_for_primary"])
    ]
    all_counts = (
        all_active.groupby("accepted_species", as_index=False)
        .agg(n_scored_H1_traits_all=("trait_name", "nunique"))
    )
    direct_counts = (
        direct_active.groupby("accepted_species", as_index=False)
        .agg(n_scored_H1_traits_direct=("trait_name", "nunique"))
    )
    candidates = (
        island_counts.merge(gift_counts, on="accepted_species", how="inner")
        .merge(all_counts, on="accepted_species", how="inner")
        .merge(direct_counts, on="accepted_species", how="left")
        .fillna({"n_scored_H1_traits_direct": 0})
    )
    candidates = candidates.loc[
        candidates["n_native_nonendemic_islands"].ge(2)
        & candidates["n_mainland_source_entities"].ge(1)
    ].copy()
    candidates["population_pair_needed"] = True
    candidates = candidates.sort_values(
        [
            "n_native_nonendemic_islands",
            "n_mainland_source_entities",
            "n_scored_H1_traits_direct",
        ],
        ascending=False,
    ).reset_index(drop=True)

    audit = {
        "trait_unit": str(config["identifiability"]["current_trait_unit"]),
        "required_for_within_species_change": list(
            config["identifiability"]["required_for_within_species_change"]
        ),
        "scope_audits": scopes,
        "within_species_change_estimable": False,
        "within_species_change_zero": False,
        "n_candidate_repeated_lineages_for_future_population_sampling": int(
            len(candidates)
        ),
        "interpretation": (
            "The current data can identify species sorting relative to a source pool, "
            "but cannot estimate a within-species evolutionary term. The latter is "
            "missing, not zero."
        ),
    }
    return audit, candidates


def _build_assignment_matrices(
    *,
    assignments: pd.DataFrame,
    islands: list[str],
    entities: list[int],
) -> dict[str, sparse.csr_matrix]:
    island_index = {island: i for i, island in enumerate(islands)}
    entity_index = {entity: i for i, entity in enumerate(entities)}
    work = assignments.loc[
        assignments["source_mode"].astype(str).isin(SOURCE_MODES)
        & assignments["island_id"].astype(str).isin(island_index),
        ["island_id", "source_mode", "entity_ID"],
    ].drop_duplicates().copy()
    work["island_id"] = work["island_id"].astype(str)
    work["entity_ID"] = pd.to_numeric(work["entity_ID"], errors="coerce")
    work = work.dropna(subset=["entity_ID"])
    work["entity_ID"] = work["entity_ID"].astype(int)
    work = work.loc[work["entity_ID"].isin(entity_index)]
    result: dict[str, sparse.csr_matrix] = {}
    for mode in SOURCE_MODES:
        part = work.loc[work["source_mode"].astype(str).eq(mode)]
        rows = part["island_id"].map(island_index).to_numpy(int)
        cols = part["entity_ID"].map(entity_index).to_numpy(int)
        result[mode] = sparse.csr_matrix(
            (np.ones(len(part), dtype=np.int8), (rows, cols)),
            shape=(len(islands), len(entities)),
        )
    return result


def _outcome_matrices(
    *,
    states: pd.DataFrame,
    outcome: str,
    matched_gift: pd.DataFrame,
    status_flora: pd.DataFrame,
    islands: list[str],
    floristic_status: str,
    assignments: pd.DataFrame,
) -> tuple[
    np.ndarray,
    sparse.csr_matrix,
    sparse.csr_matrix,
    dict[str, sparse.csr_matrix],
    list[str],
]:
    focal = states.loc[
        states["outcome"].astype(str).eq(outcome),
        ["accepted_species", "state"],
    ].drop_duplicates("accepted_species")
    source = matched_gift.merge(
        focal[["accepted_species"]],
        on="accepted_species",
        how="inner",
    ).drop_duplicates(["entity_ID", "accepted_species"])
    if source.empty:
        return (
            np.array([], dtype=float),
            sparse.csr_matrix((0, 0)),
            sparse.csr_matrix((0, 0)),
            {},
            [],
        )

    species = sorted(source["accepted_species"].astype(str).unique())
    species_index = {sp: i for i, sp in enumerate(species)}
    state_map = focal.set_index("accepted_species")["state"]
    states_array = np.array([float(state_map.loc[sp]) for sp in species], dtype=float)

    assigned_entities = pd.to_numeric(
        assignments.loc[
            assignments["island_id"].astype(str).isin(islands), "entity_ID"
        ],
        errors="coerce",
    ).dropna().astype(int)
    entities = sorted(set(source["entity_ID"].astype(int)) | set(assigned_entities))
    entity_index = {entity: i for i, entity in enumerate(entities)}

    src = source.loc[source["entity_ID"].astype(int).isin(entity_index)].copy()
    srows = src["entity_ID"].astype(int).map(entity_index).to_numpy(int)
    scols = src["accepted_species"].astype(str).map(species_index).to_numpy(int)
    source_presence = sparse.csr_matrix(
        (np.ones(len(src), dtype=np.int8), (srows, scols)),
        shape=(len(entities), len(species)),
    )

    island_index = {island: i for i, island in enumerate(islands)}
    flora = status_flora.loc[
        status_flora["floristic_status"].astype(str).eq(floristic_status)
        & status_flora["island_id"].astype(str).isin(island_index)
        & status_flora["accepted_species"].astype(str).isin(species_index),
        ["island_id", "accepted_species"],
    ].drop_duplicates()
    orows = flora["island_id"].astype(str).map(island_index).to_numpy(int)
    ocols = flora["accepted_species"].astype(str).map(species_index).to_numpy(int)
    observed = sparse.csr_matrix(
        (np.ones(len(flora), dtype=np.int8), (orows, ocols)),
        shape=(len(islands), len(species)),
    )

    assignment_matrices = _build_assignment_matrices(
        assignments=assignments,
        islands=islands,
        entities=entities,
    )
    return states_array, source_presence, observed, assignment_matrices, species


def _decompose_row(
    *,
    prevalence_row: sparse.csr_matrix,
    observed_row: sparse.csr_matrix,
    states: np.ndarray,
    n_species: int,
) -> dict[str, float | int] | None:
    if prevalence_row.nnz == 0 or observed_row.nnz == 0:
        return None
    candidate_idx = prevalence_row.indices
    prevalence = prevalence_row.data.astype(int)
    obs_idx = observed_row.indices

    buffer = np.zeros(n_species, dtype=np.int16)
    buffer[candidate_idx] = prevalence
    obs_prev = buffer[obs_idx]
    keep = obs_prev > 0
    if not np.any(keep):
        return None
    obs_source_idx = obs_idx[keep]
    obs_prev = obs_prev[keep]

    raw = float(np.mean(states[obs_source_idx]))
    max_prev = int(prevalence.max())
    source_count = np.bincount(prevalence, minlength=max_prev + 1)
    source_positive = np.bincount(
        prevalence,
        weights=states[candidate_idx],
        minlength=max_prev + 1,
    )
    class_mean = np.divide(
        source_positive,
        source_count,
        out=np.zeros_like(source_positive, dtype=float),
        where=source_count > 0,
    )
    expected = float(np.mean(class_mean[obs_prev]))
    sorting = raw - expected
    total_observed = int(observed_row.nnz)
    return {
        "n_source_candidate_species": int(len(candidate_idx)),
        "n_observed_trait_species": total_observed,
        "n_observed_source_candidate_species": int(len(obs_source_idx)),
        "source_overlap_fraction": float(len(obs_source_idx) / total_observed),
        "raw_h1_mean": raw,
        "source_species_expectation": expected,
        "species_sorting_enrichment": sorting,
        "identity_error": float(raw - expected - sorting),
    }


def build_species_sorting_scores(
    *,
    states: pd.DataFrame,
    matched_gift: pd.DataFrame,
    status_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    context: str,
    floristic_status: str,
    evidence_scope: str,
    outcomes: list[str],
) -> pd.DataFrame:
    islands = sorted(
        covariates.loc[
            covariates["analysis_regime"].astype(str).eq(context), "island_id"
        ].astype(str).unique()
    )
    if not islands:
        return pd.DataFrame()
    parts: list[pd.DataFrame] = []
    for outcome in outcomes:
        (
            state_array,
            source_presence,
            observed,
            assignment_matrices,
            species,
        ) = _outcome_matrices(
            states=states,
            outcome=outcome,
            matched_gift=matched_gift,
            status_flora=status_flora,
            islands=islands,
            floristic_status=floristic_status,
            assignments=assignments,
        )
        if not species:
            continue
        for mode in SOURCE_MODES:
            assignment = assignment_matrices.get(mode)
            if assignment is None:
                continue
            prevalence = (assignment @ source_presence).tocsr()
            rows: list[dict[str, Any]] = []
            for i, island in enumerate(islands):
                result = _decompose_row(
                    prevalence_row=prevalence.getrow(i),
                    observed_row=observed.getrow(i),
                    states=state_array,
                    n_species=len(species),
                )
                if result is None:
                    continue
                rows.append(
                    {
                        "island_id": island,
                        "context": context,
                        "evidence_scope": evidence_scope,
                        "outcome": outcome,
                        "source_mode": mode,
                        **result,
                    }
                )
            if rows:
                parts.append(pd.DataFrame(rows))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fit_components(
    *,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    context: str,
    evidence_scope: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    atomic_parts: list[pd.DataFrame] = []
    omnibus_rows: list[dict[str, Any]] = []
    minimum = int(config["model"]["minimum_islands_per_outcome"])
    for mode in SOURCE_MODES:
        for component in COMPONENTS:
            atomic, omnibus = fit_joint_source_vector(
                scores,
                covariates,
                source_mode=mode,
                response=component,
                minimum_islands_per_outcome=minimum,
            )
            if not atomic.empty:
                atomic.insert(0, "context", context)
                atomic.insert(0, "evidence_scope", evidence_scope)
                atomic_parts.append(atomic)
            omnibus.update(
                {
                    "context": context,
                    "evidence_scope": evidence_scope,
                    "component": component,
                }
            )
            omnibus_rows.append(omnibus)
    atomic = (
        pd.concat(atomic_parts, ignore_index=True)
        if atomic_parts
        else pd.DataFrame()
    )
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["q_across_source_modes"] = np.nan
    for (_, _, component), index in omnibus.groupby(
        ["evidence_scope", "context", "component"]
    ).groups.items():
        omnibus.loc[index, "q_across_source_modes"] = _bh(
            omnibus.loc[index, "p_value"]
        )
    return atomic, omnibus


def summarize_routes(omnibus: pd.DataFrame, alpha: float) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (scope, context, component), part in omnibus.groupby(
        ["evidence_scope", "context", "component"], sort=False
    ):
        fit = part.loc[part["status"].astype(str).eq("fit")].copy()
        rows.append(
            {
                "evidence_scope": scope,
                "context": context,
                "component": component,
                "n_source_modes_fit": int(len(fit)),
                "all_four_modes_testable": bool(len(fit) == 4),
                "all_four_modes_fdr_supported": bool(
                    len(fit) == 4
                    and fit["q_across_source_modes"].le(alpha).all()
                ),
                "all_four_modes_positive_mean": bool(
                    len(fit) == 4 and fit["mean_slope"].gt(0).all()
                ),
                "all_four_modes_negative_mean": bool(
                    len(fit) == 4 and fit["mean_slope"].lt(0).all()
                ),
                "min_mean_slope": (
                    float(fit["mean_slope"].min()) if len(fit) else np.nan
                ),
                "max_mean_slope": (
                    float(fit["mean_slope"].max()) if len(fit) else np.nan
                ),
                "max_q": (
                    float(fit["q_across_source_modes"].max())
                    if len(fit)
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


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
    trait_ledger_all_csv: Path | None = typer.Option(None),
    trait_ledger_direct_csv: Path | None = typer.Option(None),
) -> None:
    config = load_config(analysis_config_path)
    probability_config = yaml.safe_load(
        probability_config_path.read_text(encoding="utf-8")
    )
    all_audit = pd.read_csv(state_audit_all_csv)
    direct_audit = pd.read_csv(state_audit_direct_csv)
    status_flora = pd.read_csv(status_flora_csv)
    gift_flora = pd.read_csv(gift_flora_csv)
    assignments = pd.read_csv(source_assignments_csv)
    covariates = pd.read_csv(corrected_covariates_csv)
    trait_all = (
        pd.read_csv(trait_ledger_all_csv)
        if trait_ledger_all_csv is not None and trait_ledger_all_csv.exists()
        else None
    )
    trait_direct = (
        pd.read_csv(trait_ledger_direct_csv)
        if trait_ledger_direct_csv is not None and trait_ledger_direct_csv.exists()
        else None
    )

    all_states = build_outcome_states(all_audit, probability_config)
    matched_gift = match_gift_species(
        gift_flora,
        all_states["accepted_species"],
    )

    identifiability, candidates = audit_identifiability(
        state_audit_all=all_audit,
        state_audit_direct=direct_audit,
        trait_ledger_all=trait_all,
        trait_ledger_direct=trait_direct,
        status_flora=status_flora,
        matched_gift=matched_gift,
        config=config,
    )

    all_scores: list[pd.DataFrame] = []
    all_atomic: list[pd.DataFrame] = []
    all_omnibus: list[pd.DataFrame] = []
    outcomes = [str(x) for x in config["atomic_outcomes"]]
    floristic_status = str(config["scope"]["floristic_status"])
    for scope, audit in (
        ("all_analysis_eligible", all_audit),
        ("direct_only", direct_audit),
    ):
        states = build_outcome_states(audit, probability_config)
        scope_matched = matched_gift.loc[
            matched_gift["accepted_species"].astype(str).isin(
                set(states["accepted_species"].astype(str))
            )
        ].copy()
        for context in [str(x) for x in config["scope"]["contexts"]]:
            scores = build_species_sorting_scores(
                states=states,
                matched_gift=scope_matched,
                status_flora=status_flora,
                assignments=assignments,
                covariates=covariates,
                context=context,
                floristic_status=floristic_status,
                evidence_scope=scope,
                outcomes=outcomes,
            )
            if scores.empty:
                continue
            max_error = float(scores["identity_error"].abs().max())
            tolerance = float(config["additive_decomposition"]["identity_tolerance"])
            if max_error > tolerance:
                raise RuntimeError(
                    f"species sorting identity failed for {scope}/{context}: "
                    f"{max_error} > {tolerance}"
                )
            all_scores.append(scores)
            atomic, omnibus = fit_components(
                scores=scores,
                covariates=covariates.loc[
                    covariates["analysis_regime"].astype(str).eq(context)
                ].copy(),
                context=context,
                evidence_scope=scope,
                config=config,
            )
            if not atomic.empty:
                all_atomic.append(atomic)
            all_omnibus.append(omnibus)

    scores = (
        pd.concat(all_scores, ignore_index=True)
        if all_scores
        else pd.DataFrame()
    )
    atomic = (
        pd.concat(all_atomic, ignore_index=True)
        if all_atomic
        else pd.DataFrame()
    )
    omnibus = (
        pd.concat(all_omnibus, ignore_index=True)
        if all_omnibus
        else pd.DataFrame()
    )
    routes = summarize_routes(
        omnibus,
        alpha=float(config["multiplicity"]["alpha"]),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(output_dir / "species_sorting_scores.csv", index=False)
    atomic.to_csv(output_dir / "atomic_slopes.csv", index=False)
    omnibus.to_csv(output_dir / "omnibus.csv", index=False)
    routes.to_csv(output_dir / "route_summary.csv", index=False)
    candidates.to_csv(output_dir / "candidate_repeated_lineages.csv", index=False)
    (output_dir / "identifiability_audit.json").write_text(
        json.dumps(identifiability, indent=2) + "\n",
        encoding="utf-8",
    )

    result = {
        "contract": config["contract"],
        "status": "complete",
        "n_score_rows": int(len(scores)),
        "n_score_islands": (
            int(scores["island_id"].nunique()) if not scores.empty else 0
        ),
        "max_identity_error": (
            float(scores["identity_error"].abs().max()) if not scores.empty else np.nan
        ),
        "identifiability": identifiability,
        "claim_boundary": (
            "Species sorting is identifiable on source-evaluable support. "
            "Within-species change is not estimable from the current species-level "
            "trait state and must be treated as missing rather than zero."
        ),
    }
    (output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(routes.to_csv(index=False))


if __name__ == "__main__":
    app()
