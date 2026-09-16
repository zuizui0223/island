"""Species-level reproductive-assurance moderation of the global GloPL distance gradient.

The preflight joins only GloPL species/geography/publication metadata to the frozen
species-direct reproductive-assurance ledger. GloPL effect columns are read only after
that overlap/support audit is written. Associations do not identify causal mediation.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_global_distance import (
    CONTEXTS,
    MEASUREMENT_COLUMNS,
    _clustered_wls,
    _context_dummies,
    _lower,
    _measurement_dummies,
    _p2,
    _truthy,
    build_preflight_table,
    load_config as load_parent_config,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_glopl_reproductive_assurance_moderation_v1"
PARENT_CONFIG_PATH = Path("config/chapter1_h5_glopl_global_distance_v1.yml")
PRIMARY_CONTEXTS = ("northern_midlatitude", "tropical", "southern_extratropical")


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected GloPL reproductive-assurance moderation contract")
    return config


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_species_name(value: object) -> str:
    text = str(value or "").replace("_", " ").strip().casefold()
    return re.sub(r"\s+", " ", text)


def build_trait_assignments(ledger: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    needed = {
        "accepted_species",
        "axis",
        "trait_name",
        "normalized_value",
        "resolution_status",
        "quality",
    }
    if missing := needed - set(ledger.columns):
        raise typer.BadParameter(f"trait ledger missing columns: {sorted(missing)}")
    work = ledger.loc[
        ledger["axis"].astype(str).eq("reproductive_assurance")
        & ledger["resolution_status"].astype(str).eq("resolved")
        & ledger["quality"].astype(str).isin(["high", "medium"])
    ].copy()
    rows: list[dict[str, Any]] = []
    for trait, spec in config["trait_recodes"].items():
        source_trait = str(spec["source_trait"])
        protected = {str(x) for x in spec["protected_states"]}
        unprotected = {str(x) for x in spec["unprotected_states"]}
        part = work.loc[work["trait_name"].astype(str).eq(source_trait)].copy()
        for row in part.itertuples(index=False):
            value = str(row.normalized_value)
            if value in protected:
                state = 1
            elif value in unprotected:
                state = 0
            else:
                continue
            species_key = normalize_species_name(row.accepted_species)
            if not species_key:
                continue
            rows.append(
                {
                    "species_key": species_key,
                    "accepted_species": str(row.accepted_species),
                    "trait": str(trait),
                    "trait_state": int(state),
                    "trait_source_value": value,
                    "quality": str(row.quality),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(
            columns=[
                "species_key",
                "accepted_species",
                "trait",
                "trait_state",
                "trait_source_value",
                "quality",
            ]
        )
    conflict = out.groupby(["species_key", "trait"])["trait_state"].nunique()
    if bool((conflict > 1).any()):
        bad = conflict.loc[conflict > 1].index.tolist()[:10]
        raise typer.BadParameter(f"conflicting reproductive-assurance assignments: {bad}")
    return out.drop_duplicates(["species_key", "trait"]).reset_index(drop=True)


def _support_counts(part: pd.DataFrame) -> dict[str, int]:
    return {
        "n_rows": int(len(part)),
        "n_species": int(part["species_key"].nunique()),
        "n_sites": int(part["site_key"].nunique()),
        "n_publications": int(part["study_key"].nunique()),
        "n_offshore_sites": int(
            part.loc[
                pd.to_numeric(part["distance_to_major_continent_km"], errors="coerce") > 1e-9,
                "site_key",
            ].nunique()
        ),
    }


def summarize_trait_support(rows: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    gate = config["preflight"]["support_gates"]
    summaries: dict[str, Any] = {}
    for trait in config["trait_recodes"]:
        part = rows.loc[rows["trait"].astype(str).eq(str(trait))].copy()
        classes = {
            str(state): _support_counts(part.loc[pd.to_numeric(part["trait_state"], errors="coerce").eq(state)])
            for state in (0, 1)
        }
        n_species = int(part["species_key"].nunique())
        global_admitted = n_species >= int(gate["global_min_matched_species_per_trait"]) and all(
            classes[str(state)]["n_publications"] >= int(gate["global_min_publications_per_class"])
            and classes[str(state)]["n_sites"] >= int(gate["global_min_sites_per_class"])
            and classes[str(state)]["n_offshore_sites"] >= int(gate["global_min_offshore_sites_per_class"])
            for state in (0, 1)
        )
        contexts: dict[str, Any] = {}
        context_admitted = True
        for context in ("northern_midlatitude", "tropical"):
            context_part = part.loc[part["analysis_regime"].astype(str).eq(context)]
            state_counts = {
                str(state): _support_counts(
                    context_part.loc[pd.to_numeric(context_part["trait_state"], errors="coerce").eq(state)]
                )
                for state in (0, 1)
            }
            admitted = all(
                state_counts[str(state)]["n_publications"]
                >= int(gate["context_min_publications_per_class"])
                and state_counts[str(state)]["n_sites"] >= int(gate["context_min_sites_per_class"])
                and state_counts[str(state)]["n_offshore_sites"]
                >= int(gate["context_min_offshore_sites_per_class"])
                for state in (0, 1)
            )
            contexts[context] = {"classes": state_counts, "admitted": bool(admitted)}
            context_admitted = context_admitted and bool(admitted)
        summaries[str(trait)] = {
            "n_matched_species": n_species,
            "classes": classes,
            "global_admitted": bool(global_admitted),
            "contexts": contexts,
            "north_tropical_context_diagnostic_admitted": bool(context_admitted),
        }
    return summaries


def aggregate_species_measurement_cells(rows: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame()
    group_cols = [
        "study_key",
        "site_key",
        "species_key",
        "analysis_regime",
        "z_distance",
        "trait",
        "trait_state",
        *MEASUREMENT_COLUMNS,
    ]
    if missing := set(group_cols + ["PL_Effect_Size"]) - set(rows.columns):
        raise typer.BadParameter(f"species measurement aggregation missing: {sorted(missing)}")
    out = (
        rows.groupby(group_cols, as_index=False, dropna=False)
        .agg(PL_Effect_Size=("PL_Effect_Size", "mean"), n_effect_rows=("PL_Effect_Size", "size"))
        .reset_index(drop=True)
    )
    n_cells = out.groupby("study_key")["site_key"].transform("size").astype(float)
    out["analysis_weight"] = float(config["analysis"]["publication_total_weight"]) / n_cells
    return out


def _build_global_moderation_design(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    distance = pd.to_numeric(frame["z_distance"], errors="coerce").to_numpy(float)
    trait_state = pd.to_numeric(frame["trait_state"], errors="coerce").to_numpy(float)
    context_cols, context_names, _ = _context_dummies(frame)
    measure_cols, measure_names = _measurement_dummies(frame)
    columns = [
        np.ones(len(frame)),
        *context_cols,
        distance,
        trait_state,
        distance * trait_state,
        *measure_cols,
    ]
    names = [
        "intercept",
        *context_names,
        "z_distance",
        "protected_trait",
        "z_distance:protected_trait",
        *measure_names,
    ]
    return np.column_stack(columns), names


def fit_global_moderation(cells: pd.DataFrame) -> dict[str, Any]:
    if cells.empty or cells["trait_state"].nunique() < 2:
        return {"evaluable": False, "reason": "both_trait_states_not_available"}
    X, names = _build_global_moderation_design(cells)
    fit = _clustered_wls(cells, X, names)
    if not fit.get("evaluable"):
        return fit
    index = {name: i for i, name in enumerate(fit["names"])}
    d = index["z_distance"]
    interaction = index["z_distance:protected_trait"]
    unprotected = float(fit["beta"][d])
    inter = float(fit["beta"][interaction])
    inter_se = float(fit["se"][interaction])
    z = inter / inter_se if inter_se > 0 else float("nan")
    return {
        "evaluable": True,
        "unprotected_distance_slope": unprotected,
        "protected_distance_slope": unprotected + inter,
        "distance_by_trait_interaction": inter,
        "interaction_se": inter_se,
        "interaction_two_sided_p": _p2(z),
        "interaction_one_sided_negative_p": _lower(z),
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
        "n_species": int(cells["species_key"].nunique()),
    }


def fit_context_diagnostic(cells: pd.DataFrame) -> dict[str, Any]:
    pair = cells.loc[cells["analysis_regime"].isin(["northern_midlatitude", "tropical"])].copy()
    if pair.empty or pair["trait_state"].nunique() < 2 or pair["analysis_regime"].nunique() < 2:
        return {"evaluable": False, "reason": "context_or_trait_state_support_failed"}
    distance = pd.to_numeric(pair["z_distance"], errors="coerce").to_numpy(float)
    trait_state = pd.to_numeric(pair["trait_state"], errors="coerce").to_numpy(float)
    tropical = pair["analysis_regime"].astype(str).eq("tropical").to_numpy(float)
    measure_cols, measure_names = _measurement_dummies(pair)
    columns = [
        np.ones(len(pair)),
        tropical,
        distance,
        trait_state,
        distance * trait_state,
        distance * tropical,
        trait_state * tropical,
        distance * trait_state * tropical,
        *measure_cols,
    ]
    names = [
        "intercept",
        "tropical",
        "z_distance",
        "protected_trait",
        "z_distance:protected_trait",
        "z_distance:tropical",
        "protected_trait:tropical",
        "z_distance:protected_trait:tropical",
        *measure_names,
    ]
    fit = _clustered_wls(pair, np.column_stack(columns), names)
    if not fit.get("evaluable"):
        return fit
    index = {name: i for i, name in enumerate(fit["names"])}
    north_idx = index["z_distance:protected_trait"]
    three_idx = index["z_distance:protected_trait:tropical"]
    north_buffer = float(fit["beta"][north_idx])
    diff = float(fit["beta"][three_idx])
    diff_se = float(fit["se"][three_idx])
    z = diff / diff_se if diff_se > 0 else float("nan")
    return {
        "evaluable": True,
        "northern_distance_by_trait_interaction": north_buffer,
        "tropical_distance_by_trait_interaction": north_buffer + diff,
        "tropical_minus_northern_buffering_interaction": diff,
        "interaction_two_sided_p": _p2(z),
        "n_cells": fit["n_cells"],
        "n_publications": fit["n_publications"],
        "n_sites": fit["n_sites"],
        "n_species": int(pair["species_key"].nunique()),
    }


def classify_trait_result(
    primary: dict[str, Any], sensitivities: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    sensitivity_direction = bool(sensitivities) and all(
        float(result.get("distance_by_trait_interaction", float("nan"))) < 0
        for result in sensitivities.values()
    )
    supported = bool(
        primary.get("evaluable")
        and float(primary.get("unprotected_distance_slope", float("nan"))) > 0
        and float(primary.get("distance_by_trait_interaction", float("nan"))) < 0
        and float(primary.get("interaction_one_sided_negative_p", 1.0)) <= 0.05
        and sensitivity_direction
    )
    return {
        "supported": supported,
        "sensitivity_direction_retained": bool(sensitivity_direction),
    }


def _read_glopl_preflight(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    columns = [str(x) for x in config["preflight"]["glopl_allowed_columns"]]
    return pd.read_csv(path, usecols=columns, dtype=str, encoding="latin-1").fillna("")


def _read_trait_ledger(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    columns = [str(x) for x in config["preflight"]["trait_allowed_columns"]]
    return pd.read_csv(path, usecols=columns, dtype=str).fillna("")


def _read_effects(path: Path) -> pd.DataFrame:
    columns = ["PL_Effect_Size", *MEASUREMENT_COLUMNS]
    frame = pd.read_csv(path, usecols=columns, dtype=str, encoding="latin-1").fillna("")
    frame.insert(0, "row_id", np.arange(len(frame), dtype=int))
    return frame


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")


def _prepare_effect_rows(
    effects: pd.DataFrame,
    overlap: pd.DataFrame,
    *,
    distance_mean: float,
    distance_sd: float,
) -> pd.DataFrame:
    work = overlap.merge(effects, on="row_id", how="left", validate="many_to_one")
    work["PL_Effect_Size"] = pd.to_numeric(work["PL_Effect_Size"], errors="coerce")
    work = work.loc[np.isfinite(work["PL_Effect_Size"].to_numpy(float))].copy()
    log_distance = pd.to_numeric(work["log1p_distance_to_major_continent_km"], errors="coerce")
    work["z_distance"] = (log_distance - float(distance_mean)) / float(distance_sd)
    for column in MEASUREMENT_COLUMNS:
        work[column] = work[column].fillna("").astype(str).str.strip()
    work["Constant_added_bool"] = _truthy(work["Constant_added"])
    return work.reset_index(drop=True)


@app.command("preflight")
def preflight(
    glopl_csv: Path = typer.Option(..., exists=True),
    trait_ledger: Path = typer.Option(..., exists=True),
    land_geojson: Path = typer.Option(..., exists=True),
    parent_result_json: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    if _sha256(glopl_csv) != str(config["sources"]["glopl"]["sha256"]):
        raise typer.BadParameter("GloPL SHA-256 mismatch")
    parent_result = json.loads(parent_result_json.read_text(encoding="utf-8"))
    if parent_result.get("contract") != "chapter1_h5_glopl_global_distance_v1":
        raise typer.BadParameter("unexpected parent GloPL distance result")
    parent_config = load_parent_config(PARENT_CONFIG_PATH)
    metadata = _read_glopl_preflight(glopl_csv, config)
    land = gpd.read_file(land_geojson)
    geography = build_preflight_table(metadata, land, parent_config)
    geography["species_key"] = geography["Species_accepted_names"].map(normalize_species_name)
    ledger = _read_trait_ledger(trait_ledger, config)
    assignments = build_trait_assignments(ledger, config)
    overlap = geography.merge(assignments, on="species_key", how="inner", validate="many_to_many")
    keep = [
        "row_id",
        "Species_accepted_names",
        "species_key",
        "trait",
        "trait_state",
        "trait_source_value",
        "quality",
        "site_key",
        "study_key",
        "analysis_regime",
        "distance_to_major_continent_km",
        "log1p_distance_to_major_continent_km",
        "valid_coordinate",
    ]
    overlap = overlap[keep].copy()
    overlap = overlap.loc[
        overlap["valid_coordinate"].astype(bool)
        & overlap["site_key"].astype(str).ne("")
        & overlap["study_key"].astype(str).ne("")
        & overlap["analysis_regime"].isin(CONTEXTS)
    ].reset_index(drop=True)
    support = summarize_trait_support(overlap, config)
    distance_scaling = parent_result.get("distance_standardization", {})
    result = {
        "contract": CONTRACT,
        "status": "completed_outcome_blind_trait_overlap_preflight",
        "glopl_sha256": _sha256(glopl_csv),
        "trait_ledger_sha256": _sha256(trait_ledger),
        "geography_sha256": _sha256(land_geojson),
        "n_glopl_rows": int(len(metadata)),
        "n_trait_assignments": int(len(assignments)),
        "n_overlap_rows": int(len(overlap)),
        "n_overlap_species": int(overlap["species_key"].nunique()),
        "trait_support": support,
        "admitted_traits": [trait for trait, item in support.items() if item["global_admitted"]],
        "distance_standardization": {
            "mean": float(distance_scaling["mean"]),
            "sd": float(distance_scaling["sd"]),
            "source": "frozen_parent_global_GloPL_result",
        },
        "forbidden_outcome_columns_materialized": False,
        "claim_ceiling": config["claim_ceiling"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    overlap.to_csv(output_dir / "TRAIT_OVERLAP_PREFLIGHT.csv.gz", index=False, compression="gzip")
    _write_json(output_dir / "PREFLIGHT.json", result)
    lines = [
        "# GloPL reproductive-assurance overlap preflight",
        "",
        f"- GloPL rows: {len(metadata)}",
        f"- overlap species: {result['n_overlap_species']}",
        f"- admitted traits: {', '.join(result['admitted_traits']) or 'none'}",
    ]
    for trait, item in support.items():
        lines.append(
            f"- {trait}: species={item['n_matched_species']}, global_gate={item['global_admitted']}, "
            f"N-T_gate={item['north_tropical_context_diagnostic_admitted']}"
        )
    (output_dir / "PREFLIGHT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2))


@app.command("analyse")
def analyse(
    glopl_csv: Path = typer.Option(..., exists=True),
    preflight_csv: Path = typer.Option(..., exists=True),
    preflight_json: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    preflight_result = json.loads(preflight_json.read_text(encoding="utf-8"))
    admitted = [str(x) for x in preflight_result.get("admitted_traits", [])]
    output_dir.mkdir(parents=True, exist_ok=True)
    if not admitted:
        result = {
            "contract": CONTRACT,
            "status": "completed_not_evaluable_no_trait_passed_preflight",
            "preflight": preflight_result,
            "trait_results": {},
            "decision": {
                "supported_trait_count": 0,
                "route_A_family_supported": False,
                "classification": "reproductive_assurance_buffering_not_evaluable",
            },
            "claim_ceiling": config["claim_ceiling"],
        }
        _write_json(output_dir / "RESULT.json", result)
        (output_dir / "RESULT_SUMMARY.md").write_text(
            "# GloPL reproductive-assurance moderation\n\n- no trait passed the frozen preflight gate\n",
            encoding="utf-8",
        )
        typer.echo(json.dumps(result, indent=2))
        return
    effects = _read_effects(glopl_csv)
    overlap = pd.read_csv(preflight_csv, dtype=str).fillna("")
    overlap["row_id"] = pd.to_numeric(overlap["row_id"], errors="raise").astype(int)
    scaling = preflight_result["distance_standardization"]
    rows = _prepare_effect_rows(
        effects,
        overlap.loc[overlap["trait"].isin(admitted)].copy(),
        distance_mean=float(scaling["mean"]),
        distance_sd=float(scaling["sd"]),
    )
    trait_results: dict[str, Any] = {}
    for trait in config["trait_recodes"]:
        support = preflight_result["trait_support"][trait]
        if trait not in admitted:
            trait_results[trait] = {"evaluable": False, "reason": "preflight_support_failed"}
            continue
        trait_rows = rows.loc[rows["trait"].eq(trait)].copy()
        primary_cells = aggregate_species_measurement_cells(trait_rows, config)
        primary = fit_global_moderation(primary_cells)
        sensitivities: dict[str, Any] = {}
        sensitivity_masks = {
            "supplemental_only": trait_rows["PL_Effect_Size_Type2"].astype(str).eq("Sup"),
            "no_zero_constant": ~trait_rows["Constant_added_bool"].astype(bool),
        }
        for name, mask in sensitivity_masks.items():
            sens_cells = aggregate_species_measurement_cells(trait_rows.loc[mask].copy(), config)
            sensitivities[name] = fit_global_moderation(sens_cells)
        context = (
            fit_context_diagnostic(primary_cells)
            if support["north_tropical_context_diagnostic_admitted"]
            else {"evaluable": False, "reason": "frozen_context_preflight_support_failed"}
        )
        classification = classify_trait_result(primary, sensitivities)
        trait_results[trait] = {
            "evaluable": bool(primary.get("evaluable")),
            "primary": primary,
            "sensitivities": sensitivities,
            "context_diagnostic": context,
            "decision": classification,
        }
    supported_count = sum(
        int(bool(item.get("decision", {}).get("supported"))) for item in trait_results.values()
    )
    family_supported = supported_count >= 2
    classification = (
        "reproductive_assurance_buffering_route_supported"
        if family_supported
        else "trait_specific_buffering_only"
        if supported_count == 1
        else "reproductive_assurance_buffering_not_supported"
    )
    result = {
        "contract": CONTRACT,
        "status": "completed",
        "preflight": preflight_result,
        "n_finite_matched_effect_rows": int(len(rows)),
        "trait_results": trait_results,
        "decision": {
            "supported_trait_count": int(supported_count),
            "route_A_family_supported": bool(family_supported),
            "classification": classification,
            "causal_mediation_identified": False,
        },
        "claim_ceiling": config["claim_ceiling"],
    }
    _write_json(output_dir / "RESULT.json", result)
    rows.to_csv(output_dir / "MATCHED_EFFECT_ROWS.csv.gz", index=False, compression="gzip")
    lines = [
        "# GloPL reproductive-assurance moderation result",
        "",
        f"- matched finite effect rows: {len(rows)}",
        f"- supported traits: {supported_count}/3",
        f"- route-A family supported: {family_supported}",
        f"- classification: {classification}",
    ]
    for trait, item in trait_results.items():
        primary = item.get("primary", {})
        lines.append(
            f"- {trait}: interaction={primary.get('distance_by_trait_interaction')}, "
            f"one-sided p={primary.get('interaction_one_sided_negative_p')}, "
            f"supported={item.get('decision', {}).get('supported', False)}"
        )
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    app()
