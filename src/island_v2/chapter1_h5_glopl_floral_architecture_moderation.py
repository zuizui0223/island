"""Species-level floral-architecture moderation of the global GloPL distance gradient.

The outcome-blind preflight exact-matches GloPL plant names to three frozen atomic
floral-architecture contrasts. GloPL effect columns are read only after overlap/support
is written. The analysis is a functional association test, not causal mediation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_global_distance import (
    CONTEXTS,
    build_preflight_table,
    load_config as load_parent_config,
)
from island_v2.chapter1_h5_glopl_reproductive_assurance_moderation import (
    _prepare_effect_rows,
    _read_effects,
    aggregate_species_measurement_cells,
    fit_context_diagnostic as _fit_context_trait_diagnostic,
    fit_global_moderation as _fit_global_trait_moderation,
    normalize_species_name,
    summarize_trait_support,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_glopl_floral_architecture_moderation_v1"
PARENT_CONFIG_PATH = Path("config/chapter1_h5_glopl_global_distance_v1.yml")


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected GloPL floral-architecture moderation contract")
    return config


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
        ledger["axis"].astype(str).eq("floral_structural_complexity")
        & ledger["resolution_status"].astype(str).eq("resolved")
        & ledger["quality"].astype(str).isin(["high", "medium"])
    ].copy()
    rows: list[dict[str, Any]] = []
    for trait, spec in config["trait_recodes"].items():
        source_trait = str(spec["source_trait"])
        accessible = {str(x) for x in spec["protected_states"]}
        restricted = {str(x) for x in spec["unprotected_states"]}
        part = work.loc[work["trait_name"].astype(str).eq(source_trait)].copy()
        for row in part.itertuples(index=False):
            value = str(row.normalized_value)
            if value in accessible:
                state = 1
            elif value in restricted:
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
        raise typer.BadParameter(f"conflicting floral-architecture assignments: {bad}")
    return out.drop_duplicates(["species_key", "trait"]).reset_index(drop=True)


def fit_global_moderation(cells: pd.DataFrame) -> dict[str, Any]:
    raw = _fit_global_trait_moderation(cells)
    if not raw.get("evaluable"):
        return raw
    return {
        "evaluable": True,
        "restricted_distance_slope": raw["unprotected_distance_slope"],
        "accessible_distance_slope": raw["protected_distance_slope"],
        "distance_by_architecture_interaction": raw["distance_by_trait_interaction"],
        "interaction_se": raw["interaction_se"],
        "interaction_two_sided_p": raw["interaction_two_sided_p"],
        "interaction_one_sided_negative_p": raw["interaction_one_sided_negative_p"],
        "n_cells": raw["n_cells"],
        "n_publications": raw["n_publications"],
        "n_sites": raw["n_sites"],
        "n_species": raw["n_species"],
    }


def fit_context_diagnostic(cells: pd.DataFrame) -> dict[str, Any]:
    raw = _fit_context_trait_diagnostic(cells)
    if not raw.get("evaluable"):
        return raw
    return {
        "evaluable": True,
        "northern_distance_by_architecture_interaction": raw[
            "northern_distance_by_trait_interaction"
        ],
        "tropical_distance_by_architecture_interaction": raw[
            "tropical_distance_by_trait_interaction"
        ],
        "tropical_minus_northern_buffering_interaction": raw[
            "tropical_minus_northern_buffering_interaction"
        ],
        "interaction_two_sided_p": raw["interaction_two_sided_p"],
        "n_cells": raw["n_cells"],
        "n_publications": raw["n_publications"],
        "n_sites": raw["n_sites"],
        "n_species": raw["n_species"],
    }


def classify_trait_result(
    primary: dict[str, Any], sensitivities: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    sensitivity_direction = bool(sensitivities) and all(
        float(result.get("distance_by_architecture_interaction", float("nan"))) < 0
        for result in sensitivities.values()
    )
    supported = bool(
        primary.get("evaluable")
        and float(primary.get("restricted_distance_slope", float("nan"))) > 0
        and float(primary.get("distance_by_architecture_interaction", float("nan"))) < 0
        and float(primary.get("interaction_one_sided_negative_p", 1.0)) <= 0.05
        and sensitivity_direction
    )
    return {
        "supported": supported,
        "sensitivity_direction_retained": bool(sensitivity_direction),
    }


def classify_family_result(supported: list[bool]) -> dict[str, Any]:
    count = int(sum(bool(x) for x in supported))
    return {
        "supported_trait_count": count,
        "route_B_family_supported": count >= 2,
        "classification": (
            "floral_architecture_buffering_route_supported"
            if count >= 2
            else "trait_specific_architecture_buffering_only"
            if count == 1
            else "floral_architecture_buffering_not_supported"
        ),
    }


def _read_glopl_preflight(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    columns = [str(x) for x in config["preflight"]["glopl_allowed_columns"]]
    return pd.read_csv(path, usecols=columns, dtype=str, encoding="latin-1").fillna("")


def _read_trait_ledger(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    columns = [str(x) for x in config["preflight"]["trait_allowed_columns"]]
    return pd.read_csv(path, usecols=columns, dtype=str).fillna("")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8")


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
        "status": "completed_outcome_blind_architecture_overlap_preflight",
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
        "# GloPL floral-architecture overlap preflight",
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
                "route_B_family_supported": False,
                "classification": "floral_architecture_buffering_not_evaluable",
            },
            "claim_ceiling": config["claim_ceiling"],
        }
        _write_json(output_dir / "RESULT.json", result)
        (output_dir / "RESULT_SUMMARY.md").write_text(
            "# GloPL floral-architecture moderation\n\n- no trait passed the frozen preflight gate\n",
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
    support_flags: list[bool] = []
    for trait in config["trait_recodes"]:
        support = preflight_result["trait_support"][trait]
        if trait not in admitted:
            trait_results[trait] = {"evaluable": False, "reason": "preflight_support_failed"}
            support_flags.append(False)
            continue
        trait_rows = rows.loc[rows["trait"].eq(trait)].copy()
        primary_cells = aggregate_species_measurement_cells(trait_rows, config)
        primary = fit_global_moderation(primary_cells)
        sensitivities: dict[str, Any] = {}
        for name, mask in {
            "supplemental_only": trait_rows["PL_Effect_Size_Type2"].astype(str).eq("Sup"),
            "no_zero_constant": ~trait_rows["Constant_added_bool"].astype(bool),
        }.items():
            sens_cells = aggregate_species_measurement_cells(trait_rows.loc[mask].copy(), config)
            sensitivities[name] = fit_global_moderation(sens_cells)
        context = (
            fit_context_diagnostic(primary_cells)
            if support["north_tropical_context_diagnostic_admitted"]
            else {"evaluable": False, "reason": "frozen_context_preflight_support_failed"}
        )
        decision = classify_trait_result(primary, sensitivities)
        support_flags.append(bool(decision["supported"]))
        trait_results[trait] = {
            "evaluable": bool(primary.get("evaluable")),
            "primary": primary,
            "sensitivities": sensitivities,
            "context_diagnostic": context,
            "decision": decision,
        }

    family = classify_family_result(support_flags)
    result = {
        "contract": CONTRACT,
        "status": "completed",
        "preflight": preflight_result,
        "n_finite_matched_effect_rows": int(len(rows)),
        "trait_results": trait_results,
        "decision": {
            **family,
            "causal_mediation_identified": False,
        },
        "claim_ceiling": config["claim_ceiling"],
    }
    _write_json(output_dir / "RESULT.json", result)
    rows.to_csv(output_dir / "MATCHED_EFFECT_ROWS.csv.gz", index=False, compression="gzip")
    lines = [
        "# GloPL floral-architecture moderation result",
        "",
        f"- matched finite effect rows: {len(rows)}",
        f"- supported traits: {family['supported_trait_count']}/3",
        f"- route-B family supported: {family['route_B_family_supported']}",
        f"- classification: {family['classification']}",
    ]
    for trait, item in trait_results.items():
        primary = item.get("primary", {})
        lines.append(
            f"- {trait}: interaction={primary.get('distance_by_architecture_interaction')}, "
            f"one-sided p={primary.get('interaction_one_sided_negative_p')}, "
            f"supported={item.get('decision', {}).get('supported', False)}"
        )
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    app()
