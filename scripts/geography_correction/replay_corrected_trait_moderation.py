"""Replay frozen GloPL trait-moderation tests on corrected geography.

This is a measurement-repair replay, not a new trait-interaction search. It reuses
the exact matched effect rows, trait states, support decisions, model forms, family
rules, and sensitivities from the two 2026-09-16 frozen moderation analyses.
Only the geographic exposure and its standardization are replaced by the
24 September 2026 corrected spherical-coastline distances.

The script first reproduces the frozen original-geography estimates from the
archived matched rows. Corrected results are emitted only if that replay gate
passes to numerical tolerance.
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

from island_v2 import chapter1_h5_glopl_floral_architecture_moderation as fa
from island_v2 import chapter1_h5_glopl_reproductive_assurance_moderation as ra

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_yaml(path: Path) -> dict[str, Any]:
    obj = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise typer.BadParameter(f"invalid YAML object: {path}")
    return obj


def _load_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise typer.BadParameter(f"invalid JSON object: {path}")
    return obj


def _finite(value: object) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"non-finite value: {value}")
    return out


def _replace_distance(
    rows: pd.DataFrame,
    corrected_sites: pd.DataFrame,
    *,
    mean: float,
    sd: float,
) -> pd.DataFrame:
    work = rows.copy()
    distance = (
        corrected_sites[["site_key", "spherical_distance_km"]]
        .drop_duplicates("site_key")
        .set_index("site_key")["spherical_distance_km"]
    )
    mapped = work["site_key"].astype(str).map(distance)
    if mapped.isna().any():
        missing = work.loc[mapped.isna(), "site_key"].astype(str).drop_duplicates().tolist()[:10]
        raise typer.BadParameter(f"corrected distance missing for site keys: {missing}")
    work["distance_to_major_continent_km"] = pd.to_numeric(mapped, errors="raise")
    work["log1p_distance_to_major_continent_km"] = np.log1p(
        work["distance_to_major_continent_km"].to_numpy(float)
    )
    work["z_distance"] = (
        work["log1p_distance_to_major_continent_km"].to_numpy(float) - float(mean)
    ) / float(sd)
    return work


def _ra_replay(
    rows: pd.DataFrame,
    config: dict[str, Any],
    frozen: dict[str, Any],
) -> dict[str, Any]:
    admitted = set(frozen["outcome_blind_preflight"]["admitted_traits"])
    trait_results: dict[str, Any] = {}
    for trait in config["trait_recodes"]:
        if trait not in admitted:
            trait_results[trait] = {
                "evaluable": False,
                "reason": "frozen_preflight_support_failed",
            }
            continue
        part = rows.loc[rows["trait"].astype(str).eq(trait)].copy()
        primary_cells = ra.aggregate_species_measurement_cells(part, config)
        primary = ra.fit_global_moderation(primary_cells)
        sensitivities: dict[str, Any] = {}
        masks = {
            "supplemental_only": part["PL_Effect_Size_Type2"].astype(str).eq("Sup"),
            "no_zero_constant": ~part["Constant_added_bool"].astype(bool),
        }
        for name, mask in masks.items():
            cells = ra.aggregate_species_measurement_cells(part.loc[mask].copy(), config)
            sensitivities[name] = ra.fit_global_moderation(cells)
        support = frozen["outcome_blind_preflight"][trait]
        context = (
            ra.fit_context_diagnostic(primary_cells)
            if bool(support.get("north_tropical_context_diagnostic_admitted"))
            else {"evaluable": False, "reason": "frozen_context_preflight_support_failed"}
        )
        decision = ra.classify_trait_result(primary, sensitivities)
        trait_results[trait] = {
            "evaluable": bool(primary.get("evaluable")),
            "primary": primary,
            "sensitivities": sensitivities,
            "context_diagnostic": context,
            "decision": decision,
        }
    supported = sum(
        int(bool(x.get("decision", {}).get("supported"))) for x in trait_results.values()
    )
    return {
        "trait_results": trait_results,
        "decision": {
            "supported_trait_count": int(supported),
            "route_A_family_supported": bool(supported >= 2),
            "classification": (
                "reproductive_assurance_buffering_route_supported"
                if supported >= 2
                else "trait_specific_reproductive_assurance_buffering_only"
                if supported == 1
                else "reproductive_assurance_buffering_not_supported"
            ),
        },
    }


def _fa_replay(
    rows: pd.DataFrame,
    config: dict[str, Any],
    frozen: dict[str, Any],
) -> dict[str, Any]:
    admitted = set(frozen["outcome_blind_preflight"]["admitted_traits"])
    trait_results: dict[str, Any] = {}
    support_flags: list[bool] = []
    for trait in config["trait_recodes"]:
        if trait not in admitted:
            trait_results[trait] = {
                "evaluable": False,
                "reason": "frozen_preflight_support_failed",
            }
            support_flags.append(False)
            continue
        part = rows.loc[rows["trait"].astype(str).eq(trait)].copy()
        primary_cells = ra.aggregate_species_measurement_cells(part, config)
        primary = fa.fit_global_moderation(primary_cells)
        sensitivities: dict[str, Any] = {}
        masks = {
            "supplemental_only": part["PL_Effect_Size_Type2"].astype(str).eq("Sup"),
            "no_zero_constant": ~part["Constant_added_bool"].astype(bool),
        }
        for name, mask in masks.items():
            cells = ra.aggregate_species_measurement_cells(part.loc[mask].copy(), config)
            sensitivities[name] = fa.fit_global_moderation(cells)
        support = frozen["outcome_blind_preflight"][trait]
        context = (
            fa.fit_context_diagnostic(primary_cells)
            if bool(support.get("north_tropical_context_diagnostic_admitted"))
            else {"evaluable": False, "reason": "frozen_context_preflight_support_failed"}
        )
        decision = fa.classify_trait_result(primary, sensitivities)
        support_flags.append(bool(decision["supported"]))
        trait_results[trait] = {
            "evaluable": bool(primary.get("evaluable")),
            "primary": primary,
            "sensitivities": sensitivities,
            "context_diagnostic": context,
            "decision": decision,
        }
    family = fa.classify_family_result(support_flags)
    return {"trait_results": trait_results, "decision": family}


def _assert_close(observed: float, expected: float, *, label: str, atol: float = 1e-10) -> None:
    if not np.isclose(_finite(observed), _finite(expected), atol=atol, rtol=1e-8):
        raise AssertionError(f"{label}: observed={observed}, expected={expected}")


def _validate_original_ra(replayed: dict[str, Any], frozen: dict[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    for trait in frozen["outcome_blind_preflight"]["admitted_traits"]:
        got = replayed["trait_results"][trait]
        old = frozen["analysis"][trait]
        pairs = [
            ("interaction", got["primary"]["distance_by_trait_interaction"], old["distance_by_trait_interaction"]),
            ("interaction_se", got["primary"]["interaction_se"], old["interaction_se"]),
            ("two_sided_p", got["primary"]["interaction_two_sided_p"], old["interaction_two_sided_p"]),
            (
                "supplemental_only",
                got["sensitivities"]["supplemental_only"]["distance_by_trait_interaction"],
                old["supplemental_only_interaction"],
            ),
            (
                "no_zero_constant",
                got["sensitivities"]["no_zero_constant"]["distance_by_trait_interaction"],
                old["no_zero_constant_interaction"],
            ),
        ]
        for suffix, observed, expected in pairs:
            _assert_close(observed, expected, label=f"RA {trait} {suffix}")
            checks.append({"trait": trait, "quantity": suffix, "absolute_difference": abs(float(observed) - float(expected))})
    return {"status": "pass", "checks": checks}


def _validate_original_fa(replayed: dict[str, Any], frozen: dict[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    for trait in frozen["outcome_blind_preflight"]["admitted_traits"]:
        got = replayed["trait_results"][trait]
        old = frozen["analysis"][trait]
        primary_old = old["primary"]
        pairs = [
            (
                "interaction",
                got["primary"]["distance_by_architecture_interaction"],
                primary_old["distance_by_architecture_interaction"],
            ),
            ("interaction_se", got["primary"]["interaction_se"], primary_old["interaction_se"]),
            ("two_sided_p", got["primary"]["interaction_two_sided_p"], primary_old["interaction_two_sided_p"]),
            (
                "supplemental_only",
                got["sensitivities"]["supplemental_only"]["distance_by_architecture_interaction"],
                old["sensitivities"]["supplemental_only_interaction"],
            ),
            (
                "no_zero_constant",
                got["sensitivities"]["no_zero_constant"]["distance_by_architecture_interaction"],
                old["sensitivities"]["no_zero_constant_interaction"],
            ),
        ]
        for suffix, observed, expected in pairs:
            _assert_close(observed, expected, label=f"FA {trait} {suffix}")
            checks.append({"trait": trait, "quantity": suffix, "absolute_difference": abs(float(observed) - float(expected))})
    return {"status": "pass", "checks": checks}


def _compact_rows(
    family: str,
    original: dict[str, Any],
    corrected: dict[str, Any],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for trait, current in corrected["trait_results"].items():
        if not current.get("evaluable"):
            rows.append(
                {
                    "family": family,
                    "trait": trait,
                    "evaluable": False,
                    "reason": current.get("reason", "not_evaluable"),
                }
            )
            continue
        old = original["trait_results"][trait]
        key = (
            "distance_by_trait_interaction"
            if family == "reproductive_assurance"
            else "distance_by_architecture_interaction"
        )
        rows.append(
            {
                "family": family,
                "trait": trait,
                "evaluable": True,
                "original_interaction": old["primary"][key],
                "corrected_interaction": current["primary"][key],
                "corrected_interaction_se": current["primary"]["interaction_se"],
                "corrected_two_sided_p": current["primary"]["interaction_two_sided_p"],
                "corrected_one_sided_buffering_p": current["primary"][
                    "interaction_one_sided_negative_p"
                ],
                "corrected_supplemental_only_interaction": current["sensitivities"][
                    "supplemental_only"
                ][key],
                "corrected_no_zero_constant_interaction": current["sensitivities"][
                    "no_zero_constant"
                ][key],
                "supported_under_frozen_rule": current["decision"]["supported"],
            }
        )
    return rows


@app.command("run")
def run(
    ra_rows_csv: Path = typer.Option(..., exists=True),
    fa_rows_csv: Path = typer.Option(..., exists=True),
    corrected_site_distances_csv: Path = typer.Option(..., exists=True),
    h3_comparison_json: Path = typer.Option(..., exists=True),
    ra_config_path: Path = typer.Option(..., exists=True),
    fa_config_path: Path = typer.Option(..., exists=True),
    ra_result_lock_path: Path = typer.Option(..., exists=True),
    fa_result_lock_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    ra_config = _load_yaml(ra_config_path)
    fa_config = _load_yaml(fa_config_path)
    ra_lock = _load_json(ra_result_lock_path)
    fa_lock = _load_json(fa_result_lock_path)
    h3 = _load_json(h3_comparison_json)
    scaling = h3["corrected"]["distance_standardization"]
    mean = _finite(scaling["mean"])
    sd = _finite(scaling["sd"])

    ra_rows = pd.read_csv(ra_rows_csv)
    fa_rows = pd.read_csv(fa_rows_csv)
    sites = pd.read_csv(corrected_site_distances_csv)
    if sites["site_key"].duplicated().any():
        raise typer.BadParameter("corrected site-distance table contains duplicate site_key")

    original_ra = _ra_replay(ra_rows, ra_config, ra_lock)
    original_fa = _fa_replay(fa_rows, fa_config, fa_lock)
    gate = {
        "reproductive_assurance": _validate_original_ra(original_ra, ra_lock),
        "floral_architecture": _validate_original_fa(original_fa, fa_lock),
    }

    corrected_ra_rows = _replace_distance(ra_rows, sites, mean=mean, sd=sd)
    corrected_fa_rows = _replace_distance(fa_rows, sites, mean=mean, sd=sd)
    corrected_ra = _ra_replay(corrected_ra_rows, ra_config, ra_lock)
    corrected_fa = _fa_replay(corrected_fa_rows, fa_config, fa_lock)

    compact = pd.DataFrame(
        _compact_rows("reproductive_assurance", original_ra, corrected_ra)
        + _compact_rows("floral_architecture", original_fa, corrected_fa)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    compact.to_csv(output_dir / "corrected_moderation_summary.csv", index=False)

    result = {
        "contract": "chapter1_corrected_trait_moderation_replay_v1",
        "inferential_role": "measurement_repair_replay_of_frozen_20260916_moderation_tests",
        "rules": [
            "no_new_trait_contrasts",
            "no_support_threshold_changes",
            "no_model_changes",
            "same_matched_effect_rows",
            "same_trait_states",
            "same_sensitivities",
            "only_distance_exposure_and_standardization_replaced",
            "historical_mediation_not_identified",
        ],
        "corrected_distance_standardization": {"mean": mean, "sd": sd},
        "original_reproduction_gate": gate,
        "reproductive_assurance": {
            "original": original_ra,
            "corrected": corrected_ra,
        },
        "floral_architecture": {
            "original": original_fa,
            "corrected": corrected_fa,
        },
        "decision": {
            "reproductive_assurance_family_supported": bool(
                corrected_ra["decision"]["route_A_family_supported"]
            ),
            "floral_architecture_family_supported": bool(
                corrected_fa["decision"]["route_B_family_supported"]
            ),
            "causal_mediation_identified": False,
        },
    }
    (output_dir / "RESULT.json").write_text(
        json.dumps(result, indent=2, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    app()
