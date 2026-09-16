"""Post-hoc functional triangulation for the Chapter 1 v13 synthesis.

This module does not redefine the frozen Route A/B moderation tests. It reuses their
exact-species matched effect rows and frozen trait recodes to ask a different, explicitly
post-hoc question: are the predeclared response states associated with lower current
experimental pollen limitation after the same distance, context, measurement and
publication-weight controls?

The analysis is functional compatibility evidence, not historical causal mediation.
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

from island_v2.chapter1_h5_glopl_global_distance import (
    MEASUREMENT_COLUMNS,
    _clustered_wls,
    _context_dummies,
    _lower,
    _measurement_dummies,
    _p2,
    _truthy,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_v13_functional_bridge_v1"


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected v13 functional-bridge contract")
    if value.get("inferential_role") != "posthoc_functional_triangulation":
        raise typer.BadParameter("v13 functional bridge must remain post-hoc")
    return value


def aggregate_species_measurement_cells(
    rows: pd.DataFrame,
    *,
    publication_total_weight: float = 1.0,
) -> pd.DataFrame:
    """Aggregate duplicate GloPL rows exactly at the frozen Route A/B cell grain."""

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
    required = set(group_cols + ["PL_Effect_Size"])
    if missing := required - set(rows.columns):
        raise typer.BadParameter(f"matched effect rows missing columns: {sorted(missing)}")
    if rows.empty:
        return pd.DataFrame(columns=[*group_cols, "PL_Effect_Size", "n_effect_rows"])

    out = (
        rows.groupby(group_cols, as_index=False, dropna=False)
        .agg(
            PL_Effect_Size=("PL_Effect_Size", "mean"),
            n_effect_rows=("PL_Effect_Size", "size"),
        )
        .reset_index(drop=True)
    )
    counts = out.groupby("study_key")["site_key"].transform("size").astype(float)
    out["analysis_weight"] = float(publication_total_weight) / counts
    return out


def _trait_result_from_fit(
    cells: pd.DataFrame,
    fit: dict[str, Any],
    *,
    trait_name: str = "trait_state",
) -> dict[str, Any]:
    if not fit.get("evaluable"):
        return dict(fit)
    index = {name: i for i, name in enumerate(fit["names"])}
    if trait_name not in index:
        return {"evaluable": False, "reason": "trait_state_not_estimable"}
    i = index[trait_name]
    estimate = float(fit["beta"][i])
    se = float(fit["se"][i])
    z = estimate / se if se > 0 else float("nan")
    return {
        "evaluable": True,
        "trait_state_estimate": estimate,
        "trait_state_se": se,
        "trait_state_z": z,
        "trait_state_two_sided_p": _p2(z),
        "trait_state_one_sided_negative_p": _lower(z),
        "n_cells": int(len(cells)),
        "n_publications": int(cells["study_key"].nunique()),
        "n_sites": int(cells["site_key"].nunique()),
        "n_species": int(cells["species_key"].nunique()),
    }


def fit_global_trait_level(cells: pd.DataFrame) -> dict[str, Any]:
    """Fit PL ~ trait + distance + context + measurement fixed effects."""

    if cells.empty or cells["trait_state"].nunique() < 2:
        return {"evaluable": False, "reason": "both_trait_states_not_available"}
    distance = pd.to_numeric(cells["z_distance"], errors="coerce").to_numpy(float)
    trait_state = pd.to_numeric(cells["trait_state"], errors="coerce").to_numpy(float)
    context_cols, context_names, _ = _context_dummies(cells)
    measure_cols, measure_names = _measurement_dummies(cells)
    columns = [
        np.ones(len(cells), dtype=float),
        *context_cols,
        distance,
        trait_state,
        *measure_cols,
    ]
    names = [
        "intercept",
        *context_names,
        "z_distance",
        "trait_state",
        *measure_names,
    ]
    fit = _clustered_wls(cells, np.column_stack(columns), names)
    return _trait_result_from_fit(cells, fit)


def _mixed_group_frame(cells: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    if missing := set(group_columns) - set(cells.columns):
        raise typer.BadParameter(f"within-group columns missing: {sorted(missing)}")
    counts = (
        cells.groupby(group_columns, dropna=False)["trait_state"]
        .nunique()
        .rename("n_trait_states")
        .reset_index()
    )
    eligible = counts.loc[counts["n_trait_states"].ge(2), group_columns]
    if eligible.empty:
        return cells.iloc[0:0].copy()
    return cells.merge(eligible, on=group_columns, how="inner", validate="many_to_one")


def _weighted_demean(
    work: pd.DataFrame,
    values: np.ndarray,
    group_columns: list[str],
) -> np.ndarray:
    weights = pd.to_numeric(work["analysis_weight"], errors="coerce").to_numpy(float)
    if values.ndim == 1:
        values = values[:, None]
    out = values.astype(float, copy=True)
    group_indices = work.groupby(group_columns, dropna=False, sort=False).indices
    for indices in group_indices.values():
        idx = np.asarray(indices, dtype=int)
        local_weights = weights[idx]
        total_weight = float(local_weights.sum())
        if not np.isfinite(total_weight) or total_weight <= 0:
            raise typer.BadParameter("invalid within-group analysis weight")
        mean = (local_weights[:, None] * values[idx]).sum(axis=0) / total_weight
        out[idx] -= mean
    return out


def fit_within_group_trait_level(
    cells: pd.DataFrame,
    group_columns: list[str],
) -> dict[str, Any]:
    """Estimate the trait contrast after weighted fixed-effect demeaning."""

    work = _mixed_group_frame(cells, group_columns)
    if work.empty:
        return {"evaluable": False, "reason": "no_groups_with_both_trait_states"}

    distance = pd.to_numeric(work["z_distance"], errors="coerce").to_numpy(float)
    trait_state = pd.to_numeric(work["trait_state"], errors="coerce").to_numpy(float)
    context_cols, context_names, _ = _context_dummies(work)
    measure_cols, measure_names = _measurement_dummies(work)
    raw_columns = [*context_cols, distance, trait_state, *measure_cols]
    raw_names = [*context_names, "z_distance", "trait_state", *measure_names]
    raw_design = np.column_stack(raw_columns)
    design = _weighted_demean(work, raw_design, group_columns)
    y = _weighted_demean(
        work,
        pd.to_numeric(work["PL_Effect_Size"], errors="coerce").to_numpy(float),
        group_columns,
    )[:, 0]

    keep = np.nanstd(design, axis=0) > 1e-12
    names = [name for name, retained in zip(raw_names, keep, strict=True) if retained]
    design = design[:, keep]
    if "trait_state" not in names or design.size == 0:
        return {"evaluable": False, "reason": "trait_state_not_estimable_after_fixed_effects"}

    fit_frame = work.copy()
    fit_frame["PL_Effect_Size"] = y
    fit = _clustered_wls(fit_frame, design, names)
    result = _trait_result_from_fit(fit_frame, fit)
    result["n_fixed_effect_groups"] = int(
        work[group_columns].drop_duplicates().shape[0]
    )
    result["fixed_effect_group_columns"] = list(group_columns)
    return result


def _result_row(
    *,
    family: str,
    trait: str,
    analysis: str,
    result: dict[str, Any],
) -> dict[str, Any]:
    return {
        "family": family,
        "trait": trait,
        "analysis": analysis,
        "inferential_role": "posthoc_functional_triangulation",
        **result,
    }


def analyse_trait_rows(
    rows: pd.DataFrame,
    *,
    family: str,
    trait: str,
    publication_total_weight: float,
) -> list[dict[str, Any]]:
    part = rows.loc[rows["trait"].astype(str).eq(trait)].copy()
    if part.empty:
        return [
            _result_row(
                family=family,
                trait=trait,
                analysis="primary",
                result={"evaluable": False, "reason": "trait_absent_from_parent_rows"},
            )
        ]

    cells = aggregate_species_measurement_cells(
        part,
        publication_total_weight=publication_total_weight,
    )
    rows_out = [
        _result_row(
            family=family,
            trait=trait,
            analysis="primary",
            result=fit_global_trait_level(cells),
        )
    ]

    supplemental = part.loc[part["PL_Effect_Size_Type2"].astype(str).eq("Sup")].copy()
    supplemental_cells = aggregate_species_measurement_cells(
        supplemental,
        publication_total_weight=publication_total_weight,
    )
    rows_out.append(
        _result_row(
            family=family,
            trait=trait,
            analysis="supplemental_only",
            result=fit_global_trait_level(supplemental_cells),
        )
    )

    no_zero = part.loc[~_truthy(part["Constant_added"])].copy()
    no_zero_cells = aggregate_species_measurement_cells(
        no_zero,
        publication_total_weight=publication_total_weight,
    )
    rows_out.append(
        _result_row(
            family=family,
            trait=trait,
            analysis="no_zero_constant",
            result=fit_global_trait_level(no_zero_cells),
        )
    )

    rows_out.append(
        _result_row(
            family=family,
            trait=trait,
            analysis="within_publication",
            result=fit_within_group_trait_level(cells, ["study_key"]),
        )
    )
    rows_out.append(
        _result_row(
            family=family,
            trait=trait,
            analysis="within_publication_site",
            result=fit_within_group_trait_level(cells, ["study_key", "site_key"]),
        )
    )
    return rows_out


def run_functional_bridge(
    reproductive_assurance_rows: pd.DataFrame,
    floral_architecture_rows: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    publication_total_weight = float(config["analysis"]["publication_total_weight"])
    parent_frames = {
        "reproductive_assurance": reproductive_assurance_rows,
        "floral_architecture": floral_architecture_rows,
    }
    result_rows: list[dict[str, Any]] = []
    for family, trait_spec in config["traits"].items():
        parent = parent_frames[family]
        for trait in [str(x) for x in trait_spec["evaluable"]]:
            result_rows.extend(
                analyse_trait_rows(
                    parent,
                    family=family,
                    trait=trait,
                    publication_total_weight=publication_total_weight,
                )
            )
        for trait in [str(x) for x in trait_spec["support_failed_parent_traits"]]:
            result_rows.append(
                _result_row(
                    family=family,
                    trait=trait,
                    analysis="parent_support_gate",
                    result={
                        "evaluable": False,
                        "reason": "parent_support_gate_failed",
                    },
                )
            )

    results = pd.DataFrame(result_rows)
    primary = results.loc[results["analysis"].eq("primary")].copy()
    primary_supported = primary.loc[
        primary["evaluable"].astype(bool)
        & pd.to_numeric(
            primary["trait_state_one_sided_negative_p"], errors="coerce"
        ).le(0.05)
        & pd.to_numeric(primary["trait_state_estimate"], errors="coerce").lt(0)
    ]
    manifest = {
        "contract": CONTRACT,
        "inferential_role": "posthoc_functional_triangulation",
        "n_result_rows": int(len(results)),
        "n_primary_traits_evaluable": int(primary["evaluable"].astype(bool).sum()),
        "n_primary_traits_negative_directional_p_le_0_05": int(len(primary_supported)),
        "historical_pollen_limitation_selected_traits": False,
        "trait_mediates_isolation_effect": False,
        "failed_parent_moderation_tests_reclassified": False,
    }
    return results, manifest


def _summary_markdown(results: pd.DataFrame, manifest: dict[str, Any]) -> str:
    lines = [
        "# Chapter 1 v13 functional bridge",
        "",
        "Inferential role: **post-hoc functional triangulation**.",
        "",
        "| family | trait | analysis | estimate | one-sided negative p | evaluable |",
        "|---|---|---|---:|---:|---|",
    ]
    for row in results.to_dict("records"):
        estimate = row.get("trait_state_estimate")
        p_value = row.get("trait_state_one_sided_negative_p")
        estimate_text = "" if pd.isna(estimate) else f"{float(estimate):.6f}"
        p_text = "" if pd.isna(p_value) else f"{float(p_value):.6g}"
        lines.append(
            "| {family} | {trait} | {analysis} | {estimate} | {p_value} | {evaluable} |".format(
                family=row.get("family", ""),
                trait=row.get("trait", ""),
                analysis=row.get("analysis", ""),
                estimate=estimate_text,
                p_value=p_text,
                evaluable=bool(row.get("evaluable", False)),
            )
        )
    lines.extend(
        [
            "",
            f"Primary evaluable traits: {manifest['n_primary_traits_evaluable']}",
            (
                "Primary negative directional associations at p<=0.05: "
                f"{manifest['n_primary_traits_negative_directional_p_le_0_05']}"
            ),
            "",
            "These associations do not identify historical trait selection or causal mediation.",
        ]
    )
    return "\n".join(lines) + "\n"


@app.command("analyse")
def analyse(
    reproductive_assurance_rows: Path = typer.Option(..., exists=True, dir_okay=False),
    floral_architecture_rows: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    ra = pd.read_csv(reproductive_assurance_rows)
    architecture = pd.read_csv(floral_architecture_rows)
    results, manifest = run_functional_bridge(ra, architecture, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "functional_bridge_trait_results.csv", index=False)
    (output_dir / "functional_bridge_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "functional_bridge_summary.md").write_text(
        _summary_markdown(results, manifest),
        encoding="utf-8",
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
