"""P1a audit of paired support across observed/family/genus taxonomic stages.

This module does not refit any biological model. It audits the frozen Chapter 1
V2 taxonomic-depth outputs and fails closed if stage-specific sample changes or
algebraic inconsistencies are detected.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

STAGES = ("observed_score", "after_family_residual", "after_genus_residual")
SCOPES = ("direct_only", "all_analysis_eligible")


class P1AuditError(ValueError):
    """Raised when frozen P1a support invariants fail."""


def _read_csv_any(root: Path, stem: str) -> pd.DataFrame:
    for suffix in (".csv", ".csv.gz"):
        path = root / f"{stem}{suffix}"
        if path.is_file():
            return pd.read_csv(path)
    raise FileNotFoundError(root / f"{stem}.csv[.gz]")


def load_config(path: Path) -> dict[str, Any]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or cfg.get("contract") != "chapter1_p1_assembly_depth_defense_v1":
        raise P1AuditError("unexpected P1 defense contract")
    return cfg


def audit_scope(root: Path, scope: str, cfg: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    scope_root = root / "taxonomic-depth" / scope
    decomposition = _read_csv_any(scope_root, "decomposition")
    long_scores = _read_csv_any(scope_root, "long_scores")

    required_decomp = {
        "island_id",
        "syndrome",
        "stratum",
        "source_mode",
        "observed_score",
        "family_expected",
        "genus_expected",
        "after_family_residual",
        "after_genus_residual",
        "n_species",
    }
    missing = required_decomp - set(decomposition.columns)
    if missing:
        raise P1AuditError(f"{scope} decomposition missing columns: {sorted(missing)}")

    required_long = {
        "island_id",
        "syndrome",
        "stratum",
        "source_mode",
        "n_species",
        "taxonomic_stage",
        "syndrome_score",
    }
    missing = required_long - set(long_scores.columns)
    if missing:
        raise P1AuditError(f"{scope} long_scores missing columns: {sorted(missing)}")

    for column in (
        "observed_score",
        "family_expected",
        "genus_expected",
        "after_family_residual",
        "after_genus_residual",
        "n_species",
    ):
        decomposition[column] = pd.to_numeric(decomposition[column], errors="coerce")
    if decomposition[list(required_decomp - {"island_id", "syndrome", "stratum", "source_mode"})].isna().any().any():
        raise P1AuditError(f"{scope} decomposition contains missing retained-stage values")

    family_error = np.abs(
        decomposition["after_family_residual"]
        - (decomposition["observed_score"] - decomposition["family_expected"])
    )
    genus_error = np.abs(
        decomposition["after_genus_residual"]
        - (decomposition["observed_score"] - decomposition["genus_expected"])
    )
    tolerance = 1e-10
    if float(family_error.max()) > tolerance or float(genus_error.max()) > tolerance:
        raise P1AuditError(f"{scope} residual algebra does not reproduce frozen decomposition")

    expected_stages = set(STAGES)
    long_scores["taxonomic_stage"] = long_scores["taxonomic_stage"].astype(str)
    if not set(long_scores["taxonomic_stage"]).issubset(expected_stages):
        unexpected = sorted(set(long_scores["taxonomic_stage"]) - expected_stages)
        raise P1AuditError(f"{scope} unexpected taxonomic stages: {unexpected}")

    key = ["island_id", "syndrome", "stratum", "source_mode"]
    decomp_keys = decomposition[key + ["n_species"]].copy()
    decomp_keys["island_id"] = decomp_keys["island_id"].astype(str)
    decomp_keys["syndrome"] = decomp_keys["syndrome"].astype(str)
    decomp_keys["stratum"] = decomp_keys["stratum"].astype(str)
    decomp_keys["source_mode"] = decomp_keys["source_mode"].astype(str)
    if decomp_keys.duplicated(key).any():
        raise P1AuditError(f"{scope} duplicate decomposition keys")

    stage_sets: dict[str, set[tuple[str, str, str, str]]] = {}
    stage_species_weights: dict[str, pd.Series] = {}
    for stage in STAGES:
        part = long_scores.loc[long_scores["taxonomic_stage"].eq(stage)].copy()
        for col in key:
            part[col] = part[col].astype(str)
        if part.duplicated(key).any():
            raise P1AuditError(f"{scope} duplicate long-score keys at {stage}")
        stage_sets[stage] = set(part[key].itertuples(index=False, name=None))
        stage_species_weights[stage] = part.set_index(key)["n_species"].sort_index()

    decomp_set = set(decomp_keys[key].itertuples(index=False, name=None))
    for stage in STAGES:
        if stage_sets[stage] != decomp_set:
            raise P1AuditError(f"{scope} stage-specific support differs at {stage}")

    ref_weight = stage_species_weights[STAGES[0]]
    for stage in STAGES[1:]:
        current = stage_species_weights[stage]
        if not ref_weight.index.equals(current.index):
            raise P1AuditError(f"{scope} n_species index differs at {stage}")
        if not np.allclose(
            pd.to_numeric(ref_weight, errors="coerce"),
            pd.to_numeric(current, errors="coerce"),
            equal_nan=False,
        ):
            raise P1AuditError(f"{scope} n_species weights differ at {stage}")

    target = cfg["primary_context"]
    axes = {str(x) for x in target["response_axes"]}
    strata = {str(x) for x in target["strata"]}
    source_modes = {str(x) for x in target["source_modes"]}
    focal = decomposition.loc[
        decomposition["syndrome"].astype(str).isin(axes)
        & decomposition["stratum"].astype(str).isin(strata)
        & decomposition["source_mode"].astype(str).isin(source_modes)
    ].copy()

    rows: list[dict[str, Any]] = []
    group_cols = ["stratum", "source_mode", "syndrome"]
    for group_key, group in focal.groupby(group_cols, sort=True):
        rows.append(
            {
                "evidence_scope": scope,
                **dict(zip(group_cols, group_key, strict=True)),
                "n_rows": int(len(group)),
                "n_islands": int(group["island_id"].astype(str).nunique()),
                "median_n_species": float(pd.to_numeric(group["n_species"]).median()),
                "min_n_species": int(pd.to_numeric(group["n_species"]).min()),
                "max_family_algebra_error": float(family_error.loc[group.index].max()),
                "max_genus_algebra_error": float(genus_error.loc[group.index].max()),
                "same_support_all_three_stages": True,
                "same_n_species_weight_all_three_stages": True,
            }
        )

    summary = {
        "evidence_scope": scope,
        "n_decomposition_rows": int(len(decomposition)),
        "n_unique_decomposition_keys": int(len(decomp_set)),
        "same_support_all_three_stages": True,
        "same_n_species_weight_all_three_stages": True,
        "max_family_algebra_error": float(family_error.max()),
        "max_genus_algebra_error": float(genus_error.max()),
        "family_and_genus_are_grouping_not_trait_imputation": True,
    }
    return pd.DataFrame(rows), summary


def run_audit(artifact_root: Path, config_path: Path, output_dir: Path) -> dict[str, Any]:
    cfg = load_config(config_path)
    tables: list[pd.DataFrame] = []
    summaries: list[dict[str, Any]] = []
    for scope in SCOPES:
        table, summary = audit_scope(artifact_root, scope, cfg)
        tables.append(table)
        summaries.append(summary)

    output_dir.mkdir(parents=True, exist_ok=True)
    combined = pd.concat(tables, ignore_index=True)
    combined.to_csv(output_dir / "p1a_paired_support_cells.csv", index=False)

    direct = next(x for x in summaries if x["evidence_scope"] == "direct_only")
    manifest = {
        "contract": "chapter1_p1a_paired_support_audit_v1",
        "source_contract": cfg["contract"],
        "source_workflow_run_id": int(cfg["pinned_primary_artifact"]["workflow_run_id"]),
        "source_artifact_id": int(cfg["pinned_primary_artifact"]["artifact_id"]),
        "source_artifact_digest": str(cfg["pinned_primary_artifact"]["artifact_digest"]),
        "scope_summaries": summaries,
        "primary_direct_only_pass": bool(
            direct["same_support_all_three_stages"]
            and direct["same_n_species_weight_all_three_stages"]
            and direct["max_family_algebra_error"] <= 1e-10
            and direct["max_genus_algebra_error"] <= 1e-10
        ),
        "new_biological_model_fitted": False,
        "new_p_value_generated": False,
        "claim_boundary": (
            "This audit can exclude stage-specific support/weight changes within the frozen V2 "
            "decomposition. It does not establish that true genus boundaries outperform matched "
            "arbitrary partitions; that is P1c."
        ),
    }
    (output_dir / "chapter1_p1a_paired_support_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command()
def main(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(run_audit(artifact_root, config_path, output_dir), indent=2))


if __name__ == "__main__":
    app()
