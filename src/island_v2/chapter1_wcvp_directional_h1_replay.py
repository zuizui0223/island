"""Replay the frozen one-dimensional H1 on WCVP regional-native-compatible flora.

This module changes only the floristic provenance surface. It reuses:
- the frozen seven H1 indicators and weights;
- the corrected geography and covariates;
- the final finite-cluster t and Rademacher sign-flip inference;
- the H1a Paule-Mandel / modified Hartung-Knapp synthesis.

WCVP compatibility is regional (TDWG level 3), not exact island nativeness.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_h1_final_directional import run_final_directional_h1
from island_v2.chapter1_wcvp_native_compatibility import (
    build_regional_native_compatible_flora,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def regional_native_h1_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = copy.deepcopy(config)
    cfg["strata"] = ["all_native"]
    cfg["flora_roles"]["broad_primary"] = "all_native"
    cfg["flora_roles"]["status_sensitivities"] = []
    return cfg


def _relabel_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "stratum" in out.columns:
        out["stratum"] = out["stratum"].replace(
            {"all_native": "regional_native_compatible"}
        )
    return out


def _relabel_mapping(value: dict[str, Any]) -> dict[str, Any]:
    out = dict(value)
    if out.get("stratum") == "all_native":
        out["stratum"] = "regional_native_compatible"
    return out


def run_replay(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    config: dict[str, Any],
    *,
    evidence_scope: str,
) -> dict[str, Any]:
    upgraded, audit = build_regional_native_compatible_flora(
        status_flora,
        wcvp_ranges,
        island_tdwg,
    )
    result = run_final_directional_h1(
        upgraded,
        state_audit,
        covariates,
        regional_native_h1_config(config),
        evidence_scope=evidence_scope,
    )
    frames = {
        key: _relabel_frame(value)
        for key, value in result.items()
        if isinstance(value, pd.DataFrame)
    }
    mappings = {
        key: _relabel_mapping(value)
        for key, value in result.items()
        if isinstance(value, dict)
    }
    return {
        **frames,
        **mappings,
        "coverage": audit,
    }


@app.command()
def main(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    wcvp_ranges_csv: Path = typer.Option(..., exists=True),
    island_tdwg_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    result = run_replay(
        pd.read_csv(status_flora_csv, dtype=str).fillna(""),
        pd.read_csv(state_audit_csv, dtype=str).fillna(""),
        pd.read_csv(covariates_csv),
        pd.read_csv(wcvp_ranges_csv, dtype=str).fillna(""),
        pd.read_csv(island_tdwg_csv, dtype=str).fillna(""),
        config,
        evidence_scope=evidence_scope,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    result["results"].to_csv(
        output_dir / "regional_directional_score_results.csv",
        index=False,
    )
    result["components"].to_csv(
        output_dir / "regional_directional_score_components.csv",
        index=False,
    )
    pd.DataFrame([result["coverage"]]).to_csv(
        output_dir / "regional_native_coverage.csv",
        index=False,
    )
    for key in (
        "iut",
        "synthesis",
        "top_cluster_iut",
        "top_cluster_synthesis",
    ):
        (output_dir / f"{key}.json").write_text(
            json.dumps(result[key], indent=2) + "\n",
            encoding="utf-8",
        )
    result["top_cluster_sensitivity"].to_csv(
        output_dir / "top_cluster_leaveout_results.csv",
        index=False,
    )

    summary = {
        "evidence_scope": evidence_scope,
        "coverage": result["coverage"],
        "regional_results": result["results"].loc[
            result["results"]["stratum"].eq(
                "regional_native_compatible"
            ),
            [
                "context",
                "status",
                "estimate",
                "cluster_robust_se",
                "p_one_sided_t",
                "p_one_sided_wild",
                "n_unique_islands",
                "n_clusters",
                "effective_clusters",
            ],
        ].to_dict(orient="records"),
        "strict_four_region_iut": result["iut"],
        "global_synthesis": result["synthesis"],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
