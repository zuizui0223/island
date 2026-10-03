"""Replay the final one-dimensional H1 on the frozen WCVP regional-native-compatible flora.

This is an exploratory provenance replay. It does not redefine the confirmatory H1.
It reuses:
- the exact seven frozen H1 indicators and weights;
- the exact final finite-cluster directional-score engine;
- the frozen WCVP regional-native compatibility rule.

Only the floristic sampling frame is changed.
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

EXPECTED_ROWS = 513_320
EXPECTED_ISLANDS = 2_372
STRATUM = "regional_native_compatible"


def build_replay_config(parent: dict[str, Any]) -> dict[str, Any]:
    cfg = copy.deepcopy(parent)
    cfg["strata"] = ["all_native"]
    cfg["flora_roles"] = {
        "broad_primary": "all_native",
        "status_sensitivities": [],
    }
    return cfg


def relabel_stratum(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "stratum" in out.columns:
        out.loc[out["stratum"].astype(str).eq("all_native"), "stratum"] = STRATUM
    return out


def run_replay(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    parent_config: dict[str, Any],
    *,
    evidence_scope: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    upgraded, audit = build_regional_native_compatible_flora(
        status_flora,
        wcvp_ranges,
        island_tdwg,
    )
    if int(audit["n_regional_native_rows"]) != EXPECTED_ROWS:
        raise typer.BadParameter(
            "WCVP regional-native row count drift: "
            f'{audit["n_regional_native_rows"]} != {EXPECTED_ROWS}'
        )
    if int(audit["n_regional_native_islands"]) != EXPECTED_ISLANDS:
        raise typer.BadParameter(
            "WCVP regional-native island count drift: "
            f'{audit["n_regional_native_islands"]} != {EXPECTED_ISLANDS}'
        )

    cfg = build_replay_config(parent_config)
    result = run_final_directional_h1(
        upgraded,
        state_audit,
        covariates,
        cfg,
        evidence_scope=evidence_scope,
    )

    for key in (
        "results",
        "components",
        "counts",
        "top_cluster_sensitivity",
        "weight_sensitivity",
        "weight_sensitivity_summary",
    ):
        if isinstance(result.get(key), pd.DataFrame):
            result[key] = relabel_stratum(result[key])

    for key in (
        "iut",
        "meta",
        "synthesis",
        "top_cluster_iut",
        "top_cluster_meta",
        "top_cluster_synthesis",
    ):
        if isinstance(result.get(key), dict):
            item = dict(result[key])
            if item.get("stratum") == "all_native":
                item["stratum"] = STRATUM
            result[key] = item

    return result, audit


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
    parent = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    result, audit = run_replay(
        pd.read_csv(status_flora_csv, dtype=str).fillna(""),
        pd.read_csv(state_audit_csv, dtype=str).fillna(""),
        pd.read_csv(covariates_csv),
        pd.read_csv(wcvp_ranges_csv),
        pd.read_csv(island_tdwg_csv),
        parent,
        evidence_scope=evidence_scope,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    result["results"].to_csv(output_dir / "directional_score_results.csv", index=False)
    result["components"].to_csv(output_dir / "directional_score_components.csv", index=False)
    result["counts"].to_csv(
        output_dir / "directional_score_counts.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    result["top_cluster_sensitivity"].to_csv(
        output_dir / "top_cluster_leaveout_results.csv", index=False
    )
    result["weight_sensitivity"].to_csv(
        output_dir / "weight_sensitivity_results.csv", index=False
    )
    result["weight_sensitivity_summary"].to_csv(
        output_dir / "weight_sensitivity_summary.csv", index=False
    )

    for key, filename in (
        ("iut", "directional_score_iut.json"),
        ("meta", "directional_score_meta.json"),
        ("synthesis", "directional_score_synthesis.json"),
        ("top_cluster_iut", "top_cluster_leaveout_iut.json"),
        ("top_cluster_meta", "top_cluster_leaveout_meta.json"),
        ("top_cluster_synthesis", "top_cluster_leaveout_synthesis.json"),
    ):
        (output_dir / filename).write_text(
            json.dumps(result[key], indent=2) + "\n",
            encoding="utf-8",
        )

    summary = {
        "analysis": "final H1 directional-score replay on WCVP regional-native-compatible flora",
        "inferential_role": "exploratory_provenance_replay",
        "evidence_scope": evidence_scope,
        "wcvp_audit": audit,
        "regional_results": relabel_stratum(result["results"]).to_dict(orient="records"),
        "iut": result["iut"],
        "synthesis": result["synthesis"],
    }
    (output_dir / "SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
