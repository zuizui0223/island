"""Execute the frozen Chapter 1 N1 model on standardized channel evidence.

This module is intentionally an execution wrapper around the already frozen statistical
implementation in :mod:`island_v2.chapter1_nee_n1_model`.  It does not alter model
formulas, gates, channel definitions, or robustness rules.  Its only data adaptation is
an audited rename of the canonical Chapter 1 geography/covariate artifact into the
column names required by the N1 contract.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.chapter1_nee_n1_model import (
    deletion_robustness,
    fit_n1_primary,
    gate_receipt,
    load_config,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

RAW_COVARIATE_COLUMNS = {
    "island_id",
    "distance_to_continent_km",
    "log_distance_to_continent_km",
    "area_km2",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
    "spatial_block",
}
N1_COVARIATE_COLUMNS = [
    "island_id",
    "log1p_distance_to_continent_km",
    "log_area",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
    "spatial_block",
]


@app.callback()
def _root() -> None:
    """Prospective N1 execution utilities."""


def prepare_canonical_covariates(
    raw: pd.DataFrame,
    *,
    expected_islands: int | None = 8265,
    atol: float = 1e-10,
) -> pd.DataFrame:
    """Audit and rename the frozen Chapter 1 geography/covariate table.

    The historical canonical artifact names the two log-transformed geography columns
    ``log_distance_to_continent_km`` and ``log_island_area_km2``.  N1 uses the explicit
    names ``log1p_distance_to_continent_km`` and ``log_area``.  We verify that the
    historical columns are exactly the declared ``log1p`` transforms before renaming;
    no geography is recomputed or substituted.
    """
    missing = RAW_COVARIATE_COLUMNS.difference(raw.columns)
    if missing:
        raise ValueError(f"canonical Chapter 1 covariates missing columns: {sorted(missing)}")

    work = raw.copy()
    work["island_id"] = work["island_id"].fillna("").astype(str).str.strip()
    if work["island_id"].eq("").any():
        raise ValueError("canonical Chapter 1 covariates contain blank island_id")
    if work["island_id"].duplicated().any():
        raise ValueError("canonical Chapter 1 covariates must contain one row per island_id")
    if expected_islands is not None and len(work) != int(expected_islands):
        raise ValueError(
            f"canonical Chapter 1 covariate universe changed: expected {expected_islands}, "
            f"found {len(work)}"
        )

    distance = pd.to_numeric(work["distance_to_continent_km"], errors="coerce")
    distance_log = pd.to_numeric(work["log_distance_to_continent_km"], errors="coerce")
    area = pd.to_numeric(work["area_km2"], errors="coerce")
    area_log = pd.to_numeric(work["log_island_area_km2"], errors="coerce")

    if distance.isna().any() or distance_log.isna().any():
        raise ValueError("canonical distance columns contain nonnumeric/missing values")
    if area.isna().any() or area_log.isna().any():
        raise ValueError("canonical area columns contain nonnumeric/missing values")
    if (distance < 0).any():
        raise ValueError("distance_to_continent_km must be nonnegative")
    if (area < 0).any():
        raise ValueError("area_km2 must be nonnegative")

    expected_distance = np.log1p(distance.to_numpy(float))
    expected_area = np.log1p(area.to_numpy(float))
    if not np.allclose(distance_log.to_numpy(float), expected_distance, rtol=0.0, atol=atol):
        raise ValueError("historical log_distance_to_continent_km is not log1p(raw distance)")
    if not np.allclose(area_log.to_numpy(float), expected_area, rtol=0.0, atol=atol):
        raise ValueError("historical log_island_area_km2 is not log1p(raw area)")

    result = pd.DataFrame(
        {
            "island_id": work["island_id"],
            "log1p_distance_to_continent_km": distance_log,
            "log_area": area_log,
            "climate_pc1": pd.to_numeric(work["climate_pc1"], errors="coerce"),
            "climate_pc2": pd.to_numeric(work["climate_pc2"], errors="coerce"),
            "climate_pc3": pd.to_numeric(work["climate_pc3"], errors="coerce"),
            "climate_pc4": pd.to_numeric(work["climate_pc4"], errors="coerce"),
            "spatial_block": work["spatial_block"].fillna("").astype(str).str.strip(),
        },
        columns=N1_COVARIATE_COLUMNS,
    )
    if result["spatial_block"].eq("").any():
        raise ValueError("canonical Chapter 1 covariates contain blank spatial_block")
    return result


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _empty_gate(config: dict[str, Any], metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        "contract": config["contract"],
        "N1_pass": False,
        "status": metadata.get("status", "N1_not_evaluable_for_NEE_gate"),
        "confirmatory_channels": list(metadata.get("confirmatory_channels", [])),
        "n_confirmatory_channels": int(metadata.get("n_confirmatory_channels", 0)),
        "required_confirmatory_channels": int(
            metadata.get(
                "required_confirmatory_channels",
                config["channel_gate"]["min_confirmatory_channels_for_primary_N1"],
            )
        ),
        "failure_action": "stop_before_N2_and_keep_frozen_Chapter1",
        "claim_ceiling": "N1_not_promoted",
    }


def run_n1(
    projected: pd.DataFrame,
    qualification: pd.DataFrame,
    canonical_covariates_raw: pd.DataFrame,
    config: dict[str, Any],
    *,
    expected_islands: int | None = 8265,
) -> dict[str, Any]:
    """Run the frozen primary N1 model and required deletion robustness checks."""
    covariates = prepare_canonical_covariates(
        canonical_covariates_raw, expected_islands=expected_islands
    )
    frame, coefficients, global_test, result = fit_n1_primary(
        projected, qualification, covariates, config
    )
    if result.get("status") != "qualified":
        return {
            "frame": frame,
            "coefficients": pd.DataFrame(
                columns=["predictor", "estimate_log_odds", "cluster_robust_se"]
            ),
            "slopes": pd.DataFrame(
                columns=[
                    "channel_id",
                    "is_reference",
                    "isolation_slope_log_odds_per_sd",
                    "cluster_robust_se",
                ]
            ),
            "global_test": {},
            "model_comparison": {},
            "block_deletions": pd.DataFrame(columns=["deleted", "status", "p_value"]),
            "source_deletions": pd.DataFrame(columns=["deleted", "status", "p_value"]),
            "gate": _empty_gate(config, result),
            "metadata": result,
        }

    channels = list(result["confirmatory_channels"])
    reference = str(result["reference_channel"])
    block_deletions = deletion_robustness(
        frame, channels, reference, "spatial_block"
    )
    source_deletions = deletion_robustness(
        frame, channels, reference, "source_region_id"
    )
    comparison = dict(result["model_comparison"])
    gate = gate_receipt(
        global_test,
        comparison,
        block_deletions,
        source_deletions,
        len(channels),
        config,
    )
    return {
        "frame": frame,
        "coefficients": coefficients,
        "slopes": result["slopes"],
        "global_test": global_test,
        "model_comparison": comparison,
        "block_deletions": block_deletions,
        "source_deletions": source_deletions,
        "gate": gate,
        "metadata": {key: value for key, value in result.items() if key != "slopes"},
    }


def write_outputs(result: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    result["frame"].to_csv(output_dir / "N1_model_support.csv", index=False)
    result["coefficients"].to_csv(output_dir / "N1_coefficients.csv", index=False)
    result["slopes"].to_csv(output_dir / "N1_channel_slopes.csv", index=False)
    result["block_deletions"].to_csv(
        output_dir / "N1_leave_one_spatial_block_out.csv", index=False
    )
    result["source_deletions"].to_csv(
        output_dir / "N1_leave_one_source_region_out.csv", index=False
    )
    _write_json(output_dir / "N1_global_heterogeneity_test.json", result["global_test"])
    _write_json(output_dir / "N1_model_comparison.json", result["model_comparison"])
    _write_json(output_dir / "N1_gate_receipt.json", result["gate"])
    _write_json(output_dir / "N1_run_metadata.json", result["metadata"])


@app.command("run")
def run_command(
    projected_csv: Path = typer.Option(..., exists=True),
    qualification_csv: Path = typer.Option(..., exists=True),
    canonical_covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_nee_n1_model.yml"), exists=True),
    expected_islands: int = typer.Option(8265, min=1),
) -> None:
    config = load_config(config_path)
    result = run_n1(
        pd.read_csv(projected_csv, dtype=str).fillna(""),
        pd.read_csv(qualification_csv, dtype=str).fillna(""),
        pd.read_csv(canonical_covariates_csv),
        config,
        expected_islands=expected_islands,
    )
    write_outputs(result, output_dir)
    typer.echo(json.dumps(result["gate"], indent=2))


if __name__ == "__main__":
    app()
