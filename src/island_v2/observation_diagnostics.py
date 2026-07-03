"""Island observation-process diagnostics and realised-coverage inclusion.

Stage 3.5 of the v2 pipeline: before any bumblebee predictor or trait model,
audit *how* each island flora was observed and decide which islands carry enough
realised data to support a stable island-level trait proportion. This implements
the v2 rule that analysis inclusion follows realised coverage, not island area
(area stays an analytical covariate), and that observation-process concerns
(single-source floras, weak voucher support, cultivation, coarse coordinates,
stale records) are stated as explicit audit flags rather than silent filters.

Input is the collector's `island_observation_effort.csv`; the pure functions
here are unit-testable without GBIF.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_EFFORT_COLUMNS = {
    "island_id",
    "n_records",
    "n_species",
    "n_datasets",
    "n_preserved_specimen",
    "n_establishment_cultivated",
    "median_coordinate_uncertainty_m",
    "year_max",
}


@app.callback()
def main() -> None:
    """Audit island observation process and realised-coverage inclusion."""


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "flags" not in config:
        raise typer.BadParameter("diagnostics config must contain a top-level 'flags' mapping")
    return config


def _safe_fraction(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    denom = denominator.where(denominator > 0)
    return (numerator / denom).fillna(0.0)


def compute_island_diagnostics(effort: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Per-island observation-process metrics, concern flags, and inclusion.

    `effort` is the collector's island_observation_effort table. Returns one row
    per island with derived fractions, a pipe-joined `data_process_flags`
    string, and a boolean `analysis_included` from the realised-coverage gate.
    """
    missing = REQUIRED_EFFORT_COLUMNS.difference(effort.columns)
    if missing:
        raise typer.BadParameter(f"effort table missing columns: {', '.join(sorted(missing))}")
    if effort.empty:
        return pd.DataFrame(columns=[*effort.columns, "specimen_fraction", "cultivated_fraction",
                                     "data_process_flags", "analysis_included"])

    flags_config = config["flags"]
    min_species = int(config["min_species_for_proportions"])
    min_records = int(config["min_records"])

    out = effort.copy()
    for column in ("n_records", "n_species", "n_datasets", "n_preserved_specimen",
                   "n_establishment_cultivated", "year_max"):
        out[column] = pd.to_numeric(out[column], errors="coerce")
    uncertainty = pd.to_numeric(out["median_coordinate_uncertainty_m"], errors="coerce")

    out["specimen_fraction"] = _safe_fraction(out["n_preserved_specimen"], out["n_records"])
    out["cultivated_fraction"] = _safe_fraction(out["n_establishment_cultivated"], out["n_records"])

    single_dataset = out["n_datasets"] <= int(flags_config["single_dataset_max"])
    low_specimen = out["specimen_fraction"] < float(flags_config["low_specimen_support_fraction"])
    high_cultivation = out["cultivated_fraction"] > float(flags_config["high_cultivation_fraction"])
    coarse_coords = uncertainty > float(flags_config["coarse_coordinate_uncertainty_m"])
    stale = out["year_max"] < int(flags_config["stale_last_record_year"])
    sparse = out["n_species"] < min_species

    flag_columns = {
        "single_dataset": single_dataset,
        "low_specimen_support": low_specimen,
        "high_cultivation": high_cultivation,
        "coarse_coordinates": coarse_coords,
        "stale_records": stale,
        "sparse_flora": sparse,
    }
    out["data_process_flags"] = [
        "|".join(name for name, series in flag_columns.items() if bool(series.iloc[i]))
        for i in range(len(out))
    ]
    out["analysis_included"] = (out["n_species"] >= min_species) & (out["n_records"] >= min_records)
    return out.reset_index(drop=True)


def coverage_sensitivity(effort: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """How many islands stay included at each species threshold in the config."""
    min_records = int(config["min_records"])
    species = pd.to_numeric(effort.get("n_species", pd.Series(dtype=float)), errors="coerce")
    records = pd.to_numeric(effort.get("n_records", pd.Series(dtype=float)), errors="coerce")
    rows = []
    for threshold in config.get("species_threshold_sensitivity", [config["min_species_for_proportions"]]):
        included = (species >= int(threshold)) & (records >= min_records)
        rows.append({"min_species": int(threshold), "n_islands_included": int(included.sum())})
    return pd.DataFrame(rows)


@app.command("run")
def run(
    effort_csv: Path = typer.Option(..., exists=True, help="Collector island_observation_effort.csv."),
    output_dir: Path = typer.Option(..., help="Directory for diagnostics outputs."),
    config_path: Path = typer.Option(Path("config/observation_diagnostics.yml")),
) -> None:
    """Write per-island diagnostics, a coverage-sensitivity table, and a summary."""
    config = load_config(config_path)
    effort = pd.read_csv(effort_csv, dtype=str).fillna("")
    diagnostics = compute_island_diagnostics(effort, config)
    sensitivity = coverage_sensitivity(effort, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics.to_csv(output_dir / "island_data_process_diagnostics.csv", index=False)
    sensitivity.to_csv(output_dir / "coverage_inclusion_sensitivity.csv", index=False)

    n_islands = int(len(diagnostics))
    n_included = int(diagnostics["analysis_included"].sum()) if n_islands else 0
    flag_counts: dict[str, int] = {}
    for value in diagnostics.get("data_process_flags", pd.Series(dtype=str)):
        for flag in str(value).split("|"):
            if flag:
                flag_counts[flag] = flag_counts.get(flag, 0) + 1
    summary = {
        "config_version": config.get("version"),
        "min_species_for_proportions": int(config["min_species_for_proportions"]),
        "min_records": int(config["min_records"]),
        "n_islands_with_records": n_islands,
        "n_islands_analysis_included": n_included,
        "data_process_flag_counts": flag_counts,
    }
    (output_dir / "diagnostics_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(
        f"Diagnosed {n_islands} islands; {n_included} pass the realised-coverage inclusion gate "
        f"(>= {summary['min_species_for_proportions']} species and >= {summary['min_records']} records)."
    )


if __name__ == "__main__":
    app()
