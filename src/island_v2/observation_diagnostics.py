"""Island observation-process diagnostics and realised-coverage inclusion.

Stage 3.5 audits *how* each island flora was observed and decides which islands
have enough reviewed **angiosperm** data to support a stable island-level floral
trait proportion. Raw Tracheophyta effort remains useful for observation-process
flags, but non-angiosperms and unreviewed taxa cannot inflate primary floral
analysis denominators.
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
REQUIRED_ANGIOSPERM_COVERAGE_COLUMNS = {
    "island_id",
    "n_accepted_angiosperm_species",
    "n_accepted_angiosperm_records",
    "accepted_angiosperm_species_fraction",
    "accepted_angiosperm_record_fraction",
}
RAW_SCREENING_SCOPE = "raw_tracheophyta_screening_only"
PRIMARY_SCOPE = "accepted_angiosperms"


@app.callback()
def main() -> None:
    """Audit island observation process and realised angiosperm coverage."""


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "flags" not in config:
        raise typer.BadParameter("diagnostics config must contain a top-level 'flags' mapping")
    return config


def _safe_fraction(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    denom = denominator.where(denominator > 0)
    return (numerator / denom).fillna(0.0)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _validate_unique_islands(table: pd.DataFrame, label: str) -> None:
    values = _text(table["island_id"])
    if values.eq("").any():
        raise typer.BadParameter(f"{label} contains blank island_id values")
    duplicate_ids = values.loc[values.duplicated()].unique().tolist()
    if duplicate_ids:
        raise typer.BadParameter(f"{label} contains duplicate island_id values: {duplicate_ids[:5]}")


def _prepare_angiosperm_coverage(coverage: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_ANGIOSPERM_COVERAGE_COLUMNS.difference(coverage.columns)
    if missing:
        raise typer.BadParameter(
            "angiosperm coverage table missing columns: " + ", ".join(sorted(missing))
        )
    _validate_unique_islands(coverage, "angiosperm coverage table")
    ordered = ["island_id", *sorted(REQUIRED_ANGIOSPERM_COVERAGE_COLUMNS - {"island_id"})]
    result = coverage[ordered].copy()
    result["island_id"] = _text(result["island_id"])
    for column in REQUIRED_ANGIOSPERM_COVERAGE_COLUMNS - {"island_id"}:
        result[column] = pd.to_numeric(result[column], errors="coerce").fillna(0)
    return result


def compute_island_diagnostics(
    effort: pd.DataFrame,
    config: dict[str, Any],
    angiosperm_coverage: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return data-process flags and a taxonomically explicit inclusion gate.

    With ``angiosperm_coverage``, inclusion uses accepted angiosperm species and
    record counts. Without it, the raw coverage diagnostic remains available but
    is labelled clearly as a non-primary Tracheophyta screen.
    """
    missing = REQUIRED_EFFORT_COLUMNS.difference(effort.columns)
    if missing:
        raise typer.BadParameter(f"effort table missing columns: {', '.join(sorted(missing))}")
    if effort.empty:
        return pd.DataFrame(
            columns=[
                *effort.columns,
                "specimen_fraction",
                "cultivated_fraction",
                "data_process_flags",
                "taxonomic_scope_for_inclusion",
                "analysis_included",
            ]
        )

    _validate_unique_islands(effort, "effort table")
    flags_config = config["flags"]
    min_species = int(config["min_species_for_proportions"])
    min_records = int(config["min_records"])

    out = effort.copy()
    out["island_id"] = _text(out["island_id"])
    for column in (
        "n_records",
        "n_species",
        "n_datasets",
        "n_preserved_specimen",
        "n_establishment_cultivated",
        "year_max",
    ):
        out[column] = pd.to_numeric(out[column], errors="coerce")
    uncertainty = pd.to_numeric(out["median_coordinate_uncertainty_m"], errors="coerce")
    out["specimen_fraction"] = _safe_fraction(out["n_preserved_specimen"], out["n_records"])
    out["cultivated_fraction"] = _safe_fraction(out["n_establishment_cultivated"], out["n_records"])

    if angiosperm_coverage is None:
        gate_species = out["n_species"]
        gate_records = out["n_records"]
        scope_name = RAW_SCREENING_SCOPE
    else:
        coverage = _prepare_angiosperm_coverage(angiosperm_coverage)
        out = out.merge(coverage, on="island_id", how="left", validate="one_to_one")
        for column in REQUIRED_ANGIOSPERM_COVERAGE_COLUMNS - {"island_id"}:
            out[column] = pd.to_numeric(out[column], errors="coerce").fillna(0)
        gate_species = out["n_accepted_angiosperm_species"]
        gate_records = out["n_accepted_angiosperm_records"]
        scope_name = PRIMARY_SCOPE

    flag_columns = {
        "single_dataset": out["n_datasets"] <= int(flags_config["single_dataset_max"]),
        "low_specimen_support": out["specimen_fraction"] < float(flags_config["low_specimen_support_fraction"]),
        "high_cultivation": out["cultivated_fraction"] > float(flags_config["high_cultivation_fraction"]),
        "coarse_coordinates": uncertainty > float(flags_config["coarse_coordinate_uncertainty_m"]),
        "stale_records": out["year_max"] < int(flags_config["stale_last_record_year"]),
        "sparse_flora": gate_species < min_species,
    }
    out["data_process_flags"] = [
        "|".join(name for name, series in flag_columns.items() if bool(series.iloc[index]))
        for index in range(len(out))
    ]
    out["taxonomic_scope_for_inclusion"] = scope_name
    out["analysis_included"] = (gate_species >= min_species) & (gate_records >= min_records)
    return out.reset_index(drop=True)


def coverage_sensitivity(
    diagnostics: pd.DataFrame,
    config: dict[str, Any],
    taxonomic_scope: str = RAW_SCREENING_SCOPE,
) -> pd.DataFrame:
    """Count islands retained at predeclared thresholds in one declared scope."""
    min_records = int(config["min_records"])
    if taxonomic_scope == PRIMARY_SCOPE:
        species = pd.to_numeric(
            diagnostics.get("n_accepted_angiosperm_species", pd.Series(dtype=float)), errors="coerce"
        )
        records = pd.to_numeric(
            diagnostics.get("n_accepted_angiosperm_records", pd.Series(dtype=float)), errors="coerce"
        )
    else:
        species = pd.to_numeric(diagnostics.get("n_species", pd.Series(dtype=float)), errors="coerce")
        records = pd.to_numeric(diagnostics.get("n_records", pd.Series(dtype=float)), errors="coerce")
    rows = []
    for threshold in config.get("species_threshold_sensitivity", [config["min_species_for_proportions"]]):
        included = (species >= int(threshold)) & (records >= min_records)
        rows.append(
            {
                "taxonomic_scope": taxonomic_scope,
                "min_species": int(threshold),
                "n_islands_included": int(included.sum()),
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    effort_csv: Path = typer.Option(..., exists=True, help="Collector island_observation_effort.csv."),
    output_dir: Path = typer.Option(..., help="Directory for diagnostics outputs."),
    angiosperm_coverage_csv: Path | None = typer.Option(
        None,
        help="Reviewed island_angiosperm_coverage.csv; required for confirmatory floral analyses.",
    ),
    config_path: Path = typer.Option(Path("config/observation_diagnostics.yml")),
) -> None:
    """Write data-process diagnostics and the taxonomically explicit coverage gate."""
    config = load_config(config_path)
    effort = pd.read_csv(effort_csv, dtype=str).fillna("")
    coverage = (
        pd.read_csv(angiosperm_coverage_csv, dtype=str).fillna("")
        if angiosperm_coverage_csv is not None
        else None
    )
    diagnostics = compute_island_diagnostics(effort, config, coverage)
    scope = PRIMARY_SCOPE if coverage is not None else RAW_SCREENING_SCOPE
    sensitivity = coverage_sensitivity(diagnostics, config, scope)

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
        "taxonomic_scope_for_inclusion": scope,
        "min_species_for_proportions": int(config["min_species_for_proportions"]),
        "min_records": int(config["min_records"]),
        "n_islands_with_records": n_islands,
        "n_islands_analysis_included": n_included,
        "data_process_flag_counts": flag_counts,
    }
    (output_dir / "diagnostics_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    typer.echo(
        f"Diagnosed {n_islands} islands; {n_included} pass the {scope} coverage gate "
        f"(>= {summary['min_species_for_proportions']} species and >= {summary['min_records']} records)."
    )


if __name__ == "__main__":
    app()
