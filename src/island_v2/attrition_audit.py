"""Phase 1 island coverage and attrition audit for v2.

The audit is deliberately upstream of outcome modelling. It joins realised flora
coverage, outcome-blind Bombus applicability, effort-aware Bombus observation
evidence, and reviewed trait evidence to show exactly which islands can enter
each domain-specific analysis.

No absent trait record is interpreted as a biological absence. Trait coverage
measures only whether accepted evidence exists for a taxon in a defined evidence
track.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_FLORA_COLUMNS = {"island_id", "analysis_included"}
REQUIRED_APPLICABILITY_COLUMNS = {"island_id", "applicability"}
REQUIRED_BOMBUS_COLUMNS = {"island_id", "bombus_occurrence_evidence"}
REQUIRED_ISLAND_SPECIES_COLUMNS = {"island_id", "species"}
REQUIRED_EVIDENCE_COLUMNS = {
    "accepted_species",
    "trait_name",
    "evidence_scope",
    "source_reliability",
    "review_status",
}


@app.callback()
def main() -> None:
    """Build trait-coverage and island-attrition outputs before v2 outcome models."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the versioned attrition-audit configuration."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {"primary_evidence_track", "trait_evidence_tracks", "domain_rules", "attrition_rules"}
    if not isinstance(config, dict) or not required.issubset(config):
        raise typer.BadParameter(
            "attrition config must contain primary_evidence_track, trait_evidence_tracks, "
            "domain_rules, and attrition_rules"
        )
    if config["primary_evidence_track"] not in config["trait_evidence_tracks"]:
        raise typer.BadParameter("primary_evidence_track is not defined in trait_evidence_tracks")
    return config


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _coerce_bool(series: pd.Series) -> pd.Series:
    return _text(series).str.lower().isin({"true", "1", "yes", "y"})


def _require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise typer.BadParameter(f"{label} missing columns: {', '.join(sorted(missing))}")


def _unique_by_island(table: pd.DataFrame, label: str) -> pd.DataFrame:
    if table["island_id"].duplicated().any():
        duplicates = sorted(_text(table.loc[table["island_id"].duplicated(), "island_id"]).unique())
        raise typer.BadParameter(f"{label} has duplicate island_id values: {duplicates[:5]}")
    return table.copy()


def _accepted_species_for_domain(
    evidence: pd.DataFrame, track: dict[str, Any], trait_names: list[str]
) -> set[str]:
    """Return species with accepted direct evidence in one track/domain."""
    accepted = evidence.loc[
        _text(evidence["review_status"]).str.lower().eq("accepted")
        & _text(evidence["evidence_scope"]).isin(set(track["accepted_scopes"]))
        & _text(evidence["source_reliability"]).isin(set(track["accepted_reliability"]))
        & _text(evidence["trait_name"]).isin(set(trait_names))
    ]
    return set(_text(accepted["accepted_species"]).loc[lambda values: values.ne("")])


def compute_trait_coverage(
    island_species: pd.DataFrame, accepted_evidence: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """Compute island-by-domain reviewed-trait coverage for each direct evidence track."""
    _require_columns(island_species, REQUIRED_ISLAND_SPECIES_COLUMNS, "island species table")
    _require_columns(accepted_evidence, REQUIRED_EVIDENCE_COLUMNS, "accepted evidence table")

    species_table = island_species[["island_id", "species"]].copy()
    species_table["island_id"] = _text(species_table["island_id"])
    species_table["species"] = _text(species_table["species"])
    species_table = species_table.loc[species_table["island_id"].ne("") & species_table["species"].ne("")]
    species_table = species_table.drop_duplicates(["island_id", "species"])

    evidence = accepted_evidence.copy()
    for column in REQUIRED_EVIDENCE_COLUMNS:
        evidence[column] = _text(evidence[column])

    rows: list[dict[str, Any]] = []
    for track_name, track in config["trait_evidence_tracks"].items():
        for domain, rule in config["domain_rules"].items():
            covered_species = _accepted_species_for_domain(evidence, track, list(rule["trait_names"]))
            for island_id, group in species_table.groupby("island_id", sort=True):
                total = int(group["species"].nunique())
                covered = int(group["species"].isin(covered_species).sum())
                fraction = covered / total if total else 0.0
                eligible = (
                    covered >= int(rule["min_covered_species"])
                    and fraction >= float(rule["min_coverage_fraction"])
                )
                rows.append(
                    {
                        "island_id": island_id,
                        "evidence_track": track_name,
                        "domain": domain,
                        "n_species_total": total,
                        "n_species_covered": covered,
                        "coverage_fraction": fraction,
                        "domain_eligible": eligible,
                    }
                )
    return pd.DataFrame(
        rows,
        columns=[
            "island_id",
            "evidence_track",
            "domain",
            "n_species_total",
            "n_species_covered",
            "coverage_fraction",
            "domain_eligible",
        ],
    )


def build_attrition_table(
    flora_diagnostics: pd.DataFrame,
    applicability: pd.DataFrame,
    bombus_diagnostics: pd.DataFrame,
    trait_coverage: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Join Phase 1 audit components into one island-level eligibility table."""
    _require_columns(flora_diagnostics, REQUIRED_FLORA_COLUMNS, "flora diagnostics")
    _require_columns(applicability, REQUIRED_APPLICABILITY_COLUMNS, "Bombus applicability table")
    _require_columns(bombus_diagnostics, REQUIRED_BOMBUS_COLUMNS, "Bombus diagnostics table")
    _unique_by_island(flora_diagnostics, "flora diagnostics")
    _unique_by_island(applicability, "Bombus applicability table")
    _unique_by_island(bombus_diagnostics, "Bombus diagnostics table")

    base = flora_diagnostics[["island_id", "analysis_included"]].copy()
    base["island_id"] = _text(base["island_id"])
    base["flora_eligible"] = _coerce_bool(base["analysis_included"])
    base = base.drop(columns=["analysis_included"])

    app_table = applicability[["island_id", "applicability"]].copy()
    app_table["island_id"] = _text(app_table["island_id"])
    app_table["applicability"] = _text(app_table["applicability"])
    base = base.merge(app_table, on="island_id", how="left", validate="one_to_one")
    base["applicability"] = base["applicability"].fillna("unresolved")
    base["applicability_missing"] = base["applicability"].eq("unresolved")

    bombus_columns = ["island_id", "bombus_occurrence_evidence"]
    if "bombus_channel_evaluable" in bombus_diagnostics.columns:
        bombus_columns.append("bombus_channel_evaluable")
    bombus_table = bombus_diagnostics[bombus_columns].copy()
    bombus_table["island_id"] = _text(bombus_table["island_id"])
    bombus_table["bombus_occurrence_evidence"] = _text(bombus_table["bombus_occurrence_evidence"])
    base = base.merge(bombus_table, on="island_id", how="left", validate="one_to_one")
    base["bombus_occurrence_evidence"] = base["bombus_occurrence_evidence"].fillna("unresolved")
    evaluable_states = set(config["attrition_rules"]["evaluable_bombus_evidence"])
    if "bombus_channel_evaluable" in base.columns:
        base["bombus_evaluable"] = _coerce_bool(base["bombus_channel_evaluable"])
    else:
        base["bombus_evaluable"] = base["bombus_occurrence_evidence"].isin(evaluable_states)
    base["bombus_diagnostics_missing"] = base["bombus_occurrence_evidence"].eq("unresolved")

    primary_track = config["primary_evidence_track"]
    primary = trait_coverage.loc[trait_coverage["evidence_track"].eq(primary_track)].copy()
    if primary.empty:
        coverage_pivot = pd.DataFrame(index=base["island_id"])
    else:
        primary["domain_eligible"] = _coerce_bool(primary["domain_eligible"])
        coverage_pivot = primary.pivot(index="island_id", columns="domain", values="domain_eligible")
    coverage_pivot = coverage_pivot.reindex(base["island_id"]).fillna(False)
    for domain in config["domain_rules"]:
        source = coverage_pivot[domain] if domain in coverage_pivot.columns else pd.Series(False, index=base.index)
        base[f"{domain}_eligible"] = source.to_numpy(dtype=bool)

    confirmatory_applicability = str(config["attrition_rules"]["confirmatory_applicability"])
    base["core_confirmatory_eligible"] = (
        base["flora_eligible"]
        & base["applicability"].eq(confirmatory_applicability)
        & base["bombus_evaluable"]
    )
    base["m1_eligible"] = (
        base["core_confirmatory_eligible"]
        & base["floral_signal_eligible"]
        & base["floral_architecture_eligible"]
    )
    base["m2_m3_eligible"] = base["m1_eligible"] & base["reproductive_assurance_eligible"]
    base["whole_flora_pollen_vector_eligible"] = base["flora_eligible"] & base["pollen_vector_eligible"]
    return base.sort_values("island_id").reset_index(drop=True)


def attrition_stage_summary(attrition: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Summarise sequential Phase 1 attrition for the predeclared primary track."""
    stages = [
        ("all_islands", pd.Series(True, index=attrition.index)),
        ("flora_eligible", attrition["flora_eligible"]),
        (
            "applicable_core_islands",
            attrition["flora_eligible"]
            & attrition["applicability"].eq(config["attrition_rules"]["confirmatory_applicability"]),
        ),
        ("bombus_observation_evaluable", attrition["core_confirmatory_eligible"]),
        (
            "floral_signal_eligible",
            attrition["core_confirmatory_eligible"] & attrition["floral_signal_eligible"],
        ),
        (
            "floral_architecture_eligible",
            attrition["core_confirmatory_eligible"]
            & attrition["floral_signal_eligible"]
            & attrition["floral_architecture_eligible"],
        ),
        ("m1_eligible", attrition["m1_eligible"]),
        ("reproductive_assurance_eligible", attrition["m2_m3_eligible"]),
        ("m2_m3_eligible", attrition["m2_m3_eligible"]),
    ]
    return pd.DataFrame(
        {
            "stage": [name for name, _ in stages],
            "n_islands": [int(mask.fillna(False).sum()) for _, mask in stages],
        }
    )


def continuation_summary(attrition: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    """Report predeclared feasibility status without inspecting trait outcomes."""
    thresholds = config.get("continuation_thresholds", {})
    n_m1 = int(attrition["m1_eligible"].sum())
    n_m2 = int(attrition["m2_m3_eligible"].sum())
    m1_pilot = int(thresholds.get("m1_pilot_min_islands", 30))
    m1_confirmatory = int(thresholds.get("m1_confirmatory_target_islands", 50))
    m2_exploratory = int(thresholds.get("m2_m3_exploratory_min_islands", 30))
    m2_confirmatory = int(thresholds.get("m2_m3_confirmatory_target_islands", 50))

    def status(value: int, exploratory: int, confirmatory: int) -> str:
        if value < exploratory:
            return "below_exploratory_minimum"
        if value < confirmatory:
            return "exploratory_or_pilot_only"
        return "confirmatory_target_met"

    return {
        "config_version": config.get("version"),
        "n_m1_eligible_islands": n_m1,
        "n_m2_m3_eligible_islands": n_m2,
        "m1_status": status(n_m1, m1_pilot, m1_confirmatory),
        "m2_m3_status": status(n_m2, m2_exploratory, m2_confirmatory),
        "thresholds": thresholds,
    }


@app.command("run")
def run(
    flora_diagnostics_csv: Path = typer.Option(..., exists=True, help="island_data_process_diagnostics.csv"),
    applicability_csv: Path = typer.Option(..., exists=True, help="Outcome-blind Bombus applicability table."),
    bombus_diagnostics_csv: Path = typer.Option(..., exists=True, help="Bombus observation diagnostics CSV."),
    island_species_csv: Path = typer.Option(..., exists=True, help="island_species_occurrences.csv"),
    accepted_evidence_csv: Path = typer.Option(..., exists=True, help="Human-reviewed accepted trait evidence CSV."),
    output_dir: Path = typer.Option(..., help="Directory for Phase 1 audit outputs."),
    config_path: Path = typer.Option(Path("config/attrition_audit.yml")),
) -> None:
    """Write reviewed trait coverage, island attrition, and feasibility summaries."""
    config = load_config(config_path)
    flora = pd.read_csv(flora_diagnostics_csv, dtype=str).fillna("")
    applicability = pd.read_csv(applicability_csv, dtype=str).fillna("")
    bombus = pd.read_csv(bombus_diagnostics_csv, dtype=str).fillna("")
    island_species = pd.read_csv(island_species_csv, dtype=str).fillna("")
    evidence = pd.read_csv(accepted_evidence_csv, dtype=str).fillna("")

    coverage = compute_trait_coverage(island_species, evidence, config)
    attrition = build_attrition_table(flora, applicability, bombus, coverage, config)
    stages = attrition_stage_summary(attrition, config)
    summary = continuation_summary(attrition, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(output_dir / "trait_coverage_by_island.csv", index=False)
    attrition.to_csv(output_dir / "island_analysis_attrition.csv", index=False)
    stages.to_csv(output_dir / "attrition_stage_summary.csv", index=False)
    (output_dir / "attrition_continuation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Audited {len(attrition)} islands: {summary['n_m1_eligible_islands']} M1-eligible, "
        f"{summary['n_m2_m3_eligible_islands']} M2/M3-eligible."
    )


if __name__ == "__main__":
    app()
