"""Build an outcome-blind island-level Bombus applicability registry.

The Bombus channel is a regional biological hypothesis, not a global absence
variable. This module makes the Phase 0 classification reproducible by joining:

1. a reviewed source-region registry;
2. reviewed island-to-source-region assignments; and
3. source-region evidence records.

The inputs intentionally exclude island climate, occurrence/non-detection,
floral traits, reproductive traits, and model results. Those are downstream
predictors or outcomes and must not be used to decide whether the Bombus
hypothesis is applicable.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_MANIFEST_COLUMNS = {"island_id"}
REQUIRED_REGION_COLUMNS = {
    "source_region_id",
    "biogeographic_region",
    "floristic_source_pool",
    "plausible_dispersal_connection_or_gateway",
    "native_bombus_status",
    "applicability",
    "applicability_rationale",
    "source_region_evidence_id",
    "review_status",
    "regional_analysis_only",
    "reviewer",
    "decision_date",
}
REQUIRED_ASSIGNMENT_COLUMNS = {
    "island_id",
    "source_region_id",
    "assignment_rationale",
    "assignment_evidence_id",
    "review_status",
    "reviewer",
    "decision_date",
}
REQUIRED_EVIDENCE_COLUMNS = {
    "source_region_evidence_id",
    "source_region_id",
    "evidence_type",
    "source_citation",
    "source_url",
    "evidence_excerpt",
    "review_status",
    "reviewer",
    "decision_date",
}
VALID_APPLICABILITY = {"applicable", "structurally_not_applicable", "unresolved"}
VALID_NATIVE_STATUS = {"native_present", "native_absent", "uncertain"}
PROHIBITED_INPUT_COLUMNS = {
    "island_climate",
    "island_environmental_compatibility_score",
    "island_Bombus_occurrence_or_non_detection",
    "floral_traits",
    "reproductive_traits",
    "model_results",
}
OUTPUT_COLUMNS = [
    "island_id",
    "source_region_id",
    "biogeographic_region",
    "floristic_source_pool",
    "plausible_dispersal_connection_or_gateway",
    "native_bombus_status",
    "applicability",
    "applicability_rationale",
    "source_region_evidence_id",
    "regional_analysis_only",
    "reviewer",
    "decision_date",
    "source_region_review_status",
    "assignment_review_status",
    "source_region_assignment_missing",
]


@app.callback()
def main() -> None:
    """Build a frozen Bombus applicability registry before outcome analysis."""


def load_config(path: Path) -> dict[str, Any]:
    """Load and validate the governing applicability policy."""
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or "categories" not in config:
        raise typer.BadParameter("Bombus applicability config must contain a categories mapping")
    policy_categories = set(config["categories"])
    if policy_categories != VALID_APPLICABILITY:
        raise typer.BadParameter(
            "Bombus applicability config categories must be exactly "
            f"{sorted(VALID_APPLICABILITY)}"
        )
    return config


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _as_bool(series: pd.Series) -> pd.Series:
    values = _text(series).str.lower()
    invalid = values.loc[~values.isin({"true", "false", "1", "0", "yes", "no", "y", "n", ""})]
    if not invalid.empty:
        raise typer.BadParameter(
            "regional_analysis_only must contain only true/false values; found "
            f"{sorted(invalid.unique())[:5]}"
        )
    return values.isin({"true", "1", "yes", "y"})


def _require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise typer.BadParameter(f"{label} missing columns: {', '.join(sorted(missing))}")
    prohibited = PROHIBITED_INPUT_COLUMNS.intersection(table.columns)
    if prohibited:
        raise typer.BadParameter(
            f"{label} contains prohibited outcome/predictor columns: {', '.join(sorted(prohibited))}"
        )


def _require_unique(table: pd.DataFrame, column: str, label: str) -> None:
    values = _text(table[column])
    if values.eq("").any():
        raise typer.BadParameter(f"{label} has blank {column} values")
    duplicates = values.loc[values.duplicated()].unique().tolist()
    if duplicates:
        raise typer.BadParameter(f"{label} has duplicate {column} values: {duplicates[:5]}")


def _reviewed(table: pd.DataFrame, column: str = "review_status") -> pd.Series:
    return _text(table[column]).str.lower().eq("accepted")


def validate_region_registry(regions: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Validate reviewed source-region entries and their declared evidence links."""
    _require_columns(regions, REQUIRED_REGION_COLUMNS, "source-region registry")
    _require_columns(evidence, REQUIRED_EVIDENCE_COLUMNS, "source-region evidence")
    _require_unique(regions, "source_region_id", "source-region registry")
    _require_unique(evidence, "source_region_evidence_id", "source-region evidence")

    result = regions.copy()
    for column in REQUIRED_REGION_COLUMNS:
        result[column] = _text(result[column])
    result["regional_analysis_only"] = _as_bool(result["regional_analysis_only"])

    invalid_applicability = set(result["applicability"]).difference(VALID_APPLICABILITY)
    if invalid_applicability:
        raise typer.BadParameter(f"Invalid applicability values: {sorted(invalid_applicability)}")
    invalid_native = set(result["native_bombus_status"]).difference(VALID_NATIVE_STATUS)
    if invalid_native:
        raise typer.BadParameter(f"Invalid native_bombus_status values: {sorted(invalid_native)}")

    inconsistent_present = result.loc[
        result["applicability"].eq("applicable")
        & ~result["native_bombus_status"].eq("native_present"),
        "source_region_id",
    ].tolist()
    if inconsistent_present:
        raise typer.BadParameter(
            "applicable source regions require native_bombus_status=native_present: "
            f"{inconsistent_present[:5]}"
        )
    inconsistent_absent = result.loc[
        result["applicability"].eq("structurally_not_applicable")
        & ~result["native_bombus_status"].eq("native_absent"),
        "source_region_id",
    ].tolist()
    if inconsistent_absent:
        raise typer.BadParameter(
            "structurally_not_applicable regions require native_bombus_status=native_absent: "
            f"{inconsistent_absent[:5]}"
        )

    evidence_table = evidence.copy()
    for column in REQUIRED_EVIDENCE_COLUMNS:
        evidence_table[column] = _text(evidence_table[column])
    evidence_lookup = evidence_table.set_index("source_region_evidence_id", drop=False)
    for row in result.itertuples(index=False):
        evidence_id = str(row.source_region_evidence_id)
        if evidence_id not in evidence_lookup.index:
            raise typer.BadParameter(
                f"source-region registry references missing source_region_evidence_id={evidence_id!r}"
            )
        evidence_row = evidence_lookup.loc[evidence_id]
        if str(evidence_row["source_region_id"]) != str(row.source_region_id):
            raise typer.BadParameter(
                "source-region evidence must reference the same source_region_id: "
                f"{evidence_id!r}"
            )
    return result


def build_applicability_table(
    island_manifest: pd.DataFrame,
    assignments: pd.DataFrame,
    regions: pd.DataFrame,
    evidence: pd.DataFrame,
) -> pd.DataFrame:
    """Build one classification row for every island in the frozen universe.

    Missing or not-yet-reviewed source-region assignments remain explicitly
    `unresolved`; they are never silently dropped and never classified from
    island occurrence or trait data.
    """
    _require_columns(island_manifest, REQUIRED_MANIFEST_COLUMNS, "island manifest")
    _require_columns(assignments, REQUIRED_ASSIGNMENT_COLUMNS, "island source-region assignments")
    _require_columns(regions, REQUIRED_REGION_COLUMNS, "source-region registry")
    _require_columns(evidence, REQUIRED_EVIDENCE_COLUMNS, "source-region evidence")
    _require_unique(island_manifest, "island_id", "island manifest")
    _require_unique(assignments, "island_id", "island source-region assignments")
    validated_regions = validate_region_registry(regions, evidence)

    islands = island_manifest[["island_id"]].copy()
    islands["island_id"] = _text(islands["island_id"])
    assignment_table = assignments.copy()
    for column in REQUIRED_ASSIGNMENT_COLUMNS:
        assignment_table[column] = _text(assignment_table[column])
    unknown_islands = sorted(set(assignment_table["island_id"]).difference(set(islands["island_id"])))
    if unknown_islands:
        raise typer.BadParameter(
            "island source-region assignments contain island IDs absent from the frozen manifest: "
            f"{unknown_islands[:5]}"
        )

    merged = islands.merge(
        assignment_table,
        on="island_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_assignment"),
    ).merge(
        validated_regions,
        on="source_region_id",
        how="left",
        validate="many_to_one",
        suffixes=("_assignment", ""),
    )

    rows: list[dict[str, Any]] = []
    for row in merged.to_dict("records"):
        island_id = str(row["island_id"])
        missing_assignment = not str(row.get("source_region_id", "")).strip()
        assignment_accepted = str(row.get("review_status_assignment", "")).strip().lower() == "accepted"
        region_accepted = str(row.get("review_status", "")).strip().lower() == "accepted"
        source_region_exists = bool(str(row.get("biogeographic_region", "")).strip())

        if missing_assignment:
            applicability = "unresolved"
            rationale = "No source-region assignment has been recorded for this island."
        elif not assignment_accepted:
            applicability = "unresolved"
            rationale = "The source-region assignment is not yet accepted for analysis."
        elif not source_region_exists:
            applicability = "unresolved"
            rationale = "The assigned source_region_id is absent from the reviewed source-region registry."
        elif not region_accepted:
            applicability = "unresolved"
            rationale = "The assigned source-region registry entry is not yet accepted for analysis."
        else:
            applicability = str(row["applicability"])
            rationale = str(row["applicability_rationale"])
            assignment_rationale = str(row.get("assignment_rationale", "")).strip()
            if assignment_rationale:
                rationale = f"{rationale} Island assignment: {assignment_rationale}"

        rows.append(
            {
                "island_id": island_id,
                "source_region_id": str(row.get("source_region_id", "")),
                "biogeographic_region": str(row.get("biogeographic_region", "")),
                "floristic_source_pool": str(row.get("floristic_source_pool", "")),
                "plausible_dispersal_connection_or_gateway": str(
                    row.get("plausible_dispersal_connection_or_gateway", "")
                ),
                "native_bombus_status": str(row.get("native_bombus_status", "")),
                "applicability": applicability,
                "applicability_rationale": rationale,
                "source_region_evidence_id": str(row.get("source_region_evidence_id", "")),
                "regional_analysis_only": bool(row.get("regional_analysis_only", False)),
                "reviewer": str(row.get("reviewer_assignment", "")) or str(row.get("reviewer", "")),
                "decision_date": str(row.get("decision_date_assignment", "")) or str(row.get("decision_date", "")),
                "source_region_review_status": str(row.get("review_status", "")),
                "assignment_review_status": str(row.get("review_status_assignment", "")),
                "source_region_assignment_missing": missing_assignment,
            }
        )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS).sort_values("island_id").reset_index(drop=True)


def applicability_summary(table: pd.DataFrame) -> dict[str, Any]:
    """Summarise coverage and the predeclared Core/contingency/out-of-domain strata."""
    category_counts = table["applicability"].value_counts().sort_index().to_dict()
    core = table["applicability"].eq("applicable") & ~table["regional_analysis_only"]
    contingency = table["applicability"].eq("applicable") & table["regional_analysis_only"]
    return {
        "n_islands": int(len(table)),
        "applicability_counts": {str(key): int(value) for key, value in category_counts.items()},
        "n_core_confirmatory_islands": int(core.sum()),
        "n_regional_contingency_islands": int(contingency.sum()),
        "n_out_of_domain_islands": int(table["applicability"].eq("structurally_not_applicable").sum()),
        "n_unresolved_islands": int(table["applicability"].eq("unresolved").sum()),
        "n_missing_source_region_assignments": int(table["source_region_assignment_missing"].sum()),
    }


@app.command("build")
def build(
    island_manifest_csv: Path = typer.Option(..., exists=True, help="Frozen island manifest with island_id."),
    island_source_assignment_csv: Path = typer.Option(..., exists=True, help="Reviewed island-to-source-region table."),
    source_region_registry_csv: Path = typer.Option(..., exists=True, help="Reviewed source-region Bombus registry."),
    source_region_evidence_csv: Path = typer.Option(..., exists=True, help="Source-region evidence table."),
    output_dir: Path = typer.Option(..., help="Directory for applicability outputs."),
    config_path: Path = typer.Option(Path("config/bombus_applicability.yml")),
) -> None:
    """Write the frozen island-level Bombus applicability registry and strata."""
    load_config(config_path)
    island_manifest = pd.read_csv(island_manifest_csv, dtype=str).fillna("")
    assignments = pd.read_csv(island_source_assignment_csv, dtype=str).fillna("")
    regions = pd.read_csv(source_region_registry_csv, dtype=str).fillna("")
    evidence = pd.read_csv(source_region_evidence_csv, dtype=str).fillna("")

    table = build_applicability_table(island_manifest, assignments, regions, evidence)
    summary = applicability_summary(table)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "bombus_applicability_by_island.csv", index=False)
    table.loc[
        table["applicability"].eq("applicable") & ~table["regional_analysis_only"]
    ].to_csv(output_dir / "core_confirmatory_islands.csv", index=False)
    table.loc[
        table["applicability"].eq("applicable") & table["regional_analysis_only"]
    ].to_csv(output_dir / "regional_contingency_islands.csv", index=False)
    table.loc[table["applicability"].eq("structurally_not_applicable")].to_csv(
        output_dir / "out_of_domain_islands.csv", index=False
    )
    table.loc[table["applicability"].eq("unresolved")].to_csv(
        output_dir / "unresolved_applicability_islands.csv", index=False
    )
    (output_dir / "bombus_applicability_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Built applicability registry for {summary['n_islands']} islands: "
        f"{summary['n_core_confirmatory_islands']} Core, "
        f"{summary['n_regional_contingency_islands']} contingency, "
        f"{summary['n_out_of_domain_islands']} out-of-domain, "
        f"{summary['n_unresolved_islands']} unresolved."
    )


if __name__ == "__main__":
    app()
