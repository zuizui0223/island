"""Join Bombus environmental compatibility to the frozen applicability regime."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

VALID_APPLICABILITY = {"applicable", "structurally_not_applicable", "unresolved"}


def classify_regime(applicability: str, compatibility: float | None) -> str:
    """Classify interpretation without using compatibility to define applicability."""
    if applicability == "structurally_not_applicable":
        return "structural_bombus_absence_comparison"
    if applicability == "unresolved":
        return "applicability_unresolved"
    if compatibility is None or pd.isna(compatibility):
        return "applicable_environment_unresolved"
    if float(compatibility) >= 0.5:
        return "applicable_environment_compatible"
    return "applicable_environment_mismatch"


@app.command("run")
def run(
    compatibility_csv: Path = typer.Option(..., exists=True),
    applicability_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    """Attach the predeclared regional applicability regime to island×Bombus scores."""
    scores = pd.read_csv(compatibility_csv)
    applicability = pd.read_csv(applicability_csv, dtype=str).fillna("")
    required_scores = {"island_id", "environmental_compatibility"}
    required_applicability = {"island_id", "applicability", "native_bombus_status"}
    if missing := required_scores.difference(scores.columns):
        raise typer.BadParameter(f"compatibility table missing columns: {sorted(missing)}")
    if missing := required_applicability.difference(applicability.columns):
        raise typer.BadParameter(f"applicability table missing columns: {sorted(missing)}")
    invalid = set(applicability["applicability"]).difference(VALID_APPLICABILITY)
    if invalid:
        raise typer.BadParameter(f"invalid applicability values: {sorted(invalid)}")

    registry = applicability.drop_duplicates("island_id")
    merged = scores.merge(registry, on="island_id", how="left", validate="many_to_one")
    merged["applicability"] = merged["applicability"].replace("", pd.NA).fillna("unresolved")
    merged["bombus_regime"] = [
        classify_regime(app, comp)
        for app, comp in zip(
            merged["applicability"], merged["environmental_compatibility"], strict=True
        )
    ]
    merged["primary_bombus_channel_analysis"] = merged["applicability"].eq("applicable")
    merged["global_comparison_only"] = ~merged["primary_bombus_channel_analysis"]

    output_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_dir / "bombus_environmental_compatibility_with_regime.csv.gz", index=False, compression="gzip")
    summary = {
        "contract": "bombus_regime_join_v1",
        "n_rows": int(len(merged)),
        "n_islands": int(merged["island_id"].nunique()),
        "regime_counts": {str(k): int(v) for k, v in merged["bombus_regime"].value_counts().items()},
        "interpretation": "Applicability is frozen upstream and is never inferred from climate, occurrence, or plant traits.",
    }
    (output_dir / "bombus_regime_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
