"""Aggregate species-level Bombus niche and regime outputs to island analysis inputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def aggregate_islands(scores: pd.DataFrame) -> pd.DataFrame:
    required = {
        "island_id",
        "bombus_species",
        "environmental_compatibility",
        "environmental_extrapolation",
        "bombus_regime",
        "primary_bombus_channel_analysis",
        "global_comparison_only",
    }
    missing = required.difference(scores.columns)
    if missing:
        raise typer.BadParameter(f"regime score table missing columns: {sorted(missing)}")

    work = scores.copy()
    work["environmental_compatibility"] = pd.to_numeric(
        work["environmental_compatibility"], errors="coerce"
    )
    work["environmental_extrapolation"] = (
        work["environmental_extrapolation"]
        .fillna(False)
        .astype(str)
        .str.lower()
        .isin({"true", "1", "yes"})
    )

    rows: list[dict[str, object]] = []
    for island_id, group in work.groupby("island_id", sort=True):
        scored = group["environmental_compatibility"].dropna().astype(float)
        regimes = sorted(set(group["bombus_regime"].dropna().astype(str)))
        rows.append(
            {
                "island_id": island_id,
                "n_bombus_species_total": int(group["bombus_species"].nunique()),
                "n_bombus_species_scored": int(scored.size),
                "n_bombus_species_unscored": int(group["bombus_species"].nunique() - scored.size),
                "n_environmentally_compatible_species": int((scored >= 0.5).sum()),
                "share_environmentally_compatible_species": (
                    float((scored >= 0.5).mean()) if not scored.empty else pd.NA
                ),
                "environmental_compatibility_max": float(scored.max()) if not scored.empty else pd.NA,
                "environmental_compatibility_mean": float(scored.mean()) if not scored.empty else pd.NA,
                "environmental_extrapolation_share": float(group["environmental_extrapolation"].mean()),
                "bombus_regime_set": "|".join(regimes),
                "primary_bombus_channel_analysis": bool(
                    group["primary_bombus_channel_analysis"].fillna(False).astype(bool).any()
                ),
                "global_comparison_only": bool(
                    group["global_comparison_only"].fillna(False).astype(bool).all()
                ),
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    regime_scores_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    scores = pd.read_csv(regime_scores_csv)
    result = aggregate_islands(scores)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "bombus_island_analysis_input.csv", index=False)
    summary = {
        "contract": "bombus_island_analysis_input_v1",
        "n_islands": int(result["island_id"].nunique()),
        "n_primary_analysis_islands": int(result["primary_bombus_channel_analysis"].sum()),
        "n_global_comparison_only_islands": int(result["global_comparison_only"].sum()),
        "interpretation": (
            "Island-level summaries preserve the distinction between environmental compatibility, "
            "regional applicability, and extrapolation."
        ),
    }
    (output_dir / "bombus_island_analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
