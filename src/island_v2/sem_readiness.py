"""Audit whether the preregistered M0-M4 SEM ladder is estimable under current gates."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def audit_readiness(
    island_metrics: pd.DataFrame,
    applicability: pd.DataFrame,
    minimum_complete_islands: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    required_metrics = {"island_id", "environmental_compatibility_max"}
    required_app = {"island_id", "applicability", "regional_analysis_only"}
    missing_metrics = required_metrics.difference(island_metrics.columns)
    missing_app = required_app.difference(applicability.columns)
    if missing_metrics:
        raise typer.BadParameter(f"island metrics missing columns: {sorted(missing_metrics)}")
    if missing_app:
        raise typer.BadParameter(f"applicability table missing columns: {sorted(missing_app)}")

    metrics = island_metrics.drop_duplicates("island_id").copy()
    app_table = applicability.drop_duplicates("island_id").copy()
    data = metrics.merge(app_table, on="island_id", how="left", validate="one_to_one")
    data["applicability"] = data["applicability"].fillna("unresolved").astype(str)
    regional = data["regional_analysis_only"].fillna(False)
    if not pd.api.types.is_bool_dtype(regional):
        regional = regional.astype(str).str.lower().isin({"true", "1", "yes", "y"})
    compatibility = pd.to_numeric(data["environmental_compatibility_max"], errors="coerce")

    resolved = data["applicability"].isin({"applicable", "structurally_not_applicable"})
    primary = data["applicability"].eq("applicable") & ~regional
    weak_primary = primary & compatibility.lt(0.5)
    structural = data["applicability"].eq("structurally_not_applicable")
    m4_scope = structural | weak_primary

    scope_masks = {
        "M0": resolved,
        "M1": primary,
        "M2": primary,
        "M3": primary,
        "M4": m4_scope,
    }
    rows = []
    for model, mask in scope_masks.items():
        n = int(data.loc[mask, "island_id"].nunique())
        rows.append(
            {
                "model": model,
                "n_islands_scope": n,
                "minimum_complete_islands": int(minimum_complete_islands),
                "ready": bool(n >= minimum_complete_islands),
            }
        )
    readiness = pd.DataFrame(rows)
    summary = {
        "contract": "m0_m4_sem_readiness_v1",
        "n_islands_total": int(data["island_id"].nunique()),
        "n_applicable_core": int(primary.sum()),
        "n_structurally_not_applicable": int(structural.sum()),
        "n_unresolved": int(data["applicability"].eq("unresolved").sum()),
        "n_weak_primary_bombus_channel": int(weak_primary.sum()),
        "minimum_complete_islands": int(minimum_complete_islands),
        "model_ready": {
            row.model: bool(row.ready) for row in readiness.itertuples(index=False)
        },
        "all_models_ready": bool(readiness["ready"].all()),
        "interpretation": (
            "This is a design/readiness gate, not a model result. Models below the preregistered "
            "minimum island count must not be reported as fitted evidence."
        ),
    }
    return readiness, summary


@app.command("run")
def run(
    island_metrics_csv: Path = typer.Option(..., exists=True),
    applicability_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    sem_config: Path = typer.Option(Path("config/m0_m4_sem.yml"), exists=True),
) -> None:
    config = yaml.safe_load(sem_config.read_text(encoding="utf-8"))
    minimum = int(config.get("reporting", {}).get("minimum_complete_islands", 20))
    metrics = pd.read_csv(island_metrics_csv)
    applicability = pd.read_csv(applicability_csv)
    readiness, summary = audit_readiness(metrics, applicability, minimum)
    output_dir.mkdir(parents=True, exist_ok=True)
    readiness.to_csv(output_dir / "m0_m4_sem_readiness.csv", index=False)
    (output_dir / "m0_m4_sem_readiness.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
