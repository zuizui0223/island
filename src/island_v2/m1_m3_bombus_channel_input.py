"""Prepare component-first Bombus channel inputs for M1-M3.

The M1-M3 analysis must not collapse Bombus evidence to climate compatibility alone.
This module joins effort-aware observation evidence with environmental compatibility
while preserving both components and the secondary combined score.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_OBSERVATION_COLUMNS = {
    "island_id",
    "bombus_occurrence_evidence",
    "bombus_channel_evaluable",
    "observation_availability",
    "observation_deficit",
    "observation_effort_quality",
}
REQUIRED_ENV_COLUMNS = {
    "island_id",
    "environmental_compatibility_max",
}


def build_channel_input(observation: pd.DataFrame, environmental: pd.DataFrame) -> pd.DataFrame:
    missing_obs = REQUIRED_OBSERVATION_COLUMNS.difference(observation.columns)
    if missing_obs:
        raise typer.BadParameter(f"observation table missing columns: {sorted(missing_obs)}")
    missing_env = REQUIRED_ENV_COLUMNS.difference(environmental.columns)
    if missing_env:
        raise typer.BadParameter(f"environment table missing columns: {sorted(missing_env)}")

    obs = observation.drop_duplicates("island_id").copy()
    env = environmental.drop_duplicates("island_id").copy()
    out = obs.merge(env, on="island_id", how="outer", validate="one_to_one")

    out["environmental_compatibility"] = pd.to_numeric(
        out["environmental_compatibility_max"], errors="coerce"
    )
    out["observation_availability"] = pd.to_numeric(
        out["observation_availability"], errors="coerce"
    )
    out["observation_effort_quality"] = pd.to_numeric(
        out["observation_effort_quality"], errors="coerce"
    )

    both = out["environmental_compatibility"].notna() & out["observation_availability"].notna()
    out["bombus_channel_availability"] = pd.NA
    out.loc[both, "bombus_channel_availability"] = (
        out.loc[both, "environmental_compatibility"]
        * out.loc[both, "observation_availability"]
    ) ** 0.5
    out["bombus_channel_deficit"] = pd.to_numeric(
        1.0 - pd.to_numeric(out["bombus_channel_availability"], errors="coerce"),
        errors="coerce",
    )

    # Primary M1-M3 inference requires interpretable observation evidence and a
    # fitted environmental score. Sparse observation coverage remains visible but
    # is not silently treated as absence.
    out["m1_m3_channel_evaluable"] = (
        out["bombus_channel_evaluable"].fillna(False).astype(bool)
        & out["environmental_compatibility"].notna()
    )
    out["channel_component_status"] = "complete"
    out.loc[out["observation_availability"].isna(), "channel_component_status"] = "missing_observation"
    out.loc[out["environmental_compatibility"].isna(), "channel_component_status"] = "missing_environment"
    out.loc[
        out["observation_availability"].isna() & out["environmental_compatibility"].isna(),
        "channel_component_status",
    ] = "missing_both"
    return out


@app.command("run")
def run(
    observation_scores_csv: Path = typer.Option(..., exists=True),
    environmental_island_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    observation = pd.read_csv(observation_scores_csv)
    environmental = pd.read_csv(environmental_island_csv)
    result = build_channel_input(observation, environmental)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "m1_m3_bombus_channel_input.csv", index=False)

    summary = {
        "contract": "m1_m3_component_first_bombus_channel_v1",
        "n_islands": int(result["island_id"].nunique()),
        "n_m1_m3_channel_evaluable": int(result["m1_m3_channel_evaluable"].sum()),
        "occurrence_evidence_counts": {
            str(k): int(v)
            for k, v in result["bombus_occurrence_evidence"].fillna("missing").value_counts().items()
        },
        "component_status_counts": {
            str(k): int(v) for k, v in result["channel_component_status"].value_counts().items()
        },
        "policy": (
            "Observation availability and environmental compatibility remain separate primary components; "
            "their geometric-mean availability is a secondary summary only."
        ),
    }
    (output_dir / "m1_m3_bombus_channel_manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary, indent=2))


if __name__ == "__main__":
    app()
