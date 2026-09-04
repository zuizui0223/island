"""Information-weight sensitivity for the PR138 eight-axis observed syndrome.

Each island's observed broad-trait share is held fixed while its effective trial count
is varied. This checks whether highly trait-resolved islands are necessary for the
northern classic-syndrome projection or the north-vs-tropical vector difference.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_information_weight_sensitivity import reweight_composition
from island_v2.chapter1_pr138_biogeographic_pattern import run_observed_pattern

app = typer.Typer(add_completion=False, no_args_is_help=True)


def run_pr138_information_weight(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    modes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    configured_strata = {str(x) for x in pattern_config.get("strata", [])}
    headline_strata = [
        x for x in ("all_native", "native_nonendemic") if x in configured_strata
    ]
    for mode in modes:
        weighted = reweight_composition(counts, mode)
        _, _, within, between = run_observed_pattern(weighted, covariates, pattern_config)
        if not within.empty:
            within = within.copy()
            within.insert(0, "information_weight_mode", mode)
            within_parts.append(within)
        if not between.empty:
            between = between.copy()
            between.insert(0, "information_weight_mode", mode)
            between_parts.append(between)

        for stratum in headline_strata:
            w = within.loc[
                within["stratum"].eq(stratum)
                & within["support_tier"].eq("confirmatory")
            ] if not within.empty else pd.DataFrame()
            north = w.loc[w["context"].eq("northern_midlatitude")] if not w.empty else pd.DataFrame()
            tropical = w.loc[w["context"].eq("tropical")] if not w.empty else pd.DataFrame()
            b = between.loc[
                between["stratum"].eq(stratum)
                & between["support_tier"].eq("confirmatory")
                & between["context_a"].eq("northern_midlatitude")
                & between["context_b"].eq("tropical")
            ] if not between.empty else pd.DataFrame()
            summary_rows.append(
                {
                    "information_weight_mode": mode,
                    "stratum": stratum,
                    "north_testable": bool(not north.empty and north.iloc[0]["status"] == "fit"),
                    "north_vector_q": float(north.iloc[0].get("q_joint_within_stratum_tier", np.nan)) if not north.empty else np.nan,
                    "north_classic_projection": float(north.iloc[0].get("classic_projection_estimate", np.nan)) if not north.empty else np.nan,
                    "north_classic_projection_q": float(north.iloc[0].get("q_classic_projection_within_stratum_tier", np.nan)) if not north.empty else np.nan,
                    "north_classic_supported": bool(not north.empty and north.iloc[0].get("classic_direction_supported", False)),
                    "tropical_classic_projection": float(tropical.iloc[0].get("classic_projection_estimate", np.nan)) if not tropical.empty else np.nan,
                    "tropical_classic_projection_q": float(tropical.iloc[0].get("q_classic_projection_within_stratum_tier", np.nan)) if not tropical.empty else np.nan,
                    "north_tropical_difference_q": float(b.iloc[0].get("q_joint_within_stratum_tier", np.nan)) if not b.empty else np.nan,
                    "north_tropical_difference_supported": bool(not b.empty and b.iloc[0].get("biogeographic_difference_supported", False)),
                }
            )
    return (
        pd.concat(within_parts, ignore_index=True) if within_parts else pd.DataFrame(),
        pd.concat(between_parts, ignore_index=True) if between_parts else pd.DataFrame(),
        pd.DataFrame(summary_rows),
    )


@app.command("run")
def run(
    counts_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    pattern_config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True
    ),
    modes: str = typer.Option("canonical,cap_100,cap_50,cap_20,equal_island"),
) -> None:
    config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    mode_list = [x.strip() for x in modes.split(",") if x.strip()]
    within, between, summary = run_pr138_information_weight(
        pd.read_csv(counts_csv),
        pd.read_csv(covariates_csv),
        config,
        mode_list,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "pr138_information_weight_within.csv", index=False)
    between.to_csv(output_dir / "pr138_information_weight_between.csv", index=False)
    summary.to_csv(output_dir / "pr138_information_weight_summary.csv", index=False)
    manifest = {
        "contract": "pr138_eight_axis_information_weight_v1",
        "modes": mode_list,
        "trait_shares_fixed_within_island": True,
        "interpretation": "Pseudo-likelihood sensitivity; not an abundance model.",
    }
    (output_dir / "pr138_information_weight_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
