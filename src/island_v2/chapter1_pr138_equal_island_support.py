"""Equal-island support-threshold sensitivity for PR138.

GBIF-derived observed floras are opportunistic samples, so the number of directly
trait-resolved species is not assumed to be a literal binomial sampling effort. This
analysis gives each eligible island equal information weight, while separately raising
the minimum number of direct species required for each island x syndrome axis.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_information_weight_sensitivity import reweight_composition
from island_v2.chapter1_pr138_biogeographic_pattern import run_observed_pattern

app = typer.Typer(add_completion=False, no_args_is_help=True)


def run_equal_island_support(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
    min_trials_values: list[int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    within_slope_parts: list[pd.DataFrame] = []
    between_slope_parts: list[pd.DataFrame] = []
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    for minimum in min_trials_values:
        filtered = counts.loc[
            pd.to_numeric(counts["trials"], errors="coerce").ge(int(minimum))
        ].copy()
        weighted = reweight_composition(filtered, "equal_island")
        within_slopes, between_slopes, within, between = run_observed_pattern(
            weighted, covariates, config
        )
        for frame, parts in (
            (within_slopes, within_slope_parts),
            (between_slopes, between_slope_parts),
            (within, within_parts),
            (between, between_parts),
        ):
            if not frame.empty:
                frame = frame.copy()
                frame.insert(0, "minimum_direct_species_per_island_axis", int(minimum))
                parts.append(frame)

        for stratum in ("all_native", "native_nonendemic"):
            if stratum not in {str(x) for x in config.get("strata", [])}:
                continue
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
                    "minimum_direct_species_per_island_axis": int(minimum),
                    "stratum": stratum,
                    "north_status": str(north.iloc[0]["status"]) if not north.empty else "not_testable",
                    "north_n_retained_outcomes": int(north.iloc[0].get("n_retained_outcomes", 0)) if not north.empty else 0,
                    "north_n_unique_islands": float(north.iloc[0].get("n_unique_islands", float("nan"))) if not north.empty else float("nan"),
                    "north_vector_q": float(north.iloc[0].get("q_joint_within_stratum_tier", float("nan"))) if not north.empty else float("nan"),
                    "north_classic_projection": float(north.iloc[0].get("classic_projection_estimate", float("nan"))) if not north.empty else float("nan"),
                    "north_classic_projection_q": float(north.iloc[0].get("q_classic_projection_within_stratum_tier", float("nan"))) if not north.empty else float("nan"),
                    "north_classic_supported": bool(not north.empty and north.iloc[0].get("classic_direction_supported", False)),
                    "tropical_status": str(tropical.iloc[0]["status"]) if not tropical.empty else "not_testable",
                    "tropical_classic_projection": float(tropical.iloc[0].get("classic_projection_estimate", float("nan"))) if not tropical.empty else float("nan"),
                    "tropical_classic_projection_q": float(tropical.iloc[0].get("q_classic_projection_within_stratum_tier", float("nan"))) if not tropical.empty else float("nan"),
                    "north_tropical_status": str(b.iloc[0]["status"]) if not b.empty else "not_testable",
                    "north_tropical_difference_q": float(b.iloc[0].get("q_joint_within_stratum_tier", float("nan"))) if not b.empty else float("nan"),
                    "north_tropical_difference_supported": bool(not b.empty and b.iloc[0].get("biogeographic_difference_supported", False)),
                }
            )

    concat = lambda parts: pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return (
        concat(within_slope_parts),
        concat(between_slope_parts),
        concat(within_parts),
        concat(between_parts),
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
    minimum_direct_species: str = typer.Option("1,3,5,10"),
) -> None:
    config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    thresholds = [int(x.strip()) for x in minimum_direct_species.split(",") if x.strip()]
    within_slopes, between_slopes, within, between, summary = run_equal_island_support(
        pd.read_csv(counts_csv), pd.read_csv(covariates_csv), config, thresholds
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within_slopes.to_csv(output_dir / "equal_island_support_within_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "equal_island_support_between_slopes.csv", index=False)
    within.to_csv(output_dir / "equal_island_support_within.csv", index=False)
    between.to_csv(output_dir / "equal_island_support_between.csv", index=False)
    summary.to_csv(output_dir / "equal_island_support_summary.csv", index=False)
    manifest = {
        "contract": "pr138_equal_island_support_v1",
        "minimum_direct_species_per_island_axis": thresholds,
        "effective_information_weight_per_eligible_island_axis": 1,
        "interpretation": (
            "Separates equal-island inference from the reliability gate imposed by "
            "minimum direct species support."
        ),
    }
    (output_dir / "equal_island_support_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
