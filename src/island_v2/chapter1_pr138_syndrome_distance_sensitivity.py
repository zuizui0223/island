"""Distance functional-form sensitivity for PR138 syndrome conclusions.

The island set is unchanged and no distance threshold is introduced. The same
canonical island syndrome scores are fitted against raw, square-root, and log1p
mainland distance, with predictors standardized inside the existing model.
"""

from __future__ import annotations

import copy
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_syndrome_analysis import (
    _between_contexts,
    _bh,
    _prepare,
    _within_context,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _fit_family(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = _prepare(island_scores, covariates, pattern_config)
    contexts = [str(x) for x in pattern_config["contexts"]]
    strata = [
        s
        for s in ("all_native", "native_nonendemic")
        if s in {str(x) for x in pattern_config["strata"]}
    ]
    threshold = int(pattern_config["support_tiers"]["confirmatory"])

    slope_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        start_w = len(within_rows)
        start_b = len(between_rows)
        for context in contexts:
            slopes, result = _within_context(
                data,
                stratum=stratum,
                context_value=context,
                support_tier="confirmatory",
                threshold=threshold,
                pattern_config=pattern_config,
                syndrome_config=syndrome_config,
            )
            if not slopes.empty:
                slope_parts.append(slopes)
            within_rows.append(result)
        for context_a, context_b in itertools.combinations(contexts, 2):
            between_rows.append(
                _between_contexts(
                    data,
                    stratum=stratum,
                    context_a=context_a,
                    context_b=context_b,
                    support_tier="confirmatory",
                    threshold=threshold,
                    pattern_config=pattern_config,
                )
            )

        # FDR is applied within the same stratum/tier family used canonically.
        w = pd.DataFrame(within_rows[start_w:])
        if "northern_classic_projection_one_sided_p" in w.columns:
            q = _bh(w["northern_classic_projection_one_sided_p"]).to_numpy()
            for row, value in zip(within_rows[start_w:], q, strict=True):
                row["q_classic_projection"] = value
        b = pd.DataFrame(between_rows[start_b:])
        if "p_value" in b.columns:
            q = _bh(b["p_value"]).to_numpy()
            for row, value in zip(between_rows[start_b:], q, strict=True):
                row["q_vector_difference"] = value

    return (
        pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame(),
        pd.DataFrame(within_rows),
        pd.DataFrame(between_rows),
    )


def _headline(
    form: str,
    slopes: pd.DataFrame,
    within: pd.DataFrame,
    between: pd.DataFrame,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if within.empty or "stratum" not in within.columns:
        return rows
    fitted_strata = sorted({str(x) for x in within["stratum"].dropna().unique()})
    for stratum in fitted_strata:
        north = within.loc[
            within["stratum"].eq(stratum)
            & within["context"].eq("northern_midlatitude")
        ]
        tropical = within.loc[
            within["stratum"].eq(stratum) & within["context"].eq("tropical")
        ]
        b = between.loc[
            between["stratum"].eq(stratum)
            & between["context_a"].eq("northern_midlatitude")
            & between["context_b"].eq("tropical")
        ]
        ns = slopes.loc[
            slopes["stratum"].eq(stratum)
            & slopes["context"].eq("northern_midlatitude")
        ]

        def slope(name: str) -> float:
            hit = ns.loc[ns["syndrome"].eq(name), "distance_slope"]
            return float(hit.iloc[0]) if not hit.empty else np.nan

        nq = float(north.iloc[0].get("q_classic_projection", np.nan)) if not north.empty else np.nan
        bq = float(b.iloc[0].get("q_vector_difference", np.nan)) if not b.empty else np.nan
        tq = float(tropical.iloc[0].get("q_classic_projection", np.nan)) if not tropical.empty else np.nan
        large_bee = slope("large_bee_like")
        generalized = slope("generalized_accessible")
        rows.append(
            {
                "distance_form": form,
                "stratum": stratum,
                "north_classic_projection": (
                    float(north.iloc[0].get("northern_classic_projection", np.nan))
                    if not north.empty
                    else np.nan
                ),
                "north_classic_q": nq,
                "north_classic_supported": bool(np.isfinite(nq) and nq <= 0.05),
                "tropical_classic_q": tq,
                "tropical_classic_supported": bool(np.isfinite(tq) and tq <= 0.05),
                "north_tropical_vector_q": bq,
                "north_tropical_vector_difference_supported": bool(
                    np.isfinite(bq) and bq <= 0.05
                ),
                "north_large_bee_slope": large_bee,
                "north_generalized_slope": generalized,
                "north_attraction_direction_supported": bool(
                    np.isfinite(large_bee)
                    and np.isfinite(generalized)
                    and large_bee < 0
                    and generalized > 0
                ),
            }
        )
    return rows


def run_distance_sensitivity(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw_column = "distance_to_continent_km"
    if raw_column not in covariates.columns:
        raise typer.BadParameter(f"covariates missing column: {raw_column}")
    cov = covariates.copy()
    raw = pd.to_numeric(cov[raw_column], errors="coerce")
    cov["pr138_distance_raw"] = raw
    cov["pr138_distance_sqrt"] = np.sqrt(raw.clip(lower=0))
    cov["pr138_distance_log1p"] = np.log1p(raw.clip(lower=0))

    slope_parts: list[pd.DataFrame] = []
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    headline_rows: list[dict[str, Any]] = []
    forms = {
        "raw": "pr138_distance_raw",
        "sqrt": "pr138_distance_sqrt",
        "log1p": "pr138_distance_log1p",
    }
    for form, column in forms.items():
        cfg = copy.deepcopy(pattern_config)
        cfg["geography_column"] = column
        slopes, within, between = _fit_family(
            island_scores,
            cov,
            cfg,
            syndrome_config,
        )
        for frame, parts in (
            (slopes, slope_parts),
            (within, within_parts),
            (between, between_parts),
        ):
            if not frame.empty:
                tagged = frame.copy()
                tagged.insert(0, "distance_form", form)
                parts.append(tagged)
        headline_rows.extend(_headline(form, slopes, within, between))

    return (
        pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame(),
        pd.concat(within_parts, ignore_index=True) if within_parts else pd.DataFrame(),
        pd.concat(between_parts, ignore_index=True) if between_parts else pd.DataFrame(),
        pd.DataFrame(headline_rows),
    )


@app.command("run")
def run(
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    syndrome_config = yaml.safe_load(syndrome_config_path.read_text(encoding="utf-8"))
    slopes, within, between, headline = run_distance_sensitivity(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pattern_config,
        syndrome_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    slopes.to_csv(output_dir / "syndrome_distance_form_slopes.csv", index=False)
    within.to_csv(output_dir / "syndrome_distance_form_within.csv", index=False)
    between.to_csv(output_dir / "syndrome_distance_form_between.csv", index=False)
    headline.to_csv(output_dir / "syndrome_distance_form_headline.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_syndrome_distance_form_v1",
        "forms": ["raw", "sqrt", "log1p"],
        "island_set_changed": False,
        "distance_threshold_used": False,
    }
    (output_dir / "syndrome_distance_form_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
