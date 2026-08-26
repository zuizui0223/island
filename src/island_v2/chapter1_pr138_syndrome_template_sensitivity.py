"""Template-robustness analysis for PR138 pollination/selfing syndrome scores.

The canonical syndrome templates are scientific hypotheses, not observed pollinator
identities. This module therefore perturbs those templates *without looking at the
outcomes* and reruns exactly the same biogeographic analysis.

Prespecified variants:
- canonical literature weights;
- equal weights;
- no flower colour;
- pollination morphology only (form/symmetry/tube; reproductive templates unchanged);
- global leave-one-trait-out for every trait represented in the templates;
- require at least three informative traits per syndrome score.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_syndrome_analysis import run_syndrome_analysis

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _all_template_traits(config: dict[str, Any]) -> list[str]:
    return sorted(
        {
            str(trait)
            for spec in config["syndromes"].values()
            for trait in spec["traits"]
        }
    )


def build_template_variants(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    variants: dict[str, dict[str, Any]] = {"canonical": copy.deepcopy(config)}

    equal = copy.deepcopy(config)
    for spec in equal["syndromes"].values():
        for trait_spec in spec["traits"].values():
            trait_spec["weight"] = 1.0
    variants["equal_weights"] = equal

    no_colour = copy.deepcopy(config)
    for spec in no_colour["syndromes"].values():
        spec["traits"].pop("flower_primary_color", None)
    variants["no_colour"] = no_colour

    morphology = copy.deepcopy(config)
    floral_keep = {"floral_form", "floral_symmetry", "tube_depth_class"}
    reproductive = {"selfing_syndrome", "selfing_core"}
    for syndrome, spec in morphology["syndromes"].items():
        if syndrome in reproductive:
            continue
        spec["traits"] = {
            trait: trait_spec
            for trait, trait_spec in spec["traits"].items()
            if trait in floral_keep
        }
    variants["pollination_morphology_only"] = morphology

    for trait in _all_template_traits(config):
        variant = copy.deepcopy(config)
        for spec in variant["syndromes"].values():
            spec["traits"].pop(trait, None)
        variants[f"drop_{trait}"] = variant

    min3 = copy.deepcopy(config)
    min3["score_definition"]["minimum_informative_traits"] = 3
    variants["minimum_three_traits"] = min3
    return variants


def _headline(
    variant: str,
    slopes: pd.DataFrame,
    within: pd.DataFrame,
    between: pd.DataFrame,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stratum in ("all_native", "native_nonendemic"):
        w = within.loc[
            within["stratum"].eq(stratum)
            & within["support_tier"].eq("confirmatory")
        ]
        north = w.loc[w["context"].eq("northern_midlatitude")]
        trop = w.loc[w["context"].eq("tropical")]
        b = between.loc[
            between["stratum"].eq(stratum)
            & between["support_tier"].eq("confirmatory")
            & between["context_a"].eq("northern_midlatitude")
            & between["context_b"].eq("tropical")
        ]
        slope_subset = slopes.loc[
            slopes["stratum"].eq(stratum)
            & slopes["support_tier"].eq("confirmatory")
        ]

        def slope(context: str, syndrome: str) -> float:
            hit = slope_subset.loc[
                slope_subset["context"].eq(context)
                & slope_subset["syndrome"].eq(syndrome),
                "distance_slope",
            ]
            return float(hit.iloc[0]) if not hit.empty else float("nan")

        rows.append(
            {
                "variant": variant,
                "stratum": stratum,
                "north_testable": bool(not north.empty and north.iloc[0]["status"] == "fit"),
                "north_classic_projection": float(north.iloc[0].get("northern_classic_projection", float("nan"))) if not north.empty else float("nan"),
                "north_classic_q": float(north.iloc[0].get("q_classic_projection", float("nan"))) if not north.empty else float("nan"),
                "north_classic_supported": bool(not north.empty and north.iloc[0].get("classic_projection_supported", False)),
                "north_large_bee_slope": slope("northern_midlatitude", "large_bee_like"),
                "north_generalized_slope": slope("northern_midlatitude", "generalized_accessible"),
                "north_selfing_slope": slope("northern_midlatitude", "selfing_syndrome"),
                "tropical_butterfly_slope": slope("tropical", "butterfly_like"),
                "tropical_bird_slope": slope("tropical", "bird_like"),
                "tropical_large_bee_slope": slope("tropical", "large_bee_like"),
                "tropical_generalized_slope": slope("tropical", "generalized_accessible"),
                "tropical_selfing_slope": slope("tropical", "selfing_syndrome"),
                "north_tropical_vector_q": float(b.iloc[0].get("q_within_stratum_tier", float("nan"))) if not b.empty else float("nan"),
                "north_tropical_vector_difference_supported": bool(not b.empty and b.iloc[0].get("regional_syndrome_vector_difference_supported", False)),
                "tropical_classic_supported": bool(not trop.empty and trop.iloc[0].get("classic_projection_supported", False)),
            }
        )
    return rows


def run_template_sensitivity(
    trait_ledger: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    slope_parts: list[pd.DataFrame] = []
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    headline_rows: list[dict[str, Any]] = []

    for name, variant in build_template_variants(syndrome_config).items():
        _, _, slopes, within, between = run_syndrome_analysis(
            trait_ledger,
            status_flora,
            covariates,
            pattern_config,
            variant,
        )
        for frame, parts in (
            (slopes, slope_parts),
            (within, within_parts),
            (between, between_parts),
        ):
            if not frame.empty:
                tagged = frame.copy()
                tagged.insert(0, "variant", name)
                parts.append(tagged)
        headline_rows.extend(_headline(name, slopes, within, between))

    return (
        pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame(),
        pd.concat(within_parts, ignore_index=True) if within_parts else pd.DataFrame(),
        pd.concat(between_parts, ignore_index=True) if between_parts else pd.DataFrame(),
        pd.DataFrame(headline_rows),
    )


@app.command("run")
def run(
    trait_ledger_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    syndrome_config = yaml.safe_load(syndrome_config_path.read_text(encoding="utf-8"))
    slopes, within, between, headline = run_template_sensitivity(
        pd.read_csv(trait_ledger_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        pattern_config,
        syndrome_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    slopes.to_csv(output_dir / "template_sensitivity_slopes.csv", index=False)
    within.to_csv(output_dir / "template_sensitivity_within.csv", index=False)
    between.to_csv(output_dir / "template_sensitivity_between.csv", index=False)
    headline.to_csv(output_dir / "template_sensitivity_headline.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_syndrome_template_sensitivity_v1",
        "variants": list(build_template_variants(syndrome_config)),
        "outcome_blind_template_perturbations": True,
        "pollinator_identity_inferred": False,
    }
    (output_dir / "template_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
