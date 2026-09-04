"""Outcome-blind template robustness for PR138 pollination/selfing syndromes.

Signed species x syndrome x trait concordances are computed once from the canonical
preferred/opposed state sets. Prespecified sensitivity variants then change only trait
inclusion, weights, or the minimum informative-trait requirement. This cache changes
runtime only, never the biological state mapping or estimand.
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
    _trait_concordance,
    _within_context,
    build_island_syndrome_scores,
)

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


def precompute_trait_concordances(
    trait_ledger: pd.DataFrame,
    canonical_config: dict[str, Any],
) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "normalized_value"}
    missing = required - set(trait_ledger.columns)
    if missing:
        raise typer.BadParameter(f"trait ledger missing columns: {sorted(missing)}")

    ledger = trait_ledger[list(required)].copy()
    for column in required:
        ledger[column] = ledger[column].fillna("").astype(str).str.strip()
    ledger = ledger.loc[
        ledger["accepted_species"].ne("") & ledger["trait_name"].ne("")
    ].drop_duplicates(["accepted_species", "trait_name"], keep="first")

    parts: list[pd.DataFrame] = []
    for syndrome, syndrome_spec in canonical_config["syndromes"].items():
        for trait_name, trait_spec in syndrome_spec["traits"].items():
            subset = ledger.loc[
                ledger["trait_name"].eq(str(trait_name)),
                ["accepted_species", "normalized_value"],
            ].copy()
            if subset.empty:
                continue
            preferred = {str(x) for x in trait_spec.get("preferred", [])}
            opposed = {str(x) for x in trait_spec.get("opposed", [])}
            subset["trait_concordance"] = [
                _trait_concordance(value, preferred, opposed)
                for value in subset["normalized_value"]
            ]
            subset = subset.loc[
                np.isfinite(pd.to_numeric(subset["trait_concordance"], errors="coerce"))
            ].copy()
            if subset.empty:
                continue
            subset["syndrome"] = str(syndrome)
            subset["trait_name"] = str(trait_name)
            parts.append(
                subset[["accepted_species", "syndrome", "trait_name", "trait_concordance"]]
            )

    if not parts:
        return pd.DataFrame(
            columns=["accepted_species", "syndrome", "trait_name", "trait_concordance"]
        )
    return pd.concat(parts, ignore_index=True)


def score_variant_from_cache(
    concordance_cache: pd.DataFrame,
    variant_config: dict[str, Any],
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    minimum = int(variant_config["score_definition"]["minimum_informative_traits"])

    for syndrome, spec in variant_config["syndromes"].items():
        weight_map = {
            str(trait): float(trait_spec.get("weight", 1.0))
            for trait, trait_spec in spec["traits"].items()
        }
        if not weight_map:
            continue
        work = concordance_cache.loc[
            concordance_cache["syndrome"].eq(str(syndrome))
            & concordance_cache["trait_name"].isin(weight_map),
            ["accepted_species", "trait_name", "trait_concordance"],
        ].copy()
        if work.empty:
            continue
        work["weight"] = work["trait_name"].map(weight_map).astype(float)
        work["weighted"] = work["trait_concordance"].astype(float) * work["weight"]
        grouped = (
            work.groupby("accepted_species", as_index=False)
            .agg(
                weighted_sum=("weighted", "sum"),
                weight_total=("weight", "sum"),
                n_informative_traits=("trait_name", "nunique"),
                informative_traits=(
                    "trait_name",
                    lambda x: "|".join(sorted(set(str(v) for v in x))),
                ),
            )
        )
        grouped = grouped.loc[
            grouped["n_informative_traits"].ge(minimum) & grouped["weight_total"].gt(0)
        ].copy()
        if grouped.empty:
            continue
        grouped["syndrome_concordance"] = grouped["weighted_sum"] / grouped["weight_total"]
        grouped["soft_membership"] = (grouped["syndrome_concordance"] + 1.0) / 2.0
        grouped["syndrome"] = str(syndrome)
        frames.append(
            grouped[
                [
                    "accepted_species",
                    "syndrome",
                    "syndrome_concordance",
                    "soft_membership",
                    "n_informative_traits",
                    "informative_traits",
                ]
            ]
        )

    if not frames:
        return pd.DataFrame(
            columns=[
                "accepted_species",
                "syndrome",
                "syndrome_concordance",
                "soft_membership",
                "n_informative_traits",
                "informative_traits",
            ]
        )
    return pd.concat(frames, ignore_index=True)


def analyze_cached_species_scores(
    species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    strata = [str(x) for x in pattern_config["strata"]]
    island_scores = build_island_syndrome_scores(status_flora, species_scores, strata)
    data = _prepare(island_scores, covariates, pattern_config)
    contexts = [str(x) for x in pattern_config["contexts"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}

    slope_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for tier, threshold in tiers.items():
            for context_value in contexts:
                slopes, result = _within_context(
                    data,
                    stratum=stratum,
                    context_value=context_value,
                    support_tier=tier,
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
                        support_tier=tier,
                        threshold=threshold,
                        pattern_config=pattern_config,
                    )
                )

    slopes = pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame()
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)

    if "p_value" in within.columns:
        within["q_within_stratum_tier"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        within["syndrome_vector_supported"] = (
            within["q_within_stratum_tier"].le(0.05).fillna(False)
        )
    if "northern_classic_projection_one_sided_p" in within.columns:
        within["q_classic_projection"] = within.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["northern_classic_projection_one_sided_p"].apply(_bh)
        within["classic_projection_supported"] = (
            within["q_classic_projection"].le(0.05).fillna(False)
        )
    if "p_value" in between.columns:
        between["q_within_stratum_tier"] = between.groupby(
            ["stratum", "support_tier"], group_keys=False
        )["p_value"].apply(_bh)
        between["regional_syndrome_vector_difference_supported"] = (
            between["q_within_stratum_tier"].le(0.05).fillna(False)
        )
    return slopes, within, between


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
                "north_classic_projection": (
                    float(north.iloc[0].get("northern_classic_projection", float("nan")))
                    if not north.empty
                    else float("nan")
                ),
                "north_classic_q": (
                    float(north.iloc[0].get("q_classic_projection", float("nan")))
                    if not north.empty
                    else float("nan")
                ),
                "north_classic_supported": bool(
                    not north.empty and north.iloc[0].get("classic_projection_supported", False)
                ),
                "north_large_bee_slope": slope("northern_midlatitude", "large_bee_like"),
                "north_generalized_slope": slope(
                    "northern_midlatitude", "generalized_accessible"
                ),
                "north_selfing_slope": slope("northern_midlatitude", "selfing_syndrome"),
                "tropical_butterfly_slope": slope("tropical", "butterfly_like"),
                "tropical_bird_slope": slope("tropical", "bird_like"),
                "tropical_large_bee_slope": slope("tropical", "large_bee_like"),
                "tropical_generalized_slope": slope("tropical", "generalized_accessible"),
                "tropical_selfing_slope": slope("tropical", "selfing_syndrome"),
                "north_tropical_vector_q": (
                    float(b.iloc[0].get("q_within_stratum_tier", float("nan")))
                    if not b.empty
                    else float("nan")
                ),
                "north_tropical_vector_difference_supported": bool(
                    not b.empty
                    and b.iloc[0].get("regional_syndrome_vector_difference_supported", False)
                ),
                "tropical_classic_supported": bool(
                    not trop.empty and trop.iloc[0].get("classic_projection_supported", False)
                ),
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

    cache = precompute_trait_concordances(trait_ledger, syndrome_config)
    for name, variant in build_template_variants(syndrome_config).items():
        species_scores = score_variant_from_cache(cache, variant)
        slopes, within, between = analyze_cached_species_scores(
            species_scores,
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
    ledger = pd.read_csv(trait_ledger_csv)
    slopes, within, between, headline = run_template_sensitivity(
        ledger,
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
        "contract": "chapter1_pr138_syndrome_template_sensitivity_v2_cached",
        "variants": list(build_template_variants(syndrome_config)),
        "outcome_blind_template_perturbations": True,
        "trait_concordance_cache_reused_across_variants": True,
        "trait_concordance_cache_changes_estimand": False,
        "pollinator_identity_inferred": False,
    }
    (output_dir / "template_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
