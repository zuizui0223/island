"""Run the PR138 attraction-vs-selfing decomposition after external source-pool adjustment.

This sensitivity asks whether the attraction/accessibility distance association survives after
subtracting the same prespecified GIFT mainland expectation used for the syndrome analysis,
and whether the strict reproductive selfing core behaves differently. It remains a conditional
decomposition, not causal mediation and not a historical source reconstruction.
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from island_v2.chapter1_pr138_pathway_decomposition import run_pathway_decomposition
from island_v2.chapter1_pr138_source_pool_sensitivity import (
    build_adjusted_island_scores,
    build_mainland_syndrome_scores,
    build_source_expectations,
    match_gift_species_to_frozen_scores,
)


def run_source_adjusted_pathway(
    species_scores: pd.DataFrame,
    island_scores: pd.DataFrame,
    gift_flora: pd.DataFrame,
    source_assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    source_config: dict[str, Any],
    realm_assignment: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    min_species = int(source_config["source_region_scores"]["minimum_trait_scored_species_per_region_syndrome"])
    min_sources = int(source_config["response"]["source_expectation_requires_minimum_source_regions"])
    strata = [str(x) for x in source_config["strata"]]
    modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]

    matched, match_audit = match_gift_species_to_frozen_scores(gift_flora, species_scores)
    mainland_scores = build_mainland_syndrome_scores(matched, min_species=min_species)
    expectations = build_source_expectations(
        source_assignments,
        mainland_scores,
        min_source_regions=min_sources,
    )
    adjusted = build_adjusted_island_scores(island_scores, expectations, strata)

    regime_parts: list[pd.DataFrame] = []
    realm_parts: list[pd.DataFrame] = []
    for mode in modes:
        mode_scores = adjusted.loc[adjusted["source_mode"].eq(mode)].copy()
        config = deepcopy(pattern_config)
        config["strata"] = strata
        result = run_pathway_decomposition(mode_scores, covariates, config)
        result.insert(0, "source_mode", mode)
        result.insert(0, "context_layer", "analysis_regime")
        regime_parts.append(result)

        if realm_assignment is not None and not realm_assignment.empty:
            realm_cov = covariates.merge(
                realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
                on="island_id",
                how="left",
                validate="one_to_one",
            )
            realm_config = deepcopy(config)
            realm_config["context_column"] = "biogeographic_realm"
            realm_config["contexts"] = sorted(
                realm_cov["biogeographic_realm"].dropna().astype(str).loc[lambda s: s.ne("")].unique()
            )
            realm_result = run_pathway_decomposition(mode_scores, realm_cov, realm_config)
            realm_result.insert(0, "source_mode", mode)
            realm_result.insert(0, "context_layer", "biogeographic_realm")
            realm_parts.append(realm_result)

    regime = pd.concat(regime_parts, ignore_index=True) if regime_parts else pd.DataFrame()
    realm = pd.concat(realm_parts, ignore_index=True) if realm_parts else pd.DataFrame()
    manifest = {
        "contract": "chapter1_pr138_source_adjusted_pathway_v1",
        "source_modes": modes,
        "strata": strata,
        "attraction_shift": "(-large_bee_like + generalized_accessible) / 2",
        "selfing_core": "SC + mating system + autonomous selfing; floral size excluded",
        "source_adjustment_applied_before_pathway_decomposition": True,
        "causal_mediation_claimed": False,
        "source_assignment_uses_island_traits": False,
        "claim_boundary": "Differential persistence after source adjustment separates measured response components; it does not identify pollinator loss or causal selection.",
    }
    return {
        "match_audit": match_audit,
        "mainland_scores": mainland_scores,
        "expectations": expectations,
        "adjusted_scores": adjusted,
        "regime": regime,
        "realm": realm,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-scores-csv", type=Path, required=True)
    parser.add_argument("--island-scores-csv", type=Path, required=True)
    parser.add_argument("--gift-flora-csv", type=Path, required=True)
    parser.add_argument("--source-assignments-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--source-config-path", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pattern_config = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    source_config = yaml.safe_load(args.source_config_path.read_text(encoding="utf-8"))
    realm = pd.read_csv(args.realm_assignment_csv) if args.realm_assignment_csv else None
    outputs = run_source_adjusted_pathway(
        pd.read_csv(args.species_scores_csv),
        pd.read_csv(args.island_scores_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pattern_config,
        source_config,
        realm,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in {
        "match_audit": "gift_pathway_species_match_audit.csv",
        "mainland_scores": "mainland_pathway_syndrome_scores.csv.gz",
        "expectations": "island_pathway_source_expectations.csv.gz",
        "adjusted_scores": "source_adjusted_pathway_island_scores.csv.gz",
        "regime": "source_adjusted_pathway_regime.csv",
        "realm": "source_adjusted_pathway_realm.csv",
    }.items():
        frame = outputs[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(args.output_dir / filename, index=False)
    manifest = outputs["manifest"]
    assert isinstance(manifest, dict)
    (args.output_dir / "source_adjusted_pathway_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
