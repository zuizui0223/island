"""Prospective V4 shared-architecture and lineage decomposition.

The sampled large-bee, butterfly, and bird syndrome templates overlap strongly in
their floral-trait definitions.  V4 therefore decomposes them into one common
architecture coordinate and three template-specific residuals.  The coordinate is
estimated only among species matched to the frozen GIFT mainland native flora.  No
island geography, regional label, distance response, or pollinator observation enters
the factor fit.

The outputs remain plant-trait concordances.  A residual labelled for a sampled
template is not evidence for the corresponding pollinator, and persistence after the
genus/source-pool sensitivities does not identify functional replacement.
"""

from __future__ import annotations

import argparse
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_global_branching import run_global_branching
from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _availability_matrices,
    _genus,
    _representation_count_matrix,
    _source_assignment_matrix,
    broad_source_availability,
    compute_island_enrichment,
    match_gift_species,
)
from island_v2.chapter1_pr138_source_pool_sensitivity import (
    build_adjusted_island_scores,
    build_mainland_syndrome_scores,
    build_source_expectations,
    match_gift_species_to_frozen_scores,
)
from island_v2.chapter1_pr138_syndrome_analysis import (
    _bh,
    build_island_syndrome_scores,
)

DEFAULT_TEMPLATES = ("large_bee_like", "butterfly_like", "bird_like")
SHARED_FACTOR = "shared_architecture_factor"
RESIDUAL_NAMES = {
    "large_bee_like": "large_bee_like_residual",
    "butterfly_like": "butterfly_like_residual",
    "bird_like": "bird_like_residual",
}
LINEAGE_OUTCOMES = (
    "entry_enrichment",
    "loading_increment",
    "beyond_genus_residual",
)


def _validate_species_scores(
    species_scores: pd.DataFrame, templates: tuple[str, ...]
) -> pd.DataFrame:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    if missing := required - set(species_scores.columns):
        raise ValueError(f"species scores missing columns: {sorted(missing)}")
    work = species_scores.loc[
        species_scores["syndrome"].astype(str).isin(templates),
        ["accepted_species", "syndrome", "syndrome_concordance"],
    ].copy()
    work["accepted_species"] = work["accepted_species"].fillna("").astype(str)
    work["syndrome_concordance"] = pd.to_numeric(
        work["syndrome_concordance"], errors="coerce"
    )
    work = work.dropna(subset=["syndrome_concordance"])
    duplicate = work.duplicated(["accepted_species", "syndrome"], keep=False)
    if duplicate.any():
        examples = work.loc[
            duplicate, ["accepted_species", "syndrome"]
        ].head(5)
        raise ValueError(f"duplicate species-template scores: {examples.to_dict('records')}")
    wide = work.pivot(
        index="accepted_species", columns="syndrome", values="syndrome_concordance"
    )
    missing_templates = set(templates) - set(wide.columns)
    if missing_templates:
        raise ValueError(f"sampled templates absent: {sorted(missing_templates)}")
    return wide.loc[:, list(templates)].dropna().sort_index()


def fit_and_project_source_factor(
    species_scores: pd.DataFrame,
    gift_flora: pd.DataFrame,
    *,
    templates: tuple[str, ...] = DEFAULT_TEMPLATES,
    minimum_complete_source_species: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Fit PCA1 on GIFT-matched complete cases, then project all complete cases."""

    wide = _validate_species_scores(species_scores, templates)
    species_universe = wide.reset_index()[["accepted_species"]]
    matched = match_gift_species(gift_flora, species_universe)
    source_species = sorted(
        set(matched["accepted_species"].astype(str)) & set(wide.index.astype(str))
    )
    if len(source_species) < int(minimum_complete_source_species):
        raise ValueError(
            "too few complete GIFT source species for frozen architecture factor: "
            f"{len(source_species)} < {minimum_complete_source_species}"
        )
    source = wide.loc[source_species, list(templates)].astype(float)
    means = source.mean(axis=0)
    scales = source.std(axis=0, ddof=0)
    if not np.isfinite(scales.to_numpy(float)).all() or scales.le(0).any():
        raise ValueError("constant or invalid sampled template in source training set")
    source_z = (source - means) / scales
    correlation = np.cov(source_z.to_numpy(float), rowvar=False, ddof=0)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    order = np.argsort(eigenvalues)[::-1]
    first_value = float(eigenvalues[order[0]])
    loadings = eigenvectors[:, order[0]].astype(float)
    source_factor = source_z.to_numpy(float) @ loadings
    correlations = np.array(
        [
            np.corrcoef(source_factor, source_z.iloc[:, index].to_numpy(float))[0, 1]
            for index in range(len(templates))
        ]
    )
    if float(np.nansum(correlations)) < 0:
        loadings *= -1.0
        source_factor *= -1.0
        correlations *= -1.0

    all_z = (wide.astype(float) - means) / scales
    factor = all_z.to_numpy(float) @ loadings
    residual = all_z.to_numpy(float) - np.outer(factor, loadings)
    component_wide = pd.DataFrame(
        {SHARED_FACTOR: factor}, index=wide.index.astype(str)
    )
    for index, template in enumerate(templates):
        component_wide[RESIDUAL_NAMES[template]] = residual[:, index]

    projected = (
        component_wide.rename_axis("accepted_species")
        .reset_index()
        .melt(
            id_vars="accepted_species",
            var_name="syndrome",
            value_name="syndrome_concordance",
        )
    )
    projected["soft_membership"] = np.nan
    projected["n_informative_traits"] = len(templates)
    projected["informative_traits"] = "|".join(templates)

    model_rows: list[dict[str, Any]] = []
    for index, template in enumerate(templates):
        model_rows.append(
            {
                "template": template,
                "source_mean": float(means[template]),
                "source_population_sd": float(scales[template]),
                "loading": float(loadings[index]),
                "source_factor_correlation": float(correlations[index]),
                "first_eigenvalue": first_value,
                "variance_fraction": first_value / float(np.sum(eigenvalues)),
            }
        )
    model = pd.DataFrame(model_rows)

    source_projected = component_wide.loc[source_species]
    orthogonality = {
        residual_name: float(
            np.cov(
                source_projected[SHARED_FACTOR].to_numpy(float),
                source_projected[residual_name].to_numpy(float),
                ddof=0,
            )[0, 1]
        )
        for residual_name in RESIDUAL_NAMES.values()
    }
    audit = {
        "n_complete_scored_species": len(wide),
        "n_complete_gift_source_species": len(source_species),
        "n_gift_entity_species_memberships": len(matched),
        "first_eigenvalue": first_value,
        "variance_fraction": first_value / float(np.sum(eigenvalues)),
        "sum_source_factor_correlations": float(np.sum(correlations)),
        "source_factor_residual_covariances": orthogonality,
        "factor_fit_used_island_status": False,
        "factor_fit_used_island_distance": False,
        "factor_fit_used_context_or_realm": False,
        "missing_template_imputation": False,
    }
    return projected, model, audit


def architecture_branching_config(
    branching_config: dict[str, Any], components: list[str]
) -> dict[str, Any]:
    """Create the frozen V4 multivariate family while retaining context layers."""

    config = deepcopy(branching_config)
    config["contract"] = "chapter1_v4_source_trained_architecture_decomposition_v1"
    config["hypothesis_name"] = "V4 shared pollination-architecture decomposition"
    config["branch_axes"] = {
        component: {"components": {component: 1.0}} for component in components
    }
    config["axis_sets"] = {
        "shared_architecture_decomposition": {
            "axes": components,
            "role": "secondary_discussion_concordance_decomposition",
            "classify": False,
        }
    }
    return config


def _run_source_adjusted_branching(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    architecture_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    outputs: dict[str, list[pd.DataFrame]] = {
        "branch_scores": [],
        "slopes": [],
        "within": [],
        "between": [],
    }
    for source_mode, subset in adjusted_scores.groupby("source_mode", sort=True):
        branch, slopes, within, between, _ = run_global_branching(
            subset,
            covariates,
            realm_assignment,
            pattern_config,
            architecture_config,
        )
        for key, frame in zip(
            outputs, (branch, slopes, within, between), strict=True
        ):
            frame = frame.copy()
            frame.insert(0, "source_mode", str(source_mode))
            outputs[key].append(frame)
    combined = {
        key: pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        for key, parts in outputs.items()
    }
    correction_specs = {
        "slopes": (
            ["context_layer", "axis_set", "stratum", "support_tier"],
            "p_value",
            "q_V4_across_axes_contexts_source_modes",
        ),
        "within": (
            ["context_layer", "axis_set", "stratum", "support_tier"],
            "p_value",
            "q_V4_across_contexts_source_modes",
        ),
        "between": (
            ["context_layer", "axis_set", "stratum", "support_tier"],
            "p_value",
            "q_V4_across_context_pairs_source_modes",
        ),
    }
    for key, (family, p_column, q_column) in correction_specs.items():
        frame = combined[key]
        if frame.empty:
            continue
        frame[q_column] = frame.groupby(family, group_keys=False)[p_column].transform(
            _bh
        )
        frame["V4_source_mode_family_supported"] = frame[q_column].le(0.05).fillna(
            False
        )
    return combined


def build_genus_decomposition(
    projected_species: pd.DataFrame,
    island_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    source_modes: list[str],
    strata: list[str],
    minimum_represented_genera: int,
) -> pd.DataFrame:
    """Separate source-matched genus entry, loading, and within-genus residual."""

    components = sorted(projected_species["syndrome"].astype(str).unique())
    wide = projected_species.pivot(
        index="accepted_species", columns="syndrome", values="syndrome_concordance"
    ).dropna(subset=components)
    functional = wide.reset_index()
    functional["accepted_species"] = functional["accepted_species"].astype(str)
    functional["genus"] = functional["accepted_species"].map(_genus)
    functional = functional.loc[functional["genus"].ne("")].copy()
    matched = match_gift_species(
        gift_flora, functional[["accepted_species"]]
    )
    source_species = matched.merge(
        functional, on="accepted_species", how="inner", validate="many_to_one"
    ).drop_duplicates("accepted_species")
    position_table = source_species.groupby("genus", as_index=True)[components].mean()
    genera = sorted(position_table.index.astype(str))
    genus_index = {genus: index for index, genus in enumerate(genera)}
    availability = broad_source_availability(gift_flora, set(genera))

    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: index for index, island in enumerate(islands)}
    entities = sorted(
        set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
        | set(pd.to_numeric(availability["entity_ID"], errors="coerce").dropna().astype(int))
    )
    entity_index = {entity: index for index, entity in enumerate(entities)}
    presence, richness = _availability_matrices(
        availability, entity_index, genus_index
    )
    actual_lookup = island_scores.set_index(
        ["stratum", "island_id", "syndrome"]
    )["syndrome_score"]
    species_filter = set(functional["accepted_species"].astype(str))

    parts: list[pd.DataFrame] = []
    for stratum in strata:
        counts = _representation_count_matrix(
            status_flora,
            island_index,
            genus_index,
            stratum=stratum,
            species_filter=species_filter,
        ).toarray()
        for source_mode in source_modes:
            assignment_matrix = _source_assignment_matrix(
                assignments,
                island_index,
                entity_index,
                source_mode=source_mode,
            )
            prevalence = (assignment_matrix @ presence).toarray().astype(np.int16)
            source_richness = (assignment_matrix @ richness).toarray().astype(np.float32)
            rows: list[dict[str, Any]] = []
            for component in components:
                positions = position_table.loc[genera, component].to_numpy(float)
                for island_position, island_id in enumerate(islands):
                    result = compute_island_enrichment(
                        prevalence[island_position],
                        source_richness[island_position],
                        counts[island_position],
                        positions,
                        matching="prevalence_richness",
                        minimum_represented_genera=minimum_represented_genera,
                    )
                    if result is None:
                        continue
                    actual_key = (stratum, island_id, component)
                    actual = (
                        float(actual_lookup.loc[actual_key])
                        if actual_key in actual_lookup.index
                        else float("nan")
                    )
                    result["beyond_genus_residual"] = (
                        actual - float(result["observed_species_mean"])
                        if math.isfinite(actual)
                        else float("nan")
                    )
                    result.update(
                        {
                            "island_id": island_id,
                            "stratum": stratum,
                            "source_mode": source_mode,
                            "source_matching": "prevalence_richness",
                            "minimum_represented_genera": minimum_represented_genera,
                            "architecture_component": component,
                        }
                    )
                    rows.append(result)
            if rows:
                parts.append(pd.DataFrame(rows))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _run_lineage_branching(
    lineage_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    architecture_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    outputs: dict[str, list[pd.DataFrame]] = {
        "slopes": [],
        "within": [],
        "between": [],
    }
    for (source_mode, outcome), subset in lineage_scores.groupby(
        ["source_mode", "lineage_outcome"], sort=True
    ):
        frame = subset.rename(
            columns={
                "architecture_component": "syndrome",
                "lineage_value": "syndrome_score",
                "n_represented_species": "n_species",
            }
        )
        _, slopes, within, between, _ = run_global_branching(
            frame,
            covariates,
            realm_assignment,
            pattern_config,
            architecture_config,
        )
        for key, table in zip(outputs, (slopes, within, between), strict=True):
            table = table.copy()
            table.insert(0, "source_mode", str(source_mode))
            table.insert(1, "lineage_outcome", str(outcome))
            outputs[key].append(table)
    combined = {
        key: pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        for key, parts in outputs.items()
    }
    correction_specs = {
        "slopes": (
            [
                "lineage_outcome",
                "context_layer",
                "stratum",
                "support_tier",
                "context",
            ],
            "p_value",
            "q_V4_across_axes_source_modes",
        ),
        "within": (
            [
                "lineage_outcome",
                "context_layer",
                "stratum",
                "support_tier",
                "context",
            ],
            "p_value",
            "q_V4_across_source_modes",
        ),
        "between": (
            [
                "lineage_outcome",
                "context_layer",
                "stratum",
                "support_tier",
                "context_a",
                "context_b",
            ],
            "p_value",
            "q_V4_across_source_modes",
        ),
    }
    for key, (family, p_column, q_column) in correction_specs.items():
        frame = combined[key]
        if frame.empty:
            continue
        frame[q_column] = frame.groupby(family, group_keys=False)[p_column].transform(
            _bh
        )
        frame["V4_source_mode_family_supported"] = frame[q_column].le(0.05).fillna(
            False
        )
    return combined


def run_v4(
    all_species_scores: pd.DataFrame,
    direct_species_scores: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    source_config: dict[str, Any],
    explanation_config: dict[str, Any],
) -> dict[str, Any]:
    spec = explanation_config["validations"][
        "V4_secondary_guild_architecture_decomposition"
    ]["frozen_implementation"]
    templates = tuple(str(x) for x in spec["sampled_templates"])
    minimum_source = int(spec["minimum_complete_source_species"])
    strata = [str(x) for x in spec["strata"]]
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]
    minimum_genera = int(spec["genus_decomposition"]["minimum_represented_genera"])
    scopes = {
        "all_analysis_eligible": all_species_scores,
        "direct_only": direct_species_scores,
    }
    all_outputs: dict[str, Any] = {}
    scope_manifests: list[dict[str, Any]] = []
    taxonomic_rows: list[dict[str, Any]] = []
    for evidence_scope, scores in scopes.items():
        projected, factor_model, factor_audit = fit_and_project_source_factor(
            scores,
            gift_flora,
            templates=templates,
            minimum_complete_source_species=minimum_source,
        )
        components = [SHARED_FACTOR, *RESIDUAL_NAMES.values()]
        architecture_config = architecture_branching_config(
            branching_config, components
        )
        island_scores = build_island_syndrome_scores(
            status_flora, projected, strata
        )
        branch, slopes, within, between, _ = run_global_branching(
            island_scores,
            covariates,
            realm_assignment,
            pattern_config,
            architecture_config,
        )

        matched_scores, match_audit = match_gift_species_to_frozen_scores(
            gift_flora, projected
        )
        source_scores = build_mainland_syndrome_scores(
            matched_scores,
            min_species=int(
                spec["source_adjustment"][
                    "minimum_scored_species_per_source_component"
                ]
            ),
        )
        expectations = build_source_expectations(
            assignments,
            source_scores,
            min_source_regions=int(
                spec["source_adjustment"]["minimum_source_regions"]
            ),
        )
        adjusted = build_adjusted_island_scores(
            island_scores, expectations, strata
        )
        adjusted_results = _run_source_adjusted_branching(
            adjusted,
            covariates,
            realm_assignment,
            pattern_config,
            architecture_config,
        )

        lineage = build_genus_decomposition(
            projected,
            island_scores,
            status_flora,
            gift_flora,
            assignments,
            covariates,
            source_modes=source_modes,
            strata=strata,
            minimum_represented_genera=minimum_genera,
        )
        lineage_long = lineage.melt(
            id_vars=[
                "island_id",
                "stratum",
                "source_mode",
                "source_matching",
                "minimum_represented_genera",
                "architecture_component",
                "n_represented_genera",
                "n_represented_species",
                "n_candidate_genera",
            ],
            value_vars=list(LINEAGE_OUTCOMES),
            var_name="lineage_outcome",
            value_name="lineage_value",
        ).dropna(subset=["lineage_value"])
        lineage_results = _run_lineage_branching(
            lineage_long,
            covariates,
            realm_assignment,
            pattern_config,
            architecture_config,
        )

        for component in components:
            taxonomic_rows.extend(
                [
                    {
                        "evidence_scope": evidence_scope,
                        "architecture_component": component,
                        "taxonomic_depth": "family",
                        "status": "not_evaluable_pending_V2_taxonomy_contract",
                        "reason": "formal_workflow_has_no_frozen_comprehensive_species_to_family_mapping",
                    },
                    {
                        "evidence_scope": evidence_scope,
                        "architecture_component": component,
                        "taxonomic_depth": "genus",
                        "status": "evaluated",
                        "reason": "source_matched_entry_loading_and_beyond_genus_residual",
                    },
                ]
            )
        scope_manifests.append(
            {
                "evidence_scope": evidence_scope,
                **factor_audit,
                "n_projected_species_component_rows": len(projected),
                "n_island_component_rows": len(island_scores),
                "n_source_adjusted_rows": len(adjusted),
                "n_genus_decomposition_rows": len(lineage),
                "gift_match_audit": match_audit.to_dict("records"),
            }
        )
        all_outputs[evidence_scope] = {
            "projected_species": projected,
            "factor_model": factor_model,
            "island_scores": island_scores,
            "observed_branch_scores": branch,
            "observed_slopes": slopes,
            "observed_within": within,
            "observed_between": between,
            "source_scores": source_scores,
            "source_expectations": expectations,
            "source_adjusted_scores": adjusted,
            "source_adjusted": adjusted_results,
            "genus_scores": lineage,
            "genus_long": lineage_long,
            "genus_results": lineage_results,
        }

    manifest = {
        "contract": "chapter1_v4_source_trained_architecture_decomposition_v1",
        "factor_training_population": "frozen_GIFT_mainland_native_matched_species_only",
        "templates": list(templates),
        "components": [SHARED_FACTOR, *RESIDUAL_NAMES.values()],
        "complete_case_no_imputation": True,
        "source_pool_adjusted": True,
        "genus_entry_loading_and_beyond_genus_residual": True,
        "family_decomposition": "not_evaluable_pending_V2_taxonomy_contract",
        "pollinator_identity_observed": False,
        "effective_pollination_service_observed": False,
        "functional_replacement_identified": False,
        "template_residual_is_pollinator_identity": False,
        "scope_manifests": scope_manifests,
        "claim_boundary": (
            "Shared or residual plant-trait architecture only. Neither a guild-labelled "
            "residual nor persistence after source/genus adjustment identifies pollinator "
            "presence, loss, effectiveness, or functional replacement."
        ),
    }
    return {
        "scopes": all_outputs,
        "taxonomic_status": pd.DataFrame(taxonomic_rows),
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-species-scores-csv", type=Path, required=True)
    parser.add_argument("--direct-species-scores-csv", type=Path, required=True)
    parser.add_argument("--status-flora-csv", type=Path, required=True)
    parser.add_argument("--gift-flora-csv", type=Path, required=True)
    parser.add_argument("--source-assignments-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--branching-config-path", type=Path, required=True)
    parser.add_argument("--source-config-path", type=Path, required=True)
    parser.add_argument("--explanation-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    outputs = run_v4(
        pd.read_csv(args.all_species_scores_csv),
        pd.read_csv(args.direct_species_scores_csv),
        pd.read_csv(args.status_flora_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.branching_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.source_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.explanation_config_path.read_text(encoding="utf-8")),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for scope, tables in outputs["scopes"].items():
        root = args.output_dir / scope
        root.mkdir(parents=True, exist_ok=True)
        for name, table in tables.items():
            if isinstance(table, pd.DataFrame):
                suffix = ".csv.gz" if len(table) > 100_000 else ".csv"
                table.to_csv(root / f"{name}{suffix}", index=False)
            elif isinstance(table, dict):
                for child_name, child_table in table.items():
                    child_table.to_csv(
                        root / f"{name}_{child_name}.csv", index=False
                    )
    outputs["taxonomic_status"].to_csv(
        args.output_dir / "V4_taxonomic_depth_status.csv", index=False
    )
    (args.output_dir / "chapter1_v4_architecture_manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outputs["manifest"], indent=2))


if __name__ == "__main__":
    main()
