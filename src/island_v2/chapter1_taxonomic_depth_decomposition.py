"""Prospective V2 family-to-genus depth decomposition for Chapter 1 H3.

Family and genus are grouping variables only.  They never fill a missing trait.
Source-derived group positions are estimated from scored species matched to the frozen
GIFT mainland flora, while group availability is built outcome-blind from all GIFT
species that match the frozen taxonomy table.  Family and genus residuals are evaluated
on exactly the same observed island species.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from island_v2.chapter1_global_branching import run_global_branching
from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _source_assignment_matrix,
    match_gift_species,
)
from island_v2.chapter1_pr138_syndrome_analysis import _bh

STAGES = ("observed_score", "after_family_residual", "after_genus_residual")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_taxonomy(
    taxonomy: pd.DataFrame, expected_sha256: str, observed_sha256: str
) -> pd.DataFrame:
    required = {"accepted_species", "genus", "family"}
    if missing := required - set(taxonomy.columns):
        raise ValueError(f"taxonomy input missing columns: {sorted(missing)}")
    if observed_sha256 != expected_sha256:
        raise ValueError(
            "taxonomy SHA-256 differs from the frozen V2 contract: "
            f"{observed_sha256} != {expected_sha256}"
        )
    work = taxonomy[["accepted_species", "genus", "family"]].copy().fillna("")
    for column in work:
        work[column] = work[column].astype(str).str.strip()
    work = work.loc[work["accepted_species"].ne("")].copy()
    if work["accepted_species"].duplicated().any():
        raise ValueError("taxonomy input contains duplicate accepted_species")
    return work


def _validate_scores(
    species_scores: pd.DataFrame, axes: list[str]
) -> pd.DataFrame:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    if missing := required - set(species_scores.columns):
        raise ValueError(f"species scores missing columns: {sorted(missing)}")
    work = species_scores.loc[
        species_scores["syndrome"].astype(str).isin(axes),
        ["accepted_species", "syndrome", "syndrome_concordance"],
    ].copy()
    work["accepted_species"] = work["accepted_species"].fillna("").astype(str)
    work["syndrome"] = work["syndrome"].astype(str)
    work["syndrome_concordance"] = pd.to_numeric(
        work["syndrome_concordance"], errors="coerce"
    )
    work = work.dropna(subset=["syndrome_concordance"])
    duplicated = work.duplicated(["accepted_species", "syndrome"], keep=False)
    if duplicated.any():
        raise ValueError("duplicate species response scores in V2 input")
    missing_axes = set(axes) - set(work["syndrome"])
    if missing_axes:
        raise ValueError(f"V2 response axes absent: {sorted(missing_axes)}")
    return work


def build_source_group_contract(
    species_scores: pd.DataFrame,
    taxonomy: pd.DataFrame,
    gift_flora: pd.DataFrame,
    *,
    axes: list[str],
    minimum_source_scored_species: int,
    matched_taxa: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame, dict[str, Any]]:
    """Build outcome-blind availability and scored family/genus positions."""

    scores = _validate_scores(species_scores, axes)
    matched = (
        matched_taxa.copy()
        if matched_taxa is not None
        else match_gift_species(gift_flora, taxonomy[["accepted_species"]])
    )
    matched_identity = matched.merge(
        taxonomy, on="accepted_species", how="left", validate="many_to_one"
    )
    source_species = (
        matched_identity[["accepted_species", "genus", "family"]]
        .drop_duplicates("accepted_species")
        .merge(scores, on="accepted_species", how="inner", validate="one_to_many")
    )
    position_parts: list[pd.DataFrame] = []
    availability: dict[str, pd.DataFrame] = {}
    for level in ("family", "genus"):
        eligible_source = source_species.loc[source_species[level].ne("")].copy()
        positions = (
            eligible_source.groupby([level, "syndrome"], as_index=False)
            .agg(
                source_group_position=("syndrome_concordance", "mean"),
                n_source_scored_species=("accepted_species", "nunique"),
            )
        )
        positions = positions.loc[
            positions["n_source_scored_species"].ge(
                int(minimum_source_scored_species)
            )
        ].copy()
        positions = positions.rename(columns={level: "taxon_group"})
        positions["taxonomic_level"] = level
        position_parts.append(positions)

        eligible_groups = set(positions["taxon_group"].astype(str))
        group_availability = matched_identity.loc[
            matched_identity[level].astype(str).isin(eligible_groups),
            ["entity_ID", level],
        ].drop_duplicates()
        availability[level] = group_availability.rename(
            columns={level: "taxon_group"}
        )
    position_table = pd.concat(position_parts, ignore_index=True)
    audit = {
        "n_taxonomy_species": len(taxonomy),
        "n_gift_taxonomy_memberships": len(matched),
        "n_unique_gift_taxonomy_species": int(
            matched["accepted_species"].astype(str).nunique()
        ),
        "n_scored_source_species": int(
            source_species["accepted_species"].astype(str).nunique()
        ),
        "minimum_source_scored_species_per_group": int(
            minimum_source_scored_species
        ),
        "family_is_trait_inference": False,
        "availability_uses_trait_scores": False,
        "source_assignment_uses_island_outcomes": False,
    }
    return position_table, availability, matched, audit


def _availability_matrix(
    availability: pd.DataFrame,
    entity_index: dict[int, int],
    group_index: dict[str, int],
) -> sparse.csr_matrix:
    work = availability.copy()
    work["entity_ID"] = pd.to_numeric(work["entity_ID"], errors="coerce")
    work = work.dropna(subset=["entity_ID"])
    work["entity_ID"] = work["entity_ID"].astype(int)
    work["taxon_group"] = work["taxon_group"].fillna("").astype(str)
    work = work.loc[
        work["entity_ID"].isin(entity_index)
        & work["taxon_group"].isin(group_index)
    ].drop_duplicates(["entity_ID", "taxon_group"])
    rows = work["entity_ID"].map(entity_index).to_numpy(int)
    columns = work["taxon_group"].map(group_index).to_numpy(int)
    return sparse.csr_matrix(
        (np.ones(len(work), dtype=np.int8), (rows, columns)),
        shape=(len(entity_index), len(group_index)),
    )


def _stratum_mask(frame: pd.DataFrame, stratum: str) -> pd.Series:
    if stratum == "all_native":
        return frame["origin_status"].astype(str).eq("native")
    if stratum == "native_nonendemic":
        return frame["floristic_status"].astype(str).eq("native_nonendemic")
    raise ValueError(f"unsupported V2 stratum: {stratum}")


def build_taxonomic_decomposition(
    species_scores: pd.DataFrame,
    taxonomy: pd.DataFrame,
    status_flora: pd.DataFrame,
    positions: pd.DataFrame,
    availability: dict[str, pd.DataFrame],
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    axes: list[str],
    source_modes: list[str],
    strata: list[str],
    minimum_species: int,
    minimum_families: int,
    minimum_genera: int,
) -> pd.DataFrame:
    scores = _validate_scores(species_scores, axes)
    required_status = {
        "island_id",
        "accepted_species",
        "origin_status",
        "floristic_status",
    }
    if missing := required_status - set(status_flora.columns):
        raise ValueError(f"status flora missing V2 columns: {sorted(missing)}")

    position_lookup: dict[str, pd.Series] = {}
    group_indices: dict[str, dict[str, int]] = {}
    presence: dict[str, sparse.csr_matrix] = {}
    entities = sorted(
        set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
        | {
            int(value)
            for table in availability.values()
            for value in pd.to_numeric(table["entity_ID"], errors="coerce")
            .dropna()
            .astype(int)
        }
    )
    entity_index = {entity: index for index, entity in enumerate(entities)}
    for level in ("family", "genus"):
        table = positions.loc[positions["taxonomic_level"].eq(level)].copy()
        groups = sorted(table["taxon_group"].astype(str).unique())
        group_indices[level] = {group: index for index, group in enumerate(groups)}
        presence[level] = _availability_matrix(
            availability[level], entity_index, group_indices[level]
        )
        position_lookup[level] = table.set_index(
            ["taxon_group", "syndrome"]
        )["source_group_position"]

    islands = sorted(covariates["island_id"].astype(str).unique())
    island_index = {island: index for index, island in enumerate(islands)}
    flora = status_flora[list(required_status)].drop_duplicates(
        ["island_id", "accepted_species"]
    )
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.merge(taxonomy, on="accepted_species", how="left", validate="many_to_one")
    flora = flora.merge(scores, on="accepted_species", how="inner", validate="many_to_many")
    flora = flora.loc[flora["island_id"].isin(island_index)].copy()
    flora["island_position"] = flora["island_id"].map(island_index).astype(int)
    for level in ("family", "genus"):
        flora[f"{level}_position"] = [
            position_lookup[level].get((str(group), str(axis)), np.nan)
            for group, axis in zip(flora[level], flora["syndrome"], strict=True)
        ]
        flora[f"{level}_index"] = flora[level].map(group_indices[level])

    parts: list[pd.DataFrame] = []
    for source_mode in source_modes:
        assignment = _source_assignment_matrix(
            assignments,
            island_index,
            entity_index,
            source_mode=source_mode,
        )
        available_by_level = {
            level: (assignment @ presence[level]).tocsr()
            for level in ("family", "genus")
        }
        for stratum in strata:
            work = flora.loc[_stratum_mask(flora, stratum)].copy()
            complete_group = work["family_index"].notna() & work["genus_index"].notna()
            work = work.loc[complete_group].copy()
            work["family_index"] = work["family_index"].astype(int)
            work["genus_index"] = work["genus_index"].astype(int)
            family_present = np.asarray(
                available_by_level["family"][
                    work["island_position"].to_numpy(int),
                    work["family_index"].to_numpy(int),
                ]
            ).reshape(-1)
            genus_present = np.asarray(
                available_by_level["genus"][
                    work["island_position"].to_numpy(int),
                    work["genus_index"].to_numpy(int),
                ]
            ).reshape(-1)
            work = work.loc[(family_present > 0) & (genus_present > 0)].copy()
            if work.empty:
                continue
            summary = (
                work.groupby(["island_id", "syndrome"], as_index=False)
                .agg(
                    observed_score=("syndrome_concordance", "mean"),
                    family_expected=("family_position", "mean"),
                    genus_expected=("genus_position", "mean"),
                    n_species=("accepted_species", "nunique"),
                    n_families=("family", "nunique"),
                    n_genera=("genus", "nunique"),
                )
            )
            summary = summary.loc[
                summary["n_species"].ge(int(minimum_species))
                & summary["n_families"].ge(int(minimum_families))
                & summary["n_genera"].ge(int(minimum_genera))
            ].copy()
            summary["after_family_residual"] = (
                summary["observed_score"] - summary["family_expected"]
            )
            summary["family_to_genus_increment"] = (
                summary["genus_expected"] - summary["family_expected"]
            )
            summary["after_genus_residual"] = (
                summary["observed_score"] - summary["genus_expected"]
            )
            summary["source_mode"] = source_mode
            summary["stratum"] = stratum
            parts.append(summary)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def taxonomic_branching_config(
    branching_config: dict[str, Any], axes: list[str]
) -> dict[str, Any]:
    config = deepcopy(branching_config)
    config["contract"] = "chapter1_v2_source_matched_taxonomic_depth_v1"
    config["hypothesis_name"] = "V2 family-to-genus taxonomic-depth decomposition"
    branch_axes: dict[str, Any] = {}
    axis_sets: dict[str, Any] = {}
    for stage in STAGES:
        names = []
        for axis in axes:
            name = f"{stage}__{axis}"
            names.append(name)
            branch_axes[name] = {"components": {name: 1.0}}
        axis_sets[f"taxonomic_stage__{stage}"] = {
            "axes": names,
            "role": "prospective_H3_taxonomic_depth_falsification",
            "classify": False,
        }
    config["branch_axes"] = branch_axes
    config["axis_sets"] = axis_sets
    return config


def run_taxonomic_models(
    decomposition: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    taxonomic_config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    long = decomposition.melt(
        id_vars=[
            "island_id",
            "syndrome",
            "stratum",
            "source_mode",
            "n_species",
        ],
        value_vars=list(STAGES),
        var_name="taxonomic_stage",
        value_name="syndrome_score",
    )
    long["syndrome"] = long["taxonomic_stage"] + "__" + long["syndrome"]
    outputs: dict[str, list[pd.DataFrame]] = {
        "slopes": [],
        "within": [],
        "between": [],
    }
    for source_mode, subset in long.groupby("source_mode", sort=True):
        _, slopes, within, between, _ = run_global_branching(
            subset,
            covariates,
            realm_assignment,
            pattern_config,
            taxonomic_config,
        )
        for key, table in zip(outputs, (slopes, within, between), strict=True):
            table = table.copy()
            table.insert(0, "source_mode", str(source_mode))
            outputs[key].append(table)
    combined = {
        key: pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        for key, parts in outputs.items()
    }
    specs = {
        "slopes": (
            ["context_layer", "axis_set", "stratum", "support_tier"],
            "q_V2_across_axes_contexts_source_modes",
        ),
        "within": (
            ["context_layer", "axis_set", "stratum", "support_tier"],
            "q_V2_across_contexts_source_modes",
        ),
        "between": (
            ["context_layer", "axis_set", "stratum", "support_tier"],
            "q_V2_across_context_pairs_source_modes",
        ),
    }
    for key, (family, q_column) in specs.items():
        frame = combined[key]
        if frame.empty:
            continue
        frame[q_column] = frame.groupby(family, group_keys=False)["p_value"].transform(
            _bh
        )
        frame["V2_source_mode_family_supported"] = frame[q_column].le(0.05).fillna(
            False
        )
    combined["long_scores"] = long
    return combined


def classify_depth(
    within: pd.DataFrame, source_modes: list[str], strata: list[str]
) -> pd.DataFrame:
    confirmatory = within.loc[
        within["support_tier"].eq("confirmatory")
        & within["stratum"].astype(str).isin(strata)
    ].copy()
    keys = ["context_layer", "stratum", "context"]
    rows: list[dict[str, Any]] = []
    expected_modes = len(source_modes)
    for key, group in confirmatory.groupby(keys, sort=True):
        meta = dict(zip(keys, key, strict=True))
        support: dict[str, int] = {}
        for stage in STAGES:
            axis_set = f"taxonomic_stage__{stage}"
            subset = group.loc[group["axis_set"].eq(axis_set)]
            support[stage] = int(
                subset.loc[
                    subset["V2_source_mode_family_supported"].astype(bool),
                    "source_mode",
                ].nunique()
            )
        if support["observed_score"] != expected_modes:
            classification = "observed_vector_not_robust_across_source_modes"
        elif support["after_family_residual"] == 0:
            classification = "compatible_with_deep_family_sorting"
        elif support["after_family_residual"] == expected_modes and support[
            "after_genus_residual"
        ] == 0:
            classification = "compatible_with_genus_level_assembly_beyond_family"
        elif support["after_genus_residual"] == expected_modes:
            classification = "within_genus_or_non_taxonomic_residual_retained"
        else:
            classification = "source_definition_sensitive_taxonomic_depth"
        rows.append(
            {
                **meta,
                "n_source_modes_expected": expected_modes,
                "observed_supported_modes": support["observed_score"],
                "after_family_supported_modes": support["after_family_residual"],
                "after_genus_supported_modes": support["after_genus_residual"],
                "V2_taxonomic_depth_classification": classification,
            }
        )
    return pd.DataFrame(rows)


def run_v2(
    all_species_scores: pd.DataFrame,
    direct_species_scores: pd.DataFrame,
    taxonomy: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    source_config: dict[str, Any],
    explanation_config: dict[str, Any],
    *,
    taxonomy_sha256: str,
) -> dict[str, Any]:
    spec = explanation_config["validations"]["V2_H3_taxonomic_depth"][
        "frozen_implementation"
    ]
    taxonomy = validate_taxonomy(
        taxonomy, str(spec["taxonomy_input"]["sha256"]), taxonomy_sha256
    )
    axes = [str(x) for x in spec["plant_response_axes"]]
    source_modes = [str(x) for x in source_config["source_assignment"]["primary_modes"]]
    strata = [str(x) for x in spec["strata"]]
    minimum = spec["minimum_support"]
    matched_taxa = match_gift_species(gift_flora, taxonomy[["accepted_species"]])
    taxonomic_config = taxonomic_branching_config(branching_config, axes)
    scopes = {
        "all_analysis_eligible": all_species_scores,
        "direct_only": direct_species_scores,
    }
    scope_outputs: dict[str, Any] = {}
    audits: list[dict[str, Any]] = []
    for evidence_scope, scores in scopes.items():
        positions, availability, _, audit = build_source_group_contract(
            scores,
            taxonomy,
            gift_flora,
            axes=axes,
            minimum_source_scored_species=2,
            matched_taxa=matched_taxa,
        )
        decomposition = build_taxonomic_decomposition(
            scores,
            taxonomy,
            status_flora,
            positions,
            availability,
            assignments,
            covariates,
            axes=axes,
            source_modes=source_modes,
            strata=strata,
            minimum_species=int(minimum["observed_species_per_island_response"]),
            minimum_families=int(minimum["represented_families"]),
            minimum_genera=int(minimum["represented_genera"]),
        )
        models = run_taxonomic_models(
            decomposition,
            covariates,
            realm_assignment,
            pattern_config,
            taxonomic_config,
        )
        classification = classify_depth(models["within"], source_modes, strata)
        audits.append(
            {
                "evidence_scope": evidence_scope,
                **audit,
                "n_source_position_rows": len(positions),
                "n_decomposition_rows": len(decomposition),
                "n_depth_classification_rows": len(classification),
            }
        )
        scope_outputs[evidence_scope] = {
            "source_group_positions": positions,
            "decomposition": decomposition,
            **models,
            "classification": classification,
        }
    manifest = {
        "contract": "chapter1_v2_source_matched_taxonomic_depth_v1",
        "taxonomy_sha256": taxonomy_sha256,
        "n_taxonomy_species": len(taxonomy),
        "n_family_resolved_taxonomy_species": int(taxonomy["family"].ne("").sum()),
        "plant_response_axes": axes,
        "formal_vector_stages": list(STAGES),
        "family_used_for_trait_inference": False,
        "family_and_genus_nulls_use_common_observed_species": True,
        "source_assignment_uses_island_outcomes": False,
        "pollinator_mechanism_identified": False,
        "scope_audits": audits,
        "claim_boundary": (
            "Taxonomic assembly depth only. Residual depth does not identify historical "
            "colonisation, in-situ evolution, pollinator loss, or effective service."
        ),
    }
    return {"scopes": scope_outputs, "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-species-scores-csv", type=Path, required=True)
    parser.add_argument("--direct-species-scores-csv", type=Path, required=True)
    parser.add_argument("--taxonomy-csv", type=Path, required=True)
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

    taxonomy_hash = file_sha256(args.taxonomy_csv)
    outputs = run_v2(
        pd.read_csv(args.all_species_scores_csv),
        pd.read_csv(args.direct_species_scores_csv),
        pd.read_csv(args.taxonomy_csv),
        pd.read_csv(args.status_flora_csv),
        pd.read_csv(args.gift_flora_csv),
        pd.read_csv(args.source_assignments_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.branching_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.source_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.explanation_config_path.read_text(encoding="utf-8")),
        taxonomy_sha256=taxonomy_hash,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for scope, tables in outputs["scopes"].items():
        root = args.output_dir / scope
        root.mkdir(parents=True, exist_ok=True)
        for name, table in tables.items():
            suffix = ".csv.gz" if len(table) > 100_000 else ".csv"
            table.to_csv(root / f"{name}{suffix}", index=False)
    (args.output_dir / "chapter1_v2_taxonomic_depth_manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outputs["manifest"], indent=2))


if __name__ == "__main__":
    main()
