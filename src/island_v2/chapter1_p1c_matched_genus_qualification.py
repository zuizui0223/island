"""Qualification gate for the P1c matched-complexity genus null.

This module MUST NOT report a pseudo-genus attenuation statistic.  It verifies only
that (1) a lightweight Palearctic refit reproduces the frozen true-genus slopes,
(2) one deterministic qualification partition preserves within-family genus
complexity exactly, and (3) all eight direct-only primary profiles remain fit-able.
"""
from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_lineage_representation_bridge import match_gift_species
from island_v2.chapter1_pr138_syndrome_analysis import _prepare, _within_context
from island_v2.chapter1_taxonomic_depth_decomposition import (
    build_source_group_contract,
    build_taxonomic_decomposition,
    file_sha256,
    validate_taxonomy,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


class P1CQualificationError(ValueError):
    """Raised when a P1c qualification invariant is violated."""


def load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise P1CQualificationError(f"expected YAML mapping: {path}")
    return value


def make_matched_pseudo_taxonomy(
    taxonomy: pd.DataFrame,
    *,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Randomize species into pseudo-genera within families preserving exact group sizes."""

    required = {"accepted_species", "family", "genus"}
    missing = required - set(taxonomy.columns)
    if missing:
        raise P1CQualificationError(f"taxonomy missing columns: {sorted(missing)}")
    work = taxonomy.copy()
    for column in required:
        work[column] = work[column].fillna("").astype(str).str.strip()
    if work["accepted_species"].duplicated().any():
        raise P1CQualificationError("duplicate accepted_species in taxonomy")

    rng = np.random.default_rng(int(seed))
    pseudo = pd.Series("", index=work.index, dtype=object)
    n_families = 0
    n_groups = 0
    size_mismatch_families: list[str] = []

    for family, family_frame in work.loc[
        work["family"].ne("") & work["genus"].ne("")
    ].groupby("family", sort=True):
        n_families += 1
        original_counts = (
            family_frame.groupby("genus")["accepted_species"].size().sort_index()
        )
        sizes = original_counts.to_numpy(int)
        indices = family_frame.index.to_numpy(copy=True)
        rng.shuffle(indices)
        offset = 0
        generated_counts: list[int] = []
        for group_number, size in enumerate(sizes):
            chosen = indices[offset : offset + int(size)]
            label = f"PSEUDO::{family}::{group_number:05d}"
            pseudo.loc[chosen] = label
            generated_counts.append(len(chosen))
            offset += int(size)
            n_groups += 1
        if offset != len(indices):
            raise P1CQualificationError(f"pseudo partition did not consume family {family}")
        if sorted(generated_counts) != sorted(sizes.tolist()):
            size_mismatch_families.append(str(family))

    out = work.copy()
    eligible = out["family"].ne("") & out["genus"].ne("")
    out.loc[eligible, "genus"] = pseudo.loc[eligible]
    if out.loc[eligible, "genus"].eq("").any():
        raise P1CQualificationError("eligible species missing pseudo-genus assignment")

    audit = {
        "seed": int(seed),
        "n_families_randomized": int(n_families),
        "n_pseudo_genera": int(n_groups),
        "n_size_mismatch_families": int(len(size_mismatch_families)),
        "group_size_multiset_preserved_all_families": not size_mismatch_families,
    }
    return out, audit


def fit_palearctic_profiles(
    decomposition: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    *,
    stage: str,
    source_modes: list[str],
    strata: list[str],
    axes: list[str],
) -> pd.DataFrame:
    """Fit only the frozen direct-only Palearctic profiles using the primary estimator."""

    required = {"island_id", "syndrome", "stratum", "source_mode", "n_species", stage}
    missing = required - set(decomposition.columns)
    if missing:
        raise P1CQualificationError(f"decomposition missing columns: {sorted(missing)}")

    realm_column = "biogeographic_realm"
    cov = covariates.copy()
    if realm_column not in cov.columns:
        realm_required = {"island_id", realm_column}
        realm_missing = realm_required - set(realm_assignment.columns)
        if realm_missing:
            raise P1CQualificationError(f"realm assignment missing: {sorted(realm_missing)}")
        cov = cov.merge(
            realm_assignment[["island_id", realm_column]].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="one_to_one",
        )

    layer_pattern = deepcopy(pattern_config)
    layer_pattern["context_column"] = realm_column
    layer_pattern["contexts"] = ["Palearctic"]
    threshold = int(pattern_config["support_tiers"]["confirmatory"])

    rows: list[pd.DataFrame] = []
    for source_mode in source_modes:
        subset = decomposition.loc[
            decomposition["source_mode"].astype(str).eq(source_mode)
            & decomposition["syndrome"].astype(str).isin(axes)
        ].copy()
        scores = subset[
            ["island_id", "stratum", "syndrome", stage, "n_species"]
        ].rename(columns={stage: "syndrome_score"})
        data = _prepare(scores, cov, layer_pattern)
        for stratum in strata:
            slopes, result = _within_context(
                data,
                stratum=stratum,
                context_value="Palearctic",
                support_tier="confirmatory",
                threshold=threshold,
                pattern_config=layer_pattern,
                syndrome_config={},
            )
            if str(result.get("status", "")) != "fit":
                continue
            fitted_axes = set(slopes["syndrome"].astype(str)) if not slopes.empty else set()
            if fitted_axes != set(axes):
                continue
            slopes = slopes.copy()
            slopes.insert(0, "source_mode", source_mode)
            slopes.insert(1, "taxonomic_stage", stage)
            rows.append(slopes)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def compare_to_frozen_slopes(
    reproduced: pd.DataFrame,
    frozen_slopes: pd.DataFrame,
    *,
    stage: str,
    source_modes: list[str],
    strata: list[str],
    axes: list[str],
) -> dict[str, Any]:
    frozen = frozen_slopes.loc[
        frozen_slopes["context_layer"].astype(str).eq("biogeographic_realm")
        & frozen_slopes["axis_set"].astype(str).eq(f"taxonomic_stage__{stage}")
        & frozen_slopes["stratum"].astype(str).isin(strata)
        & frozen_slopes["support_tier"].astype(str).eq("confirmatory")
        & frozen_slopes["context"].astype(str).eq("Palearctic")
        & frozen_slopes["source_mode"].astype(str).isin(source_modes)
    ].copy()
    frozen["raw_axis"] = frozen["syndrome"].astype(str).str.removeprefix(f"{stage}__")
    frozen = frozen.loc[frozen["raw_axis"].isin(axes)].copy()

    rep = reproduced.copy()
    keys_rep = ["source_mode", "stratum", "syndrome"]
    keys_frozen = ["source_mode", "stratum", "raw_axis"]
    if len(rep) != len(source_modes) * len(strata) * len(axes):
        raise P1CQualificationError("true lightweight fit did not produce all expected slope rows")
    if len(frozen) != len(rep):
        raise P1CQualificationError("frozen slope table does not contain expected focal rows")

    joined = rep.merge(
        frozen,
        left_on=keys_rep,
        right_on=keys_frozen,
        how="inner",
        suffixes=("_reproduced", "_frozen"),
        validate="one_to_one",
    )
    if len(joined) != len(rep):
        raise P1CQualificationError("true reproduced and frozen slope keys do not match")
    difference = np.abs(
        pd.to_numeric(joined["distance_slope_reproduced"])
        - pd.to_numeric(joined["distance_slope_frozen"])
    )
    return {
        "stage": stage,
        "n_compared_slope_rows": int(len(joined)),
        "max_absolute_slope_difference": float(difference.max()),
    }


def run_qualification(
    *,
    artifact_root: Path,
    taxonomy_path: Path,
    p1_config_path: Path,
    pattern_config_path: Path,
    source_config_path: Path,
    explanation_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    p1 = load_yaml(p1_config_path)
    pattern = load_yaml(pattern_config_path)
    source_config = load_yaml(source_config_path)
    explanation = load_yaml(explanation_config_path)
    if p1.get("contract") != "chapter1_p1_assembly_depth_defense_v1":
        raise P1CQualificationError("unexpected P1 contract")

    v2_spec = explanation["validations"]["V2_H3_taxonomic_depth"]["frozen_implementation"]
    taxonomy_hash = file_sha256(taxonomy_path)
    taxonomy = validate_taxonomy(
        pd.read_csv(taxonomy_path),
        str(v2_spec["taxonomy_input"]["sha256_newline_canonicalized"]),
        taxonomy_hash,
    )

    direct_scores = pd.read_csv(
        artifact_root / "syndrome/direct/species_syndrome_concordance.csv.gz"
    )
    status_flora = pd.read_csv(
        artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz"
    )
    gift_flora = pd.read_csv(
        artifact_root / "fixed/source/source/gift_native_mainland_flora.csv.gz"
    )
    assignments = pd.read_csv(
        artifact_root / "fixed/source/result/island_source_assignments.csv.gz"
    )
    covariates = pd.read_csv(
        artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv"
    )
    realm_assignment = pd.read_csv(
        artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv"
    )
    true_decomposition = pd.read_csv(
        artifact_root / "taxonomic-depth/direct_only/decomposition.csv"
    )
    frozen_slopes = pd.read_csv(
        artifact_root / "taxonomic-depth/direct_only/slopes.csv"
    )

    primary = p1["primary_context"]
    axes = [str(x) for x in primary["response_axes"]]
    source_modes = [str(x) for x in primary["source_modes"]]
    strata = [str(x) for x in primary["strata"]]
    qualification = p1["p1c_matched_complexity_genus_null"]["qualification_gate"]
    tolerance = float(qualification["slope_reproduction_absolute_tolerance"])

    true_checks = []
    for stage in ("after_family_residual", "after_genus_residual"):
        reproduced = fit_palearctic_profiles(
            true_decomposition,
            covariates,
            realm_assignment,
            pattern,
            stage=stage,
            source_modes=source_modes,
            strata=strata,
            axes=axes,
        )
        true_checks.append(
            compare_to_frozen_slopes(
                reproduced,
                frozen_slopes,
                stage=stage,
                source_modes=source_modes,
                strata=strata,
                axes=axes,
            )
        )
    true_reproduction_pass = all(
        check["max_absolute_slope_difference"] <= tolerance for check in true_checks
    )

    pseudo_taxonomy, partition_audit = make_matched_pseudo_taxonomy(
        taxonomy,
        seed=int(qualification["qualification_seed"]),
    )
    matched_taxa = match_gift_species(gift_flora, taxonomy[["accepted_species"]])
    positions, availability, _, _ = build_source_group_contract(
        direct_scores,
        pseudo_taxonomy,
        gift_flora,
        axes=axes,
        minimum_source_scored_species=2,
        matched_taxa=matched_taxa,
    )
    minimum = v2_spec["minimum_support"]
    pseudo_decomposition = build_taxonomic_decomposition(
        direct_scores,
        pseudo_taxonomy,
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
    pseudo_slopes = fit_palearctic_profiles(
        pseudo_decomposition,
        covariates,
        realm_assignment,
        pattern,
        stage="after_genus_residual",
        source_modes=source_modes,
        strata=strata,
        axes=axes,
    )
    expected_rows = len(source_modes) * len(strata) * len(axes)
    pseudo_fit_pass = len(pseudo_slopes) == expected_rows
    n_profiles_fit = 0
    if not pseudo_slopes.empty:
        n_profiles_fit = int(
            pseudo_slopes[["source_mode", "stratum"]].drop_duplicates().shape[0]
        )

    qualified = bool(
        true_reproduction_pass
        and partition_audit["group_size_multiset_preserved_all_families"]
        and pseudo_fit_pass
        and n_profiles_fit == len(source_modes) * len(strata)
    )
    manifest = {
        "contract": "chapter1_p1c_matched_genus_qualification_v1",
        "source_p1_contract": str(p1["contract"]),
        "source_workflow_run_id": int(p1["pinned_primary_artifact"]["workflow_run_id"]),
        "source_artifact_id": int(p1["pinned_primary_artifact"]["artifact_id"]),
        "source_artifact_digest": str(p1["pinned_primary_artifact"]["artifact_digest"]),
        "taxonomy_sha256_newline_canonicalized": taxonomy_hash,
        "qualification_seed": int(qualification["qualification_seed"]),
        "true_genus_slope_reproduction": true_checks,
        "true_reproduction_pass": true_reproduction_pass,
        "pseudo_partition_audit": partition_audit,
        "pseudo_decomposition_rows": int(len(pseudo_decomposition)),
        "pseudo_primary_profiles_fit": n_profiles_fit,
        "pseudo_expected_primary_profiles": int(len(source_modes) * len(strata)),
        "pseudo_two_axis_slope_rows_fit": int(len(pseudo_slopes)),
        "qualified_for_full_2000_permutation_null": qualified,
        "pseudo_attenuation_reported_or_stored": False,
        "new_biological_p_value_generated": False,
        "claim_boundary": (
            "Qualification establishes implementation fidelity and feasibility only. It does "
            "not compare true genus attenuation with the pseudo-genus null distribution."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "chapter1_p1c_qualification_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command()
def main(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    taxonomy_path: Path = typer.Option(..., exists=True, dir_okay=False),
    p1_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    pattern_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    source_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    explanation_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    result = run_qualification(
        artifact_root=artifact_root,
        taxonomy_path=taxonomy_path,
        p1_config_path=p1_config_path,
        pattern_config_path=pattern_config_path,
        source_config_path=source_config_path,
        explanation_config_path=explanation_config_path,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    app()
