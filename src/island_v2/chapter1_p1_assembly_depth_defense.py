"""Post-baseline P1 defense of Chapter 1 taxonomic assembly depth.

P1a reconstructs the frozen family/genus decomposition from the pinned PR142 artifact
and verifies exact common island-species support against the frozen decomposition.
P1b estimates paired spatial-block bootstrap uncertainty for attenuation itself.

This module cannot change the historical H3 classification and does not implement the
matched-complexity pseudo-genus null (P1c), which opens only after P1a/P1b are frozen.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml
from scipy import sparse

from island_v2.chapter1_pr138_lineage_representation_bridge import (
    _source_assignment_matrix,
    match_gift_species,
)
from island_v2.chapter1_taxonomic_depth_decomposition import (
    _availability_matrix,
    _stratum_mask,
    _validate_scores,
    build_source_group_contract,
    file_sha256,
    validate_taxonomy,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

STAGES = ("observed_score", "after_family_residual", "after_genus_residual")
AXES = ("generalized_accessible", "selfing_core")
SCOPES = ("direct_only", "all_analysis_eligible")


class P1ContractError(ValueError):
    """Raised when the frozen P1 contract or artifact does not match expectations."""


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise P1ContractError(f"expected mapping in {path}")
    return value


def _hash_rows(frame: pd.DataFrame, columns: list[str]) -> str:
    ordered = frame[columns].copy()
    for column in columns:
        ordered[column] = ordered[column].fillna("").astype(str)
    ordered = ordered.sort_values(columns, kind="stable").reset_index(drop=True)
    payload = ordered.to_csv(index=False, lineterminator="\n")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _find_table(root: Path, stem: str) -> Path:
    for suffix in (".csv", ".csv.gz"):
        path = root / f"{stem}{suffix}"
        if path.is_file():
            return path
    raise FileNotFoundError(root / stem)


def _membership_and_summary(
    species_scores: pd.DataFrame,
    taxonomy: pd.DataFrame,
    status_flora: pd.DataFrame,
    gift_flora: pd.DataFrame,
    assignments: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    axes: list[str],
    source_modes: list[str],
    strata: list[str],
    minimum_species: int,
    minimum_families: int,
    minimum_genera: int,
    minimum_source_scored_species: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reconstruct exact species membership used by the frozen V2 decomposition."""
    scores = _validate_scores(species_scores, axes)
    matched_taxa = match_gift_species(gift_flora, taxonomy[["accepted_species"]])
    positions, availability, _, _ = build_source_group_contract(
        scores,
        taxonomy,
        gift_flora,
        axes=axes,
        minimum_source_scored_species=minimum_source_scored_species,
        matched_taxa=matched_taxa,
    )

    position_lookup: dict[str, pd.Series] = {}
    group_indices: dict[str, dict[str, int]] = {}
    presence: dict[str, sparse.csr_matrix] = {}
    entities = sorted(
        set(pd.to_numeric(assignments["entity_ID"], errors="coerce").dropna().astype(int))
        | {
            int(value)
            for table in availability.values()
            for value in pd.to_numeric(table["entity_ID"], errors="coerce").dropna().astype(int)
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
    required_status = ["island_id", "accepted_species", "origin_status", "floristic_status"]
    flora = status_flora[required_status].drop_duplicates(["island_id", "accepted_species"])
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

    membership_parts: list[pd.DataFrame] = []
    summary_parts: list[pd.DataFrame] = []
    for source_mode in source_modes:
        assignment = _source_assignment_matrix(
            assignments,
            island_index,
            entity_index,
            source_mode=source_mode,
        )
        available_by_level = {
            level: (assignment @ presence[level]).tocsr() for level in ("family", "genus")
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
            support = (
                work.groupby(["island_id", "syndrome"], as_index=False)
                .agg(
                    n_species=("accepted_species", "nunique"),
                    n_families=("family", "nunique"),
                    n_genera=("genus", "nunique"),
                )
            )
            support = support.loc[
                support["n_species"].ge(minimum_species)
                & support["n_families"].ge(minimum_families)
                & support["n_genera"].ge(minimum_genera)
            ].copy()
            if support.empty:
                continue
            work = work.merge(
                support[["island_id", "syndrome"]],
                on=["island_id", "syndrome"],
                how="inner",
                validate="many_to_one",
            )
            work = work.drop_duplicates(["island_id", "accepted_species", "syndrome"])
            member = work[
                ["island_id", "accepted_species", "syndrome", "family", "genus"]
            ].copy()
            member["source_mode"] = source_mode
            member["stratum"] = stratum
            membership_parts.append(member)

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
            summary["after_family_residual"] = summary["observed_score"] - summary["family_expected"]
            summary["family_to_genus_increment"] = summary["genus_expected"] - summary["family_expected"]
            summary["after_genus_residual"] = summary["observed_score"] - summary["genus_expected"]
            summary["source_mode"] = source_mode
            summary["stratum"] = stratum
            summary_parts.append(summary)

    membership = pd.concat(membership_parts, ignore_index=True) if membership_parts else pd.DataFrame()
    summary = pd.concat(summary_parts, ignore_index=True) if summary_parts else pd.DataFrame()
    return membership, summary


def _compare_reconstructed_to_frozen(
    reconstructed: pd.DataFrame, frozen: pd.DataFrame
) -> tuple[float, bool]:
    keys = ["island_id", "syndrome", "source_mode", "stratum"]
    value_cols = [
        "observed_score",
        "family_expected",
        "genus_expected",
        "n_species",
        "n_families",
        "n_genera",
        "after_family_residual",
        "family_to_genus_increment",
        "after_genus_residual",
    ]
    left = reconstructed.sort_values(keys).reset_index(drop=True)
    right = frozen.sort_values(keys).reset_index(drop=True)
    same_keys = left[keys].astype(str).equals(right[keys].astype(str))
    if not same_keys or len(left) != len(right):
        return float("inf"), False
    max_diff = 0.0
    for column in value_cols:
        a = pd.to_numeric(left[column], errors="coerce").to_numpy(float)
        b = pd.to_numeric(right[column], errors="coerce").to_numpy(float)
        if not np.isfinite(a).all() or not np.isfinite(b).all():
            return float("inf"), False
        max_diff = max(max_diff, float(np.max(np.abs(a - b))))
    return max_diff, bool(max_diff <= 1.0e-12)


def _weighted_standardized_slope(
    frame: pd.DataFrame,
    response: str,
    block_counts: dict[str, int],
    baseline: list[str],
    geography: str,
) -> float:
    weights = frame["spatial_block"].astype(str).map(block_counts).fillna(0).to_numpy(float)
    keep = weights > 0
    if int(keep.sum()) < len(baseline) + 3:
        return float("nan")
    data = frame.loc[keep]
    weights = weights[keep]
    columns: list[np.ndarray] = []
    for name in [*baseline, geography]:
        values = pd.to_numeric(data[name], errors="coerce").to_numpy(float)
        if not np.isfinite(values).all():
            return float("nan")
        mean = float(np.average(values, weights=weights))
        variance = float(np.average((values - mean) ** 2, weights=weights))
        if not math.isfinite(variance) or variance <= 0:
            return float("nan")
        columns.append((values - mean) / math.sqrt(variance))
    design = np.column_stack([np.ones(len(data)), *columns])
    outcome = pd.to_numeric(data[response], errors="coerce").to_numpy(float)
    if not np.isfinite(outcome).all():
        return float("nan")
    root_w = np.sqrt(weights)
    beta, *_ = np.linalg.lstsq(design * root_w[:, None], outcome * root_w, rcond=None)
    return float(beta[-1])


def _prepare_profile(
    decomposition: pd.DataFrame,
    covariates: pd.DataFrame,
    realm: pd.DataFrame,
    *,
    source_mode: str,
    stratum: str,
    target_context: str,
    geography: str,
    baseline: list[str],
) -> dict[str, pd.DataFrame]:
    realm_small = realm[["island_id", "biogeographic_realm"]].drop_duplicates("island_id")
    cov = covariates[["island_id", "spatial_block", geography, *baseline]].drop_duplicates("island_id")
    cov = cov.merge(realm_small, on="island_id", how="left", validate="one_to_one")
    subset = decomposition.loc[
        decomposition["source_mode"].astype(str).eq(source_mode)
        & decomposition["stratum"].astype(str).eq(stratum)
        & decomposition["syndrome"].astype(str).isin(AXES)
    ].copy()
    prepared: dict[str, pd.DataFrame] = {}
    for axis in AXES:
        frame = subset.loc[subset["syndrome"].astype(str).eq(axis)].merge(
            cov, on="island_id", how="inner", validate="one_to_one"
        )
        frame = frame.loc[frame["biogeographic_realm"].astype(str).eq(target_context)].copy()
        frame = frame.dropna(subset=[*STAGES, geography, *baseline, "spatial_block"])
        prepared[axis] = frame
    return prepared


def _stage_vector_norm(
    prepared: dict[str, pd.DataFrame],
    stage: str,
    block_counts: dict[str, int],
    baseline: list[str],
    geography: str,
) -> float:
    slopes = np.asarray(
        [
            _weighted_standardized_slope(
                prepared[axis], stage, block_counts, baseline, geography
            )
            for axis in AXES
        ],
        dtype=float,
    )
    return float(np.linalg.norm(slopes)) if np.isfinite(slopes).all() else float("nan")


def _attenuation_from_norms(observed: float, family: float, genus: float, floor: float) -> dict[str, float]:
    if not all(math.isfinite(x) for x in (observed, family, genus)) or observed <= floor or family <= floor:
        return {key: float("nan") for key in (
            "family_attenuation",
            "genus_attenuation",
            "family_to_genus_extra_attenuation",
            "conditional_genus_attenuation",
        )}
    return {
        "family_attenuation": 1.0 - family / observed,
        "genus_attenuation": 1.0 - genus / observed,
        "family_to_genus_extra_attenuation": (family - genus) / observed,
        "conditional_genus_attenuation": 1.0 - genus / family,
    }


def _bootstrap_scope(
    decomposition: pd.DataFrame,
    covariates: pd.DataFrame,
    realm: pd.DataFrame,
    *,
    scope: str,
    source_modes: list[str],
    strata: list[str],
    baseline: list[str],
    geography: str,
    target_context: str,
    draws: int,
    seed: int,
    floor: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    draw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    seed_sequence = np.random.SeedSequence(seed)
    child_seeds = seed_sequence.spawn(len(source_modes) * len(strata))
    profile_index = 0
    for source_mode in source_modes:
        for stratum in strata:
            prepared = _prepare_profile(
                decomposition,
                covariates,
                realm,
                source_mode=source_mode,
                stratum=stratum,
                target_context=target_context,
                geography=geography,
                baseline=baseline,
            )
            blocks = sorted(
                set().union(
                    *[set(frame["spatial_block"].astype(str)) for frame in prepared.values()]
                )
            )
            if len(blocks) < 2:
                raise P1ContractError(f"too few spatial blocks for {scope}/{source_mode}/{stratum}")
            full_counts = {block: 1 for block in blocks}
            full_norms = [
                _stage_vector_norm(prepared, stage, full_counts, baseline, geography)
                for stage in STAGES
            ]
            full_metrics = _attenuation_from_norms(*full_norms, floor)
            rng = np.random.default_rng(child_seeds[profile_index])
            profile_index += 1
            for draw in range(draws):
                sampled = rng.choice(np.asarray(blocks, dtype=object), size=len(blocks), replace=True)
                values, counts = np.unique(sampled, return_counts=True)
                block_counts = {str(value): int(count) for value, count in zip(values, counts, strict=True)}
                norms = [
                    _stage_vector_norm(prepared, stage, block_counts, baseline, geography)
                    for stage in STAGES
                ]
                metrics = _attenuation_from_norms(*norms, floor)
                draw_rows.append(
                    {
                        "evidence_scope": scope,
                        "source_mode": source_mode,
                        "stratum": stratum,
                        "draw": draw,
                        "n_spatial_blocks": len(blocks),
                        **metrics,
                    }
                )
            draw_frame = pd.DataFrame(
                [
                    row
                    for row in draw_rows
                    if row["evidence_scope"] == scope
                    and row["source_mode"] == source_mode
                    and row["stratum"] == stratum
                ]
            )
            rec: dict[str, Any] = {
                "evidence_scope": scope,
                "source_mode": source_mode,
                "stratum": stratum,
                "n_spatial_blocks": len(blocks),
                "n_valid_draws": int(draw_frame["family_to_genus_extra_attenuation"].notna().sum()),
                **{f"observed_{key}": value for key, value in full_metrics.items()},
            }
            for metric in full_metrics:
                values = pd.to_numeric(draw_frame[metric], errors="coerce").dropna().to_numpy(float)
                if len(values) == 0:
                    rec[f"{metric}_ci_low"] = np.nan
                    rec[f"{metric}_median"] = np.nan
                    rec[f"{metric}_ci_high"] = np.nan
                else:
                    low, med, high = np.quantile(values, [0.025, 0.5, 0.975])
                    rec[f"{metric}_ci_low"] = float(low)
                    rec[f"{metric}_median"] = float(med)
                    rec[f"{metric}_ci_high"] = float(high)
            rec["extra_attenuation_ci_above_zero"] = bool(
                pd.notna(rec["family_to_genus_extra_attenuation_ci_low"])
                and float(rec["family_to_genus_extra_attenuation_ci_low"]) > 0
            )
            summary_rows.append(rec)
    return pd.DataFrame(draw_rows), pd.DataFrame(summary_rows)


def run_p1ab(
    *,
    artifact_root: Path,
    taxonomy_csv: Path,
    pattern_config_path: Path,
    source_config_path: Path,
    explanation_config_path: Path,
    p1_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    p1 = _load_yaml(p1_config_path)
    if p1.get("contract") != "chapter1_p1_assembly_depth_defense_v1":
        raise P1ContractError("unexpected P1 contract")
    pattern = _load_yaml(pattern_config_path)
    source = _load_yaml(source_config_path)
    explanation = _load_yaml(explanation_config_path)
    v2 = explanation["validations"]["V2_H3_taxonomic_depth"]["frozen_implementation"]

    manifest_path = artifact_root / "taxonomic-depth/chapter1_v2_taxonomic_depth_manifest.json"
    frozen_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if frozen_manifest.get("family_and_genus_nulls_use_common_observed_species") is not True:
        raise P1ContractError("frozen V2 manifest does not certify common observed species")

    taxonomy_hash = file_sha256(taxonomy_csv)
    taxonomy = validate_taxonomy(
        pd.read_csv(taxonomy_csv),
        str(v2["taxonomy_input"]["sha256_newline_canonicalized"]),
        taxonomy_hash,
    )
    status_flora = pd.read_csv(artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz")
    gift_flora = pd.read_csv(artifact_root / "fixed/source/source/gift_native_mainland_flora.csv.gz")
    assignments = pd.read_csv(artifact_root / "fixed/source/result/island_source_assignments.csv.gz")
    covariates = pd.read_csv(artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv")
    realm = pd.read_csv(artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv")

    source_modes = [str(x) for x in source["source_assignment"]["primary_modes"]]
    strata = [str(x) for x in v2["strata"]]
    axes = [str(x) for x in v2["plant_response_axes"]]
    minimum = v2["minimum_support"]

    audit_rows: list[dict[str, Any]] = []
    decompositions: dict[str, pd.DataFrame] = {}
    for scope, score_path in (
        ("all_analysis_eligible", artifact_root / "syndrome/all/species_syndrome_concordance.csv.gz"),
        ("direct_only", artifact_root / "syndrome/direct/species_syndrome_concordance.csv.gz"),
    ):
        scores = pd.read_csv(score_path)
        membership, reconstructed = _membership_and_summary(
            scores,
            taxonomy,
            status_flora,
            gift_flora,
            assignments,
            covariates,
            axes=axes,
            source_modes=source_modes,
            strata=strata,
            minimum_species=int(minimum["observed_species_per_island_response"]),
            minimum_families=int(minimum["represented_families"]),
            minimum_genera=int(minimum["represented_genera"]),
            minimum_source_scored_species=2,
        )
        frozen = pd.read_csv(_find_table(artifact_root / f"taxonomic-depth/{scope}", "decomposition"))
        max_diff, exact_match = _compare_reconstructed_to_frozen(reconstructed, frozen)
        membership_hash = _hash_rows(
            membership,
            ["island_id", "accepted_species", "syndrome", "family", "genus", "source_mode", "stratum"],
        )
        key_hash = _hash_rows(frozen, ["island_id", "syndrome", "source_mode", "stratum"])
        audit_rows.append(
            {
                "evidence_scope": scope,
                "n_membership_rows": len(membership),
                "n_unique_species": int(membership["accepted_species"].astype(str).nunique()),
                "n_decomposition_rows": len(frozen),
                "membership_sha256": membership_hash,
                "decomposition_key_sha256": key_hash,
                "reconstructed_frozen_max_abs_difference": max_diff,
                "exact_reconstruction_pass": exact_match,
                "common_species_manifest_pass": True,
            }
        )
        if not exact_match:
            raise P1ContractError(f"P1a reconstruction differs from frozen decomposition: {scope}")
        decompositions[scope] = frozen
        output_dir.mkdir(parents=True, exist_ok=True)
        membership.to_csv(
            output_dir / f"p1a_{scope}_common_species_membership.csv.gz",
            index=False,
            compression="gzip",
        )

    baseline = [str(x) for x in pattern["baseline_covariates"]]
    geography = str(pattern["geography_column"])
    boot = p1["p1b_paired_spatial_block_uncertainty"]
    draws = int(boot["draws"])
    seed = int(boot["seed"])
    floor = float(p1.get("zero_floor", 1.0e-12))
    draw_parts: list[pd.DataFrame] = []
    summary_parts: list[pd.DataFrame] = []
    for offset, scope in enumerate(SCOPES):
        draw_frame, summary = _bootstrap_scope(
            decompositions[scope],
            covariates,
            realm,
            scope=scope,
            source_modes=source_modes,
            strata=strata,
            baseline=baseline,
            geography=geography,
            target_context="Palearctic",
            draws=draws,
            seed=seed + offset,
            floor=floor,
        )
        draw_parts.append(draw_frame)
        summary_parts.append(summary)
    draws_out = pd.concat(draw_parts, ignore_index=True)
    summary_out = pd.concat(summary_parts, ignore_index=True)
    audits = pd.DataFrame(audit_rows)

    direct = summary_out.loc[summary_out["evidence_scope"].eq("direct_only")]
    p1b_primary_pass = bool(
        len(direct) == len(source_modes) * len(strata)
        and direct["extra_attenuation_ci_above_zero"].astype(bool).all()
    )
    p1a_pass = bool(audits["exact_reconstruction_pass"].astype(bool).all())
    if not p1a_pass:
        status = "p1a_failed_stop"
    elif not p1b_primary_pass:
        status = "p1a_pass_p1b_narrows_claim_before_p1c"
    else:
        status = "p1a_p1b_pass_p1c_may_open"

    output_dir.mkdir(parents=True, exist_ok=True)
    audits.to_csv(output_dir / "p1a_exact_pair_audit.csv", index=False)
    summary_out.to_csv(output_dir / "p1b_paired_bootstrap_summary.csv", index=False)
    draws_out.to_csv(
        output_dir / "p1b_paired_bootstrap_draws.csv.gz", index=False, compression="gzip"
    )
    manifest = {
        "contract": p1["contract"],
        "status": status,
        "pinned_workflow_run_id": int(p1["pinned_empirical_input"]["workflow_run_id"]),
        "pinned_artifact_id": int(p1["pinned_empirical_input"]["artifact_id"]),
        "pinned_artifact_digest": str(p1["pinned_empirical_input"]["digest"]),
        "p1a_exact_reconstruction_pass": p1a_pass,
        "p1b_direct_only_all_profiles_ci_above_zero": p1b_primary_pass,
        "n_direct_only_profiles": int(len(direct)),
        "n_direct_only_profiles_ci_above_zero": int(
            direct["extra_attenuation_ci_above_zero"].astype(bool).sum()
        ),
        "bootstrap_draws": draws,
        "bootstrap_seed": seed,
        "p1c_opened": False,
        "historical_H3_reclassified": False,
        "claim_boundary": (
            "P1a/P1b are post-baseline robustness extensions. Failure narrows the assembly-depth "
            "interpretation and cannot be repaired by retuning bootstrap or support rules."
        ),
    }
    (output_dir / "chapter1_p1ab_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    taxonomy_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    pattern_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    source_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    explanation_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    p1_config_path: Path = typer.Option(..., exists=True, dir_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    manifest = run_p1ab(
        artifact_root=artifact_root,
        taxonomy_csv=taxonomy_csv,
        pattern_config_path=pattern_config_path,
        source_config_path=source_config_path,
        explanation_config_path=explanation_config_path,
        p1_config_path=p1_config_path,
        output_dir=output_dir,
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
