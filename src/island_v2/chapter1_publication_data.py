from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import typer

APP = typer.Typer(add_completion=False)

TRAIT_FAMILIES = {
    "generalized_form": "floral_accessibility",
    "actinomorphic_symmetry": "floral_accessibility",
    "shallow_open_tube": "floral_accessibility",
    "self_compatibility": "reproductive_assurance",
    "selfing_mating_system": "reproductive_assurance",
    "autonomous_selfing": "reproductive_assurance",
}

INTERNAL_PUBLICATION_TERMS = ("v13", "workflow", "artifact", "main figure")

PINNED_ARTIFACTS = {
    "isolation": {
        "run_id": 29228212586,
        "name": "purpose-shortest-distance-regime-29228212586",
        "digest": "sha256:695f35b97bae07e81b05deab537dd73fa687b2d99e9efdb1cb3babd2fa12dfb6",
    },
    "atomic": {
        "run_id": 34961775336,
        "name": "chapter1-all-data-probability-34961775336",
        "digest": "sha256:af7ad0cc68c5e475fd13cc1be309dbb118f21f84b4a1cb11277da36daa11d87b",
    },
    "glopl_overlap": {
        "run_id": 35085378171,
        "name": "chapter1-h5-glopl-island-overlap-preflight-35085378171",
        "digest": "sha256:5ab7aecfbeb4fac2db26bcdc13edd4b47e953b1d1e6f35255330a0b2f70fe041",
    },
    "glopl_global": {
        "run_id": 35090599662,
        "name": "chapter1-h5-glopl-global-distance-35090599662",
        "digest": "sha256:f11046d441fba6cdc3a2842a5a19fdc0bed3661b9ce258625eb93acc879ce623",
    },
    "glopl_shape": {
        "run_id": 35091385147,
        "name": "chapter1-h5-glopl-global-shape-audit-35091385147",
        "digest": "sha256:ee7700a03dcbe9b0d6184195c482e3ec7915a5e7a9e8f1c6a85d0e4fe2de8e44",
    },
    "functional_bridge": {
        "run_id": 35141624253,
        "name": "chapter1-v13-functional-bridge-35141624253",
        "digest": "sha256:50b6ca693ec989690854a9e9aadb5cd8a672b5fbcf3ff84e5aabbab06a6abc8e",
    },
}


@dataclass(frozen=True)
class PublicationInputs:
    isolation_dir: Path
    atomic_dir: Path
    glopl_overlap_dir: Path
    glopl_global_dir: Path
    glopl_shape_dir: Path
    functional_bridge_dir: Path
    unified_lock: Path
    functional_lock: Path
    world_geojson: Path | None = None


def publication_text_is_clean(text: str) -> bool:
    lower = str(text).lower()
    return not any(term in lower for term in INTERNAL_PUBLICATION_TERMS)


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    return path


def _ci_frame(frame: pd.DataFrame, estimate: str, se: str) -> pd.DataFrame:
    result = frame.copy()
    est = pd.to_numeric(result[estimate], errors="coerce")
    stderr = pd.to_numeric(result[se], errors="coerce")
    result["ci_low"] = est - 1.96 * stderr
    result["ci_high"] = est + 1.96 * stderr
    return result


def _materialize_map_data(inputs: PublicationInputs, output_dir: Path) -> tuple[Path, Path, Path]:
    islands = pd.read_csv(
        _require(inputs.isolation_dir / "results/purpose_shortest_island_data.csv")
    )
    islands["trait_informed"] = (
        pd.to_numeric(islands["n_trait_species"], errors="coerce").fillna(0) > 0
    )
    if islands["island_id"].nunique() != 8265:
        raise ValueError("frozen island universe must contain 8,265 unique islands")
    if islands.loc[islands["trait_informed"], "island_id"].nunique() != 4453:
        raise ValueError("trait-informed input set must contain 4,453 unique islands")

    island_cols = [
        "island_id",
        "island_latitude",
        "island_longitude",
        "analysis_regime",
        "n_trait_species",
        "trait_informed",
        "distance_to_continent_km",
        "area_km2",
    ]
    island_out = output_dir / "figure1_islands.csv"
    islands[island_cols].to_csv(island_out, index=False)

    overlap = pd.read_csv(
        _require(inputs.glopl_overlap_dir / "GLOPL_METADATA_ISLAND_MATCH.csv.gz")
    )
    sites = (
        overlap.sort_values("island_id", na_position="last")
        .drop_duplicates("site_key")
        [["site_key", "Latitude", "Longitude", "island_id", "analysis_regime"]]
        .copy()
    )
    sites["is_frozen_island_site"] = sites["island_id"].notna()
    if sites["site_key"].nunique() != 1248:
        raise ValueError("GloPL sampling frame must contain 1,248 unique sites")
    if sites.loc[sites["is_frozen_island_site"], "site_key"].nunique() != 197:
        raise ValueError("exact-island GloPL overlap must contain 197 unique sites")
    sites_out = output_dir / "figure1_glopl_sites.csv"
    sites.to_csv(sites_out, index=False)

    matched = overlap.loc[overlap["island_id"].notna()].copy()
    counts = (
        matched.groupby("island_id")
        .agg(
            n_glopl_sites=("site_key", "nunique"),
            n_glopl_rows=("row_id", "size"),
            n_glopl_studies=("study_key", "nunique"),
        )
        .reset_index()
    )
    tested = islands[island_cols].merge(counts, on="island_id", how="inner", validate="one_to_one")
    if tested["island_id"].nunique() != 37:
        raise ValueError("exact GloPL overlap must contain 37 frozen islands")
    if not tested["trait_informed"].all():
        raise ValueError("all GloPL-tested frozen islands must be trait-informed inputs")
    tested_out = output_dir / "extended_data_table6_glopl_tested_islands.csv"
    tested.sort_values(["analysis_regime", "island_id"]).to_csv(tested_out, index=False)
    return island_out, sites_out, tested_out


def _materialize_plant_results(
    inputs: PublicationInputs, output_dir: Path, unified: dict[str, object]
) -> tuple[Path, Path, Path]:
    atomic_parts: list[pd.DataFrame] = []
    omnibus_parts: list[pd.DataFrame] = []
    for directory, scope in (("all", "all_analysis"), ("direct", "direct_only")):
        slopes = pd.read_csv(
            _require(inputs.atomic_dir / directory / "beta_binomial_within_slopes.csv")
        )
        slopes = slopes.loc[slopes["stratum"].eq("all_observed")].copy()
        slopes["evidence_scope"] = scope
        slopes["family"] = slopes["outcome"].map(TRAIT_FAMILIES)
        slopes = slopes.rename(
            columns={
                "geography_slope_log_odds": "estimate",
                "cluster_robust_se": "se",
                "p_value": "two_sided_p",
            }
        )
        slopes = _ci_frame(slopes, "estimate", "se")
        atomic_parts.append(slopes)

        omnibus = pd.read_csv(
            _require(inputs.atomic_dir / directory / "beta_binomial_within_omnibus.csv")
        )
        omnibus = omnibus.loc[omnibus["stratum"].eq("all_observed")].copy()
        omnibus["evidence_scope"] = scope
        omnibus_parts.append(omnibus)

    atomic = pd.concat(atomic_parts, ignore_index=True)
    omnibus = pd.concat(omnibus_parts, ignore_index=True)
    if len(atomic) != 48:
        raise ValueError("publication atomic table must contain 48 rows")
    if int(omnibus.loc[omnibus["evidence_scope"].eq("all_analysis"), "n_unique_islands"].sum()) != 4334:
        raise ValueError("primary complete-model support must sum to 4,334 context-specific islands")
    if int(omnibus.loc[omnibus["evidence_scope"].eq("direct_only"), "n_unique_islands"].sum()) != 4256:
        raise ValueError("Direct-only complete-model support must sum to 4,256 context-specific islands")

    atomic_out = output_dir / "figure2_atomic_coefficients.csv"
    omnibus_out = output_dir / "figure2_omnibus_support.csv"
    atomic.to_csv(atomic_out, index=False)
    omnibus.to_csv(omnibus_out, index=False)

    h1 = unified["H1_global_recurrent_island_syndrome"]
    h3 = unified["H3_dual_plant_response_pathways"]
    summary_rows: list[dict[str, object]] = []
    for scope_key, scope_label in (("all_analysis", "all_analysis"), ("direct_only", "direct_only")):
        for context, values in h1[scope_key].items():
            summary_rows.append(
                {
                    "evidence_scope": scope_label,
                    "context": context,
                    "summary": "classic_orientation",
                    "estimate": values["descriptive_mean_classic_slope"],
                    "inferential_role": "descriptive",
                }
            )
        family_key = (
            "all_analysis_descriptive_family_mean_slopes"
            if scope_key == "all_analysis"
            else "direct_only_descriptive_family_mean_slopes"
        )
        for family, family_block in (
            ("reproductive_assurance", h3["reproductive_assurance"]),
            ("floral_accessibility", h3["generalized_accessibility"]),
        ):
            for context, estimate in family_block[family_key].items():
                summary_rows.append(
                    {
                        "evidence_scope": scope_label,
                        "context": context,
                        "summary": family,
                        "estimate": estimate,
                        "inferential_role": "descriptive",
                    }
                )
    summaries = pd.DataFrame(summary_rows)
    summary_out = output_dir / "figure2_descriptive_summaries.csv"
    summaries.to_csv(summary_out, index=False)
    return atomic_out, omnibus_out, summary_out


def _global_effect_rows(global_result: dict[str, object], shape_result: dict[str, object]) -> pd.DataFrame:
    global_gradient = global_result["global_gradient"]
    rows = [
        {
            "row_id": "global_distance",
            "section": "global_distance",
            "label": "Global distance effect",
            "estimate": global_gradient["distance_slope"],
            "se": global_gradient["distance_slope_se"],
            "two_sided_p": global_gradient["two_sided_p"],
            "one_sided_p": global_gradient["one_sided_positive_p"],
            "n_cells": global_gradient["n_cells"],
            "n_sites": global_gradient["n_sites"],
            "n_publications": global_gradient["n_publications"],
            "inferential_role": "primary",
            "evaluable": True,
        }
    ]
    for key, label in (
        ("supplemental_only", "Supplemental-only"),
        ("no_zero_constant", "No-zero-constant"),
    ):
        value = global_result["sensitivities"][key]["global"]
        rows.append(
            {
                "row_id": f"global_{key}",
                "section": "global_distance",
                "label": label,
                "estimate": value["distance_slope"],
                "se": value["distance_slope_se"],
                "two_sided_p": value["two_sided_p"],
                "one_sided_p": value["one_sided_positive_p"],
                "n_cells": value["n_cells"],
                "n_sites": value["n_sites"],
                "n_publications": value["n_publications"],
                "inferential_role": "sensitivity",
                "evaluable": True,
            }
        )
    two_part = shape_result["two_part_global"]
    rows.extend(
        [
            {
                "row_id": "mainland_to_offshore_step",
                "section": "shape_diagnostic",
                "label": "Mainland-to-offshore step",
                "estimate": two_part["offshore_step"],
                "se": two_part["offshore_step_se"],
                "two_sided_p": two_part["offshore_step_two_sided_p"],
                "one_sided_p": np.nan,
                "n_cells": two_part["n_cells"],
                "n_sites": two_part["n_sites"],
                "n_publications": two_part["n_publications"],
                "inferential_role": "posthoc_shape_diagnostic",
                "evaluable": True,
            },
            {
                "row_id": "within_offshore_gradient",
                "section": "shape_diagnostic",
                "label": "Within-offshore gradient",
                "estimate": two_part["within_offshore_slope"],
                "se": two_part["within_offshore_slope_se"],
                "two_sided_p": two_part["within_offshore_slope_two_sided_p"],
                "one_sided_p": np.nan,
                "n_cells": two_part["n_cells"],
                "n_sites": two_part["n_sites"],
                "n_publications": two_part["n_publications"],
                "inferential_role": "posthoc_shape_diagnostic",
                "evaluable": True,
            },
        ]
    )
    offshore = shape_result["offshore_only"]
    rows.append(
        {
            "row_id": "offshore_only_gradient",
            "section": "shape_diagnostic",
            "label": "Offshore-only gradient",
            "estimate": offshore["distance_slope"],
            "se": offshore["distance_slope_se"],
            "two_sided_p": offshore["distance_slope_two_sided_p"],
            "one_sided_p": np.nan,
            "n_cells": offshore["n_cells"],
            "n_sites": offshore["n_sites"],
            "n_publications": offshore["n_publications"],
            "inferential_role": "posthoc_shape_diagnostic",
            "evaluable": True,
        }
    )
    return _ci_frame(pd.DataFrame(rows), "estimate", "se")


def _materialize_glopl_and_bridge(
    inputs: PublicationInputs,
    output_dir: Path,
    unified: dict[str, object],
) -> tuple[Path, Path, Path]:
    global_result = _read_json(_require(inputs.glopl_global_dir / "out/RESULT.json"))
    shape_result = _read_json(_require(inputs.glopl_shape_dir / "RESULT.json"))
    global_rows = _global_effect_rows(global_result, shape_result)

    bridge = pd.read_csv(
        _require(inputs.functional_bridge_dir / "functional_bridge_trait_results.csv")
    )
    bridge["row_id"] = "trait_" + bridge["trait"].astype(str) + "_" + bridge["analysis"].astype(str)
    bridge["section"] = "functional_bridge"
    bridge["label"] = bridge["trait"].astype(str).str.replace("_", " ")
    bridge["estimate"] = pd.to_numeric(bridge["estimate"], errors="coerce")
    bridge["se"] = pd.to_numeric(bridge["se"], errors="coerce")
    bridge["two_sided_p"] = pd.to_numeric(bridge["two_sided_p"], errors="coerce")
    bridge["one_sided_p"] = pd.to_numeric(
        bridge["one_sided_negative_p"], errors="coerce"
    )
    bridge = _ci_frame(bridge, "estimate", "se")

    moderation = unified["H4_functional_bridge"]["parent_distance_by_trait_moderation"]
    boundary = pd.DataFrame(
        [
            {
                "row_id": "parent_reproductive_assurance_moderation",
                "section": "negative_result_boundary",
                "label": "Distance × reproductive assurance",
                "estimate": np.nan,
                "se": np.nan,
                "two_sided_p": np.nan,
                "one_sided_p": np.nan,
                "inferential_role": "frozen_negative_result",
                "evaluable": True,
                "reason": "not supported",
                "family": "reproductive_assurance",
                "analysis": "parent_moderation",
            },
            {
                "row_id": "parent_floral_architecture_moderation",
                "section": "negative_result_boundary",
                "label": "Distance × floral architecture",
                "estimate": np.nan,
                "se": np.nan,
                "two_sided_p": np.nan,
                "one_sided_p": np.nan,
                "inferential_role": "frozen_negative_result",
                "evaluable": True,
                "reason": "not supported",
                "family": "floral_architecture",
                "analysis": "parent_moderation",
            },
        ]
    )
    if moderation["reproductive_assurance_family_supported"]:
        raise ValueError("frozen reproductive-assurance parent moderation must remain unsupported")
    if moderation["floral_architecture_family_supported"]:
        raise ValueError("frozen floral-architecture parent moderation must remain unsupported")

    common_cols = sorted(set(global_rows.columns) | set(bridge.columns) | set(boundary.columns))
    figure3 = pd.concat(
        [
            global_rows.reindex(columns=common_cols),
            bridge.reindex(columns=common_cols),
            boundary.reindex(columns=common_cols),
        ],
        ignore_index=True,
    )
    figure3_out = output_dir / "figure3_source_data.csv"
    figure3.to_csv(figure3_out, index=False)

    glopl_out = output_dir / "extended_data_table4_glopl_estimates.csv"
    global_rows.to_csv(glopl_out, index=False)
    bridge_out = output_dir / "extended_data_table5_functional_bridge.csv"
    bridge.to_csv(bridge_out, index=False)
    return figure3_out, glopl_out, bridge_out


def _materialize_data_layers(
    output_dir: Path,
    map_islands: Path,
    glopl_sites: Path,
    atomic: Path,
    figure3: Path,
) -> Path:
    islands = pd.read_csv(map_islands)
    sites = pd.read_csv(glopl_sites)
    atomic_frame = pd.read_csv(atomic)
    figure3_frame = pd.read_csv(figure3)
    primary_atomic = atomic_frame.loc[atomic_frame["evidence_scope"].eq("all_analysis")]
    global_row = figure3_frame.loc[figure3_frame["row_id"].eq("global_distance")].iloc[0]
    data = pd.DataFrame(
        [
            {
                "evidence_layer": "Geographic universe",
                "estimand": "Frozen island sampling frame",
                "observations": len(islands),
                "sites_or_islands": islands["island_id"].nunique(),
                "publications": np.nan,
                "inferential_role": "sampling frame",
            },
            {
                "evidence_layer": "Plant traits",
                "estimand": "Trait-informed island inputs",
                "observations": int(islands["trait_informed"].sum()),
                "sites_or_islands": islands.loc[islands["trait_informed"], "island_id"].nunique(),
                "publications": np.nan,
                "inferential_role": "input support",
            },
            {
                "evidence_layer": "Plant response",
                "estimand": "Six-atomic isolation response",
                "observations": len(primary_atomic),
                "sites_or_islands": 4334,
                "publications": np.nan,
                "inferential_role": "primary",
            },
            {
                "evidence_layer": "Experimental pollen limitation",
                "estimand": "Global distance effect",
                "observations": 2969,
                "sites_or_islands": sites["site_key"].nunique(),
                "publications": int(global_row["n_publications"]),
                "inferential_role": "primary",
            },
            {
                "evidence_layer": "Exact-island overlap",
                "estimand": "GloPL sites within frozen islands",
                "observations": int(sites["is_frozen_island_site"].sum()),
                "sites_or_islands": int(sites.loc[sites["is_frozen_island_site"], "site_key"].nunique()),
                "publications": np.nan,
                "inferential_role": "geographic overlap",
            },
        ]
    )
    path = output_dir / "extended_data_table1_data_layers.csv"
    data.to_csv(path, index=False)
    return path


def build_publication_source_data(
    inputs: PublicationInputs, output_dir: Path
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    unified = _read_json(_require(inputs.unified_lock))
    functional = _read_json(_require(inputs.functional_lock))
    if unified.get("status") != "submission_candidate_global_only_frozen":
        raise ValueError("unexpected unified result-lock status")
    if functional.get("status") != "reproduced_and_locked":
        raise ValueError("unexpected functional-bridge lock status")

    figure1_islands, figure1_glopl_sites, tested_islands = _materialize_map_data(
        inputs, output_dir
    )
    figure2_atomic, figure2_omnibus, figure2_summaries = _materialize_plant_results(
        inputs, output_dir, unified
    )
    figure3, table4, table5 = _materialize_glopl_and_bridge(
        inputs, output_dir, unified
    )

    table2 = output_dir / "extended_data_table2_atomic_coefficients.csv"
    shutil.copyfile(figure2_atomic, table2)
    table3 = output_dir / "extended_data_table3_descriptive_summaries.csv"
    shutil.copyfile(figure2_summaries, table3)
    table1 = _materialize_data_layers(
        output_dir, figure1_islands, figure1_glopl_sites, figure2_atomic, figure3
    )

    if inputs.world_geojson and inputs.world_geojson.is_file():
        shutil.copyfile(inputs.world_geojson, output_dir / "world_land.geojson")

    paths = {
        "figure1_islands": figure1_islands,
        "figure1_glopl_sites": figure1_glopl_sites,
        "figure2_atomic": figure2_atomic,
        "figure2_omnibus": figure2_omnibus,
        "figure2_summaries": figure2_summaries,
        "figure3": figure3,
        "extended_data_table1": table1,
        "extended_data_table2": table2,
        "extended_data_table3": table3,
        "extended_data_table4": table4,
        "extended_data_table5": table5,
        "extended_data_table6": tested_islands,
    }
    manifest = {
        "contract": "chapter1_publication_figure_bundle_v1",
        "counts": {
            "frozen_islands": 8265,
            "trait_informed_islands": 4453,
            "primary_complete_model_context_sum": 4334,
            "direct_complete_model_context_sum": 4256,
            "glopl_effect_rows": 2969,
            "glopl_sites": 1248,
            "glopl_island_sites": 197,
            "glopl_tested_islands": 37,
        },
        "pinned_inputs": PINNED_ARTIFACTS,
        "source_files": {
            key: {"path": path.name, "sha256": _sha256(path)} for key, path in paths.items()
        },
    }
    manifest_path = output_dir / "publication_figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    paths["manifest"] = manifest_path
    return paths


@APP.command()
def cli(
    input_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    unified_lock: Path = typer.Option(
        Path("config/chapter1_v13_unified_island_syndrome_result_lock.json")
    ),
    functional_lock: Path = typer.Option(
        Path("config/chapter1_v13_functional_bridge_result_lock.json")
    ),
) -> None:
    world_geojson = input_root / "natural_earth/ne_110m_land.geojson"
    inputs = PublicationInputs(
        isolation_dir=input_root / "isolation",
        atomic_dir=input_root / "atomic",
        glopl_overlap_dir=input_root / "glopl_overlap",
        glopl_global_dir=input_root / "glopl_global",
        glopl_shape_dir=input_root / "glopl_shape",
        functional_bridge_dir=input_root / "functional_bridge",
        unified_lock=unified_lock,
        functional_lock=functional_lock,
        world_geojson=world_geojson if world_geojson.is_file() else None,
    )
    build_publication_source_data(inputs, output_dir)


if __name__ == "__main__":
    APP()
