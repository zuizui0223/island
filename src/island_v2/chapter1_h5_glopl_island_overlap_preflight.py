"""Outcome-blind GloPL x frozen-island overlap preflight for Chapter 1 H5.

This module deliberately reads only geographic/study metadata from GloPL. Pollen-
limitation effect sizes and treatment outcomes are never materialized here. The only
question is whether exact GSHHG-island overlap contains enough independent islands and
studies in both primary biogeographic contexts to justify a separately frozen service-
limitation analysis.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT = "chapter1_h5_glopl_island_overlap_preflight_v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected GloPL island-overlap preflight contract")
    return config


def read_glopl_metadata(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    """Read only predeclared non-outcome columns from the public GloPL CSV."""
    allowed = [str(x) for x in config["allowed_columns"]]
    forbidden = {str(x) for x in config["forbidden_outcome_columns"]}
    if forbidden.intersection(allowed):
        raise typer.BadParameter("forbidden GloPL outcomes overlap allowed preflight columns")
    # The archived file contains non-ASCII author strings; latin-1 is lossless for bytes
    # and is sufficient because preflight identity/geography does not parse prose.
    frame = pd.read_csv(path, usecols=allowed, dtype=str, encoding="latin-1").fillna("")
    if list(frame.columns) != allowed:
        frame = frame[allowed]
    return frame


def build_study_key(frame: pd.DataFrame) -> pd.Series:
    required = {"DOI", "Author", "Year"}
    if missing := required - set(frame.columns):
        raise ValueError(f"study-key columns missing: {sorted(missing)}")
    doi = frame["DOI"].fillna("").astype(str).str.strip().str.casefold()
    author = frame["Author"].fillna("").astype(str).str.strip()
    year = frame["Year"].fillna("").astype(str).str.strip()
    out = pd.Series("", index=frame.index, dtype=object)
    has_doi = doi.ne("")
    out.loc[has_doi] = "doi:" + doi.loc[has_doi]
    fallback = ~has_doi & author.ne("") & year.ne("")
    out.loc[fallback] = "author_year:" + author.loc[fallback] + "|" + year.loc[fallback]
    return out


def _site_key(latitude: pd.Series, longitude: pd.Series, valid: pd.Series) -> pd.Series:
    out = pd.Series("", index=latitude.index, dtype=object)
    for idx in latitude.index[valid]:
        out.loc[idx] = f"{float(latitude.loc[idx]):.15g}|{float(longitude.loc[idx]):.15g}"
    return out


def _unique_island_matches(joined: gpd.GeoDataFrame) -> dict[int, list[str]]:
    result: dict[int, list[str]] = {}
    if joined.empty:
        return result
    valid = joined.dropna(subset=["island_id"])
    for row_id, part in valid.groupby("row_id", sort=False):
        result[int(row_id)] = sorted(set(part["island_id"].astype(str)))
    return result


def match_glopl_to_islands(
    glopl: pd.DataFrame,
    islands: gpd.GeoDataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Match GloPL coordinates to frozen islands using strict unique `within` first."""
    required = set(config["allowed_columns"])
    if missing := required - set(glopl.columns):
        raise typer.BadParameter(f"GloPL metadata missing columns: {sorted(missing)}")
    if "island_id" not in islands.columns:
        raise typer.BadParameter("island geometry lacks island_id")

    work = glopl[list(config["allowed_columns"])].copy().reset_index(drop=True)
    work["row_id"] = range(len(work))
    lat = pd.to_numeric(work["Latitude"], errors="coerce")
    lon = pd.to_numeric(work["Longitude"], errors="coerce")
    limits = config["spatial_matching"]["coordinate_validity"]
    valid = (
        lat.between(float(limits["latitude_min"]), float(limits["latitude_max"]), inclusive="both")
        & lon.between(float(limits["longitude_min"]), float(limits["longitude_max"]), inclusive="both")
    )
    work["coordinate_valid"] = valid.astype(bool)
    work["site_key"] = _site_key(lat, lon, valid)
    work["study_key"] = build_study_key(work)
    work["island_id"] = pd.Series(pd.NA, index=work.index, dtype="string")
    work["boundary_only_match"] = False
    work["multi_polygon_match"] = False

    island_frame = islands[["island_id", "geometry"]].copy()
    if island_frame["island_id"].astype(str).duplicated().any():
        raise typer.BadParameter("island geometry contains duplicate island_id")
    target_crs = str(config["spatial_matching"]["crs"])
    if island_frame.crs is None:
        raise typer.BadParameter("island geometry CRS is missing")
    island_frame = island_frame.to_crs(target_crs)
    island_frame["island_id"] = island_frame["island_id"].astype(str)

    points = gpd.GeoDataFrame(
        work.loc[valid, ["row_id"]].copy(),
        geometry=gpd.points_from_xy(lon.loc[valid], lat.loc[valid]),
        crs=target_crs,
    )
    within = gpd.sjoin(
        points,
        island_frame,
        how="left",
        predicate=str(config["spatial_matching"]["primary_predicate"]),
    )
    within_matches = _unique_island_matches(within)
    for row_id, matches in within_matches.items():
        if len(matches) == 1:
            work.loc[work["row_id"].eq(row_id), "island_id"] = matches[0]
        elif len(matches) > 1:
            work.loc[work["row_id"].eq(row_id), "multi_polygon_match"] = True

    unmatched_ids = work.loc[
        work["coordinate_valid"] & work["island_id"].isna() & ~work["multi_polygon_match"],
        "row_id",
    ].astype(int)
    if len(unmatched_ids):
        sensitivity_points = points.loc[points["row_id"].isin(unmatched_ids)].copy()
        intersects = gpd.sjoin(
            sensitivity_points,
            island_frame,
            how="left",
            predicate=str(config["spatial_matching"]["boundary_sensitivity_predicate"]),
        )
        sensitivity_matches = _unique_island_matches(intersects)
        for row_id, matches in sensitivity_matches.items():
            if matches:
                # Boundary/intersection evidence is intentionally never promoted to
                # the primary island_id assignment.
                work.loc[work["row_id"].eq(row_id), "boundary_only_match"] = True

    audit = {
        "n_rows_total": int(len(work)),
        "n_rows_valid_coordinates": int(work["coordinate_valid"].sum()),
        "n_unique_coordinate_sites": int(work.loc[work["coordinate_valid"], "site_key"].nunique()),
        "n_rows_within_any_frozen_island": int(work["island_id"].notna().sum()),
        "n_unique_sites_within_any_frozen_island": int(
            work.loc[work["island_id"].notna(), "site_key"].nunique()
        ),
        "boundary_only_match_count": int(work["boundary_only_match"].sum()),
        "multi_polygon_match_count": int(work["multi_polygon_match"].sum()),
    }
    return work, audit


def _count_nonblank_unique(series: pd.Series) -> int:
    clean = series.fillna("").astype(str).str.strip()
    return int(clean.loc[clean.ne("")].nunique())


def summarize_overlap(matched: pd.DataFrame, config: dict[str, Any]) -> dict[str, Any]:
    context_column = str(config["context_source"]["context_column"])
    primary_contexts = [str(x) for x in config["context_source"]["primary_contexts"]]
    primary = matched.loc[matched["island_id"].notna()].copy()
    coordinate_valid = (
        matched["coordinate_valid"].astype(bool)
        if "coordinate_valid" in matched.columns
        else matched.get("site_key", pd.Series("", index=matched.index)).fillna("").astype(str).ne("")
    )
    summary: dict[str, Any] = {
        "n_rows_total": int(len(matched)),
        "n_rows_valid_coordinates": int(coordinate_valid.sum()),
        "n_unique_coordinate_sites": _count_nonblank_unique(matched.get("site_key", pd.Series(dtype=object))),
        "n_rows_within_any_frozen_island": int(len(primary)),
        "n_unique_sites_within_any_frozen_island": _count_nonblank_unique(primary.get("site_key", pd.Series(dtype=object))),
        "n_unique_islands_with_GloPL": _count_nonblank_unique(primary["island_id"]),
        "n_unique_studies_within_any_frozen_island": _count_nonblank_unique(primary.get("study_key", pd.Series(dtype=object))),
        "n_unique_species_within_any_frozen_island": _count_nonblank_unique(primary.get("Species_accepted_names", pd.Series(dtype=object))),
        "boundary_only_match_count": int(matched.get("boundary_only_match", pd.Series(False, index=matched.index)).fillna(False).astype(bool).sum()),
        "multi_polygon_match_count": int(matched.get("multi_polygon_match", pd.Series(False, index=matched.index)).fillna(False).astype(bool).sum()),
        "by_analysis_regime": {},
    }
    for context in primary_contexts:
        if context_column in primary.columns:
            part = primary.loc[primary[context_column].fillna("").astype(str).eq(context)]
        else:
            part = primary.iloc[0:0]
        summary["by_analysis_regime"][context] = {
            "n_rows": int(len(part)),
            "n_sites": _count_nonblank_unique(part.get("site_key", pd.Series(dtype=object))),
            "n_islands": _count_nonblank_unique(part["island_id"]) if "island_id" in part.columns else 0,
            "n_studies": _count_nonblank_unique(part.get("study_key", pd.Series(dtype=object))),
            "n_species": _count_nonblank_unique(part.get("Species_accepted_names", pd.Series(dtype=object))),
        }
    return summary


def evaluate_admission_gate(summary: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    gate = config["admission_gate"]
    contexts = [str(x) for x in config["context_source"]["primary_contexts"]]
    minimum_islands = int(gate["minimum_unique_islands_per_primary_context"])
    minimum_studies = int(gate["minimum_unique_studies_per_primary_context"])
    context_pass: dict[str, bool] = {}
    for context in contexts:
        counts = summary["by_analysis_regime"].get(context, {})
        context_pass[context] = bool(
            int(counts.get("n_islands", 0)) >= minimum_islands
            and int(counts.get("n_studies", 0)) >= minimum_studies
        )
    admitted = all(context_pass.values()) if bool(gate["both_primary_contexts_must_pass"]) else any(context_pass.values())
    return {
        "contract": CONTRACT,
        "admitted": bool(admitted),
        "context_pass": context_pass,
        "minimum_unique_islands_per_primary_context": minimum_islands,
        "minimum_unique_studies_per_primary_context": minimum_studies,
        "effect_size_analysis_allowed": bool(admitted),
        "threshold_relaxation_after_overlap_result": False,
        "classification": (
            "GloPL_exact_island_overlap_gate_passed"
            if admitted
            else "GloPL_exact_island_overlap_gate_failed"
        ),
    }


def _unique_context(covariates: pd.DataFrame, context_column: str) -> pd.DataFrame:
    required = {"island_id", context_column}
    if missing := required - set(covariates.columns):
        raise typer.BadParameter(f"context table missing columns: {sorted(missing)}")
    work = covariates[["island_id", context_column]].copy()
    multiplicity = work.groupby("island_id", dropna=False)[context_column].nunique(dropna=False)
    if bool(multiplicity.gt(1).any()):
        examples = [str(x) for x in multiplicity.index[multiplicity.gt(1)][:5]]
        raise typer.BadParameter(f"conflicting analysis_regime rows for island_id: {examples}")
    return work.drop_duplicates("island_id")


@app.command("run")
def run(
    glopl_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = load_config(config_path)
    glopl = read_glopl_metadata(glopl_csv, config)
    islands = gpd.read_file(
        islands_gpkg,
        layer=str(config["island_universe"]["geometry_layer"]),
    )
    expected_islands = int(config["island_universe"]["n_islands"])
    if len(islands) != expected_islands or islands["island_id"].astype(str).nunique() != expected_islands:
        raise typer.BadParameter(
            f"frozen island universe mismatch: expected {expected_islands}, got {len(islands)}"
        )
    matched, spatial_audit = match_glopl_to_islands(glopl, islands, config)

    context_column = str(config["context_source"]["context_column"])
    contexts = _unique_context(pd.read_csv(covariates_csv, dtype=str).fillna(""), context_column)
    matched = matched.merge(contexts, on="island_id", how="left", validate="many_to_one")
    summary = summarize_overlap(matched, config)
    # Preserve the spatial audit's explicit coordinate/match counts.
    summary.update(spatial_audit)
    decision = evaluate_admission_gate(summary, config)
    source_sha256 = _sha256(glopl_csv)

    output_dir.mkdir(parents=True, exist_ok=True)
    matched.to_csv(output_dir / "GLOPL_METADATA_ISLAND_MATCH.csv.gz", index=False, compression="gzip")
    (output_dir / "GLOPL_OVERLAP_SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "GLOPL_PREFLIGHT_DECISION.json").write_text(
        json.dumps(decision, indent=2) + "\n", encoding="utf-8"
    )
    receipt = {
        "contract": CONTRACT,
        "status": "completed_outcome_blind_preflight",
        "source_repository": config["source"]["repository"],
        "source_git_commit": config["source"]["git_commit"],
        "source_file": config["source"]["file"],
        "source_file_sha256": source_sha256,
        "allowed_columns": list(config["allowed_columns"]),
        "forbidden_outcome_columns_materialized": False,
        "island_universe_n": expected_islands,
        "summary": summary,
        "decision": decision,
        "claim_ceiling": config["claim_ceiling"],
    }
    (output_dir / "PREFLIGHT_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# H5 GloPL exact-island overlap preflight",
        "",
        f"- GloPL SHA-256: `{source_sha256}`",
        f"- valid-coordinate rows/sites: {summary['n_rows_valid_coordinates']} / {summary['n_unique_coordinate_sites']}",
        f"- exact-island rows/sites/islands: {summary['n_rows_within_any_frozen_island']} / {summary['n_unique_sites_within_any_frozen_island']} / {summary['n_unique_islands_with_GloPL']}",
        f"- island studies/species: {summary['n_unique_studies_within_any_frozen_island']} / {summary['n_unique_species_within_any_frozen_island']}",
    ]
    for context in config["context_source"]["primary_contexts"]:
        counts = summary["by_analysis_regime"][str(context)]
        lines.append(
            f"- {context}: {counts['n_islands']} islands / {counts['n_studies']} studies / {counts['n_sites']} sites / {counts['n_species']} species"
        )
    lines.extend(
        [
            f"- boundary-only / multi-polygon: {summary['boundary_only_match_count']} / {summary['multi_polygon_match_count']}",
            "",
            f"**Classification:** {decision['classification']}",
            "",
            "No pollen-limitation effect size or treatment outcome was read by this preflight.",
        ]
    )
    (output_dir / "RESULT_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    typer.echo(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    app()
