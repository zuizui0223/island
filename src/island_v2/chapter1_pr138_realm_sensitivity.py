"""Formal biogeographic-realm sensitivity for PR138 syndrome results.

The primary PR138 contrast uses broad northern/tropical/southern contexts. This module
asks whether that result is merely an artefact of those latitude-based labels by
reassigning island points to published terrestrial biogeographic realms and rerunning
the same predeclared syndrome analysis.

Realm assignment is fail-closed: only point-in-polygon intersections are accepted.
Unresolved island points are retained in the assignment audit and are never assigned to
the nearest realm. The realm analysis is a sensitivity analysis, not a replacement for
the full-island primary analysis.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_syndrome_analysis import run_syndrome_analysis

app = typer.Typer(add_completion=False, no_args_is_help=True)

REALM_ALIASES = {
    "nearctic": "Nearctic",
    "palearctic": "Palearctic",
    "palaearctic": "Palearctic",
    "neotropical": "Neotropical",
    "neotropic": "Neotropical",
    "afrotropical": "Afrotropical",
    "afrotropic": "Afrotropical",
    "indo-malay": "Indomalayan",
    "indo-malay realm": "Indomalayan",
    "indomalayan": "Indomalayan",
    "australasia": "Australasian",
    "australasian": "Australasian",
    "oceania": "Oceanian",
    "oceanian": "Oceanian",
    "antarctic": "Antarctic",
}


def _realm_column(frame: gpd.GeoDataFrame) -> str:
    exact = ["REALM", "Realm", "realm", "REALM_NAME", "realm_name"]
    for column in exact:
        if column in frame.columns:
            return column
    candidates = [c for c in frame.columns if "realm" in str(c).lower()]
    if not candidates:
        raise typer.BadParameter("ecoregion layer has no realm attribute")
    return str(candidates[0])


def _normalize_realm(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    normalized = REALM_ALIASES.get(text.lower())
    return normalized if normalized else text


def assign_biogeographic_realms(
    covariates: pd.DataFrame,
    ecoregions: gpd.GeoDataFrame,
) -> pd.DataFrame:
    required = {"island_id", "island_latitude", "island_longitude"}
    missing = required - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"island covariates missing columns: {sorted(missing)}")
    realm_column = _realm_column(ecoregions)

    islands = covariates[["island_id", "island_latitude", "island_longitude"]].copy()
    islands["island_latitude"] = pd.to_numeric(islands["island_latitude"], errors="coerce")
    islands["island_longitude"] = pd.to_numeric(islands["island_longitude"], errors="coerce")
    islands = islands.drop_duplicates("island_id")
    valid = islands.dropna(subset=["island_latitude", "island_longitude"]).copy()
    points = gpd.GeoDataFrame(
        valid,
        geometry=gpd.points_from_xy(valid["island_longitude"], valid["island_latitude"]),
        crs="EPSG:4326",
    )

    realms = ecoregions[[realm_column, "geometry"]].copy()
    if realms.crs is None:
        realms = realms.set_crs("EPSG:4326")
    else:
        realms = realms.to_crs("EPSG:4326")
    realms["biogeographic_realm"] = realms[realm_column].map(_normalize_realm)
    realms = realms.loc[realms["biogeographic_realm"].ne("")].copy()

    joined = gpd.sjoin(
        points,
        realms[["biogeographic_realm", "geometry"]],
        how="left",
        predicate="intersects",
    ).drop(columns=["index_right"], errors="ignore")

    rows: list[dict[str, Any]] = []
    for island_id, group in joined.groupby("island_id", sort=True):
        realm_values = sorted({str(x) for x in group["biogeographic_realm"].dropna() if str(x)})
        if len(realm_values) == 1:
            realm = realm_values[0]
            status = "point_intersection"
        elif len(realm_values) > 1:
            realm = ""
            status = "realm_boundary_conflict"
        else:
            realm = ""
            status = "unresolved_no_intersection"
        first = group.iloc[0]
        rows.append(
            {
                "island_id": str(island_id),
                "island_latitude": float(first["island_latitude"]),
                "island_longitude": float(first["island_longitude"]),
                "biogeographic_realm": realm,
                "realm_assignment_status": status,
                "n_intersected_realms": len(realm_values),
            }
        )

    assigned = pd.DataFrame(rows)
    missing_ids = set(islands["island_id"].astype(str)) - set(assigned["island_id"].astype(str))
    if missing_ids:
        missing_rows = islands.loc[islands["island_id"].astype(str).isin(missing_ids)].copy()
        missing_rows["biogeographic_realm"] = ""
        missing_rows["realm_assignment_status"] = "unresolved_missing_coordinates"
        missing_rows["n_intersected_realms"] = 0
        assigned = pd.concat(
            [
                assigned,
                missing_rows.rename(columns={
                    "island_latitude": "island_latitude",
                    "island_longitude": "island_longitude",
                })[
                    [
                        "island_id",
                        "island_latitude",
                        "island_longitude",
                        "biogeographic_realm",
                        "realm_assignment_status",
                        "n_intersected_realms",
                    ]
                ],
            ],
            ignore_index=True,
        )
    return assigned.sort_values("island_id").reset_index(drop=True)


def run_realm_sensitivity(
    trait_ledger: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    ecoregions: gpd.GeoDataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
    realm_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assignment = assign_biogeographic_realms(covariates, ecoregions)
    augmented = covariates.merge(
        assignment[["island_id", "biogeographic_realm", "realm_assignment_status"]],
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    sensitivity_config = dict(pattern_config)
    sensitivity_config["context_column"] = "biogeographic_realm"
    sensitivity_config["contexts"] = [str(x) for x in realm_config["realms"]]
    species, islands, slopes, within, between = run_syndrome_analysis(
        trait_ledger,
        status_flora,
        augmented,
        sensitivity_config,
        syndrome_config,
    )
    return assignment, species, islands, slopes, within, between


@app.command("run")
def run(
    trait_ledger_csv: Path = typer.Option(..., exists=True),
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    ecoregions_path: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    realm_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    syndrome_config = yaml.safe_load(syndrome_config_path.read_text(encoding="utf-8"))
    realm_config = yaml.safe_load(realm_config_path.read_text(encoding="utf-8"))
    assignment, species, islands, slopes, within, between = run_realm_sensitivity(
        pd.read_csv(trait_ledger_csv),
        pd.read_csv(status_flora_csv),
        pd.read_csv(covariates_csv),
        gpd.read_file(ecoregions_path),
        pattern_config,
        syndrome_config,
        realm_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    assignment.to_csv(output_dir / "island_biogeographic_realm_assignment.csv", index=False)
    species.to_csv(output_dir / "realm_species_syndrome_concordance.csv.gz", index=False)
    islands.to_csv(output_dir / "realm_island_syndrome_scores.csv.gz", index=False)
    slopes.to_csv(output_dir / "realm_syndrome_distance_slopes.csv", index=False)
    within.to_csv(output_dir / "realm_syndrome_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "realm_syndrome_between_omnibus.csv", index=False)

    resolved = assignment["biogeographic_realm"].fillna("").ne("")
    manifest = {
        "contract": "chapter1_pr138_realm_sensitivity_v1",
        "realm_source": realm_config["dataset"],
        "assignment_rule": "point_intersection_only",
        "nearest_realm_fallback": False,
        "n_islands_total": int(len(assignment)),
        "n_islands_realm_resolved": int(resolved.sum()),
        "n_islands_realm_unresolved": int((~resolved).sum()),
        "sensitivity_not_primary_replacement": True,
    }
    (output_dir / "realm_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
