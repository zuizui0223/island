from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, help="Build and validate Island Plant Database bundles.")

CORE_TABLES = ("islands", "taxa", "island_taxa", "evidence")
OPTIONAL_TABLES = ("traits",)
RIGHTS_COLUMNS = {
    "object_type",
    "object_id",
    "source_id",
    "source_url",
    "source_license",
    "rights_status",
    "rights_evidence",
    "redistributable_value",
    "redistributable_provenance",
    "review_status",
}


def load_contract(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported Island Plant Database schema_version")
    if payload.get("database", {}).get("primary_object") != "island_x_taxon":
        raise ValueError("database.primary_object must be island_x_taxon")
    tables = payload.get("canonical_tables", {})
    missing = sorted(set(CORE_TABLES + OPTIONAL_TABLES) - set(tables))
    if missing:
        raise ValueError(f"contract missing canonical tables: {missing}")
    return payload


def _read_table(bundle_dir: Path, name: str) -> pd.DataFrame:
    path = bundle_dir / f"{name}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, dtype=str).fillna("")


def _validate_table(name: str, frame: pd.DataFrame, spec: dict[str, Any]) -> None:
    required = list(spec.get("required_columns", []))
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{name} missing required columns: {missing}")

    keys = list(spec.get("key", []))
    if not keys:
        raise ValueError(f"{name} has no declared key")
    if frame[keys].eq("").any(axis=None):
        raise ValueError(f"{name} has blank primary-key values")
    if frame.duplicated(keys).any():
        raise ValueError(f"{name} has duplicate primary keys: {keys}")

    controlled = spec.get("controlled_values", {})
    for column, allowed_values in controlled.items():
        observed = set(frame[column].astype(str))
        invalid = sorted(observed - set(allowed_values))
        if invalid:
            raise ValueError(f"{name}.{column} contains unsupported values: {invalid}")


def _require_nonnegative_integer(
    frame: pd.DataFrame,
    column: str,
    *,
    allow_blank: bool,
) -> None:
    values = frame[column].astype(str)
    if allow_blank:
        values = values[values.ne("")]
    elif values.eq("").any():
        raise ValueError(f"{column} contains blank values")
    if values.empty:
        return
    parsed = pd.to_numeric(values, errors="coerce")
    if parsed.isna().any() or (parsed < 0).any() or ((parsed % 1) != 0).any():
        raise ValueError(f"{column} must contain non-negative integers")


def _geometry_sha256(geometry: Any) -> str:
    return hashlib.sha256(geometry.wkb).hexdigest()


def build_islands_core(
    islands_gpkg: Path,
    source_policy_json: Path,
    contract_path: Path,
) -> pd.DataFrame:
    contract = load_contract(contract_path)
    frame = gpd.read_file(islands_gpkg, layer="islands")
    policy = json.loads(source_policy_json.read_text(encoding="utf-8"))

    required = {
        "island_id",
        "source_label",
        "parent_feature_id",
        "island_name",
        "area_km2",
        "geometry",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"prepared island geometry missing columns: {missing}")
    if frame.crs is None:
        raise ValueError("prepared island geometry has no CRS")
    if frame["island_id"].duplicated().any():
        raise ValueError("prepared island geometry has duplicate island_id values")

    backend = str(policy.get("source_backend", ""))
    if backend == "gshhg":
        source_license = "LGPL-3.0-or-later"
        geometry_source = "GSHHG"
    elif backend == "natural_earth_10m_fallback":
        source_license = "PUBLIC-DOMAIN"
        geometry_source = "Natural Earth 10m"
    else:
        raise ValueError(f"unsupported island source backend: {backend!r}")

    projected = frame.to_crs(6933)
    centroids = gpd.GeoSeries(projected.geometry.centroid, crs=6933).to_crs(4326)

    source_labels = frame["source_label"].astype(str)
    source_feature_ids = frame["parent_feature_id"].astype(str)
    output = pd.DataFrame(
        {
            "island_id": frame["island_id"].astype(str),
            "source_island_id": source_labels + ":" + source_feature_ids,
            "island_name": frame["island_name"].fillna("").astype(str),
            "archipelago": "",
            "country_or_territory": "",
            "area_km2": frame["area_km2"].astype(float),
            "centroid_lat": centroids.y,
            "centroid_lon": centroids.x,
            "geometry_source": geometry_source,
            "geometry_version": source_labels,
            "geometry_sha256": [_geometry_sha256(value) for value in frame.geometry],
            "release_status": "redistributable",
            "source_license": source_license,
        }
    ).sort_values("island_id", kind="stable")
    output = output.reset_index(drop=True)
    _validate_table("islands", output, contract["canonical_tables"]["islands"])
    if output["geometry_sha256"].str.fullmatch(r"[0-9a-f]{64}").eq(False).any():
        raise ValueError("islands.geometry_sha256 must contain full SHA-256 digests")
    if (output["area_km2"] <= 0).any():
        raise ValueError("islands.area_km2 must be positive")
    return output


def validate_bundle(bundle_dir: Path, contract_path: Path) -> dict[str, Any]:
    contract = load_contract(contract_path)
    specs = contract["canonical_tables"]

    tables: dict[str, pd.DataFrame] = {}
    for name in CORE_TABLES:
        frame = _read_table(bundle_dir, name)
        _validate_table(name, frame, specs[name])
        tables[name] = frame

    traits_path = bundle_dir / "traits.csv"
    if traits_path.exists():
        traits = pd.read_csv(traits_path, dtype=str).fillna("")
        _validate_table("traits", traits, specs["traits"])
        tables["traits"] = traits

    islands = tables["islands"]
    taxa = tables["taxa"]
    island_taxa = tables["island_taxa"]
    evidence = tables["evidence"]

    island_ids = set(islands["island_id"])
    taxon_ids = set(taxa["taxon_id"])
    island_taxon_keys = set(
        zip(island_taxa["island_id"], island_taxa["taxon_id"], strict=True)
    )

    missing_islands = sorted(set(island_taxa["island_id"]) - island_ids)
    missing_taxa = sorted(set(island_taxa["taxon_id"]) - taxon_ids)
    if missing_islands:
        raise ValueError(f"island_taxa references unknown islands: {missing_islands[:10]}")
    if missing_taxa:
        raise ValueError(f"island_taxa references unknown taxa: {missing_taxa[:10]}")

    _require_nonnegative_integer(island_taxa, "occurrence_record_count", allow_blank=False)
    _require_nonnegative_integer(island_taxa, "specimen_record_count", allow_blank=True)
    _require_nonnegative_integer(island_taxa, "evidence_count", allow_blank=False)

    evidence_islands = set(evidence.loc[evidence["island_id"].ne(""), "island_id"])
    evidence_taxa = set(evidence.loc[evidence["taxon_id"].ne(""), "taxon_id"])
    if evidence_islands - island_ids:
        raise ValueError(
            f"evidence references unknown islands: {sorted(evidence_islands - island_ids)[:10]}"
        )
    if evidence_taxa - taxon_ids:
        raise ValueError(
            f"evidence references unknown taxa: {sorted(evidence_taxa - taxon_ids)[:10]}"
        )

    pair_evidence = evidence[evidence["island_id"].ne("") & evidence["taxon_id"].ne("")]
    unknown_pairs = sorted(
        set(zip(pair_evidence["island_id"], pair_evidence["taxon_id"], strict=True))
        - island_taxon_keys
    )
    if unknown_pairs:
        raise ValueError(f"evidence references unknown island × taxon pairs: {unknown_pairs[:10]}")

    if "traits" in tables:
        traits = tables["traits"]
        trait_taxa = set(traits["taxon_id"])
        if trait_taxa - taxon_ids:
            raise ValueError(f"traits reference unknown taxa: {sorted(trait_taxa - taxon_ids)[:10]}")
        trait_islands = set(traits.loc[traits["island_id"].ne(""), "island_id"])
        if trait_islands - island_ids:
            raise ValueError(
                f"traits reference unknown islands: {sorted(trait_islands - island_ids)[:10]}"
            )
        evidence_ids = set(evidence["evidence_id"])
        missing_evidence = sorted(set(traits["evidence_id"]) - evidence_ids)
        if missing_evidence:
            raise ValueError(f"traits reference unknown evidence: {missing_evidence[:10]}")

    rights_path = bundle_dir / "RIGHTS_LEDGER.csv"
    if rights_path.exists():
        rights = pd.read_csv(rights_path, dtype=str).fillna("")
        missing_rights = sorted(RIGHTS_COLUMNS - set(rights.columns))
        if missing_rights:
            raise ValueError(f"RIGHTS_LEDGER missing columns: {missing_rights}")
        allowed = {"redistributable", "reference_only", "review_required"}
        invalid = sorted(set(rights["rights_status"]) - allowed)
        if invalid:
            raise ValueError(f"RIGHTS_LEDGER contains unsupported rights_status: {invalid}")

    return {
        "database_id": contract["database"]["database_id"],
        "version": contract["database"]["version"],
        "islands": int(len(islands)),
        "taxa": int(len(taxa)),
        "island_taxa": int(len(island_taxa)),
        "evidence": int(len(evidence)),
        "traits": int(len(tables.get("traits", []))),
        "rights_ledger_present": rights_path.exists(),
        "valid": True,
    }


@app.command("export-islands")
def export_islands_command(
    islands_gpkg: Path = typer.Option(..., "--islands-gpkg", exists=True, dir_okay=False),
    source_policy_json: Path = typer.Option(
        ...,
        "--source-policy",
        exists=True,
        dir_okay=False,
    ),
    output_csv: Path = typer.Option(..., "--output-csv"),
    contract_path: Path = typer.Option(
        Path("config/island_plant_database_v2.yml"),
        "--contract",
        exists=True,
        dir_okay=False,
    ),
) -> None:
    frame = build_islands_core(islands_gpkg, source_policy_json, contract_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_csv, index=False)
    typer.echo(f"Wrote {len(frame)} Database 2.0 island rows to {output_csv}")


@app.command("validate")
def validate_command(
    bundle_dir: Path = typer.Option(..., "--bundle-dir", exists=True, file_okay=False),
    contract_path: Path = typer.Option(
        Path("config/island_plant_database_v2.yml"),
        "--contract",
        exists=True,
        dir_okay=False,
    ),
) -> None:
    report = validate_bundle(bundle_dir, contract_path)
    typer.echo(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
