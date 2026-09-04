"""Build an outcome-blind GIFT mainland source-pool universe for PR138.

The builder retrieves complete native Tracheophyta checklists from GIFT, restricts the
source universe to suitable Mainland entities, and attaches representative coordinates
and the same CHELSA/elevation variables used by the frozen island climate PCA.

No island floral outcome is read by this module.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests
from shapely.geometry import shape

API = "https://gift.uni-goettingen.de/api/extended/index3.2.php"
GEOJSON = "https://gift.uni-goettingen.de/geojson/geojson_smaller3.2/{entity_id}.geojson"
CHELSA_LAYERS = {
    "bio1": "CHELSA_bio1_1981-2010_V.2.1",
    "bio4": "CHELSA_bio4_1981-2010_V.2.1",
    "bio5": "CHELSA_bio5_1981-2010_V.2.1",
    "bio6": "CHELSA_bio6_1981-2010_V.2.1",
    "bio12": "CHELSA_bio12_1981-2010_V.2.1",
    "bio15": "CHELSA_bio15_1981-2010_V.2.1",
    "bio18": "CHELSA_bio18_1981-2010_V.2.1",
    "bio19": "CHELSA_bio19_1981-2010_V.2.1",
    "elevation": "mn30_grd",
}
NATIVE_COMPLETE_SUBSETS = {
    "all",
    "native",
    "native and naturalized",
    "native and historically introduced",
}


def _num(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _int(value: object) -> int | None:
    value = _num(value)
    return int(value) if math.isfinite(value) else None


def _text(value: object) -> str:
    return "" if value is None else str(value).strip()


def _request_json(url: str, *, retries: int = 4, timeout: int = 120) -> Any:
    error: Exception | None = None
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except Exception as exc:  # network boundary: retry, then fail closed upstream
            error = exc
            if attempt + 1 < retries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GIFT request failed after {retries} attempts: {url}: {error}")


def _metadata() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    lists = _request_json(f"{API}?query=lists")
    taxonomy = _request_json(f"{API}?query=taxonomy")
    if not isinstance(lists, list) or not lists:
        raise RuntimeError("GIFT lists metadata missing or malformed")
    if not isinstance(taxonomy, list) or not taxonomy:
        raise RuntimeError("GIFT taxonomy metadata missing or malformed")
    return lists, taxonomy


def _tracheophyta_contract(taxonomy: list[dict[str, Any]]) -> tuple[int, set[int]]:
    target = [row for row in taxonomy if _text(row.get("taxon_name")) == "Tracheophyta"]
    if len(target) != 1:
        raise RuntimeError(f"Expected exactly one Tracheophyta taxon, found {len(target)}")
    target_id = _int(target[0].get("taxon_ID"))
    target_lft = _num(target[0].get("lft"))
    target_rgt = _num(target[0].get("rgt"))
    if target_id is None or not all(math.isfinite(x) for x in (target_lft, target_rgt)):
        raise RuntimeError("Tracheophyta taxonomy row lacks ID or nested-set coordinates")
    complete_ids: set[int] = set()
    for row in taxonomy:
        taxon_id = _int(row.get("taxon_ID"))
        lft = _num(row.get("lft"))
        rgt = _num(row.get("rgt"))
        if taxon_id is None or not all(math.isfinite(x) for x in (lft, rgt)):
            continue
        if lft <= target_lft and rgt >= target_rgt:
            complete_ids.add(taxon_id)
    return target_id, complete_ids


def eligible_mainland_lists(
    lists: list[dict[str, Any]], taxonomy: list[dict[str, Any]]
) -> tuple[pd.DataFrame, int]:
    target_id, complete_taxon_ids = _tracheophyta_contract(taxonomy)
    rows: list[dict[str, Any]] = []
    for row in lists:
        if _text(row.get("entity_class")) != "Mainland":
            continue
        if _int(row.get("suit_geo")) != 1:
            continue
        if _int(row.get("native_indicated")) != 1:
            continue
        if _int(row.get("restricted")) not in (0, None):
            continue
        if _int(row.get("taxon_ID")) not in complete_taxon_ids:
            continue
        if _text(row.get("subset")).lower() not in NATIVE_COMPLETE_SUBSETS:
            continue
        list_id = _int(row.get("list_ID"))
        entity_id = _int(row.get("entity_ID"))
        if list_id is None or entity_id is None:
            continue
        rows.append(
            {
                "list_ID": list_id,
                "entity_ID": entity_id,
                "geo_entity": _text(row.get("geo_entity")),
                "ref_ID": _int(row.get("ref_ID")),
                "reference_type": _text(row.get("type")),
                "subset": _text(row.get("subset")),
                "taxon_ID": _int(row.get("taxon_ID")),
            }
        )
    out = pd.DataFrame(rows).drop_duplicates("list_ID")
    if out.empty:
        raise RuntimeError("No GIFT mainland lists satisfy the frozen source-pool contract")
    return out.sort_values(["entity_ID", "list_ID"]).reset_index(drop=True), target_id


def _fetch_one_checklist(list_id: int, target_id: int) -> tuple[int, list[dict[str, Any]]]:
    url = (
        f"{API}?query=checklists&listid={int(list_id)}&taxonid={int(target_id)}"
        "&namesmatched=0&filter=native"
    )
    payload = _request_json(url)
    if payload is None:
        payload = []
    if not isinstance(payload, list):
        raise RuntimeError(f"Malformed checklist payload for list {list_id}")
    return int(list_id), payload


def fetch_native_flora(eligible: pd.DataFrame, target_id: int, *, workers: int = 8) -> tuple[pd.DataFrame, pd.DataFrame]:
    list_to_entity = eligible.set_index("list_ID")["entity_ID"].astype(int).to_dict()
    rows: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_fetch_one_checklist, int(list_id), target_id): int(list_id)
            for list_id in eligible["list_ID"].astype(int).unique()
        }
        for future in as_completed(futures):
            list_id = futures[future]
            try:
                _, payload = future.result()
                kept = 0
                for record in payload:
                    if _int(record.get("native")) != 1:
                        continue
                    if _int(record.get("questionable")) == 1:
                        continue
                    if _int(record.get("quest_native")) == 1:
                        continue
                    species = _text(record.get("work_species"))
                    if not species:
                        continue
                    rows.append(
                        {
                            "entity_ID": int(list_to_entity[list_id]),
                            "list_ID": list_id,
                            "work_species": species,
                            "work_ID": _int(record.get("work_ID")),
                        }
                    )
                    kept += 1
                audit.append({"list_ID": list_id, "status": "ok", "n_native_species": kept, "error": ""})
            except Exception as exc:
                audit.append({"list_ID": list_id, "status": "failed", "n_native_species": 0, "error": str(exc)[:500]})
    flora = pd.DataFrame(rows)
    if flora.empty:
        raise RuntimeError("All GIFT source checklists failed or yielded no strict native species")
    flora = flora.drop_duplicates(["entity_ID", "work_species"]).sort_values(["entity_ID", "work_species"])
    return flora.reset_index(drop=True), pd.DataFrame(audit).sort_values("list_ID").reset_index(drop=True)


def _fetch_one_region(entity_id: int) -> dict[str, Any]:
    payload = _request_json(GEOJSON.format(entity_id=int(entity_id)))
    features = payload.get("features", []) if isinstance(payload, dict) else []
    if not features:
        raise RuntimeError(f"No GeoJSON feature for GIFT entity {entity_id}")
    feature = features[0]
    properties = feature.get("properties") or {}
    lon = _num(properties.get("point_x"))
    lat = _num(properties.get("point_y"))
    if not (math.isfinite(lon) and math.isfinite(lat)):
        geometry = shape(feature["geometry"])
        centroid = geometry.centroid
        lon, lat = float(centroid.x), float(centroid.y)
    return {
        "entity_ID": int(entity_id),
        "source_longitude": lon,
        "source_latitude": lat,
        "source_area_km2": _num(properties.get("area")),
        "source_entity_type": _text(properties.get("entity_type")),
    }


def fetch_region_points(entity_ids: list[int], *, workers: int = 8) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_fetch_one_region, entity_id): entity_id for entity_id in entity_ids}
        for future in as_completed(futures):
            entity_id = futures[future]
            try:
                rows.append(future.result())
            except Exception as exc:
                failures.append({"entity_ID": entity_id, "error": str(exc)[:500]})
    return pd.DataFrame(rows).sort_values("entity_ID"), pd.DataFrame(failures)


def _fetch_environment_layer(variable: str, layer_name: str) -> pd.DataFrame:
    url = f"{API}?query=geoentities_env_raster&layername={layer_name}&sumstat=mean"
    payload = _request_json(url)
    if not isinstance(payload, list):
        raise RuntimeError(f"Malformed environment payload for {layer_name}")
    rows = []
    for row in payload:
        entity_id = _int(row.get("entity_ID"))
        value = _num(row.get("mean"))
        if entity_id is not None and math.isfinite(value):
            rows.append({"entity_ID": entity_id, variable: value})
    return pd.DataFrame(rows).drop_duplicates("entity_ID")


def fetch_environment() -> pd.DataFrame:
    merged: pd.DataFrame | None = None
    for variable, layer_name in CHELSA_LAYERS.items():
        frame = _fetch_environment_layer(variable, layer_name)
        merged = frame if merged is None else merged.merge(frame, on="entity_ID", how="outer")
    return merged if merged is not None else pd.DataFrame(columns=["entity_ID"])


def project_frozen_climate_pcs(environment: pd.DataFrame, pca_loadings: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"variable", "mean", "sd", "climate_pc1_loading", "climate_pc2_loading", "climate_pc3_loading", "climate_pc4_loading"}
    missing = required - set(pca_loadings.columns)
    if missing:
        raise RuntimeError(f"Frozen PCA loadings missing columns: {sorted(missing)}")
    load = pca_loadings.loc[pca_loadings["variable"].isin(CHELSA_LAYERS)].copy()
    if set(load["variable"].astype(str)) != set(CHELSA_LAYERS):
        raise RuntimeError("Frozen PCA loadings do not contain all nine source-pool variables")
    out = environment.copy()
    audit_rows: list[dict[str, Any]] = []
    z_columns: dict[str, np.ndarray] = {}
    for row in load.itertuples(index=False):
        variable = str(row.variable)
        values = pd.to_numeric(out[variable], errors="coerce").to_numpy(float)
        mean = float(row.mean)
        sd = float(row.sd)
        z = (values - mean) / sd
        z_columns[variable] = z
        finite = z[np.isfinite(z)]
        audit_rows.append(
            {
                "variable": variable,
                "n_source_regions": int(len(finite)),
                "source_value_median": float(np.nanmedian(values)) if np.isfinite(values).any() else np.nan,
                "frozen_island_mean": mean,
                "frozen_island_sd": sd,
                "median_abs_frozen_z": float(np.median(np.abs(finite))) if len(finite) else np.nan,
                "p95_abs_frozen_z": float(np.quantile(np.abs(finite), 0.95)) if len(finite) else np.nan,
            }
        )
    for pc in range(1, 5):
        score = np.zeros(len(out), dtype=float)
        complete = np.ones(len(out), dtype=bool)
        for row in load.itertuples(index=False):
            variable = str(row.variable)
            z = z_columns[variable]
            complete &= np.isfinite(z)
            score += np.nan_to_num(z, nan=0.0) * float(getattr(row, f"climate_pc{pc}_loading"))
        score[~complete] = np.nan
        out[f"climate_pc{pc}"] = score
    return out, pd.DataFrame(audit_rows)


def run(output_dir: Path, pca_loadings_csv: Path, *, workers: int = 8) -> dict[str, Any]:
    lists, taxonomy = _metadata()
    eligible, target_id = eligible_mainland_lists(lists, taxonomy)
    flora, checklist_audit = fetch_native_flora(eligible, target_id, workers=workers)
    entity_ids = sorted(flora["entity_ID"].astype(int).unique())
    points, geometry_failures = fetch_region_points(entity_ids, workers=workers)
    environment = fetch_environment()
    environment, climate_audit = project_frozen_climate_pcs(environment, pd.read_csv(pca_loadings_csv))

    entity_meta = (
        eligible.groupby("entity_ID", as_index=False)
        .agg(geo_entity=("geo_entity", "first"), n_eligible_lists=("list_ID", "nunique"))
    )
    native_counts = flora.groupby("entity_ID")["work_species"].nunique().rename("n_native_species").reset_index()
    regions = entity_meta.merge(native_counts, on="entity_ID", how="inner")
    regions = regions.merge(points, on="entity_ID", how="left").merge(environment, on="entity_ID", how="left")

    output_dir.mkdir(parents=True, exist_ok=True)
    eligible.to_csv(output_dir / "gift_eligible_mainland_lists.csv", index=False)
    flora.to_csv(output_dir / "gift_native_mainland_flora.csv.gz", index=False)
    checklist_audit.to_csv(output_dir / "gift_checklist_fetch_audit.csv", index=False)
    geometry_failures.to_csv(output_dir / "gift_geometry_fetch_failures.csv", index=False)
    climate_audit.to_csv(output_dir / "gift_frozen_pca_scale_audit.csv", index=False)
    regions.to_csv(output_dir / "gift_mainland_source_regions.csv", index=False)

    n_failed_lists = int((checklist_audit["status"] == "failed").sum())
    n_regions_coords = int(regions[["source_longitude", "source_latitude"]].notna().all(axis=1).sum())
    n_regions_climate = int(regions[[f"climate_pc{i}" for i in range(1, 5)]].notna().all(axis=1).sum())
    manifest = {
        "contract": "chapter1_pr138_gift_source_pool_v1",
        "gift_version": "3.2",
        "target_taxon": "Tracheophyta",
        "entity_class": "Mainland",
        "n_eligible_lists": int(eligible["list_ID"].nunique()),
        "n_failed_list_requests": n_failed_lists,
        "n_source_entities_with_native_flora": int(len(entity_ids)),
        "n_unique_native_species_rows": int(len(flora)),
        "n_source_entities_with_coordinates": n_regions_coords,
        "n_source_entities_with_complete_frozen_climate_pcs": n_regions_climate,
        "outcome_blind": True,
        "island_trait_outcomes_used": False,
        "claim_boundary": "External mainland source universe only; no island is assigned a source in this stage.",
    }
    if n_regions_coords < 100:
        raise RuntimeError(f"Only {n_regions_coords} source regions have coordinates; source proxy fails closed")
    (output_dir / "gift_source_pool_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pca-loadings-csv", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir, args.pca_loadings_csv, workers=args.workers), indent=2))


if __name__ == "__main__":
    main()
