"""Match source-backed GIFT island floristic status to the frozen GSHHG flora.

Primary matching is intentionally conservative:
- only GIFT regions whose entity_class is exactly ``Island`` are eligible;
- geometry matching requires strong bidirectional overlap and an unambiguous best
  GSHHG polygon;
- species matching to the frozen island flora is exact after whitespace cleanup;
- GIFT taxonomic/status uncertainty fails closed to ``unresolved``;
- focal-island endemism uses GIFT ``endemic_list`` only when ``end_list`` says
  list-level endemism is explicitly indicated. ``endemic_ref`` is retained for
  audit but never substituted for focal-island endemism.

The output is a source-level status ledger. Multiple independent GIFT lists are
preserved; downstream ``flora_status_support.collapse_status_ledger`` resolves
agreement and converts conflicts to unresolved.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import typer

from island_v2.flora_status_support import _text, normalize_island_species

app = typer.Typer(add_completion=False, no_args_is_help=True)

EQUAL_AREA_CRS = "EPSG:6933"


def _flag(value: object) -> bool:
    text = str(value if value is not None else "").strip().casefold()
    return text in {"1", "1.0", "true", "yes", "y"}


def _known_zero(value: object) -> bool:
    text = str(value if value is not None else "").strip().casefold()
    return text in {"0", "0.0", "false", "no", "n"}


def _column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise typer.BadParameter(f"none of the required columns are present: {list(candidates)}")


def prepare_list_metadata(list_meta: pd.DataFrame) -> pd.DataFrame:
    required = {"list_ID", "entity_ID", "entity_class"}
    missing = required - set(list_meta.columns)
    if missing:
        raise typer.BadParameter(f"GIFT list metadata missing columns: {sorted(missing)}")
    frame = list_meta.copy()
    frame["entity_class"] = _text(frame["entity_class"])
    frame = frame.loc[frame["entity_class"].eq("Island")].copy()
    if "suit_geo" in frame.columns:
        frame = frame.loc[frame["suit_geo"].map(_flag)].copy()
    if "restricted" in frame.columns:
        frame = frame.loc[~frame["restricted"].map(_flag)].copy()
    frame["list_ID"] = _text(frame["list_ID"])
    frame["entity_ID"] = _text(frame["entity_ID"])
    if "ref_ID" in frame.columns:
        frame["ref_ID"] = _text(frame["ref_ID"])
    else:
        frame["ref_ID"] = ""
    frame = frame.loc[frame["list_ID"].ne("") & frame["entity_ID"].ne("")]
    return frame.drop_duplicates(["list_ID", "entity_ID", "ref_ID"]).reset_index(drop=True)


def match_gift_islands_to_gshhg(
    gift_shapes: gpd.GeoDataFrame,
    gshhg_islands: gpd.GeoDataFrame,
    *,
    min_bidirectional_overlap: float = 0.70,
    min_iou: float = 0.50,
    ambiguity_margin: float = 0.20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return accepted GIFT->GSHHG matches plus all spatial candidates for audit."""
    if "entity_ID" not in gift_shapes.columns:
        raise typer.BadParameter("GIFT shapes missing entity_ID")
    if "island_id" not in gshhg_islands.columns:
        raise typer.BadParameter("GSHHG islands missing island_id")
    if gift_shapes.crs is None or gshhg_islands.crs is None:
        raise typer.BadParameter("both spatial inputs require CRS metadata")

    gift = gift_shapes.copy()
    if "entity_class" in gift.columns:
        gift = gift.loc[_text(gift["entity_class"]).eq("Island")].copy()
    gift["entity_ID"] = _text(gift["entity_ID"])
    gift = gift.loc[gift["entity_ID"].ne("") & gift.geometry.notna()].copy()
    gshhg = gshhg_islands[["island_id", "geometry"]].copy()
    gshhg["island_id"] = _text(gshhg["island_id"])
    gshhg = gshhg.loc[gshhg["island_id"].ne("") & gshhg.geometry.notna()].copy()

    gift = gift.to_crs(EQUAL_AREA_CRS)
    gshhg = gshhg.to_crs(EQUAL_AREA_CRS)
    gift["geometry"] = gift.geometry.make_valid()
    gshhg["geometry"] = gshhg.geometry.make_valid()
    sindex = gshhg.sindex

    candidate_rows: list[dict[str, Any]] = []
    match_rows: list[dict[str, Any]] = []
    for row in gift.itertuples(index=False):
        entity_id = str(row.entity_ID)
        geom = row.geometry
        gift_area = float(geom.area)
        if not np.isfinite(gift_area) or gift_area <= 0:
            match_rows.append(
                {
                    "entity_ID": entity_id,
                    "island_id": "",
                    "spatial_match_status": "invalid_gift_geometry",
                }
            )
            continue
        candidate_idx = list(sindex.query(geom, predicate="intersects"))
        scored: list[dict[str, Any]] = []
        for idx in candidate_idx:
            target = gshhg.iloc[int(idx)]
            inter = geom.intersection(target.geometry)
            inter_area = float(inter.area)
            if inter_area <= 0:
                continue
            target_area = float(target.geometry.area)
            union_area = gift_area + target_area - inter_area
            gift_coverage = inter_area / gift_area
            gshhg_coverage = inter_area / target_area if target_area > 0 else 0.0
            iou = inter_area / union_area if union_area > 0 else 0.0
            record = {
                "entity_ID": entity_id,
                "island_id": str(target.island_id),
                "gift_coverage": gift_coverage,
                "gshhg_coverage": gshhg_coverage,
                "bidirectional_overlap": min(gift_coverage, gshhg_coverage),
                "iou": iou,
            }
            scored.append(record)
            candidate_rows.append(record)
        scored.sort(key=lambda item: (item["iou"], item["bidirectional_overlap"]), reverse=True)
        if not scored:
            match_rows.append(
                {
                    "entity_ID": entity_id,
                    "island_id": "",
                    "spatial_match_status": "no_intersection",
                }
            )
            continue
        top = scored[0]
        second_iou = scored[1]["iou"] if len(scored) > 1 else 0.0
        strong = (
            top["bidirectional_overlap"] >= min_bidirectional_overlap
            and top["iou"] >= min_iou
        )
        ambiguous = len(scored) > 1 and (top["iou"] - second_iou) < ambiguity_margin
        status = "accepted" if strong and not ambiguous else (
            "ambiguous" if ambiguous else "insufficient_overlap"
        )
        match_rows.append(
            {
                **top,
                "second_best_iou": second_iou,
                "spatial_match_status": status,
            }
        )
    return pd.DataFrame(match_rows), pd.DataFrame(candidate_rows)


def _taxon_safe(row: pd.Series) -> tuple[bool, str]:
    if "questionable" in row.index and _flag(row.get("questionable")):
        return False, "questionable_occurrence"
    if "matched" in row.index and str(row.get("matched", "")).strip() and not _flag(row.get("matched")):
        return False, "taxon_not_matched"
    if "resolved" in row.index and str(row.get("resolved", "")).strip() and not _flag(row.get("resolved")):
        return False, "taxon_not_resolved"
    for column in ("cf_genus", "cf_species", "aff_species"):
        if column in row.index and _flag(row.get(column)):
            return False, "uncertain_source_name"
    if "subtaxon" in row.index and str(row.get("subtaxon", "")).strip():
        return False, "infraspecific_source_scope"
    return True, ""


def classify_gift_status(row: pd.Series) -> tuple[str, str, str]:
    """Return origin, focal-island endemism, and an audit note."""
    safe, reason = _taxon_safe(row)
    if not safe:
        return "unresolved", "unresolved", reason

    native_available = _flag(row.get("native_indicated"))
    naturalized_available = _flag(row.get("natural_indicated"))
    native = native_available and _flag(row.get("native")) and not _flag(row.get("quest_native"))
    introduced = naturalized_available and _flag(row.get("naturalized"))
    if native and introduced:
        return "unresolved", "unresolved", "native_naturalized_conflict"
    if native:
        origin = "native"
    elif introduced:
        origin = "introduced"
    else:
        origin = "unresolved"

    endemic = "unresolved"
    note = ""
    if origin == "native" and _flag(row.get("end_list")) and not _flag(row.get("quest_end_list")):
        if _flag(row.get("endemic_list")):
            endemic = "endemic"
        elif _known_zero(row.get("endemic_list")):
            endemic = "nonendemic"
        else:
            note = "endemism_not_explicit"
    elif origin == "native" and not _flag(row.get("end_list")):
        note = "list_level_endemism_not_indicated"
    return origin, endemic, note


def build_status_ledger(
    checklist_rows: pd.DataFrame,
    list_meta: pd.DataFrame,
    spatial_matches: pd.DataFrame,
    island_species: pd.DataFrame,
    *,
    gift_version: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build source-level status rows only for exact frozen-flora memberships."""
    meta = prepare_list_metadata(list_meta)
    accepted_matches = spatial_matches.loc[
        spatial_matches["spatial_match_status"].eq("accepted")
        & spatial_matches["island_id"].astype(str).ne("")
    ][["entity_ID", "island_id"]].drop_duplicates("entity_ID")
    accepted_matches["entity_ID"] = _text(accepted_matches["entity_ID"])
    if accepted_matches.empty:
        return pd.DataFrame(), pd.DataFrame()

    rows = checklist_rows.copy()
    rows["list_ID"] = _text(rows["list_ID"])
    if "ref_ID" in rows.columns:
        rows["ref_ID"] = _text(rows["ref_ID"])
    else:
        rows["ref_ID"] = ""
    species_column = _column(rows, ("work_species", "species"))
    rows["accepted_species"] = _text(rows[species_column])
    rows = rows.loc[rows["accepted_species"].ne("") & rows["list_ID"].ne("")].copy()

    # List metadata are authoritative for entity_ID and for whether status fields
    # were actually stated by the source/list.
    merge_columns = ["list_ID", "entity_ID", "ref_ID"]
    for column in ("native_indicated", "natural_indicated", "end_ref", "end_list", "geo_entity"):
        if column in meta.columns:
            merge_columns.append(column)
    meta_for_merge = meta[merge_columns].drop_duplicates("list_ID")
    rows = rows.drop(columns=[column for column in merge_columns if column in rows.columns and column != "list_ID"])
    rows = rows.merge(meta_for_merge, on="list_ID", how="inner", validate="many_to_one")
    rows["entity_ID"] = _text(rows["entity_ID"])
    rows = rows.merge(accepted_matches, on="entity_ID", how="inner", validate="many_to_one")

    flora = normalize_island_species(island_species)
    rows = rows.merge(
        flora.assign(_in_frozen_flora=True),
        on=["island_id", "accepted_species"],
        how="left",
        validate="many_to_one",
    )
    audit_rows: list[dict[str, Any]] = []
    ledger_rows: list[dict[str, Any]] = []
    for _, row in rows.iterrows():
        in_flora = bool(row.get("_in_frozen_flora") is True)
        origin, endemic, note = classify_gift_status(row)
        ref_id = str(row.get("ref_ID", "")).strip()
        list_id = str(row.get("list_ID", "")).strip()
        entity_id = str(row.get("entity_ID", "")).strip()
        evidence_key = "|".join(
            [gift_version, ref_id, list_id, entity_id, str(row["accepted_species"])]
        )
        evidence_id = "GIFT_STATUS_" + hashlib.sha256(evidence_key.encode("utf-8")).hexdigest()[:20]
        audit_rows.append(
            {
                "island_id": row["island_id"],
                "accepted_species": row["accepted_species"],
                "entity_ID": entity_id,
                "list_ID": list_id,
                "ref_ID": ref_id,
                "origin_status": origin,
                "endemic_status": endemic,
                "status_note": note,
                "in_frozen_flora": in_flora,
            }
        )
        if not in_flora:
            continue
        ledger_rows.append(
            {
                "island_id": row["island_id"],
                "accepted_species": row["accepted_species"],
                "origin_status": origin,
                "endemic_status": endemic,
                "status_source": f"GIFT_{gift_version}",
                "status_reference": f"ref_ID={ref_id};list_ID={list_id};entity_ID={entity_id}",
                "status_evidence_id": evidence_id,
            }
        )
    return pd.DataFrame(ledger_rows), pd.DataFrame(audit_rows)


@app.command("run")
def run(
    gift_checklists_csv: Path = typer.Option(..., exists=True),
    gift_list_metadata_csv: Path = typer.Option(..., exists=True),
    gift_shapes_gpkg: Path = typer.Option(..., exists=True),
    gshhg_islands_gpkg: Path = typer.Option(..., exists=True),
    island_species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    gift_version: str = typer.Option("latest"),
    min_bidirectional_overlap: float = typer.Option(0.70, min=0.0, max=1.0),
    min_iou: float = typer.Option(0.50, min=0.0, max=1.0),
    ambiguity_margin: float = typer.Option(0.20, min=0.0, max=1.0),
) -> None:
    gift_shapes = gpd.read_file(gift_shapes_gpkg)
    gshhg = gpd.read_file(gshhg_islands_gpkg)
    matches, candidates = match_gift_islands_to_gshhg(
        gift_shapes,
        gshhg,
        min_bidirectional_overlap=min_bidirectional_overlap,
        min_iou=min_iou,
        ambiguity_margin=ambiguity_margin,
    )
    checklist_rows = pd.read_csv(gift_checklists_csv, dtype=str).fillna("")
    list_meta = pd.read_csv(gift_list_metadata_csv, dtype=str).fillna("")
    island_species = pd.read_csv(island_species_csv)
    ledger, audit = build_status_ledger(
        checklist_rows,
        list_meta,
        matches,
        island_species,
        gift_version=gift_version,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    matches.to_csv(output_dir / "gift_to_gshhg_island_matches.csv", index=False)
    candidates.to_csv(output_dir / "gift_to_gshhg_spatial_candidates.csv.gz", index=False)
    ledger.to_csv(output_dir / "gift_floristic_status_ledger.csv.gz", index=False)
    audit.to_csv(output_dir / "gift_floristic_status_audit.csv.gz", index=False)
    manifest = {
        "contract": "gift_floristic_status_match_v1",
        "gift_version": gift_version,
        "n_gift_entities": int(matches["entity_ID"].nunique()) if len(matches) else 0,
        "spatial_match_status": (
            {str(k): int(v) for k, v in matches["spatial_match_status"].value_counts().items()}
            if len(matches)
            else {}
        ),
        "n_status_ledger_rows": int(len(ledger)),
        "n_status_species": int(ledger["accepted_species"].nunique()) if len(ledger) else 0,
        "n_status_islands": int(ledger["island_id"].nunique()) if len(ledger) else 0,
        "status_policy": (
            "GIFT source/list status only; exact Island geometry only; ambiguous spatial/taxonomic "
            "matches unresolved; endemic_list only for focal-island endemism."
        ),
    }
    (output_dir / "gift_floristic_status_match_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
