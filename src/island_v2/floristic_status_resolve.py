"""Resolve native/introduced status independently from conservative endemism.

GIFT source/list endemism is retained only as a source claim. GIFT explicitly
warns that a taxon stated as endemic in one checklist may occur in another
checklist, so it is not sufficient for focal-island endemism.

Final endemic status is therefore corroborated only when:
1. GIFT reports the species as endemic in a source list for the focal island;
2. WCVP has an exact accepted-species match and exactly one non-doubtful,
   non-extinct native TDWG level-3 region; and
3. the focal frozen island maps unambiguously to that same TDWG level-3 region.

A WCVP range spanning >1 native level-3 region is sufficient to classify the
species as non-endemic at this conservative regional scale. Everything else
remains unresolved.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

ORIGIN_VALUES = {"native", "introduced", "unresolved"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def collapse_gift_origin(raw: pd.DataFrame) -> pd.DataFrame:
    """Collapse origin independently; endemism disagreement never erases origin."""
    required = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required - set(raw.columns)
    if missing:
        raise typer.BadParameter(f"GIFT status ledger missing columns: {sorted(missing)}")

    rows: list[dict[str, Any]] = []
    for (island_id, species), group in raw.groupby(["island_id", "accepted_species"], sort=True):
        origins = {
            _text(value).lower()
            for value in group["origin_status"]
            if _text(value).lower() in {"native", "introduced"}
        }
        origin = next(iter(origins)) if len(origins) == 1 else "unresolved"
        gift_endemic_claim = group["endemic_status"].astype(str).str.lower().eq("endemic").any()
        gift_nonendemic_claim = group["endemic_status"].astype(str).str.lower().eq("nonendemic").any()
        refs = "|".join(sorted({_text(v) for v in group.get("status_reference", pd.Series(dtype=str)) if _text(v)}))
        rows.append(
            {
                "island_id": _text(island_id),
                "accepted_species": _text(species),
                "origin_status": origin,
                "origin_conflict": len(origins) > 1,
                "gift_endemic_claim": bool(gift_endemic_claim),
                "gift_nonendemic_claim": bool(gift_nonendemic_claim),
                "gift_status_references": refs,
            }
        )
    return pd.DataFrame(rows)


def map_islands_to_tdwg_l3(
    islands: gpd.GeoDataFrame, tdwg: gpd.GeoDataFrame
) -> pd.DataFrame:
    """Map an island to one TDWG-L3 region only when its representative point is unambiguous."""
    if "island_id" not in islands.columns:
        raise typer.BadParameter("island geometry missing island_id")
    if "LEVEL3_COD" not in tdwg.columns:
        raise typer.BadParameter("TDWG geometry missing LEVEL3_COD")
    if islands.crs is None or tdwg.crs is None:
        raise typer.BadParameter("island and TDWG geometries require CRS metadata")

    isl = islands[["island_id", "geometry"]].copy().to_crs(4326)
    geo = tdwg[["LEVEL3_COD", "LEVEL3_NAM", "geometry"]].copy().to_crs(4326)
    isl["geometry"] = isl.geometry.representative_point()
    joined = gpd.sjoin(isl, geo, how="left", predicate="within")

    rows: list[dict[str, Any]] = []
    for island_id, group in joined.groupby("island_id", sort=True):
        codes = sorted({_text(v) for v in group["LEVEL3_COD"] if _text(v)})
        names = sorted({_text(v) for v in group["LEVEL3_NAM"] if _text(v)})
        rows.append(
            {
                "island_id": _text(island_id),
                "tdwg_l3_code": codes[0] if len(codes) == 1 else "",
                "tdwg_l3_name": names[0] if len(names) == 1 else "",
                "tdwg_match_status": "accepted" if len(codes) == 1 else (
                    "no_match" if len(codes) == 0 else "ambiguous"
                ),
            }
        )
    return pd.DataFrame(rows)


def resolve_endemism(
    gift_origin: pd.DataFrame,
    wcvp: pd.DataFrame,
    island_tdwg: pd.DataFrame,
) -> pd.DataFrame:
    required_wcvp = {"accepted_species", "n_native_l3", "native_l3_codes"}
    missing = required_wcvp - set(wcvp.columns)
    if missing:
        raise typer.BadParameter(f"WCVP summary missing columns: {sorted(missing)}")

    work = gift_origin.merge(
        wcvp[list(required_wcvp)].drop_duplicates("accepted_species"),
        on="accepted_species",
        how="left",
        validate="many_to_one",
    ).merge(
        island_tdwg.drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    work["n_native_l3"] = pd.to_numeric(work["n_native_l3"], errors="coerce")

    endemic_status: list[str] = []
    basis: list[str] = []
    scope: list[str] = []
    for row in work.itertuples(index=False):
        if row.origin_status != "native":
            endemic_status.append("unresolved")
            basis.append("origin_not_confirmed_native")
            scope.append("")
            continue
        n_native = row.n_native_l3
        codes = [v for v in _text(row.native_l3_codes).split("|") if v]
        island_code = _text(getattr(row, "tdwg_l3_code", ""))
        match_status = _text(getattr(row, "tdwg_match_status", ""))
        if pd.notna(n_native) and float(n_native) > 1:
            endemic_status.append("nonendemic")
            basis.append("wcvp_multiple_native_tdwg_l3_regions")
            scope.append("tdwg_l3")
        elif (
            pd.notna(n_native)
            and int(n_native) == 1
            and len(codes) == 1
            and match_status == "accepted"
            and island_code == codes[0]
            and bool(row.gift_endemic_claim)
        ):
            endemic_status.append("endemic")
            basis.append("gift_endemic_claim_plus_wcvp_single_matching_native_tdwg_l3")
            scope.append("tdwg_l3_corroborated")
        else:
            endemic_status.append("unresolved")
            basis.append("insufficient_global_endemism_corroboration")
            scope.append("")

    work["endemic_status"] = endemic_status
    work["endemism_basis"] = basis
    work["endemism_scope"] = scope
    work["status_source"] = "GIFT_origin+WCVP_endemism"
    return work


@app.command("run")
def run(
    gift_status_ledger_csv: Path = typer.Option(..., exists=True),
    wcvp_native_range_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    tdwg_gpkg: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    raw = pd.read_csv(gift_status_ledger_csv)
    wcvp = pd.read_csv(wcvp_native_range_csv)
    islands = gpd.read_file(islands_gpkg)
    tdwg = gpd.read_file(tdwg_gpkg)

    origin = collapse_gift_origin(raw)
    mapping = map_islands_to_tdwg_l3(islands, tdwg)
    resolved = resolve_endemism(origin, wcvp, mapping)

    output_dir.mkdir(parents=True, exist_ok=True)
    origin.to_csv(output_dir / "gift_origin_collapsed.csv.gz", index=False)
    mapping.to_csv(output_dir / "island_tdwg_l3_mapping.csv", index=False)
    resolved.to_csv(output_dir / "resolved_floristic_status_ledger.csv.gz", index=False)

    manifest = {
        "contract": "floristic_status_resolve_v1",
        "n_rows": int(len(resolved)),
        "n_islands": int(resolved["island_id"].nunique()),
        "n_species": int(resolved["accepted_species"].nunique()),
        "origin_status_counts": {
            str(k): int(v) for k, v in resolved["origin_status"].value_counts(dropna=False).items()
        },
        "endemic_status_counts": {
            str(k): int(v) for k, v in resolved["endemic_status"].value_counts(dropna=False).items()
        },
        "policy": (
            "GIFT list endemism is never sufficient alone; endemic requires a GIFT source claim "
            "plus one matching WCVP native TDWG-L3 region. WCVP >1 native regions => nonendemic."
        ),
    }
    (output_dir / "floristic_status_resolve_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
