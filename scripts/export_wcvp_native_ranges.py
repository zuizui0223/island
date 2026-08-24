#!/usr/bin/env python3
"""Export exact-name WCVP native TDWG-L3 ranges from the official bulk ZIP."""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import pandas as pd


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def export_ranges(island_taxa_csv: Path, wcvp_zip: Path, output_dir: Path) -> dict[str, object]:
    target = pd.read_csv(island_taxa_csv, usecols=["accepted_species"])
    target["accepted_species"] = _text(target["accepted_species"])
    target_names = set(target.loc[target["accepted_species"].ne(""), "accepted_species"])

    output_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(wcvp_zip) as archive:
        if "wcvp_names.csv" not in archive.namelist() or "wcvp_distribution.csv" not in archive.namelist():
            raise RuntimeError("official WCVP ZIP lacks expected names/distribution tables")

        matched_chunks: list[pd.DataFrame] = []
        with archive.open("wcvp_names.csv") as handle:
            for chunk in pd.read_csv(
                handle,
                sep="|",
                usecols=["plant_name_id", "taxon_rank", "taxon_status", "taxon_name"],
                dtype={"plant_name_id": "Int64", "taxon_rank": str, "taxon_status": str, "taxon_name": str},
                chunksize=250_000,
                low_memory=False,
            ):
                chunk["taxon_name"] = _text(chunk["taxon_name"])
                keep = (
                    _text(chunk["taxon_rank"]).str.casefold().eq("species")
                    & _text(chunk["taxon_status"]).str.casefold().eq("accepted")
                    & chunk["taxon_name"].isin(target_names)
                )
                selected = chunk.loc[keep, ["plant_name_id", "taxon_name"]].copy()
                if not selected.empty:
                    selected = selected.rename(columns={"taxon_name": "accepted_species"})
                    matched_chunks.append(selected)
        matched = (
            pd.concat(matched_chunks, ignore_index=True).dropna(subset=["plant_name_id"]).drop_duplicates()
            if matched_chunks
            else pd.DataFrame(columns=["plant_name_id", "accepted_species"])
        )
        if matched.empty:
            raise RuntimeError("no fixed-master species exactly matched accepted WCVP species")
        matched["plant_name_id"] = matched["plant_name_id"].astype("int64")
        ids = set(matched["plant_name_id"].tolist())
        id_to_species = matched.set_index("plant_name_id")["accepted_species"].to_dict()

        distribution_chunks: list[pd.DataFrame] = []
        with archive.open("wcvp_distribution.csv") as handle:
            for chunk in pd.read_csv(
                handle,
                sep="|",
                usecols=[
                    "plant_name_id",
                    "area_code_l3",
                    "introduced",
                    "extinct",
                    "location_doubtful",
                ],
                dtype={"plant_name_id": "Int64", "area_code_l3": str},
                chunksize=500_000,
                low_memory=False,
            ):
                chunk = chunk.loc[chunk["plant_name_id"].isin(ids)].copy()
                if chunk.empty:
                    continue
                for column in ("introduced", "extinct", "location_doubtful"):
                    chunk[column] = pd.to_numeric(chunk[column], errors="coerce")
                chunk["area_code_l3"] = _text(chunk["area_code_l3"])
                chunk = chunk.loc[
                    chunk["introduced"].eq(0)
                    & chunk["extinct"].eq(0)
                    & chunk["location_doubtful"].eq(0)
                    & chunk["area_code_l3"].ne("")
                ].copy()
                if not chunk.empty:
                    chunk["accepted_species"] = chunk["plant_name_id"].map(id_to_species)
                    distribution_chunks.append(
                        chunk[["accepted_species", "plant_name_id", "area_code_l3"]]
                    )

    native = (
        pd.concat(distribution_chunks, ignore_index=True).drop_duplicates()
        if distribution_chunks
        else pd.DataFrame(columns=["accepted_species", "plant_name_id", "area_code_l3"])
    )
    if native.empty:
        summary = matched.copy()
        summary["n_native_l3"] = 0
        summary["native_l3_codes"] = ""
    else:
        grouped = (
            native.groupby(["accepted_species", "plant_name_id"], as_index=False)
            .agg(
                n_native_l3=("area_code_l3", "nunique"),
                native_l3_codes=("area_code_l3", lambda x: "|".join(sorted(set(x)))),
            )
        )
        summary = matched.merge(
            grouped,
            on=["accepted_species", "plant_name_id"],
            how="left",
            validate="one_to_one",
        )
        summary["n_native_l3"] = summary["n_native_l3"].fillna(0).astype(int)
        summary["native_l3_codes"] = summary["native_l3_codes"].fillna("")

    summary.sort_values("accepted_species").to_csv(
        output_dir / "wcvp_native_range_summary.csv.gz", index=False
    )
    manifest = {
        "contract": "official_wcvp_bulk_exact_name_native_range_v1",
        "n_target_species": len(target_names),
        "n_exact_accepted_species": int(summary["accepted_species"].nunique()),
        "n_species_with_native_range": int(summary["n_native_l3"].gt(0).sum()),
        "n_single_native_l3_species": int(summary["n_native_l3"].eq(1).sum()),
        "source_files": ["wcvp_names.csv", "wcvp_distribution.csv"],
        "distribution_filters": {
            "introduced": 0,
            "extinct": 0,
            "location_doubtful": 0,
        },
        "policy": (
            "Exact accepted species only. A single native TDWG-L3 range is corroboration "
            "only; final endemic status also requires an independent GIFT endemic claim "
            "and focal-island membership in the same TDWG region."
        ),
    }
    (output_dir / "wcvp_native_range_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("island_taxa_csv", type=Path)
    parser.add_argument("wcvp_zip", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(export_ranges(args.island_taxa_csv, args.wcvp_zip, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
