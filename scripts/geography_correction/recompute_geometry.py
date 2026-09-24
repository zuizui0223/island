"""Recompute locked-island distances from the exact parent GSHHG archive."""

import argparse
import csv
import hashlib
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely import make_valid
from shapely.geometry import Point
from shapely.ops import unary_union

from island_v2.spherical_coast_distance import CoastIndex, geometry_segments

SHA = "8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
BAD = "gshhg_2.3.7_h_8b13189234b949ee1ff6"
SEEDS = {
    "africa": (20, 0),
    "eurasia": (80, 50),
    "north_america": (-100, 40),
    "south_america": (-60, -15),
    "australia": (135, -25),
    "antarctica": (0, -80),
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--islands", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--sites", type=Path, help="Optional locked site table with site_key, latitude, longitude"
    )
    args = p.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest() != SHA:
        raise ValueError("GSHHG archive SHA256 mismatch")
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Geometry output must be new or empty")
    args.output.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.archive) as archive:
        for name in archive.namelist():
            if Path(name).name.startswith(("GSHHS_h_L1.", "GSHHS_h_L5.")):
                (args.output / Path(name).name).write_bytes(archive.read(name))
    land = gpd.read_file(args.output / "GSHHS_h_L1.shp")
    ice = gpd.read_file(args.output / "GSHHS_h_L5.shp")
    if (
        hashlib.sha256(args.islands.read_bytes()).hexdigest()
        != "0d57c81501bc50e36bec0f3875b47d934899543ec2687289bb3a7a99f20dc5b8"
    ):
        raise ValueError("Locked island GPKG SHA256 mismatch")
    islands = gpd.read_file(args.islands)
    if len(islands) != 8265 or not islands.island_id.is_unique:
        raise ValueError("Expected unique locked 8265-unit universe")
    raw = land[land.id.eq("0-W")]
    bad = islands[islands.island_id.eq(BAD)]
    if len(raw) != 1 or len(bad) != 1 or not bad.geometry.iloc[0].equals(raw.geometry.iloc[0]):
        raise ValueError("Excluded component does not match exact source geometry")
    rows = []
    for name, (lon, lat) in SEEDS.items():
        layer = ice if name == "antarctica" else land
        hits = layer[layer.geometry.covers(Point(lon, lat))]
        if len(hits) == 0 or hits.sibling_id.nunique() != 1:
            raise ValueError(f"Ambiguous continental seed: {name}")
        sibling = int(hits.sibling_id.iloc[0])
        rows.append(
            {
                "continent": name,
                "sibling_id": sibling,
                "geometry": make_valid(unary_union(layer[layer.sibling_id.eq(sibling)].geometry)),
            }
        )
    coast = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
    coast.to_file(args.output / "gshhg_continents.gpkg", driver="GPKG")
    parts = [geometry_segments(g) for g in coast.geometry]
    index = CoastIndex(np.concatenate([a for a, b in parts]), np.concatenate([b for a, b in parts]))
    with (args.output / "gshhg_spherical_distances_all.csv").open(
        "w", newline="", encoding="utf8"
    ) as f:
        writer = csv.writer(f)
        writer.writerow(["island_id", "spherical_coast_distance_km"])
        for n, row in enumerate(islands[~islands.island_id.eq(BAD)].itertuples(), 1):
            distance = index.distance(*geometry_segments(row.geometry))
            if not np.isfinite(distance) or distance <= 0:
                raise ValueError(f"Unresolved island distance: {row.island_id}")
            writer.writerow([row.island_id, distance])
            if n % 100 == 0:
                print("Islands:", n, flush=True)
    if args.sites:
        import pandas as pd

        sites = pd.read_csv(args.sites)
        if sites.site_key.duplicated().any():
            raise ValueError("Duplicate site identity")
        continental_land = unary_union(coast.geometry)
        sites["spherical_distance_km"] = [
            index.point_distance(r.longitude, r.latitude, continental_land)
            for r in sites.itertuples()
        ]
        sites.to_csv(args.output / "glopl_corrected_site_distances.csv", index=False)


if __name__ == "__main__":
    main()
