#!/usr/bin/env python3
"""Island x bird incidence from BirdLife BOTW range maps (contract: docs/RANGE_MAP_NECTARIVORE_CONTRACT_20261125.md).

Needs the optional geo extras: ``pip install -e .[geo]`` (shapely, pyproj, pyogrio).

    python scripts/range_map_nectarivore_incidence.py \
        --ranges BOTW_2024.gdb --layer All_Species \
        --islands gshhg_gift_islands.gpkg --island-id-col gshhg_id \
        --elton BirdFuncDat.txt --synonyms botw_to_elton_synonyms.csv \
        --out work/range_map_incidence

BirdLife maps are licensed for research use only. Everything written to --out is derived from them:
keep it in protected storage, never commit it, and export only the anonymous aggregates produced by
scripts/range_map_nectarivore_flora_run.py.
"""
import argparse
import hashlib
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import sparse

from fd.range_map_nectarivore import (
    CLADES,
    SEASONAL_PRIMARY,
    SEASONAL_WITH_NONBREEDING,
    build_species_table,
    crosswalk_names,
    island_incidence_streaming,
    nectarivore_clade,
    normalise_binomial,
)

EQUAL_AREA = "EPSG:6933"
VARIANTS = {
    # name: (seasonal codes, island buffer in metres)
    "primary": (SEASONAL_PRIMARY, 0.0),
    "nonbreeding": (SEASONAL_WITH_NONBREEDING, 0.0),
    "buffer1km": (SEASONAL_PRIMARY, 1000.0),
}
COVER_THRESHOLDS = (0.5, 0.1, 0.0)


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    if path.is_dir():  # file geodatabase
        for f in sorted(p for p in path.rglob("*") if p.is_file()):
            h.update(f.name.encode())
            with f.open("rb") as fh:
                for block in iter(lambda: fh.read(1 << 20), b""):
                    h.update(block)
    else:
        with path.open("rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
    return h.hexdigest()


def projector(src_crs):
    import shapely
    from pyproj import Transformer

    tr = Transformer.from_crs(src_crs, EQUAL_AREA, always_xy=True)

    def project(geoms):
        return shapely.transform(geoms, lambda xy: np.column_stack(tr.transform(xy[:, 0], xy[:, 1])))

    return project


def read_islands(path, layer, id_col):
    import pyogrio
    import shapely

    meta, _, wkb, fields = pyogrio.raw.read(path, layer=layer, columns=[id_col])
    ids = fields[0]
    geoms = projector(meta["crs"])(shapely.from_wkb(wkb))
    # islands split at the antimeridian come as several rows with one id: dissolve them
    frame = pd.DataFrame({"id": ids, "g": list(geoms)})
    merged = frame.groupby("id")["g"].agg(lambda s: shapely.union_all(list(s)))
    return tuple(merged.index), np.asarray(merged.to_list(), dtype=object)


def range_chunks(path, layer, seasonal, chunk, name_col="sci_name"):
    import pyogrio
    import shapely

    info = pyogrio.read_info(path, layer=layer)
    cols = {c.lower(): c for c in info["fields"]}
    for need in (name_col, "presence", "origin", "seasonal"):
        if need not in cols:
            raise KeyError(f"range layer lacks {need}; fields are {list(info['fields'])}")
    where = (
        f"{cols['presence']} IN (1,2) AND {cols['origin']} IN (1,2) "
        f"AND {cols['seasonal']} IN ({','.join(str(s) for s in sorted(seasonal))})"
    )
    project = projector(info["crs"])
    start = 0
    while True:
        meta, _, wkb, fields = pyogrio.raw.read(
            path, layer=layer, columns=[cols[name_col]], where=where, skip_features=start, max_features=chunk
        )
        if len(wkb) == 0:
            return
        geoms = shapely.make_valid(project(shapely.from_wkb(wkb)))
        yield list(fields[0]), geoms
        start += len(wkb)


def read_elton(path):
    t = pd.read_csv(path, sep="\t", encoding="latin-1", low_memory=False)
    need = {"Scientific", "Diet-Nect", "BLFamilyLatin"}
    if need - set(t.columns):
        raise KeyError(f"not an EltonTraits bird table: missing {sorted(need - set(t.columns))}")
    t = t[t["Scientific"].notna() & t["Diet-Nect"].notna()]
    return t.rename(columns={"Scientific": "scientific", "Diet-Nect": "diet_nect", "BLFamilyLatin": "family"})[
        ["scientific", "diet_nect", "family"]
    ]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ranges", type=pathlib.Path, required=True)
    p.add_argument("--layer", default=None)
    p.add_argument("--islands", type=pathlib.Path, required=True)
    p.add_argument("--islands-layer", default=None)
    p.add_argument("--island-id-col", required=True)
    p.add_argument("--elton", type=pathlib.Path, required=True)
    p.add_argument("--synonyms", type=pathlib.Path, default=None)
    p.add_argument("--chunk", type=int, default=2000)
    p.add_argument("--out", type=pathlib.Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    island_ids, island_geoms = read_islands(a.islands, a.islands_layer, a.island_id_col)
    import shapely

    area_km2 = shapely.area(island_geoms) / 1e6
    pd.DataFrame({"island": island_ids, "area_km2": area_km2}).to_csv(a.out / "islands.csv", index=False)
    traits = read_elton(a.elton)
    synonyms = pd.read_csv(a.synonyms) if a.synonyms else None
    report = {
        "inputs": {
            "ranges_sha256": sha256(a.ranges),
            "islands_sha256": sha256(a.islands),
            "elton_sha256": sha256(a.elton),
            "synonyms_sha256": sha256(a.synonyms) if a.synonyms else None,
        },
        "islands": len(island_ids),
        "islands_ge_1km2": int((area_km2 >= 1).sum()),
        "variants": {},
    }
    qualifying = area_km2 >= 1
    for name, (seasonal, buf) in VARIANTS.items():
        inc = island_incidence_streaming(
            island_ids, island_geoms, range_chunks(a.ranges, a.layer, seasonal, a.chunk), min_cover=0.0, island_buffer=buf
        )
        cw = crosswalk_names(inc.species, traits["scientific"], synonyms)
        sub, table, rep = build_species_table(inc, traits, cw, island_mask=qualifying)
        sparse.save_npz(a.out / f"cover_{name}.npz", sub.cover.tocsr())
        fam_name = traits.assign(key=traits["scientific"].map(normalise_binomial)).set_index("key")["family"]
        species = pd.DataFrame(
            {
                "species": table.species,
                "diet_class": table.diet,
                "genus_block": table.genus,
                "family_block": table.family,
                "clade": [
                    nectarivore_clade(str(fam_name.get(cw.mapping[s], "")), s.split(" ")[0]) for s in table.species
                ],
            }
        )
        species.to_csv(a.out / f"species_{name}.csv", index=False)
        rep["unmatched_species"] = list(cw.unmatched)
        rep["species_by_clade"] = {
            c: int(((species["clade"] == k) & (species["diet_class"] == 2)).sum()) for k, c in enumerate(CLADES)
        }
        rep["incidence_pairs_by_cover"] = {str(t): int((sub.cover.data >= t if t > 0 else sub.cover.data > 0).sum()) for t in COVER_THRESHOLDS}
        report["variants"][name] = rep
        print(name, json.dumps({k: v for k, v in rep.items() if k != "unmatched_species"}), flush=True)
    (a.out / "incidence_report.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
