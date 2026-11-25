"""End-to-end run of both range-map scripts on synthetic GeoPackages (no real data)."""
import json
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

pyogrio = pytest.importorskip("pyogrio")
shapely = pytest.importorskip("shapely")
pytest.importorskip("pyproj")
from shapely.geometry import box  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]


def alpha(i):
    return "abcdefghijklmnopqrstuvwxyz"[i // 26] + "abcdefghijklmnopqrstuvwxyz"[i % 26]


def write_layer(path, geoms, fields):
    pyogrio.raw.write(
        str(path), shapely.to_wkb(np.asarray(geoms, dtype=object)), geometry_type="Polygon",
        field_data=[np.asarray(v) for v in fields.values()], fields=list(fields), crs="EPSG:4326", driver="GPKG",
    )


@pytest.fixture
def world(tmp_path):
    rng = np.random.default_rng(7)
    n_groups, per = 34, 7
    island_rows, geoms = [], []
    for g in range(n_groups):
        for k in range(per):
            lon, lat = -170 + g * 9.5 + k * 0.6, -20 + (g % 8) * 5 + k * 0.2
            geoms.append(box(lon, lat, lon + 0.1, lat + 0.1))  # ~110 km2
            island_rows.append(dict(island=f"I{g}_{k}", group=f"A{g}", region="tropics",
                                    dist_km=float(rng.uniform(50, 3000)), lat=lat + 0.05, lon=lon + 0.05))
    islands = pd.DataFrame(island_rows)
    write_layer(tmp_path / "islands.gpkg", geoms, {"gshhg_id": islands["island"].to_numpy()})

    # birds: 3 specialist genera, partial and other species, each ranging over random islands
    names, diet, fams, rgeoms, codes = [], [], [], [], []
    for j in range(60):
        kind = "spec" if j < 15 else "part" if j < 30 else "oth"
        genus = "Gen" + alpha(j // 3)
        name = f"{genus} sp{alpha(j)}"
        names.append(name)
        diet.append({"spec": 70, "part": 20, "oth": 0}[kind])
        fams.append("Trochilidae" if kind == "spec" and j < 6 else "Meliphagidae" if kind == "spec" else "Turdidae")
        for idx in rng.choice(len(geoms), size=int(rng.integers(20, 90)), replace=False):
            x0, y0, x1, y1 = geoms[idx].bounds
            rgeoms.append(box(x0 - 0.01, y0 - 0.01, x1 + 0.01, y1 + 0.01))
            codes.append((name, 1, 1, 1))
    # a non-native polygon that must be ignored, and a non-breeding one for the sensitivity only
    rgeoms += [box(-180, -30, 180, 30), box(-180, -30, 180, 30)]
    codes += [("Genaa spaa", 1, 3, 1), ("Genzz spzz", 1, 1, 3)]
    names.append("Genzz spzz"); diet.append(0); fams.append("Turdidae")
    c = pd.DataFrame(codes, columns=["sci_name", "presence", "origin", "seasonal"])
    write_layer(tmp_path / "ranges.gpkg", rgeoms, {k: c[k].to_numpy() for k in c})
    pd.DataFrame({"Scientific": names, "Diet-Nect": diet, "BLFamilyLatin": fams}).to_csv(
        tmp_path / "elton.txt", sep="\t", index=False
    )
    islands.to_csv(tmp_path / "islands.csv", index=False)
    return tmp_path, islands


def run(script, *args):
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-3000:]
    return r.stdout


def test_incidence_then_flora_run(world):
    tmp, islands = world
    out = tmp / "inc"
    run("range_map_nectarivore_incidence.py", "--ranges", tmp / "ranges.gpkg", "--islands", tmp / "islands.gpkg",
        "--island-id-col", "gshhg_id", "--elton", tmp / "elton.txt", "--out", out, "--chunk", "500")
    rep = json.loads((out / "incidence_report.json").read_text())
    assert rep["islands"] == len(islands) and rep["islands_ge_1km2"] == len(islands)
    prim, nonb = rep["variants"]["primary"], rep["variants"]["nonbreeding"]
    assert prim["unmatched_occurring"] == 0
    assert prim["diet_class_counts"] == {"other": 30, "partial": 15, "specialist": 15}
    assert nonb["species_kept"] == prim["species_kept"] + 1  # non-breeding range only counts in sensitivity
    sp = pd.read_csv(out / "species_primary.csv")
    assert set(sp.loc[sp.diet_class == 2, "clade"]) == {0, 1}  # trochilidae, meliphagidae

    # plants: tubular share rises with specialist count through genus composition only (co-occurrence)
    from scipy import sparse

    cover = sparse.load_npz(out / "cover_primary.npz").tocsr()
    inc_isl = pd.read_csv(out / "islands.csv")
    spec = np.asarray((cover[:, (sp.diet_class == 2).to_numpy()] >= 0.5).sum(axis=1)).ravel()
    spec_of = dict(zip(inc_isl["island"], spec))
    rng = np.random.default_rng(11)
    rows = []
    for isl in islands["island"]:
        p_tub_genus = 0.1 + 0.08 * spec_of[isl]
        for k in range(int(rng.integers(8, 25))):
            tub_genus = rng.random() < p_tub_genus
            g = f"T{rng.integers(0, 4)}" if tub_genus else f"O{rng.integers(0, 20)}"
            s = f"{g}_{rng.integers(0, 6)}"
            rows.append(dict(island=isl, species=s, genus=g, grade="high", status_scope="native",
                             tubular=int(tub_genus), deep=int(tub_genus), extended_access=int(tub_genus),
                             red_extended=0))
    plants = pd.DataFrame(rows).drop_duplicates(["island", "species"])
    plants.to_csv(tmp / "plants.csv", index=False)
    summary = tmp / "summary.json"
    stdout = run("range_map_nectarivore_flora_run.py", "--incidence", out, "--islands", tmp / "islands.csv",
                 "--plants", tmp / "plants.csv", "--n-boot", "99", "--n-perm", "49", "--out", summary)
    s = json.loads(summary.read_text())
    prim = s["primary"]
    assert prim["status"] == "ok" and prim["n_informative_groups"] >= 30
    assert prim["decision"] == "lineage_cooccurrence", stdout
    assert set(prim["clade_replication"]) == {"trochilidae", "meliphagidae", "nectariniidae", "drepanidini", "other"}
    red = [x for x in s["secondary"] if x["endpoint"] == "red_extended"][0]
    assert red["decision"] == "not_estimable"  # state absent everywhere: not a null result
    assert prim["bird_genera_with_mixed_diet_present"] == 0
    assert prim["permutation"]["share"]["genus"]["p_one_sided"] == 1.0  # no within-genus diet variation
    assert len(s["secondary"]) == 5 and all("holm_p_residual_genus" in x for x in s["secondary"])
    assert len(s["sensitivity"]) == 6
    text = summary.read_text()
    assert "I0_0" not in text and "Genaa" not in text  # no identifiers leave the protected run
