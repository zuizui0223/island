"""Geological-origin audit and LGM-connectivity crosswalk for fd.

Uses independent, published Weigelt 2013 island-centroid/area/GMMC data
and the Zizka 2022 supplementary *missing-geology additions only*.
Neither great-circle distance nor GMMC=0 is geological proof of oceanic origin.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

import numpy as np

RADIUS_KM = 6371.0088
WEIGELT = "https://zenodo.org/records/4971536/files/Weigelt_etal_2013_PNAS_islanddata.csv?download=1"
ZIZKA = "https://zenodo.org/records/6325640/files/input.zip?download=1"
ROEBLE = ("https://media.springernature.com/original/springer-static/"
          "esm/art%3A10.1038%2Fs41467-024-51556-7/MediaObjects/"
          "41467_2024_51556_MOESM6_ESM.xlsx")
GSHHG_COMMIT = "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c"
GSHHG_URL = ("https://raw.githubusercontent.com/zuizui0223/island/" + GSHHG_COMMIT
             + "/results/geography_20260924/corrected_geography_covariates.csv")


def download(url: str, dest: Path, max_bytes: int) -> str:
    sha = hashlib.sha256()
    size = 0
    request = urllib.request.Request(url, headers={"User-Agent": "fd-geology-audit/1.0"})
    with urllib.request.urlopen(request, timeout=120) as source, dest.open("wb") as out:
        while chunk := source.read(1 << 20):
            size += len(chunk)
            if size > max_bytes:
                raise ValueError("Source exceeds declared maximum byte count")
            sha.update(chunk)
            out.write(chunk)
    if size == 0:
        raise ValueError("Downloaded empty source")
    return sha.hexdigest()


def _number(value) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def _csv(path: Path, *, encoding: str = "utf-8-sig"):
    with path.open(encoding=encoding, newline="") as handle:
        yield from csv.DictReader(handle)


def xyz(lat, lon) -> np.ndarray:
    a, b = np.radians(np.asarray(lat, dtype=float)), np.radians(np.asarray(lon, dtype=float))
    return np.column_stack((np.cos(a) * np.cos(b), np.cos(a) * np.sin(b), np.sin(a)))


def _valid_features(lat: float, lon: float, area: float) -> bool:
    return bool(np.isfinite([lat, lon, area]).all() and -90 <= lat <= 90
                and -180 <= lon <= 180 and area > 0)


def parse_gshhg(path: Path) -> list[dict]:
    rows = []
    for r in _csv(path):
        lat, lon, area = (_number(r.get(k)) for k in
                          ("island_latitude", "island_longitude", "area_km2"))
        if _valid_features(lat, lon, area):
            rows.append({"id": r["island_id"], "lat": lat, "lon": lon, "area": area,
                         "regime": r.get("analysis_regime", ""),
                         "continent_distance_km": _number(r.get("distance_to_continent_km")),
                         "spatial_block": r.get("spatial_block", ""),
                         **{f"climate_pc{j}": _number(r.get(f"climate_pc{j}"))
                            for j in range(1, 5)}})
    if not rows or len({r["id"] for r in rows}) != len(rows):
        raise ValueError("No valid unique geographic islands")
    return rows


def parse_weigelt(path: Path) -> list[dict]:
    rows = []
    for r in _csv(path, encoding="latin1"):
        lat, lon, area = (_number(r.get(k)) for k in ("Lat", "Long", "Area"))
        if _valid_features(lat, lon, area) and r.get("GMMC", "") in ("0", "1"):
            rows.append({"id": str(r["ID"]), "lat": lat, "lon": lon, "area": area,
                         "connection": "lgm_connected" if r["GMMC"] == "1"
                         else "lgm_unconnected", "name": r.get("Island", ""),
                         "archipelago": r.get("Archip", "")})
    return rows


def parse_dbf(buf: bytes) -> list[dict]:
    if len(buf) < 32:
        raise ValueError("Invalid DBF")
    n = int.from_bytes(buf[4:8], "little")
    header = int.from_bytes(buf[8:10], "little")
    record = int.from_bytes(buf[10:12], "little")
    fields = []
    offset = 1
    for position in range(32, header - 1, 32):
        desc = buf[position:position + 32]
        if not desc or desc[0] == 13:
            break
        name = desc[:11].split(bytes([0]), 1)[0].decode("ascii", errors="replace")
        fields.append((name, offset, desc[16]))
        offset += desc[16]
    if offset > record:
        raise ValueError("DBF field offsets exceed record size")
    result = []
    for i in range(n):
        row_bytes = buf[header + i * record:header + (i + 1) * record]
        if len(row_bytes) != record:
            raise ValueError("Truncated DBF rows")
        if row_bytes[0:1] == b"*":
            continue
        row = {name: row_bytes[start:start + width].decode("latin1").strip()
               for name, start, width in fields}
        result.append(row)
    return result


def parse_verified_geology(archive: Path) -> list[dict]:
    with zipfile.ZipFile(archive) as z:
        with z.open("input/entities_with_missing_geology_added.csv") as handle:
            manual = list(csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8-sig")))
        entities = {r["entt_ID"]: r for r in parse_dbf(z.read("input/geoentities_simple.dbf"))}
    rows = []
    for r in manual:
        entity = entities.get(r["entity_ID"])
        if not entity:
            continue
        label = r.get("geology_simple_add", "").strip().lower()
        # Preserve the geological distinction, with mixed and unknown not
        # automatically forced into a binary oceanic/continental class.
        if label in ("oceanic", "volcanic", "atoll", "raised ocean floor"):
            kind = "oceanic"
        elif label in ("continental", "continental shelf", "shelf", "continental fragment"):
            kind = "continental"
        elif label == "mixed":
            kind = "mixed"
        else:
            kind = "unresolved"
        lat, lon, area = (_number(entity.get(k)) for k in ("point_y", "point_x", "area"))
        if _valid_features(lat, lon, area):
            rows.append({"id": r["entity_ID"], "lat": lat, "lon": lon, "area": area,
                         "geological_origin": kind,
                         "geological_source_value": label,
                         "geological_source_name": r.get("geo_entity", ""),
                         "geological_reference": r.get("reference", "")})
    return rows



def parse_roeble_xlsx(path: Path) -> list[dict]:
    """Use Roeble et al 2024 Supplementary Data 3 geology for 510 islands.

    Standard-library XLSX parser reads values from the first worksheet;
    the publicly released table is one worksheet with inline string cells.
    It does not redistribute the supplementary workbook.
    """
    namespace = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in si.findall(".//s:t", namespace))
                      for si in root.findall("s:si", namespace)]
        root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in root.findall(".//s:sheetData/s:row", namespace):
        values = {}
        for cell in row.findall("s:c", namespace):
            address = cell.attrib.get("r", "")
            col = 0
            for char in address:
                if not char.isalpha():
                    break
                col = 26 * col + ord(char.upper()) - 64
            value_node = cell.find("s:v", namespace)
            inline = cell.find("s:is", namespace)
            value = value_node.text if value_node is not None else None
            if cell.attrib.get("t") == "s" and value is not None:
                value = shared[int(value)]
            elif inline is not None:
                value = "".join(t.text or "" for t in inline.findall(".//s:t", namespace))
            values[col-1] = value
        rows.append(values)
    if not rows:
        raise ValueError("Empty Roeble supplementary spreadsheet")
    header = rows[0]
    lookup = {str(value): column for column, value in header.items() if value}
    required = ("entity_ID", "longitude", "latitude", "area", "island_type", "geology")
    if any(k not in lookup for k in required):
        raise ValueError(f"Roeble geological table missing columns: {required}")
    result = []
    for cells in rows[1:]:
        get = lambda field: cells.get(lookup.get(field, -1))
        label = (get("island_type") or "").strip().lower()
        if label not in ("continental", "oceanic"):
            continue
        area = _number(get("area"))
        lat = _number(get("latitude"))
        lon = _number(get("longitude"))
        if not _valid_features(lat, lon, area):
            continue
        result.append({
            "id": str(get("entity_ID")), "lat": lat, "lon": lon, "area": area,
            "geological_origin": label,
            "geological_source_value": str(get("geology") or ""),
            "geological_source_name": str(get("geo_entity") or ""),
            "geological_reference": str(get("geology_ref") or ""),
        })
    if len(result) < 300 or len({r["id"] for r in result}) != len(result):
        raise ValueError(f"Unexpected verified Roeble geology support: {len(result)}")
    return result


def crosswalk(target: list[dict], source: list[dict], *, k: int = 24) -> tuple[dict, dict]:
    """Conservative many-to-one ambiguity-screened spherical centroid crosswalk.

    A candidate must match area within a factor 1.8 and centroid distance
    under a predeclared size-adaptive bound. Reused source IDs rejected.
    """
    if not source:
        return {}, {"matched": 0, "reason": "empty source"}
    features = xyz([r["lat"] for r in source], [r["lon"] for r in source])
    result = {}
    uncertain = {"no_candidate": 0, "ambiguous": 0, "source_reused": 0}
    for begin in range(0, len(target), 80):
        subset = target[begin: begin + 80]
        query = xyz([r["lat"] for r in subset], [r["lon"] for r in subset])
        similarities = query @ features.T
        topn = min(k, len(source))
        indices = np.argpartition(-similarities, topn - 1, axis=1)[:, :topn]
        for local, row in enumerate(subset):
            candidates = []
            limit = max(9.0, min(70.0, 0.32 * math.sqrt(row["area"])))
            for index in indices[local]:
                other = source[int(index)]
                d = RADIUS_KM * math.acos(float(np.clip(similarities[local, index], -1.0, 1.0)))
                ratio = max(row["area"], other["area"]) / min(row["area"], other["area"])
                if d <= limit and ratio <= 1.8:
                    score = d / limit + abs(math.log(row["area"] / other["area"])) / math.log(1.8)
                    candidates.append((score, d, ratio, other))
            if not candidates:
                uncertain["no_candidate"] += 1
                continue
            candidates.sort(key=lambda x: (x[0], x[1], x[3]["id"]))
            best = candidates[0]
            # A second plausible island is problematic even if the best has
            # a slightly better score; do not force an arbitrary name join.
            if len(candidates) > 1 and candidates[1][0] - best[0] < 0.20:
                uncertain["ambiguous"] += 1
                continue
            result[row["id"]] = {"source": best[3], "km": best[1],
                                 "area_ratio": best[2], "match_score": best[0],
                                 "candidate_count": len(candidates)}
    reverse = {}
    for island_id, info in result.items():
        reverse.setdefault(info["source"]["id"], []).append(island_id)
    for ids in reverse.values():
        if len(ids) > 1:
            for island_id in ids:
                del result[island_id]
                uncertain["source_reused"] += 1
    audit = {"target_total": len(target), "source_total": len(source),
             "matched": len(result), **uncertain,
             "max_area_ratio": 1.8, "min_score_gap_second_best": 0.20,
             "max_centroid_km": "max(9,min(70,0.32*sqrt(area_km2)))",
             "note": "conservative centroid+area join, NOT verified polygon identity"}
    return result, audit


def run_crosswalk(destination: Path) -> dict:
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        root = Path(work)
        files = {"gshhg": root / "gshhg.csv", "weigelt": root / "weigelt.csv",
                 "zizka": root / "zizka.zip", "roeble": root / "roeble.xlsx"}
        hashes = {"gshhg": download(GSHHG_URL, files["gshhg"], 80_000_000),
                  "weigelt": download(WEIGELT, files["weigelt"], 15_000_000),
                  "zizka": download(ZIZKA, files["zizka"], 100_000_000),
                  "roeble": download(ROEBLE, files["roeble"], 8_000_000)}
        geo = parse_gshhg(files["gshhg"])
        conn = parse_weigelt(files["weigelt"])
        verified = parse_verified_geology(files["zizka"])
        roeble = parse_roeble_xlsx(files["roeble"])
    conn_match, conn_audit = crosswalk(geo, conn)
    geo_match, geo_audit = crosswalk(geo, verified)
    roeble_match, roeble_audit = crosswalk(geo, roeble)
    source = {r["id"]: r for r in geo}
    rows = []
    for island in geo:
        matched = conn_match.get(island["id"])
        zmatch = geo_match.get(island["id"])
        rmatch = roeble_match.get(island["id"])
        conflict = bool(zmatch and rmatch and
                        zmatch["source"]["geological_origin"] != rmatch["source"]["geological_origin"])
        strict = None if conflict else (rmatch or zmatch)
        direct = strict["source"] if strict else {}
        connected = matched["source"] if matched else {}
        label = direct.get("geological_origin", "unresolved")
        if label not in ("continental", "oceanic", "mixed"):
            label = "unresolved"
        rows.append({
            "island_id": island["id"],
            "analysis_regime": island["regime"],
            "area_km2": island["area"],
            "spatial_block": island["spatial_block"],
            **{f"climate_pc{j}": island[f"climate_pc{j}"] for j in range(1, 5)},
            "distance_to_continent_km": island["continent_distance_km"],
            "island_latitude": island["lat"], "island_longitude": island["lon"],
            "geological_origin": label,
            "geological_class_source": ("Roeble_2024_510_island_geology" if rmatch and not conflict
                                         else "Zizka_2022_missing_geology_additions" if zmatch and not conflict
                                         else "conflicted_sources" if conflict else ""),
            "geological_record_id": direct.get("id", ""),
            "geological_reference": direct.get("geological_reference", ""),
            "geological_match_km": strict["km"] if strict else "",
            "geological_area_ratio": strict["area_ratio"] if strict else "",
            "LGM_connectivity": connected.get("connection", "unresolved"),
            "weigelt_record_id": connected.get("id", ""),
            "weigelt_island": connected.get("name", ""),
            "weigelt_archipelago": connected.get("archipelago", ""),
            "weigelt_match_km": matched["km"] if matched else "",
            "weigelt_area_ratio": matched["area_ratio"] if matched else "",
        })
    columns = list(rows[0])
    with (destination / "island_origin_crosswalk.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    from collections import Counter
    geology_counts = dict(Counter(r["geological_origin"] for r in rows))
    connectivity_counts = dict(Counter(r["LGM_connectivity"] for r in rows))
    contradictions = [r["island_id"] for r in rows if r["geological_origin"] == "oceanic"
                      and r["LGM_connectivity"] == "lgm_connected"]
    audit = {
        "status": "geological_origin_partial_and_lgm_proxy_separate",
        "target_gshhg": len(geo), "weigelt_source_records": len(conn),
        "zizka_verified_additional_records": len(verified),
        "roeble_2024_geological_records": len(roeble),
        "roeble_match_audit": roeble_audit,
        "roeble_zizka_disagreement": sum(bool(geo_match.get(r["id"]) and roeble_match.get(r["id"]) and
                geo_match[r["id"]]["source"]["geological_origin"] != roeble_match[r["id"]]["source"]["geological_origin"])
                for r in geo),
        "geological_counts": geology_counts, "LGM_counts": connectivity_counts,
        "geological_match_audit": geo_audit, "LGM_match_audit": conn_audit,
        "origin_LGM_inconsistency_count": len(contradictions),
        "origin_LGM_inconsistency_ids": contradictions[:20],
        "source_sha256": hashes,
        "interpretation_limit": "LGM-unconnected does not imply geologically oceanic",
        "source_links": {"gshhg": GSHHG_URL, "weigelt": WEIGELT, "zizka": ZIZKA,
                         "roeble": ROEBLE},
    }
    (destination / "crosswalk_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"status": audit["status"], "geological_counts": geology_counts,
                      "LGM_counts": connectivity_counts,
                      "geology_matched": len(geo_match), "roeble_matched": len(roeble_match),
                      "LGM_matched": len(conn_match),
                      "conflicting": len(contradictions)}, ensure_ascii=False))
    return audit
