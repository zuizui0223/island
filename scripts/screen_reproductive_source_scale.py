from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import tempfile
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
REPRO_COL = re.compile(
    r"outcross|selfing|self[_ .-]?compat|self[_ .-]?incompat|breeding|mating|"
    r"autofertil|autogam|cleistog|(^|[_.])si($|[_.])|isi|mean[_.]?tm|^tm$",
    re.I,
)

SOURCES = [
    {
        "source_id": "goodwillie_kalisz_eckert_2012_dryad",
        "doi": "10.5061/dryad.292q34fp",
        "status": "integrated_control",
        "source_type": "structured_xls",
        "note": "Already source-scale processed on 2026-08-21 (469 source rows); retain only as a control.",
    },
    {
        "source_id": "goodwillie_etal_2010_floral_display_dryad",
        "doi": "10.5061/dryad.bn5v3cm5",
        "status": "candidate",
        "source_type": "structured_xls",
        "note": "Outcrossing-rate plus floral-display tables; only reproductive fields are screened here.",
    },
    {
        "source_id": "busch_prior_2021_dryad",
        "doi": "10.5061/dryad.3j9kd51jw",
        "status": "candidate",
        "source_type": "structured_csv_xlsx",
        "note": "Population selfing-rate compilation; species-level aggregation requires later review.",
    },
    {
        "source_id": "moeller_etal_2017_dryad",
        "doi": "10.5061/dryad.577q1",
        "status": "integrated_control",
        "source_type": "structured_csv",
        "note": "Already source-scale processed in the repository; retained as a control for the screen.",
    },
    {
        "source_id": "razanajatovo_2018_dryad",
        "doi": "10.5061/dryad.89pd54f",
        "status": "integrated_control",
        "source_type": "structured_csv",
        "note": "Already integrated self-compatibility/autofertility source; retained as a control.",
    },
]


def _norm_col(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().casefold()).strip("_")


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str).fillna("")
    if suffix == ".xlsx":
        return pd.read_excel(path, dtype=str, engine="openpyxl").fillna("")
    if suffix == ".xls":
        return pd.read_excel(path, dtype=str, engine="xlrd").fillna("")
    raise ValueError(suffix)


def _species_from_table(frame: pd.DataFrame) -> set[str]:
    if frame.empty:
        return set()
    columns = {_norm_col(c): c for c in frame.columns}
    direct_names = [
        "accepted_species",
        "scientific_name",
        "scientificname",
        "species_name",
        "speciesname",
        "taxon_name",
        "taxon",
        "species",
    ]
    out: set[str] = set()
    for key in direct_names:
        if key not in columns:
            continue
        for raw in frame[columns[key]].astype(str):
            value = " ".join(raw.replace("_", " ").split())
            if BINOMIAL.fullmatch(value):
                out.add(value)
    if "genus" in columns and "species" in columns:
        for genus, epithet in zip(
            frame[columns["genus"]].astype(str),
            frame[columns["species"]].astype(str),
            strict=False,
        ):
            g = " ".join(genus.split())
            e = " ".join(epithet.split())
            candidate = f"{g} {e}"
            if BINOMIAL.fullmatch(candidate):
                out.add(candidate)
    return out


def _reproductive_columns(frame: pd.DataFrame) -> list[str]:
    return [str(c) for c in frame.columns if REPRO_COL.search(_norm_col(c))]


def _download_dryad(doi: str, target: Path) -> Path:
    encoded = urllib.parse.quote(f"doi:{doi}", safe="")
    url = f"https://datadryad.org/api/v2/datasets/{encoded}/download"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "island-source-scale-screen/1.0"},
    )
    with urllib.request.urlopen(request, timeout=180) as response, target.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return target


def _unpack_recursive(archive: Path, root: Path) -> list[Path]:
    pending = [archive]
    seen: set[Path] = set()
    tables: list[Path] = []
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        suffix = path.suffix.casefold()
        if suffix in {".csv", ".xls", ".xlsx"}:
            tables.append(path)
            continue
        if suffix == ".zip" or zipfile.is_zipfile(path):
            out = root / f"unzipped_{len(seen)}"
            out.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(path) as zf:
                zf.extractall(out)
            pending.extend(p for p in out.rglob("*") if p.is_file())
    return tables


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    reproductive = coverage.loc[coverage["axis"].eq(AXIS)].copy()
    unresolved = set(reproductive.loc[reproductive["quality"].eq(""), "accepted_species"])
    universe = set(reproductive["accepted_species"])
    if len(universe) != 106295:
        raise ValueError(f"expected 106295 reproductive species, got {len(universe)}")

    source_rows: list[dict[str, object]] = []
    table_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, str]] = []

    with tempfile.TemporaryDirectory(prefix="island-source-screen-") as tmp:
        root = Path(tmp)
        for source in SOURCES:
            source_root = root / str(source["source_id"])
            source_root.mkdir(parents=True, exist_ok=True)
            error = ""
            total_rows = 0
            all_taxa: set[str] = set()
            relevant_taxa: set[str] = set()
            relevant_tables = 0
            tables: list[Path] = []
            try:
                archive = _download_dryad(str(source["doi"]), source_root / "dataset.zip")
                tables = _unpack_recursive(archive, source_root)
                for table in sorted(tables):
                    try:
                        frame = _read_table(table)
                    except Exception as exc:
                        table_rows.append(
                            {
                                "source_id": source["source_id"],
                                "file": table.name,
                                "rows": 0,
                                "taxa": 0,
                                "reproductive_columns": "",
                                "parse_status": f"error:{type(exc).__name__}",
                            }
                        )
                        continue
                    taxa = _species_from_table(frame)
                    columns = _reproductive_columns(frame)
                    total_rows += len(frame)
                    all_taxa.update(taxa)
                    if columns:
                        relevant_tables += 1
                        relevant_taxa.update(taxa)
                    table_rows.append(
                        {
                            "source_id": source["source_id"],
                            "file": table.name,
                            "rows": len(frame),
                            "taxa": len(taxa),
                            "reproductive_columns": "|".join(columns),
                            "parse_status": "ok",
                        }
                    )
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"

            exact_overlap = relevant_taxa & unresolved
            in_universe = relevant_taxa & universe
            estimated_minutes = 20.0 + 2.0 * relevant_tables + 0.30 * len(exact_overlap)
            cells_per_hour = (len(exact_overlap) / (estimated_minutes / 60.0)) if exact_overlap else 0.0
            source_rows.append(
                {
                    **source,
                    "download_status": "error" if error else "ok",
                    "error": error,
                    "tables_found": len(tables),
                    "reproductive_tables": relevant_tables,
                    "source_rows_scanned": total_rows,
                    "taxa_detected": len(all_taxa),
                    "reproductive_taxa_detected": len(relevant_taxa),
                    "reproductive_taxa_in_fixed_universe": len(in_universe),
                    "exact_unresolved_reproductive_overlap": len(exact_overlap),
                    "estimated_review_minutes": round(estimated_minutes, 1),
                    "expected_cells_per_hour_upper_bound": round(cells_per_hour, 2),
                    "formal_gain": 0,
                    "promotion_allowed": False,
                }
            )
            for species in sorted(exact_overlap):
                overlap_rows.append(
                    {
                        "source_id": str(source["source_id"]),
                        "accepted_species": species,
                        "axis": AXIS,
                        "status": "unresolved_overlap_only_not_evidence",
                    }
                )

    ranking = pd.DataFrame(source_rows)
    ranking["priority_class"] = ranking.apply(
        lambda row: (
            "candidate_high_roi"
            if row["status"] == "candidate" and int(row["exact_unresolved_reproductive_overlap"]) > 0
            else ("candidate_no_exact_overlap" if row["status"] == "candidate" else "integrated_control")
        ),
        axis=1,
    )
    ranking = ranking.sort_values(
        ["status", "exact_unresolved_reproductive_overlap", "expected_cells_per_hour_upper_bound"],
        ascending=[True, False, False],
    )
    ranking.to_csv(args.output / "reproductive_source_roi.csv", index=False)
    pd.DataFrame(table_rows).to_csv(args.output / "reproductive_source_table_inventory.csv", index=False)
    pd.DataFrame(overlap_rows, columns=["source_id", "accepted_species", "axis", "status"]).to_csv(
        args.output / "reproductive_source_unresolved_overlap.csv", index=False
    )

    summary = {
        "contract": "reproductive_source_scale_roi_screen_v1",
        "fixed_species": 106295,
        "current_reproductive_unresolved": len(unresolved),
        "sources_screened": len(SOURCES),
        "candidate_sources": sum(s["status"] == "candidate" for s in SOURCES),
        "formal_gain": 0,
        "selection_only": True,
        "top_candidate": (
            ranking.loc[ranking["status"].eq("candidate"), "source_id"].iloc[0]
            if (ranking["status"] == "candidate").any()
            else None
        ),
    }
    (args.output / "reproductive_source_roi_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
