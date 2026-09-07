from __future__ import annotations

import argparse
import json
import re
import zipfile
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
SOURCE_ID = "monnet_etal_2025_zenodo_selfing"
SOURCE_DOI = "10.5281/zenodo.15288584"
SOURCE_URL = "https://zenodo.org/records/15288584"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
SELFING_COL = re.compile(r"self(?:ing)?(?:_?rate|_?estimate|_?estimates)?$|selfing|self_rate", re.I)


def norm_col(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().casefold()).strip("_")


def normalize_species(value: object) -> str:
    text = " ".join(str(value).replace("_", " ").split())
    return text if BINOMIAL.fullmatch(text) else ""


def read_delimited(path: Path) -> pd.DataFrame:
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=str, sep=None, engine="python").fillna("")
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, dtype=str, sep="\t").fillna("")
    if suffix == ".txt":
        return pd.read_csv(path, dtype=str, sep=None, engine="python").fillna("")
    if suffix == ".xlsx":
        return pd.read_excel(path, dtype=str, engine="openpyxl").fillna("")
    if suffix == ".xls":
        return pd.read_excel(path, dtype=str, engine="xlrd").fillna("")
    raise ValueError(f"unsupported table: {path.name}")


def species_columns(frame: pd.DataFrame) -> list[str]:
    preferred = []
    for column in frame.columns:
        key = norm_col(column)
        if key in {"species", "species_name", "speciesname", "scientific_name", "scientificname", "taxon", "taxon_name"}:
            preferred.append(str(column))
        elif "species" in key and not any(token in key for token in {"number", "count", "pair", "n_species"}):
            preferred.append(str(column))
    return list(dict.fromkeys(preferred))


def selfing_columns(frame: pd.DataFrame) -> list[str]:
    out = []
    for column in frame.columns:
        key = norm_col(column)
        if SELFING_COL.search(key) and not any(token in key for token in {"method", "source", "reference", "citation"}):
            out.append(str(column))
    return out


def numeric_rate(value: object) -> tuple[float | None, str]:
    text = str(value).strip()
    if not text:
        return None, "empty"
    try:
        raw = float(text)
    except ValueError:
        return None, "non_numeric"
    if 0.0 <= raw <= 1.0:
        return raw, "fraction"
    if 1.0 < raw <= 100.0:
        return raw / 100.0, "percent_scaled"
    return None, "out_of_range"


def mating_state(rate: float) -> str:
    if rate < 0.2:
        return "predominantly_outcrossing"
    if rate > 0.8:
        return "predominantly_selfing"
    return "mixed_mating"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    reproductive = coverage.loc[coverage["axis"].eq(AXIS)].copy()
    universe = set(reproductive["accepted_species"])
    unresolved = set(reproductive.loc[reproductive["quality"].eq(""), "accepted_species"])
    if len(universe) != 106295:
        raise ValueError(f"expected 106295 fixed species, got {len(universe)}")

    extract_root = args.output / "source_files"
    extract_root.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.source_zip) as zf:
        zf.extractall(extract_root)

    table_inventory: list[dict[str, object]] = []
    row_candidates: list[dict[str, object]] = []
    unmatched_names: set[str] = set()
    files = sorted(p for p in extract_root.rglob("*") if p.is_file())

    for path in files:
        relative = str(path.relative_to(extract_root))
        if path.suffix.casefold() not in {".csv", ".tsv", ".tab", ".txt", ".xlsx", ".xls"}:
            table_inventory.append({
                "file": relative, "parse_status": "non_table", "rows": "", "columns": "",
                "species_columns": "", "selfing_columns": "", "candidate_numeric_rows": 0,
            })
            continue
        try:
            frame = read_delimited(path)
        except Exception as exc:
            table_inventory.append({
                "file": relative, "parse_status": f"error:{type(exc).__name__}", "rows": 0,
                "columns": "", "species_columns": "", "selfing_columns": "",
                "candidate_numeric_rows": 0,
            })
            continue
        sp_cols = species_columns(frame)
        sf_cols = selfing_columns(frame)
        n_candidates_before = len(row_candidates)
        for row_index, row in frame.iterrows():
            species = ""
            source_species = ""
            species_col = ""
            for column in sp_cols:
                candidate = normalize_species(row[column])
                if candidate:
                    species = candidate
                    source_species = " ".join(str(row[column]).replace("_", " ").split())
                    species_col = column
                    break
            if not species:
                continue
            if species not in universe:
                unmatched_names.add(species)
                continue
            for column in sf_cols:
                rate, normalization = numeric_rate(row[column])
                if rate is None:
                    continue
                row_candidates.append({
                    "source_id": SOURCE_ID,
                    "source_doi": SOURCE_DOI,
                    "source_url": SOURCE_URL,
                    "source_file": relative,
                    "source_row": int(row_index) + 2,
                    "source_species_name": source_species,
                    "accepted_species": species,
                    "species_column": species_col,
                    "selfing_column": column,
                    "raw_selfing_value": str(row[column]),
                    "selfing_rate": rate,
                    "rate_normalization": normalization,
                    "axis": AXIS,
                    "trait_name": "mating_system",
                    "proposed_value": mating_state(rate),
                    "public_restart_cell_status": "unresolved" if species in unresolved else "already_resolved",
                    "evidence_scope": "source_dataset_numeric_species_row",
                    "review_status": "source_scale_extracted_pending_method_and_identity_review",
                    "quality": "unreviewed",
                    "promotion_allowed": False,
                    "genus_rule_training_allowed": False,
                })
        table_inventory.append({
            "file": relative,
            "parse_status": "ok",
            "rows": len(frame),
            "columns": "|".join(str(c) for c in frame.columns),
            "species_columns": "|".join(sp_cols),
            "selfing_columns": "|".join(sf_cols),
            "candidate_numeric_rows": len(row_candidates) - n_candidates_before,
        })

    candidates = pd.DataFrame(row_candidates)
    candidate_columns = [
        "source_id", "source_doi", "source_url", "source_file", "source_row",
        "source_species_name", "accepted_species", "species_column", "selfing_column",
        "raw_selfing_value", "selfing_rate", "rate_normalization", "axis", "trait_name",
        "proposed_value", "public_restart_cell_status", "evidence_scope", "review_status",
        "quality", "promotion_allowed", "genus_rule_training_allowed",
    ]
    if candidates.empty:
        candidates = pd.DataFrame(columns=candidate_columns)
    else:
        candidates = candidates[candidate_columns].drop_duplicates()
    candidates.to_csv(args.output / "zenodo_selfing_row_candidates.csv", index=False)
    pd.DataFrame(table_inventory).to_csv(args.output / "zenodo_source_table_inventory.csv", index=False)
    pd.DataFrame({"source_species_name": sorted(unmatched_names)}).to_csv(
        args.output / "zenodo_unmatched_species_names.csv", index=False
    )

    aggregate_rows: list[dict[str, object]] = []
    if len(candidates):
        for species, group in candidates.groupby("accepted_species", sort=True):
            rates = sorted(set(float(v) for v in group["selfing_rate"]))
            states = sorted(set(str(v) for v in group["proposed_value"]))
            aggregate_rows.append({
                "accepted_species": species,
                "axis": AXIS,
                "trait_name": "mating_system",
                "n_numeric_source_rows": len(group),
                "n_unique_selfing_rates": len(rates),
                "selfing_rate_min": min(rates),
                "selfing_rate_max": max(rates),
                "selfing_rate_median": float(pd.Series(rates).median()),
                "proposed_state_set": "|".join(states),
                "single_state_across_source_rows": len(states) == 1,
                "proposed_value": states[0] if len(states) == 1 else "",
                "public_restart_cell_status": (
                    "unresolved" if species in unresolved else "already_resolved"
                ),
                "review_status": "pending_method_identity_and_duplicate_source_review",
                "quality": "unreviewed",
                "promotion_allowed": False,
                "genus_rule_training_allowed": False,
            })
    aggregate = pd.DataFrame(aggregate_rows)
    if aggregate.empty:
        aggregate = pd.DataFrame(columns=[
            "accepted_species", "axis", "trait_name", "n_numeric_source_rows",
            "n_unique_selfing_rates", "selfing_rate_min", "selfing_rate_max",
            "selfing_rate_median", "proposed_state_set", "single_state_across_source_rows",
            "proposed_value", "public_restart_cell_status", "review_status", "quality",
            "promotion_allowed", "genus_rule_training_allowed",
        ])
    aggregate.to_csv(args.output / "zenodo_selfing_species_candidates.csv", index=False)

    unresolved_candidates = aggregate.loc[
        aggregate["public_restart_cell_status"].eq("unresolved")
    ] if len(aggregate) else aggregate
    summary = {
        "contract": "monnet_2025_zenodo_selfing_source_scale_v1",
        "source_doi": SOURCE_DOI,
        "source_archive_complete": True,
        "files_in_archive": len(files),
        "parsed_tables": sum(row["parse_status"] == "ok" for row in table_inventory),
        "numeric_species_rows_in_fixed_universe": len(candidates),
        "fixed_universe_species_with_numeric_selfing": int(aggregate["accepted_species"].nunique()) if len(aggregate) else 0,
        "exact_unresolved_reproductive_species_overlap": int(unresolved_candidates["accepted_species"].nunique()) if len(unresolved_candidates) else 0,
        "unresolved_species_single_proposed_mating_state": int(
            unresolved_candidates["single_state_across_source_rows"].astype(bool).sum()
        ) if len(unresolved_candidates) else 0,
        "unresolved_species_variable_across_state_thresholds": int(
            (~unresolved_candidates["single_state_across_source_rows"].astype(bool)).sum()
        ) if len(unresolved_candidates) else 0,
        "formal_gain": 0,
        "promotion_allowed": False,
        "mapping_note": "Selfing rate is screened only toward mating_system. It is never used to infer self-incompatibility or autonomous selfing. Values <0.2, 0.2-0.8, >0.8 are planning bins matching the prior outcrossing-rate convention after inversion; all rows remain unreviewed until source methods, identity, duplicate source lineage and species-level interpretation are checked.",
    }
    (args.output / "zenodo_selfing_source_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
