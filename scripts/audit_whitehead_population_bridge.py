"""Audit the Whitehead population-level mating-system bridge to Chapter 1.

This is a diagnostic bridge, not a within-lineage evolutionary test. Whitehead et al.
(2018) retain population-level multilocus outcrossing estimates (tm), whereas the active
Chapter 1 trait ledger collapses mating system to one accepted-species state. This audit
asks how many of the Chapter 1 repeated island/source lineages already have population-
level tm data and whether the Whitehead population table itself contains enough geography
to label observations as island versus mainland/source.

If explicit site geography is absent, the bridge is classified as recoverable-but-not-yet-
identified: original study site metadata must be recovered before any isolation effect is
fitted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

EXPECTED_ZIP_SHA256 = "7067064bf2bc37ade116b33dd3dc710f2938a36fc069265b0628499f887658f8"
SPECIES_FILE = "SpeciesTmData_Whitehead_etal.csv"
POP_FILE = "PopTmData_Whitehead_etal.csv"

EXPLICIT_GEO_COLUMNS = {
    "latitude",
    "longitude",
    "lat",
    "lon",
    "long",
    "country",
    "region",
    "state",
    "province",
    "locality",
    "location",
    "site",
    "site_name",
    "population_name",
    "island",
    "island_name",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean_species(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).replace("_", " ").strip().split())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-zip", type=Path, required=True)
    parser.add_argument("--candidate-lineages", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    observed_hash = sha256(args.source_zip)
    if observed_hash != EXPECTED_ZIP_SHA256:
        raise ValueError(f"Whitehead source hash mismatch: {observed_hash}")

    with ZipFile(args.source_zip) as archive:
        with archive.open(SPECIES_FILE) as handle:
            species = pd.read_csv(handle, dtype=str).fillna("")
        with archive.open(POP_FILE) as handle:
            populations = pd.read_csv(handle, dtype=str).fillna("")

    required_species = {"species", "mean_tm", "var_tm", "pops"}
    required_pop = {"species_name", "study_number", "pop_code", "Tm"}
    if missing := required_species - set(species.columns):
        raise ValueError(f"Whitehead species table missing columns: {sorted(missing)}")
    if missing := required_pop - set(populations.columns):
        raise ValueError(f"Whitehead population table missing columns: {sorted(missing)}")

    candidates = pd.read_csv(args.candidate_lineages, dtype=str).fillna("")
    if "accepted_species" not in candidates.columns:
        raise ValueError("candidate lineage table lacks accepted_species")

    candidate_set = set(candidates["accepted_species"].map(clean_species))
    species = species.copy()
    species["accepted_species"] = species["species"].map(clean_species)
    species["mean_tm_numeric"] = pd.to_numeric(species["mean_tm"], errors="coerce")
    species["var_tm_numeric"] = pd.to_numeric(species["var_tm"], errors="coerce")
    species["pops_numeric"] = pd.to_numeric(species["pops"], errors="coerce")
    species = species.loc[
        species["accepted_species"].ne("")
        & species["mean_tm_numeric"].notna()
        & species["pops_numeric"].notna()
    ].drop_duplicates("accepted_species")

    overlap_species = species.loc[
        species["accepted_species"].isin(candidate_set)
    ].copy()
    overlap_species = overlap_species.merge(
        candidates,
        on="accepted_species",
        how="left",
        validate="one_to_one",
    )
    overlap_species = overlap_species.sort_values(
        ["var_tm_numeric", "n_native_nonendemic_islands"],
        ascending=[False, False],
    )

    populations = populations.copy()
    populations["accepted_species"] = populations["species_name"].map(clean_species)
    populations["tm_numeric"] = pd.to_numeric(populations["Tm"], errors="coerce")
    overlap_populations = populations.loc[
        populations["accepted_species"].isin(
            set(overlap_species["accepted_species"])
        )
    ].copy()

    population_column_map = {str(c).strip().casefold(): str(c) for c in populations.columns}
    geography_columns = sorted(
        population_column_map[key]
        for key in EXPLICIT_GEO_COLUMNS
        if key in population_column_map
    )
    has_explicit_geography = bool(geography_columns)
    n_valid_tm = int(overlap_populations["tm_numeric"].notna().sum())

    per_species_population_counts = (
        overlap_populations.groupby("accepted_species", as_index=False)
        .agg(
            n_population_rows=("pop_code", "size"),
            n_valid_tm=("tm_numeric", lambda x: int(x.notna().sum())),
            n_studies=("study_number", "nunique"),
            tm_min=("tm_numeric", "min"),
            tm_max=("tm_numeric", "max"),
        )
    )
    overlap_species = overlap_species.merge(
        per_species_population_counts,
        on="accepted_species",
        how="left",
        validate="one_to_one",
    )
    overlap_species["tm_range"] = (
        overlap_species["tm_max"] - overlap_species["tm_min"]
    )
    overlap_species = overlap_species.sort_values(
        ["tm_range", "var_tm_numeric", "n_native_nonendemic_islands"],
        ascending=[False, False, False],
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    overlap_species.to_csv(
        args.output_dir / "whitehead_chapter1_overlap_species.csv",
        index=False,
    )
    overlap_populations.to_csv(
        args.output_dir / "whitehead_chapter1_overlap_population_rows.csv",
        index=False,
    )

    summary = {
        "contract": "whitehead_population_bridge_audit_v1",
        "source_zip_sha256": observed_hash,
        "source_population_rows": int(len(populations)),
        "source_species_with_numeric_mean_tm": int(len(species)),
        "chapter1_repeated_lineage_candidates": int(len(candidates)),
        "overlap_species": int(len(overlap_species)),
        "overlap_population_rows": int(len(overlap_populations)),
        "overlap_population_rows_with_numeric_tm": n_valid_tm,
        "population_table_columns": [str(x) for x in populations.columns],
        "explicit_geography_columns": geography_columns,
        "has_explicit_population_geography": has_explicit_geography,
        "within_species_island_source_contrast_estimable_directly": has_explicit_geography,
        "next_gate": (
            "fit population-level island/source contrast"
            if has_explicit_geography
            else "recover study-specific population locality metadata from original sources"
        ),
        "claim_boundary": (
            "Population-level mating-system variation is available for overlapping lineages, "
            "but without explicit population geography it cannot yet be attributed to island "
            "isolation or called within-lineage island evolution."
        ),
    }
    (args.output_dir / "whitehead_population_bridge_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    print("\nTop overlap species by observed tm range:")
    print(
        overlap_species[
            [
                "accepted_species",
                "n_population_rows",
                "tm_min",
                "tm_max",
                "tm_range",
                "var_tm_numeric",
                "n_native_nonendemic_islands",
                "n_mainland_source_entities",
            ]
        ].head(40).to_string(index=False)
    )


if __name__ == "__main__":
    main()
