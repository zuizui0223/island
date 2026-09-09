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
REQUIRED_COLUMNS = {"species", "individual", "selfing_rate", "logLik"}


def normalize_species(value: object) -> str:
    text = " ".join(str(value).replace("_", " ").split())
    return text if BINOMIAL.fullmatch(text) else ""


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
    files = sorted(p for p in extract_root.rglob("*") if p.is_file())

    inventory: list[dict[str, object]] = []
    individual_rows: list[dict[str, object]] = []
    unmatched: set[str] = set()
    ambiguous_mle_individuals = 0

    for path in files:
        relative = str(path.relative_to(extract_root))
        try:
            frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        except Exception as exc:
            inventory.append({
                "file": relative, "parse_status": f"error:{type(exc).__name__}",
                "rows": 0, "species": 0, "individuals": 0, "mle_individuals": 0,
            })
            continue
        if not REQUIRED_COLUMNS.issubset(frame.columns):
            inventory.append({
                "file": relative, "parse_status": "unsupported_schema",
                "rows": len(frame), "species": 0, "individuals": 0, "mle_individuals": 0,
                "columns": "|".join(frame.columns),
            })
            continue
        frame["selfing_rate_num"] = pd.to_numeric(frame["selfing_rate"], errors="coerce")
        frame["logLik_num"] = pd.to_numeric(frame["logLik"], errors="coerce")
        frame = frame.loc[
            frame["selfing_rate_num"].between(0, 1, inclusive="both")
            & frame["logLik_num"].notna()
        ].copy()
        mle_count = 0
        for (source_species, individual), group in frame.groupby(["species", "individual"], sort=False):
            maximum = group["logLik_num"].max()
            best = group.loc[group["logLik_num"].eq(maximum), "selfing_rate_num"].dropna().unique()
            species = normalize_species(source_species)
            if not species:
                continue
            if species not in universe:
                unmatched.add(species)
                continue
            if len(best) != 1:
                ambiguous_mle_individuals += 1
                continue
            rate = float(best[0])
            mle_count += 1
            individual_rows.append({
                "source_id": SOURCE_ID,
                "source_doi": SOURCE_DOI,
                "source_url": SOURCE_URL,
                "source_file": relative,
                "source_species_name": " ".join(str(source_species).replace("_", " ").split()),
                "accepted_species": species,
                "individual": individual,
                "mle_selfing_rate": rate,
                "max_logLik": float(maximum),
                "mle_mating_state": mating_state(rate),
                "axis": AXIS,
                "trait_name": "mating_system",
                "public_restart_cell_status": "unresolved" if species in unresolved else "already_resolved",
                "evidence_scope": "source_dataset_individual_selfing_mle",
                "review_status": "source_scale_extracted_pending_method_identity_and_species_aggregation_review",
                "quality": "unreviewed",
                "promotion_allowed": False,
                "genus_rule_training_allowed": False,
            })
        inventory.append({
            "file": relative,
            "parse_status": "ok",
            "rows": len(frame),
            "species": frame["species"].nunique(),
            "individuals": frame[["species", "individual"]].drop_duplicates().shape[0],
            "mle_individuals_in_fixed_universe": mle_count,
            "columns": "|".join(frame.columns),
        })

    individuals = pd.DataFrame(individual_rows)
    individual_columns = [
        "source_id", "source_doi", "source_url", "source_file", "source_species_name",
        "accepted_species", "individual", "mle_selfing_rate", "max_logLik", "mle_mating_state",
        "axis", "trait_name", "public_restart_cell_status", "evidence_scope", "review_status",
        "quality", "promotion_allowed", "genus_rule_training_allowed",
    ]
    if individuals.empty:
        individuals = pd.DataFrame(columns=individual_columns)
    else:
        individuals = individuals[individual_columns].drop_duplicates(
            ["accepted_species", "individual", "source_file"]
        )
    individuals.to_csv(args.output / "zenodo_selfing_individual_mle.csv", index=False)
    pd.DataFrame(inventory).to_csv(args.output / "zenodo_source_table_inventory.csv", index=False)
    pd.DataFrame({"source_species_name": sorted(unmatched)}).to_csv(
        args.output / "zenodo_unmatched_species_names.csv", index=False
    )

    species_rows: list[dict[str, object]] = []
    for species, group in individuals.groupby("accepted_species", sort=True):
        rates = group["mle_selfing_rate"].astype(float)
        states = sorted(set(group["mle_mating_state"].astype(str)))
        species_rows.append({
            "accepted_species": species,
            "axis": AXIS,
            "trait_name": "mating_system",
            "n_individual_mle": len(group),
            "selfing_mle_min": float(rates.min()),
            "selfing_mle_max": float(rates.max()),
            "selfing_mle_mean": float(rates.mean()),
            "selfing_mle_median": float(rates.median()),
            "individual_state_set": "|".join(states),
            "single_state_across_individual_mles": len(states) == 1,
            "proposed_value": states[0] if len(states) == 1 else "",
            "public_restart_cell_status": "unresolved" if species in unresolved else "already_resolved",
            "review_status": "pending_method_identity_duplicate_lineage_and_species_aggregation_review",
            "quality": "unreviewed",
            "promotion_allowed": False,
            "genus_rule_training_allowed": False,
        })
    species_candidates = pd.DataFrame(species_rows)
    if species_candidates.empty:
        species_candidates = pd.DataFrame(columns=[
            "accepted_species", "axis", "trait_name", "n_individual_mle", "selfing_mle_min",
            "selfing_mle_max", "selfing_mle_mean", "selfing_mle_median", "individual_state_set",
            "single_state_across_individual_mles", "proposed_value", "public_restart_cell_status",
            "review_status", "quality", "promotion_allowed", "genus_rule_training_allowed",
        ])
    species_candidates.to_csv(args.output / "zenodo_selfing_species_candidates.csv", index=False)

    unresolved_candidates = species_candidates.loc[
        species_candidates["public_restart_cell_status"].eq("unresolved")
    ] if len(species_candidates) else species_candidates
    stable_unresolved = unresolved_candidates.loc[
        unresolved_candidates["single_state_across_individual_mles"].astype(bool)
    ] if len(unresolved_candidates) else unresolved_candidates

    summary = {
        "contract": "monnet_2025_zenodo_selfing_source_scale_v2",
        "source_doi": SOURCE_DOI,
        "source_archive_complete": True,
        "source_semantics": "Each file contains a selfing-rate likelihood grid per individual. The extracted estimate is the unique selfing_rate at maximum logLik for each individual; grid rows are never treated as observations.",
        "files_in_archive": len(files),
        "parsed_likelihood_files": sum(row["parse_status"] == "ok" for row in inventory),
        "ambiguous_mle_individuals_dropped": ambiguous_mle_individuals,
        "individual_mles_in_fixed_universe": len(individuals),
        "fixed_universe_species_with_individual_mles": int(species_candidates["accepted_species"].nunique()),
        "exact_unresolved_reproductive_species_overlap": int(unresolved_candidates["accepted_species"].nunique()),
        "unresolved_species_single_mating_state_across_individual_mles": int(stable_unresolved["accepted_species"].nunique()),
        "unresolved_species_variable_across_individual_mle_states": int(len(unresolved_candidates) - len(stable_unresolved)),
        "stable_unresolved_state_counts": stable_unresolved["proposed_value"].value_counts().sort_index().astype(int).to_dict() if len(stable_unresolved) else {},
        "formal_gain": 0,
        "promotion_allowed": False,
        "mapping_note": "Selfing estimates are screened only toward mating_system; never SI or autonomous selfing. The <0.2 / 0.2-0.8 / >0.8 bins are planning bins and do not become accepted direct evidence until the paper methods, taxonomic identity, original-source duplication and species-level aggregation rule are reviewed.",
    }
    (args.output / "zenodo_selfing_source_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
