#!/usr/bin/env python3
"""Reproduce the separate floral-profile exploration from frozen source ledgers.

The upstream population is retained in prepared CSVs, including names with no
eligible traits. The analysis applies its lexical taxon-name gate separately.
All four evidence/single-state scopes use the restored All reference space.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys

import pandas as pd


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from island_v2.floral_space import CORE_TRAITS, analyze  # noqa: E402
from island_v2.floral_space_data import TRAITS, TRAIT_SETS, load_profiles  # noqa: E402


REFERENCE_SHA256 = {
    "raw_direct": "7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
    "species_axis_coverage": "a6eef8d731b2730a99e388ddd0683d7ac9b2af30d43f1f9569c857f89f664b2a",
    "snapshot_all": "2d95136fc7df1dcce92aab738ce0364cdb1434746e1f7b794c529900192cbf4d",
    "snapshot_direct": "d7d6a0372ce7f3ebe958047a222ae0504f01564e364f93bab360273235c6661b",
    "flora": "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
UPSTREAM_REFERENCE_ARTIFACTS = {
    "trait_source": {
        "repository": "zuizui0223/island",
        "run_id": 34191508045,
        "artifact_id": 10042377593,
        "name": "source-scale-batch-integration-34191508045",
        "run_url": "https://github.com/zuizui0223/island/actions/runs/34191508045",
        "archive_sha256": "5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd",
    },
    "progressive": {
        "repository": "zuizui0223/island",
        "run_id": 34232450884,
        "artifact_id": 10058653212,
        "name": "chapter1-progressive-analysis-34232450884",
        "run_url": "https://github.com/zuizui0223/island/actions/runs/34232450884",
        "archive_sha256": "b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999",
    },
}


def file_record(path: Path) -> dict:
    """Hash the actual bytes, so alternate inputs are never attributed to a reference."""
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "sha256": digest}


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def source_paths(input_root: Path, geography: Path, ontology: Path) -> dict[str, Path]:
    return {
        "raw_direct": input_root / "trait_source/direct_species_trait_ledger.csv.gz",
        "species_axis_coverage": input_root / "trait_source/species_axis_coverage.csv.gz",
        "snapshot_all": input_root / "progressive/snapshot/chapter1_trait_ledger_all_analysis.csv.gz",
        "snapshot_direct": input_root / "progressive/snapshot/chapter1_trait_ledger_direct_only.csv.gz",
        "flora": input_root / "progressive/fixed/canonical/input/chapter1_status_flora.csv.gz",
        "geography": geography,
        "ontology": ontology,
    }


def prepare_profiles(paths: dict[str, Path], *, expected_species: int) -> tuple[dict, dict]:
    """Validate the frozen denominator and prepare four scopes without dropping rows."""
    population_table = pd.read_csv(paths["species_axis_coverage"], dtype=str).fillna("")
    names = population_table["accepted_species"].str.strip()
    if names.eq("").any():
        raise ValueError("Frozen species coverage has blank accepted_species")
    population_table["accepted_species"] = names
    if "axis" in population_table and population_table.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("Frozen species coverage has duplicate species x axis keys")
    population = pd.Index(sorted(names.unique()), name="accepted_species")
    if len(population) != expected_species:
        raise ValueError(f"Expected {expected_species} frozen name units; found {len(population)}")
    all_profiles, direct_profiles, audit = load_profiles(
        paths["raw_direct"], paths["snapshot_all"], paths["snapshot_direct"], paths["ontology"]
    )
    columns = [*TRAITS, *[f"{trait}__quality" for trait in TRAITS]]
    profiles = {}
    for scope, frame in (("all", all_profiles), ("direct", direct_profiles)):
        outside = sorted(set(frame["accepted_species"]) - set(population))
        if outside:
            raise ValueError(f"{scope}: {len(outside)} ledger names outside frozen population: {outside[:5]}")
        if frame["accepted_species"].duplicated().any():
            raise ValueError(f"{scope}: duplicate profile species key")
        profiles[scope] = frame.set_index("accepted_species")[columns].reindex(population).reset_index()
    for scope in ("all", "direct"):
        single = profiles[scope].copy()
        multistate = pd.Series(False, index=single.index)
        for trait in CORE_TRAITS:
            multistate |= single[trait].fillna("").str.contains("|", regex=False)
        single.loc[multistate, list(CORE_TRAITS)] = pd.NA
        profiles[f"{scope}_single_state"] = single
    audit["population_species"] = len(population)
    audit["population_validation"] = {
        "coverage_axis_rows": len(population_table),
        "ledger_name_union": len(all_profiles),
        "names_outside_population": 0,
        "population_names_without_ledger_rows": len(population) - len(all_profiles),
        "analysis_unit": "taxon_name_unit",
        "taxonomic_exclusion_stage": "analyze; prepared profiles retain the frozen population",
    }
    return profiles, audit


def provenance_manifest(
    paths: dict[str, Path], profile_paths: dict[str, Path], *, input_root: Path,
    population_species: int, parameters: dict,
) -> dict:
    """Record actual inputs/code/software and distinguish them from reference releases."""
    sources = {label: file_record(path) for label, path in paths.items()}
    for label, digest in REFERENCE_SHA256.items():
        sources[label]["reference_sha256"] = digest
        sources[label]["matches_reference_sha256"] = sources[label]["sha256"] == digest
    references = {label: dict(reference) for label, reference in UPSTREAM_REFERENCE_ARTIFACTS.items()}
    for label, reference in references.items():
        archive = input_root / f"{label}.zip"
        if archive.is_file():
            reference["local_archive"] = file_record(archive)
            reference["local_archive"]["matches_reference_sha256"] = (
                reference["local_archive"]["sha256"] == reference["archive_sha256"]
            )
    git = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPOSITORY, capture_output=True, text=True,
    )
    return {
        "contract": "floral_space_exploration_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "prepared",
        "analysis_unit": "taxon_name_unit",
        "population_species": population_species,
        "parameters": parameters,
        "sources": sources,
        "upstream_reference_artifacts": references,
        "prepared_profiles": {scope: file_record(path) for scope, path in profile_paths.items()},
        "source_code": {label: file_record(REPOSITORY / path) for label, path in {
            "runner": "scripts/run_floral_space_exploration.py",
            "data_adapter": "src/island_v2/floral_space_data.py",
            "analysis": "src/island_v2/floral_space.py",
            "geometry": "src/island_v2/floral_space_geometry.py",
        }.items()},
        "repository_head": git.stdout.strip() if git.returncode == 0 else None,
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "libraries": {name: version(name) for name in ("numpy", "pandas", "scipy", "PyYAML")},
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=REPOSITORY / "work/floral_space_inputs")
    parser.add_argument("--output", type=Path, default=REPOSITORY / "results/floral_space_20261106")
    parser.add_argument("--geography", type=Path,
                        default=REPOSITORY / "results/geography_20260924/corrected_geography_covariates.csv")
    parser.add_argument("--ontology", type=Path, default=REPOSITORY / "config/trait_ontology.yml")
    parser.add_argument("--expected-species", type=int, default=106295)
    parser.add_argument("--min-species", type=int, default=50)
    parser.add_argument("--rarefaction-size", type=int, default=50)
    parser.add_argument("--pair-count", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.expected_species < 1 or args.min_species < 1 or args.rarefaction_size < 1 or args.pair_count < 0:
        parser.error("Population, minimum species and rarefaction size must be positive; pairs must be nonnegative")
    paths = source_paths(args.input_root, args.geography, args.ontology)
    for path in paths.values():
        if not path.is_file():
            parser.error(f"Required input not found: {path}")
    profiles, audit = prepare_profiles(paths, expected_species=args.expected_species)
    flora = pd.read_csv(paths["flora"])
    geography = pd.read_csv(paths["geography"])
    args.output.mkdir(parents=True, exist_ok=True)
    profile_paths = {}
    for scope, frame in profiles.items():
        path = args.input_root / f"{scope}_profiles.csv.gz"
        frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
        profile_paths[scope] = path
    _write_json(args.output / "data_quality_audit.json", audit)
    population_species = audit["population_species"]
    pd.DataFrame([
        {"scope": scope, "trait_set": trait_set, "n_species": counts[trait_set],
         "population_species": population_species, "fraction": counts[trait_set] / population_species}
        for scope, counts in audit["coverage"].items() for trait_set in TRAIT_SETS
    ]).to_csv(args.output / "joint_coverage.csv", index=False)
    analysis_parameters = {name: getattr(args, name) for name in
                           ("min_species", "rarefaction_size", "pair_count", "seed")}
    manifest = provenance_manifest(
        paths, profile_paths, input_root=args.input_root, population_species=population_species,
        parameters={
            **analysis_parameters,
            "expected_species": args.expected_species,
            "traits": list(CORE_TRAITS),
            "space_reference": "all",
            "single_state_policy": "mask both core traits if either contains multiple states; retain population rows",
        },
    )
    manifest_path = args.output / "provenance_manifest.json"
    _write_json(manifest_path, manifest)
    summaries = {}
    for scope, frame in profiles.items():
        print(f"Analyzing {scope} in the common All reference space", flush=True)
        scope_output = args.output / scope
        summaries[scope] = analyze(
            frame, flora, geography, scope_output, space_profiles=profiles["all"], **analysis_parameters,
        )
        scope_inputs = {"profiles": profile_paths[scope], "space_profiles": profile_paths["all"],
                        "flora": paths["flora"], "geography": paths["geography"], "ontology": paths["ontology"]}
        _write_json(scope_output / "input_manifest.json", {
            label: file_record(path) for label, path in scope_inputs.items()
        })
    manifest["status"] = "complete"
    manifest["summaries"] = summaries
    _write_json(manifest_path, manifest)
    print(json.dumps({"output": str(args.output), "population_species": population_species,
                      "scopes": list(summaries)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
