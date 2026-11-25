"""Run the public-aggregate first stage of fd's floral functional-diversity line.

No original species-to-trait records are loaded, copied or released.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

from .core import ProfileSpace, csv_rows, load_carriers, load_island_metadata, load_profiles


SCOPES = ("all", "direct", "all_single_state", "direct_single_state")
UPSTREAM_COMMIT = "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c"
PAIR_INPUT = Path("floral_space_pairs_20261106") / "pairwise_observed_and_null.csv.gz"
REF_INPUT = Path("floral_space_20261106") / "all" / "profile_coordinates.csv"
ISLAND_COLUMNS = (
    "scope", "island_id", "analysis_regime", "n_names", "n_profiles",
    "rao_q", "rao_q_without_replacement", "rao_q_profile_equal",
    "joint_coverage", "area_km2", "distance_to_continent_km",
)
PAIR_COLUMNS = (
    "scope", "encoding", "analysis_regime", "island_a", "island_b",
    "n_names_a", "n_names_b", "n_profiles_a", "n_profiles_b",
    "profile_jaccard", "nn_ab", "nn_ba", "functional_nn_dissimilarity",
    "mean_cross_profile_distance",
)


def file_sha256(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            sha.update(chunk)
    return sha.hexdigest()


def _write_csv(path: Path, rows: list[dict], columns: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(input_root: Path, output_root: Path, *, min_names: int = 50) -> dict:
    if isinstance(min_names, bool) or min_names < 1:
        raise ValueError("min_names must be positive")
    input_root = Path(input_root)
    output_root = Path(output_root)
    reference = input_root / REF_INPUT
    pairs = input_root / PAIR_INPUT
    sources = {"reference": reference, "pairs": pairs}
    for scope in SCOPES:
        for name in ("island_metrics.csv", "island_profile_carriers.csv.gz"):
            sources[f"{scope}/{name}"] = input_root / "floral_space_20261106" / scope / name
    missing = [str(p) for p in sources.values() if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing required public-aggregate files: " + "; ".join(missing))

    profile_reference = load_profiles(reference)
    space = ProfileSpace(profile_reference)
    output_root.mkdir(parents=True, exist_ok=True)
    summaries = {}
    staged: dict[str, tuple[dict, dict]] = {}

    for scope in SCOPES:
        base = input_root / "floral_space_20261106" / scope
        carriers = load_carriers(base / "island_profile_carriers.csv.gz")
        metadata = load_island_metadata(base / "island_metrics.csv")
        if set(carriers) != set(metadata):
            raise ValueError(f"{scope}: carrier and metadata island populations disagree")
        rows = []
        eligible = {}
        for island in sorted(metadata):
            c = carriers[island]
            metric = space.island(c)
            meta = metadata[island]
            if metric["n_names"] != int(meta["n_species"]):
                raise ValueError(f"{scope}: name count mismatch for {island}")
            if metric["n_profiles"] != int(meta["n_profiles"]):
                raise ValueError(f"{scope}: profile count mismatch for {island}")
            if metric["n_names"] < min_names:
                continue
            eligible[island] = c
            rows.append({
                "scope": scope, "island_id": island,
                "analysis_regime": meta["analysis_regime"],
                **metric,
                "joint_coverage": meta.get("joint_coverage", ""),
                "area_km2": meta.get("area_km2", ""),
                "distance_to_continent_km": meta.get("distance_to_continent_km", ""),
            })
        staged[scope] = (eligible, metadata)
        _write_csv(output_root / f"{scope}_islands.csv", rows, ISLAND_COLUMNS)
        summaries[scope] = {
            "n_islands_eligible": len(rows),
            "n_islands_with_any_complete_profile": len(metadata),
            "n_records_eligible": sum(int(r["n_names"]) for r in rows),
            "n_pairs": 0,
        }

    seen: dict[str, set[tuple[str, str]]] = {scope: set() for scope in SCOPES}
    pair_rows: dict[str, list[dict]] = {scope: [] for scope in SCOPES}
    for source in csv_rows(pairs):
        scope = source["evidence_scope"]
        if scope not in SCOPES:
            raise ValueError(f"Unknown evidence scope in matched pairs: {scope}")
        a, b = source["island_a"], source["island_b"]
        encoding, regime = source["encoding"], source["analysis_regime"]
        if a == b:
            raise ValueError("Self-pair in PR321 matched-cohort input")
        key = tuple(sorted((a, b)))
        if key in seen[scope]:
            raise ValueError(f"Duplicate {scope} pair: {key}")
        seen[scope].add(key)
        eligible, meta = staged[scope]
        if a not in eligible or b not in eligible:
            # Explicit support gate: do not impute absent endpoint profiles.
            continue
        if (meta[a]["analysis_regime"], meta[b]["analysis_regime"]) != (regime, regime):
            raise ValueError(f"Pair and island regimes disagree for {a}/{b}")
        scores = space.pair(eligible[a], eligible[b])
        pair_rows[scope].append({
            "scope": scope, "encoding": encoding, "analysis_regime": regime,
            "island_a": a, "island_b": b, **scores,
        })

    for scope in SCOPES:
        _write_csv(output_root / f"{scope}_pairs.csv", pair_rows[scope], PAIR_COLUMNS)
        summaries[scope]["n_pairs"] = len(pair_rows[scope])

    manifest = {
        "schema_version": "fd_flower_functional_diversity_v1",
        "status": "computed_descriptive_not_mechanistic",
        "upstream_repository": "zuizui0223/island",
        "upstream_commit": UPSTREAM_COMMIT,
        "upstream_pull_request": 321,
        "input_source": "public aggregate island-profile carriers and matching PR321 pair IDs",
        "data_scope": "all recorded origins; name units, not validated botanical species",
        "method": {
            "traits": ["flower_primary_color", "floral_form"],
            "distance": "equal-weight mean state-set Jaccard dissimilarity, ORIGINAL distance",
            "raos_q": "species-name-carrier-weighted mean within-island pairwise distance",
            "pair_nn": "symmetrized carrier-weighted nearest-profile distance, whole-island",
            "shared_taxa_removed": False,
            "genus_conditioning": False,
            "significance_claim": False,
            "min_names": min_names,
        },
        "dependencies": {"python": platform.python_version()},
        "source_files_sha256": {name: file_sha256(path) for name, path in sources.items()},
        "reference_profiles": len(profile_reference),
        "summaries": summaries,
        "release": "aggregate-only outputs; do not add original taxon x trait ledgers",
    }
    (output_root / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-names", type=int, default=50)
    args = parser.parse_args()
    result = run(args.input, args.output, min_names=args.min_names)
    print(json.dumps({"status": result["status"], "reference_profiles": result["reference_profiles"],
                      "summaries": result["summaries"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
