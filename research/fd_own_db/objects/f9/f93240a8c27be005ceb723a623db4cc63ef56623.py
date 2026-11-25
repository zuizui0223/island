"""Rights-restricted geographic cross-fitting coverage audit for floral N2.

Assesses whether a genus/native-mainland trait predictor can be honestly
evaluated without reusing target taxon names. This is an INPUT COVERAGE
DIAGNOSTIC, not a floral packing p-value, ecological response or mainland
ancestry reconstruction. Outputs ONLY grouped counts and digests.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .genus_source_null import load_private_maps, validate_island_counts, _sha
from .pipeline import _write_csv
from .core import load_profiles

SCOPES = ("all", "direct")
FIELDS = (
    "scope", "fold", "n_heldout_islands", "n_spatial_blocks",
    "n_heldout_unique_names", "n_global_donor_names_after_holdout",
    "fraction_focal_carriers_with_global_genus_donor",
    "fraction_focal_carriers_with_local_source_genus_donor",
    "status",
)
UNION_FIELDS = (
    "scope", "eligible_islands", "eligible_distinct_target_names",
    "original_complete_native_GIFT_donor_names",
    "donor_names_after_all_eligible_island_names_removed",
    "fraction_native_donor_names_retained",
    "status",
)


def balanced_block_folds(eligible: pd.DataFrame, *, folds: int = 5) -> dict[str, int]:
    """Deterministic nonoverlapping block assignment, balancing island counts."""
    if not 2 <= folds <= 10:
        raise ValueError("Supported number of cross-fit folds is 2..10")
    if not {"island_id", "spatial_block"}.issubset(eligible):
        raise ValueError("Missing geographic fold columns")
    if eligible.island_id.duplicated().any():
        raise ValueError("Duplicate eligible island")
    if eligible.spatial_block.isna().any() or eligible.spatial_block.astype(str).str.strip().eq("").any():
        raise ValueError("Missing source geographic block")
    counts = eligible.groupby("spatial_block", sort=True).size()
    if len(counts) < folds:
        raise ValueError("Fewer distinct source blocks than folds")
    # stable tie order via deterministic two-key sorting
    order = sorted(counts.items(), key=lambda x: (-x[1], str(x[0])))
    load = [0] * folds
    mapping = {}
    for block, n in order:
        f = min(range(folds), key=lambda i: (load[i], i))
        mapping[str(block)] = f
        load[f] += int(n)
    return mapping


def retained_donor_names(all_donor_names: set[str],
                         heldout_names: set[str]) -> set[str]:
    """Remove the SAME taxon name from all GIFT source entities globally."""
    if not heldout_names:
        raise ValueError("No heldout target names")
    return all_donor_names - heldout_names


def geographic_fold_coverage(
    membership: pd.DataFrame,
    donor_names: set[str],
    entity_names: dict[str, set[str]],
    assignments: dict[str, set[str]],
    block_map: dict[str, int], *, min_islands: int = 5,
) -> list[dict]:
    """Fold-disjoint genus support, both global and geographic source-specific.

    No island names/blocks/source keys leave this function's aggregated rows.
    Columns in membership: island_id, accepted_species, genus, spatial_block.
    """
    required = {"island_id", "accepted_species", "genus", "spatial_block"}
    if not required.issubset(membership):
        raise ValueError("Missing matched island-identity membership fields")
    if membership.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Duplicate island/accepted name in frozen flora")
    working = membership.assign(
        _fold=membership.spatial_block.astype(str).map(block_map))
    if working._fold.isna().any():
        raise ValueError("Unassigned spatial block")
    result = []
    for fold, part in working.groupby("_fold", sort=True):
        targets = set(part.accepted_species)
        trained = retained_donor_names(donor_names, targets)
        global_genera = {name.split()[0] for name in trained}
        global_supported = int(part.genus.isin(global_genera).sum())
        local_supported = 0
        for island, one in part.groupby("island_id", sort=False):
            source_ids = assignments.get(island, set())
            if not source_ids:
                raise ValueError("No potential source assignment for an eligible island")
            candidates = set().union(*(entity_names.get(src, set()) for src in source_ids))
            candidates.intersection_update(trained)
            local_genera = {name.split()[0] for name in candidates}
            local_supported += int(one.genus.isin(local_genera).sum())
        n_islands = part.island_id.nunique()
        report = {
            "fold": int(fold), "n_heldout_islands": n_islands,
            "n_spatial_blocks": int(part.spatial_block.nunique()),
            "n_heldout_unique_names": len(targets),
            "n_global_donor_names_after_holdout": len(trained),
            "status": "heldout_crossfit_donor_support" if n_islands >= min_islands
                      else "small_cell_suppressed",
        }
        if n_islands >= min_islands:
            report.update({
                "fraction_focal_carriers_with_global_genus_donor":
                    global_supported / len(part),
                "fraction_focal_carriers_with_local_source_genus_donor":
                    local_supported / len(part),
            })
        result.append(report)
    return result


def run(private_all: Path, private_direct: Path, flora: Path,
        g0_root: Path, geography: Path, mainland_native: Path,
        assignments: Path, reference: Path, output: Path) -> dict:
    input_paths = list(map(Path, (private_all, private_direct, flora, g0_root,
                                   geography, mainland_native, assignments,
                                   reference)))
    private_all, private_direct, flora, g0_root, geography, mainland_native, assignments, reference = input_paths
    output = Path(output)
    space = load_profiles(reference)
    if len(space) != 496:
        raise ValueError("Frozen original reference not 496 profiles")
    profile_maps = {
        scope: load_private_maps(path, scope, set(space), strict=True)
        for scope, path in (("all", private_all), ("direct", private_direct))
    }
    geography_table = pd.read_csv(geography, usecols=["island_id", "spatial_block"],
                                  dtype=str)
    if geography_table.island_id.duplicated().any():
        raise ValueError("Duplicate island geographic block metadata")
    g0 = {
        scope: pd.read_csv(g0_root / f"{scope}_islands.csv", dtype={"island_id": str})
        for scope in SCOPES
    }
    assign = pd.read_csv(assignments,
                         usecols=["island_id", "source_mode", "entity_ID"],
                         dtype=str)
    assign = assign.loc[assign.source_mode.eq("geo_k5")].drop_duplicates(
        ["island_id", "entity_ID"])
    source_lookup = assign.groupby("island_id").entity_ID.agg(set).to_dict()
    native = pd.read_csv(mainland_native,
                         usecols=["entity_ID", "work_species"], dtype=str)
    native = native.drop_duplicates(["entity_ID", "work_species"])
    flora_table = pd.read_csv(flora,
        usecols=["island_id", "accepted_species"], dtype=str).drop_duplicates()
    if flora_table.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Original island/accepted species duplicate")
    all_g0 = g0["all"].merge(geography_table, on="island_id", validate="one_to_one")
    if len(all_g0) != len(g0["all"]):
        raise ValueError("Incomplete eligible All geographic block coverage")
    block_to_fold = balanced_block_folds(all_g0)
    results, unions, validations = [], [], {}
    for scope in SCOPES:
        selected = g0[scope][["island_id", "n_names", "n_profiles"]].merge(
            geography_table, on="island_id", validate="one_to_one")
        sample = flora_table.loc[
            flora_table.island_id.isin(selected.island_id)
        ].merge(profile_maps[scope], on="accepted_species",
                how="inner", validate="many_to_one")
        validate_island_counts(sample, selected, scope)
        if sample.island_id.nunique() != len(selected):
            raise ValueError("Missing eligible test island membership")
        sample = sample.merge(
            selected[["island_id", "spatial_block"]],
            on="island_id", validate="many_to_one")
        target_union = set(sample.accepted_species)
        donor_names = set(native.work_species) & set(profile_maps[scope].accepted_species)
        remaining = retained_donor_names(donor_names, target_union)
        unions.append({
            "scope": scope, "eligible_islands": len(selected),
            "eligible_distinct_target_names": len(target_union),
            "original_complete_native_GIFT_donor_names": len(donor_names),
            "donor_names_after_all_eligible_island_names_removed": len(remaining),
            "fraction_native_donor_names_retained": len(remaining) / len(donor_names),
            "status": "global_holdout_not_viable_for_full_island_population",
        })
        relevant = native.loc[native.work_species.isin(donor_names)]
        by_entity = relevant.groupby("entity_ID").work_species.agg(set).to_dict()
        fold_rows = geographic_fold_coverage(
            sample, donor_names, by_entity, source_lookup, block_to_fold)
        results += [{"scope": scope, **row} for row in fold_rows]
        validations[scope] = {
            "source_scoped_profile_names": len(profile_maps[scope]),
            "matched_eligible_islands": len(selected),
            "private_public_name_counts_verified": True,
        }
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "n2_five_geographic_folds_donor_support.csv", results, FIELDS)
    _write_csv(output / "n2_all_eligible_global_holdout_loss.csv", unions, UNION_FIELDS)
    manifest = {
        "schema": "fd_n2_private_geographic_fold_source_support_v1",
        "status": "authorized_private_input_aggregate_only_coverage_audit",
        "folds": 5, "spatial_partition": "All eligible islands, spatial blocks balanced by island counts",
        "donor_holdout": "all complete taxon-name carriers in the heldout spatial fold removed from EVERY GIFT native source and trait donor",
        "not_a_n2_effect_test": True,
        "validation": validations,
        "rights": "only support counts and hashes released; no species, genus, source, island or block identifiers",
        "input_sha256": {
            "private_all": _sha(private_all), "private_direct": _sha(private_direct),
            "flora": _sha(flora), "geography": _sha(geography),
            "mainland_native": _sha(mainland_native),
            "source_assignments": _sha(assignments), "profile_reference": _sha(reference),
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"cohorts": unions, "folds": results, "status": manifest["status"]}))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("private-all", "private-direct", "flora", "g0",
                 "geography", "mainland-native", "assignments",
                 "reference", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    run(a.private_all, a.private_direct, a.flora, a.g0, a.geography,
        a.mainland_native, a.assignments, a.reference, a.output)


if __name__ == "__main__":
    main()
