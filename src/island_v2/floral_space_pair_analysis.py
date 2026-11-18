"""Matched between-island profile sharing under joint regional genus permutations.

Membership-only reservoir selection fixes identical All/Direct pairs within each
encoding. Shared names are removed by the kernel; complete regional permutation
pools and the saved-metric validation gates remain those of floral_space_null.
"""
from __future__ import annotations

import argparse
from collections import Counter
import importlib.metadata
import itertools
import json
from pathlib import Path
import platform
import subprocess

import numpy as np
import pandas as pd

from island_v2 import floral_space_pairs
from island_v2.floral_space_null import (
    REGIMES,
    SCOPES,
    RegionPermutationDesign,
    _file_record,
    bh_adjust_fixed,
    build_membership,
    prepare_region,
)
from island_v2.floral_space_pairs import prepare_pair_design, run_pair_comparison

ENCODINGS = (("full_state_sets", "all", "direct"),
             ("single_state", "all_single_state", "direct_single_state"))
DRAW_METRICS = ("profile_jaccard", "profile_intersection", "profile_union",
                "n_profiles_a", "n_profiles_b")
PAIR_COLUMNS = ["encoding", "evidence_scope", "analysis_regime", "island_a", "island_b",
                "n_taxa_a", "n_taxa_b", "n_shared_taxa", "n_unique_a", "n_unique_b",
                "n_profiles_a", "n_profiles_b", "profile_intersection", "profile_union",
                "profile_jaccard", "null_mean_jaccard", "null_sd_jaccard"]
DRAW_COLUMNS = ["encoding", "evidence_scope", "analysis_regime", "permutation",
                *[f"median_{metric}" for metric in DRAW_METRICS]]
INFERENCE_COLUMNS = [
    "observed_median_jaccard", "p_lower", "p_upper", "p_two_sided", "null_median",
    "null_q025", "null_q975", "null_sd", "observed_minus_null_median",
    "observed_to_null_median", "null_min", "null_max", "null_n_distinct_sampled_values",
    *[f"{kind}_median_{metric}" for metric in DRAW_METRICS[1:]
      for kind in ("observed", "null")],
]


def _positive_integer(value, label: str) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{label} must be a positive integer")


def _members(design: RegionPermutationDesign) -> dict[str, set[int]]:
    """Return incidence only: profile values never enter cohort selection."""
    members = {island: set() for island in design.island_ids}
    for island, taxon in zip(design.member_island_indices, design.member_taxon_indices, strict=True):
        members[design.island_ids[island]].add(int(taxon))
    return members


def select_matched_pairs(
    all_design: RegionPermutationDesign, direct_design: RegionPermutationDesign, *,
    min_unique: int = 20, max_pairs: int = 2000, seed=20261106,
) -> tuple[pd.DataFrame, dict]:
    """Uniform reservoir over sorted jointly supported unordered island pairs.

    All and Direct can have nonnested incidence and separate taxon-ID mappings.
    The eligibility intersection precedes both scope-specific support filters.
    Memory is bounded by incidence plus the retained reservoir and endpoint counts.
    """
    _positive_integer(min_unique, "Minimum nonshared support")
    _positive_integer(max_pairs, "Maximum pair count")
    if all_design.analysis_regime != direct_design.analysis_regime:
        raise ValueError("Matched designs must refer to the same regime")
    all_members, direct_members = _members(all_design), _members(direct_design)
    islands = sorted(set(all_design.island_ids) & set(direct_design.island_ids))
    rng = np.random.default_rng(seed)
    reservoir: list[tuple[str, str]] = []
    n_all = n_direct = n_joint = n_fail_both = 0
    supported_degrees: Counter = Counter()
    for a, b in itertools.combinations(islands, 2):
        all_shared = len(all_members[a] & all_members[b])
        direct_shared = len(direct_members[a] & direct_members[b])
        all_ok = min(len(all_members[a]), len(all_members[b])) - all_shared >= min_unique
        direct_ok = min(len(direct_members[a]), len(direct_members[b])) - direct_shared >= min_unique
        n_all += int(all_ok)
        n_direct += int(direct_ok)
        n_fail_both += int(not all_ok and not direct_ok)
        if not (all_ok and direct_ok):
            continue
        n_joint += 1
        supported_degrees.update((a, b))
        if len(reservoir) < max_pairs:
            reservoir.append((a, b))
        else:
            index = int(rng.integers(n_joint))
            if index < max_pairs:
                reservoir[index] = (a, b)
    pairs = pd.DataFrame(sorted(reservoir), columns=["island_a", "island_b"])
    sampled_degrees = Counter(endpoint for pair in reservoir for endpoint in pair)
    candidates = len(islands) * (len(islands) - 1) // 2
    audit = {
        "analysis_regime": all_design.analysis_regime,
        "min_unique": min_unique, "max_pairs": max_pairs,
        "n_eligible_islands_all": len(all_design.island_ids),
        "n_eligible_islands_direct": len(direct_design.island_ids),
        "n_matched_eligible_islands": len(islands),
        "n_candidate_pairs": candidates,
        "n_supported_pairs_all": n_all, "n_supported_pairs_direct": n_direct,
        "n_excluded_support_all": candidates - n_all,
        "n_excluded_support_direct": candidates - n_direct,
        "n_excluded_support_both": n_fail_both,
        "n_joint_supported_pairs": n_joint, "n_excluded_pairs": candidates - n_joint,
        "n_islands_with_supported_pairs": len(supported_degrees),
        "n_sampled_pairs": len(reservoir), "n_sampled_islands": len(sampled_degrees),
        "n_supported_pairs_not_sampled": n_joint - len(reservoir),
        "sampled_partners_min": min(sampled_degrees.values(), default=0),
        "sampled_partners_median": float(np.median(list(sampled_degrees.values()))) if sampled_degrees else 0.0,
        "sampled_partners_max": max(sampled_degrees.values(), default=0),
    }
    return pairs, audit


def _empty_summary(design: RegionPermutationDesign) -> dict:
    """Retain an empty hypothesis slot with explicit, undefined inference."""
    return {
        "analysis_regime": design.analysis_regime, "status": "no_eligible_pairs",
        "n_pairs": 0, "n_islands": 0, "n_permutations": 0,
        "n_taxa_pool": len(design.taxon_names), "n_genera": len(design.genus_indices),
        "n_permutable_genera": len(design.permutable_genus_indices),
        **dict.fromkeys(INFERENCE_COLUMNS, np.nan),
        "null_median_is_zero": None, "ratio_status": "no_eligible_pairs",
    }


def _parent_commit() -> str:
    process = subprocess.run(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[2],
                             capture_output=True, text=True, check=True)
    return process.stdout.strip()


def run_suite(
    results: Path, flora_path: Path, geography_path: Path, output: Path, *,
    permutations: int = 999, min_species: int = 50, min_unique: int = 20,
    max_pairs: int = 2000, seed: int = 20261106,
) -> dict:
    """Validate sources, run all 16 fixed slots, and save aggregate-only outputs."""
    parameters = dict(permutations=permutations, min_species=min_species,
                      min_unique=min_unique, max_pairs=max_pairs, seed=seed)
    for label, value in parameters.items():
        if label != "seed":
            _positive_integer(value, label)
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    results, flora_path, geography_path, output = map(
        lambda path: Path(path).resolve(), (results, flora_path, geography_path, output))
    flora = pd.read_csv(flora_path, usecols=["island_id", "accepted_species"], dtype=str)
    geography = pd.read_csv(geography_path, usecols=["island_id", "analysis_regime"], dtype=str)
    sources = {"flora": _file_record(flora_path), "geography": _file_record(geography_path)}
    source_provenance = {}
    # Hash provenance documents without copying arbitrary source content (which
    # could include restricted species-to-profile records) into public outputs.
    provenance_paths = [results / name for name in (
        "provenance_manifest.json", "permutation_manifest.json", "data_quality_audit.json")]
    provenance_paths.extend(results / scope / name for scope in SCOPES
                            for name in ("input_manifest.json", "summary.json"))
    for path in provenance_paths:
        if path.is_file():
            source_provenance[str(path.relative_to(results))] = _file_record(path)
    designs = {}
    for scope in SCOPES:
        profile_path = results / scope / "species_profiles.csv.gz"
        metric_path = results / scope / "island_metrics.csv"
        profiles = pd.read_csv(profile_path, usecols=["accepted_species", "genus", "profile_id"],
                               dtype={"accepted_species": str, "genus": str})
        metrics = pd.read_csv(metric_path, dtype={"island_id": str, "analysis_regime": str})
        sources[scope] = {"profiles": _file_record(profile_path), "island_metrics": _file_record(metric_path)}
        membership = build_membership(profiles, flora, geography)
        for regime in REGIMES:
            designs[scope, regime] = prepare_region(membership, metrics, regime, min_species=min_species)
    cohorts, audits, selection_streams = {}, [], []
    for encoding_index, (encoding, all_scope, direct_scope) in enumerate(ENCODINGS):
        for regime_index, regime in enumerate(REGIMES):
            entropy = [seed, 17, encoding_index, regime_index]
            pairs, audit = select_matched_pairs(
                designs[all_scope, regime], designs[direct_scope, regime],
                min_unique=min_unique, max_pairs=max_pairs, seed=np.random.SeedSequence(entropy))
            cohorts[encoding, regime] = pairs
            audits.append({"encoding": encoding, "all_scope": all_scope, "direct_scope": direct_scope,
                           "min_species": min_species, **audit})
            selection_streams.append({"encoding": encoding, "analysis_regime": regime, "entropy": entropy})
    summaries, draw_tables, pair_tables, permutation_streams = [], [], [], []
    scope_encodings = {scope: encoding for encoding, a, d in ENCODINGS for scope in (a, d)}
    for scope_index, scope in enumerate(SCOPES):
        encoding = scope_encodings[scope]
        for regime_index, regime in enumerate(REGIMES):
            entropy = [seed, 29, scope_index, regime_index]
            permutation_streams.append({"evidence_scope": scope, "analysis_regime": regime, "entropy": entropy})
            design, pairs = designs[scope, regime], cohorts[encoding, regime]
            if pairs.empty:
                summary = _empty_summary(design)
            else:
                pair_design = prepare_pair_design(design, pairs, min_unique=min_unique)
                summary, draws, pair_rows = run_pair_comparison(
                    pair_design, permutations=permutations, seed=np.random.SeedSequence(entropy))
                metadata = dict(encoding=encoding, evidence_scope=scope, analysis_regime=regime)
                draw_tables.append(draws.assign(**metadata))
                pair_tables.append(pair_rows.assign(**metadata))
            summaries.append({"encoding": encoding, "evidence_scope": scope, **summary,
                              "n_eligible_islands_scope": len(design.island_ids),
                              "maximum_expectation_difference": design.maximum_expectation_difference})
    summary_table = pd.DataFrame(summaries)
    assert len(summary_table) == len(SCOPES) * len(REGIMES) == 16
    summary_table["q_bh_16"] = bh_adjust_fixed(summary_table.p_two_sided)
    summary_table["q_by_16"] = np.minimum(1, summary_table.q_bh_16 * np.sum(1 / np.arange(1, 17)))
    draw_table = (pd.concat(draw_tables, ignore_index=True).reindex(columns=DRAW_COLUMNS)
                  if draw_tables else pd.DataFrame(columns=DRAW_COLUMNS))
    pair_table = (pd.concat(pair_tables, ignore_index=True).reindex(columns=PAIR_COLUMNS)
                  if pair_tables else pd.DataFrame(columns=PAIR_COLUMNS))
    audit_table = pd.DataFrame(audits)
    parent_commit = _parent_commit()
    output.mkdir(parents=True, exist_ok=True)
    summary_table.to_csv(output / "permutation_summary.csv", index=False)
    audit_table.to_csv(output / "pair_cohort_audit.csv", index=False)
    for name, table in (("permutation_null_draws.csv.gz", draw_table),
                        ("pairwise_observed_and_null.csv.gz", pair_table)):
        table.to_csv(output / name, index=False, compression={"method": "gzip", "mtime": 0})
    manifest = {
        "schema_version": "floral_profile_matched_pair_joint_genus_permutation_v1",
        "status": "complete", "design_status": "exploratory conditional assemblage comparison",
        "parent_commit": parent_commit,
        "fixed_comparison_family_size": 16, "evidence_scopes": list(SCOPES),
        "analysis_regimes": list(REGIMES),
        "encodings": [{"encoding": enc, "all_scope": a, "direct_scope": d} for enc, a, d in ENCODINGS],
        "parameters": parameters,
        "cli_arguments": {"results": str(results), "flora": str(flora_path),
                          "geography": str(geography_path), "output": str(output), **parameters},
        "software": {"python": platform.python_version(),
                     "libraries": {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy")}},
        "seed_streams": {"pair_selection": selection_streams, "permutation": permutation_streams},
        "rng": "NumPy default_rng with explicit SeedSequence entropy; no rejected identities or duplicate assignments",
        "pair_selection": "uniform reservoir without replacement over sorted jointly supported unordered pairs; membership only, no profile overlap",
        "statistic": "equal-pair-weighted median Jaccard of profile sets after removing shared taxon names before profile lookup",
        "conditioning": ["fixed observed island-by-taxon membership",
                         "full complete-taxon regional pool, including ineligible and unselected islands and pair-shared taxa",
                         "within-genus regional taxon profile multisets; one joint taxon assignment reused across all pairs"],
        "regional_interpretation": "separate conditional references, not reconstructed globally coherent species traits",
        "tail_definition": "inclusive finite-corrected tails (1+tail_count)/(B+1); two-sided=min(1,2*min(p_lower,p_upper))",
        "tie_tolerance": "100 * machine_epsilon * max(1, abs(observed)); numerical ties included in both tails",
        "multiple_testing": "BH fixed family of 16; BY=BH*H_16 capped at 1; missing comparisons enter as p=1 and retain NaN q",
        "null_range": "2.5 and 97.5 percentiles of permutation statistics, not an effect confidence interval",
        "pair_null_sd": "population standard deviation (ddof=0)",
        "claim_limits": ["conditional descriptive sharing or differentiation, not evolutionary convergence or causal selection",
                         "source-dependent observational coverage and profile assignments remain part of the inputs",
                         "taxon range-size associations can contribute; frequencies fixed per taxon rather than island occurrence",
                         "overlapping scopes and encodings are not independent replication",
                         "pair weighting is not equal-island weighting; support filters alter the target population"],
        "actual_counts": {"summary_slots": len(summary_table),
                          "nonempty_comparisons": int(summary_table.status.eq("ok").sum()),
                          "null_draw_rows": len(draw_table), "pair_rows": len(pair_table),
                          "cohort_audit_rows": len(audit_table),
                          "candidate_pairs": int(audit_table.n_candidate_pairs.sum()),
                          "joint_supported_pairs": int(audit_table.n_joint_supported_pairs.sum()),
                          "sampled_pairs_across_encoding_regimes": int(audit_table.n_sampled_pairs.sum())},
        "sources": sources, "source_result_provenance": source_provenance,
        "code": {"floral_space_pairs": _file_record(Path(floral_space_pairs.__file__).resolve()),
                 "floral_space_pair_analysis": _file_record(Path(__file__).resolve())},
        "outputs": {name: _file_record(output / name) for name in (
            "pair_cohort_audit.csv", "pairwise_observed_and_null.csv.gz",
            "permutation_summary.csv", "permutation_null_draws.csv.gz")},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("results", "flora", "geography", "output"):
        parser.add_argument("--" + argument, type=Path, required=True)
    for argument, default in (("permutations", 999), ("min-species", 50), ("min-unique", 20),
                              ("max-pairs", 2000), ("seed", 20261106)):
        parser.add_argument("--" + argument, type=int, default=default)
    args = parser.parse_args()
    run_suite(args.results, args.flora, args.geography, args.output,
              permutations=args.permutations, min_species=args.min_species,
              min_unique=args.min_unique, max_pairs=args.max_pairs, seed=args.seed)


if __name__ == "__main__":
    main()
