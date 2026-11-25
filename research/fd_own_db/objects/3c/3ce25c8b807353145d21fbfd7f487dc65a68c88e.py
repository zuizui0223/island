"""Geographic source-assignment placebo with genus/donor support matching.

PRIVATE INPUT ONLY: PR321 name-to-flower profiles, island floras, GIFT native
mainland source lists and assignments. ALL raw island, name, genus, fold-block,
and source region identifiers stay in memory. Publish grouped diagnostics only.

This is post-hoc, observational prediction robustness, not historical
colonization, causal selection or geological-origin inference.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_carriers, load_profiles
from .genus_source_null import load_private_maps, validate_island_counts, _sha
from .heldout_mainland_prediction import DonorModel, DEFAULT_SHRINK
from .n2_source_coverage_audit import balanced_block_folds
from .n2_source_distance import donor_training, _validate_private_against_public
from .n2_crossfit_macro import block_bootstrap_mean
from .pipeline import _write_csv
from .proper_scoring import heldout_probability_scores, inspect_kernel
from .trait_kernels import KERNEL_NAMES, kernel_quality, heldout_axis_source_gains
from .marginal_fixed_coupling import marginal_fixed_placebo_gain, profile_axis_codes

SCOPES = ("all", "direct")
POLICIES = ("genus_coverage_only", "genus_coverage_plus_donor_size_factor2")
MAX_ABS_COVERAGE_DIFF = 0.02
MAX_DONOR_SIZE_RATIO = 2.0
MIN_CONTROLS = 2
COUPLING_FIELDS = (
    "scope", "policy", "region", "n_islands", "n_spatial_blocks", "status",
    "mean_joint_before_marginal_fixed", "mean_joint_marginal_fixed",
    "joint_marginal_fixed_block_q025", "joint_marginal_fixed_block_q975",
    "mean_brier_marginal_fixed", "brier_marginal_fixed_block_q025",
    "brier_marginal_fixed_block_q975", "largest_ipf_marginal_error",
    "mean_ipf_iterations", "fraction_joint_marginal_fixed_positive",
)
TRAIT_PLACEBO_FIELDS = (
    "scope", "policy", "region", "n_islands", "n_spatial_blocks", "status",
    *tuple(field for trait in KERNEL_NAMES
            for field in (f"mean_{trait}_gain", f"{trait}_block_q025",
                          f"{trait}_block_q975")),
    "mean_form_minus_colour_gain",
    "form_minus_colour_block_q025",
    "form_minus_colour_block_q975",
)

REPORT_COLUMNS = (
    "scope", "policy", "region", "n_islands", "n_spatial_blocks",
    "status", "median_eligible_placebo_assignments",
    "mean_focal_minus_control_genus_coverage",
    "mean_control_to_focal_donor_name_count_ratio",
    "median_control_to_focal_donor_name_count_ratio",
    "mean_geographic_source_entity_jaccard", "fraction_actual_better_kernel",
    "mean_kernel_gain", "kernel_block_q025", "kernel_block_q975",
    "mean_brier_gain", "brier_block_q025", "brier_block_q975",
    "mean_coverage_gain", "coverage_block_q025", "coverage_block_q975",
    "leave_one_block_out_kernel_gain_min", "leave_one_block_out_kernel_gain_max",
)
FOLD_COLUMNS = (
    "scope", "policy", "fold", "n_islands",
    "mean_kernel_gain", "mean_brier_gain",
    "mean_coverage_difference", "mean_donor_name_ratio",
)


def matching_candidates(
    focal_id: str, focal_region: str, focal_block: str,
    focal_genus_counts: dict[str, int],
    catalog: dict[str, dict], policy: str,
) -> tuple[list[dict], int]:
    """Pick two closest coverage-matched *different* source assignments.

    Match on traits-free observed genus carrier counts and independent
    donor support, never observed floral profiles, Brier, or K-score.
    Controls are within the focal cross-fitting fold and broad region, but
    from a different geographic block; identical source entity sets are
    excluded to ensure the placebo actually changes source assignment.
    """
    if policy not in POLICIES or focal_id not in catalog:
        raise ValueError("Invalid focal source or control policy")
    if not focal_genus_counts or min(focal_genus_counts.values()) < 1:
        raise ValueError("Invalid genus-carrier counts")
    focal = catalog[focal_id]
    if (focal["region"], focal["spatial_block"]) != (focal_region, focal_block):
        raise ValueError("Focal geography disagrees with test membership")
    if not focal["source_entities"] or focal["n_donor_names"] < 1:
        raise ValueError("Focal source has no independent donor evidence")

    total = sum(focal_genus_counts.values())
    def coverage(source):
        donor_genera = source["genera"]
        return sum(count for genus, count in focal_genus_counts.items()
                   if genus in donor_genera) / total

    actual = coverage(focal)
    candidates = []
    for source_id, candidate in catalog.items():
        if (source_id == focal_id or candidate["region"] != focal_region
                or candidate["spatial_block"] == focal_block
                or candidate["source_entities"] == focal["source_entities"]
                or candidate["n_donor_names"] < 1):
            continue
        control = coverage(candidate)
        diff = actual - control
        if abs(diff) > MAX_ABS_COVERAGE_DIFF + 1e-10:
            continue
        ratio = candidate["n_donor_names"] / focal["n_donor_names"]
        if (policy == POLICIES[1] and
                not 1 / MAX_DONOR_SIZE_RATIO <= ratio <= MAX_DONOR_SIZE_RATIO):
            continue
        a, b = focal["source_entities"], candidate["source_entities"]
        candidates.append({
            "source_island": source_id, "support": control,
            "coverage_diff": diff, "donor_ratio": ratio,
            "entity_jaccard": len(a & b) / len(a | b),
        })
    # Sorting uses ONLY support covariates, with a stable anonymous tie
    # breaker; never select a more favorable floral outcome.
    selected = sorted(
        candidates, key=lambda c: (
            abs(c["coverage_diff"]), abs(np.log(c["donor_ratio"])),
            c["source_island"]))[:MIN_CONTROLS]
    return (selected if len(selected) == MIN_CONTROLS else []), len(candidates)


def summarize_placebo(records: list[dict], *,
                      draws: int = 4999, seed: int = 20261110
                      ) -> tuple[list[dict], list[dict]]:
    """Block clustered summaries; no exact source/island identity emitted."""
    if draws < 99 or not records:
        raise ValueError("No placebo evidence or too few uncertainty draws")
    sample = pd.DataFrame(records)
    grouped = []
    for (scope, policy), grp in sample.groupby(["scope", "policy"], sort=True):
        grouped.append((scope, policy, "all", grp))
        grouped += [(scope, policy, regime, sub)
                    for regime, sub in grp.groupby("analysis_regime", sort=True)]
    reports = []
    for scope, policy, regime, group in grouped:
        result = {
            "scope": scope, "policy": policy, "region": regime,
            "n_islands": len(group),
            "n_spatial_blocks": group.spatial_block.nunique(),
            "status": "support_matched_posthoc_prediction"
                      if len(group) >= 5 and group.spatial_block.nunique() >= 5
                      else "insufficient_geographic_support",
        }
        if result["status"] == "support_matched_posthoc_prediction":
            result.update({
                "median_eligible_placebo_assignments": float(
                    group.n_eligible_control_assignments.median()),
                "mean_focal_minus_control_genus_coverage": float(
                    group.coverage_difference.mean()),
                "mean_control_to_focal_donor_name_count_ratio": float(
                    group.control_source_donor_count_ratio.mean()),
                "median_control_to_focal_donor_name_count_ratio": float(
                    group.control_source_donor_count_ratio.median()),
                "mean_geographic_source_entity_jaccard": float(
                    group.mean_overlap_source_entity_jaccard.mean()),
                "fraction_actual_better_kernel": float(
                    (group.matched_kernel_gain > 0).mean()),
            })
            for col, dest in (
                ("matched_kernel_gain", "kernel"),
                ("matched_brier_gain", "brier"),
                ("coverage_difference", "coverage"),
            ):
                mean, low, high = block_bootstrap_mean(
                    group, response=col, draws=draws, seed=seed)
                result[f"mean_{dest}_gain"] = mean
                result[f"{dest}_block_q025"] = low
                result[f"{dest}_block_q975"] = high
            block = group.groupby("spatial_block").matched_kernel_gain.agg(
                ["sum", "size"])
            total = group.matched_kernel_gain.sum()
            remaining = ((total - block["sum"]) /
                         (len(group) - block["size"]))
            result["leave_one_block_out_kernel_gain_min"] = float(remaining.min())
            result["leave_one_block_out_kernel_gain_max"] = float(remaining.max())
        reports.append(result)

    folds = []
    for (scope, policy, fold), group in sample.groupby(
            ["scope", "policy", "fold"], sort=True):
        result = {"scope": scope, "policy": policy, "fold": int(fold),
                  "n_islands": len(group)}
        if len(group) >= 5:
            result.update({
                "mean_kernel_gain": float(group.matched_kernel_gain.mean()),
                "mean_brier_gain": float(group.matched_brier_gain.mean()),
                "mean_coverage_difference": float(
                    group.coverage_difference.mean()),
                "mean_donor_name_ratio": float(
                    group.control_source_donor_count_ratio.mean()),
            })
        folds.append(result)
    return reports, folds


def summarize_trait_placebo(records: list[dict], *, draws: int = 4999,
                            seed: int = 20261110) -> list[dict]:
    """Post-hoc marginal vs joint source-placebo scores on SAME matched islands."""
    if draws < 99 or not records:
        raise ValueError("No source-matched test evidence")
    data = pd.DataFrame(records)
    output = []
    for (scope, policy), group in data.groupby(["scope", "policy"], sort=True):
        collection = [(scope, policy, "all", group)]
        collection += [(scope, policy, regime, sub)
                       for regime, sub in group.groupby("analysis_regime", sort=True)]
        for _, _, regime, test in collection:
            record = {"scope": scope, "policy": policy, "region": regime,
                      "n_islands": len(test),
                      "n_spatial_blocks": test.spatial_block.nunique(),
                      "status": "posthoc_marginal_joint_source_placebo"}
            if len(test) < 5 or test.spatial_block.nunique() < 5:
                record["status"] = "insufficient_geographic_support"
            else:
                for name in KERNEL_NAMES:
                    avg, lo, hi = block_bootstrap_mean(
                        test, response="matched_" + name + "_gain",
                        draws=draws, seed=seed)
                    record["mean_" + name + "_gain"] = avg
                    record[name + "_block_q025"] = lo
                    record[name + "_block_q975"] = hi
                if not np.isclose(
                    record["mean_additive_gain"],
                    .5*(record["mean_colour_gain"] +
                        record["mean_form_gain"]), atol=1e-9):
                    raise AssertionError("Additive source placebo is not the mean of marginals")
                # A paired difference on identical islands, NOT a comparison
                # of two independent bootstrap intervals. Sign and magnitude
                # remain post-hoc and Jaccard scales differ by trait.
                matched = test.assign(
                    _form_minus_colour=test["matched_form_gain"] -
                    test["matched_colour_gain"])
                delta, lo, hi = block_bootstrap_mean(
                    matched, response="_form_minus_colour",
                    draws=draws, seed=seed)
                record["mean_form_minus_colour_gain"] = delta
                record["form_minus_colour_block_q025"] = lo
                record["form_minus_colour_block_q975"] = hi
            output.append(record)
    return output


def summarize_marginal_fixed(records: list[dict], *, draws: int = 4999,
                             seed: int = 20261110) -> list[dict]:
    """Aggregation of the actual-vs-placebo source gain with BOTH margins fixed."""
    data = pd.DataFrame(records)
    if data.empty:
        return []
    out = []
    for (scope, policy), frame in data.groupby(["scope", "policy"], sort=True):
        strata = [(scope, policy, "all", frame)]
        strata += [(scope, policy, region, sub) for region, sub
                   in frame.groupby("analysis_regime", sort=True)]
        for _, _, region, sample in strata:
            item = {"scope": scope, "policy": policy, "region": region,
                    "n_islands": len(sample),
                    "n_spatial_blocks": sample.spatial_block.nunique(),
                    "status": "posthoc_marginal_fixed_source_association_prediction"}
            if len(sample) < 5 or sample.spatial_block.nunique() < 5:
                item["status"] = "insufficient_independent_geographic_support"
            else:
                item["mean_joint_before_marginal_fixed"] = float(
                    sample.matched_joint_product_gain.mean())
                for col, name in (
                    ("coupling_fixed_joint_gain", "joint_marginal_fixed"),
                    ("coupling_fixed_brier_gain", "brier_marginal_fixed"),
                ):
                    mean, lo, hi = block_bootstrap_mean(
                        sample, response=col, draws=draws, seed=seed)
                    item["mean_" + name] = mean
                    item[name + "_block_q025"] = lo
                    item[name + "_block_q975"] = hi
                item["largest_ipf_marginal_error"] = float(
                    sample.max_ipf_error.max())
                item["mean_ipf_iterations"] = float(
                    sample.mean_ipf_iterations.mean())
                item["fraction_joint_marginal_fixed_positive"] = float(
                    (sample.coupling_fixed_joint_gain > 0).mean())
            out.append(item)
    return out


def run(private_all: Path, private_direct: Path, flora: Path,
        mainland_native: Path, assignments: Path, g0_root: Path,
        geography: Path, reference: Path, input_root: Path, output: Path,
        *, shrink: float = DEFAULT_SHRINK,
        bootstrap_draws: int = 4999,
        evaluate_marginal_fixed: bool = False) -> dict:
    """Source-grade five-fold controlled comparison; aggregate-only results."""
    if not np.isfinite(shrink) or shrink <= 0:
        raise ValueError("Require positive frozen source shrinkage")
    private_all, private_direct, flora, mainland_native, assignments = map(
        Path, (private_all, private_direct, flora, mainland_native, assignments))
    g0_root, geography, reference, input_root, output = map(
        Path, (g0_root, geography, reference, input_root, output))
    space = ProfileSpace(load_profiles(reference))
    if len(space.ids) != 496 or space.ids != list(range(496)):
        raise ValueError("Frozen colour/form profile reference changed")
    kernel_evidence = inspect_kernel(space.matrix)
    similarity = kernel_evidence.pop("kernel")
    reference_profiles = load_profiles(reference)
    trait_k, trait_audit = kernel_quality(reference_profiles)
    colour_codes, form_codes = profile_axis_codes(reference_profiles)
    if not np.allclose(trait_k["additive"], similarity, atol=1e-12):
        raise AssertionError("Original D additive kernel identity failed")
    geo = pd.read_csv(
        geography, usecols=["island_id", "spatial_block"], dtype=str)
    if geo.island_id.duplicated().any():
        raise ValueError("Duplicated frozen island blocks")
    g0 = {scope: pd.read_csv(
        g0_root / f"{scope}_islands.csv", dtype={"island_id": str})
          for scope in SCOPES}
    all_islands = g0["all"].merge(geo, on="island_id", validate="one_to_one")
    if len(all_islands) != len(g0["all"]):
        raise ValueError("Missing All source geographic blocks")
    fold_lookup = balanced_block_folds(all_islands)
    profile = {
        scope: load_private_maps(p, scope, set(space.ids))
        for scope, p in (("all", private_all), ("direct", private_direct))
    }
    flora_names = pd.read_csv(
        flora, usecols=["island_id", "accepted_species"], dtype=str)
    if flora_names.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Original island×name records repeated")
    native = pd.read_csv(
        mainland_native, usecols=["entity_ID", "work_species"],
        dtype=str).drop_duplicates()
    source = pd.read_csv(
        assignments, usecols=["island_id", "source_mode", "entity_ID"],
        dtype=str)
    source = source.loc[
        source.source_mode.eq("geo_k5"),
        ["island_id", "entity_ID"]].drop_duplicates()
    all_records = []
    fixed_records = []
    support = {}
    for scope in SCOPES:
        selected = g0[scope].merge(geo, on="island_id", validate="one_to_one")
        selected["fold"] = selected.spatial_block.map(fold_lookup)
        if selected.fold.isna().any():
            raise ValueError("Unassigned eligible geographic block")
        carriers = load_carriers(
            input_root / "floral_space_20261106" / scope /
            "island_profile_carriers.csv.gz")
        members = flora_names.loc[
            flora_names.island_id.isin(selected.island_id)].merge(
                profile[scope], on="accepted_species",
                how="inner", validate="many_to_one")
        validate_island_counts(members, g0[scope], scope)
        _validate_private_against_public(members, g0[scope], carriers, scope)
        members = members.merge(
            selected[["island_id", "analysis_regime", "spatial_block", "fold"]],
            on="island_id", validate="many_to_one")
        assigned = source.loc[source.island_id.isin(selected.island_id)]
        if assigned.island_id.nunique() != len(selected):
            raise ValueError("Potential source entity missing from eligible island")
        support[scope] = {"n_islands": len(selected), "n_scoreable": 0}
        for fold, batch in members.groupby("fold", sort=True):
            heldout = set(batch.accepted_species)
            donors = donor_training(native, profile[scope], heldout)
            fold_ids = set(batch.island_id)
            fold_source = assigned.loc[assigned.island_id.isin(fold_ids)]
            model = DonorModel(donors, fold_source, len(space.ids), shrink=shrink)
            entity_by_island = (fold_source.groupby("island_id").entity_ID
                                .agg(lambda x: frozenset(x)).to_dict())
            by_island = {}
            for island, group in batch.groupby("island_id", sort=True):
                local_total = model.local_total_counts.get(island)
                if local_total is None or local_total.sum() <= 0:
                    continue
                local_genera = {genus for key, genus
                               in model.local_genus_counts
                               if key == island}
                by_island[island] = {
                    "region": group.analysis_regime.iloc[0],
                    "spatial_block": group.spatial_block.iloc[0],
                    "genera": local_genera,
                    "n_donor_names": int(local_total.sum()),
                    "source_entities": entity_by_island[island],
                }
            for island, test in batch.groupby("island_id", sort=True):
                if island not in by_island:
                    continue
                focal_genus = {
                    str(genus): int(count)
                    for genus, count in test.genus.value_counts().items()}
                actual = heldout_probability_scores(
                    test, model, island, similarity)
                actual_axis = heldout_axis_source_gains(
                    test, model, island, trait_k)
                support[scope]["n_scoreable"] += 1
                for policy in POLICIES:
                    control, n_available = matching_candidates(
                        island, test.analysis_regime.iloc[0],
                        test.spatial_block.iloc[0], focal_genus,
                        by_island, policy)
                    if len(control) != MIN_CONTROLS:
                        continue
                    alternate = [heldout_probability_scores(
                        test, model, c["source_island"], similarity)
                        for c in control]
                    alt_axis = [heldout_axis_source_gains(
                        test, model, c["source_island"], trait_k)
                        for c in control]
                    actual_cov = sum(
                        count for genus, count in focal_genus.items()
                        if genus in by_island[island]["genera"]) / len(test)
                    alt_cover = np.mean([c["support"] for c in control])
                    all_records.append({
                        "scope": scope, "fold": int(fold),
                        "policy": policy,
                        "analysis_regime": test.analysis_regime.iloc[0],
                        "spatial_block": test.spatial_block.iloc[0],
                        "matched_kernel_gain": float(
                            np.mean([a["local_genus_kernel"] for a in alternate])
                            - actual["local_genus_kernel"]),
                        "matched_brier_gain": float(
                            np.mean([a["local_genus_brier"] for a in alternate])
                            - actual["local_genus_brier"]),
                        "coverage_difference": float(actual_cov - alt_cover),
                        "control_source_donor_count_ratio": float(np.mean([
                            c["donor_ratio"] for c in control])),
                        "mean_overlap_source_entity_jaccard": float(np.mean([
                            c["entity_jaccard"] for c in control])),
                        "n_eligible_control_assignments": n_available,
                        **{
                            "matched_" + name + "_gain": float(
                                actual_axis["source_" + name] -
                                np.mean([a["source_" + name] for a in alt_axis]))
                            for name in KERNEL_NAMES
                        },
                    })
                    if abs(all_records[-1]["matched_additive_gain"] -
                           all_records[-1]["matched_kernel_gain"]) > 1e-9:
                        raise AssertionError("Additive K gain must equal original source-placebo score")
                    if evaluate_marginal_fixed and policy == POLICIES[1]:
                        association = marginal_fixed_placebo_gain(
                            test, model, island,
                            [c["source_island"] for c in control],
                            reference_profiles,
                            colour_codes=colour_codes, form_codes=form_codes,
                            joint_kernel=trait_k["joint_product"])
                        fixed_records.append({
                            **all_records[-1],
                            "coupling_fixed_joint_gain": association[
                                "actual_minus_placebo_marginal_fixed_joint"],
                            "coupling_fixed_brier_gain": association[
                                "actual_minus_placebo_marginal_fixed_brier"],
                            "max_ipf_error": association["maximum_ipf_error"],
                            "mean_ipf_iterations": association[
                                "mean_ipf_iterations"],
                        })
    reports, folds = summarize_placebo(
        all_records, draws=bootstrap_draws)
    trait_scores = summarize_trait_placebo(
        all_records, draws=bootstrap_draws)
    fixed_scores = summarize_marginal_fixed(
        fixed_records, draws=bootstrap_draws)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "n2_source_placebo_matched_groups.csv",
               reports, REPORT_COLUMNS)
    _write_csv(output / "n2_source_placebo_matched_folds.csv",
               folds, FOLD_COLUMNS)
    _write_csv(output / "n2_source_placebo_trait_kernels.csv",
               trait_scores, TRAIT_PLACEBO_FIELDS)
    if evaluate_marginal_fixed:
        if not fixed_scores:
            raise ValueError("No eligible marginal-fixed support-matched source cohort")
        _write_csv(output / "n2_source_placebo_marginal_fixed.csv",
                   fixed_scores, COUPLING_FIELDS)
    manifest = {
        "schema": "fd_heldout_geo_source_placebo_support_match_v1",
        "status": "private_input_posthoc_coverage_matched_source_prediction",
        "frozen_profiles": len(space.ids),
        "folds": 5,
        "control": ("same fold and region, different geographic block and "
                    "different source entity set; choose two controls with "
                    "nearest observed genus donor coverage, |difference|<=0.02;"
                    " add factor-two total source donor name caliper"),
        "source_selected_without_floral_outcome": True,
        "source_prediction": ("focal genus identities known; proper "
                              "K-score with PSD K=1-original D and Brier"),
        "comparisons": list(POLICIES),
        "noncausal": ("Neither true ancestral sources, colonization, "
                      "phylogenetic independence nor pollination adaptation"),
        "missing_support": support,
        "kernel": kernel_evidence,
        "marginal_and_joint_trait_kernels": trait_audit,
        "marginal_fixed_coupling_computed": evaluate_marginal_fixed,
        "marginal_fixed_method": "IPF project local source and each matched alternate onto global genus colour-set/form-set marginals, preserve association odds; compare strictly proper joint and Brier loss on exact same recorded profile types",
        "marginal_fixed_output": ("n2_source_placebo_marginal_fixed.csv"
                                  if evaluate_marginal_fixed else None),
        "trait_kernel_output": "n2_source_placebo_trait_kernels.csv; post-hoc colour, form, additive, product-joint source placebo on same island/control groups",
        "interaction_limit": "joint-product prediction can detect paired colour/form associations but is not by itself a pure interaction effect",
        "privacy": ("Only group/fold summary and hashes persisted; "
                    "no island, source, genus or taxon-level outputs"),
        "source_sha256": {
            "private_all": _sha(private_all),
            "private_direct": _sha(private_direct),
            "flora": _sha(flora),
            "mainland_native": _sha(mainland_native),
            "assignments": _sha(assignments),
            "geography": _sha(geography),
            "profile_reference": _sha(reference),
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": manifest["status"],
        "n_aggregate_rows": len(reports),
        "n_fold_rows": len(folds),
        "eligible_island_counts": support,
        "strict_support": [r for r in reports
                           if r["region"] == "all"
                           and r["policy"] == POLICIES[1]],
    }))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for path in ("private-all", "private-direct", "flora", "mainland-native",
                 "assignments", "g0", "geography", "reference",
                 "input", "output"):
        p.add_argument("--" + path, type=Path, required=True)
    p.add_argument("--shrink", type=float, default=DEFAULT_SHRINK)
    p.add_argument("--bootstrap-draws", type=int, default=4999)
    p.add_argument("--evaluate-marginal-fixed", action="store_true",
                   help="Private data, strict support-matched cohort only; IPF to identical colour and form marginals")
    a = p.parse_args()
    run(a.private_all, a.private_direct, a.flora, a.mainland_native,
        a.assignments, a.g0, a.geography, a.reference,
        a.input, a.output, shrink=a.shrink,
        bootstrap_draws=a.bootstrap_draws,
        evaluate_marginal_fixed=a.evaluate_marginal_fixed)


if __name__ == "__main__":
    main()
