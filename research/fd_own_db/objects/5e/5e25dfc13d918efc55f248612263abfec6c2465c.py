"""Geography-blocked held-out taxon-name flower-profile prediction (private input).

This evaluates the *predictive adequacy* of global genus and geographically
local native-mainland genus donors. It does NOT test geological selection,
causal filtering, phylogenetically independent islands or actual colonization.
All individual island, taxon, genus, GIFT-source and spatial-block identifiers
remain in memory. Only non-identifying group summaries and hashes are written.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_profiles, load_carriers
from .genus_source_null import load_private_maps, validate_island_counts, _sha
from .n2_source_distance import (
    donor_training, heldout_name_calibration, source_coverage,
    _validate_private_against_public,
)
from .n2_source_coverage_audit import balanced_block_folds
from .heldout_mainland_prediction import DonorModel, DEFAULT_SHRINK
from .pipeline import _write_csv
from .proper_scoring import (
    inspect_kernel, heldout_probability_scores, heldout_name_source_gains,
)
from .trait_kernels import KERNEL_NAMES, kernel_quality, heldout_axis_source_gains

SCOPES = ("all", "direct")
AXIS_METRICS = tuple(
    f"{prefix}_{name}" for prefix in ("genus", "source")
    for name in KERNEL_NAMES
)
AXIS_FIELDS = (
    "scope", "region", "n_islands", "n_spatial_blocks", "status",
    *tuple(field for metric in AXIS_METRICS
           for field in (f"mean_{metric}", f"{metric}_block_q025",
                         f"{metric}_block_q975")),
)
SOURCE_NAME_METRICS = (
    "source_gain_logloss", "source_gain_trait_distance_loss",
    "source_gain_brier", "source_gain_kernel",
)
NAME_BALANCE_FIELDS = (
    "scope", "n_unique_names", "n_unique_genera", "n_name_fold_units",
    "mean_local_genus_donor_support",
    *tuple(v for k in SOURCE_NAME_METRICS
           for v in (f"mean_{k}", f"{k}_genus_boot_q025",
                     f"{k}_genus_boot_q975")),
)
METRICS = (
    "genus_gain", "source_gain", "source_distance_gain",
    "genus_brier_gain", "source_brier_gain",
    "genus_kernel_gain", "source_kernel_gain",
)
FIELDS = (
    "scope", "region", "n_islands", "n_spatial_blocks", "status",
    "mean_global_genus_coverage", "mean_local_source_genus_coverage",
    "mean_prior_logloss", "mean_global_genus_logloss",
    "mean_local_source_genus_logloss", "fraction_islands_local_better",
    "mean_genus_gain", "genus_gain_block_q025", "genus_gain_block_q975",
    "mean_source_gain", "source_gain_block_q025", "source_gain_block_q975",
    "mean_source_distance_gain", "source_distance_gain_block_q025",
    "source_distance_gain_block_q975",
    "mean_genus_brier_gain", "genus_brier_gain_block_q025", "genus_brier_gain_block_q975",
    "mean_source_brier_gain", "source_brier_gain_block_q025", "source_brier_gain_block_q975",
    "mean_genus_kernel_gain", "genus_kernel_gain_block_q025", "genus_kernel_gain_block_q975",
    "mean_source_kernel_gain", "source_kernel_gain_block_q025", "source_kernel_gain_block_q975",
)


def _seed(seed: int, scope: str, region: str) -> int:
    payload = f"{seed}:{scope}:{region}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def block_bootstrap_mean(
    data: pd.DataFrame, *, response: str, clusters: str = "spatial_block",
    draws: int = 999, seed: int = 1,
) -> tuple[float, float, float]:
    """Equal-island mean, whole-geographic-block bootstrap confidence limits."""
    if draws < 99 or len(data) < 5 or response not in data or clusters not in data:
        raise ValueError("Insufficient block bootstrap support")
    if data[clusters].isna().any() or not np.isfinite(data[response]).all():
        raise ValueError("Missing/invalid geographic block or response")
    samples = data.groupby(clusters, sort=True).agg(
        size=(response, "size"), response_sum=(response, "sum"))
    n = samples["size"].to_numpy(dtype=int)
    total = samples["response_sum"].to_numpy(dtype=float)
    if len(n) < 5:
        raise ValueError("At least five separate geographic blocks required")
    rng = np.random.default_rng(seed)
    indices = rng.integers(len(n), size=(draws, len(n)))
    estimates = total[indices].sum(axis=1) / n[indices].sum(axis=1)
    lo, hi = np.quantile(estimates, [.025, .975])
    return float(data[response].mean()), float(lo), float(hi)


def summarize(rows: list[dict], *, draws: int = 999,
              seed: int = 20261109, min_cell: int = 5) -> list[dict]:
    """Only grouped scores; never serialize internal target rows."""
    data = pd.DataFrame(rows)
    groups = []
    for scope, frame in data.groupby("scope", sort=True):
        groups.append((scope, "all", frame))
        groups += [(scope, reg, part) for reg, part in
                   frame.groupby("analysis_regime", sort=True)]
    output = []
    for scope, region, frame in groups:
        row = {"scope": scope, "region": region,
               "n_islands": len(frame),
               "n_spatial_blocks": frame.spatial_block.nunique(),
               "status": "geography_block_heldout_exploratory_prediction"}
        if len(frame) < min_cell or frame.spatial_block.nunique() < 5:
            row["status"] = "insufficient_independent_geographic_support"
        else:
            for field, label in (
                ("global_coverage", "mean_global_genus_coverage"),
                ("local_coverage", "mean_local_source_genus_coverage"),
                ("prior_logloss", "mean_prior_logloss"),
                ("global_logloss", "mean_global_genus_logloss"),
                ("local_logloss", "mean_local_source_genus_logloss"),
            ):
                row[label] = float(frame[field].mean())
            row["fraction_islands_local_better"] = float(
                (frame.source_gain > 0).mean())
            for metric in METRICS:
                mu, low, high = block_bootstrap_mean(
                    frame, response=metric, draws=draws,
                    seed=_seed(seed, scope, region + metric))
                row[f"mean_{metric}"] = mu
                row[f"{metric}_block_q025"] = low
                row[f"{metric}_block_q975"] = high
        output.append(row)
    return output


def taxon_balanced_source_summary(
    name_fold: list[pd.DataFrame], *, draws: int = 999,
    seed: int = 20261110,
) -> list[dict]:
    """One taxon name = one unit; bootstrap entire recorded genera.

    Requires upstream aggregation of repeated island observations WITHIN
    each geographic fold, then equal weighting of distinct fold predictions
    within each taxon. Genus blocks are NOT phylogenetically independent
    families; this is a strong sensitivity, not a causal bootstrap.
    """
    if draws < 99 or not name_fold:
        raise ValueError("Missing private name-fold prediction evidence")
    data = pd.concat(name_fold, ignore_index=True)
    keys = ["scope", "fold", "accepted_species"]
    if data.duplicated(keys).any():
        raise ValueError("Name-fold outcomes were not deduplicated")
    results = []
    for scope, part in data.groupby("scope", sort=True):
        taxon = part.groupby(["accepted_species", "genus"], sort=False)[
            [*SOURCE_NAME_METRICS, "has_local_genus_donor"]
        ].mean().reset_index()
        row = {
            "scope": scope, "n_unique_names": len(taxon),
            "n_unique_genera": taxon.genus.nunique(),
            "n_name_fold_units": len(part),
            "mean_local_genus_donor_support": float(
                taxon.has_local_genus_donor.mean()),
        }
        for metric in SOURCE_NAME_METRICS:
            group = taxon.groupby("genus", sort=True)[metric].agg(
                ["sum", "size"]).to_numpy(dtype=float)
            if len(group) < 20:
                raise ValueError("Insufficient taxonomically distinct genus clusters")
            rng = np.random.default_rng(_seed(seed, scope, metric))
            estimates = []
            for chunk in range(0, draws, 199):
                nrep = min(199, draws-chunk)
                idx = rng.integers(len(group), size=(nrep, len(group)))
                estimates.append(
                    group[idx, 0].sum(axis=1) /
                    group[idx, 1].sum(axis=1))
            values = np.concatenate(estimates)
            row[f"mean_{metric}"] = float(taxon[metric].mean())
            row[f"{metric}_genus_boot_q025"], row[
                f"{metric}_genus_boot_q975"] = map(
                    float, np.quantile(values, [.025, .975]))
        results.append(row)
    return results


def summarize_trait_axes(rows: list[dict], *, draws: int = 999,
                         seed: int = 20261110) -> list[dict]:
    """One original island vote for each post-hoc trait kernel.

    Retains the same geographically held-out island population across all
    kernels. No source IDs, species, genera or individual island outcomes.
    """
    sample = pd.DataFrame(rows)
    out = []
    for scope, by_scope in sample.groupby("scope", sort=True):
        groups = [(scope, "all", by_scope)]
        groups += [(scope, region, sub) for region, sub
                   in by_scope.groupby("analysis_regime", sort=True)]
        for _, region, frame in groups:
            row = {"scope": scope, "region": region,
                   "n_islands": len(frame),
                   "n_spatial_blocks": frame.spatial_block.nunique(),
                   "status": "posthoc_heldout_trait_kernel_decomposition"}
            if len(frame) < 5 or frame.spatial_block.nunique() < 5:
                row["status"] = "insufficient_independent_geographic_support"
            else:
                for metric in AXIS_METRICS:
                    avg, lo, hi = block_bootstrap_mean(
                        frame, response=metric, draws=draws,
                        seed=_seed(seed, scope, region + metric))
                    row["mean_" + metric] = avg
                    row[metric + "_block_q025"] = lo
                    row[metric + "_block_q975"] = hi
                if not np.isclose(
                        row["mean_source_additive"],
                        .5*(row["mean_source_colour"] + row["mean_source_form"]),
                        atol=1e-9):
                    raise AssertionError("Additive marginal kernel identity broken")
            out.append(row)
    return out


def run(profile_all: Path, profile_direct: Path, flora: Path,
        mainland_native: Path, assignments: Path, g0_root: Path,
        geographic_blocks: Path, reference: Path, input_root: Path,
        output: Path, *, source_mode: str = "geo_k5",
        shrink: float = DEFAULT_SHRINK, seed: int = 20261109,
        bootstrap_draws: int = 999) -> dict:
    if source_mode != "geo_k5" or not np.isfinite(shrink) or shrink <= 0:
        raise ValueError("Frozen geo_k5 source assignment and positive shrink required")
    paths = list(map(Path, (profile_all, profile_direct, flora,
                            mainland_native, assignments, geographic_blocks,
                            reference, input_root)))
    profile_all, profile_direct, flora, mainland_native, assignments, geographic_blocks, reference, input_root = paths
    g0_root, output = Path(g0_root), Path(output)
    space = ProfileSpace(load_profiles(reference))
    if len(space.ids) != 496 or space.ids != list(range(496)):
        raise ValueError("Frozen 496 categorical profiles and IDs changed")
    kernel_audit = inspect_kernel(space.matrix)
    similarity = kernel_audit.pop("kernel")
    trait_k, trait_audit = kernel_quality(load_profiles(reference))
    if not np.allclose(trait_k["additive"], similarity, atol=1e-12):
        raise AssertionError("Original D and additive colour/form kernel mismatch")
    profiles = {
        "all": load_private_maps(profile_all, "all", set(space.ids)),
        "direct": load_private_maps(profile_direct, "direct", set(space.ids)),
    }
    g0 = {
        scope: pd.read_csv(g0_root / f"{scope}_islands.csv", dtype={"island_id": str})
        for scope in SCOPES
    }
    geography = pd.read_csv(geographic_blocks,
        usecols=["island_id", "spatial_block"], dtype=str)
    if geography.island_id.duplicated().any():
        raise ValueError("Duplicate geographic block island")
    all_eligible = g0["all"].merge(geography, on="island_id",
                                    validate="one_to_one")
    if len(all_eligible) != len(g0["all"]):
        raise ValueError("Incomplete All geography")
    folds = balanced_block_folds(all_eligible)
    native = pd.read_csv(mainland_native,
        usecols=["entity_ID", "work_species"], dtype=str).drop_duplicates()
    source = pd.read_csv(assignments,
        usecols=["island_id", "source_mode", "entity_ID"], dtype=str)
    source = source.loc[source.source_mode.eq(source_mode),
                        ["island_id", "entity_ID"]].drop_duplicates()
    flora_table = pd.read_csv(flora,
        usecols=["island_id", "accepted_species"], dtype=str).drop_duplicates()
    if flora_table.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Duplicated name-bearing island records")
    internal = []
    axis_islands = []
    private_name_fold = []
    support = {}
    for scope in SCOPES:
        eligible = g0[scope].merge(geography, on="island_id",
                                    validate="one_to_one")
        eligible["fold"] = eligible.spatial_block.map(folds)
        if eligible.fold.isna().any():
            raise ValueError("Missing geographic fold assignment")
        member = flora_table.loc[flora_table.island_id.isin(eligible.island_id)].merge(
            profiles[scope], on="accepted_species", validate="many_to_one")
        validate_island_counts(member, g0[scope], scope)
        carriers = load_carriers(input_root / "floral_space_20261106" /
                                scope / "island_profile_carriers.csv.gz")
        _validate_private_against_public(member, g0[scope], carriers, scope)
        member = member.merge(
            eligible[["island_id", "analysis_regime", "spatial_block", "fold"]],
            on="island_id", validate="many_to_one")
        source_scoped = source.loc[source.island_id.isin(eligible.island_id)]
        if source_scoped.island_id.nunique() != len(eligible):
            raise ValueError("Missing independent source entity assignment")
        support[scope] = {"n_islands": len(eligible), "scored": 0}
        for fold, batch in member.groupby("fold", sort=True):
            named_parts = []
            heldout = set(batch.accepted_species)
            donor = donor_training(native, profiles[scope], heldout)
            selected_assign = source_scoped.loc[
                source_scoped.island_id.isin(batch.island_id)]
            model = DonorModel(donor, selected_assign, len(space.ids),
                               shrink=shrink)
            for island, part in batch.groupby("island_id", sort=True):
                audit = heldout_name_calibration(part, model, island, space)
                proper = heldout_probability_scores(part, model, island, similarity)
                axis_gains = heldout_axis_source_gains(
                    part, model, island, trait_k)
                if abs(axis_gains["source_additive"] - (
                        proper["global_genus_kernel"] -
                        proper["local_genus_kernel"])) > 1e-9:
                    raise AssertionError("Trait-axis score disagrees with source proper score")
                axis_islands.append({
                    "scope": scope,
                    "analysis_regime": part.analysis_regime.iloc[0],
                    "spatial_block": part.spatial_block.iloc[0],
                    **axis_gains,
                })
                names = heldout_name_source_gains(part, model, island, similarity)
                for name_metric, check_metric in (
                    ("source_gain_logloss", "logloss"),
                    ("source_gain_trait_distance_loss", "trait_distance_loss"),
                    ("source_gain_brier", "brier"),
                    ("source_gain_kernel", "kernel"),
                ):
                    check = (
                        proper["global_genus_" + check_metric] -
                        proper["local_genus_" + check_metric])
                    if abs(float(names[name_metric].mean())-check) > 1e-9:
                        raise AssertionError("Per-name/whole-island scoring disagreement")
                named_parts.append(names)
                if (abs(proper["global_genus_logloss"] - audit["log_loss_global_genus"]) > 1e-9 or
                        abs(proper["local_genus_logloss"] - audit["log_loss_local_genus"]) > 1e-9):
                    raise AssertionError("Two implementations of donor heldout log loss disagree")
                frac_g, frac_l = source_coverage(
                    model, island,
                    {str(g): int(n) for g, n in part.genus.value_counts().items()})
                row = {
                    "scope": scope, "analysis_regime": part.analysis_regime.iloc[0],
                    "spatial_block": part.spatial_block.iloc[0],
                    "global_coverage": frac_g, "local_coverage": frac_l,
                    "prior_logloss": audit["log_loss_global_prior"],
                    "global_logloss": audit["log_loss_global_genus"],
                    "local_logloss": audit["log_loss_local_genus"],
                    "genus_gain": (
                        audit["log_loss_global_prior"] -
                        audit["log_loss_global_genus"]),
                    "source_gain": (
                        audit["log_loss_global_genus"] -
                        audit["log_loss_local_genus"]),
                    "source_distance_gain": (
                        audit["trait_distance_global_genus"] -
                        audit["trait_distance_local_genus"]),
                    "genus_brier_gain": (
                        proper["global_prior_brier"] - proper["global_genus_brier"]),
                    "source_brier_gain": (
                        proper["global_genus_brier"] - proper["local_genus_brier"]),
                    "genus_kernel_gain": (
                        proper["global_prior_kernel"] - proper["global_genus_kernel"]),
                    "source_kernel_gain": (
                        proper["global_genus_kernel"] - proper["local_genus_kernel"]),
                }
                internal.append(row)
                support[scope]["scored"] += 1
            # Keep only one private row per (taxon name, heldout fold), despite
            # the same taxon occurring on many focal islands.
            fold_name = pd.concat(named_parts, ignore_index=True).groupby(
                ["accepted_species", "genus"], sort=False
            )[[*SOURCE_NAME_METRICS, "has_local_genus_donor"]].mean().reset_index()
            fold_name["scope"] = scope
            fold_name["fold"] = int(fold)
            private_name_fold.append(fold_name)
    summarized = summarize(internal, draws=bootstrap_draws, seed=seed)
    axis_summary = summarize_trait_axes(
        axis_islands, draws=bootstrap_draws, seed=seed)
    names_summary = taxon_balanced_source_summary(
        private_name_fold, draws=bootstrap_draws, seed=seed)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "n2_geographic_crossfit_prediction.csv",
               summarized, FIELDS)
    _write_csv(output / "n2_trait_component_crossfit.csv",
               axis_summary, AXIS_FIELDS)
    _write_csv(output / "n2_taxon_balanced_source_scores.csv",
               names_summary, NAME_BALANCE_FIELDS)
    manifest = {
        "schema": "fd_n2_5_spatial_block_fold_heldout_prediction_v1",
        "status": "private_input_only_exploratory_predictive_calibration",
        "models": ["global prior", "known global genus", "known local source genus"],
        "frozen_source": source_mode, "shrinkage": shrink,
        "fold_assignment": "five disjoint source spatial blocks balanced by All island counts; Direct inherits All block assignment",
        "global_name_holdout": "all names in the entire heldout geographic fold removed from every mainland source entity before trait donor training",
        "proper_name_score": "categorical negative log probability of complete colour/form profile; equal-island grouped average",
        "proper_probability_scores": "Brier (strictly proper) and K=1-original-D kernel quadratic score (proper but not strictly proper at rank deficient K); equal-island scoring",
        "kernel_audit": kernel_audit,
        "trait_axis_kernels": trait_audit,
        "trait_axis_inference": "post-hoc: separate colour and form marginal similarity, additive identity, product joint association; product signal is not uniquely pure interaction",
        "trait_axis_output": "n2_trait_component_crossfit.csv, grouped island-equal geographic-block intervals only",
        "taxon_balanced_sensitivity": "within heldout fold average each accepted name over islands, then across folds equally by name; whole-genus bootstrap, not independent family phylogeny",
        "taxon_balanced_output": "n2_taxon_balanced_source_scores.csv; two aggregate scope rows, no names or genus keys",
        "distance_loss": "expected original mean state-set Jaccard distance; not a proper probabilistic calibration score",
        "bootstrap": "999 whole-spatial-block resamplings within each group; exploratory confidence limits; not source phylogenetic independence",
        "support": support,
        "not_estimated": "geological origin effect, dispersal, plant selection, causal trait convergence, pollinator fitness",
        "privacy": "no per-island, per-name, genus, source or spatial block identity written",
        "source_sha256": {
            "private_all": _sha(profile_all), "private_direct": _sha(profile_direct),
            "flora": _sha(flora), "native_mainland": _sha(mainland_native),
            "source_assignments": _sha(assignments),
            "geography_blocks": _sha(geographic_blocks),
            "public_reference": _sha(reference),
        },
    }
    (output / "manifest.json").write_text(json.dumps(
        manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"support": support, "n_aggregate_rows": len(summarized),
                      "status": manifest["status"]}))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("profile-all", "profile-direct", "flora",
                "mainland-native", "assignments", "g0",
                "geographic-blocks", "reference", "input", "output"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--shrink", type=float, default=DEFAULT_SHRINK)
    p.add_argument("--bootstrap", type=int, default=999)
    args = p.parse_args()
    run(args.profile_all, args.profile_direct, args.flora,
        args.mainland_native, args.assignments, args.g0,
        args.geographic_blocks, args.reference, args.input,
        args.output, shrink=args.shrink, bootstrap_draws=args.bootstrap)


if __name__ == "__main__":
    main()
