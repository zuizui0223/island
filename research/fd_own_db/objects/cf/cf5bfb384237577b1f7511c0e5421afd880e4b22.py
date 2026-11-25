"""Private-input-only genus/source-conditional floral distance reference (N2).

Hold out ALL evaluated island taxon names from every native GIFT source entity
before donor training. Condition post hoc on recorded island genus/carrier counts
and geography-specific putative mainland donor floral profiles. N2 is a
model-based fixed-distinct-profile support benchmark, NOT an exact
genus-preserving permutation, pollination niche, dispersal or causal null.

NEVER commit or emit taxon names, island IDs, raw floras, GIFT source identities,
per-island metric rows or matched-pair membership. Public CI tests only
synthetic helper behavior. Real inference requires licensed local inputs.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_carriers, load_profiles
from .distance_ses import (
    ENDPOINTS, MODES, SUMMARY_FIELDS, _score, _seed, _summary,
    distance_statistics, exact_sign_flip, regional_leave_one_out_pool,
    regional_profile_counts, sample_null,
)
from .geography_independent import _build_edges, independent_pair_matching
from .g1_fixed_margins import _bh
from .genus_source_null import load_private_maps, validate_island_counts, _sha
from .heldout_mainland_prediction import DonorModel, VALID_MODES, DEFAULT_SHRINK
from .pipeline import REF_INPUT, _write_csv
from .repertoire import match_origin_islands
from .spatial_clusters import connected_block_components

SCOPES = ("all", "direct")
NULLS = ("N1b_regional_islandfreq_LOO", "N2_global_genus", "N2_local_genus")
DONOR_VARIANTS = {"N2_global_genus": 2, "N2_local_genus": 3}
PAIR_COHORTS = ("original_block_components", "geography_both_verified")
PAIR_COLUMNS = (
    "cohort", "scope", "null", "mode", "endpoint", "status",
    "n_available_pairs", "n_independent_sign_groups",
    "observed_mean_SES_oceanic_minus_continental",
    "median_pair_SES_oceanic_minus_continental",
    "exact_two_sided_sign_p", "sign_null_q025", "sign_null_q975",
    "n_sign_configurations", "q_bh_posthoc4", "q_by_posthoc4",
)
CALIBRATION_FIELDS = (
    "scope", "geological_origin", "n_islands", "status",
    "mean_global_prior_name_log_loss", "mean_global_genus_name_log_loss",
    "mean_local_genus_name_log_loss",
    "mean_local_vs_global_genus_logloss_improvement",
    "fraction_islands_local_genus_better_than_global_genus",
    "mean_global_prior_trait_distance", "mean_global_genus_trait_distance",
    "mean_local_genus_trait_distance",
    "mean_local_vs_global_genus_trait_distance_improvement",
)

AUDIT_COLUMNS = (
    "scope", "geological_origin", "n_islands", "median_names",
    "median_observed_distinct_profiles", "mean_global_genus_donor_fraction",
    "mean_local_genus_donor_fraction", "median_N2_global_expected_unique",
    "median_N2_local_expected_unique", "median_local_expected_to_observed_ratio",
    "fraction_local_genus_support_below_025",
    "median_generative_richness", "median_generative_equal30_valid_fraction",
    "fraction_observed_richness_outside_generative_95pct", "status",
)


def donor_training(native: pd.DataFrame, profiles: pd.DataFrame,
                   heldout_names: set[str]) -> pd.DataFrame:
    """Remove heldout names before merging with private trait donors.

    A species in 20 source entities is still globally removed from all 20.
    """
    if not {"entity_ID", "work_species"}.issubset(native):
        raise ValueError("GIFT source entity/name columns missing")
    if not {"accepted_species", "genus", "profile_id"}.issubset(profiles):
        raise ValueError("Private flower profile map fields missing")
    if not heldout_names:
        raise ValueError("No evaluation names for global donor holdout")
    eligible = native.loc[~native.work_species.isin(heldout_names),
                          ["entity_ID", "work_species"]].drop_duplicates()
    if set(eligible.work_species) & heldout_names:
        raise AssertionError("Heldout target taxon still in mainland source feature")
    donor = eligible.merge(
        profiles[["accepted_species", "genus", "profile_id"]],
        left_on="work_species", right_on="accepted_species",
        how="inner", validate="many_to_one",
    ).drop(columns=["work_species"])
    if donor.empty or set(donor.accepted_species) & heldout_names:
        raise ValueError("Empty/leaky complete native GIFT floral donors")
    if donor.duplicated(["entity_ID", "accepted_species"]).any():
        raise AssertionError("Source-by-name donor duplicated after join")
    return donor


def expected_occupancy_weights(
    genus_counts: dict[str, int], model: DonorModel,
    island: str, variant: str,
) -> tuple[np.ndarray, float]:
    """Approximate probability of EACH category being present among carriers.

    For name carriers with independent donor profile probabilities p_gp,
    q_p = 1 - product_g(1-p_gp)^count_g. Conditional fixed-n category
    sampling using weights q_p is a **working reference**, not the exact
    conditional distribution of genus-level trait assignments given n
    different profiles: correlations between occupancy events are ignored.
    """
    if variant not in DONOR_VARIANTS or not genus_counts:
        raise ValueError("Unknown donor model or empty observed genus counts")
    log_absent = np.zeros(model.n_profiles, dtype=float)
    for genus, carriers in genus_counts.items():
        if not isinstance(carriers, (int, np.integer)) or carriers < 1:
            raise ValueError("Genus carrier count must be positive integer")
        prob = np.asarray(model.probabilities(island, genus)[DONOR_VARIANTS[variant]])
        if prob.shape != (model.n_profiles,) or not np.isfinite(prob).all():
            raise ValueError("Invalid donor probability dimension")
        if np.any(prob < 0) or not np.isclose(prob.sum(), 1., atol=1e-9):
            raise ValueError("Donor profile distribution not normalized")
        log_absent += carriers * np.log1p(-np.clip(prob, 0, 1-1e-12))
    inclusion = -np.expm1(log_absent)
    expected_unique = float(inclusion.sum())
    if not np.isfinite(expected_unique) or expected_unique <= 0:
        raise ValueError("Donor model predicts no finite profile occupancy")
    return inclusion / inclusion.sum(), expected_unique



def simulate_genus_traits(
    genus_counts: dict[str, int], model: DonorModel, island: str,
    distance: np.ndarray, *, variant: str, draws: int = 199,
    equal_n: int = 30, seed: int = 20261109,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate EACH held-out island name's trait from the donor genus model.

    Unlike the fixed-n inclusion-weight sample, the number of distinct
    profiles is generated rather than forced to the observed number.
    If >=equal_n types emerge, rarefy to equal_n. Output is a model adequacy
    diagnostic; never treat it as an additional independent origin test.
    """
    if variant not in DONOR_VARIANTS or draws < 20 or equal_n < 2:
        raise ValueError("Unsupported donor generator or simulation support")
    if not genus_counts or distance.shape != (model.n_profiles, model.n_profiles):
        raise ValueError("Empty model or incorrect categorical D dimension")
    parameters = []
    for genus, count in sorted(genus_counts.items()):
        if not isinstance(count, (int, np.integer)) or count < 1:
            raise ValueError("Invalid number of recorded genus-name carriers")
        probs = np.asarray(
            model.probabilities(island, genus)[DONOR_VARIANTS[variant]],
            dtype=float,
        )
        if probs.shape != (model.n_profiles,) or not np.isfinite(probs).all():
            raise ValueError("Invalid genus-specific model probabilities")
        if np.any(probs < 0) or not np.isclose(probs.sum(), 1., atol=1e-9):
            raise ValueError("Unnormalized source genus profile probability")
        parameters.append((count, probs))
    rng = np.random.default_rng(seed)
    richness = np.empty(draws, dtype=int)
    equal_scores = []
    for i in range(draws):
        realized = np.unique(np.concatenate([
            rng.choice(model.n_profiles, size=count, replace=True, p=p)
            for count, p in parameters
        ]))
        richness[i] = len(realized)
        if len(realized) >= equal_n:
            fixed = rng.choice(realized, size=equal_n, replace=False)
            equal_scores.append(distance_statistics(distance, fixed))
    return richness, np.asarray(equal_scores, dtype=float).reshape(-1, len(ENDPOINTS))


def source_coverage(model: DonorModel, island: str,
                    genus_counts: dict[str, int]) -> tuple[float, float]:
    """Fraction of recorded name carriers with global/local donor genus support."""
    n = sum(genus_counts.values())
    global_supported = local_supported = 0
    for genus, count in genus_counts.items():
        glob, local = model.genus_support(island, genus)
        global_supported += int(glob) * count
        local_supported += int(local) * count
    return global_supported/n, local_supported/n


def _validate_private_against_public(
    membership: pd.DataFrame, expected_g0: pd.DataFrame,
    carriers: dict[str, dict[int, int]], scope: str,
) -> None:
    validate_island_counts(membership, expected_g0, scope)
    for island, records in membership.groupby("island_id", sort=False):
        published = carriers.get(island)
        if published is None:
            raise ValueError("Private island missing from public carrier evidence")
        actual = {int(k): int(v) for k, v in
                  records.profile_id.value_counts().items()}
        if actual != published:
            raise ValueError(f"{scope}: private profile frequencies do not reproduce pinned public aggregate")


def heldout_name_calibration(part: pd.DataFrame, model: DonorModel,
                             island: str, space: ProfileSpace) -> dict:
    """Score unseen test name profiles; NEVER choose shrinkage using these scores.

    Per name: exact categorical negative log probability and expected original
    colour/form Jaccard distance. Evaluation names were removed from all
    source entities BEFORE DonorModel training. These are observational
    diagnostic losses rather than an independent evolutionary replication.
    """
    if part.empty or part.accepted_species.duplicated().any():
        raise ValueError("Invalid private heldout evaluation island")
    labels = ("global_prior", "global_genus", "local_genus")
    loss = np.zeros((len(labels), 2), dtype=float)
    for genus, profile in part[["genus", "profile_id"]].itertuples(index=False, name=None):
        pid = space.lookup[int(profile)]
        p = model.probabilities(island, str(genus))[[0, 2, 3]]
        if p.shape != (len(labels), len(space.ids)) or not np.isfinite(p).all():
            raise ValueError("Malformed donor predictive probabilities")
        loss[:, 0] -= np.log(np.clip(p[:, pid], 1e-15, 1.))
        loss[:, 1] += p @ space.matrix[:, pid]
    loss /= len(part)
    return {f"{metric}_{label}": float(loss[j, k])
            for j, label in enumerate(labels)
            for k, metric in enumerate(("log_loss", "trait_distance"))}


def _calibration_summaries(audits: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for row in audits:
        groups[(row["scope"], row["geological_origin"])].append(row)
    output = []
    for (scope, origin), samples in sorted(groups.items()):
        entry = {
            "scope": scope, "geological_origin": origin,
            "n_islands": len(samples),
            "status": "heldout_taxon_name_prediction_diagnostic"
                      if len(samples) >= 5 else "small_cell_suppressed",
        }
        if len(samples) >= 5:
            for label in ("global_prior", "global_genus", "local_genus"):
                for metric in ("name_log_loss", "trait_distance"):
                    field = ("log_loss_" + label if metric == "name_log_loss"
                             else "trait_distance_" + label)
                    entry[f"mean_{label}_{metric}"] = float(
                        np.mean([sample[field] for sample in samples]))
            improvements = np.array([
                x["log_loss_global_genus"] - x["log_loss_local_genus"]
                for x in samples
            ])
            distance_gain = np.array([
                x["trait_distance_global_genus"] -
                x["trait_distance_local_genus"] for x in samples
            ])
            entry["mean_local_vs_global_genus_logloss_improvement"] = float(
                np.mean(improvements))
            entry["fraction_islands_local_genus_better_than_global_genus"] = float(
                np.mean(improvements > 0))
            entry["mean_local_vs_global_genus_trait_distance_improvement"] = float(
                np.mean(distance_gain))
        output.append(entry)
    return output


def paired_primary(scores: list[dict], cohorts: dict[str, list[dict]],
                   geo: dict[str, dict]) -> list[dict]:
    """Four predeclared post-hoc tests: local-genus equal30 nearest-D, 2×2.

    Scope overlap + geographic limitations prohibit independent confirmations.
    Other nulls, full support, and pairwise D remain purely descriptive.
    """
    look = {(r["scope"], r["island_id"]): r["ses"] for r in scores
            if r["null"] == "N2_local_genus" and
            r["mode"] == "equal30" and
            r["endpoint"] == "mean_nearest_profile_distance"}
    out = []
    for cohort in PAIR_COHORTS:
        selected = cohorts[cohort]
        for scope in SCOPES:
            eligible = [p for p in selected if
                        (scope, p["continental_id"]) in look and
                        (scope, p["oceanic_id"]) in look]
            r = {
                "cohort": cohort, "scope": scope, "null": "N2_local_genus",
                "mode": "equal30", "endpoint": "mean_nearest_profile_distance",
                "n_available_pairs": len(eligible),
            }
            if len(eligible) < 8:
                r["status"] = "insufficient_complete_matched_pairs"
            else:
                diff = np.asarray([
                    look[(scope, p["oceanic_id"])] -
                    look[(scope, p["continental_id"])]
                    for p in eligible
                ])
                comps = (connected_block_components(eligible, geo)
                         if cohort == PAIR_COHORTS[0]
                         else [[i] for i in range(len(eligible))])
                r.update(exact_sign_flip(diff, comps))
                r["n_independent_sign_groups"] = len(comps)
                r["status"] = ("conservative_spatial_component_stress_test"
                               if cohort == PAIR_COHORTS[0] else
                               "posthoc_source_archipelago_disjoint_sensitivity")
            out.append(r)
    if len(out) != 4:
        raise AssertionError("Exactly four primary N2 scalar slots required")
    q = _bh([r.get("exact_two_sided_sign_p", 1.) for r in out])
    for r, adj in zip(out, q):
        if "exact_two_sided_sign_p" in r:
            r["q_bh_posthoc4"] = float(adj)
            r["q_by_posthoc4"] = float(min(1., adj * (1+1/2+1/3+1/4)))
    return out


def _coverage_summaries(audits: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for x in audits:
        groups[x["scope"], x["geological_origin"]].append(x)
    output = []
    for (scope, origin), rows in sorted(groups.items()):
        r = {"scope": scope, "geological_origin": origin,
             "n_islands": len(rows), "status": "observational_private_donor_support"}
        if len(rows) < 5:
            r["status"] = "small_cell_suppressed"
        else:
            for src, dst, fun in (
                ("n_names", "median_names", np.median),
                ("n_profiles", "median_observed_distinct_profiles", np.median),
                ("global_fraction", "mean_global_genus_donor_fraction", np.mean),
                ("local_fraction", "mean_local_genus_donor_fraction", np.mean),
                ("global_expected", "median_N2_global_expected_unique", np.median),
                ("local_expected", "median_N2_local_expected_unique", np.median),
                ("expected_ratio", "median_local_expected_to_observed_ratio", np.median),
                ("generated_richness_median", "median_generative_richness", np.median),
                ("generated_equal30_valid_fraction", "median_generative_equal30_valid_fraction", np.median),
            ):
                r[dst] = float(fun([x[src] for x in rows]))
            r["fraction_local_genus_support_below_025"] = float(
                np.mean([x["local_fraction"] < 0.25 for x in rows]))
            r["fraction_observed_richness_outside_generative_95pct"] = float(
                np.mean([x["generated_obs_outside95"] for x in rows]))
        output.append(r)
    return output


def run(private_all: Path, private_direct: Path, flora: Path,
        source_assignments: Path, mainland_native: Path,
        crosswalk: Path, g0_root: Path, reference: Path, input_root: Path,
        output: Path, *, source_mode: str = "geo_k5",
        shrink: float = DEFAULT_SHRINK,
        draws: int = 999, thinnings: int = 99,
        generative_draws: int = 199,
        seed: int = 20261109, strict: bool = True) -> dict:
    if source_mode not in VALID_MODES or not np.isfinite(shrink) or shrink <= 0:
        raise ValueError("Unrecognized GIFT source mode/shrinkage")
    if draws < 20 or thinnings < 20 or generative_draws < 20:
        raise ValueError("Need >=20 support reference and observed thinning draws")
    private_all, private_direct, flora, source_assignments, mainland_native = map(
        Path, (private_all, private_direct, flora, source_assignments, mainland_native))
    crosswalk, g0_root, reference, input_root, output = map(
        Path, (crosswalk, g0_root, reference, input_root, output))
    space = ProfileSpace(load_profiles(reference))
    if len(space.ids) != 496 or space.ids != list(range(496)):
        raise ValueError("Frozen profile reference 0..495 changed")
    g0 = {scope: pd.read_csv(g0_root/f"{scope}_islands.csv",
                             dtype={"island_id": str}) for scope in SCOPES}
    geo_frame = pd.read_csv(crosswalk, dtype={"island_id": str})
    if geo_frame.island_id.duplicated().any():
        raise ValueError("Duplicate geographic identity")
    geo = geo_frame.set_index("island_id").to_dict("index")
    original, _ = match_origin_islands(
        g0["all"].to_dict("records"), g0["direct"].to_dict("records"), geo)
    edges, _ = _build_edges(
        g0["all"].to_dict("records"), g0["direct"].to_dict("records"), geo)
    strict_pairs, verified = independent_pair_matching(edges, "both_verified")
    if strict and (len(original) != 37 or len(strict_pairs) != 10):
        raise ValueError("Frozen geographic cohort counts changed")
    cohorts = dict(zip(PAIR_COHORTS, (original, strict_pairs)))
    focal_ids = {p[side] for pairs in cohorts.values()
                 for p in pairs for side in ("continental_id", "oceanic_id")}
    reference_maps = {
        scope: load_private_maps(path, scope, set(space.ids), strict=strict)
        for scope, path in (("all", private_all), ("direct", private_direct))
    }
    all_flora = pd.read_csv(
        flora, usecols=["island_id", "accepted_species", "origin_status"],
        dtype=str)
    if all_flora.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Duplicate island/taxon membership in private checklist")
    all_flora = all_flora.loc[all_flora.island_id.isin(focal_ids)].copy()
    assignments = pd.read_csv(
        source_assignments,
        usecols=["island_id", "source_mode", "entity_ID"], dtype=str)
    assignments = assignments.loc[
        assignments.island_id.isin(focal_ids) &
        assignments.source_mode.eq(source_mode),
        ["island_id", "entity_ID"]].drop_duplicates()
    if assignments.island_id.nunique() != len(focal_ids):
        raise ValueError("At least one tested island has no authorized GIFT source assignment")
    mainland = pd.read_csv(
        mainland_native, usecols=["entity_ID", "work_species"], dtype=str)
    if mainland.empty:
        raise ValueError("Empty native mainland GIFT source")
    rows, audits, source_hashes = [], [], {}
    carrier_files = {}
    for scope in SCOPES:
        carrier_path = input_root / "floral_space_20261106" / scope / "island_profile_carriers.csv.gz"
        carrier_files[scope] = carrier_path
        carriers = load_carriers(carrier_path)
        group = g0[scope].loc[g0[scope].island_id.isin(focal_ids)]
        sample = all_flora.merge(
            reference_maps[scope], on="accepted_species", how="inner",
            validate="many_to_one")
        sample = sample.loc[sample.island_id.isin(group.island_id)]
        _validate_private_against_public(sample, group, carriers, scope)
        if sample.island_id.nunique() != len(focal_ids):
            raise ValueError("One matched island missing in private complete-profile sample")
        heldout = set(sample.accepted_species)
        train = donor_training(mainland, reference_maps[scope], heldout)
        model = DonorModel(train, assignments, len(space.ids), shrink=shrink)
        source_hashes[scope] = {
            "n_global_target_names_heldout": len(heldout),
            "n_complete_native_donor_names_after_holdout": model.n_donor_names,
            "n_assigned_sources": model.n_source_entities,
        }
        regional = regional_profile_counts(g0[scope].to_dict("records"), carriers)
        g0_lookup = group.set_index("island_id").to_dict("index")
        for island, part in sample.groupby("island_id", sort=True):
            row = g0_lookup[island]
            origin = geo[island].get("geological_origin")
            if origin not in ("continental", "oceanic"):
                raise ValueError("Target origin outside geological cohort")
            profile_support = np.asarray(
                sorted({space.lookup[int(pid)] for pid in part.profile_id}),
                dtype=int)
            n = len(profile_support)
            if n < 2:
                continue
            genus_counts = {str(g): int(v) for g, v in
                            part.genus.value_counts().items()}
            gl, local = source_coverage(model, island, genus_counts)
            weights = {}
            expected = {}
            for variant in DONOR_VARIANTS:
                weights[variant], expected[variant] = expected_occupancy_weights(
                    genus_counts, model, island, variant)
            generated_richness, generated_equal30 = simulate_genus_traits(
                genus_counts, model, island, space.matrix,
                variant="N2_local_genus", draws=generative_draws,
                equal_n=30, seed=_seed(seed, scope, island, "N2g"))
            richness_q025, richness_q975 = np.quantile(
                generated_richness, [.025, .975])
            audits.append({
                "scope": scope, "geological_origin": origin,
                "n_names": len(part), "n_profiles": n,
                "global_fraction": gl, "local_fraction": local,
                "global_expected": expected["N2_global_genus"],
                "local_expected": expected["N2_local_genus"],
                "expected_ratio": expected["N2_local_genus"]/n,
                "generated_richness_median": float(np.median(generated_richness)),
                "generated_equal30_valid_fraction": float(
                    len(generated_equal30)/generative_draws),
                "generated_obs_outside95": bool(
                    n < richness_q025 or n > richness_q975),
                **heldout_name_calibration(part, model, island, space),
            })
            rng = np.random.default_rng(_seed(seed, scope, island, "N2_thinning"))
            thin = (np.asarray([
                distance_statistics(space.matrix, rng.choice(profile_support, 30, replace=False))
                for _ in range(thinnings)
            ]) if n >= 30 else None)
            regional_pool = regional_leave_one_out_pool(
                regional, row["analysis_regime"], set(carriers[island]), space.lookup)
            regional_weights = np.asarray([
                regional[row["analysis_regime"]][space.ids[int(idx)]]
                - int(space.ids[int(idx)] in carriers[island])
                for idx in regional_pool], dtype=float)
            for mode in MODES:
                if mode == "equal30" and thin is None:
                    continue
                nn = n if mode == "full" else 30
                observed = (
                    distance_statistics(space.matrix, profile_support)
                    if mode == "full" else thin.mean(axis=0))
                for variant in NULLS:
                    if variant == "N1b_regional_islandfreq_LOO":
                        simulated = sample_null(
                            space.matrix, regional_pool, nn, draws=draws,
                            seed=_seed(seed, "N2compareN1b", scope, island, mode),
                            weights=regional_weights,
                        )
                    else:
                        simulated = sample_null(
                            space.matrix, np.arange(len(space.ids)), nn,
                            draws=draws,
                            seed=_seed(seed, variant, scope, island, mode),
                            weights=weights[variant])
                    if simulated is None:
                        continue
                    scored = _score(
                        observed, simulated,
                        thinning=thin if mode == "equal30" else None)
                    for item in scored:
                        if item is not None:
                            rows.append({
                                "scope": scope, "island_id": island,
                                "null": variant, "mode": mode,
                                "endpoint": item["endpoint"],
                                "geological_origin": origin,
                                "region": row["analysis_regime"],
                                "n_profiles": n, **item,
                            })
                if mode == "equal30" and len(generated_equal30) >= 20:
                    for item in _score(observed, generated_equal30, thinning=thin):
                        if item is not None:
                            rows.append({
                                "scope": scope, "island_id": island,
                                "null": "N2g_local_genus_taxon_simulation",
                                "mode": mode, "endpoint": item["endpoint"],
                                "geological_origin": origin,
                                "region": row["analysis_regime"],
                                "n_profiles": n, **item,
                            })
    summary = _summary(rows)
    pair_rows = paired_primary(rows, cohorts, geo)
    coverage = _coverage_summaries(audits)
    calibration = _calibration_summaries(audits)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output/"n2_distance_ses_summary.csv", summary, SUMMARY_FIELDS)
    _write_csv(output/"n2_source_donor_coverage.csv", coverage, AUDIT_COLUMNS)
    _write_csv(output/"n2_heldout_prediction_calibration.csv", calibration, CALIBRATION_FIELDS)
    _write_csv(output/"n2_paired_spatial_sensitivity.csv", pair_rows, PAIR_COLUMNS)
    manifest = {
        "schema": "fd_private_heldout_genus_source_distance_SES_v1",
        "status": "private_input_only_exploratory_model_reference",
        "source_mode": source_mode, "shrinkage": shrink,
        "n_profiles_reference": len(space.ids),
        "n_unique_focal_islands": len(focal_ids),
        "cohorts": {
            "original_pairs": len(original),
            "original_spatial_components": len(connected_block_components(original, geo)),
            "source_block_disjoint_verified_pairs": len(strict_pairs),
            "source_unknown_archipelago_endpoints": verified["n_unknown_archipelago_endpoints"],
        },
        "draws": draws, "thinnings": thinnings,
        "generative_draws": generative_draws, "seed": seed,
        "generative_model_check": (
            "Independently sample one profile per observed genus name, then "
            "evaluate simulated number of distinct types and equally thin "
            "only simulations with >=30 distinct types. No observed-n "
            "conditioning; no extra geology tests. A low valid fraction "
            "or observed richness outside 95% simulation range warns "
            "against using fixed-n N2 as an ecological null."),
        "calibration_gate": "heldout individual-name categorical log loss and expected original trait distance, averaged across island units, for global prior/global genus/local genus; positive local-vs-global gains are necessary but not sufficient",
        "calibration_output": "n2_heldout_prediction_calibration.csv; grouped means only, no test name/island/source identities",
        "reference_definition": (
            "Model 1-probability-absent-over-taxon-carriers from GIFT native donors "
            "conditioned on OBSERVED genera, then sample n distinct profiles "
            "without replacement weighted by their expected inclusion. "
            "Not an exact genus-fixed permutation or fully generative null."),
        "global_holdout": "All tested name units excluded from EVERY native GIFT region before donor construction, per evidence scope",
        "donor_models": {
            "N2_global_genus": "global GIFT native donors in recorded genus with shrinkage to global trait prior",
            "N2_local_genus": "candidate five-nearest-mainland GIFT donors in recorded genus, shrunk to global genus posterior",
            "N1b": "focal-island-excluded same-regime island profile-frequency reference",
        },
        "critical_limitations": [
            "Profile inclusion indicators are dependent; weighted fixed-n draws are an approximation, not an exact conditional genus randomization",
            "Genus composition is an observed oracle, not predicted colonization or true mainland ancestry",
            "Every genus without independent source donor falls back to a smoothed global/genus prior: audit coverage before inference",
            "Native status refers to GIFT mainland donors, not necessarily explicitly native status of focal island flora",
            "Geography, shared genera, island-level observation/detection biases and phylogeny are not randomized",
            "No unrestricted trait files or local island/genus/source identities in the outputs",
        ],
        "donor_support_audit": source_hashes,
        "input_sha256": {
            "private_all": _sha(private_all),
            "private_direct": _sha(private_direct),
            "flora": _sha(flora),
            "mainland_native": _sha(mainland_native),
            "source_assignments": _sha(source_assignments),
            "geology_crosswalk": _sha(crosswalk),
            "profile_reference": _sha(reference),
            **{f"{scope}_carrier_aggregate": _sha(path)
               for scope, path in carrier_files.items()},
        },
    }
    (output/"manifest.json").write_text(json.dumps(
        manifest, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({
        "status": manifest["status"],
        "focal_islands": len(focal_ids),
        "n_aggregated_summary": len(summary),
        "n_paired_slots": len(pair_rows),
        "donor_support": source_hashes,
        "n_min_p": min((r["exact_two_sided_sign_p"] for r in pair_rows
                         if "exact_two_sided_sign_p" in r), default=None),
        "interpretation": "No source-independent selection inference",
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("private-all", "private-direct", "flora", "source-assignments",
                "mainland-native", "crosswalk", "g0", "reference",
                "input", "output"):
        parser.add_argument("--"+arg, type=Path, required=True)
    parser.add_argument("--source-mode", choices=VALID_MODES, default="geo_k5")
    parser.add_argument("--shrink", type=float, default=DEFAULT_SHRINK)
    parser.add_argument("--draws", type=int, default=999)
    parser.add_argument("--thinnings", type=int, default=99)
    parser.add_argument("--generative-draws", type=int, default=199)
    parser.add_argument("--allow-synthetic", action="store_true")
    args = parser.parse_args()
    run(args.private_all, args.private_direct, args.flora,
        args.source_assignments, args.mainland_native,
        args.crosswalk, args.g0, args.reference, args.input, args.output,
        source_mode=args.source_mode, shrink=args.shrink,
        draws=args.draws, thinnings=args.thinnings,
        generative_draws=args.generative_draws,
        strict=not args.allow_synthetic)


if __name__ == "__main__":
    main()
