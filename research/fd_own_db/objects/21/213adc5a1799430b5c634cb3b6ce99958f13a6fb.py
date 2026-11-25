"""Taxon-name-held-out GIFT mainland prediction of island floral repertoires.

PRIVATE source inputs stay local to the authorized researcher. Only aggregate
model scores, model-comparison uncertainty, source support, and hashes are
written. No species names, island identifiers, source entities, exact pairs,
taxon-to-trait rows, or restricted source text are released.

This is conditional prediction GIVEN the observed island genus composition.
It is not a colonization model, an ancestral mainland reconstruction, an
independent species-level evolutionary test, or a conditional-genus permutation.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_profiles
from .genus_source_null import load_private_maps, validate_island_counts, _sha
from .repertoire import (
    REGIONS, _paired_swap_statistic, _upper_tail, match_origin_islands,
)
from .g1_fixed_margins import _bh

MODELS = ("global_mainland", "local_mainland", "global_genus", "local_genus")
COMPARISONS = (
    ("global_genus", "global_mainland"),
    ("local_genus", "global_genus"),
    ("local_mainland", "global_mainland"),
    ("local_genus", "local_mainland"),
)
LOSSES = ("log_loss", "trait_distance_loss",
          "repertoire_L2_error", "repertoire_kernel_error")
NULL_ENDPOINTS = ("exact_profile_L2", "categorical_similarity_kernel")
PROVENANCE_COMMIT = "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c"
VALID_MODES = ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10")
PSEUDO_GLOBAL_PER_PROFILE = 0.02
DEFAULT_SHRINK = 8.0
SCORE_COLUMNS = (
    "scope", "model", "endpoint", "n_islands", "mean_loss", "median_loss",
)
IMPROVEMENT_COLUMNS = (
    "scope", "better_model", "worse_model", "endpoint", "n_pairs",
    "mean_improvement", "pair_bootstrap_q025", "pair_bootstrap_q975",
)
COVERAGE_COLUMNS = (
    "scope", "origin", "n_islands", "median_complete_name_count",
    "mean_global_genus_supported_fraction", "mean_local_genus_supported_fraction",
)
RESIDUAL_COLUMNS = (
    "scope", "source_model", "response", "endpoint", "n_pairs",
    "observed_statistic", "swap_null_median", "swap_null_q025",
    "swap_null_q975", "p_swap_upper", "q_bh_primary4", "q_by_primary4",
)


def _write_csv(path: Path, rows: list[dict], columns: tuple[str, ...]) -> None:
    pd.DataFrame(rows).reindex(columns=columns).to_csv(path, index=False)


class DonorModel:
    """Four probabilistic models estimated only from non-target mainland taxa.

    global_mainland: no genus and no island-specific source entity
    local_mainland: source entity union, no genus
    global_genus: global mainland genus, shrink to global reference
    local_genus: island source genus, shrink to global genus posterior
    """

    def __init__(self, donors: pd.DataFrame, assignments: pd.DataFrame,
                 n_profiles: int, *, shrink: float = DEFAULT_SHRINK):
        if shrink <= 0 or not math.isfinite(shrink):
            raise ValueError("shrink must be a positive finite pseudo-count")
        required = {"entity_ID", "accepted_species", "genus", "profile_id"}
        if not required.issubset(donors.columns):
            raise ValueError("Missing donor columns")
        if not {"island_id", "entity_ID"}.issubset(assignments.columns):
            raise ValueError("Missing source-assignment fields")
        if donors.empty:
            raise ValueError("No training mainland traits")
        donors = donors.copy()
        donors.profile_id = donors.profile_id.astype(int)
        if donors.profile_id.min() < 0 or donors.profile_id.max() >= n_profiles:
            raise ValueError("Donor profile outside reference")
        self.n_profiles = n_profiles
        self.shrink = float(shrink)
        global_species = donors.drop_duplicates("accepted_species")
        counts = np.bincount(global_species.profile_id, minlength=n_profiles).astype(float)
        self.global_prior = (
            (counts + PSEUDO_GLOBAL_PER_PROFILE) /
            (counts.sum() + n_profiles * PSEUDO_GLOBAL_PER_PROFILE)
        )
        self.global_genus_counts = {}
        for genus, part in global_species.groupby("genus", sort=False):
            self.global_genus_counts[genus] = np.bincount(
                part.profile_id, minlength=n_profiles).astype(float)
        local = donors.merge(
            assignments.drop_duplicates(["island_id", "entity_ID"]),
            on="entity_ID", how="inner",
        ).drop_duplicates(["island_id", "accepted_species"])
        self.local_total_counts = {}
        self.local_genus_counts = {}
        for island, part in local.groupby("island_id", sort=False):
            self.local_total_counts[island] = np.bincount(
                part.profile_id, minlength=n_profiles).astype(float)
            for genus, members in part.groupby("genus", sort=False):
                self.local_genus_counts[(island, genus)] = np.bincount(
                    members.profile_id, minlength=n_profiles).astype(float)
        self.global_genus_cache = {}
        self.local_genus_cache = {}
        self.local_total_cache = {}
        self.n_donor_names = int(global_species.accepted_species.nunique())
        self.n_donor_genera = len(self.global_genus_counts)
        self.n_local_source_names = len(local)
        self.n_source_entities = int(assignments.entity_ID.nunique())
        self.n_local_source_islands = len(self.local_total_counts)

    def _global_genus(self, genus: str) -> np.ndarray:
        if genus not in self.global_genus_cache:
            counts = self.global_genus_counts.get(genus)
            if counts is None:
                self.global_genus_cache[genus] = self.global_prior
            else:
                self.global_genus_cache[genus] = (
                    (counts + self.shrink * self.global_prior) /
                    (counts.sum() + self.shrink)
                )
        return self.global_genus_cache[genus]

    def probabilities(self, island: str, genus: str) -> np.ndarray:
        if island not in self.local_total_counts:
            raise ValueError("Matched island without eligible GIFT source donors")
        if island not in self.local_total_cache:
            counts = self.local_total_counts[island]
            self.local_total_cache[island] = (
                (counts + self.shrink * self.global_prior) /
                (counts.sum() + self.shrink)
            )
        key = (island, genus)
        if key not in self.local_genus_cache:
            base = self._global_genus(genus)
            counts = self.local_genus_counts.get(key)
            if counts is None:
                self.local_genus_cache[key] = base
            else:
                self.local_genus_cache[key] = (
                    (counts + self.shrink * base) /
                    (counts.sum() + self.shrink)
                )
        return np.stack((
            self.global_prior, self.local_total_cache[island],
            self._global_genus(genus), self.local_genus_cache[key],
        ))

    def genus_support(self, island: str, genus: str) -> tuple[bool, bool]:
        return (genus in self.global_genus_counts,
                (island, genus) in self.local_genus_counts)


def predict_one_island(observed: pd.DataFrame, model: DonorModel,
                       island: str, space: ProfileSpace) -> tuple[dict, np.ndarray, np.ndarray]:
    """Proper name-level scores and approximate normalized expected occupancy.

    Approximation: normalize expected binary category incidence; this is NOT
    the same as the expectation of normalized incidence.
    """
    if observed.empty or observed.accepted_species.duplicated().any():
        raise ValueError("Island must have distinct complete taxon names")
    n_profiles = len(space.ids)
    losses = np.zeros((len(MODELS), 2), dtype=float)
    log_absent = np.zeros((len(MODELS), n_profiles))
    global_n = local_n = 0
    ids = []
    for genus, profile in observed[["genus", "profile_id"]].itertuples(index=False, name=None):
        profile = int(profile)
        if profile not in space.lookup:
            raise ValueError("Unrecognized profile ID in held-out flora")
        index = space.lookup[profile]
        ids.append(index)
        glob, local = model.genus_support(island, genus)
        global_n += int(glob)
        local_n += int(local)
        probabilities = model.probabilities(island, genus)
        if not np.allclose(probabilities.sum(axis=1), 1, atol=1e-10):
            raise ValueError("Prediction probabilities failed normalization")
        losses[:, 0] -= np.log(np.clip(probabilities[:, index], 1e-13, None))
        losses[:, 1] += probabilities @ space.matrix[index]
        log_absent += np.log1p(-np.clip(probabilities, 0, 1 - 1e-12))
    expected_incidence = -np.expm1(log_absent)
    expected_normalized = (
        expected_incidence / expected_incidence.sum(axis=1, keepdims=True)
    )
    observed_normalized = np.zeros(n_profiles, dtype=float)
    distinct = np.unique(ids)
    observed_normalized[distinct] = 1 / len(distinct)
    kernel = 1 - space.matrix
    metrics = {}
    for i, key in enumerate(MODELS):
        err = expected_normalized[i] - observed_normalized
        metrics[key] = {
            "log_loss": float(losses[i, 0] / len(ids)),
            "trait_distance_loss": float(losses[i, 1] / len(ids)),
            "repertoire_L2_error": float(err @ err),
            "repertoire_kernel_error": float(err @ kernel @ err),
            "expected_number_profiles": float(expected_incidence[i].sum()),
        }
    audit = {
        "n_names": len(ids), "n_profiles": len(distinct),
        "global_genus_support": global_n / len(ids),
        "local_genus_support": local_n / len(ids),
    }
    return {"metrics": metrics, "audit": audit}, observed_normalized, expected_normalized


def _paired_uncertainty(values: np.ndarray, *, seed: int, draws: int) -> tuple[float, float]:
    if len(values) < 2 or draws < 99:
        raise ValueError("Too few independent matched pairs or bootstrap draws")
    rng = np.random.default_rng(seed)
    ids = rng.integers(len(values), size=(draws, len(values)))
    means = values[ids].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def _origin_residual_records(
    scope: str, matches: list[dict], observed: dict,
    predicted: dict, space: ProfileSpace, *,
    swaps: int, seed: int,
) -> list[dict]:
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1., 1.], size=(swaps, len(matches)))
    delta_obs = np.stack([
        observed[p["oceanic_id"]] - observed[p["continental_id"]]
        for p in matches
    ])
    kernel = 1 - space.matrix
    matrices = {
        "exact_profile_L2": np.eye(len(space.ids)),
        "categorical_similarity_kernel": kernel,
    }
    results = []
    for idx, model in enumerate(MODELS):
        delta_pred = np.stack([
            predicted[p["oceanic_id"]][idx] -
            predicted[p["continental_id"]][idx]
            for p in matches
        ])
        for response, delta in (
            ("observed", delta_obs),
            ("predicted", delta_pred),
            ("source_adjusted_residual", delta_obs - delta_pred),
        ):
            for endpoint, km in matrices.items():
                obs, reference = _paired_swap_statistic(delta, km, signs)
                results.append({
                    "scope": scope, "source_model": model,
                    "response": response, "endpoint": endpoint,
                    "n_pairs": len(matches), "observed_statistic": float(obs),
                    "swap_null_median": float(np.median(reference)),
                    "swap_null_q025": float(np.quantile(reference, .025)),
                    "swap_null_q975": float(np.quantile(reference, .975)),
                    "p_swap_upper": float(_upper_tail(obs, reference)),
                })
    return results


def run(private_all: Path, private_direct: Path, flora: Path,
        source_assignments: Path, mainland_native: Path,
        crosswalk: Path, g0_root: Path, reference: Path, output: Path,
        *, source_mode: str = "geo_k5", shrink: float = DEFAULT_SHRINK,
        native_only: bool = False, min_native_names: int = 20,
        bootstrap_draws: int = 4999, permutations: int = 9999,
        seed: int = 20261108, strict: bool = True) -> dict:
    if source_mode not in VALID_MODES:
        raise ValueError("Unknown GIFT source assignment mode")
    if bootstrap_draws < 99 or permutations < 999 or min_native_names < 1:
        raise ValueError("Analysis support or draw count outside allowed range")
    private_all, private_direct, flora, source_assignments, mainland_native = map(
        Path, (private_all, private_direct, flora, source_assignments, mainland_native))
    crosswalk, g0_root, reference, output = map(
        Path, (crosswalk, g0_root, reference, output))
    output.mkdir(parents=True, exist_ok=True)
    space = ProfileSpace(load_profiles(reference))
    profiles = {
        "all": load_private_maps(private_all, "all", set(space.ids), strict=strict),
        "direct": load_private_maps(private_direct, "direct", set(space.ids), strict=strict),
    }
    geology = pd.read_csv(crosswalk, dtype={"island_id": str})
    if geology.island_id.duplicated().any():
        raise ValueError("Duplicate geological geography island ID")
    g0 = {
        scope: pd.read_csv(g0_root / f"{scope}_islands.csv", dtype={"island_id": str})
        for scope in ("all", "direct")
    }
    pairs, matching = match_origin_islands(
        g0["all"].to_dict("records"), g0["direct"].to_dict("records"),
        geology.set_index("island_id").to_dict("index"),
    )
    if strict and (len(pairs) != 37 or
                   sum(p["region"] == "tropical" for p in pairs) != 14):
        raise ValueError("Reproduced origin-matched cohort does not equal frozen 37 pairs")
    all_ids = {p[key] for p in pairs for key in ("continental_id", "oceanic_id")}
    source_status = pd.read_csv(
        flora, usecols=["island_id", "accepted_species", "origin_status"], dtype=str)
    membership = {}
    for scope in ("all", "direct"):
        support = source_status.loc[source_status.island_id.isin(all_ids)]
        joined = support.merge(profiles[scope], on="accepted_species",
                               how="inner", validate="many_to_one")
        if joined.duplicated(["island_id", "accepted_species"]).any():
            raise ValueError("Duplicate private island taxon memberships")
        membership[scope] = joined
        # Before native/status filtering, exact original PR321 complete-name
        # and distinct-profile counts must agree on every selected island.
        # This is a fail-closed validation of the restricted floristic join.
        subset_expected = g0[scope].loc[
            g0[scope].island_id.isin(all_ids),
            ["island_id", "n_names", "n_profiles"],
        ]
        check = validate_island_counts(joined, subset_expected, scope)
        if check["eligible_islands_checked"] != len(all_ids):
            raise ValueError("Some matched island profile counts were not audited")
    if native_only:
        native = {s: part.loc[part.origin_status == "native"].copy()
                  for s, part in membership.items()}
        counts = {s: p.groupby("island_id").size().to_dict()
                  for s, p in native.items()}
        pairs = [p for p in pairs if all(
            counts[s].get(p[key], 0) >= min_native_names
            for s in counts for key in ("continental_id", "oceanic_id")
        )]
        if strict and (len(pairs) != 28 or
                       sum(p["region"] == "tropical" for p in pairs) != 8):
            raise ValueError("Native-only original matched cohort did not reproduce 28 pairs")
        if len(pairs) < 8:
            raise ValueError("Insufficient explicitly native floras")
        membership = native
    matched_ids = {p[key] for p in pairs for key in ("continental_id", "oceanic_id")}
    assignments = pd.read_csv(
        source_assignments,
        usecols=["island_id", "source_mode", "entity_ID"],
        dtype={"island_id": str, "entity_ID": str, "source_mode": str})
    assignments = assignments.loc[
        assignments.island_id.isin(matched_ids) &
        assignments.source_mode.eq(source_mode),
        ["island_id", "entity_ID"],
    ].drop_duplicates()
    if assignments.island_id.nunique() != len(matched_ids):
        raise ValueError("One or more matched islands lack the selected mainland source")
    source_native = pd.read_csv(
        mainland_native, usecols=["entity_ID", "work_species"],
        dtype={"entity_ID": str, "work_species": str})
    aggregate_rows, improvements, coverage, residuals = [], [], [], []
    source_audits = {}
    for scope in ("all", "direct"):
        sample = membership[scope].loc[
            membership[scope].island_id.isin(matched_ids)].copy()
        if sample.island_id.nunique() != len(matched_ids):
            raise ValueError("Some matched islands lack complete taxon profiles")
        heldout = set(sample.accepted_species)
        train = source_native.merge(
            profiles[scope], left_on="work_species", right_on="accepted_species",
            how="inner").drop(columns=["work_species"])
        train = train.loc[~train.accepted_species.isin(heldout)]
        if set(train.accepted_species) & heldout:
            raise AssertionError("Taxon-name holdout leaked target into donor training")
        if train.empty:
            raise ValueError("No remaining donor references after whole-cohort taxon holdout")
        predictor = DonorModel(train, assignments, len(space.ids), shrink=shrink)
        source_audits[scope] = {
            "n_heldout_taxon_names": len(heldout),
            "n_test_taxon_island_records": len(sample),
            "n_reference_native_taxon_names": predictor.n_donor_names,
            "n_global_donor_genera": predictor.n_donor_genera,
            "n_nearby_source_taxon_island_records": predictor.n_local_source_names,
            "n_assigned_mainland_entities": predictor.n_source_entities,
        }
        observed_vectors, predicted_vectors = {}, {}
        scores, audit_rows = {}, []
        for island, part in sample.groupby("island_id", sort=True):
            summary, obs, pred = predict_one_island(part, predictor, island, space)
            observed_vectors[island] = obs
            predicted_vectors[island] = pred
            scores[island] = summary["metrics"]
            origin = ("continental" if island in
                      {p["continental_id"] for p in pairs} else "oceanic")
            audit_rows.append({"origin": origin, **summary["audit"]})
        for origin in ("continental", "oceanic"):
            rows = [r for r in audit_rows if r["origin"] == origin]
            coverage.append({
                "scope": scope, "origin": origin,
                "n_islands": len(rows),
                "median_complete_name_count": float(np.median([r["n_names"] for r in rows])),
                "mean_global_genus_supported_fraction": float(
                    np.mean([r["global_genus_support"] for r in rows])),
                "mean_local_genus_supported_fraction": float(
                    np.mean([r["local_genus_support"] for r in rows])),
            })
        for model in MODELS:
            for endpoint in LOSSES:
                values = [scores[island][model][endpoint] for island in sorted(scores)]
                aggregate_rows.append({
                    "scope": scope, "model": model, "endpoint": endpoint,
                    "n_islands": len(values), "mean_loss": float(np.mean(values)),
                    "median_loss": float(np.median(values)),
                })
        for better, worse in COMPARISONS:
            for endpoint in LOSSES:
                differences = np.asarray([
                    np.mean([
                        scores[p[key]][worse][endpoint] -
                        scores[p[key]][better][endpoint]
                        for key in ("continental_id", "oceanic_id")
                    ]) for p in pairs
                ])
                low, high = _paired_uncertainty(
                    differences, draws=bootstrap_draws, seed=seed)
                improvements.append({
                    "scope": scope, "better_model": better, "worse_model": worse,
                    "endpoint": endpoint, "n_pairs": len(pairs),
                    "mean_improvement": float(np.mean(differences)),
                    "pair_bootstrap_q025": low, "pair_bootstrap_q975": high,
                })
        residuals.extend(_origin_residual_records(
            scope, pairs, observed_vectors, predicted_vectors, space,
            swaps=permutations, seed=seed))
    # Exploratory residual family for a stated genus+nearest-source predictor
    # (All/Direct × exact-category and category-distance outcomes).
    primary = [r for r in residuals if (
        r["source_model"] == "local_genus" and
        r["response"] == "source_adjusted_residual")]
    if len(primary) != 4:
        raise AssertionError("Residual primary four-hypothesis family changed")
    q = _bh([r["p_swap_upper"] for r in primary])
    harmonic4 = 1 + .5 + 1/3 + .25
    for r, qv in zip(primary, q):
        r["q_bh_primary4"] = float(qv)
        r["q_by_primary4"] = float(min(1., qv * harmonic4))
    _write_csv(output / "heldout_model_scores.csv", aggregate_rows, SCORE_COLUMNS)
    _write_csv(output / "paired_prediction_improvements.csv", improvements, IMPROVEMENT_COLUMNS)
    _write_csv(output / "mainland_genus_support.csv", coverage, COVERAGE_COLUMNS)
    _write_csv(output / "origin_repertoire_residuals.csv", residuals, RESIDUAL_COLUMNS)
    manifest = {
        "schema_version": "fd_heldout_gift_mainland_genus_v1",
        "status": "completed_private_input_aggregate_only_exploratory",
        "source_mode": source_mode, "native_only": native_only,
        "n_matched_pairs": len(pairs), "n_test_islands": len(matched_ids),
        "regions": {region: sum(p["region"] == region for p in pairs)
                    for region in REGIONS},
        "source_audits": source_audits,
        "source_holdout": "Every exact name in any matched tested island is removed globally from every GIFT donor pool within that evidence scope",
        "models": {
            "global_mainland": "single global GIFT native profile prior",
            "local_mainland": "same profile prior updated with sampled potential mainland source taxa",
            "global_genus": "native GIFT genus profile counts updated from global prior",
            "local_genus": "local GIFT source genus profile counts updated from global genus model",
        },
        "parameters": {"global_prior_per_profile": PSEUDO_GLOBAL_PER_PROFILE,
                       "hierarchical_shrinkage": shrink,
                       "matched_pair_bootstrap_draws": bootstrap_draws,
                       "paired_label_swaps": permutations, "seed": seed},
        "metrics": {
            "log_loss": "true held-out name-level categorical log score, each name occurrence counted within its island",
            "trait_distance_loss": "expected untransformed colour/form set-Jaccard distance from predicted profiles to recorded held-out profile",
            "repertoire_error": "observed unique profile incidence versus normalized expectation of category occupancy, approximate",
            "residual_pair_test": "unadjusted origin differences remaining after model-predicted normalized repertoire subtraction, paired sign-swap reference",
        },
        "claim_limits": [
            "The known island genus is an input; this does not forecast which genera colonize",
            "Complete name-level holdout is not phylogenetic or historical source-region holdout",
            "Actual island source ancestors and native statuses not fully reconstructed",
            "Observed name occurrences may repeat across matched islands, so pair bootstrap CIs are conditional on these pairs",
            "Nearest-source information can be sparse and strongly influenced by smoothing hyperparameters",
            "Residual paired origin association is not a fully genus/source-conditioned permutation null",
            "No observed pollinator function, adaptation or evolutionary mechanism identified",
        ],
        "data_release": "aggregate output only, no island IDs, mainland source entity IDs, species or genus names, taxon-to-trait rows, or pair identities",
        "source_sha256": {
            "private_all": _sha(private_all), "private_direct": _sha(private_direct),
            "private_flora": _sha(flora),
            "GIFT_native": _sha(mainland_native),
            "GIFT_assignments": _sha(source_assignments),
            "geology_crosswalk": _sha(crosswalk), "public_profile_reference": _sha(reference),
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "n_matched_pairs": len(pairs), "native_only": native_only,
        "global_genus_log_loss": [
            {"scope": r["scope"], "improvement": r["mean_improvement"]}
            for r in improvements if r["better_model"] == "global_genus"
            and r["endpoint"] == "log_loss"
        ],
        "primary_residuals": [
            {"scope": r["scope"], "endpoint": r["endpoint"],
             "p": r["p_swap_upper"], "BH": r.get("q_bh_primary4")}
            for r in primary
        ],
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("private-all", "private-direct", "flora", "source-assignments",
                 "mainland-native", "crosswalk", "g0", "reference", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--source-mode", choices=VALID_MODES, default="geo_k5")
    p.add_argument("--shrink", type=float, default=DEFAULT_SHRINK)
    p.add_argument("--bootstrap-draws", type=int, default=4999)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--native-only", action="store_true")
    p.add_argument("--min-native-names", type=int, default=20)
    p.add_argument("--allow-synthetic", action="store_true")
    args = p.parse_args()
    run(
        args.private_all, args.private_direct, args.flora,
        args.source_assignments, args.mainland_native, args.crosswalk,
        args.g0, args.reference, args.output,
        source_mode=args.source_mode, shrink=args.shrink,
        bootstrap_draws=args.bootstrap_draws, permutations=args.permutations,
        native_only=args.native_only, min_native_names=args.min_native_names,
        strict=not args.allow_synthetic,
    )


if __name__ == "__main__":
    main()
