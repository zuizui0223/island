"""Two-stage mainland source-genus opportunity to floral repertoire prediction.

Authorized PRIVATE input: complete taxon→flower-profile table, matched island
taxon membership, original GIFT native mainland floras and source assignment
identities. Training removes ALL 74 matched-island taxon names globally.
Public output: aggregate model error, pair-bootstrap improvements and digest
manifest ONLY. Original species, genus, island and donor-entity IDs stay local.

Source-only predictors do NOT use known island genera. A separately labeled
'known_genus_oracle' benchmark does, to bound missing lineage-assembly skill.
Unrecorded genera are NOT presumed biologically absent.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_profiles
from .heldout_mainland_prediction import _paired_uncertainty
from .source_genus_opportunity import (_native_mainland_source_data, _sha, exclude_heldout_species)

MODELS = ("global_prior", "source_genus_density",
          "source_genus_prevalence", "known_genus_oracle")
ENDPOINTS = ("name_log_loss", "name_trait_distance",
             "frequency_L2", "frequency_kernel")
PSEUDO_GLOBAL = 0.02
SHRINK = 8.0


def source_genus_weights(
    assignment: list[tuple[str, float]],
    source_maps: dict[str, dict[str, int]],
    source_sizes: dict[str, int],
    global_prevalence: dict[str, int],
) -> tuple[list[str], dict[str, np.ndarray]]:
    rows = sorted(
        [(ent, float(km)) for ent, km in assignment
         if ent in source_maps and np.isfinite(km) and km >= 0],
        key=lambda x: (x[1], x[0]),
    )
    if len(rows) < 3:
        raise ValueError("Fewer than 3 usable mainland GIFT source regions")
    genera = sorted(set().union(*(set(source_maps[entity]) for entity, _ in rows)))
    if len(genera) < 5:
        raise ValueError("Too few source-eligible genera")
    region_weights = np.asarray(
        [1.0 / np.sqrt(1.0 + km / 500.0) for _, km in rows], dtype=float)
    region_weights /= region_weights.sum()
    local_density = np.asarray([
        sum(w * source_maps[entity].get(genus, 0) / max(source_sizes[entity], 1)
            for (entity, _), w in zip(rows, region_weights))
        for genus in genera
    ], dtype=float)
    global_source_prevalence = np.asarray(
        [global_prevalence.get(genus, 0) for genus in genera], dtype=float)
    if local_density.sum() <= 0 or global_source_prevalence.sum() <= 0:
        raise ValueError("All potential donor lineage scores were zero")
    return genera, {
        "source_genus_density": local_density / local_density.sum(),
        "source_genus_prevalence": (
            global_source_prevalence / global_source_prevalence.sum()),
    }


def genus_flower_model(donors: pd.DataFrame, n_profiles: int,
                       *, shrinkage: float = SHRINK
                       ) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Globally coherent per-genus profile posterior, from held-out donors."""
    if donors.empty or shrinkage <= 0 or not np.isfinite(shrinkage):
        raise ValueError("Need non-empty donors and positive shrinkage")
    donor = donors.drop_duplicates("accepted_species")
    if donor.profile_id.min() < 0 or donor.profile_id.max() >= n_profiles:
        raise ValueError("Donor profile outside reference")
    counts = np.bincount(donor.profile_id, minlength=n_profiles).astype(float)
    base = (counts + PSEUDO_GLOBAL) / (counts.sum() + n_profiles * PSEUDO_GLOBAL)
    by_genus = {}
    for genus, frame in donor.groupby("genus", sort=False):
        counts = np.bincount(frame.profile_id, minlength=n_profiles).astype(float)
        by_genus[genus] = (counts + shrinkage * base) / (
            counts.sum() + shrinkage)
    return base, by_genus


def predict_island_flora(
    island_membership: pd.DataFrame, assignment: list[tuple[str,float]],
    source_maps: dict[str,dict[str,int]], source_sizes: dict[str,int],
    global_prevalence: dict[str,int], prior: np.ndarray,
    genus_pdfs: dict[str,np.ndarray],
) -> tuple[dict[str,dict[str,float]],dict]:
    """Flower *frequency* distribution, not normalized unique-profile count.

    All actual island genus records enter only observed evaluation and the
    explicitly post-hoc oracle benchmark, not source-only model selection.
    """
    if island_membership.empty or island_membership.accepted_species.duplicated().any():
        raise ValueError("Missing or duplicate island taxon records")
    candidates, weights = source_genus_weights(
        assignment, source_maps, source_sizes, global_prevalence)
    source_profiles = np.stack([genus_pdfs.get(g, prior) for g in candidates])
    lineage_counts = island_membership.groupby("genus").size()
    oracle = np.stack([genus_pdfs.get(g, prior) for g in lineage_counts.index])
    oracle_weights = lineage_counts.to_numpy(dtype=float)
    oracle_weights /= oracle_weights.sum()
    estimates = {
        "global_prior": prior,
        "source_genus_density": weights["source_genus_density"] @ source_profiles,
        "source_genus_prevalence": weights["source_genus_prevalence"] @ source_profiles,
        "known_genus_oracle": oracle_weights @ oracle,
    }
    observed_counts = np.bincount(island_membership.profile_id, minlength=len(prior))
    observed = observed_counts / observed_counts.sum()
    scores = {}
    for model, probs in estimates.items():
        if not np.isclose(probs.sum(), 1, atol=1e-10):
            raise ValueError("Predicted mixture does not normalize")
        probs = np.maximum(probs, 1e-15)
        probs /= probs.sum()
        scores[model] = (observed, probs)
    audit = {
        "n_test_names": len(island_membership),
        "n_observed_genera": int(lineage_counts.size),
        "n_source_candidate_genera": len(candidates),
        "n_source_candidate_with_flower_donor": sum(g in genus_pdfs for g in candidates),
    }
    return scores, audit


def score_flora(observed: np.ndarray, predicted: np.ndarray,
                distance: np.ndarray) -> dict[str,float]:
    if not np.isfinite(predicted).all() or np.any(predicted <= 0):
        raise ValueError("Nonpositive predicted category probability")
    residual = observed - predicted
    similarity = 1 - distance
    return {
        "name_log_loss": float(-observed @ np.log(predicted)),
        "name_trait_distance": float(observed @ distance @ predicted),
        "frequency_L2": float(residual @ residual),
        "frequency_kernel": float(residual @ similarity @ residual),
    }


def run(profile_all: Path, profile_direct: Path,
        membership_all: Path, membership_direct: Path,
        mainland_native: Path, assignments: Path, matched_pairs: Path,
        reference: Path, output: Path, *,
        shrinkage: float = SHRINK, source_mode: str = "geo_k5",
        bootstrap_draws: int = 4999, seed: int = 20261108) -> dict:
    paths = list(map(Path, (profile_all, profile_direct, membership_all,
                            membership_direct, mainland_native,
                            assignments, matched_pairs, reference)))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if bootstrap_draws < 99:
        raise ValueError("Insufficient bootstrap draws")
    matched = pd.read_csv(matched_pairs, dtype=str)
    if len(matched) != 37 or not {"continental_id", "oceanic_id", "region"}.issubset(matched):
        raise ValueError("Must use frozen 37-pair geological matching table")
    matched_islands = set(matched.continental_id) | set(matched.oceanic_id)
    if len(matched_islands) != 74:
        raise ValueError("Paired island endpoint reuse or missing IDs")
    reference_space = ProfileSpace(load_profiles(reference))
    n_profiles = len(reference_space.ids)
    if n_profiles != 496 or reference_space.ids != list(range(n_profiles)):
        raise ValueError("Changed profile ID reference: must be frozen 0..495")
    gift = pd.read_csv(
        mainland_native, usecols=["entity_ID", "work_species"], dtype=str)
    source_maps, source_sizes, global_prevalence, source_native = (
        _native_mainland_source_data(gift))
    assign = pd.read_csv(
        assignments,
        usecols=["island_id","entity_ID","source_mode","source_distance_km"],
        dtype=str)
    assign = assign.loc[
        assign.island_id.isin(matched_islands) & assign.source_mode.eq(source_mode)]
    nearby = {}
    for island, frame in assign.groupby("island_id"):
        nearby[island] = []
        for row in frame.itertuples(index=False):
            try:
                nearby[island].append(
                    (row.entity_ID, float(row.source_distance_km)))
            except (ValueError, TypeError):
                continue
    if set(nearby) != matched_islands:
        raise ValueError("One matched island has no nearest mainland reference")
    metrics = []
    improvements = []
    source_audit = []
    for scope, profpath, memberpath in (
        ("all", profile_all, membership_all),
        ("direct", profile_direct, membership_direct),
    ):
        catalog = pd.read_csv(
            profpath, usecols=["accepted_species","genus","profile_id"],
            dtype={"accepted_species":str, "genus":str,"profile_id":int})
        membership = pd.read_csv(
            memberpath,
            usecols=["island_id","accepted_species","genus","profile_id"],
            dtype={"island_id":str, "accepted_species":str,"genus":str,"profile_id":int})
        membership = membership.loc[membership.island_id.isin(matched_islands)].copy()
        if membership.island_id.nunique() != 74:
            raise ValueError("One matched island lacks floral profiles")
        if membership.duplicated(["island_id","accepted_species"]).any():
            raise ValueError("Private flora duplicated island/name")
        if (membership.genus != membership.accepted_species.str.split().str[0]).any():
            raise ValueError("Private flora has genus/name mismatch")
        if catalog.accepted_species.duplicated().any():
            raise ValueError("Profile catalog repeats taxon names")
        names_heldout = set(membership.accepted_species)
        # Strict independent taxon-name holdout applies also to *source genus
        # opportunity* features, not merely to floral donor profile counts.
        filtered_native = exclude_heldout_species(source_native, names_heldout)
        filtered_source_map, filtered_source_sizes, filtered_global_prevalence, _ = (
            _native_mainland_source_data(filtered_native))
        donors = filtered_native.merge(
            catalog[["accepted_species","profile_id"]], left_on="work_species", right_on="accepted_species",
            how="inner", validate="many_to_one")
        donors = donors.loc[~donors.accepted_species.isin(names_heldout)].copy()
        if set(donors.accepted_species) & names_heldout:
            raise AssertionError("Test taxa leaked into mainland profile donors")
        donor_names = donors.drop_duplicates("accepted_species")
        prior, genus_pdfs = genus_flower_model(
            donor_names, n_profiles, shrinkage=shrinkage)
        source_audit.append({
            "scope": scope, "n_heldout_names": len(names_heldout),
            "n_heldout_name_island_records": len(membership),
            "n_GIFT_source_entity_name_records_after_holdout": len(filtered_native),
            "n_training_donor_names": len(donor_names),
            "n_donor_genera": len(genus_pdfs),
            "n_GIFT_source_regions": int(source_native.entity_ID.nunique()),
        })
        scores = {}
        for island, observed in membership.groupby("island_id"):
            preds, audit = predict_island_flora(
                observed, nearby[island], filtered_source_map, filtered_source_sizes,
                filtered_global_prevalence, prior, genus_pdfs)
            for model, (actual, predicted) in preds.items():
                values = score_flora(actual, predicted, reference_space.matrix)
                scores[island, model] = values
                metrics.append({"scope": scope, "model": model, **values,
                                "n_source_candidates": audit["n_source_candidate_genera"],
                                "n_source_genus_trait_donors": audit[
                                    "n_source_candidate_with_flower_donor"]})
        for model in MODELS:
            subset = [r for r in metrics if r["scope"] == scope and r["model"] == model]
            # Collapse to island means; no taxon or genus name is output.
            if len(subset) != 74:
                raise ValueError("Expected 74 evaluated island floras")
        for model in MODELS[1:]:
            for endpoint in ENDPOINTS:
                diffs = np.asarray([
                    np.mean([
                        scores[getattr(row, key), "global_prior"][endpoint] -
                        scores[getattr(row, key), model][endpoint]
                        for key in ("continental_id", "oceanic_id")
                    ])
                    for row in matched.itertuples(index=False)
                ], dtype=float)
                lo, hi = _paired_uncertainty(
                    diffs, seed=seed, draws=bootstrap_draws)
                improvements.append({
                    "scope": scope, "model": model, "baseline": "global_prior",
                    "endpoint": endpoint, "n_matched_pairs": len(matched),
                    "mean_improvement": float(diffs.mean()),
                    "pair_bootstrap_q025": float(lo),
                    "pair_bootstrap_q975": float(hi),
                })
    raw_metrics = pd.DataFrame(metrics)
    rows = (raw_metrics.groupby(["scope","model"],as_index=False)
            .agg(n_islands=("name_log_loss","size"),
                 mean_name_log_loss=("name_log_loss","mean"),
                 mean_name_trait_distance=("name_trait_distance","mean"),
                 mean_frequency_L2=("frequency_L2","mean"),
                 mean_frequency_kernel=("frequency_kernel","mean"),
                 median_source_candidate_genera=("n_source_candidates","median"),
                 median_source_genus_trait_donors=("n_source_genus_trait_donors","median")))
    rows.to_csv(output / "two_stage_source_prediction_scores.csv", index=False)
    pd.DataFrame(improvements).to_csv(
        output / "two_stage_paired_improvements.csv", index=False)
    manifest = {
        "schema": "fd_two_stage_source_genus_flora_prediction_v1",
        "status": "conditional_recorded_flora_prediction_not_colonization",
        "matched_island_pairs": len(matched),
        "source_mode": source_mode, "shrinkage": shrinkage,
        "bootstrap_draws": bootstrap_draws,
        "n_reference_profiles": n_profiles,
        "source_audit": source_audit,
        "source_only_models": list(MODELS[:3]),
        "oracle_model_uses_known_island_genus": MODELS[3],
        "source_name_holdout": "all target-island taxon names globally removed from GIFT mainland source genus opportunity features AND floral donors before prediction",
        "prediction_unit": "observed taxon-name carrier-weighted recorded flower-category frequencies, conditional on evaluation island flora",
        "not_prediction_of": ["which genus truly colonizes",
                              "new arrivals", "true absences",
                              "taxon-specific dispersal",
                              "pollination function"],
        "source_hashes": {
            name: _sha(path) for name, path in zip(
                ("all_profiles","direct_profiles","all_membership",
                 "direct_membership","GIFT_native","mainland_assignments",
                 "matched_pairs","reference"), paths)
        },
        "privacy": "aggregate-only outputs without taxon, genus, island or source IDs",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "status": manifest["status"], "n_pairs": len(matched),
        "results": rows.to_dict("records"),
        "source_genus_improvements": [r for r in improvements
                                      if r["model"] == "source_genus_density"
                                      and r["endpoint"] == "frequency_kernel"],
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("profile-all", "profile-direct", "membership-all",
                 "membership-direct", "mainland-native", "assignments",
                 "matched-pairs", "reference", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--source-mode", default="geo_k5")
    parser.add_argument("--shrinkage", type=float, default=SHRINK)
    parser.add_argument("--bootstrap-draws", type=int, default=4999)
    args = parser.parse_args()
    run(args.profile_all, args.profile_direct,
        args.membership_all, args.membership_direct,
        args.mainland_native, args.assignments, args.matched_pairs,
        args.reference, args.output, source_mode=args.source_mode,
        shrinkage=args.shrinkage, bootstrap_draws=args.bootstrap_draws)


if __name__ == "__main__":
    main()
