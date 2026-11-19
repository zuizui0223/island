"""Original 4,505-GIFT flora geography->genus->Direct floral-category source pilot.

This runner needs the SHA-pinned original archival GIFT source restored in the
island GitHub Actions runner. It neither prints nor exports accepted taxon
names, original islands, candidate source IDs or geographic coordinates.

Stage 1 re-fits botanical genus frequencies using only OTHER eligible island
floras: unlike the previous geographic-only genus run it removes any accepted
plant name from the full focal 10-degree block OR focal first-source group,
including names recorded again on distant donor islands. Stage 2 independently
learns P(Direct flower category | genus) from disjoint, held-out external
region species under the same entire-full-flora name ban. Tests are on
entire target island Direct carrier sets; all four stage-1 models are scored
on identical real plants. This is floristic prediction, not measured dispersal,
local adaptation, pollinator selection or a mainland ancestor study.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from fd_island_functional_beta import strict_traits
from fd_global_genus_assembly_geometry import source_endpoints
from fd_global_original_oracle_free_genus_mix import refit_no_oracle_regional_genus_source
from fd_global_geo_source_training_holdout import (
    original_full_flora_group_blacklists, source_training_ban, _source_group,
)
from fd_interregional_genus_name_exclusive_transfer import (
    genus, donor_source_region_ownership, sha, REGIONS,
    SCOPES, MIN_DIRECT,
)
from fd_original_environmental_genus_composition import (
    _floats, _xyz, _jaccard, _nearest_probabilities,
    K, DELTA_UNKNOWN, SOURCE_WEIGHT, FEATURES,
)
from fd.two_stage_floral_forecast import (
    REQUIRED_MODELS, UNKNOWN_GENUS, score_focal_flora,
    paired_cluster_summary, dual_loss_dual_group_positive,
)

PILOT_ENDPOINTS = ("floral_form", "flower_primary_color__floral_form")
PRIMARY_REGIONS = ("northern_midlatitude", "tropical", "southern_extratropical")
MODEL_DISTANCES = ("geographic_nearest", "environment_nearest",
                   "environment_source_nearest")
NBOOT = 499
MIN_SOURCE_ISLANDS = 30
MIN_SOURCE_BLOCKS = 10
MIN_SOURCE_SPECIES = 10
MIN_SOURCE_GENERA = 5
SOURCE_PINS = {
    "flora": "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    "assign": "8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
    "traits": "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
}
SUMMARY_GAIN_COLUMNS = tuple(
    f"{m}_gain_vs_{base}_{loss}"
    for m in REQUIRED_MODELS
    for base in ("external_blind", "regional_blind", "regional_genus_prior")
    if m != base
    for loss in ("Brier", "logloss")
) + tuple(
    f"{name}_{loss}" for name in
    ("source_increment", "source_over_geo", "source_over_environment")
    for loss in ("Brier", "logloss")
)


EXPANDED_GAIN_COLUMNS = tuple(
    f"expanded_{m}_gain_vs_{base}_{loss}"
    for m in REQUIRED_MODELS
    for base in ("regional_blind", "external_blind", "regional_genus_prior",
                 "constrained")
    if m != base
    for loss in ("Brier", "logloss")
)


def extend_external_only_congener_profiles(
    constrained_profiles: dict[str, np.ndarray],
    exclusive_external_source: set[str],
    forbidden_full_flora_names: set[str],
    states: dict[str, str],
    indices: dict[str, int],
    blind: np.ndarray,
    *, alpha: float = 5.0, max_donors: int = 20,
) -> tuple[dict[str, np.ndarray], int]:
    """Add external trait genera NOT represented by local Direct donor genera.

    Keep every existing local-congener-supported genus profile exactly
    unchanged, including its original equal-local/external-name sample cap.
    Hence the paired sensitivity isolates the EFFECT OF THE PREVIOUS
    LOCAL-GENUS-SUPPORT GATE rather than changing the learned profiles of
    genera already permitted. Candidate species names are original
    source-private; return only anonymous counts/profile probabilities.
    """
    if alpha <= 0 or max_donors < 1:
        raise ValueError("External genus shrinkage settings invalid")
    if (not np.isclose(np.sum(blind), 1)
            or any(not np.isclose(np.sum(p), 1)
                   for p in constrained_profiles.values())):
        raise ValueError("External genus blind/profile simplex invalid")
    accepted = (set(exclusive_external_source) - forbidden_full_flora_names) & set(states)
    if accepted & forbidden_full_flora_names:
        raise AssertionError("Forbidden accepted plant in external congener model")
    extra = defaultdict(list)
    for name in accepted:
        g = genus(name)
        if g not in constrained_profiles:
            extra[g].append(name)
    out = dict(constrained_profiles)
    for g, names in extra.items():
        names = sorted(names, key=lambda sp: hashlib.sha256(
            ("fd:20261119:expanded-external:" + sp).encode()).hexdigest())
        counts = np.zeros(len(blind), dtype=float)
        for name in names[:max_donors]:
            counts[indices[states[name]]] += 1.0
        n = min(len(names), max_donors)
        p = (counts + alpha * blind) / (n + alpha)
        if not np.isclose(p.sum(), 1):
            raise AssertionError("Expanded source genus conditional not normalized")
        out[g] = p
    return out, len(extra)


def filtered_genus_matrix(
    donor_islands: list[str],
    island_species: dict[str, set[str]],
    forbidden: set[str],
    vocab: dict[str, int],
    *, min_species: int = MIN_SOURCE_SPECIES,
    min_genera: int = MIN_SOURCE_GENERA,
) -> tuple[list[str], csr_matrix, set[str], int]:
    """Construct genuine name-disjoint per-island genus frequency rows."""
    rows, cols, vals, kept = [], [], [], []
    all_donor_names: set[str] = set()
    lost = 0
    for island in donor_islands:
        names = island_species[island] - forbidden
        counts = Counter(genus(sp) for sp in names)
        if len(names) < min_species or len(counts) < min_genera:
            lost += 1
            continue
        row = len(kept)
        kept.append(island)
        all_donor_names.update(names)
        for g, n in counts.items():
            rows.append(row)
            cols.append(vocab[g])
            vals.append(n / len(names))
    q = csr_matrix((vals, (rows, cols)),
                   shape=(len(kept), len(vocab)), dtype=float)
    if len(kept) and not np.allclose(np.asarray(q.sum(axis=1)).ravel(), 1.0):
        raise AssertionError("Filtered donor genus profiles do not sum to one")
    if all_donor_names & forbidden:
        raise AssertionError("Focal/block/first-source plant leaked into genus matrix")
    return kept, q, all_donor_names, lost


def weights_for_flower(p: np.ndarray, names: list[str]) -> dict[str, float]:
    """Keep only genuinely positive source genera; never use the target roster."""
    p = np.asarray(p, float)
    if (p.ndim != 1 or len(p) != len(names)
            or (p < 0).any() or not np.isfinite(p).all()
            or not np.isclose(p.sum(), 1.0)):
        raise ValueError("Stage-1 genus weights invalid")
    out = {names[j]: float((1-DELTA_UNKNOWN)*p[j])
           for j in np.flatnonzero(p)}
    out[UNKNOWN_GENUS] = DELTA_UNKNOWN
    return out


def _group_source_features(ids: list[str], geography: pd.DataFrame,
                           assign: pd.DataFrame) -> dict:
    geo = geography.set_index("island_id").reindex(ids)
    features = _floats(geo)
    lat = pd.to_numeric(geo.island_latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(geo.island_longitude, errors="coerce").to_numpy(float)
    valid_coords = (np.isfinite(lat) & np.isfinite(lon)
                    & (np.abs(lat) <= 90) & (np.abs(lon) <= 180))
    xyz = _xyz(np.where(valid_coords, lat, 0), np.where(valid_coords, lon, 0))
    region = geo.analysis_regime.fillna("geography_unresolved").astype(str).to_numpy()
    blocks = np.array([_source_group(v) for v in geo.spatial_block], dtype=object)
    ranked = assign.loc[assign.source_mode.eq("geo_k5")].copy()
    ranked.source_rank = pd.to_numeric(ranked.source_rank, errors="raise")
    if ranked.duplicated(["island_id", "source_rank"]).any():
        raise ValueError("Original candidate-source rankings not unique")
    first = ranked.loc[ranked.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    five = ranked.loc[ranked.source_rank.le(5)].groupby("island_id").entity_ID.agg(set).to_dict()
    first_keys = np.array([_source_group(first.get(i)) for i in ids], dtype=object)
    valid = valid_coords & np.isfinite(features).all(axis=1)
    valid &= np.array([v is not None for v in blocks])
    valid &= np.array([v is not None for v in first_keys])
    return dict(region=region, block=blocks, source=first_keys, valid=valid,
                xyz=xyz, features=features,
                topfive=[five.get(i, set()) for i in ids])


def _distance_predictions(focal: int, donor: np.ndarray, Q: csr_matrix,
                          geo: dict) -> dict[str, np.ndarray]:
    """All candidate stage-1 models use the same name-disjoint donor corpus."""
    uniform = np.asarray(Q.mean(axis=0)).ravel()
    if not np.isclose(uniform.sum(), 1.0):
        raise AssertionError("Stage-1 source genus category baseline invalid")
    z, xyz, sets = geo["features"], geo["xyz"], geo["topfive"]
    env = np.linalg.norm(z[donor] - z[focal], axis=1)
    spatial = np.linalg.norm(xyz[donor] - xyz[focal], axis=1)
    source = env + SOURCE_WEIGHT*np.array(
        [1.-_jaccard(sets[focal], sets[int(k)]) for k in donor], float)
    return dict(
        regional_genus_prior=uniform,
        geographic_nearest=_nearest_probabilities(Q, np.arange(len(donor)),
                                                  spatial, uniform, k=K),
        environment_nearest=_nearest_probabilities(Q, np.arange(len(donor)),
                                                   env, uniform, k=K),
        environment_source_nearest=_nearest_probabilities(Q, np.arange(len(donor)),
                                                          source, uniform, k=K),
    )


def run_analysis(floras: pd.DataFrame, geography: pd.DataFrame,
                 assignments: pd.DataFrame, direct_states: dict,
                 *, endpoints: tuple[str, ...] = PILOT_ENDPOINTS,
                 scopes: tuple[str, ...] = SCOPES,
                 regions: tuple[str, ...] = PRIMARY_REGIONS,
                 nboot: int = NBOOT) -> dict:
    ids = sorted(floras.island_id.unique())
    if len(ids) != 4505 or floras.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Original 4505 accepted botanical island population changed")
    if geography.island_id.nunique() != 8264:
        raise ValueError("Independent GSHHG 8264 island source geography changed")
    full = floras.groupby("island_id").accepted_species.agg(set).to_dict()
    native = floras.loc[floras.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    geo = _group_source_features(ids, geography, assignments)
    block_by = dict(zip(ids, geo["block"]))
    source_by = dict(zip(ids, geo["source"]))
    blacklists = original_full_flora_group_blacklists(full, block_by, source_by)
    original = source_endpoints(direct_states)
    if len(original) != 13 or not set(endpoints) <= original.keys():
        raise ValueError("Original 13 true Direct endpoints missing")
    indices_by_endpoint = {
        ep: {value: k for k, value in
             enumerate(sorted(set(original[ep].values()), key=repr))}
        for ep in endpoints
    }
    scope_floras = {"all_recorded": full, "original_native": native}
    all_rows = []
    for scope in scopes:
        membership = scope_floras[scope]
        # Every single candidate donor species has a global genus dictionary;
        # only strictly filtered source island names determine probabilities.
        vocab_names = sorted({genus(n) for names in membership.values() for n in names})
        vocab = {name: k for k, name in enumerate(vocab_names)}
        by_reg = {reg: set().union(*(
            membership.get(i, set())
            for i, r in zip(ids, geo["region"]) if r == reg))
            for reg in REGIONS}
        # No Direct state is needed for geographic genus source eligibility.
        eligible_stage1 = np.array([
            geo["valid"][j] and len(membership.get(i, set())) >= MIN_SOURCE_SPECIES
            and len({genus(n) for n in membership.get(i, set())}) >= MIN_SOURCE_GENERA
            for j, i in enumerate(ids)], bool)
        for reg in regions:
            regional = [j for j, i in enumerate(ids)
                        if eligible_stage1[j] and geo["region"][j] == reg]
            # The unqualified islands remain in the declared 4505 parent.
            region_targets = [j for j, i in enumerate(ids) if geo["region"][j] == reg]
            grouped: dict[tuple[str, str], list[int]] = defaultdict(list)
            for focal in regional:
                grouped[(geo["block"][focal], geo["source"][focal])].append(focal)
            source_banks = {}
            for ep in endpoints:
                known = set(original[ep])
                scope_sets = {r: names & known for r, names in by_reg.items()}
                source_banks[ep] = donor_source_region_ownership(scope_sets, reg)
            data_by_ep = {ep: [] for ep in endpoints}
            support = {ep: Counter() for ep in endpoints}
            for (block, first_source), focals in grouped.items():
                rep = focals[0]
                forbidden = source_training_ban(
                    ids[rep], full[ids[rep]], block, first_source,
                    blacklists, "block_and_first_source_out")
                if forbidden is None:
                    for ep in endpoints:
                        support[ep]["missing_block_source_identity"] += len(focals)
                    continue
                potential = [
                    j for j in regional
                    if geo["block"][j] != block and geo["source"][j] != first_source
                ]
                if len(potential) < MIN_SOURCE_ISLANDS:
                    for ep in endpoints:
                        support[ep]["too_few_original_group_donors"] += len(focals)
                    continue
                # Refit genus matrix from true permitted botanical species only.
                source_ids, Q, kept_names, failed_donors = filtered_genus_matrix(
                    [ids[j] for j in potential], membership, forbidden, vocab)
                donor_ids = np.array([ids.index(i) for i in source_ids], dtype=int)
                if (len(source_ids) < MIN_SOURCE_ISLANDS
                    or len(set(geo["block"][donor_ids])) < MIN_SOURCE_BLOCKS
                    or len(set(geo["source"][donor_ids])) < MIN_SOURCE_BLOCKS):
                    for ep in endpoints:
                        support[ep]["too_few_name_disjoint_source_donors"] += len(focals)
                    continue
                for ep in endpoints:
                    indices = indices_by_endpoint[ep]
                    local, external = source_banks[ep]
                    source_model = refit_no_oracle_regional_genus_source(
                        local, external, forbidden, original[ep], indices)
                    if source_model is None:
                        support[ep]["no_disjoint_direct_trait_source"] += len(focals)
                        continue
                    if (source_model["n_distinct_local_direct_source_species"] < 20
                        or source_model["n_distinct_external_direct_source_species"] < 20
                        or source_model["n_distinct_genus_profiles_with_external_direct_congeners"] < 2
                        or source_model["local_source_species_fraction_with_external_same_genus_donor"] < .10):
                        support[ep]["low_external_same_genus_coverage"] += len(focals)
                        continue
                    if ((local | external) - forbidden) & forbidden:
                        raise AssertionError("External Direct species leaked")
                    genus_profiles = source_model["external_conditional_genus_profiles"]
                    for focal in focals:
                        direct = sorted(membership.get(ids[focal], set()) & set(original[ep]))
                        if len(direct) < MIN_DIRECT:
                            support[ep]["fewer_than_five_Direct_carriers"] += 1
                            continue
                        preds = _distance_predictions(focal, donor_ids, Q, geo)
                        scores = score_focal_flora(
                            observed_codes=[indices[original[ep][name]] for name in direct],
                            weights_by_model={
                                m: weights_for_flower(preds[m], vocab_names) for m in REQUIRED_MODELS},
                            external_genus_profiles=genus_profiles,
                            external_blind=source_model["external_genus_blind_prior"],
                            regional_blind=source_model["regional_local_genus_blind_prior"],
                            observed_focal_genera_for_ORACLE_only=[genus(name) for name in direct],
                            min_direct=MIN_DIRECT,
                        )
                        scores["source_increment_Brier"] = (
                            scores["environment_nearest_Brier"]-
                            scores["environment_source_nearest_Brier"])
                        scores["source_increment_logloss"] = (
                            scores["environment_nearest_logloss"]-
                            scores["environment_source_nearest_logloss"])
                        scores["source_over_geo_Brier"] = (
                            scores["geographic_nearest_Brier"]-
                            scores["environment_source_nearest_Brier"])
                        scores["source_over_geo_logloss"] = (
                            scores["geographic_nearest_logloss"]-
                            scores["environment_source_nearest_logloss"])
                        scores["source_over_environment_Brier"] = scores["source_increment_Brier"]
                        scores["source_over_environment_logloss"] = scores["source_increment_logloss"]
                        # True local genus identities appear only in AFTER-PREDICTION audit.
                        truth_genus = Counter(genus(n) for n in membership.get(ids[focal], set()))
                        weights = preds["environment_source_nearest"]
                        known = np.array([truth_genus.get(g, 0) for g in vocab_names], float)
                        known /= max(known.sum(), 1.)
                        scores["true_focal_genus_species_fraction_outside_source_support"] = float(
                            known[(np.asarray(Q.sum(axis=0)).ravel() == 0)].sum())
                        scores["spatial_block"] = block
                        scores["first_source"] = first_source
                        data_by_ep[ep].append(scores)
            for ep in endpoints:
                scored = pd.DataFrame(data_by_ep[ep])
                rec = dict(region=reg, scope=scope, trait_endpoint=ep,
                           original_parent_botanical_floras=4505,
                           n_original_region_flora=len(region_targets),
                           n_stage1_original_region_geography_eligible=len(regional),
                           n_original_target_islands_scored=len(scored),
                           support_attrition=dict(support[ep]),
                           original_focal_genus_oracle_NOT_used_in_predictor=True,
                           original_full_accepted_name_group_ban_on_both_models=True)
                if not scored.empty:
                    rec["mean_target_Direct_carriers"] = float(scored.n_original_direct_plant_carriers.mean())
                    rec["unknown_genus_target_fraction_mean"] = float(
                        scored.true_focal_genus_species_fraction_outside_source_support.mean())
                    for k, grp in (("physical_block", "spatial_block"),
                                   ("first_source_group", "first_source")):
                        rec[k] = paired_cluster_summary(
                            scored, group=grp, columns=SUMMARY_GAIN_COLUMNS,
                            nboot=nboot, seed=int.from_bytes(hashlib.sha256(
                                (scope+reg+ep+grp).encode()).digest()[:8], "big"))
                    rec["support_both_independent_source_groups"] = bool(
                        rec["physical_block"]["supported_for_inference"] and
                        rec["first_source_group"]["supported_for_inference"])
                    rec["positive_vs_local_blind_both_scores_and_groups"] = bool(
                        dual_loss_dual_group_positive(
                            rec["physical_block"], rec["first_source_group"],
                            "environment_source_nearest", "regional_blind"))
                    rec["positive_vs_external_blind_both_scores_and_groups"] = bool(
                        dual_loss_dual_group_positive(
                            rec["physical_block"], rec["first_source_group"],
                            "environment_source_nearest", "external_blind"))
                else:
                    rec["support_both_independent_source_groups"] = False
                    rec["positive_vs_local_blind_both_scores_and_groups"] = False
                    rec["positive_vs_external_blind_both_scores_and_groups"] = False
                all_rows.append(rec)
    return {
        "status": "ORIGINAL_GIFT_STAGE2_NAME_DISJOINT_GEO_GENUS_FLORAL_PILOT",
        "original_parent_GIFT_botanical_floras": 4505,
        "original_distinct_Direct_endpoints_total": 13,
        "frozen_pilot_endpoints": list(endpoints),
        "original_group_training_ban": "whole_full_flora_same_10degree_block_and_first_source_NAMES",
        "training_uses_other_GIFT_islands_not_independent_mainland_ancestors": True,
        "focal_full_flora_species_names_used_for_CV_EXCLUSION_only_not_predictor_weight": True,
        "animal_pollinator_transition_not_inferred": True,
        "prediction_of_establishment_or_adaptation_not_inferred": True,
        "pilot_regions": list(regions),
        "scopes": list(scopes),
        "regional_source_heldout_results": all_rows,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--original-root", type=Path, required=True)
    ap.add_argument("--geography", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--endpoints", nargs="+", default=list(PILOT_ENDPOINTS))
    ap.add_argument("--nboot", type=int, default=NBOOT)
    args = ap.parse_args()
    paths = {
        "flora": args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assign": args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        "traits": args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
    }
    for name, path in paths.items():
        if sha(path) != SOURCE_PINS[name]:
            raise ValueError("Original GIFT source archive identity changed: "+name)
    flora = pd.read_csv(paths["flora"], dtype=str,
                        usecols=["island_id", "accepted_species", "origin_status"])
    ass = pd.read_csv(paths["assign"], dtype=str,
                      usecols=["island_id", "source_mode", "source_rank", "entity_ID"])
    geo = pd.read_csv(args.geography, dtype=str, usecols=[
        "island_id", "analysis_regime", "spatial_block", "island_latitude",
        "island_longitude", *FEATURES])
    traits = strict_traits(pd.read_csv(paths["traits"], dtype=str, usecols=[
        "accepted_species", "trait_name", "normalized_value",
        "evidence_scope", "evidence_quality"]))
    r = run_analysis(flora, geo, ass, traits,
                     endpoints=tuple(args.endpoints), nboot=args.nboot)
    raw = json.dumps(r, indent=2, sort_keys=True, allow_nan=False)
    for private in ('"island_id":', '"accepted_species":', '"entity_ID":',
                    '"spatial_block":', '"first_source":', '"genus_name":',
                    '"latitude":', '"longitude":', '"source_genus":'):
        if private in raw:
            raise AssertionError("Private original taxon/island/source identifier in report")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + "\n", encoding="utf-8")
    for row in r["regional_source_heldout_results"]:
        print("ANONYMOUS_GIFT_STAGE2_REGIONAL="+json.dumps(row, sort_keys=True))
    print("ANONYMOUS_GIFT_STAGE2_PILOT_SUMMARY="+json.dumps(dict(
        status=r["status"], parent=4505,
        panel_count=len(r["regional_source_heldout_results"]),
        n_qualified=sum(q["support_both_independent_source_groups"]
                        for q in r["regional_source_heldout_results"]),
        n_local_blind_positive=sum(q["positive_vs_local_blind_both_scores_and_groups"]
                                   for q in r["regional_source_heldout_results"])),
        sort_keys=True))


if __name__ == "__main__":
    main()
