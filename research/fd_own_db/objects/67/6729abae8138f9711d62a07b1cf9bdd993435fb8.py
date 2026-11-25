"""Geological-island genus opportunity rankings from potential native mainland sources.

Reads authorized identity-level GIFT mainland lists, original island floras,
source assignments, and independently audited geology locally.  Writes ONLY
aggregate rank metrics, audits and cryptographic input digests.

No unrecorded genus is assumed truly absent. This models geographic source
OPPORTUNITY for a recorded genus among a source-defined candidate set. It
does NOT model actual arrival, persistence or colonization success, and its
AUC is only recorded-vs-unrecorded label discrimination.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .g1_fixed_margins import _bh
from .repertoire import REGIONS

MODELS = ("global", "local_count", "local_distance", "local_density")
STATUSES = ("all_recorded", "explicit_native")
METRICS = ("average_precision", "roc_auc")
PAIR_MODELS = ("local_count", "local_distance", "local_density")
SOURCE_MODE = "geo_k5"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1 << 20), b""):
            digest.update(part)
    return digest.hexdigest()


def ranking_metrics(positives: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    """Tie-aware AP and Mann–Whitney AUC on recorded/unrecorded candidates."""
    y = np.asarray(positives, dtype=bool)
    s = np.asarray(scores, dtype=float)
    if y.ndim != 1 or s.shape != y.shape or len(y) < 2 or not np.isfinite(s).all():
        raise ValueError("Invalid ranked candidate labels/scores")
    pos = int(y.sum())
    neg = len(y) - pos
    if pos == 0 or neg == 0:
        raise ValueError("Must have both recorded and unrecorded candidate genera")
    order = np.argsort(-s, kind="stable")
    sy, ss = y[order], s[order]
    tp_cum = np.cumsum(sy)
    breaks = np.r_[np.flatnonzero(ss[1:] != ss[:-1]), len(ss) - 1]
    old = np.r_[0, tp_cum[breaks[:-1]]]
    increment = (tp_cum[breaks] - old) / pos
    ap = float(np.sum(increment * (tp_cum[breaks] / (breaks + 1))))
    ranks = rankdata(s, method="average")
    auc = float((ranks[y].sum() - pos * (pos + 1) / 2) / (pos * neg))
    return {"average_precision": ap, "roc_auc": auc}


def exclude_heldout_species(mainland: pd.DataFrame, names: set[str]) -> pd.DataFrame:
    """Drop target taxon names from EVERY mainland entity before feature fitting."""
    if not {"work_species", "entity_ID"}.issubset(mainland.columns):
        raise ValueError("Source records must include work_species and entity_ID")
    cleaned = mainland.loc[~mainland.work_species.isin(names)].copy()
    if set(cleaned.work_species) & names:
        raise AssertionError("Exact test taxon leaked into GIFT source features")
    if cleaned.empty:
        raise ValueError("Name holdout exhausted the source training sample")
    return cleaned


def _native_mainland_source_data(native: pd.DataFrame):
    required = {"entity_ID", "work_species"}
    if not required.issubset(native.columns):
        raise ValueError("Missing source native mainland columns")
    native = native[["entity_ID", "work_species"]].drop_duplicates().copy()
    if native.isna().any().any() or native.empty:
        raise ValueError("Missing native GIFT entity/name fields")
    native["genus"] = native.work_species.str.split().str[0]
    by_entity = native.groupby(["entity_ID", "genus"]).size().rename("n_taxa").reset_index()
    source_map = {k: p.set_index("genus").n_taxa.to_dict()
                  for k, p in by_entity.groupby("entity_ID", sort=False)}
    sizes = native.groupby("entity_ID").work_species.nunique().to_dict()
    global_genus = native.groupby("genus").entity_ID.nunique().to_dict()
    return source_map, sizes, global_genus, native


def source_candidate_scores(assignment: list[tuple[str, float]],
                            source_map: dict, source_sizes: dict,
                            global_genus: dict, n_source_regions: int
                            ) -> tuple[list[str], dict[str, np.ndarray]]:
    """Four *mainland-only* scores on a common source-eligible genus cohort."""
    if n_source_regions < 1:
        raise ValueError("At least one native source region required")
    rows = [(str(ent), float(dist)) for ent, dist in assignment
            if ent in source_map and np.isfinite(dist) and dist >= 0]
    rows.sort(key=lambda t: (t[1], t[0]))
    if len(rows) < 3:
        raise ValueError("Insufficient nearby source entities")
    candidates = sorted(set().union(*(set(source_map[ent]) for ent, _ in rows)))
    if len(candidates) < 5:
        raise ValueError("Insufficient source-eligible genera")
    weights = np.asarray([1.0 / np.sqrt(1.0 + dist / 500.0) for _, dist in rows])
    weights /= weights.sum()
    global_score = np.array([global_genus.get(g, 0) / n_source_regions
                             for g in candidates], dtype=float)
    present = np.array([sum(g in source_map[ent] for ent, _ in rows)
                        for g in candidates], dtype=float)
    local_distance = np.array([
        sum(w for (ent, _), w in zip(rows, weights) if g in source_map[ent])
        for g in candidates
    ], dtype=float)
    local_density = np.array([
        sum(w * source_map[ent].get(g, 0) / max(source_sizes[ent], 1)
            for (ent, _), w in zip(rows, weights))
        for g in candidates
    ], dtype=float)
    return candidates, {
        "global": global_score,
        "local_count": present,
        "local_distance": local_distance,
        "local_density": local_density,
    }


def _paired_uncertainty(values: np.ndarray, *, seed: int = 20261108,
                        draws: int = 1999) -> tuple[float, float, float]:
    if draws < 99 or len(values) < 3:
        raise ValueError("Insufficient matched clusters/draws")
    rng = np.random.default_rng(seed)
    selected = rng.integers(len(values), size=(draws, len(values)))
    ci = np.quantile(values[selected].mean(axis=1), [0.025, 0.975])
    signs = rng.choice([-1.0, 1.0], size=(draws, len(values)))
    obs = float(values.mean())
    ref = (signs * values).mean(axis=1)
    p = (1 + int((abs(ref) >= abs(obs) - 1e-14).sum())) / (draws + 1)
    return float(ci[0]), float(ci[1]), float(p)


def run(flora_path: Path, native_mainland_path: Path, assignment_path: Path,
        geological_crosswalk_path: Path, output: Path, *,
        matched_pairs_path: Path | None = None,
        source_mode: str = SOURCE_MODE, seed: int = 20261108,
        holdout_observed_species: bool = True) -> dict:
    flora_path, native_mainland_path, assignment_path, geological_crosswalk_path, output = (
        map(Path, (flora_path, native_mainland_path, assignment_path,
                   geological_crosswalk_path, output)))
    output.mkdir(parents=True, exist_ok=True)
    geography = pd.read_csv(geological_crosswalk_path, dtype=str)
    geography = geography.loc[
        geography.geological_origin.isin(("oceanic", "continental"))]
    geo = {r.island_id: r for r in geography.itertuples(index=False)}
    flora = pd.read_csv(flora_path,
                        usecols=["island_id", "accepted_species", "origin_status"],
                        dtype=str)
    flora = flora.loc[flora.island_id.isin(geo)].copy()
    if flora.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Duplicate island/taxon source incidence")
    flora["genus"] = flora.accepted_species.str.split().str[0]
    all_observed = flora.groupby("island_id").genus.agg(lambda x: set(x)).to_dict()
    native_observed = flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").genus.agg(lambda x: set(x)).to_dict()
    mainland = pd.read_csv(native_mainland_path,
                           usecols=["entity_ID", "work_species"], dtype=str)
    raw_mainland_rows = len(mainland)
    source_overlap_count = int(mainland.work_species.isin(
        set(flora.accepted_species)).sum())
    if holdout_observed_species:
        # Entire geological reference population, not individual target
        # island held out in isolation. Never exploit identical native
        # mainland species as a shortcut to recognizing its island record.
        mainland = exclude_heldout_species(
            mainland, set(flora.accepted_species))
    source_map, source_sizes, global_genus, native = _native_mainland_source_data(mainland)
    assignments = pd.read_csv(
        assignment_path, dtype=str,
        usecols=["island_id", "entity_ID", "source_mode", "source_distance_km"],
    )
    assignments = assignments.loc[
        assignments.island_id.isin(geo) &
        assignments.source_mode.eq(source_mode)]
    records = []
    for island, part in assignments.groupby("island_id", sort=False):
        if island not in all_observed:
            continue
        selected = []
        for r in part.itertuples(index=False):
            try:
                distance = float(r.source_distance_km)
            except (TypeError, ValueError):
                continue
            selected.append((r.entity_ID, distance))
        try:
            candidates, scoring = source_candidate_scores(
                selected, source_map, source_sizes, global_genus,
                native.entity_ID.nunique())
        except ValueError:
            continue
        row = geo[island]
        for status, observed in (
            ("all_recorded", all_observed.get(island, set())),
            ("explicit_native", native_observed.get(island, set())),
        ):
            if len(observed) < 10:
                continue
            positive = np.asarray([g in observed for g in candidates], dtype=bool)
            if sum(positive) < 5 or sum(positive) == len(positive):
                continue
            baseline_prevalence = float(positive.mean())
            coverage = float(positive.sum() / len(observed))
            for model, score in scoring.items():
                metrics = ranking_metrics(positive, score)
                records.append({
                    "island_id": island, "geological_origin": row.geological_origin,
                    "analysis_regime": row.analysis_regime,
                    "spatial_block": row.spatial_block, "status": status,
                    "model": model, **metrics,
                    "AP_lift_over_prevalence": metrics["average_precision"] -
                                               baseline_prevalence,
                    "n_candidate_genera": len(candidates),
                    "n_observed_genera": len(observed),
                    "n_observed_candidate_genera": int(positive.sum()),
                    "source_genus_coverage": coverage,
                    "candidate_positive_fraction": baseline_prevalence,
                })
    if not records:
        raise ValueError("No evaluable island source/observed genus combinations")
    frame = pd.DataFrame(records)
    summary = frame.groupby(
        ["status", "geological_origin", "model"], as_index=False
    ).agg(
        n_islands=("island_id", "nunique"),
        mean_AP=("average_precision", "mean"),
        median_AP=("average_precision", "median"),
        mean_auc=("roc_auc", "mean"),
        median_auc=("roc_auc", "median"),
        mean_AP_lift=("AP_lift_over_prevalence", "mean"),
        median_genus_coverage=("source_genus_coverage", "median"),
        median_candidate_genera=("n_candidate_genera", "median"),
        median_recorded_genera=("n_observed_genera", "median"),
    )
    summary.to_csv(output / "source_genus_ranking_summary.csv", index=False)
    original = {}
    for (status, origin), members in frame.groupby(
        ["status", "geological_origin"]):
        for endpoint in METRICS:
            piv = members.pivot(index="island_id", columns="model", values=endpoint)
            for m in PAIR_MODELS:
                delta = (piv[m] - piv["global"]).to_numpy()
                if len(delta) < 10:
                    continue
                lo, hi, _ = _paired_uncertainty(delta, draws=1999, seed=seed)
                original[(status, origin, endpoint, m)] = {
                    "status": status, "origin": origin, "endpoint": endpoint,
                    "comparison": m + " minus global", "n_islands": len(delta),
                    "mean_ranking_gain": float(delta.mean()),
                    "island_bootstrap_q025": lo, "island_bootstrap_q975": hi,
                }
    pd.DataFrame(list(original.values())).to_csv(
        output / "source_ranking_improvements.csv", index=False)
    paired_tests = []
    if matched_pairs_path is not None:
        matched_pairs_path = Path(matched_pairs_path)
        pairs = pd.read_csv(matched_pairs_path, dtype=str)
        if not {"continental_id", "oceanic_id", "region"}.issubset(pairs.columns):
            raise ValueError("Matched-pair private input is missing required columns")
        for status in STATUSES:
            for endpoint in METRICS:
                pivot = frame.loc[frame.status.eq(status)].pivot(
                    index="island_id", columns="model", values=endpoint)
                for model in PAIR_MODELS:
                    data = []
                    for p in pairs.itertuples(index=False):
                        if p.continental_id not in pivot.index or p.oceanic_id not in pivot.index:
                            continue
                        gain_o = pivot.at[p.oceanic_id, model] - pivot.at[p.oceanic_id, "global"]
                        gain_c = pivot.at[p.continental_id, model] - pivot.at[p.continental_id, "global"]
                        data.append(float(gain_o - gain_c))
                    if len(data) < 8:
                        continue
                    data = np.asarray(data)
                    lo, hi, pval = _paired_uncertainty(data, draws=9999, seed=seed)
                    paired_tests.append({
                        "status": status, "endpoint": endpoint,
                        "model": model, "n_pairs": len(data),
                        "mean_oceanic_minus_continental_gain": float(data.mean()),
                        "paired_bootstrap_q025": lo, "paired_bootstrap_q975": hi,
                        "p_pair_swap_two_sided": pval,
                    })
        if len(paired_tests) != 12:
            raise ValueError("Cannot form complete predeclared 12-test pair contrast family")
        q = _bh([r["p_pair_swap_two_sided"] for r in paired_tests])
        harmonic = float(sum(1 / k for k in range(1, len(q) + 1)))
        for i, row in enumerate(paired_tests):
            row["q_bh_12"] = float(q[i])
            row["q_by_12"] = float(min(1., q[i] * harmonic))
        pd.DataFrame(paired_tests).to_csv(
            output / "geological_source_gain_paired_contrasts.csv", index=False)
    summary_native = frame.loc[
        frame.model.eq("global")].groupby(
        ["status", "geological_origin"], as_index=False
    ).agg(
        n_islands=("island_id", "nunique"),
        median_source_genus_coverage=("source_genus_coverage", "median"),
        median_n_observed_genera=("n_observed_genera", "median"),
        median_n_candidates=("n_candidate_genera", "median"),
    )
    summary_native.to_csv(output / "source_support_by_status.csv", index=False)
    manifest = {
        "schema": "fd_gift_genus_source_opportunity_ranking_v1",
        "status": "source_only_recorded_genus_ranking_not_colonization",
        "source_mode": source_mode,
        "source_species_holdout": holdout_observed_species,
        "n_GIFT_entity_name_records_original": raw_mainland_rows,
        "n_GIFT_entity_name_records_with_test_taxa": source_overlap_count,
        "n_GIFT_entity_name_records_retained": len(mainland),
        "n_islands_evaluated": int(frame.island_id.nunique()),
        "n_geology_classified": len(geo),
        "n_GIFT_native_source_regions": int(native.entity_ID.nunique()),
        "n_GIFT_native_donor_genera": int(native.genus.nunique()),
        "n_original_island_flora_rows_geocoded": len(flora),
        "ranking_reference": "within each island, genera with >=1 native mainland source-region record",
        "models": {
            "global": "share of mainland regions with genus",
            "local_count": "count among nearest source regions with genus",
            "local_distance": "distance-weighted source-region presence",
            "local_density": "distance-weighted genus taxa/source-list taxon richness",
        },
        "candidate_negative_policy":
            "unrecorded genus is unlabeled, not confirmed biological absence",
        "prediction_type": "source-only zero-fit mainland opportunity ranking of island recorded genera",
        "not_predicted": ["establishment", "genuine absence", "colonization probability",
                          "dispersal trait adaptation", "pollination function"],
        "matched_pair_family": len(paired_tests),
        "input_hashes_sha256": {
            "island_flora": _sha(flora_path),
            "GIFT_native_mainland": _sha(native_mainland_path),
            "source_assignments": _sha(assignment_path),
            "geology": _sha(geological_crosswalk_path),
            **({"matched_pairs": _sha(matched_pairs_path)}
               if matched_pairs_path is not None else {}),
        },
        "privacy": "Only aggregate results and hash manifests released. No island IDs, "
                   "native GIFT entity IDs, species names or genus keys emitted.",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "status": manifest["status"], "n_islands": manifest["n_islands_evaluated"],
        "mean_AUC": [
            {"status": r.status, "geological_origin": r.geological_origin,
             "model": r.model, "n": int(r.n_islands), "AUC": float(r.mean_auc)}
            for r in summary.itertuples(index=False)
            if r.model in ("global", "local_density")
        ],
        "matched_pair_BH_significant": sum(
            r.get("q_bh_12", 1.0) < .05 for r in paired_tests),
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("flora", "mainland-native", "source-assignments",
                 "geology", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--matched-pairs", type=Path)
    parser.add_argument("--source-mode", default=SOURCE_MODE)
    parser.add_argument("--allow-overlapping-source-species", action="store_true",
        help="Descriptive source opportunity only; NOT independent predictive validation")
    args = parser.parse_args()
    run(args.flora, args.mainland_native, args.source_assignments,
        args.geology, args.output, matched_pairs_path=args.matched_pairs,
        source_mode=args.source_mode,
        holdout_observed_species=not args.allow_overlapping_source_species)


if __name__ == "__main__":
    main()
