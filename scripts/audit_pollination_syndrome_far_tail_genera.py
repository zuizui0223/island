"""Post-hoc genus attribution diagnostic for the >783 km yellow/orange × deep-tube signal.

This diagnostic is explicitly outside the frozen Chapter 1 submission inference.
It identifies which plant genera contribute eligible yellow/orange deep-tube
occurrences in the far tail and refits the same hinge model after removing one
genus at a time.  Pollinator-labelled architecture names are descriptive only;
this script does not infer realized pollinators.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_v13_raw_colour_coupling_audit import _species_raw_states

# Import the exact support/model helpers used by the parent diagnostic.
from audit_pollination_syndrome_nonlinearity import _fit_models, _load_windows, _prepare_focus

COMBINATION = "yellow_orange__butterfly_deep_tube_given_colour"
CONTEXTS = ("northern_midlatitude", "tropical")


def _genus(species: str) -> str:
    text = str(species).strip().strip('"').strip()
    if not text:
        return ""
    return text.split()[0].strip('"')


def _species_occurrences(
    species_axis: pd.DataFrame,
    status_flora: pd.DataFrame,
    *,
    evidence_scope: str,
) -> pd.DataFrame:
    states = _species_raw_states(species_axis, evidence_scope=evidence_scope)
    if states.empty:
        return pd.DataFrame(columns=["island_id", "accepted_species", "genus", "success"])
    eligible = states.loc[
        states["colour_states"].map(lambda x: "yellow_orange" in x)
        & states["tube_depth_class"].map(bool)
    ].copy()
    eligible["success"] = eligible["tube_depth_class"].map(lambda x: "deep" in x).astype(int)
    eligible["genus"] = eligible["accepted_species"].map(_genus)
    flora = status_flora[["island_id", "accepted_species"]].copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.drop_duplicates(["island_id", "accepted_species"])
    out = flora.merge(
        eligible[["accepted_species", "genus", "success"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    return out


def _validate_occurrence_reconstruction(occ: pd.DataFrame, work: pd.DataFrame) -> None:
    rebuilt = (
        occ.groupby("island_id", as_index=False)
        .agg(rebuilt_trials=("success", "size"), rebuilt_successes=("success", "sum"))
    )
    check = work[["island_id", "trials", "successes"]].merge(
        rebuilt, on="island_id", how="left", validate="one_to_one"
    )
    check[["rebuilt_trials", "rebuilt_successes"]] = check[
        ["rebuilt_trials", "rebuilt_successes"]
    ].fillna(0)
    if not np.array_equal(
        check["trials"].to_numpy(int), check["rebuilt_trials"].to_numpy(int)
    ):
        raise RuntimeError("species-level reconstruction does not match island trial counts")
    if not np.array_equal(
        check["successes"].to_numpy(int), check["rebuilt_successes"].to_numpy(int)
    ):
        raise RuntimeError("species-level reconstruction does not match island success counts")


def _summarize_occurrences(
    occ: pd.DataFrame,
    model_islands: pd.DataFrame,
    *,
    context: str,
    scope: str,
    knot_km: float,
    cluster_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ids = model_islands[["island_id", "distance_to_continent_km", cluster_column]].drop_duplicates("island_id")
    z = occ.merge(ids, on="island_id", how="inner", validate="many_to_one")
    z["segment"] = np.where(z["distance_to_continent_km"].gt(knot_km), "far_tail", "pre_hinge")

    rows = []
    for segment, part in z.groupby("segment", sort=False):
        total_trials = int(len(part))
        total_success = int(part["success"].sum())
        for genus, g in part.groupby("genus", sort=False):
            trials = int(len(g))
            success = int(g["success"].sum())
            rows.append({
                "evidence_scope": scope,
                "context": context,
                "segment": segment,
                "genus": genus,
                "n_trial_occurrences": trials,
                "n_deep_occurrences": success,
                "n_trial_species": int(g["accepted_species"].nunique()),
                "n_deep_species": int(g.loc[g["success"].eq(1), "accepted_species"].nunique()),
                "n_trial_islands": int(g["island_id"].nunique()),
                "n_deep_islands": int(g.loc[g["success"].eq(1), "island_id"].nunique()),
                "within_genus_deep_rate": success / trials if trials else np.nan,
                "share_of_all_trials": trials / total_trials if total_trials else np.nan,
                "share_of_all_deep_occurrences": success / total_success if total_success else np.nan,
            })
    genus_segment = pd.DataFrame(rows)

    species = (
        z.loc[z["segment"].eq("far_tail")]
        .groupby(["genus", "accepted_species"], as_index=False)
        .agg(
            n_islands=("island_id", "nunique"),
            n_deep_occurrences=("success", "sum"),
            deep=("success", "max"),
        )
    )
    species.insert(0, "context", context)
    species.insert(0, "evidence_scope", scope)
    species = species.sort_values(["n_deep_occurrences", "n_islands"], ascending=False)

    wide = genus_segment.pivot(index="genus", columns="segment", values=["n_trial_occurrences", "n_deep_occurrences"]).fillna(0)
    if not wide.empty:
        total = z.groupby("segment").agg(trials=("success", "size"), deep=("success", "sum"))
        far_trials = float(total.loc["far_tail", "trials"]) if "far_tail" in total.index else np.nan
        pre_trials = float(total.loc["pre_hinge", "trials"]) if "pre_hinge" in total.index else np.nan
        pooled = []
        for genus in wide.index:
            gf = float(wide.loc[genus].get(("n_deep_occurrences", "far_tail"), 0.0))
            gp = float(wide.loc[genus].get(("n_deep_occurrences", "pre_hinge"), 0.0))
            pooled.append({
                "evidence_scope": scope,
                "context": context,
                "genus": genus,
                "pooled_additive_far_minus_pre": (gf / far_trials if far_trials else 0.0) - (gp / pre_trials if pre_trials else 0.0),
            })
        pooled = pd.DataFrame(pooled).sort_values("pooled_additive_far_minus_pre", ascending=False)
    else:
        pooled = pd.DataFrame()
    return genus_segment, species, pooled


def _genus_island_counts(occ: pd.DataFrame, genus: str) -> pd.DataFrame:
    g = occ.loc[occ["genus"].eq(genus)].copy()
    if g.empty:
        return pd.DataFrame(columns=["island_id", "genus_trials", "genus_successes"])
    return (
        g.groupby("island_id", as_index=False)
        .agg(genus_trials=("success", "size"), genus_successes=("success", "sum"))
    )


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--status-flora", type=Path, required=True)
    parser.add_argument("--all-counts", type=Path, required=True)
    parser.add_argument("--direct-counts", type=Path, required=True)
    parser.add_argument("--all-scores", type=Path, required=True)
    parser.add_argument("--direct-scores", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    species_axis = pd.read_csv(args.species_axis)
    status_flora = pd.read_csv(args.status_flora, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    lower_km, _iqr, knot_km = _load_windows(args.windows)
    cluster = str(cfg["cluster_column"])

    args.output.mkdir(parents=True, exist_ok=True)
    segment_tables = []
    species_tables = []
    pooled_tables = []
    logo_rows = []
    headline_rows = []

    for scope, counts_path, scores_path in [
        ("all", args.all_counts, args.all_scores),
        ("direct", args.direct_counts, args.direct_scores),
    ]:
        counts = pd.read_csv(counts_path, dtype={"island_id": str})
        scores = pd.read_csv(scores_path, dtype={"island_id": str})
        occ = _species_occurrences(species_axis, status_flora, evidence_scope=scope)

        for context in CONTEXTS:
            work = _prepare_focus(
                counts,
                scores,
                cov,
                cfg,
                context=context,
                combination=COMBINATION,
                lower_km=lower_km,
            )
            _validate_occurrence_reconstruction(occ, work)
            base = _fit_models(work, cfg, knot_km=knot_km)
            tail = work.loc[work["distance_to_continent_km"].gt(knot_km)]
            headline_rows.append({
                "evidence_scope": scope,
                "context": context,
                "n_model_islands": int(work["island_id"].nunique()),
                "n_far_tail_islands": int(tail["island_id"].nunique()),
                "n_far_tail_blocks": int(tail[cluster].nunique()),
                "baseline_post_hinge_slope": base["above_slope_per_log1p_km"],
                "baseline_post_hinge_p": base["above_p"],
                "hinge_km": knot_km,
            })

            seg, sp, pooled = _summarize_occurrences(
                occ, work, context=context, scope=scope, knot_km=knot_km, cluster_column=cluster
            )
            segment_tables.append(seg)
            species_tables.append(sp)
            pooled_tables.append(pooled)

            far = seg.loc[seg["segment"].eq("far_tail")].copy()
            candidates = far.loc[
                (far["n_deep_occurrences"].ge(2)) | (far["n_trial_occurrences"].ge(5))
            ].sort_values(["n_deep_occurrences", "n_trial_occurrences"], ascending=False)
            named = {"Isoplexis", "Canarina", "Lotus", "Clermontia"}
            genera = list(dict.fromkeys(candidates["genus"].tolist() + [g for g in named if g in set(far["genus"])]))

            for genus in genera:
                gi = _genus_island_counts(occ, genus)
                adj = work.merge(gi, on="island_id", how="left", validate="one_to_one")
                adj[["genus_trials", "genus_successes"]] = adj[["genus_trials", "genus_successes"]].fillna(0)
                adj["trials"] = adj["trials"] - adj["genus_trials"]
                adj["successes"] = adj["successes"] - adj["genus_successes"]
                adj = adj.loc[adj["trials"].gt(0)].copy()
                if adj["island_id"].nunique() < 50:
                    continue
                fit = _fit_models(adj, cfg, knot_km=knot_km)
                logo_rows.append({
                    "evidence_scope": scope,
                    "context": context,
                    "genus_removed": genus,
                    "n_islands_after_removal": int(adj["island_id"].nunique()),
                    "n_far_tail_islands_after_removal": fit["n_above_knot"],
                    "post_hinge_slope_after_removal": fit["above_slope_per_log1p_km"],
                    "post_hinge_p_after_removal": fit["above_p"],
                    "baseline_post_hinge_slope": base["above_slope_per_log1p_km"],
                    "delta_post_hinge_slope_removed_minus_baseline": fit["above_slope_per_log1p_km"] - base["above_slope_per_log1p_km"],
                    "hinge_change_after_removal": fit["hinge_change"],
                    "hinge_change_p_after_removal": fit["hinge_change_p"],
                })

    segment = pd.concat(segment_tables, ignore_index=True) if segment_tables else pd.DataFrame()
    species = pd.concat(species_tables, ignore_index=True) if species_tables else pd.DataFrame()
    pooled = pd.concat(pooled_tables, ignore_index=True) if pooled_tables else pd.DataFrame()
    logo = pd.DataFrame(logo_rows)
    headline = pd.DataFrame(headline_rows)

    segment.to_csv(args.output / "far_tail_genus_segment_counts.csv", index=False)
    species.to_csv(args.output / "far_tail_species_counts.csv", index=False)
    pooled.to_csv(args.output / "far_tail_genus_pooled_contribution.csv", index=False)
    logo.to_csv(args.output / "far_tail_genus_leave_one_out.csv", index=False)
    headline.to_csv(args.output / "far_tail_genus_headline.csv", index=False)

    direct_far = segment.loc[
        segment["evidence_scope"].eq("direct") & segment["segment"].eq("far_tail")
    ].copy()
    direct_far = direct_far.sort_values(
        ["context", "n_deep_occurrences", "n_trial_occurrences"],
        ascending=[True, False, False],
    )
    direct_far.to_csv(args.output / "far_tail_genus_ranked_direct.csv", index=False)

    print(headline.to_string(index=False))
    for context in CONTEXTS:
        top = direct_far.loc[direct_far["context"].eq(context)].head(20)
        print(f"\n=== direct top genera: {context} ===")
        print(top[["genus", "n_deep_occurrences", "n_trial_occurrences", "n_deep_islands", "n_deep_species", "share_of_all_deep_occurrences"]].to_string(index=False))
        if not logo.empty:
            l = logo.loc[(logo["evidence_scope"].eq("direct")) & (logo["context"].eq(context))].copy()
            l = l.sort_values("delta_post_hinge_slope_removed_minus_baseline")
            print("\nLeave-one-genus-out (most slope-reducing removals):")
            print(l.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
