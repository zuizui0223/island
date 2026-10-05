"""Audit genus/species contributors to the post-hoc far-tail yellow/orange x deep-tube signal.

This is a descriptive diagnostic only. It reuses the exact complete-case island set
from the hinge analysis and decomposes grouped-binomial successes into species/genus
island cells. Pollinator labels are not inferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


COMBINATION = "yellow_orange__butterfly_deep_tube_given_colour"
BASELINE = ["log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
CONTEXTS = ["northern_midlatitude", "tropical"]


def parse_traits(value):
    out = {}
    if pd.isna(value):
        return out
    for part in str(value).split("|"):
        if "=" not in part:
            continue
        key, raw = part.split("=", 1)
        try:
            vals = json.loads(raw)
        except json.JSONDecodeError:
            continue
        if not isinstance(vals, list):
            vals = [vals]
        out.setdefault(key, set()).update(str(x) for x in vals)
    return out


def species_states(species_axis, scope):
    work = species_axis.loc[
        species_axis["axis"].astype(str).isin(["flower_colour", "floral_structural_complexity"])
    ].copy()
    if scope == "direct":
        work = work.loc[
            work["quality"].fillna("").astype(str).str.lower().isin(["high", "medium"])
        ].copy()
    rows = []
    for species, group in work.groupby("accepted_species", sort=False):
        colours, tubes = set(), set()
        for row in group.itertuples(index=False):
            traits = parse_traits(row.trait_composition)
            if str(row.axis) == "flower_colour":
                colours |= traits.get("flower_primary_color", set())
            else:
                tubes |= traits.get("tube_depth_class", set())
        if "yellow_orange" not in colours or not tubes:
            continue
        rows.append(
            {
                "accepted_species": str(species),
                "tube_states": "|".join(sorted(tubes)),
                "deep": "deep" in tubes,
            }
        )
    return pd.DataFrame(rows)


def selfing_core(scores):
    return (
        scores.loc[
            scores["stratum"].astype(str).eq("all_observed")
            & scores["syndrome"].astype(str).eq("selfing_core"),
            ["island_id", "syndrome_score"],
        ]
        .rename(columns={"syndrome_score": "selfing_core"})
        .drop_duplicates("island_id")
    )


def load_knots(path):
    table = pd.read_csv(path)
    common = table.loc[table["window"].eq("common_05_95")].iloc[0]
    return float(common["lower_km"]), float(common["upper_km"])


def analytic_tail(counts, scores, cov, context, lower_km, knot_km):
    base = counts.loc[
        counts["stratum"].astype(str).eq("all_observed")
        & counts["combination"].astype(str).eq(COMBINATION)
    ].copy()
    needed = ["island_id", "distance_to_continent_km", "analysis_regime", "spatial_block", *BASELINE]
    work = (
        base.merge(
            cov[needed].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="many_to_one",
        )
        .merge(selfing_core(scores), on="island_id", how="left", validate="many_to_one")
    )
    for col in ["successes", "trials", "distance_to_continent_km", "selfing_core", *BASELINE]:
        work[col] = pd.to_numeric(work[col], errors="coerce")
    work["analysis_regime"] = work["analysis_regime"].fillna("").astype(str)
    work["spatial_block"] = work["spatial_block"].fillna("").astype(str)
    work = work.loc[work["analysis_regime"].eq(context)].dropna(
        subset=[
            "successes",
            "trials",
            "distance_to_continent_km",
            "selfing_core",
            *BASELINE,
            "spatial_block",
        ]
    )
    work = work.loc[
        work["trials"].gt(0)
        & work["distance_to_continent_km"].ge(lower_km)
        & work["spatial_block"].ne("")
    ].copy()
    return work.loc[work["distance_to_continent_km"].gt(knot_km)].copy()


def genus_name(species):
    clean = str(species).strip().strip('"')
    return clean.split()[0] if clean else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-axis", type=Path, required=True)
    ap.add_argument("--status-flora", type=Path, required=True)
    ap.add_argument("--all-counts", type=Path, required=True)
    ap.add_argument("--direct-counts", type=Path, required=True)
    ap.add_argument("--all-scores", type=Path, required=True)
    ap.add_argument("--direct-scores", type=Path, required=True)
    ap.add_argument("--covariates", type=Path, required=True)
    ap.add_argument("--windows", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    axis = pd.read_csv(args.species_axis)
    flora = pd.read_csv(args.status_flora, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    lower_km, knot_km = load_knots(args.windows)

    all_rows = []
    species_rows = []
    status_rows = []
    validation_rows = []
    target_rows = []

    for scope, count_path, score_path in [
        ("all", args.all_counts, args.all_scores),
        ("direct", args.direct_counts, args.direct_scores),
    ]:
        counts = pd.read_csv(count_path, dtype={"island_id": str})
        scores = pd.read_csv(score_path, dtype={"island_id": str})
        states = species_states(axis, scope)
        joined = flora.merge(states, on="accepted_species", how="inner", validate="many_to_one")
        joined["genus"] = joined["accepted_species"].map(genus_name)

        for context in CONTEXTS:
            tail = analytic_tail(counts, scores, cov, context, lower_km, knot_km)
            ids = set(tail["island_id"].astype(str))
            cells = joined.loc[joined["island_id"].astype(str).isin(ids)].drop_duplicates(
                ["island_id", "accepted_species"]
            )

            calc = cells.groupby("island_id").agg(
                calc_trials=("accepted_species", "nunique"),
                calc_successes=("deep", "sum"),
            )
            chk = tail.set_index("island_id")[["trials", "successes"]].join(calc, how="left").fillna(0)
            chk["trial_match"] = chk["trials"].astype(int).eq(chk["calc_trials"].astype(int))
            chk["success_match"] = chk["successes"].astype(int).eq(chk["calc_successes"].astype(int))
            if not chk["trial_match"].all() or not chk["success_match"].all():
                bad = chk.loc[~(chk["trial_match"] & chk["success_match"])]
                raise RuntimeError(f"{scope} {context}: decomposition mismatch on {len(bad)} islands")

            successes = cells.loc[cells["deep"]].copy()
            total_success = len(successes)
            genus = (
                successes.groupby("genus", dropna=False)
                .agg(
                    success_cells=("accepted_species", "size"),
                    n_islands=("island_id", "nunique"),
                    n_species=("accepted_species", "nunique"),
                )
                .reset_index()
                .sort_values(
                    ["success_cells", "n_islands", "genus"],
                    ascending=[False, False, True],
                )
            )
            genus["success_cell_share"] = genus["success_cells"] / total_success if total_success else 0.0
            genus["evidence_scope"] = scope
            genus["context"] = context
            genus["tail_islands"] = len(ids)
            genus["tail_blocks"] = tail["spatial_block"].nunique()
            genus["knot_km"] = knot_km
            all_rows.extend(genus.to_dict("records"))

            sp = (
                successes.groupby(["genus", "accepted_species"], dropna=False)
                .agg(success_cells=("island_id", "size"), n_islands=("island_id", "nunique"))
                .reset_index()
                .sort_values(["success_cells", "accepted_species"], ascending=[False, True])
            )
            sp["evidence_scope"] = scope
            sp["context"] = context
            species_rows.extend(sp.to_dict("records"))

            st = (
                successes.groupby("floristic_status", dropna=False)
                .size()
                .rename("success_cells")
                .reset_index()
                .sort_values("success_cells", ascending=False)
            )
            st["share"] = st["success_cells"] / total_success if total_success else 0.0
            st["evidence_scope"] = scope
            st["context"] = context
            status_rows.extend(st.to_dict("records"))

            coords = cov.loc[
                cov["island_id"].isin(ids),
                ["island_id", "distance_to_continent_km", "spatial_block"],
            ]
            v = tail[["island_id", "trials", "successes"]].merge(coords, on="island_id", how="left")
            v["evidence_scope"] = scope
            v["context"] = context
            validation_rows.extend(v.to_dict("records"))

            targets = ["Isoplexis", "Canarina", "Lotus", "Clermontia", "Cyanea", "Lobelia"]
            for g in targets:
                part = successes.loc[successes["genus"].eq(g)]
                target_rows.append(
                    {
                        "evidence_scope": scope,
                        "context": context,
                        "genus": g,
                        "success_cells": len(part),
                        "n_islands": part["island_id"].nunique(),
                        "n_species": part["accepted_species"].nunique(),
                        "species": "; ".join(sorted(part["accepted_species"].unique())),
                    }
                )

    pd.DataFrame(all_rows).to_csv(args.output / "far_tail_genus_contributions.csv", index=False)
    pd.DataFrame(species_rows).to_csv(args.output / "far_tail_species_contributions.csv", index=False)
    pd.DataFrame(status_rows).to_csv(args.output / "far_tail_status_contributions.csv", index=False)
    pd.DataFrame(validation_rows).to_csv(args.output / "far_tail_island_validation.csv", index=False)
    pd.DataFrame(target_rows).to_csv(args.output / "far_tail_named_genera_check.csv", index=False)

    genus_df = pd.DataFrame(all_rows)
    lines = [
        "# Far-tail yellow/orange x deep-tube genus audit",
        "",
        f"Common-support far-tail knot: {knot_km:.6f} km",
        "Descriptive post-hoc decomposition only; labels do not identify realized pollinators.",
        "",
    ]
    for scope in ["direct", "all"]:
        for context in CONTEXTS:
            part = genus_df.loc[
                genus_df["evidence_scope"].eq(scope)
                & genus_df["context"].eq(context)
            ].head(20)
            lines += [
                f"## {scope} - {context}",
                "",
                part[
                    [
                        "genus",
                        "success_cells",
                        "success_cell_share",
                        "n_islands",
                        "n_species",
                        "tail_islands",
                        "tail_blocks",
                    ]
                ].to_string(index=False),
                "",
            ]
    (args.output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
