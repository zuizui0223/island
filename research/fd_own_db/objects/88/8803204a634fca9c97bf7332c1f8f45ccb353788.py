"""Outcome-blind, geographically disjoint oceanic/continental source-flora matching.

This is a conservative exploratory sensitivity of geological-origin
ASSOCIATION, not independent geological formation histories or a randomized
island treatment. An ambiguous archipelago is not assigned a made-up ID.
Florescence colour/form/symmetry outcomes are not inspected while matching.
"""
from __future__ import annotations

from itertools import product
import numpy as np
import pandas as pd

from .geography_independent import independent_pair_matching

METRICS=(
    "source_relative_marginal_rao",
    "source_relative_joint_rao",
    "source_relative_colour_architecture_pairing_departure",
    "shared_source_genus_fraction",
    "source_genus_explains_local_fraction",
    "expected_10genus_share",
)

def _float_log(x):
    n=float(x)
    return float(np.log1p(n)) if np.isfinite(n) and n>=0 else np.nan


def _features(q:pd.DataFrame)->np.ndarray:
    cols=("area_km2","distance_to_continent_km",
          "n_original_direct_species","n_source_genera",
          "climate_pc1","climate_pc2","climate_pc3","climate_pc4")
    if not set(cols).issubset(q):
        raise ValueError("Origin-controlled original source covariates missing")
    out=np.asarray([
        [_float_log(row[c]) if i<4 else float(row[c])
         for i,c in enumerate(cols)]
        for _,row in q.iterrows()],float)
    if not np.isfinite(out).all():
        raise ValueError("Source geography or candidate opportunity covariates invalid")
    return out


def matched_oceanic_continental_source_audit(
    source_rows:pd.DataFrame,*, max_pairs_report:int=20
)->list[dict]:
    """Release ONLY source-disjoint geographic pair support and conditional signs."""
    required={"island_id","origin","regime","scope",
              "spatial_block","archipelago",*METRICS}
    if not required.issubset(source_rows):
        raise ValueError("Original geological cohort source panel lacks input schema")
    res=[]
    for scope in ("all_recorded","regional_WCVP_compatible","strict_source_native"):
        frame=source_rows.loc[source_rows.scope.eq(scope)]
        for region in ("northern_midlatitude","tropical"):
            q=frame.loc[frame.regime.eq(region)].copy()
            # Missing historical climate PC values cannot be fabricated; drop
            # only for outcome-blind matched geographic eligibility.
            covariates=("area_km2","distance_to_continent_km",
                        "n_original_direct_species","n_source_genera",
                        "climate_pc1","climate_pc2","climate_pc3","climate_pc4")
            for key in covariates:
                q[key]=pd.to_numeric(q[key],errors="coerce")
            finite=np.isfinite(q[list(covariates)].to_numpy(dtype=float)).all(axis=1)
            positive=(q.area_km2.gt(0)&q.distance_to_continent_km.ge(0)&
                      q.n_original_direct_species.ge(30)&q.n_source_genera.ge(1))
            before=len(q)
            q=q.loc[finite&positive].copy()
            c=q.loc[q.origin.eq("continental")].sort_values("island_id")
            o=q.loc[q.origin.eq("oceanic")].sort_values("island_id")
            row=dict(scope=scope,region=region,
                     n_source_eligible_islands_before_geography_complete_case=before,
                     n_source_eligible_islands_complete_geo_covariates=len(q),
                     n_eligible_continental_islands=len(c),
                     n_eligible_oceanic_islands=len(o),
                     n_candidate_pairs_before_archipelago_and_source_independence=0,
                     n_spatial_and_verified_archipelago_disjoint_pairs=0,
                     status="no_outcome_blind_geology_matched_pairs")
            if c.empty or o.empty:
                res.append(row);continue
            xc=_features(c);xo=_features(o)
            X=np.vstack((xc,xo))
            sd=np.maximum(X.std(axis=0),[.25,.25,.25,.25,.1,.1,.1,.1])
            edges=[]
            for i,j in product(range(len(c)),range(len(o))):
                if (abs(xc[i,0]-xo[j,0])>1.10 or
                    abs(xc[i,1]-xo[j,1])>1.10 or
                    abs(xc[i,2]-xo[j,2])>1.10 or
                    abs(xc[i,3]-xo[j,3])>1.10):
                    continue
                clima=(xc[i,4:]-xo[j,4:])/sd[4:]
                if np.linalg.norm(clima)>2.5:continue
                ea=c.iloc[i];eb=o.iloc[j]
                edges.append({
                    "continental_id":ea.island_id,
                    "oceanic_id":eb.island_id,
                    "region":region,
                    "c_spatial_block":ea.spatial_block,
                    "o_spatial_block":eb.spatial_block,
                    "c_archipelago":ea.archipelago,
                    "o_archipelago":eb.archipelago,
                    "cost":float(np.linalg.norm((xc[i]-xo[j])/sd)),
                })
            row["n_candidate_pairs_before_archipelago_and_source_independence"]=len(edges)
            # Fixed maximum 6 opportunity-based candidate pairs PER endpoint,
            # not trait-score-driven pruning, prevents huge MILP graphs.
            if len(edges)>500:
                selected=set()
                for key in ("continental_id","oceanic_id"):
                    grouped={}
                    for pos,e in enumerate(edges):
                        grouped.setdefault(e[key],[]).append(pos)
                    for ids in grouped.values():
                        selected.update(sorted(ids,key=lambda pos:(
                            edges[pos]["cost"],
                            edges[pos]["continental_id"],
                            edges[pos]["oceanic_id"]))[:6])
                edges=[e for k,e in enumerate(edges) if k in selected]
            row["n_edges_after_outcome_blind_pair_pruning"]=len(edges)
            if not edges:
                res.append(row);continue
            eligible=[e for e in edges if
                isinstance(e["c_archipelago"],str) and e["c_archipelago"].strip()
                and isinstance(e["o_archipelago"],str) and e["o_archipelago"].strip()]
            if not eligible:
                row["status"]="no_both_verified_archipelago_source_edges"
                res.append(row);continue
            try:
                matched,audit=independent_pair_matching(
                    eligible,"both_verified",assert_optimal=True)
            except ValueError:
                row["status"]="geography_disjoint_matching_unverified"
                res.append(row);continue
            row["n_spatial_and_verified_archipelago_disjoint_pairs"]=len(matched)
            row["n_distinct_geological_pair_spatial_blocks"]=audit[
                "n_distinct_spatial_blocks_used"]
            row["n_verified_archipelago_labels_across_selected_pairs"]=audit[
                "n_distinct_labelled_archipelagos_used"]
            row["n_unknown_archipelago_endpoints"]=audit["n_unknown_archipelago_endpoints"]
            if not matched:
                res.append(row);continue
            island=c.set_index("island_id")
            ocean=o.set_index("island_id")
            for metric in METRICS:
                diffs=np.array([
                    ocean.at[e["oceanic_id"],metric]-
                    island.at[e["continental_id"],metric]
                    for e in matched],dtype=float)
                diffs=diffs[np.isfinite(diffs)]
                label="pair_oceanic_minus_continental_"+metric
                if len(diffs):
                    row[label+"_mean"]=float(diffs.mean())
                    row[label+"_median"]=float(np.median(diffs))
                    row[label+"_n"]=len(diffs)
                if 6<=len(diffs)<=max_pairs_report and len(diffs)<=16:
                    # Exact 2^n sign assignments and two-sided no-fabricated-p.
                    configs=np.indices([2]*len(diffs)).reshape(len(diffs),-1).T
                    signs=configs*2-1
                    signed=np.abs(signs@diffs)
                    row[label+"_exact_two_sided_p"]=float(np.mean(
                        signed>=abs(sum(diffs))-1e-12))
            row["status"]=(
                "descriptive_geography_disjoint_few_pairs"
                if len(matched)<6 else
                "conditional_geography_disjoint_exact_sign_test")
            res.append(row)
    return res
