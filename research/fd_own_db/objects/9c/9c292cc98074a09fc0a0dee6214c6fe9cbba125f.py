"""Conservative block-dependent origin label swap on matched floral repertoires.

A block-sharing graph connects two pairs if ANY of their four endpoints share
the same frozen geographic spatial_block. Under the most conservative
exchangeability assumption, all matched pairs in a connected component
receive one common sign reversal. Enumerate ALL 2^K component signs.

This is an intentionally severe spatial-dependence stress test. Block-grid
membership alone does not prove perfect statistical dependence, and sign
flips of geological origin remain observational exchangeability assumptions.
It cannot establish adaptation, lineage-independent selection or true arrival.

Only AGGREGATE statistics, component sizes and source hashes are released.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .pipeline import REF_INPUT, file_sha256, _write_csv
from .repertoire import match_origin_islands
from .effort_audit import read_effort, eligible_pairs


ROWS = (
    "scope", "subset", "dependency_grouping", "endpoint", "n_pairs", "n_components",
    "largest_component_pairs", "component_size_distribution",
    "n_exact_component_sign_configurations", "observed",
    "null_median", "null_q025", "null_q975", "exact_component_upper_p",
    "minimum_attainable_p_due_sign_symmetry",
    "status",
)
DROP_COLUMNS = (
    "scope", "subset", "dependency_grouping", "endpoint", "n_pairs_original",
    "n_components_original", "n_pairs_remaining", "n_components_remaining",
    "component_size_removed", "remaining_statistic",
    "ratio_norm_remaining_to_original", "cosine_with_original_mean_direction",
)


def connected_block_components(pairs: list[dict], geo: dict[str, dict],
                               *, field: str = "spatial_block"
                               ) -> list[list[int]]:
    """Return stable pair-index components under any shared geographic block.

    No independence is asserted across geographically neighboring but
    *different* block IDs; this is a conservative stress test, not a solved
    spatial covariance model.
    """
    if field not in ("spatial_block", "weigelt_archipelago"):
        raise ValueError("Invalid spatial dependency grouping")
    parent=list(range(len(pairs)))
    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index=parent[index]
        return index
    def union(left,right):
        parent[root(left)]=root(right)
    encountered={}
    for i,pair in enumerate(pairs):
        for label in ("continental_id","oceanic_id"):
            key=pair[label]
            if key not in geo:
                raise ValueError("Missing geographic identity of a matched island")
            block=geo[key].get(field,"")
            unavailable = (
                not isinstance(block, str) or
                block.strip().lower() in
                {"", "na", "n/a", "nan", "null", "none", "unknown", "unresolved", "-"}
            )
            if unavailable:
                if field == "spatial_block":
                    raise ValueError("Missing spatial block of a matched island")
                # Source table uses literal "NA" for an unknown archipelago.
                # It MUST NOT connect several otherwise unrelated pairs.
                block=f"__unlabeled_archipelago_{i}_{label}"
            if block in encountered:
                union(i,encountered[block])
            else:
                encountered[block]=i
    clusters={}
    for i in range(len(pairs)):
        clusters.setdefault(root(i),[]).append(i)
    return sorted((sorted(x) for x in clusters.values()),
                  key=lambda x:(-len(x),x[0]))


def _delta(pairs:list[dict], carriers:dict, space:ProfileSpace) -> np.ndarray:
    cache={}
    def vector(island):
        if island not in carriers:
            raise ValueError("Missing island profile")
        if island not in cache:
            c=carriers[island]
            if not c:raise ValueError("Island without profiles")
            v=np.zeros(len(space.ids),dtype=np.float64)
            v[[space.lookup[k] for k in c]]=1/len(c)
            cache[island]=v
        return cache[island]
    return np.stack([vector(p["oceanic_id"])-vector(p["continental_id"])
                     for p in pairs])


def _quad(mat: np.ndarray, kernel: np.ndarray) -> float:
    mean=mat.mean(axis=0)
    return float(max(0.0, mean @ kernel @ mean))


def cluster_permutation(
    pair_delta:np.ndarray, groups:list[list[int]], kernel:np.ndarray,
) ->dict:
    """Enumerate all groupwise flips; no Monte Carlo +1 correction needed."""
    n=len(pair_delta)
    if n < 8 or len(groups) < 1:
        raise ValueError("Too few pair observations")
    if sorted(i for group in groups for i in group) != list(range(n)):
        raise ValueError("Component IDs do not partition the whole cohort")
    if len(groups)>18:
        raise ValueError("More than 2^18 exhaustive null sign configurations")
    obs=_quad(pair_delta,kernel)
    null=[]
    for signs in itertools.product((-1.0,1.0),repeat=len(groups)):
        oriented=pair_delta.copy()
        for group,sign in zip(groups,signs):
            oriented[group]*=sign
        null.append(_quad(oriented,kernel))
    values=np.asarray(null)
    tolerance=1e-12*max(1.0,abs(obs))
    exact=float(np.mean(values >= obs-tolerance))
    symmetry_floor=2/(2**len(groups)) if len(groups)>=1 else 1.
    return {
        "observed":obs,"null_median":float(np.median(values)),
        "null_q025":float(np.quantile(values,.025)),
        "null_q975":float(np.quantile(values,.975)),
        "exact_component_upper_p":exact,
        "n_exact_component_sign_configurations":len(values),
        "minimum_attainable_p_due_sign_symmetry":min(1.,symmetry_floor),
    }


def delete_component_sensitivity(
    pair_delta:np.ndarray, groups:list[list[int]], kernel:np.ndarray,
    *, scope:str,subset:str,endpoint:str,
    dependency_grouping:str="spatial_block",
) -> list[dict]:
    """What if one spatially connected group drives most of the profile shift?"""
    center=pair_delta.mean(axis=0)
    q_all=_quad(pair_delta,kernel)
    records=[]
    for group in groups:
        mask=np.ones(len(pair_delta),dtype=bool)
        mask[group]=False
        if mask.sum()<8:
            continue
        remaining=pair_delta[mask]
        m=remaining.mean(axis=0)
        qa=_quad(remaining,kernel)
        cosine_denom=np.linalg.norm(center)*np.linalg.norm(m)
        cosine=float(center@m/cosine_denom) if cosine_denom>1e-15 else None
        records.append({
            "scope":scope,"subset":subset,
            "dependency_grouping":dependency_grouping,
            "endpoint":endpoint,
            "n_pairs_original":len(pair_delta),
            "n_components_original":len(groups),
            "n_pairs_remaining":int(mask.sum()),
            "n_components_remaining":len(groups)-1,
            "component_size_removed":len(group),
            "remaining_statistic":qa,
            "ratio_norm_remaining_to_original":float(
                np.sqrt(qa/q_all)) if q_all>1e-15 else None,
            "cosine_with_original_mean_direction":cosine,
        })
    return records


def run(input_root:Path,g0_root:Path,crosswalk_path:Path,
        effort_path:Path,output:Path)->dict:
    input_root,g0_root,crosswalk_path,effort_path,output=map(
        Path,(input_root,g0_root,crosswalk_path,effort_path,output))
    output.mkdir(parents=True,exist_ok=True)
    g0={scope:list(csv_rows(g0_root/f"{scope}_islands.csv"))
        for scope in ("all","direct")}
    geo={r["island_id"]:r for r in csv_rows(crosswalk_path)}
    matched,audit=match_origin_islands(g0["all"],g0["direct"],geo)
    if len(matched)!=37 or sum(p["region"]=="tropical" for p in matched)!=14:
        raise ValueError("Unexpected original fd PR4 geology matched cohort")
    effort=read_effort(effort_path)
    space=ProfileSpace(load_profiles(input_root/REF_INPUT))
    results=[]
    dropout=[]
    component_report=[]
    for subset in ("all_37", "effort_balanced"):
        selected = (matched if subset == "all_37" else
                    eligible_pairs(matched, effort,
                                   "min30_balanced4x_uncertainty10km")[0])
        if len(selected) < 8:
            raise ValueError("Spatial pair support insufficient for specified screen")
        for dependency_grouping in ("spatial_block", "weigelt_archipelago"):
            components = connected_block_components(
                selected, geo, field=dependency_grouping)
            component_report.append({
                "subset": subset,
                "dependency_grouping": dependency_grouping,
                "n_pairs": len(selected),
                "n_connected_components": len(components),
                "component_sizes": "|".join(map(str, sorted(
                    (len(x) for x in components), reverse=True))),
                "n_northern_midlatitude": sum(
                    p["region"]=="northern_midlatitude" for p in selected),
                "n_tropical": sum(p["region"]=="tropical" for p in selected),
            })
            for scope in ("all", "direct"):
                carriers = load_carriers(
                    input_root/"floral_space_20261106"/scope/
                    "island_profile_carriers.csv.gz")
                delta = _delta(selected, carriers, space)
                for endpoint, kernel in (
                    ("exact_profile_L2", np.eye(len(space.ids))),
                    ("categorical_similarity_kernel", 1-space.matrix),
                ):
                    test = cluster_permutation(delta, components, kernel)
                    results.append({
                        "scope": scope, "subset": subset,
                        "dependency_grouping": dependency_grouping,
                        "endpoint": endpoint, "n_pairs": len(selected),
                        "n_components": len(components),
                        "largest_component_pairs": max(len(x) for x in components),
                        "component_size_distribution": "|".join(map(
                            str,sorted((len(x) for x in components),reverse=True))),
                        "status": "conservative_exact_geography_component_permutation",
                        **test,
                    })
                    dropout.extend(delete_component_sensitivity(
                        delta, components, kernel, scope=scope, subset=subset,
                        endpoint=endpoint, dependency_grouping=dependency_grouping))
    _write_csv(output/"spatial_component_null.csv",results,ROWS)
    _write_csv(output/"leave_one_spatial_component_out.csv",dropout,DROP_COLUMNS)
    _write_csv(output/"spatial_component_support.csv",component_report,
               tuple(component_report[0]))
    manifest={
        "schema":"fd_spatial_dependence_component_signflip_v1",
        "status":"conservative_geographic_dependency_stress_test_not_spatial_occupancy_model",
        "n_original_physical_pairs":len(matched),
        "n_all_original_components":component_report[0]["n_connected_components"],
        "n_dependency_groupings":2,
        "block_graph":"Pairs linked transitively when EITHER origin partner shares spatial_block OR independently within source-defined Weigelt archipelago; missing archipelago treated as individual endpoint; each grouping independently evaluated",
        "null":"all 2^K exact component-level pair-label sign configurations, squared norm invariant to global sign reversal",
        "caveat":"This extreme graph may over-group physically distant pairs into one component. Grid co-membership alone does not prove perfect dependence. But original independent pair sign swaps should not be read as certified spatially robust.",
        "n_tests":len(results),
        "source_hashes":{
            "crosswalk":file_sha256(crosswalk_path),
            "source_effort":file_sha256(effort_path),
            "profile_reference":file_sha256(input_root/REF_INPUT),
        },
        "not_estimated":["phylogenetically independent colonization",
                         "true absence or detection","geological selection",
                         "effective number of independent island lineages"],
        "privacy":"Output has only aggregate component sizes, null statistics and source hashes. No block labels, source IDs or island-pair IDs.",
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"support":component_report,
                      "tests":[{"scope":r["scope"],"subset":r["subset"],
                                "endpoint":r["endpoint"],"n_components":r["n_components"],
                                "exact_p":r["exact_component_upper_p"]}
                               for r in results]},ensure_ascii=False))
    return manifest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("input","g0","crosswalk","effort","output"):
        ap.add_argument("--"+name,type=Path,required=True)
    args=ap.parse_args()
    run(args.input,args.g0,args.crosswalk,args.effort,args.output)


if __name__=="__main__":
    main()
