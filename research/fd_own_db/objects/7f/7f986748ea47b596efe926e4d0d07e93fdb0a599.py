"""Post-hoc equal-support rarefaction of geology-matched floral repertoires.

Scientific question: does the oceanic/continental difference between EXACT
categorical flower-colour × form profiles depend on different observed
numbers of complete taxa or category profiles per island?

No missing species are imputed here. These are sensitivity tests against
unequal observation support, NOT an occupancy/detection model, native-only
phylogenetic test or ecological colonization inference.

Original 37 pairs are selected without using floral outcomes. Every test
reads only rights-cleared island×profile carrier aggregates; it emits
ONLY aggregated statistics, no private island matches or trait provenance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .pipeline import REF_INPUT, file_sha256, _write_csv
from .repertoire import match_origin_islands, _paired_swap_statistic, _upper_tail

SCOPES = ("all", "direct")
SCREENS = (
    "original_full",
    "rarefy_30_profiles",
    "rarefy_40_name_units",
    "rarefy_to_pair_min_profiles",
    "single_state_20_profiles",
)
ENDPOINTS = ("exact_profile_L2", "categorical_similarity_kernel")
FIELDS = (
    "scope", "screen", "endpoint", "status", "n_pairs",
    "n_northern_midlatitude", "n_tropical", "n_resamples",
    "n_permutations_each", "observed_statistic_median",
    "observed_statistic_q025", "observed_statistic_q975",
    "paired_null_median", "p_median", "p_q025", "p_q975",
    "fraction_resamples_p_below_005",
)


def _seed(seed: int, *labels: str) -> int:
    payload = ":".join([str(seed), *labels]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _draw_one_carrier_support(
    carriers: dict[int,int], *, screen: str, rng: np.random.Generator,
    target_profiles: int | None = None,
) -> list[int]:
    if not carriers:
        raise ValueError("An island with zero complete traits cannot be rarefied")
    types = np.array(sorted(carriers), dtype=int)
    counts = np.asarray([carriers[int(pid)] for pid in types], dtype=int)
    if np.any(counts <= 0):
        raise ValueError("All carrier counts must be positive integers")
    if screen in ("original_full", "single_state_20_profiles"):
        return types.tolist()
    if screen == "rarefy_30_profiles":
        if len(types) < 30:
            raise ValueError("Cannot rarefy an island with <30 profiles")
        return sorted(rng.choice(types, size=30, replace=False).tolist())
    if screen == "rarefy_to_pair_min_profiles":
        if target_profiles is None or not 1 <= target_profiles <= len(types):
            raise ValueError("Invalid pair-common profile target")
        return sorted(rng.choice(types, size=target_profiles, replace=False).tolist())
    if screen == "rarefy_40_name_units":
        if counts.sum() < 40:
            raise ValueError("Cannot subsample 40 names from the island")
        sampled = rng.multivariate_hypergeometric(counts, 40)
        return types[sampled > 0].tolist()
    raise ValueError(f"Unrecognized rarefaction rule: {screen}")


def _sign_swap_scores(
    delta: np.ndarray, kernels: dict, signs: np.ndarray,
) -> dict[str, dict]:
    if delta.ndim != 2 or signs.shape[1] != len(delta):
        raise ValueError("Pair contrast dimensions do not match")
    if len(delta) < 8:
        raise ValueError("At least eight independent matched pairs are required")
    output = {}
    for endpoint, kernel in kernels.items():
        observed, reference = _paired_swap_statistic(delta, kernel, signs)
        output[endpoint] = {
            "observed": float(observed),
            "null_median": float(np.median(reference)),
            "p": float(_upper_tail(observed, reference)),
        }
    return output


def _model_one_scope(
    screen: str, scope: str, matches: list[dict], carriers: dict,
    space: ProfileSpace, kernels: dict, *, replications: int,
    permutations: int, seed: int,
) -> list[dict]:
    if screen not in SCREENS:
        raise ValueError("Invalid sensitivity screen")
    if screen == "single_state_20_profiles":
        available = [p for p in matches
                     if all(i in carriers and len(carriers[i]) >= 20
                            for i in (p["continental_id"], p["oceanic_id"]))]
    else:
        available = matches
    summary = {
        "scope": scope, "screen": screen, "n_pairs": len(available),
        "n_northern_midlatitude": sum(p["region"]=="northern_midlatitude" for p in available),
        "n_tropical": sum(p["region"]=="tropical" for p in available),
    }
    if len(available) < 8:
        return [{**summary, "endpoint": e, "status": "insufficient_pairs"}
                for e in ENDPOINTS]
    n_reps = 1 if screen in ("original_full", "single_state_20_profiles") else replications
    n_signs = 9999 if n_reps == 1 else permutations
    rng = np.random.default_rng(_seed(seed, scope, screen, "subsample"))
    swap_rng = np.random.default_rng(_seed(seed, scope, screen, "signs"))
    signs = swap_rng.choice((-1., 1.), size=(n_signs,len(available)))
    results = {e: [] for e in ENDPOINTS}
    for replicate in range(n_reps):
        pair_delta = []
        for p in available:
            continental = carriers[p["continental_id"]]
            oceanic = carriers[p["oceanic_id"]]
            target = min(len(continental),len(oceanic)) if screen == "rarefy_to_pair_min_profiles" else None
            c = _draw_one_carrier_support(
                continental, screen=screen, rng=rng, target_profiles=target)
            o = _draw_one_carrier_support(
                oceanic, screen=screen, rng=rng, target_profiles=target)
            diff = np.zeros(len(space.ids),dtype=float)
            diff[[space.lookup[k] for k in o]] += 1/len(o)
            diff[[space.lookup[k] for k in c]] -= 1/len(c)
            pair_delta.append(diff)
        stats = _sign_swap_scores(np.vstack(pair_delta), kernels, signs)
        for endpoint in ENDPOINTS:
            results[endpoint].append(stats[endpoint])
    output=[]
    for endpoint, draws in results.items():
        observed=np.array([d["observed"] for d in draws])
        pval=np.array([d["p"] for d in draws])
        refs=np.array([d["null_median"] for d in draws])
        output.append({
            **summary, "endpoint":endpoint,
            "status": "exploratory_posthoc_support_sensitivity",
            "n_resamples":n_reps, "n_permutations_each":n_signs,
            "observed_statistic_median":float(np.median(observed)),
            "observed_statistic_q025":float(np.quantile(observed,.025)),
            "observed_statistic_q975":float(np.quantile(observed,.975)),
            "paired_null_median":float(np.median(refs)),
            "p_median":float(np.median(pval)),
            "p_q025":float(np.quantile(pval,.025)),
            "p_q975":float(np.quantile(pval,.975)),
            "fraction_resamples_p_below_005":float(np.mean(pval<.05)),
        })
    return output


def run(input_root: Path,g0_root: Path,crosswalk_path: Path,output: Path,
        *,replications:int=99,permutations:int=1999,seed:int=20261109)->dict:
    if replications<20 or permutations<999:
        raise ValueError("Insufficient rarefaction/permutation support")
    input_root,g0_root,crosswalk_path,output=map(Path,(input_root,g0_root,crosswalk_path,output))
    output.mkdir(parents=True,exist_ok=True)
    g0={scope:list(csv_rows(g0_root/f"{scope}_islands.csv"))
        for scope in ("all","direct")}
    geology={r["island_id"]:r for r in csv_rows(crosswalk_path)}
    matched,matching=match_origin_islands(g0["all"],g0["direct"],geology)
    if len(matched)!=37 or sum(x["region"]=="tropical" for x in matched)!=14:
        raise ValueError("Fossil geology origin pair cohort changed relative to frozen fd PR4")
    space=ProfileSpace(load_profiles(input_root/REF_INPUT))
    kernels={
        "exact_profile_L2":np.eye(len(space.ids)),
        "categorical_similarity_kernel":1-space.matrix,
    }
    rows=[]
    validations=[]
    for scope in SCOPES:
        carriers=load_carriers(input_root/"floral_space_20261106"/scope/"island_profile_carriers.csv.gz")
        # Validate full name/profile support on every selected island against the
        # independently pinned published G0; fail closed on any mismatch.
        g0_map={x["island_id"]:x for x in g0[scope]}
        for p in matched:
            for identifier in (p["continental_id"],p["oceanic_id"]):
                if identifier not in g0_map or identifier not in carriers:
                    raise ValueError("Matched island missing from G0/support")
                c=carriers[identifier]
                r=g0_map[identifier]
                if sum(c.values())!=int(r["n_names"]) or len(c)!=int(r["n_profiles"]):
                    raise ValueError("Matched island carrier support differs from source manifest")
        counts=[len(carriers[i]) for p in matched for i in
                (p["continental_id"],p["oceanic_id"])]
        namecounts=[sum(carriers[i].values()) for p in matched for i in
                    (p["continental_id"],p["oceanic_id"])]
        validations.append({
            "scope":scope,"n_matched_islands":len(counts),
            "min_n_names":min(namecounts),"min_n_profiles":min(counts),
            "median_n_profiles":float(np.median(counts)),
            "n_islands_supporting_30_profile_rarefaction":sum(k>=30 for k in counts),
            "n_islands_supporting_40_name_rarefaction":sum(k>=40 for k in namecounts),
        })
        single_scope=f"{scope}_single_state"
        single=load_carriers(input_root/"floral_space_20261106"/single_scope/"island_profile_carriers.csv.gz")
        single_g0={r["island_id"] for r in csv_rows(g0_root/f"{single_scope}_islands.csv")}
        # Use the single-state **G0-eligible** population (>=50 complete name units)
        # rather than sneaking in islands with scant single-state support.
        single={i:v for i,v in single.items() if i in single_g0}
        for screen in SCREENS:
            reference=single if screen=="single_state_20_profiles" else carriers
            rows.extend(_model_one_scope(
                screen,scope,matched,reference,space,kernels,
                replications=replications,permutations=permutations,seed=seed))
    _write_csv(output/"support_rarefaction_results.csv",rows,FIELDS)
    manifest={
        "schema":"fd_profile_and_name_support_rarefaction_v1",
        "status":"exploratory_recorded_support_rarefaction_not_missing_species_imputation",
        "frozen_original_pairs":len(matched),
        "pair_geographic_regimes":{"northern_midlatitude":23,"tropical":14},
        "public_reference_profile_count":len(space.ids),
        "screens":list(SCREENS),"endpoints":list(ENDPOINTS),
        "sampled_profile_count_per_island":30,
        "sampled_carrier_name_count_per_island":40,
        "minimum_single_state_profiles_per_endpoint":20,
        "single_state_eligibility":"original independent G0 >=50 complete name units and >=20 distinct observed profiles",
        "stochastic_resamples":replications,"permutations_per_resample":permutations,
        "single_state_permutations":9999,
        "support_validation":validations,
        "not_claimed":[
            "true species absence","actual missing trait values","nonrandom trait recording correction",
            "colonization or establishment rate","phylogenetic source-pool balance",
            "evolutionary adaptation","independent external replication"],
        "inferential_note":"rarefaction p distributions are correlated, post-hoc sensitivity, not new independent confirmations; All/Direct share taxa and island pairs",
        "source_sha256":{
            "profile_reference":file_sha256(input_root/REF_INPUT),
            "geology_crosswalk":file_sha256(crosswalk_path),
            **{f"{scope}_carriers":file_sha256(
                input_root/"floral_space_20261106"/scope/"island_profile_carriers.csv.gz")
                for scope in (*SCOPES,"all_single_state","direct_single_state")},
        },
        "privacy":"aggregate statistics only; row-level island/geology matching never emitted",
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
        "scope_validation":validations,
        "results":[{"scope":r["scope"],"screen":r["screen"],"endpoint":r["endpoint"],
                    "n_pairs":r["n_pairs"],"p_median":r.get("p_median"),
                    "fraction_below_005":r.get("fraction_resamples_p_below_005")}
                   for r in rows]},ensure_ascii=False))
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for field in ("input","g0","crosswalk","output"):
        p.add_argument("--"+field,type=Path,required=True)
    p.add_argument("--replications",type=int,default=99)
    p.add_argument("--permutations",type=int,default=1999)
    a=p.parse_args()
    run(a.input,a.g0,a.crosswalk,a.output,
        replications=a.replications,permutations=a.permutations)


if __name__=="__main__":
    main()
