"""Public-aggregate locality-matched floral distance support benchmark N1c.

Is original-distance flower category packing unusual compared with geographically
nearby OTHER island recorded profile mixtures, rather than a broad regional
pool? The ten closest eligible same-regime islands are used as a fixed-size
observed donor set; their shared species may be unknown. This is an explicitly
post-hoc, observational GEOGRAPHY reference, NOT a genus or mainland N2 null.

No individual island IDs, coordinates, donor lists or local pair identities
are released. No geological-origin significance tests are run here.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .distance_ses import (
    ENDPOINTS, MODES, SUMMARY_FIELDS, _score, _seed, _summary,
    distance_statistics, regional_leave_one_out_pool,
    regional_profile_counts, sample_null,
)
from .pipeline import REF_INPUT, _write_csv, file_sha256

SCOPES = ("all", "direct")
NULLS = ("N1b_regional_islandfreq_LOO", "N1c_nearest10_geo_LOO")
MATCHED_FIELDS = (
    "scope", "mode", "endpoint", "region", "geological_origin",
    "n_islands_both_nulls", "status",
    "median_N1b_SES", "median_N1c_SES", "median_N1c_minus_N1b_SES",
    "fraction_negative_N1b", "fraction_negative_N1c",
    "median_N1b_observed_minus_null", "median_N1c_observed_minus_null",
    "median_N1c_minus_N1b_null_mean",
)
SUPPORT_FIELDS = (
    "scope", "geological_origin", "n_geology_islands",
    "n_valid_geographic_neighbor_sets",
    "n_n1c_full_estimable", "n_n1c_equal30_estimable",
    "median_n_profile_types_in_nearest10",
    "median_10th_neighbor_distance_km", "status",
)
RADIUS_KM = 6371.0088


def great_circle_distance_km(lat: float, lon: float,
                             candidate_lat: np.ndarray,
                             candidate_lon: np.ndarray) -> np.ndarray:
    """Vectorized haversine; no invented distance or coordinate metadata."""
    lat0, lon0 = float(lat), float(lon)
    phi = np.asarray(candidate_lat, dtype=float)
    lam = np.asarray(candidate_lon, dtype=float)
    if not np.isfinite([lat0, lon0]).all() or not (
        -90 <= lat0 <= 90 and -180 <= lon0 <= 180):
        raise ValueError("Invalid geographic reference coordinate")
    if phi.shape != lam.shape or np.any(~np.isfinite(phi)) or np.any(~np.isfinite(lam)):
        raise ValueError("Invalid geographic donor coordinate vectors")
    if np.any(abs(phi)>90) or np.any(abs(lam)>180):
        raise ValueError("Geographic donor coordinate outside Earth bounds")
    x1 = np.radians(phi-lat0)/2
    x2 = np.radians(lam-lon0)/2
    a = np.sin(x1)**2 + np.cos(np.radians(lat0)) * np.cos(np.radians(phi)) * np.sin(x2)**2
    return 2*RADIUS_KM*np.arcsin(np.sqrt(np.clip(a,0,1)))


def nearest_regional_donors(focal: str, regime: str,
                            metadata: dict[str, dict],
                            geography: dict[str, dict], *, k: int = 10
                            ) -> tuple[list[str], float]:
    """Return k nearest OTHER eligible recorded islands in the same broad regime."""
    if k < 2 or focal not in metadata or focal not in geography:
        raise ValueError("Unsupported target geographic donor query")
    lat, lon = (float(geography[focal][field])
                for field in ("island_latitude", "island_longitude"))
    candidates = []
    coords = []
    for key, row in metadata.items():
        if key == focal or row["analysis_regime"] != regime or key not in geography:
            continue
        gr = geography[key]
        try:
            y, x = float(gr["island_latitude"]), float(gr["island_longitude"])
        except (ValueError, KeyError, TypeError):
            continue
        if not np.isfinite([y,x]).all() or not (-90 <= y <= 90 and -180 <= x <= 180):
            continue
        candidates.append(key)
        coords.append((y, x))
    if len(candidates)<k:
        return [], float("nan")
    matrix=np.asarray(coords,dtype=float)
    dist=great_circle_distance_km(lat,lon,matrix[:,0],matrix[:,1])
    selected=sorted(zip(dist,candidates),key=lambda x:(x[0],x[1]))[:k]
    return [key for _,key in selected],float(selected[-1][0])


def local_pool(donors: list[str], carriers: dict[str,dict[int,int]],
               lookup: dict[int,int]) -> tuple[np.ndarray,np.ndarray]:
    """Use donor-ISLAND incidence weights, never source name abundances."""
    count=Counter()
    for island in donors:
        if island not in carriers:
            raise ValueError("Nearest neighbor island absent from public source")
        count.update(carriers[island].keys())
    ids=sorted(count)
    return (np.asarray([lookup[k] for k in ids],dtype=int),
            np.asarray([count[k] for k in ids],dtype=float))


def _support_summary(records: list[dict]) -> list[dict]:
    group=defaultdict(list)
    for r in records:
        group[(r["scope"],r["geological_origin"])].append(r)
    output=[]
    for (scope,origin),data in sorted(group.items()):
        row={
            "scope":scope,"geological_origin":origin,"n_geology_islands":len(data),
            "n_valid_geographic_neighbor_sets":sum(x["has_neighbors"] for x in data),
            "n_n1c_full_estimable":sum(x["has_full"] for x in data),
            "n_n1c_equal30_estimable":sum(x["has_equal30"] for x in data),
            "status":"descriptive_nearest10_geographic_reference",
        }
        valid=[x for x in data if x["has_neighbors"]]
        if len(valid)>=5:
            row["median_n_profile_types_in_nearest10"]=float(np.median(
                [x["n_donor_profiles"] for x in valid]))
            row["median_10th_neighbor_distance_km"]=float(np.median(
                [x["radius_km"] for x in valid]))
        else:
            row["status"]="small_cell_suppressed"
        output.append(row)
    return output


def paired_reference_summary(scores: list[dict], *, min_release_n:int=5) -> list[dict]:
    """Compare N1b and N1c ONLY on the same physical islands and same metric.

    This descriptive within-island paired change avoids interpreting unmatched
    group medians as the effect of moving from regional to local reference.
    """
    grouped=defaultdict(list)
    by_island=defaultdict(dict)
    for r in scores:
        key=(r["scope"],r["mode"],r["endpoint"],r["region"],
             r["geological_origin"],r["island_id"])
        if r["null"] in NULLS:
            if r["null"] in by_island[key]:
                raise ValueError("Duplicate island×null×metric in local reference")
            by_island[key][r["null"]]=r
    for key,variants in by_island.items():
        if len(variants)!=len(NULLS):
            continue
        n1b,n1c=[variants[null] for null in NULLS]
        if not np.isclose(n1b["observed"],n1c["observed"],rtol=1e-8,atol=1e-8):
            raise AssertionError("Same-island baseline test must hold observations fixed")
        common=(n1b["scope"],n1b["mode"],n1b["endpoint"],n1b["region"],
                n1b["geological_origin"])
        grouped[common].append((n1b,n1c))
    report=[]
    # Fold both broad regime-specific and pooled "all" summaries only after
    # evaluating each physical island on both nulls.
    folded=defaultdict(list)
    for key,records in grouped.items():
        folded[key].extend(records)
        if key[3]!="all":
            folded[(key[0],key[1],key[2],"all",key[4])].extend(records)
    for key,records in sorted(folded.items()):
        if not records:
            continue
        row=dict(zip(("scope","mode","endpoint","region","geological_origin"),key))
        row["n_islands_both_nulls"]=len(records)
        row["status"]=("descriptive_paired_reference_comparison"
                       if len(records)>=min_release_n else "small_cell_suppressed")
        if len(records)>=min_release_n:
            b=np.asarray([pair[0]["ses"] for pair in records])
            n=np.asarray([pair[1]["ses"] for pair in records])
            db=np.asarray([pair[0]["delta"] for pair in records])
            dn=np.asarray([pair[1]["delta"] for pair in records])
            row.update({
                "median_N1b_SES":float(np.median(b)),
                "median_N1c_SES":float(np.median(n)),
                "median_N1c_minus_N1b_SES":float(np.median(n-b)),
                "fraction_negative_N1b":float(np.mean(b<0)),
                "fraction_negative_N1c":float(np.mean(n<0)),
                "median_N1b_observed_minus_null":float(np.median(db)),
                "median_N1c_observed_minus_null":float(np.median(dn)),
                "median_N1c_minus_N1b_null_mean":float(np.median([
                    pair[1]["null_mean"]-pair[0]["null_mean"] for pair in records])),
            })
        report.append(row)
    return report


def run(input_root:Path,g0_root:Path,crosswalk:Path,output:Path,
        *,neighbor_k:int=10,draws:int=999,thinnings:int=99,
        seed:int=20261109)->dict:
    if neighbor_k <2 or draws<20 or thinnings<20:
        raise ValueError("Insufficient fixed neighbor/draw/sampling support")
    input_root,g0_root,crosswalk,output=map(
        Path,(input_root,g0_root,crosswalk,output))
    space=ProfileSpace(load_profiles(input_root/REF_INPUT))
    if len(space.ids)!=496:
        raise ValueError("Frozen floral-profile reference changed")
    geo={r["island_id"]:r for r in csv_rows(crosswalk)}
    score_rows=[]
    support_rows=[]
    carrier_files={}
    for scope in SCOPES:
        g0=list(csv_rows(g0_root/f"{scope}_islands.csv"))
        meta={r["island_id"]:r for r in g0}
        carrier_path=input_root/"floral_space_20261106"/scope/"island_profile_carriers.csv.gz"
        carrier_files[scope]=carrier_path
        carriers=load_carriers(carrier_path)
        if set(meta)-set(carriers):
            raise ValueError("G0 eligible islands missing public floral carriers")
        regional=regional_profile_counts(g0,carriers)
        for r in g0:
            island=r["island_id"]
            origin=geo.get(island,{}).get("geological_origin")
            if origin not in ("continental","oceanic"):
                continue
            profiles=set(carriers[island])
            selected=np.asarray([space.lookup[p] for p in sorted(profiles)],dtype=int)
            n=len(selected)
            if n<2:
                continue
            donors,radius=nearest_regional_donors(
                island,r["analysis_regime"],meta,geo,k=neighbor_k)
            near_pool,near_weights=(local_pool(donors,carriers,space.lookup)
                                    if donors else (np.asarray([],dtype=int),np.asarray([],dtype=float)))
            regional_pool=regional_leave_one_out_pool(
                regional,r["analysis_regime"],profiles,space.lookup)
            regional_weights=np.asarray([
                regional[r["analysis_regime"]][space.ids[int(idx)]]
                - int(space.ids[int(idx)] in profiles)
                for idx in regional_pool
            ],dtype=float)
            if len(regional_weights) and np.any(regional_weights<=0):
                raise AssertionError("Regional LOO incidence nonpositive")
            support_rows.append({
                "scope":scope,"geological_origin":origin,
                "has_neighbors":bool(donors),"has_full":len(near_pool)>=n,
                "has_equal30":n>=30 and len(near_pool)>=30,
                "n_donor_profiles":len(near_pool),"radius_km":radius,
            })
            thin=None
            if n>=30:
                rng=np.random.default_rng(_seed(seed,scope,island,"N1c_30_thin"))
                thin=np.asarray([
                    distance_statistics(space.matrix,rng.choice(selected,30,replace=False))
                    for _ in range(thinnings)
                ])
            for mode in MODES:
                if mode=="equal30" and thin is None:
                    continue
                n_sample=n if mode=="full" else 30
                observed=(distance_statistics(space.matrix,selected)
                          if mode=="full" else thin.mean(axis=0))
                for null in NULLS:
                    if null=="N1b_regional_islandfreq_LOO":
                        pool,weights=regional_pool,regional_weights
                    else:
                        pool,weights=near_pool,near_weights
                    if len(pool)<n_sample:
                        continue
                    sim=sample_null(
                        space.matrix,pool,n_sample,draws=draws,
                        seed=_seed(seed,null,scope,island,mode),
                        weights=weights)
                    for value in _score(observed,sim,
                                        thinning=thin if mode=="equal30" else None):
                        if value is not None:
                            score_rows.append({
                                "scope":scope,"null":null,"mode":mode,
                                "region":r["analysis_regime"],
                                "geological_origin":origin,
                                "island_id":island,"n_profiles":n,
                                **value
                            })
    summary=_summary(score_rows)
    support=_support_summary(support_rows)
    paired=paired_reference_summary(score_rows)
    output.mkdir(parents=True,exist_ok=True)
    _write_csv(output/"local_geography_ses_summary.csv",summary,SUMMARY_FIELDS)
    _write_csv(output/"local_geography_support.csv",support,SUPPORT_FIELDS)
    _write_csv(output/"local_geography_matched_reference_comparison.csv",paired,MATCHED_FIELDS)
    manifest={
        "schema":"fd_nearestK_geography_flora_fixed_support_SES_v1",
        "status":"public_aggregate_posthoc_geographic_falsification",
        "nearest_neighbors":neighbor_k,
        "nulls":{
            "N1b":"focal-excluded same-regime profile incidence in ALL other eligible recorded islands",
            "N1c":"profile incidence among k nearest OTHER same-regime eligible islands only, distance from audited island centroids",
        },
        "sampling":"fixed distinct recorded n or equal30, donor island incidence-weighted draws without replacement",
        "paired_reference_validation":"N1b vs N1c evaluated on IDENTICAL target islands; no geographic-origin hypothesis test, compare SES and original D differences jointly",
        "draws":draws,"thinnings":thinnings,"seed":seed,
        "scope":"public All/Direct overlapping evidence classes",
        "biological_limit":"geography-proxy reference only; shared taxa/genera not removed; observation processes and spatial dependence persist; no origin hypothesis test here",
        "privacy":"only n>=5 grouped scores and aggregate neighbor radius/support, no exact coordinates/island IDs/individual neighbor lists",
        "source_sha256":{
            "reference":file_sha256(input_root/REF_INPUT),
            "crosswalk":file_sha256(crosswalk),
            **{scope+"_carrier":file_sha256(path) for scope,path in carrier_files.items()},
        }
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({
        "schema":manifest["schema"],"n_group_summary":len(summary),
        "n_support_groups":len(support),
        "n_paired_reference_groups":len(paired),
        "n_internal_scores":len(score_rows),
        "n_geology_islands":len(support_rows),
        "interpretation":"only observational neighborhood reference; not genus/source conditional N2",
    }))
    return manifest


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for field in ("input","g0","crosswalk","output"):
        p.add_argument("--"+field,type=Path,required=True)
    p.add_argument("--neighbors",type=int,default=10)
    p.add_argument("--draws",type=int,default=999)
    p.add_argument("--thinnings",type=int,default=99)
    args=p.parse_args()
    run(args.input,args.g0,args.crosswalk,args.output,
        neighbor_k=args.neighbors,draws=args.draws,thinnings=args.thinnings)


if __name__=="__main__":
    main()
