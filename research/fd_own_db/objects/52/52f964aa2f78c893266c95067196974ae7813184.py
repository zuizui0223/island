"""Genus- and regional-availability-conditioned floral repertoire null.

The supplied profile maps and island flora are PRIVATE inputs reconstructed by
the upstream island PR321 runner. This module never writes taxon names, island
identifiers, pair identities, provenance tokens, or source rows to its outputs.
Only aggregate null summaries and SHA256 evidence can be published.

Genus + geographic-availability signature is a proxy for source opportunities,
NOT an ancestral mainland source pool or fully resolved phylogeny.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_profiles
from .g1_fixed_margins import _bh
from .repertoire import REGIONS, match_origin_islands

SCOPES = ("all", "direct")
POLICIES = ("genus_region_signature", "genus_global")
SCORES = ("exact_profile_L2", "categorical_similarity_kernel")
EXPECTED = {
    "all": {"complete_names": 16739, "distinct_profiles": 496,
            "flora_complete_names": 16725},
    "direct": {"complete_names": 8455, "distinct_profiles": 467,
               "flora_complete_names": 8452},
}
GEO_REGIMES = (
    "northern_midlatitude", "northern_high_latitude",
    "tropical", "southern_extratropical",
)
OUTPUT_COLS = (
    "scope", "policy", "region", "endpoint", "n_pairs", "n_islands",
    "n_complete_names_pool", "n_groups", "n_permutable_groups",
    "n_permutable_taxa", "draws", "observed", "null_median",
    "excess", "null_q025", "null_q975", "p_upper_approx",
    "q_bh_primary4", "q_by_primary4", "q_bh_secondary8",
    "q_by_secondary8",
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_private_maps(path: Path, scope: str, ids: set[int],
                      *, strict: bool = True) -> pd.DataFrame:
    """Read only the three essential columns of an authorized PR321 checkpoint.

    Never persist this DataFrame to the public project. Refuse name and profile
    duplicates before any inference.
    """
    d = pd.read_csv(path, usecols=["accepted_species", "genus", "profile_id"],
                    dtype={"accepted_species": str, "genus": str, "profile_id": int})
    if d.empty or d.accepted_species.isna().any() or d.genus.isna().any():
        raise ValueError(f"{scope}: invalid private profile mapping")
    if d.accepted_species.duplicated().any() or not set(d.profile_id).issubset(ids):
        raise ValueError(f"{scope}: repeated names or unrecognized profile IDs")
    if (d.genus != d.accepted_species.str.split().str[0]).any():
        raise ValueError(f"{scope}: genus/name identity inconsistency")
    if strict and (len(d) != EXPECTED[scope]["complete_names"]
                   or d.profile_id.nunique() != EXPECTED[scope]["distinct_profiles"]):
        raise ValueError(f"{scope}: upstream PR321 profile population changed")
    return d


def _read_membership(flora: Path, geography: pd.DataFrame, profile: pd.DataFrame,
                     *, native_only: bool = False) -> pd.DataFrame:
    f = pd.read_csv(flora, usecols=["island_id", "accepted_species", "origin_status"],
                    dtype=str)
    if f.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Private flora duplicates island/taxon-name rows")
    if native_only:
        f = f.loc[f.origin_status == "native"].copy()
    if geography.island_id.duplicated().any():
        raise ValueError("Geological crosswalk duplicated island IDs")
    merged = f[["island_id", "accepted_species"]].merge(
        geography[["island_id", "analysis_regime"]],
        on="island_id", how="inner", validate="many_to_one",
    ).merge(
        profile, on="accepted_species", how="inner", validate="many_to_one",
    )
    if merged.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Duplicate complete-profile island/taxon join")
    return merged


def validate_island_counts(membership: pd.DataFrame, g0: pd.DataFrame,
                           scope: str) -> dict:
    """Require every originally eligible island to reproduce G0 exactly."""
    counts = membership.groupby("island_id").agg(
        n_names=("accepted_species", "nunique"),
        n_profiles=("profile_id", "nunique"),
    )
    expected = g0.set_index("island_id")[["n_names", "n_profiles"]].astype(int)
    observed = expected.join(counts, how="left", rsuffix="_private")
    n = len(observed)
    if observed[["n_names_private", "n_profiles_private"]].isna().any().any():
        raise ValueError(f"{scope}: incomplete population in private flora")
    mismatches = (
        (observed.n_names != observed.n_names_private) |
        (observed.n_profiles != observed.n_profiles_private)
    )
    if mismatches.any():
        raise ValueError(f"{scope}: {int(mismatches.sum())} published island count mismatches")
    return {"eligible_islands_checked": n, "mismatches": 0}


def selected_native_pairs(matched: list[dict], all_native: pd.DataFrame,
                          direct_native: pd.DataFrame, min_names: int) -> list[dict]:
    """Taxon-status-based support gate, applied before using any trait outcomes."""
    if min_names < 1:
        raise ValueError("min_names must be positive")
    counts = (all_native.groupby("island_id").size().to_dict(),
              direct_native.groupby("island_id").size().to_dict())
    return [
        p for p in matched
        if all(c.get(p["continental_id"], 0) >= min_names and
               c.get(p["oceanic_id"], 0) >= min_names for c in counts)
    ]


def _source_input(membership: pd.DataFrame, matched: list[dict],
                  n_profiles: int) -> dict:
    # Construct four-regime incidence signature from the complete flora, not
    # from matched islands only. A taxon present in more than one regime keeps
    # one trait assignment in every randomization.
    names = membership[["accepted_species", "genus", "profile_id"]].drop_duplicates()
    if names.accepted_species.duplicated().any():
        raise ValueError("One taxon-name maps to multiple profiles")
    names = names.sort_values("accepted_species").reset_index(drop=True)
    signatures = membership.groupby("accepted_species")["analysis_regime"].agg(
        lambda x: ",".join(sorted(set(x)))
    ).to_dict()
    names["source_signature"] = names.accepted_species.map(signatures)
    row_index = {name: i for i, name in enumerate(names.accepted_species)}
    island_ids = sorted(set(id_ for pair in matched for id_ in (
        pair["continental_id"], pair["oceanic_id"]
    )))
    by_island = membership.loc[membership.island_id.isin(island_ids)].groupby(
        "island_id")["accepted_species"].agg(list).to_dict()
    if len(by_island) != len(island_ids):
        raise ValueError("Selected paired island is absent from private membership")
    incidence = [
        np.asarray([row_index[name] for name in by_island[island]], dtype=np.int32)
        for island in island_ids
    ]
    island_lookup = {key: i for i, key in enumerate(island_ids)}
    pair_indices = np.asarray([
        (island_lookup[p["oceanic_id"]], island_lookup[p["continental_id"]])
        for p in matched
    ], dtype=int)
    regional = {
        region: np.asarray([i for i,p in enumerate(matched)
                             if p["region"] == region], dtype=int)
        for region in REGIONS
    }
    if not len(names) or not len(incidence):
        raise ValueError("Insufficient private genotype/profile inputs")
    if not names.profile_id.between(0, n_profiles - 1).all():
        raise ValueError("Profile IDs outside reference")
    return {"names": names, "islands": incidence, "pairs": pair_indices,
            "regional": regional, "original": names.profile_id.to_numpy(dtype=np.int32)}


def _permutable_groups(names: pd.DataFrame, policy: str) -> tuple[list[np.ndarray],int,int]:
    if policy not in POLICIES:
        raise ValueError(f"Unrecognized null policy: {policy}")
    groups = defaultdict(list)
    for i, name in enumerate(names.itertuples(index=False)):
        key = name.genus if policy == "genus_global" else (name.genus, name.source_signature)
        groups[key].append(i)
    original = names.profile_id.to_numpy(dtype=np.int32)
    active = [np.asarray(ix,dtype=np.int32) for ix in groups.values()
              if len(ix) >= 2 and len(set(original[ix])) >= 2]
    return active, len(groups), sum(len(ix) for ix in active)


def _distance_statistics(assigned: np.ndarray, table: dict,
                         kernel: np.ndarray, n_profiles: int) -> dict[str,np.ndarray]:
    islands = np.zeros((len(table["islands"]), n_profiles), dtype=np.float64)
    for j, taxon_indices in enumerate(table["islands"]):
        islands[j, np.unique(assigned[taxon_indices])] = 1.0
    islands /= islands.sum(axis=1,keepdims=True)
    delta = islands[table["pairs"][:,0]] - islands[table["pairs"][:,1]]
    result = {}
    for name, positions in (
        ("global", np.arange(len(delta))), *list(table["regional"].items()),
    ):
        if not len(positions):
            continue
        mu = delta[positions].mean(axis=0)
        result[name] = np.asarray([float(mu @ mu), float(mu @ kernel @ mu)])
    return result


def _adjust(rows: list[dict], selector: str, bh_name: str, by_name: str) -> None:
    selected = [r for r in rows if (r["region"] == "global") == (selector == "global")]
    if not selected:
        return
    p = [r["p_upper_approx"] for r in selected]
    q = _bh(p)
    harmonic = sum(1 / i for i in range(1,len(selected)+1))
    for r, qv in zip(selected, q):
        r[bh_name] = float(qv)
        r[by_name] = float(min(1, qv * harmonic))


def infer_for_scope(scope: str, private_membership: pd.DataFrame,
                    matched: list[dict], space: ProfileSpace,
                    *, policy: str, draws: int, seed: int) -> list[dict]:
    """One global source-regional permutation shared across all matched islands."""
    if draws < 99:
        raise ValueError("Require at least 99 independent conditional draws")
    table = _source_input(private_membership, matched, len(space.ids))
    groups, n_groups, n_permutable_taxa = _permutable_groups(table["names"], policy)
    original = table["original"]
    kernel = 1.0 - space.matrix
    observed = _distance_statistics(original, table, kernel, len(space.ids))
    rng = np.random.default_rng(seed)
    null = {key: np.empty((draws, 2), dtype=float) for key in observed}
    for draw in range(draws):
        shuffled = original.copy()
        for indices in groups:
            shuffled[indices] = original[indices][rng.permutation(len(indices))]
        values = _distance_statistics(shuffled, table, kernel, len(space.ids))
        for key in null:
            null[key][draw] = values[key]
    result = []
    for region in observed:
        n_pairs = len(table["pairs"]) if region == "global" else len(table["regional"][region])
        for j, endpoint in enumerate(SCORES):
            stat = observed[region][j]
            distribution = null[region][:,j]
            tol = 100 * np.finfo(float).eps * max(1.0, abs(stat))
            p = (1 + int((distribution >= stat - tol).sum())) / (draws + 1)
            result.append({
                "scope": scope, "policy": policy, "region": region,
                "endpoint": endpoint, "n_pairs": n_pairs,
                "n_islands": len(table["islands"]),
                "n_complete_names_pool": len(table["names"]),
                "n_groups": n_groups, "n_permutable_groups": len(groups),
                "n_permutable_taxa": n_permutable_taxa, "draws": draws,
                "observed": float(stat), "null_median": float(np.median(distribution)),
                "excess": float(stat - np.median(distribution)),
                "null_q025": float(np.quantile(distribution, 0.025)),
                "null_q975": float(np.quantile(distribution, 0.975)),
                "p_upper_approx": float(p),
            })
    return result


def run(private_all: Path, private_direct: Path, flora: Path,
        crosswalk: Path, g0_root: Path, reference: Path, output: Path,
        *, draws: int = 1999, seed: int = 20261108,
        native_only: bool = False, min_native_names: int = 20,
        strict: bool = True) -> dict:
    """Compute and persist aggregate-only results. Fail on an upstream mismatch."""
    paths = (private_all, private_direct, flora, crosswalk, reference)
    private_all,private_direct,flora,crosswalk,g0_root,reference,output = map(
        Path, (*paths[:4],g0_root,reference,output))
    if len({str(p) for p in paths}) != len(paths):
        raise ValueError("Duplicate input paths")
    output.mkdir(parents=True, exist_ok=True)
    space = ProfileSpace(load_profiles(reference))
    profiles = {
        "all": load_private_maps(private_all, "all", set(space.ids), strict=strict),
        "direct": load_private_maps(private_direct, "direct", set(space.ids), strict=strict),
    }
    geology = pd.read_csv(crosswalk, dtype={"island_id": str})
    if geology.island_id.duplicated().any():
        raise ValueError("Repeated geological island IDs")
    g0 = {
        scope: pd.read_csv(g0_root / f"{scope}_islands.csv", dtype={"island_id": str})
        for scope in SCOPES
    }
    matched, matching = match_origin_islands(
        g0["all"].to_dict("records"), g0["direct"].to_dict("records"),
        geology.set_index("island_id").to_dict("index"),
    )
    if strict and (len(matched) != 37 or sum(p["region"]=="tropical" for p in matched)!=14):
        raise ValueError("Failed to recreate frozen 37-pair origin matching cohort")
    full = {scope: _read_membership(flora,geology,profiles[scope])
            for scope in SCOPES}
    validations = {scope: validate_island_counts(full[scope],g0[scope],scope)
                   for scope in SCOPES}
    if strict:
        for scope in SCOPES:
            if len(full[scope].accepted_species.unique()) != EXPECTED[scope]["flora_complete_names"]:
                raise ValueError(f"{scope}: upstream complete-flora name population changed")
    if native_only:
        native = {scope: _read_membership(flora,geology,profiles[scope],native_only=True)
                  for scope in SCOPES}
        matched = selected_native_pairs(matched,native["all"],native["direct"],
                                        min_native_names)
        if len(matched) < 8:
            raise ValueError("Insufficient paired native-only support")
        membership = native
    else:
        membership = full
    results = []
    for scope in SCOPES:
        for policy in POLICIES:
            offset = (17 if scope == "all" else 37)
            offset += 107 if policy == "genus_region_signature" else 311
            results.extend(infer_for_scope(
                scope,membership[scope],matched,space,policy=policy,
                draws=draws,seed=seed+offset,
            ))
    for policy in POLICIES:
        selected = [r for r in results if r["policy"]==policy]
        _adjust(selected,"global","q_bh_primary4","q_by_primary4")
        _adjust(selected,"regional","q_bh_secondary8","q_by_secondary8")
    rows = pd.DataFrame(results).reindex(columns=OUTPUT_COLS)
    rows.to_csv(output / "genus_conditioned_summary.csv",index=False)
    report = {
        "schema": "fd_genus_regional_availability_v1",
        "status": "computed_private_input_aggregate_only_exploratory",
        "cohort":{"n_pairs":len(matched),
                  "n_northern_midlatitude":sum(p["region"]=="northern_midlatitude" for p in matched),
                  "n_tropical":sum(p["region"]=="tropical" for p in matched),
                  "native_only":native_only,
                  "native_min_complete_names_per_endpoint":min_native_names if native_only else None},
        "verification": validations,
        "nulls":{
            "genus_region_signature":"taxon profiles permuted among names of same genus and identical observed four-regime island-flora presence signature; one assignment reused at all islands",
            "genus_global":"taxon profiles permuted within genus across all observed islands; one globally coherent taxon assignment",
        },
        "statistic":"mean across paired oceanic−continental island unique-profile-normalized incidence vectors, quadratic L2 or categorical similarity kernel",
        "conditional_on":["fixed island × taxon incidence","fixed name-level complete-colour/form support",
                          "genus profile multisets","source-regional occurrence signature for primary null",
                          "fixed geological origin matching and origin labels"],
        "NOT_conditioned_on":["true mainland colonist sources","fully resolved phylogeny",
                              "native status when native_only=false",
                              "individual trait observation errors","pollinator assemblages"],
        "main_inference":"primary genus x geography-signature null BH/BY over 4 global comparisons; region comparisons are secondary 8 family; 2 scopes nested",
        "draws":draws,"seed":seed,
        "privacy":"no private species/island/profile mapping, name keys, row provenance, source citations or pair IDs emitted",
        "input_sha256":{"private_all_profiles":_sha(private_all),
                        "private_direct_profiles":_sha(private_direct),
                        "private_flora":_sha(flora),
                        "geology_crosswalk":_sha(crosswalk),
                        "profile_reference":_sha(reference)},
    }
    (output / "manifest.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"cohort":report["cohort"],
                      "primary":[{"scope":r["scope"],"metric":r["endpoint"],
                                  "p":r["p_upper_approx"],"q":r.get("q_bh_primary4")}
                                 for r in results if r["policy"]==POLICIES[0] and r["region"]=="global"]}))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("private-all","private-direct","flora","crosswalk","g0","reference","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--draws",type=int,default=1999)
    p.add_argument("--seed",type=int,default=20261108)
    p.add_argument("--native-only",action="store_true")
    p.add_argument("--min-native-names",type=int,default=20)
    p.add_argument("--allow-synthetic",action="store_true",
                   help="Disable only fixed reference-population checks for synthetic fixtures")
    a=p.parse_args()
    run(a.private_all,a.private_direct,a.flora,a.crosswalk,a.g0,a.reference,a.output,
        draws=a.draws,seed=a.seed,native_only=a.native_only,
        min_native_names=a.min_native_names,strict=not a.allow_synthetic)


if __name__=="__main__":
    main()
