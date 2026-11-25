"""Rights-protected sufficient-statistic check for functional trait genus effects.

Input: private pair x trait positive/observed carrier counts for intact, species-
exclusive, and genus-exclusive floras. Output: coarse aggregated results only.
The main score compares observed genus-exclusive prevalence gaps with random
equal-carrier thinning of ALREADY species-exclusive records. Neither source
identities nor taxon/island ids are ever exported.

These are observational, non-independent descriptive diagnostics, NOT p-values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .functional_turnover import TRAIT_GROUPS, TRAITS

PREFIXES = ("all", "species", "genus")
KEYS = ("scope", "region", "island_a", "island_b", "trait")
COUNTS = tuple(f"{p}_{v}_{s}" for p in PREFIXES for v in ("n", "k") for s in ("a","b"))
REQUIRED = (*KEYS, *COUNTS)
DOMAINS = ("reproductive_assurance", "floral_accessibility")


def validate_counts(data: pd.DataFrame) -> pd.DataFrame:
    missing = set(REQUIRED) - set(data.columns)
    if missing:
        raise ValueError(f"Missing count fields: {sorted(missing)}")
    x = data[list(REQUIRED)].copy()
    for key in KEYS:
        x[key] = x[key].astype("string").str.strip()
        if (x[key].isna().any() or x[key].eq("").any()
                or x[key].str.lower().isin(["nan","na","unknown","none"]).any()):
            raise ValueError(f"Missing {key}")
    if not x.scope.isin(("All","Direct")).all() or not x.trait.isin(TRAITS).all():
        raise ValueError("Unknown evidence scope or frozen trait")
    if (x.island_a == x.island_b).any() or x.duplicated(list(KEYS)).any():
        raise ValueError("Self pair or duplicate pair/trait")
    for col in COUNTS:
        y = pd.to_numeric(x[col], errors="coerce")
        if y.isna().any() or not np.isfinite(y).all() or (y < 0).any() or (y != np.floor(y)).any():
            raise ValueError(f"Invalid nonnegative integer carrier count: {col}")
        x[col] = y.astype(np.int64)
    for side in ("a","b"):
        for p in PREFIXES:
            if (x[f"{p}_k_{side}"] > x[f"{p}_n_{side}"]).any():
                raise ValueError(f"Positive count exceeds observed carriers: {p}/{side}")
        if ((x[f"genus_n_{side}"] > x[f"species_n_{side}"]) |
                (x[f"species_n_{side}"] > x[f"all_n_{side}"])).any():
            raise ValueError(f"Non-nested support: {side}")
        for parent, sub in (("all","species"),("species","genus")):
            if ((x[f"{sub}_k_{side}"] > x[f"{parent}_k_{side}"]) |
                    ((x[f"{sub}_n_{side}"]-x[f"{sub}_k_{side}"]) >
                     (x[f"{parent}_n_{side}"]-x[f"{parent}_k_{side}"]))).any():
                raise ValueError(f"Non-nested positive/negative carrier counts: {side}")
    return x


def _expected_gap(n_a: np.ndarray, k_a: np.ndarray, m_a: np.ndarray,
                  n_b: np.ndarray, k_b: np.ndarray, m_b: np.ndarray,
                  *, rng: np.random.Generator, draws: int) -> np.ndarray:
    if draws < 100:
        raise ValueError("At least 100 null draws required")
    out = np.empty(len(n_a), dtype=float)
    for start in range(0, len(out), 1024):
        sl = slice(start,min(start+1024,len(out)))
        ya = rng.hypergeometric(k_a[sl],n_a[sl]-k_a[sl],m_a[sl],
                                size=(draws, len(out[sl])))
        yb = rng.hypergeometric(k_b[sl],n_b[sl]-k_b[sl],m_b[sl],
                                size=(draws, len(out[sl])))
        out[sl] = np.abs(ya/m_a[sl]-yb/m_b[sl]).mean(axis=0)
    return out


def score_pairs(data: pd.DataFrame, *, min_carriers: int = 20,
                draws: int = 399, seed: int = 20261110) -> pd.DataFrame:
    """Private per-pair scores, never directly export these data."""
    if min_carriers < 1: raise ValueError("min_carriers must be positive")
    x = validate_counts(data)
    x = x[(x.genus_n_a>=min_carriers)&(x.genus_n_b>=min_carriers)].copy()
    outcols = [*KEYS,"species_support_residual","extra_genus_gap_beyond_support"]
    if x.empty: return pd.DataFrame(columns=outcols)
    x = x.sort_values(list(KEYS),kind="mergesort").reset_index(drop=True)
    rng = np.random.default_rng(seed)
    for stage,parent in (("species","all"),("genus","species")):
        observed = np.abs(
            x[f"{stage}_k_a"].to_numpy()/x[f"{stage}_n_a"].to_numpy() -
            x[f"{stage}_k_b"].to_numpy()/x[f"{stage}_n_b"].to_numpy())
        reference = _expected_gap(
            x[f"{parent}_n_a"].to_numpy(),x[f"{parent}_k_a"].to_numpy(),
            x[f"{stage}_n_a"].to_numpy(),
            x[f"{parent}_n_b"].to_numpy(),x[f"{parent}_k_b"].to_numpy(),
            x[f"{stage}_n_b"].to_numpy(),rng=rng,draws=draws)
        x[f"{stage}_support_residual"] = observed-reference
    x["extra_genus_gap_beyond_support"] = x.genus_support_residual
    return x[outcols].copy()


def safe_aggregate(scored: pd.DataFrame, *, min_pairs: int = 20,
                   min_islands: int = 10) -> pd.DataFrame:
    """No island/name/genus identifiers in the public output."""
    if min_pairs < 20 or min_islands < 10:
        raise ValueError("Privacy minimum: >=20 pairs and >=10 different islands")
    fields = ["scope","region","measure","n_pairs","n_distinct_islands",
              "mean_excess_gap","median_excess_gap","fraction_pairs_excess_positive"]
    if scored.empty: return pd.DataFrame(columns=fields)
    rows = []
    def append_group(x: pd.DataFrame, label: str) -> None:
        if x.empty: return
        n_islands = len(set(x.island_a)|set(x.island_b))
        if len(x) < min_pairs or n_islands < min_islands: return
        v = x.extra_genus_gap_beyond_support
        rows.append({"scope":str(x.scope.iloc[0]),"region":str(x.region.iloc[0]),
                     "measure":label,"n_pairs":len(x),"n_distinct_islands":n_islands,
                     "mean_excess_gap":float(v.mean()),"median_excess_gap":float(v.median()),
                     "fraction_pairs_excess_positive":float((v>0).mean())})
    def domains(d: pd.DataFrame) -> None:
        for name in DOMAINS:
            sub = (d[d.trait.isin(TRAIT_GROUPS[name])]
                   .groupby(["scope","region","island_a","island_b"],as_index=False)
                   .agg(components=("trait","nunique"),
                        extra_genus_gap_beyond_support=("extra_genus_gap_beyond_support","mean")))
            append_group(sub[sub.components>=2], name)
    for (scope,region), chunk in scored.groupby(["scope","region"],sort=True):
        for trait in TRAITS: append_group(chunk[chunk.trait.eq(trait)],trait)
        domains(chunk)
    pooled = scored.copy()
    pooled["region"]="ALL_REGIONS"
    for scope,chunk in pooled.groupby("scope",sort=True):
        domains(chunk)
    return pd.DataFrame.from_records(rows,columns=fields)


def main(argv: list[str] | None = None) -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-counts",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--min-carriers",type=int,default=20)
    p.add_argument("--draws",type=int,default=399)
    p.add_argument("--seed",type=int,default=20261110)
    x=p.parse_args(argv)
    inputs=pd.read_csv(x.private_counts,dtype={k:str for k in KEYS})
    scored=score_pairs(inputs,min_carriers=x.min_carriers,draws=x.draws,seed=x.seed)
    summary=safe_aggregate(scored)
    x.output.mkdir(parents=True,exist_ok=True)
    summary.to_csv(x.output/"functional_genus_support_null_aggregate.csv",index=False)
    (x.output/"status.json").write_text(json.dumps({
        "status":"descriptive_sequential_thinning_not_confirmatory_inference",
        "private_input_exported":False,"n_public_group_rows":len(summary),
        "min_carriers":x.min_carriers,"null_draws":x.draws,
        "source_lineage_history_identified":False},indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
