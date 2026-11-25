"""Rights-safe *descriptive* floral trait turnover audit.

Private island x taxon x trait inputs remain in local memory. The CLI writes
ONLY group aggregates; no taxon, island or genus identifiers are exported.
No prediction, functional-service, or statistical-significance claim is made.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TRAIT_GROUPS: dict[str, tuple[str, ...]] = {
    "reproductive_assurance": (
        "self_compatibility", "selfing_mating_system", "autonomous_selfing"
    ),
    "floral_accessibility": (
        "generalized_form", "actinomorphic_symmetry", "shallow_open_tube"
    ),
    "floral_signal_proxy": ("plain_colour",),
}
TRAITS = tuple(t for values in TRAIT_GROUPS.values() for t in values)
REQUIRED = ("island_id", "taxon_name", "genus", *TRAITS)
COMPARISONS = ("all_taxa", "species_exclusive", "genus_exclusive")


def validate_trait_table(table: pd.DataFrame) -> pd.DataFrame:
    """Reject inferred zeros, conflicting taxon records, and duplicates."""
    missing = set(REQUIRED) - set(table.columns)
    if missing:
        raise ValueError(f"Required columns missing: {sorted(missing)}")
    frame = table[list(REQUIRED)].copy()
    for col in ("island_id", "taxon_name", "genus"):
        if frame[col].isna().any() or frame[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing or blank {col}")
        frame[col] = frame[col].astype(str).str.strip()
    if frame.duplicated(["island_id", "taxon_name"]).any():
        raise ValueError("Duplicate island x taxon rows")
    if frame.groupby("taxon_name")["genus"].nunique().gt(1).any():
        raise ValueError("Inconsistent genus for same taxon")
    for col in TRAITS:
        raw = frame[col]
        numeric = pd.to_numeric(raw, errors="coerce")
        nonblank = raw.notna() & raw.astype(str).str.strip().ne("")
        if (nonblank & numeric.isna()).any() or (~numeric.dropna().isin([0, 1])).any():
            raise ValueError(f"{col} must be 0, 1, or missing")
        frame[col] = numeric.astype(float)
        if frame.groupby("taxon_name")[col].nunique(dropna=True).gt(1).any():
            raise ValueError(f"Inconsistent {col} for the same taxon")
    return frame


def audit_pairs(
    traits: pd.DataFrame, pairs: pd.DataFrame, *, min_carriers: int = 20
) -> pd.DataFrame:
    """Observed taxonomic/trait distance, in three non-independent comparison scopes.

    The pairs input must contain island_a, island_b, comparison. Comparison is
    an outcome-blind predefined donor category; this function does not match pairs.
    Results contain private island keys and MUST NOT be published directly.
    """
    if min_carriers < 1:
        raise ValueError("min_carriers must be positive")
    required = {"island_a", "island_b", "comparison"}
    if not required.issubset(pairs.columns):
        raise ValueError(f"Pair columns missing: {sorted(required-set(pairs.columns))}")
    data = validate_trait_table(traits)
    by_island = {k: g.set_index("taxon_name", drop=False) for k, g in data.groupby("island_id")}
    records: list[dict] = []
    for pair in pairs.itertuples(index=False):
        a, b = str(pair.island_a), str(pair.island_b)
        if a == b:
            raise ValueError("Self-island comparison")
        if a not in by_island or b not in by_island:
            raise ValueError("Pair refers to island absent from trait table")
        aa, bb = by_island[a], by_island[b]
        sa, sb = set(aa.index), set(bb.index)
        ga, gb = set(aa["genus"]), set(bb["genus"])
        shared_taxa = sa & sb
        shared_genera = ga & gb
        union_n = len(sa | sb)
        species_jaccard = len(shared_taxa) / union_n if union_n else np.nan
        for scope in COMPARISONS:
            if scope == "all_taxa":
                fa, fb = aa, bb
            elif scope == "species_exclusive":
                fa = aa.loc[~aa.index.isin(shared_taxa)]
                fb = bb.loc[~bb.index.isin(shared_taxa)]
            else:
                fa = aa.loc[~aa["genus"].isin(shared_genera)]
                fb = bb.loc[~bb["genus"].isin(shared_genera)]
            for domain, trait_names in TRAIT_GROUPS.items():
                for t in trait_names:
                    pa, pb = fa[t].dropna(), fb[t].dropna()
                    records.append({
                        "island_a": a, "island_b": b,
                        "comparison": str(pair.comparison),
                        "trait_scope": scope, "domain": domain, "trait": t,
                        "n_a": len(pa), "n_b": len(pb),
                        "taxon_jaccard": species_jaccard,
                        "abs_prevalence_gap": (
                            abs(float(pa.mean()) - float(pb.mean()))
                            if len(pa) >= min_carriers and len(pb) >= min_carriers
                            else np.nan
                        ),
                    })
    return pd.DataFrame.from_records(records)


def safe_summary(rows: pd.DataFrame, *, minimum_pairs: int = 10) -> pd.DataFrame:
    """Suppress thin comparisons; never exports species or island identifiers."""
    if minimum_pairs < 5:
        raise ValueError("minimum_pairs must be at least 5")
    if rows.empty:
        return pd.DataFrame(columns=[
            "comparison", "trait_scope", "domain", "trait",
            "n_pairs", "median_abs_prevalence_gap"
        ])
    eligible = rows.dropna(subset=["abs_prevalence_gap"])
    result = eligible.groupby(
        ["comparison", "trait_scope", "domain", "trait"], as_index=False
    ).agg(
        n_pairs=("abs_prevalence_gap", "size"),
        median_abs_prevalence_gap=("abs_prevalence_gap", "median"),
    )
    return result.loc[result["n_pairs"] >= minimum_pairs].reset_index(drop=True)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--private-traits", type=Path, required=True)
    parser.add_argument("--private-pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-carriers", type=int, default=20)
    parser.add_argument("--minimum-pairs", type=int, default=10)
    args = parser.parse_args(argv)
    # Inputs MUST be rights-authorized, locally prepared, accepted-name and
    # evidence-scope harmonized; no public fetch or imputation takes place.
    traits = pd.read_csv(
        args.private_traits,
        dtype={"island_id": str, "taxon_name": str, "genus": str}
    )
    pairs = pd.read_csv(
        args.private_pairs, dtype={"island_a": str, "island_b": str}
    )
    rows = audit_pairs(traits, pairs, min_carriers=args.min_carriers)
    summary = safe_summary(rows, minimum_pairs=args.minimum_pairs)
    args.output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output / "functional_turnover_grouped.csv", index=False)
    (args.output / "status.json").write_text(json.dumps({
        "status": "descriptive_only_no_null_or_p_values",
        "private_input_released": False,
        "traits": list(TRAITS),
        "minimum_pairs": args.minimum_pairs,
        "minimum_carriers_per_island_and_trait": args.min_carriers,
        "released_groups": len(summary),
        "requires_additional_null_and_geo_independence_audit": True,
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
