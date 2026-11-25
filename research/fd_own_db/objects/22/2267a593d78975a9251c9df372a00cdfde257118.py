"""Separate floral marginal diversity from nonrandom cross-trait pairing.

Only DIRECT, strictly single-state species evidence may enter. Analyses are
observational morphologies, not realized guild/pollination outcomes. The
returned summaries have NO island, taxon, genus or source IDs. Input remains
local/rights-restricted. Permutations preserve recorded island trait margins;
a within-genus permutation is a POWER/identifiability diagnostic, not a
phylogenetically independent causal null.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

TRAITS = ("flower_primary_color", "floral_form", "floral_symmetry")
TEMPLATES = ("large_bee_like", "butterfly_like", "bird_like")
REQUIRED = ("island_id", "accepted_species", "genus", "analysis_regime",
            "spatial_block", "origin_status", *TRAITS)
OUTPUT_COLUMNS = (
    "floristic_scope", "analysis_regime", "n_islands", "n_spatial_blocks",
    "median_trait_complete_species", "median_marginal_rao",
    "median_color_architecture_excess_nats", "frac_islands_color_arch_positive",
    "median_form_symmetry_excess_nats", "frac_islands_form_sym_positive",
    "median_color_form_excess_nats", "frac_islands_color_form_positive",
    "median_genus_shuffled_color_arch_excess_nats",
    "median_within_genus_color_shuffle_support",
    "equal_block_mean_color_arch_excess_nats",
    "block_bootstrap_lower", "block_bootstrap_upper",
)


def _entropy(codes: np.ndarray) -> float:
    if len(codes) == 0:
        raise ValueError("Empty categorical vector")
    counts = np.bincount(codes.astype(int))
    p = counts[counts > 0] / len(codes)
    return float(-(p * np.log(p)).sum())


def marginal_rao_equal_traits(traits: np.ndarray) -> float:
    """Equal-weight Gower categorical Rao Q, WITH replacement.

    Q=(1/3) sum_j (1-sum_s p_js^2). It sees trait marginals but is
    mathematically blind to colour/form/symmetry cross-trait pairings.
    """
    arr = np.asarray(traits)
    if arr.ndim != 2 or arr.shape[1] != 3 or len(arr) < 2:
        raise ValueError("Expected n x 3 complete categorical trait matrix")
    if not np.isfinite(arr.astype(float)).all() or (arr < 0).any():
        raise ValueError("Invalid categorical codes")
    n = len(arr)
    return float(np.mean([
        1 - np.sum((np.bincount(arr[:, j].astype(int)) / n) ** 2)
        for j in range(3)
    ]))


def _excess_one_sample(
    traits: np.ndarray, genera: np.ndarray,
    rng: np.random.Generator, *, permutations: int
) -> dict:
    """Return null-corrected information in nats (NOT a p-value)."""
    if len(traits) != len(genera) or traits.shape[1] != 3:
        raise ValueError("Misaligned genus and morphology vectors")
    c, f, s = traits.T.astype(int)
    nf, ns = int(f.max()) + 1, int(s.max()) + 1
    fs = f * ns + s
    h_c, h_f, h_s, h_fs = _entropy(c), _entropy(f), _entropy(s), _entropy(fs)
    raw_ca = h_c + h_fs - _entropy(c * (nf * ns) + fs)
    raw_cf = h_c + h_f - _entropy(c * nf + f)
    raw_fs = h_f + h_s - h_fs
    groups = [np.flatnonzero(genera == g) for g in np.unique(genera)]
    groups = [idx for idx in groups if len(idx) >= 2]
    info_n = sum(len(idx) for idx in groups if np.unique(c[idx]).size > 1)
    null_ca, null_cf, null_fs, null_ca_genus = [], [], [], []
    for _ in range(permutations):
        cc = rng.permutation(c)
        ss = rng.permutation(s)
        null_ca.append(h_c + h_fs - _entropy(cc * (nf * ns) + fs))
        null_cf.append(h_c + h_f - _entropy(cc * nf + f))
        null_fs.append(h_f + h_s - _entropy(f * ns + ss))
        within = c.copy()
        for idx in groups:
            within[idx] = rng.permutation(c[idx])
        null_ca_genus.append(h_c + h_fs -
                             _entropy(within * (nf * ns) + fs))
    return {
        "rao": marginal_rao_equal_traits(traits),
        "excess_ca": float(raw_ca - np.mean(null_ca)),
        "excess_cf": float(raw_cf - np.mean(null_cf)),
        "excess_fs": float(raw_fs - np.mean(null_fs)),
        "excess_ca_genus": float(raw_ca - np.mean(null_ca_genus)),
        "fraction_genus_color_informative": float(info_n / len(c)),
    }


def _check_direct(frame: pd.DataFrame) -> pd.DataFrame:
    if not set(REQUIRED).issubset(frame.columns):
        raise ValueError("Required DIRECT island x taxon trait fields missing")
    x = frame[list(REQUIRED)].copy()
    if x.empty:
        raise ValueError("No Direct morphological evidence")
    for col in REQUIRED:
        if x[col].isna().any() or x[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Unknown/missing {col}")
        x[col] = x[col].astype(str).str.strip()
    if x.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Duplicate island x accepted species")
    if not (x.accepted_species.str.split().str[0] == x.genus).all():
        raise ValueError("Accepted species and genus mismatch")
    for col in TRAITS:
        if x[col].str.contains("|", regex=False).any():
            raise ValueError("Mixed-state evidence is NOT a single Direct state")
        if x.groupby("accepted_species")[col].nunique().gt(1).any():
            raise ValueError("Conflicting species-direct trait state")
    for col in ("analysis_regime", "spatial_block"):
        if x.groupby("island_id")[col].nunique().gt(1).any():
            raise ValueError("Inconsistent island geography")
    return x


def audit_island_flora(
    direct_complete: pd.DataFrame,
    *, min_species: int = 30, rarefactions: int = 5,
    permutations: int = 49, min_islands: int = 20,
    min_spatial_blocks: int = 5, block_draws: int = 1999,
    seed: int = 20261111
) -> pd.DataFrame:
    """Return only aggregate region × observed/native summaries.

    Inputs MUST have been joined with original Direct, strict single-state
    species trait records and source-backed floristic-status evidence.
    Regions with thin support are suppressed. Bootstrap is conditional on
    observed island/spatial-block composition and shared species histories.
    """
    if (min_species < 30 or rarefactions < 1 or permutations < 29 or
        min_islands < 20 or min_spatial_blocks < 5 or block_draws < 199):
        raise ValueError("Privacy / precision gate not met")
    x = _check_direct(direct_complete)
    for col in TRAITS:
        x[col] = pd.Categorical(x[col]).codes.astype(int)
    rng = np.random.default_rng(seed)
    records = []
    for scope in ("all_recorded", "source_native"):
        part = x if scope == "all_recorded" else x[x.origin_status.eq("native")]
        for island, group in part.groupby("island_id", sort=True):
            n = len(group)
            if n < min_species:
                continue
            arrays = group[list(TRAITS)].to_numpy(dtype=int)
            genus = group.genus.to_numpy()
            draw_scores = []
            for _ in range(rarefactions):
                idx = (np.arange(n) if n == min_species
                       else rng.choice(n, min_species, replace=False))
                draw_scores.append(_excess_one_sample(
                    arrays[idx], genus[idx], rng, permutations=permutations
                ))
            avg = {key: float(np.mean([v[key] for v in draw_scores]))
                   for key in draw_scores[0]}
            records.append(dict(scope=scope, island_id=island,
                                region=group.analysis_regime.iloc[0],
                                spatial_block=group.spatial_block.iloc[0],
                                n_species=n, **avg))
    scored = pd.DataFrame(records)
    if scored.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    out = []
    for (scope, region), g in scored.groupby(["scope", "region"], sort=True):
        b = g.groupby("spatial_block")["excess_ca"].mean()
        if len(g) < min_islands or len(b) < min_spatial_blocks:
            continue
        values = b.to_numpy()
        replicated = values[rng.integers(len(values),
                          size=(block_draws, len(values)))].mean(axis=1)
        lo, hi = np.quantile(replicated, [.025, .975])
        out.append({
            "floristic_scope": scope, "analysis_regime": region,
            "n_islands": int(len(g)), "n_spatial_blocks": int(len(b)),
            "median_trait_complete_species": float(g.n_species.median()),
            "median_marginal_rao": float(g.rao.median()),
            "median_color_architecture_excess_nats": float(g.excess_ca.median()),
            "frac_islands_color_arch_positive": float((g.excess_ca > 0).mean()),
            "median_form_symmetry_excess_nats": float(g.excess_fs.median()),
            "frac_islands_form_sym_positive": float((g.excess_fs > 0).mean()),
            "median_color_form_excess_nats": float(g.excess_cf.median()),
            "frac_islands_color_form_positive": float((g.excess_cf > 0).mean()),
            "median_genus_shuffled_color_arch_excess_nats":
                float(g.excess_ca_genus.median()),
            "median_within_genus_color_shuffle_support":
                float(g.fraction_genus_color_informative.median()),
            "equal_block_mean_color_arch_excess_nats": float(values.mean()),
            "block_bootstrap_lower": float(lo), "block_bootstrap_upper": float(hi),
        })
    return pd.DataFrame.from_records(out, columns=OUTPUT_COLUMNS)


def template_overlap(
    direct_syndrome_scores: pd.DataFrame, *,
    min_informative: int = 3, high_cutoff: float = .75,
    min_species: int = 50
) -> dict:
    """Check ambiguity in three existing *candidate* syndromes.

    This does not classify actual pollinators. Require every named candidate
    to have >=min_informative source-backed traits for a given taxon.
    """
    required = {"accepted_species", "syndrome", "soft_membership",
                "n_informative_traits"}
    if not required.issubset(direct_syndrome_scores):
        raise ValueError("Missing syndrome score fields")
    if min_informative < 3 or not (.5 <= high_cutoff <= 1) or min_species < 50:
        raise ValueError("Unsafe syndrome overlap threshold")
    d = direct_syndrome_scores[list(required)].copy()
    if d.duplicated(["accepted_species", "syndrome"]).any():
        raise ValueError("Duplicate species x syndrome score")
    d = d[d.syndrome.isin(TEMPLATES) &
          (d.n_informative_traits >= min_informative)].copy()
    if not d.soft_membership.between(0, 1).all():
        raise ValueError("Invalid template soft-membership score")
    wide = d.pivot(index="accepted_species", columns="syndrome",
                   values="soft_membership").reindex(columns=TEMPLATES).dropna()
    if len(wide) < min_species:
        raise ValueError("Too few comparable species to release")
    h = (wide >= high_cutoff).sum(axis=1)
    return {
        "n_comparable_species": int(len(wide)),
        "n_high_two_or_more": int((h >= 2).sum()),
        "n_high_all_three": int((h == 3).sum()),
        "fraction_high_two_or_more": float((h >= 2).mean()),
        "fraction_high_all_three": float((h == 3).mean()),
        "min_informative_traits_each_template": min_informative,
        "high_cutoff": high_cutoff,
        "pollinator_identity_established": False,
    }
