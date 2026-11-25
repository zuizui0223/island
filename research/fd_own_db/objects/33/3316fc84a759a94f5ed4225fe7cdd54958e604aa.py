"""Same-island Direct flower traits versus name-exclusive mainland congeners.

This is a SOURCE-TAXONOMY comparison of DIFFERENT accepted species in a shared
botanical genus. No observed seed arrival, ancestral phenotype, pollinator
service, historical shift, or island within-species adaptation is inferred.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from functools import lru_cache

import numpy as np
import pandas as pd

MORPHOLOGY = (
    "zygomorphic", "tube_or_bell", "open_radial", "deep_tube",
    "bilateral_tube_or_bell_joint",
)
SCOPES = (
    "all_recorded", "strict_native", "strict_native_nonendemic",
    "native_nonendemic_source_listed_both",
)
REGIMES = (
    "northern_midlatitude", "tropical", "southern_extratropical",
    "northern_high_latitude", "geography_unresolved",
)
MIN_ISLAND_DIRECT = 5
MIN_MODE_DIRECT_DONORS = 10
MIN_MATCHED_GENERA = 2
MIN_INDEPENDENT_ISLANDS = 30
MIN_SPATIAL_BLOCKS = 10


@lru_cache(maxsize=250000)
def genus(name: str) -> str:
    x = str(name).strip().split()
    if len(x) < 2:
        # A single unparsed name is not a trustworthy shared botanical genus.
        return "single_source_name:" + str(name)
    return x[0]


def _donor_map(names: set[str], phenotypes: dict[str, int]):
    total = Counter()
    documented = defaultdict(list)
    for sp in names:
        g = genus(sp)
        total[g] += 1
        if sp in phenotypes:
            documented[g].append((sp, int(phenotypes[sp])))
    return total, documented


def _subset_result(
    island_names: set[str], supported_genera: set[str],
    scope_names: set[str], species_states: dict[str, int],
    total_geo: Counter, direct_geo: dict,
    total_clim: Counter, direct_clim: dict,
) -> dict:
    documented = sorted(s for s in island_names
                        if genus(s) in supported_genera)
    all_g = [s for s in scope_names if genus(s) in supported_genera]
    counts = Counter(genus(s) for s in documented)
    n = len(documented)
    all_donor_geo = sum(total_geo[g] for g in counts)
    all_donor_clim = sum(total_clim[g] for g in counts)
    n_donor_geo = sum(len(direct_geo[g]) for g in counts)
    n_donor_clim = sum(len(direct_clim[g]) for g in counts)
    # Every original island accepted name was excluded from BOTH candidate
    # source pools, so no plant identity is reused as its own mainland donor.
    p_island = (float(np.mean([species_states[s] for s in documented]))
                if n else np.nan)
    p_geo = (float(sum(
        k * np.mean([z for _, z in direct_geo[g]]) for g, k in counts.items()
    ) / n) if n else np.nan)
    p_clim = (float(sum(
        k * np.mean([z for _, z in direct_clim[g]]) for g, k in counts.items()
    ) / n) if n else np.nan)
    return {
        "n_matched_island_direct": int(n),
        "n_matched_island_total": int(len(all_g)),
        "n_matched_genera": int(len(counts)),
        "n_geo5_shared_genus_candidate_direct": int(n_donor_geo),
        "n_climate10_shared_genus_candidate_direct": int(n_donor_clim),
        "n_geo5_shared_genus_candidate_total": int(all_donor_geo),
        "n_climate10_shared_genus_candidate_total": int(all_donor_clim),
        "observed_island_matched_fraction": p_island,
        "predicted_geo5_from_distinct_source_congeners": p_geo,
        "predicted_climate10_from_distinct_source_congeners": p_clim,
        "n_variable_source_genera_in_both": int(sum(
            len({v for _, v in direct_geo[g]}) >= 2 and
            len({v for _, v in direct_clim[g]}) >= 2
            for g in counts
        )),
    }


def _compare_with_original_source_maps(
    scope_names: set[str], species_states: dict[str, int],
    geo_total: Counter, geo_direct: dict,
    clim_total: Counter, clim_direct: dict,
) -> dict:
    all_scored = {sp for sp in scope_names if sp in species_states}
    eligible1 = {g for g in set(geo_direct) & set(clim_direct)
                 if geo_direct[g] and clim_direct[g]}
    eligible2 = {g for g in eligible1
                 if len(geo_direct[g]) >= 2 and len(clim_direct[g]) >= 2}
    main = _subset_result(all_scored, eligible1, scope_names, species_states,
                          geo_total, geo_direct, clim_total, clim_direct)
    strict = _subset_result(all_scored, eligible2, scope_names, species_states,
                            geo_total, geo_direct, clim_total, clim_direct)
    n_scope = len(all_scored)
    geo_vals = [v for group in geo_direct.values() for _, v in group]
    clim_vals = [v for group in clim_direct.values() for _, v in group]
    return {
        "n_scope_total": int(len(scope_names)),
        "n_scope_direct": int(n_scope),
        "unrestricted_island_direct_fraction": (
            float(np.mean([species_states[s] for s in all_scored]))
            if n_scope else np.nan
        ),
        "unrestricted_source_geo5_direct_fraction": (
            float(np.mean(geo_vals)) if geo_vals else np.nan
        ),
        "unrestricted_source_climate10_direct_fraction": (
            float(np.mean(clim_vals)) if clim_vals else np.nan
        ),
        "n_unrestricted_geo5_candidate_direct": len(geo_vals),
        "n_unrestricted_climate10_candidate_direct": len(clim_vals),
        "shared_genus_one_donor_each": main,
        "shared_genus_two_donors_each": strict,
    }


def compare_island_source_genera_scopes(
    *, full_island_names: set[str], scope_names: dict[str, set[str]],
    geo5_candidates: set[str], climate10_candidates: set[str],
    species_states: dict[str, int],
) -> dict[str, dict]:
    """Compute original genus references once, for fixed nested plant scopes.

    The exact SAME GIFT/Direct reference cohort is reused across scopes;
    each scope's *original island plant names* remain unaltered between modes.
    """
    if (geo5_candidates | climate10_candidates) & full_island_names:
        raise ValueError("leakage: donor accepted species already on focal island")
    if any(v not in (0, 1) for v in species_states.values()):
        raise ValueError("only genuinely Direct binary sector states permitted")
    if any(not names <= full_island_names for names in scope_names.values()):
        raise ValueError("scope is not part of original full botanical flora")
    gt,gd=_donor_map(geo5_candidates,species_states)
    ct,cd=_donor_map(climate10_candidates,species_states)
    return {
        scope: _compare_with_original_source_maps(
            names, species_states, gt, gd, ct, cd
        ) for scope, names in scope_names.items()
    }


def compare_island_source_genera(
    *, full_island_names: set[str], scope_names: set[str],
    geo5_candidates: set[str], climate10_candidates: set[str],
    species_states: dict[str, int],
) -> dict:
    """Single-scope adapter for tests and independent source checks."""
    return compare_island_source_genera_scopes(
        full_island_names=full_island_names,
        scope_names={"selected_original_scope": scope_names},
        geo5_candidates=geo5_candidates,
        climate10_candidates=climate10_candidates,
        species_states=species_states,
    )["selected_original_scope"]


def _support(x: pd.DataFrame, gate: str) -> pd.Series:
    return (
        x[f"{gate}__n_matched_island_direct"].ge(MIN_ISLAND_DIRECT)
        & x[f"{gate}__n_matched_genera"].ge(MIN_MATCHED_GENERA)
        & x[f"{gate}__n_geo5_shared_genus_candidate_direct"].ge(
            MIN_MODE_DIRECT_DONORS
        )
        & x[f"{gate}__n_climate10_shared_genus_candidate_direct"].ge(
            MIN_MODE_DIRECT_DONORS
        )
        & x[f"{gate}__observed_island_matched_fraction"].notna()
        & x[f"{gate}__predicted_geo5_from_distinct_source_congeners"].notna()
        & x[f"{gate}__predicted_climate10_from_distinct_source_congeners"].notna()
    )


def _cluster_summary(q: pd.DataFrame, *, group: str, key: str,
                     bootstrap: int) -> dict:
    m = q.dropna(subset=[group]).groupby(group, sort=True)
    size = m.ngroups
    out = {"n_islands": int(len(q)), "n_clusters": int(size)}
    if len(q) < MIN_INDEPENDENT_ISLANDS or size < MIN_SPATIAL_BLOCKS:
        return {**out, "status": "NOT_IDENTIFIABLE_CLUSTER_SUPPORT"}
    names = [
        "observed_minus_geo5_congener_reference",
        "observed_minus_climate10_congener_reference",
        "climate10_minus_geo5_congener_residual",
        "unrestricted_island_minus_geo5",
        "unrestricted_island_minus_climate10",
    ]
    a = m[names].mean().to_numpy(float)
    mean = np.mean(a, axis=0)
    seed = int.from_bytes(sha256(key.encode()).digest()[:8], "big")
    rng = np.random.default_rng(seed)
    boot = a[rng.integers(size, size=(bootstrap, size))].mean(axis=1)
    for ix, name in enumerate(names):
        out[name + "__mean"] = float(mean[ix])
        out[name + "__CI95"] = [float(z) for z in
                                np.quantile(boot[:, ix], [.025, .975])]
    out["status"] = "ORIGINAL_CONGENERIC_SOURCE_DESCRIPTIVE_NOT_CAUSAL"
    return out


def summarize_source_genus_frame(
    frame: pd.DataFrame, *, population: int = 4505, bootstrap: int = 999,
) -> dict:
    expected = {"island_id", "regime", "scope", "sector", "spatial_block",
                "source_signature", "n_scope_total", "n_scope_direct",
                "n_unrestricted_geo5_candidate_direct",
                "n_unrestricted_climate10_candidate_direct",
                "unrestricted_island_direct_fraction",
                "unrestricted_source_geo5_direct_fraction",
                "unrestricted_source_climate10_direct_fraction"}
    if not expected <= set(frame):
        raise ValueError(f"missing original source genus fields {expected-set(frame)}")
    for gate in ("shared_genus_one_donor_each", "shared_genus_two_donors_each"):
        for name in (
            "n_matched_island_direct", "n_matched_island_total",
            "n_matched_genera", "n_geo5_shared_genus_candidate_direct",
            "n_climate10_shared_genus_candidate_direct",
            "n_geo5_shared_genus_candidate_total",
            "n_climate10_shared_genus_candidate_total",
            "observed_island_matched_fraction",
            "predicted_geo5_from_distinct_source_congeners",
            "predicted_climate10_from_distinct_source_congeners",
        ):
            if gate + "__" + name not in frame:
                raise ValueError(f"missing genus donor field {gate}__{name}")
    if frame.island_id.nunique() != population:
        raise ValueError("the 4505 original botanical islands must remain")
    if frame.duplicated(["island_id", "scope", "sector"]).any():
        raise ValueError("same island/scope/sector duplicated")
    if len(frame) != population * len(SCOPES) * len(MORPHOLOGY):
        raise ValueError("original botanical genus rows incomplete")
    if not frame.scope.isin(SCOPES).all() or not frame.sector.isin(MORPHOLOGY).all():
        raise ValueError("strata or morphology ontology changed")
    if not frame.regime.isin(REGIMES).all():
        raise ValueError("original region support changed")
    if bootstrap < 199:
        raise ValueError("insufficient bootstrap draws")

    records = []
    for region in REGIMES:
        for scope in SCOPES:
            for sector in MORPHOLOGY:
                x = frame.loc[
                    frame.regime.eq(region) & frame.scope.eq(scope)
                    & frame.sector.eq(sector)
                ].copy()
                one = x.loc[_support(x, "shared_genus_one_donor_each")].copy()
                two = x.loc[_support(x, "shared_genus_two_donors_each")].copy()
                r = {
                    "region": region, "plant_provenance_scope": scope,
                    "floral_sector": sector,
                    "n_original_floras_in_region": int(len(x)),
                    "n_original_islands_congeneric_both_sources": int(len(one)),
                    "n_original_islands_min2_donor_species_per_genus": int(len(two)),
                    "n_original_islands_native_and_source_status_not_inferred": True,
                }
                for label, q in (("one_donor", one), ("two_donors", two)):
                    if len(q):
                        prefix = ("shared_genus_one_donor_each" if
                                  label == "one_donor" else
                                  "shared_genus_two_donors_each")
                        q["observed_minus_geo5_congener_reference"] = (
                            q[prefix + "__observed_island_matched_fraction"] -
                            q[prefix + "__predicted_geo5_from_distinct_source_congeners"]
                        )
                        q["observed_minus_climate10_congener_reference"] = (
                            q[prefix + "__observed_island_matched_fraction"] -
                            q[prefix + "__predicted_climate10_from_distinct_source_congeners"]
                        )
                        q["climate10_minus_geo5_congener_residual"] = (
                            q.observed_minus_climate10_congener_reference -
                            q.observed_minus_geo5_congener_reference
                        )
                        q["unrestricted_island_minus_geo5"] = (
                            q.unrestricted_island_direct_fraction -
                            q.unrestricted_source_geo5_direct_fraction
                        )
                        q["unrestricted_island_minus_climate10"] = (
                            q.unrestricted_island_direct_fraction -
                            q.unrestricted_source_climate10_direct_fraction
                        )
                        q["shared_genus_direct_carrier_fraction"] = (
                            q[prefix+"__n_matched_island_direct"] /
                            q.n_scope_direct.clip(lower=1)
                        )
                        q["island_genus_carrier_coverage"] = (
                            q[prefix+"__n_matched_island_direct"] /
                            q[prefix+"__n_matched_island_total"].replace(0, np.nan)
                        )
                        q["source_geo_genus_carrier_coverage"] = (
                            q[prefix+"__n_geo5_shared_genus_candidate_direct"] /
                            q[prefix+"__n_geo5_shared_genus_candidate_total"].replace(0, np.nan)
                        )
                        q["source_clim_genus_carrier_coverage"] = (
                            q[prefix+"__n_climate10_shared_genus_candidate_direct"] /
                            q[prefix+"__n_climate10_shared_genus_candidate_total"].replace(0, np.nan)
                        )
                        balanced = (
                            (q.island_genus_carrier_coverage -
                             q.source_geo_genus_carrier_coverage).abs().le(.10)
                            & (q.island_genus_carrier_coverage -
                               q.source_clim_genus_carrier_coverage).abs().le(.10)
                        )
                        r[label + "__median_island_direct_carrier_fraction"] = float(
                            q.shared_genus_direct_carrier_fraction.median()
                        )
                        r[label + "__n_coverage_gap_both_le10pp"] = int(balanced.sum())
                        r[label + "__n_positive_variable_donor_genera_both"] = int(
                            (q[prefix+"__n_variable_source_genera_in_both"] > 0).sum()
                        )
                    else:
                        balanced = pd.Series(False, index=q.index)
                    for scope_name, sub in (
                        ("full_matched_sample", q),
                        ("matched_direct_documentation_gap_le10pp", q.loc[balanced]),
                    ):
                        for group in ("spatial_block", "source_signature"):
                            r[f"{label}__{scope_name}__{group}"] = _cluster_summary(
                                sub, group=group, bootstrap=bootstrap,
                                key=f"{region}|{scope}|{sector}|{label}|{scope_name}|{group}"
                            )
                records.append(r)
    return {
        "status": "ORIGINAL_4505_GIFT_SOURCE_MODE_SHARED_GENUS_REAL_NAME_COMPARISON",
        "n_original_botanical_island_floras": int(population),
        "n_original_island_scope_morphology_rows": int(len(frame)),
        "n_regime_scope_morphology_panels": int(len(records)),
        "panels": records,
        "claim_limits": {
            "genuine_distinct_accepted_mainland_congener_species_only": True,
            "same_original_island_species_across_source_modes": True,
            "incomplete_source_genus_opportunity_is_not_no_evolution": True,
            "within_genus_not_independent_sister_ancestor": True,
            "no_plant_pollinator_causality_or_historical_guild_loss": True,
            "no_morphological_plant_population_phototypes": True,
            "phylogenetic_shared_genus_and_source_donor_overlap_uncorrected": True,
            "all_4505_botanical_islands_retained_even_when_unassessable": True,
        },
    }
