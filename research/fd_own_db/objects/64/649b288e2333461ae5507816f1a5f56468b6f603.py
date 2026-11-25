"""True same-accepted-species THREE-trait sensitivity for island floral assembly.

Radial x deep and tubular/bell x deep are predicted on EXACT same island
target plants, using EXACT same other-species mainland source donors.
Each plant has independently directly recorded symmetry, form and tube depth.
Conditioning on recorded island genera does not predict colonization, and
source botanical records are not an independent field-phenotype provider.
"""
from __future__ import annotations

from collections import Counter, defaultdict

from .source_loo_flower_placebo import heldout_source_placebo
from hashlib import sha256

import numpy as np
import pandas as pd

SECTORS = ("actinomorphic_deep", "tubular_bell_deep")
SCOPES = ("all_recorded", "strict_native", "strict_native_nonendemic")
MODES = ("geo5", "climate10")
GATES = ("any_congener", "at_least_two_congeners")
STRATA = ("triple_direct", "triple_direct_deep_only")
REGIMES = ("northern_midlatitude", "tropical", "southern_extratropical",
           "northern_high_latitude", "geography_unresolved")
MIN_TARGET = 5
MIN_GENERA = 2
MIN_SOURCE = 10
MIN_ISLANDS = 30
MIN_BLOCKS = 10


def genus(sp: str) -> str:
    parts = str(sp).split()
    return parts[0] if len(parts) >= 2 else "unparsed_" + str(sp)


def _donor_genera(source: set[str]) -> dict[str, set[str]]:
    out = defaultdict(set)
    for sp in source:
        out[genus(sp)].add(sp)
    return dict(out)


def _empty(n: int, genera: int, target_coverage: float) -> dict:
    row = dict(n_target=n, n_target_genera=genera,
               island_direct_coverage=target_coverage)
    for mode in MODES:
        row.update({f"{mode}__n_donor": 0,
                    f"{mode}__donor_direct_coverage": np.nan})
        for metric in ("gain_radial", "gain_tubular", "gain_delta",
                       "relative_gain_radial", "relative_gain_tubular",
                       "relative_gain_delta", "prevalence_radial",
                       "prevalence_tubular", "pooled_probability_radial",
                       "pooled_probability_tubular", "brier_pooled_radial",
                       "brier_pooled_tubular",
                       "island_gain_radial","island_gain_tubular",
                       "source_loo_gain_radial","source_loo_gain_tubular",
                       "island_gain_difference","source_loo_gain_difference",
                       "island_excess_over_source_loo"):
            row[f"{mode}__{metric}"] = np.nan
    return row


def triple_matched_island(
    *, full_island: set[str], scope_rosters: dict[str, set[str]],
    geo5_source: set[str], climate10_source: set[str],
    profiles: dict,
) -> list[dict]:
    """Return 12 scope x donor gate x deep stratum source-paired rows.

    The source candidate accepted species may not appear anywhere in the
    focal island flora. The strict triply Direct phenotype subset is the
    *same* accepted target and donor species for both morphology endpoints.
    Both source modes share the exact same target accepted species.
    """
    if any(not r <= full_island for r in scope_rosters.values()):
        raise ValueError("native scope not a subset of original island flora")
    if (geo5_source | climate10_source) & full_island:
        raise ValueError("island species leaked into source-only donor pool")
    if not set(SECTORS) <= set(profiles) or not set(SCOPES) <= set(scope_rosters):
        raise ValueError("true paired floral profile or provenance missing")

    radial, tubular = (profiles[k] for k in SECTORS)
    triple = set(radial["known"]) & set(tubular["known"])
    # Deep has the exact same Direct ontology on both outcomes.
    if set(radial["B"]) & triple != set(tubular["B"]) & triple:
        raise AssertionError("true directly recorded tube-depth states disagree")
    known_deep = triple & set(radial["B"])
    source_all = {"geo5": set(geo5_source), "climate10": set(climate10_source)}
    # Cache original source genus opportunities ONCE, not 24 full-source
    # scans per island across status/gate/depth/endpoint sensitivity.
    source_genus_counts = {m: Counter(map(genus,v)) for m,v in source_all.items()}
    source_triple_counts = {m: Counter(map(genus,v & triple)) for m,v in source_all.items()}
    donor_by_stratum = {}
    for stratum in STRATA:
        scored = triple if stratum == "triple_direct" else known_deep
        donors = {mode: source_all[mode] & scored for mode in MODES}
        donor_by_stratum[stratum] = (scored, {
            mode: _donor_genera(donors[mode]) for mode in MODES
        })

    scope_genus_counts = {s: Counter(map(genus,roster)) for s,roster in scope_rosters.items()}
    scope_triple_counts = {s: Counter(map(genus,roster & triple)) for s,roster in scope_rosters.items()}
    rows = []
    for scope in SCOPES:
        roster = scope_rosters[scope]
        for gate, threshold in (("any_congener", 1),
                                ("at_least_two_congeners", 2)):
            for stratum in STRATA:
                scored, donor = donor_by_stratum[stratum]
                both = {
                    g for g, members in donor["geo5"].items()
                    if len(members) >= threshold
                } & {
                    g for g, members in donor["climate10"].items()
                    if len(members) >= threshold
                }
                targets = sorted(sp for sp in (roster & scored)
                                 if genus(sp) in both)
                observed_genera = {genus(sp) for sp in targets}
                # Counts of accepted plant names in matched genera, whether
                # or not their THREE traits were directly documented.
                target_total = sum(scope_genus_counts[scope][g] for g in observed_genera)
                direct_total = sum(scope_triple_counts[scope][g] for g in observed_genera)
                target_coverage = (
                    direct_total / target_total if target_total else np.nan
                )
                row = dict(scope=scope, donor_gate=gate, stratum=stratum)
                row.update(_empty(len(targets), len(observed_genera),
                                  target_coverage))
                if not targets:
                    rows.append(row)
                    continue

                target_by_g = defaultdict(list)
                for sp in targets:
                    target_by_g[genus(sp)].append(sp)
                for mode in MODES:
                    source_by_g = donor[mode]
                    matched_source = set().union(
                        *(source_by_g[g] for g in observed_genera)
                    )
                    source_unscored = sum(source_genus_counts[mode][g] for g in observed_genera)
                    source_direct = sum(source_triple_counts[mode][g] for g in observed_genera)
                    row[f"{mode}__n_donor"] = len(matched_source)
                    row[f"{mode}__donor_direct_coverage"] = (
                        source_direct/source_unscored if source_unscored else np.nan
                    )
                    for label, prof in (("radial", radial),
                                        ("tubular", tubular)):
                        positive = set(prof["joint"])
                        pooled_p = len(matched_source & positive) / len(matched_source)
                        y = np.asarray([sp in positive for sp in targets],float)
                        genus_p = {
                            g: len(source_by_g[g] & positive) / len(source_by_g[g])
                            for g in observed_genera
                        }
                        p = np.asarray([genus_p[genus(sp)] for sp in targets],
                                       dtype=float)
                        loss_pooled = float(np.square(y-pooled_p).mean())
                        loss_genus = float(np.square(y-p).mean())
                        gain = loss_pooled-loss_genus
                        row[f"{mode}__gain_{label}"] = gain
                        row[f"{mode}__relative_gain_{label}"] = (
                            gain/loss_pooled if loss_pooled > 0 else np.nan
                        )
                        row[f"{mode}__prevalence_{label}"] = float(y.mean())
                        row[f"{mode}__pooled_probability_{label}"] = pooled_p
                        row[f"{mode}__brier_pooled_{label}"] = loss_pooled
                    row[f"{mode}__gain_delta"] = (
                        row[f"{mode}__gain_radial"] -
                        row[f"{mode}__gain_tubular"]
                    )
                    if (threshold >= 2 and len(observed_genera) >= 2
                            and len(matched_source) >= 10):
                        source_placebo = heldout_source_placebo(
                            targets=targets,
                            donors_by_genus={
                                g:source_by_g[g] for g in observed_genera
                            }, radial=radial,tubular=tubular,
                        )
                        for name,value in source_placebo.items():
                            row[f"{mode}__{name}"]=value
                    r0 = row[f"{mode}__relative_gain_radial"]
                    r1 = row[f"{mode}__relative_gain_tubular"]
                    row[f"{mode}__relative_gain_delta"] = (
                        r0-r1 if np.isfinite(r0) and np.isfinite(r1)
                        else np.nan
                    )
                rows.append(row)
    if len(rows) != len(SCOPES)*len(GATES)*len(STRATA):
        raise AssertionError("missing full botanical island status panels")
    return rows


def aggregate_triple(frame: pd.DataFrame, *, population: int=4505,
                     draws: int=999) -> dict:
    """Anonymized group-bootstrap results, status-first when unsupported."""
    needed = {"island_id","regime","spatial_block","source_signature",
              "scope","donor_gate","stratum","n_target","n_target_genera",
              "island_direct_coverage"}
    for mode in MODES:
        needed |= {f"{mode}__n_donor",
                   f"{mode}__donor_direct_coverage"}
        for col in ("gain_radial","gain_tubular","gain_delta",
                    "relative_gain_delta","prevalence_radial",
                    "prevalence_tubular","brier_pooled_radial",
                    "brier_pooled_tubular"):
            needed.add(f"{mode}__{col}")
    if missing := needed-set(frame):
        raise ValueError(f"missing true triple-source fields: {sorted(missing)}")
    keys = ["island_id","scope","donor_gate","stratum"]
    if frame.duplicated(keys).any() or len(frame) != population*12 or (
        frame.island_id.nunique() != population
    ):
        raise ValueError("original botanical population/status completeness lost")
    if not frame.regime.isin(REGIMES).all() or (
        not frame.scope.isin(SCOPES).all()
        or not frame.donor_gate.isin(GATES).all()
        or not frame.stratum.isin(STRATA).all()
    ):
        raise ValueError("unknown original geographical or cohort category")
    if draws < 199:
        raise ValueError("bootstrap not enough iterations")

    results = []
    for region in REGIMES:
        for scope in SCOPES:
            for gate in GATES:
                for stratum in STRATA:
                    pool = frame.loc[
                        frame.regime.eq(region) & frame.scope.eq(scope) &
                        frame.donor_gate.eq(gate) & frame.stratum.eq(stratum)
                    ].copy()
                    # Both source modes must meet the identical prerequired
                    # candidate donor threshold, preserving matched islands.
                    good = (pool.n_target.ge(MIN_TARGET) &
                            pool.n_target_genera.ge(MIN_GENERA) &
                            pool.geo5__n_donor.ge(MIN_SOURCE) &
                            pool.climate10__n_donor.ge(MIN_SOURCE))
                    for mode in MODES:
                        good &= np.isfinite(pool[f"{mode}__gain_delta"])
                    eligible = pool.loc[good].copy()
                    balanced = eligible.island_direct_coverage.sub(
                        eligible.geo5__donor_direct_coverage).abs().le(.10) & (
                        eligible.island_direct_coverage.sub(
                            eligible.climate10__donor_direct_coverage).abs().le(.10))
                    for support, frame2 in (("all",eligible),
                                             ("doc_gap_le10pp",eligible.loc[balanced])):
                        for mode in MODES:
                            x = frame2.dropna(
                                subset=["spatial_block","source_signature"]
                            ).copy()
                            fields = [f"{mode}__gain_radial",
                                      f"{mode}__gain_tubular",
                                      f"{mode}__gain_delta",
                                      f"{mode}__prevalence_radial",
                                      f"{mode}__prevalence_tubular",
                                      f"{mode}__brier_pooled_radial",
                                      f"{mode}__brier_pooled_tubular"]
                            blocks = x.groupby("spatial_block",sort=True)[fields].mean()
                            sources = x.groupby("source_signature",sort=True)[fields].mean()
                            item = dict(
                                region=region, scope=scope, donor_gate=gate,
                                stratum=stratum, documentation=support, mode=mode,
                                n_original_region_islands=int(len(pool)),
                                n_eligible_both_modes=int(len(eligible)),
                                n_doc_balanced=int(balanced.sum()),
                                n_supported=int(len(x)),
                                n_blocks=int(len(blocks)),
                                n_first_sources=int(len(sources)),
                                median_triply_direct_targets=(
                                    float(x.n_target.median()) if len(x) else None
                                ),
                                median_source_triple_donors=(
                                    float(x[f"{mode}__n_donor"].median())
                                    if len(x) else None
                                ),
                            )
                            if (len(x)<MIN_ISLANDS or len(blocks)<MIN_BLOCKS
                                or len(sources)<MIN_BLOCKS):
                                item["status"]="NOT_IDENTIFIABLE_TRIPLE_DIRECT_SUPPORT"
                            else:
                                item["status"]="DESCRIPTIVE_TRIPLE_DIRECT_MATCHED_SOURCE"
                                for unit, groups in (("block",blocks),
                                                     ("first_source",sources)):
                                    arr = groups.to_numpy(dtype=float)
                                    vals = arr.mean(axis=0)
                                    seed=int.from_bytes(sha256(
                                        f"{region}|{scope}|{gate}|{stratum}|{support}|{mode}|{unit}".encode()
                                    ).digest()[:8],"big")
                                    rng=np.random.default_rng(seed)
                                    bootstrap=arr[rng.integers(
                                        len(arr),size=(draws,len(arr))
                                    )].mean(axis=1)
                                    for j,k in enumerate(fields):
                                        name=k.split("__",1)[1]
                                        item[f"{unit}__{name}"]=float(vals[j])
                                        item[f"{unit}__{name}_ci95"]=[
                                            float(v) for v in np.quantile(
                                                bootstrap[:,j],[.025,.975])
                                        ]
                                    # Relative Brier gain: divide both gains by
                                    # their OWN baseline risk. This is computed
                                    # after equal-block averaging, not by
                                    # selectively dropping zero-loss islands.
                                    if vals[5]<=0 or vals[6]<=0:
                                        item[f"{unit}__normalized_status"]=(
                                            "NOT_IDENTIFIABLE_ZERO_BASELINE_RISK")
                                    else:
                                        skill_r=vals[0]/vals[5]
                                        skill_t=vals[1]/vals[6]
                                        valid=(bootstrap[:,5]>0)&(bootstrap[:,6]>0)
                                        item[f"{unit}__n_normalized_bootstrap_valid"]=int(valid.sum())
                                        if valid.sum()<int(draws*.99):
                                            item[f"{unit}__normalized_status"]=(
                                                "NOT_IDENTIFIABLE_BOOTSTRAP_ZERO_BASELINE")
                                        else:
                                            norm_r=bootstrap[valid,0]/bootstrap[valid,5]
                                            norm_t=bootstrap[valid,1]/bootstrap[valid,6]
                                            item[f"{unit}__skill_radial"]=float(skill_r)
                                            item[f"{unit}__skill_tubular"]=float(skill_t)
                                            item[f"{unit}__skill_delta"]=float(skill_r-skill_t)
                                            item[f"{unit}__skill_delta_ci95"]=[
                                                float(v) for v in np.quantile(
                                                    norm_r-norm_t,[.025,.975])
                                            ]
                                            item[f"{unit}__normalized_status"]=(
                                                "SAME_COHORT_BASELINE_RISK_NORMALIZED")
                            results.append(item)
    return {
        "status":"ORIGINAL_ISLAND_SAME_SPECIES_THREE_DIRECT_TRAIT_SENSITIVITY",
        "n_original_botanical_islands":population,
        "n_original_scope_gate_stratum_rows":int(len(frame)),
        "n_panels":len(results),
        "panels":results,
        "limits":{
            "same_target_species_for_both_endpoints":True,
            "same_source_species_for_both_endpoints_within_source_mode":True,
            "source_donor_names_excluded_from_entire_focal_island_flora":True,
            "observed_genus_membership_oracle_not_colonization_model":True,
            "provider_and_archipelago_heldout":False,
            "effective_pollinator_or_adaptation_not_identified":True,
            "triple_direct_subset_only_secondary_not_global_population_filter":True,
        },
    }
