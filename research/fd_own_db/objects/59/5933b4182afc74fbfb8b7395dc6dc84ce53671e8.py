"""Pollinator-channel *opportunity* and island floral functional-diversity audit.

The channel exposure is qualified source-available GBIF detection versus
adequate NON-detection, NOT a time series of decline. Responses are fixed
species-level recorded morphological phenotypes, NOT measured pollen service.
Never export island/taxon/source-row identifiers or label absent observations
as lost guilds. Source GIFT lists are candidate plant pools, not ancestors.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from math import comb, log
import hashlib
import numpy as np
import pandas as pd

TRAITS = ("flower_primary_color", "floral_form", "floral_symmetry")
GUILDS = ("bombus", "non_bombus_bees", "lepidoptera",
          "flower_visiting_birds", "diptera")
METRICS = ("marginal_rao", "joint_rao", "joint_hill2",
           "colour_architecture_excess_mi", "open_radial_share",
           "actinomorphic_share")
GEOG = ("island_id", "analysis_regime", "spatial_block")
CHANNEL = ("island_id", "channel_id", "channel_state", "source_state")
DIRECT = ("accepted_species", "trait_name", "normalized_value",
          "evidence_scope", "evidence_quality")


def exact_colour_architecture_excess(sample: np.ndarray) -> float:
    """I(colour;form,symmetry) minus EXACT finite-30 fixed-margin expectation.

    In a random colour-label permutation, each colour/architecture cell count
    is hypergeometric. Expectation of mutual information is additive across
    cells, so no Monte Carlo shuffles are required. Architecture pairs stay
    intact. This is NOT a genus/source-conditioned null.
    """
    n = len(sample)
    if n < 2:
        raise ValueError("At least two known complete profiles required")
    colours = Counter(sample[:, 0])
    architecture = Counter(zip(sample[:, 1], sample[:, 2]))
    cells = Counter((r[0], r[1], r[2]) for r in sample)
    observed = sum(
        (v / n) * log(v * n / (colours[c] * architecture[(f,s)]))
        for (c,f,s), v in cells.items() if v
    )
    expected = 0.0
    for a in colours.values():
        for b in architecture.values():
            den = comb(n,b)
            for k in range(max(1, a+b-n), min(a,b)+1):
                probability = comb(a,k) * comb(n-a,b-k) / den
                expected += probability * (k/n) * log(k*n/(a*b))
    return float(observed-expected)


def sample_morphology(arr: np.ndarray, *, repeats: int, rng: np.random.Generator,
                      standardized_n: int = 30) -> dict[str,float]:
    """Equal-30 species outcome, without fabricated colour×form expansions."""
    if len(arr) < standardized_n:
        raise ValueError("Insufficient fully Direct plant species")
    scores = []
    for _ in range(repeats):
        s = arr[rng.choice(len(arr), size=standardized_n, replace=False)]
        marginal = []
        for j in range(3):
            counts = np.asarray(list(Counter(s[:, j]).values()), float)
            marginal.append(1-float(np.square(counts/standardized_n).sum()))
        combos = np.asarray(list(Counter(map(tuple,s)).values()),float)
        same = float(np.square(combos/standardized_n).sum())
        scores.append((
            float(np.mean(marginal)),
            1-same,
            1/same,
            exact_colour_architecture_excess(s),
            float(np.mean(s[:,1] == "open_radial")),
            float(np.mean(s[:,2] == "actinomorphic")),
        ))
    return dict(zip(METRICS, np.mean(scores, axis=0).astype(float)))


def direct_triple(ledger: pd.DataFrame) -> pd.DataFrame:
    if not set(DIRECT).issubset(ledger):
        raise ValueError("Missing original source Direct ledger contract")
    x = ledger.loc[
        ledger.trait_name.isin(TRAITS) &
        ledger.evidence_scope.eq("species_direct") &
        ledger.evidence_quality.isin(("high","medium")) &
        ledger.normalized_value.notna() &
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False),
        list(DIRECT)].copy()
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original species/trait duplicate or conflicting evidence")
    w = x.pivot(index="accepted_species", columns="trait_name", values="normalized_value")
    if not set(TRAITS).issubset(w):
        raise ValueError("Not all source Direct morphological traits present")
    return w.dropna(subset=list(TRAITS))[list(TRAITS)]


def _private_seed(island: str, scope: str, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}:{island}:{scope}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def _block_ci(group: pd.DataFrame, key: str, *, draws: int, seed: int):
    # Entire physical blocks, NOT 928 independent island replicates.
    blocks = group.groupby("spatial_block")[key].mean().dropna().to_numpy(float)
    if len(blocks) < 5:
        return (np.nan, np.nan, np.nan)
    rng = np.random.default_rng(seed)
    bs = blocks[rng.integers(len(blocks), size=(draws,len(blocks)))].mean(axis=1)
    return float(blocks.mean()), float(np.quantile(bs,.025)), float(np.quantile(bs,.975))


def audit_original_direct_channel_opportunity(
    original_direct: pd.DataFrame, original_flora: pd.DataFrame,
    geography: pd.DataFrame, channel_states: pd.DataFrame,
    source_assignments: pd.DataFrame, mainland_species: pd.DataFrame,
    *, repeats: int = 8, bootstrap: int = 1999,
    standardized_n: int = 30, seed: int = 20261113,
) -> dict:
    """Return ONLY anonymized aggregates; private scores stay in local memory.

    An *associational*, not causal, channel-region comparison is unsupported
    unless BOTH retained and
    qualified-disrupted states have >=10 islands in >=10 spatial blocks EACH,
    AND at least 3 distinct nearest-GIFT-source signatures have both states.
    This predeclared gate is independent of measured floral diversity.
    """
    if standardized_n != 30 or repeats < 5 or bootstrap < 999:
        raise ValueError("Nonstandardized support/randomization contract")
    for tab, cols, label in (
        (original_flora,("island_id","accepted_species","origin_status"),"original flora"),
        (geography, GEOG, "source geography"),
        (channel_states, CHANNEL, "GBIF channel evidence"),
        (source_assignments,("island_id","source_mode","source_rank","entity_ID"),"GIFT source"),
        (mainland_species,("entity_ID","work_species"),"GIFT native source plants")):
        if not set(cols).issubset(tab):
            raise ValueError(label+": required source fields absent")
    x = direct_triple(original_direct)
    flora = original_flora.loc[original_flora.accepted_species.isin(x.index),
                               ["island_id","accepted_species","origin_status"]].copy()
    if flora.isna().any().any() or flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Source flora missing or duplicate plant/island")
    if not flora.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Plant-origin label not audited")
    eligible = flora.groupby("island_id").accepted_species.nunique()
    eligible = set(eligible[eligible >= 30].index)
    if len(eligible) != 928:
        raise ValueError("Original locked 928 Direct-island universe drifted")
    flora = flora.loc[flora.island_id.isin(eligible)]
    geo = geography[list(GEOG)].drop_duplicates()
    if geo.isna().any().any() or geo.duplicated("island_id").any() or not eligible.issubset(set(geo.island_id)):
        raise ValueError("Geographic block/regime unmapped")
    gmap = geo.set_index("island_id")
    states = channel_states.loc[channel_states.island_id.isin(eligible), list(CHANNEL)].copy()
    if states.isna().any().any() or states.duplicated(["island_id","channel_id"]).any():
        raise ValueError("Invalid source GBIF channel input")
    if not states.channel_id.isin(GUILDS).all() or not states.channel_state.isin(
       ("retained","disrupted","unresolved","structurally_absent")).all():
        raise ValueError("Invalid qualified GBIF guild/state")
    if ((states.channel_state.eq("disrupted")) &
        (~states.source_state.eq("available"))).any():
        raise ValueError("No qualified disruption without independent source opportunity")
    channel_table = states.pivot(index="island_id",columns="channel_id",values="channel_state")
    if set(channel_table.index)!=eligible or set(channel_table.columns)!=set(GUILDS) or channel_table.isna().any().any():
        raise ValueError("Missing one of five original island guild states")
    donors = source_assignments.loc[
        source_assignments.island_id.isin(eligible) &
        source_assignments.source_mode.eq("geo_k5") &
        source_assignments.source_rank.between(1,5),
        ["island_id","source_rank","entity_ID"]].copy()
    if donors.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicated GIFT source rank")
    source_map = donors.sort_values(["island_id","source_rank"]).groupby("island_id").entity_ID.agg(tuple).to_dict()
    if any(len(source_map.get(i,()))!=5 or len(set(source_map[i]))!=5 for i in eligible):
        raise ValueError("Incomplete fixed five-nearest GIFT source opportunities")
    mainland = mainland_species.loc[
        mainland_species.work_species.isin(x.index),
        ["entity_ID","work_species"]].dropna().drop_duplicates()
    gift_map = mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    flora_map = {i:g for i,g in flora.groupby("island_id")}
    rows = []
    for island in sorted(eligible):
        g = flora_map[island]
        raw = x.loc[g.accepted_species.tolist()].to_numpy(dtype=str)
        src = set().union(*(gift_map.get(r,set()) for r in source_map[island]))
        src &= set(x.index)
        donor_arr = x.loc[sorted(src)].to_numpy(dtype=str) if len(src)>=30 else None
        source_vals = None
        if donor_arr is not None:
            source_vals = sample_morphology(donor_arr,repeats=repeats,
                rng=np.random.default_rng(_private_seed(island,"GIFT",seed)))
        for scope in ("all_recorded","source_native"):
            subset = g if scope=="all_recorded" else g[g.origin_status.eq("native")]
            if len(subset)<30:
                continue
            scores = sample_morphology(x.loc[subset.accepted_species.tolist()].to_numpy(dtype=str),
                repeats=repeats,rng=np.random.default_rng(_private_seed(island,scope,seed)))
            d = dict(island_id=island,scope=scope,
                analysis_regime=str(gmap.loc[island,"analysis_regime"]),
                spatial_block=str(gmap.loc[island,"spatial_block"]),
                source_signature=source_map[island][0],
                n_complete_species=len(subset),n_source_complete_species=len(src),
                bombus=channel_table.loc[island,"bombus"],
                flower_visiting_birds=channel_table.loc[island,"flower_visiting_birds"])
            for guild in GUILDS:
                d[guild+"_state"] = channel_table.loc[island,guild]
            for metric in METRICS:
                d[metric]=scores[metric]
                d["source_delta_"+metric]=(scores[metric]-source_vals[metric]) if source_vals is not None else np.nan
            rows.append(d)
    private=pd.DataFrame.from_records(rows)
    # Only anonymous region/status and block support; never write 'private'.
    region=[]
    for (scope,regime),g in private.groupby(["scope","analysis_regime"]):
        if len(g)<10 or g.spatial_block.nunique()<5:
            continue
        obj=dict(scope=scope,regime=regime,n_islands=len(g),
                 n_spatial_blocks=g.spatial_block.nunique(),
                 n_nearest_mainland_source_signatures=g.source_signature.nunique(),
                 n_source_comparable=int(g["source_delta_marginal_rao"].notna().sum()))
        for j,metric in enumerate(METRICS):
            obj["median_"+metric]=float(g[metric].median())
            key="source_delta_"+metric
            sub=g.loc[g[key].notna()]
            obj["median_source_delta_"+metric]=float(sub[key].median()) if len(sub) else None
            # Two separate conservative descriptive uncertainty levels.
            # Source signature and spatial block represent overlapping, NOT
            # phylogenetically independent historical island replicates.
            if len(sub)>=10 and sub.spatial_block.nunique()>=5:
                mean,lo,hi=_block_ci(sub,key,draws=bootstrap,seed=seed+100+j)
                obj["block_equal_mean_source_delta_"+metric]=mean
                obj["block_bootstrap_lower_"+metric]=lo
                obj["block_bootstrap_upper_"+metric]=hi
            if len(sub)>=10 and sub.source_signature.nunique()>=5:
                src_means=sub.groupby("source_signature")[key].mean().to_numpy(float)
                rng=np.random.default_rng(seed+200+j)
                rs=src_means[rng.integers(len(src_means),
                    size=(bootstrap,len(src_means)))].mean(axis=1)
                obj["source_equal_mean_source_delta_"+metric]=float(src_means.mean())
                obj["source_bootstrap_lower_"+metric]=float(np.quantile(rs,.025))
                obj["source_bootstrap_upper_"+metric]=float(np.quantile(rs,.975))
        region.append(obj)
    supports=[]
    for guild in GUILDS:
        for regime,g in private.loc[private.scope.eq("all_recorded")].groupby("analysis_regime"):
            vals=g[g[guild+"_state"].isin(("retained","disrupted"))]
            a=vals[vals[guild+"_state"].eq("retained")]
            b=vals[vals[guild+"_state"].eq("disrupted")]
            # Nearest source overlap tests report only support cardinalities.
            common=set(a.source_signature)&set(b.source_signature)
            status=(len(a)>=10 and len(b)>=10 and
                    a.spatial_block.nunique()>=10 and b.spatial_block.nunique()>=10 and
                    len(common)>=3)
            supports.append(dict(guild=guild,regime=regime,
                n_retained=len(a),n_adequately_nondetected=len(b),
                retained_spatial_blocks=a.spatial_block.nunique(),
                nondetected_spatial_blocks=b.spatial_block.nunique(),
                nearest_mainland_source_signatures_with_both=len(common),
                source_comparable_retained=int(a["source_delta_marginal_rao"].notna().sum()),
                source_comparable_nondetected=int(b["source_delta_marginal_rao"].notna().sum()),
                qualified_observational_contrast_admissible=bool(status)))
    # Scientific candidate only; BOTH channels are documentary opportunity states,
    # NOT bird-pollinated flowers nor an observed historical shift.
    candidates=[]
    for scope,g in private.groupby("scope"):
        focus=g[g.bombus.eq("disrupted") & g.flower_visiting_birds.eq("retained")]
        focus=focus[focus["source_delta_marginal_rao"].notna()]
        if len(focus)<10 or focus.spatial_block.nunique()<5:
            continue
        rec=dict(scope=scope,n_candidate_islands=len(focus),
                 n_spatial_blocks=focus.spatial_block.nunique(),
                 direct_species_fixed_across_locations=True,
                 confirmed_pollinator_shifts=0)
        for j,metric in enumerate(METRICS):
            vals=focus["source_delta_"+metric]
            mean,lower,upper=_block_ci(focus,"source_delta_"+metric,
                                      draws=bootstrap,seed=seed+j)
            rec["median_source_delta_"+metric]=float(vals.median())
            rec["block_equal_mean_source_delta_"+metric]=mean
            rec["block_bootstrap_lower_"+metric]=lower
            rec["block_bootstrap_upper_"+metric]=upper
        candidates.append(rec)
    # Strict native sensitivity MUST compare the SAME ISLANDS before attempting
    # an interpretation; the 67 eligible native islands are a selected subset
    # of all 743 northern midlatitude islands, not an independent sample.
    paired_native=[]
    all_scope=private.loc[private.scope.eq("all_recorded")].set_index("island_id")
    nat_scope=private.loc[private.scope.eq("source_native")].set_index("island_id")
    for regime,native in nat_scope.groupby("analysis_regime"):
        common=native.index.intersection(all_scope.index)
        if len(common)<10 or native.spatial_block.nunique()<5:
            continue
        same_all=all_scope.loc[common]
        equal_native=native.loc[common]
        rec=dict(regime=regime,n_paired_islands=len(common),
                 n_paired_spatial_blocks=native.spatial_block.nunique(),
                 same_islands_and_fixed_30_species_each=True,
                 not_an_experiment_on_species_introduction=True)
        frame=equal_native[["spatial_block","source_signature"]].copy()
        for j,metric in enumerate(METRICS):
            rec["median_all_same_islands_source_delta_"+metric] = float(
                same_all["source_delta_"+metric].median())
            rec["median_native_same_islands_source_delta_"+metric] = float(
                equal_native["source_delta_"+metric].median())
            frame["paired_native_minus_all_"+metric]=(
                equal_native[metric]-same_all[metric]).to_numpy(float)
            mean,lo,hi=_block_ci(frame,"paired_native_minus_all_"+metric,
                                 draws=bootstrap,seed=seed+300+j)
            rec["block_mean_native_minus_all_"+metric]=mean
            rec["block_ci025_native_minus_all_"+metric]=lo
            rec["block_ci975_native_minus_all_"+metric]=hi
        paired_native.append(rec)
    return {
      "status":"OPPORTUNITY_AND_FUNCTIONAL_DIVERSITY_NOT_HISTORICAL_LOSS",
      "population_trait_evolution_tested":False,
      "real_local_visitor_shift_measured":False,
      "pollen_deposition_or_seedset_measured":False,
      "biological_guild_absence_from_unresolved":False,
      "global_source_is_not_historical_ancestor":True,
      "n_original_qualifying_islands":len(eligible),
      "n_island_scope_rows_privately_scored":len(private),
      "n_direct_triple_species":len(x),
      "n_direct_islands_with_native_30plus":int((private.scope=="source_native").sum()),
      "n_qualified_channel_region_contrasts_passed":sum(t["qualified_observational_contrast_admissible"] for t in supports),
      "support_by_guild_regime":supports,
      "floral_functional_assemblage_by_region":region,
      "source_relative_bombus_nondetected_bird_detected_candidates":candidates,
      "same_island_direct_native_status_sensitivity":paired_native,
      "sampling":{"equal_species_per_island_or_source":standardized_n,"rarefactions":repeats,
                  "exact_colour_label_shuffle_null":True,
                  "block_bootstrap_draws":bootstrap},
      "raw_island_taxon_or_source_identifiers_exported":False,
    }
