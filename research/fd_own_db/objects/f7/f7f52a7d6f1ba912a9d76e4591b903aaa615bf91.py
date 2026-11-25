"""Within-observed-network guild-loss *simulation*, not observed pollinator decline.

Only original EuPPollNet v1.3.0 flower-visitor links at exact GSHHG islands
can support this analysis. Suppress network analyses where source sampling
recorded only one order/guild. A lost DOCUMENTED link does not prove loss of
effective pollen delivery or seed set. No raw island/source IDs may be emitted.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import numpy as np
import pandas as pd
from .euppollnet_island_overlap import validate_original_visits, _guild

KEYS=("Study_id","Network_id","Year","Latitude","Longitude","island_id")
GUILDS=("bombus","all_bees")
MORPH=("flower_primary_color","floral_form","floral_symmetry")


def _q_joint(morph: pd.DataFrame, *, n: int, seed: int, draws: int=50) -> float:
    if len(morph)<n:
        raise ValueError("Insufficient directly known same-network flora")
    # No projected 17-D KDE or single-profile inflation.
    codes=morph.astype(str).to_numpy()
    rng=np.random.default_rng(seed)
    vals=[]
    for _ in range(draws):
        arr=codes[rng.choice(len(codes),size=n,replace=False)]
        counts=np.asarray(list(Counter(map(tuple,arr)).values()),dtype=float)
        p=counts/n
        vals.append(1-float((p*p).sum()))
    return float(np.mean(vals))


def _group_remove(g: pd.DataFrame, *, group: str, direct: pd.DataFrame,
                  seed: int=20261113) -> dict | None:
    """Quantify *documented* plant coverage under hypothetical guild removal."""
    if group not in GUILDS:
        raise ValueError("Unknown frozen guild stress scenario")
    if g.Plant_rank.ne("SPECIES").all():
        return None
    x=g.loc[g.Plant_rank.eq("SPECIES")].copy()
    if len(x)==0:
        return None
    x["visitor_guild"]=x.apply(_guild,axis=1)
    orders=set(x.loc[x.Pollinator_order.isin((
        "Hymenoptera","Diptera","Lepidoptera","Coleoptera")),"Pollinator_order"])
    if len(orders)<2:
        return None
    all_plants=set(x.Plant_accepted_name)
    if len(all_plants)<10:
        return None
    excluded=(x.visitor_guild.eq("bombus") if group=="bombus"
              else x.visitor_guild.isin(("bombus","other_bees")))
    lost_guild_plants=set(x.loc[excluded,"Plant_accepted_name"])
    if len(lost_guild_plants)==0:
        return None
    retained=set(x.loc[~excluded,"Plant_accepted_name"])
    at_risk=len(all_plants-retained)
    overlap=len(lost_guild_plants&retained)
    # Original sampled insect taxon identities, unlike record rows, allow a
    # qualified random VISITOR-SPECIES deletion null matching the number
    # of distinct targeted pollinator species in THIS network. This is not
    # an ecological randomization experiment and lacks visit-effort controls.
    target_visitor_names=set(x.loc[excluded,"Pollinator_accepted_name"])
    visitor_plants={
        str(k):set(z.Plant_accepted_name)
        for k,z in x.groupby("Pollinator_accepted_name")
    }
    all_visitors=sorted(visitor_plants)
    rand_risk=None
    if 0<len(target_visitor_names)<len(all_visitors) and len(all_visitors)>=3:
        fingerprint="|".join(map(str,sorted(set(zip(
            x.Study_id,x.Network_id,x.Year,x.island_id)))))
        digest=hashlib.sha256((
             "FD-random-visitor:"+str(seed)+fingerprint+group).encode()).digest()
        rng=np.random.default_rng(int.from_bytes(digest[:8],"big"))
        universe=np.asarray(all_visitors,dtype=object)
        risks=[]
        for _ in range(199):
            removed=set(rng.choice(universe,size=len(target_visitor_names),replace=False))
            surviving=set().union(*(species for visitor,species in visitor_plants.items()
                                    if visitor not in removed))
            risks.append(len(all_plants-surviving)/len(all_plants))
        rand_risk=float(np.mean(risks))
    # Shares are nonexclusive, because a flower can have multiple
    # alternative-guild documented visitor types.
    alternative_shares={}
    for alternate in ("other_bees","diptera","lepidoptera_pooled","coleoptera"):
        alternate_plants=set(x.loc[
           (~excluded)&x.visitor_guild.eq(alternate),
           "Plant_accepted_name"])
        alternative_shares[
          "share_exposed_plants_with_alternative_"+alternate
        ]=len(lost_guild_plants & alternate_plants)/len(lost_guild_plants)

    # Equal 8 Direct plants before and after, conditional on sufficient
    # exact same-island botanical source coverage. Not flower fitness.
    pre=sorted(all_plants&set(direct.index))
    post=sorted(retained&set(direct.index))
    qdelta=None
    if len(pre)>=8 and len(post)>=8:
        digest=hashlib.sha256(("FDnet:"+str(seed)+":"+str(sorted(set(
            zip(x.Study_id,x.Network_id,x.Year,x.island_id))))+group).encode()).digest()
        seed0=int.from_bytes(digest[:8],"big")
        qpre=_q_joint(direct.loc[pre,list(MORPH)],n=8,seed=seed0)
        qpost=_q_joint(direct.loc[post,list(MORPH)],n=8,seed=seed0+1)
        qdelta=float(qpost-qpre)
    return {
        "n_original_visited_plants":len(all_plants),
        "n_documented_guild_exposed_plant_species":len(lost_guild_plants),
        "n_documented_guild_plant_species_with_other_visitor_records":overlap,
        "share_guild_exposed_plant_species_with_other_visitor_records":
           overlap/len(lost_guild_plants),
        "fraction_of_all_documented_visited_plant_species_without_other_links":
           at_risk/len(all_plants),
        "n_same_original_direct_morph_species_before":len(pre),
        "n_same_original_direct_morph_species_after":len(post),
        "equal_8_direct_morph_joint_Q_post_minus_pre":qdelta,
        "n_distinct_documented_visitor_species":len(all_visitors),
        "n_removed_documented_visitor_species":len(target_visitor_names),
        "mean_random_equal_visitor_species_deletion_plant_loss":rand_risk,
        "extra_plant_loss_vs_random_equal_visitor_species_deletion":(
           at_risk/len(all_plants)-rand_risk if rand_risk is not None else None),
        **alternative_shares,
    }


def original_island_network_guild_ablation(
    visits: pd.DataFrame,
    source_site_islands: pd.DataFrame,
    original_flora: pd.DataFrame,
    original_direct: pd.DataFrame,
    island_geo: pd.DataFrame,
    metadata: pd.DataFrame,
) -> dict:
    """Aggregate only. All private network and plant data stay in process."""
    original=validate_original_visits(visits)
    site_keys=["Study_id","Network_id","Latitude","Longitude"]
    if not set(site_keys+["island_id"]).issubset(source_site_islands):
        raise ValueError("Original source network geography not available")
    if not {"island_id","accepted_species","origin_status"}.issubset(original_flora):
        raise ValueError("Original botanical flora join contract missing")
    if not {"island_id","analysis_regime","spatial_block"}.issubset(island_geo):
        raise ValueError("Original block geography missing")
    if not set(("Study_id","Year","Taxa_recorded")).issubset(metadata):
        raise ValueError("Original EuPPollNet study sampling scope metadata missing")
    source_site=source_site_islands[site_keys+["island_id"]].drop_duplicates()
    if source_site.duplicated(site_keys).any():
        raise ValueError("Original site has multiple islands")
    valid=set(original_flora.island_id)
    source_site=source_site.loc[source_site.island_id.isin(valid)]
    local=original.merge(source_site,on=site_keys,how="inner",validate="many_to_one")
    if local.empty:
        return {"status":"NO_OBSERVED_DIRECT_ISLAND_VISITOR_NETWORKS",
                "historical_pollinator_decline_measured":False}
    local=local.merge(island_geo[["island_id","analysis_regime","spatial_block"]],
                      on="island_id",how="left",validate="many_to_one")
    if local.spatial_block.isna().any():
        raise ValueError("Original spatial blocks missing")
    # Metadata describes study sampling scope; observed presence alone is NOT
    # grounds for calling a missing insect order absent.
    meta=metadata[["Study_id","Year","Taxa_recorded"]].copy()
    meta.Year=pd.to_numeric(meta.Year,errors="coerce")
    meta=meta.dropna(subset=["Year"])
    meta["Year"]=meta.Year.astype(int)
    if meta.duplicated(["Study_id","Year"]).any():
        # Study years with conflicting sampling descriptions cannot support
        # a within-study network guild deletion comparison.
        meta=meta.loc[~meta.duplicated(["Study_id","Year"],keep=False)]
    meta["declared_broad_taxa"]=meta.Taxa_recorded.fillna("").astype(str).str.lower().map(
        lambda t: ("bees" in t or "hymenoptera" in t or "bumblebees" in t) and
                  ("diptera" in t or "syrphid" in t or "hover" in t or "flies" in t))
    local=local.merge(meta[["Study_id","Year","declared_broad_taxa"]],
                      on=["Study_id","Year"],how="left",validate="many_to_one")
    local["declared_broad_taxa"]=local.declared_broad_taxa.fillna(False)
    # Strictly require independent study-level evidence that bee AND an
    # alternative pollinator order were within the reported taxa scope.
    broadly=local.loc[local.declared_broad_taxa].copy()
    n_sites=local[list(KEYS)].drop_duplicates().shape[0]
    n_studies=local.Study_id.nunique()
    n_broad_studies=broadly.Study_id.nunique()
    rows=[]
    for key,g in broadly.groupby(list(KEYS),sort=False):
        group_geo={c: g[c].iloc[0] for c in ("analysis_regime","spatial_block")}
        for scenario in GUILDS:
            d=_group_remove(g,group=scenario,direct=original_direct)
            if d is None:
                continue
            rows.append({
                "scenario":scenario,"analysis_regime":group_geo["analysis_regime"],
                "spatial_block":group_geo["spatial_block"],
                "island_id":key[5],"study_id":key[0],
                **d
            })
    sub=pd.DataFrame(rows)
    results=[]
    if not sub.empty:
        for scenario,group in sub.groupby("scenario"):
            values=group["share_guild_exposed_plant_species_with_other_visitor_records"]
            risk=group["fraction_of_all_documented_visited_plant_species_without_other_links"]
            q=group["equal_8_direct_morph_joint_Q_post_minus_pre"].dropna()
            excess=group["extra_plant_loss_vs_random_equal_visitor_species_deletion"].dropna()
            eq_island=group.groupby("island_id")[
              ["share_guild_exposed_plant_species_with_other_visitor_records",
               "fraction_of_all_documented_visited_plant_species_without_other_links"]
            ].mean()
            eq_study=group.groupby("study_id")[
              ["share_guild_exposed_plant_species_with_other_visitor_records",
               "fraction_of_all_documented_visited_plant_species_without_other_links"]
            ].mean()
            results.append({
                "scenario":scenario,
                "n_within_study_site_year_networks":int(len(group)),
                "n_original_islands":int(group.island_id.nunique()),
                "n_original_studies":int(group.study_id.nunique()),
                "n_spatial_blocks":int(group.spatial_block.nunique()),
                "median_alternative_visitor_link_coverage_among_exposed_plants":float(values.median()),
                "median_fraction_all_plants_losing_all_documented_other_visitors":float(risk.median()),
                "island_equal_mean_fraction_plant_links_lost":float(eq_island[
                  "fraction_of_all_documented_visited_plant_species_without_other_links"].mean()),
                "study_equal_mean_fraction_plant_links_lost":float(eq_study[
                  "fraction_of_all_documented_visited_plant_species_without_other_links"].mean()),
                "island_equal_mean_share_exposed_with_other_visitor_links":float(eq_island[
                  "share_guild_exposed_plant_species_with_other_visitor_records"].mean()),
                "study_equal_mean_share_exposed_with_other_visitor_links":float(eq_study[
                  "share_guild_exposed_plant_species_with_other_visitor_records"].mean()),
                "n_networks_with_equal_removed_visitor_species_null":int(len(excess)),
                "median_extra_plant_loss_vs_random_equal_visitor_species":(
                   float(excess.median()) if len(excess)>=10 else None),
                "study_equal_mean_extra_loss_vs_random_visitor_species":(
                   float(group.groupby("study_id")[
                     "extra_plant_loss_vs_random_equal_visitor_species_deletion"
                   ].mean().mean()) if len(excess)>=10 else None),
                "median_study_network_alternative_guild_coverage_by_kind":{
                   guild:float(group["share_exposed_plants_with_alternative_"+guild].median())
                   for guild in ("other_bees","diptera","lepidoptera_pooled","coleoptera")
                },
                "n_networks_with_8_matched_Direct_morph_plants_after_loss":int(len(q)),
                "median_equal_8_joint_morph_Q_change_if_supported":float(q.median()) if len(q)>=5 else None,
                "functional_diversity_claim_supported":bool(
                   len(q)>=30 and group.island_id.nunique()>=10 and
                   group.spatial_block.nunique()>=5 and group.study_id.nunique()>=3),
                "unobserved_guilds_never_called_absent":True,
                "simulated_loss_not_observed_historical_decline":True,
                "other_visitor_links_not_effective_pollination":True,
            })
    return {
        "status":"EXPLORATORY_OBSERVED_VISITOR_LINK_REMOVAL_NOT_CAUSAL_EFFECT",
        "n_all_same_island_observed_source_network_site_years":int(n_sites),
        "n_original_network_studies":int(n_studies),
        "n_studies_with_source_documented_bee_and_nonbee_sampling_scope":int(n_broad_studies),
        "n_qualified_network_simulation_scenarios":int(len(sub)),
        "scenario_aggregates":results,
        "all_networks_cannot_be_treated_as_independent_island_histories":True,
        "observed_pollinator_loss_shift_or_effective_pollen_delivery":False,
        "raw_taxon_network_site_or_source_ids_exported":False,
    }
