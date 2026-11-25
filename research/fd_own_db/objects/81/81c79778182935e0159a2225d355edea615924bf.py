"""Independent reproductive-assurance × joint-flower-profile support/assembly audit.

Original species Direct high/medium only. No low-confidence genus imputation,
no inference that SC is autonomous selfing or that reproductive states confer
effective pollen service. Island and mainland source records remain private.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import numpy as np
import pandas as pd

from .pollinator_functional_diversity import direct_triple

REPRO = {
 "self_compatibility":("self_incompatibility",{"SC"},{"SI"}),
 "autonomous_selfing":("autonomous_selfing_capacity",{"autonomous","delayed"},{"absent"}),
 "selfing_mating_system":("mating_system",{"predominantly_selfing","obligate_selfing"},{"predominantly_outcrossing"}),
}


def direct_binary_assurance(ledger:pd.DataFrame) -> pd.DataFrame:
    req={"accepted_species","trait_name","normalized_value",
         "evidence_scope","evidence_quality"}
    if not req.issubset(ledger):
        raise ValueError("Direct source evidence columns missing")
    x=ledger.loc[
        ledger.trait_name.isin([v[0] for v in REPRO.values()]) &
        ledger.evidence_scope.eq("species_direct") &
        ledger.evidence_quality.isin(("high","medium")) &
        ledger.normalized_value.notna() &
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False),
        list(req)].copy()
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Conflicting Direct reproductive trait")
    scored={}
    for trait,(original,pos,neg) in REPRO.items():
        sub=x.loc[x.trait_name.eq(original),["accepted_species","normalized_value"]].copy()
        sub=sub.loc[sub.normalized_value.isin(pos|neg)]
        if sub.accepted_species.duplicated().any():
            raise ValueError("Duplicate Direct trait after recoding")
        scored[trait]=sub.assign(state=sub.normalized_value.isin(pos).astype(int)).set_index(
            "accepted_species").state
    return pd.DataFrame(scored).astype(float)


def _joint(arr):
    n=len(arr)
    counts=np.asarray(list(Counter(map(tuple,arr)).values()),dtype=float)
    return float(1-np.square(counts/n).sum())


def _seed(name,trait,seed):
    return int.from_bytes(hashlib.sha256(f"{seed}:{name}:{trait}".encode()).digest()[:8],"big")


def _sampling(binary:np.ndarray,morph:np.ndarray,*,n:int,repeats:int,rng):
    scores=[]
    for _ in range(repeats):
        ids=rng.choice(len(binary),size=n,replace=False)
        p=float(binary[ids].mean())
        scores.append((p,2*p*(1-p),_joint(morph[ids])))
    return np.mean(scores,axis=0)


def _block_ci(frame,key,draws,seed):
    b=frame.groupby("spatial_block")[key].mean().dropna().to_numpy(float)
    if len(b)<5:return (None,None,None)
    rng=np.random.default_rng(seed)
    boot=b[rng.integers(len(b),size=(draws,len(b)))].mean(axis=1)
    return tuple(map(float,(b.mean(),np.quantile(boot,.025),np.quantile(boot,.975))))


def audit_assurance_flora(
    original_ledger:pd.DataFrame, original_flora:pd.DataFrame,
    geo:pd.DataFrame, sources:pd.DataFrame, native_mainland:pd.DataFrame,
    guilds:pd.DataFrame, *,
    seed:int=20261113, sample_size:int=20, repeats:int=20,
    draws:int=1999,
)->dict:
    if sample_size!=20 or repeats<10 or draws<999:
        raise ValueError("Reproductive trait independent support frozen contract")
    morph=direct_triple(original_ledger)
    binary=direct_binary_assurance(original_ledger)
    required={
      "flora":(original_flora,{"island_id","accepted_species","origin_status"}),
      "geography":(geo,{"island_id","analysis_regime","spatial_block"}),
      "sources":(sources,{"island_id","source_mode","source_rank","entity_ID"}),
      "GIFT":(native_mainland,{"entity_ID","work_species"}),
      "guild":(guilds,{"island_id","channel_id","channel_state","source_state"}),
    }
    for k,(t,col) in required.items():
        if not col.issubset(t):raise ValueError("Original "+k+" contract missing")
    original=original_flora.loc[original_flora.accepted_species.isin(morph.index)].copy()
    if original.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate island plant")
    if not original.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Unrecognized original origin status")
    counts=original.groupby("island_id").accepted_species.nunique()
    elig=set(counts[counts>=30].index)
    if len(elig)!=928:
        raise ValueError("Original 928 island Direct morphology universe changed")
    original=original.loc[original.island_id.isin(elig)]
    g=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates().set_index("island_id")
    if g.index.has_duplicates or not elig.issubset(g.index):
        raise ValueError("Original geography not available for all 928")
    donors=sources.loc[
       sources.island_id.isin(elig)&sources.source_mode.eq("geo_k5")&
       sources.source_rank.between(1,5),["island_id","source_rank","entity_ID"]]
    if donors.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicated fixed GIFT source")
    src=donors.sort_values(["island_id","source_rank"]).groupby("island_id").entity_ID.agg(tuple).to_dict()
    if any(len(src.get(i,()))!=5 or len(set(src[i]))!=5 for i in elig):
        raise ValueError("Missing five GIFT geography source entities")
    gift=(native_mainland.loc[native_mainland.work_species.isin(morph.index),
                              ["entity_ID","work_species"]].dropna().drop_duplicates()
          .groupby("entity_ID").work_species.agg(set).to_dict())
    ch=guilds[guilds.island_id.isin(elig)][["island_id","channel_id","channel_state","source_state"]].copy()
    if ch.duplicated(["island_id","channel_id"]).any():
        raise ValueError("Original channel IDs duplicated")
    if not (ch.loc[ch.channel_state.eq("disrupted"),"source_state"].eq("available")).all():
        raise ValueError("Qualified nondetection lacks source support")
    c=ch.pivot(index="island_id",columns="channel_id",values="channel_state")
    if not {"bombus","flower_visiting_birds"}.issubset(c.columns):
        raise ValueError("No original Bombus/bird opportunity")
    focal=set(c.index[c.bombus.eq("disrupted")&c.flower_visiting_birds.eq("retained")])
    if len(focal)!=21:
        raise ValueError("21 original focal candidate opportunity islands drifted")
    byisland={i:x for i,x in original.groupby("island_id")}
    rows=[]
    for island in sorted(elig):
        gi=byisland[island]
        src_names=set().union(*(gift.get(entity,set()) for entity in src[island]))
        src_names=src_names.intersection(morph.index)
        island_species=set(gi.accepted_species)
        for scope in ("all_recorded","source_native"):
            local=gi if scope=="all_recorded" else gi.loc[gi.origin_status.eq("native")]
            local_names=set(local.accepted_species)
            for trait in REPRO:
                values=binary[trait]
                known=set(values.dropna().index)
                isl=sorted(local_names & known)
                source=sorted(src_names & known)
                # Prevent different trait carrier count from changing effective
                # support under the metric itself. Never impute unknown states.
                if len(isl)<sample_size:
                    continue
                rand=np.random.default_rng(_seed(island,scope+trait,seed))
                local_morph=morph.loc[isl].to_numpy(str)
                local_bin=values.loc[isl].to_numpy(float)
                loc=_sampling(local_bin,local_morph,n=sample_size,repeats=repeats,rng=rand)
                src_val=None
                if len(source)>=sample_size:
                    rng=np.random.default_rng(_seed(island,"source"+trait,seed))
                    src_val=_sampling(values.loc[source].to_numpy(float),
                            morph.loc[source].to_numpy(str),
                            n=sample_size,repeats=repeats,rng=rng)
                d=dict(island_id=island,spatial_block=g.loc[island,"spatial_block"],
                    regime=g.loc[island,"analysis_regime"],source_signature=src[island][0],
                    trait=trait,scope=scope,n_island_known=len(isl),
                    n_source_known=len(source),candidate_bombus_bird=island in focal,
                    island_positive_share=float(loc[0]),island_binary_Q=float(loc[1]),
                    island_morph_joint_Q=float(loc[2]))
                for k,v in enumerate(("positive_share","binary_Q","morph_joint_Q")):
                    d["delta_"+v]=float(loc[k]-src_val[k]) if src_val is not None else np.nan
                rows.append(d)
    private=pd.DataFrame.from_records(rows)
    out=[]
    for (scope,trait,regime),q in private.groupby(["scope","trait","regime"]):
        if len(q)<10 or q.spatial_block.nunique()<5: continue
        pair=q.dropna(subset=["delta_positive_share"])
        rec=dict(scope=scope,trait=trait,regime=regime,n_islands=len(q),
            n_source_comparable=len(pair),
            n_spatial_blocks=q.spatial_block.nunique(),
            median_island_positive_share=float(q.island_positive_share.median()),
            median_island_reproductive_binary_Q=float(q.island_binary_Q.median()),
            median_island_morph_joint_Q_same_reproductive_carriers=float(q.island_morph_joint_Q.median()))
        if len(pair)>=10 and pair.spatial_block.nunique()>=5:
            for j,key in enumerate(("delta_positive_share","delta_binary_Q","delta_morph_joint_Q")):
                mean,lo,hi=_block_ci(pair,key,draws,seed+j)
                rec["block_mean_"+key]=mean
                rec["block_ci025_"+key]=lo
                rec["block_ci975_"+key]=hi
        # Post-hoc source-signature sensitivity (distinct from spatial
        # blocks, not independent evolutionary histories).
        rec["n_nearest_mainland_source_signatures"]=pair.source_signature.nunique()
        if len(pair)>=10 and pair.source_signature.nunique()>=5:
            for j,key in enumerate(("delta_positive_share","delta_binary_Q","delta_morph_joint_Q")):
                means=pair.groupby("source_signature")[key].mean().dropna().to_numpy(float)
                rng=np.random.default_rng(seed+200+j)
                sims=means[rng.integers(len(means),size=(draws,len(means)))].mean(axis=1)
                rec["source_equal_mean_"+key]=float(means.mean())
                rec["source_ci025_"+key]=float(np.quantile(sims,.025))
                rec["source_ci975_"+key]=float(np.quantile(sims,.975))
        out.append(rec)
    cand=[]
    for (scope,trait),q in private.loc[private.candidate_bombus_bird].groupby(["scope","trait"]):
        pair=q.dropna(subset=["delta_positive_share"])
        rec=dict(scope=scope,trait=trait,n_candidate_islands_with_20_known=len(q),
                 n_source_comparable=len(pair),n_spatial_blocks=q.spatial_block.nunique())
        if len(pair)>=10 and pair.spatial_block.nunique()>=5:
            for j,key in enumerate(("delta_positive_share","delta_morph_joint_Q")):
                rec["block_mean_"+key],rec["block_ci025_"+key],rec["block_ci975_"+key]=(
                    _block_ci(pair,key,draws,seed+100+j))
        cand.append(rec)
    # In exactly the SAME selected native-eligible islands, compare scope to
    # avoid attributing different geography/coverage to island origin status.
    paired=[]
    for trait,sub in private.groupby("trait"):
        all_=sub.loc[sub.scope.eq("all_recorded")].set_index("island_id")
        native=sub.loc[sub.scope.eq("source_native")].set_index("island_id")
        for regime,n in native.groupby("regime"):
            common=n.index.intersection(all_.index)
            if len(common)<10 or n.spatial_block.nunique()<5:
                continue
            a=all_.loc[common]
            n=n.loc[common]
            row=dict(trait=trait,regime=regime,n_matched_native_eligible_islands=len(common),
                n_matched_spatial_blocks=n.spatial_block.nunique(),
                compare_exact_same_island_not_independent_replication=True)
            vals=n[["spatial_block"]].copy()
            for j,key in enumerate(("island_positive_share","island_binary_Q",
                                    "island_morph_joint_Q","delta_positive_share")):
                vals["native_minus_all_"+key]=(n[key]-a[key]).to_numpy(float)
                mean,lo,hi=_block_ci(vals,"native_minus_all_"+key,draws,seed+500+j)
                row["block_mean_native_minus_all_"+key]=mean
                row["block_ci025_native_minus_all_"+key]=lo
                row["block_ci975_native_minus_all_"+key]=hi
                if key=="delta_positive_share":
                    row["median_same_island_all_source_delta"]=(
                        float(a[key].median()) if a[key].notna().any() else None)
                    row["median_same_island_native_source_delta"]=(
                        float(n[key].median()) if n[key].notna().any() else None)
            paired.append(row)
    # Post-hoc safeguard against comparing SC and mating system on different
    # trait-known species: both binary states, morphology and island SOURCE
    # comparison now use the EXACT SAME 20 accepted plant species per draw.
    shared=binary[["self_compatibility","selfing_mating_system"]].dropna()
    shared_names=set(shared.index)&set(morph.index)
    matched_rows=[]
    for island in sorted(elig):
        group=byisland[island]
        source_names=set().union(*(gift.get(entity,set()) for entity in src[island]))
        source_names=source_names&shared_names
        if len(source_names)<sample_size:
            continue
        local_all=set(group.accepted_species)
        for scope in ("all_recorded","source_native"):
            local_names=(
                local_all if scope=="all_recorded"
                else set(group.loc[group.origin_status.eq("native"),"accepted_species"])
            ) & shared_names
            if len(local_names)<sample_size:
                continue
            score={}
            for label,names_ in (("island",sorted(local_names)),
                                 ("source",sorted(source_names))):
                state=shared.loc[names_].to_numpy(float)
                rng=np.random.default_rng(_seed(island,scope+label+"joint_repro",seed))
                sampled=[]
                for _ in range(repeats):
                    ids=rng.choice(len(names_),size=sample_size,replace=False)
                    sampled.append(state[ids].mean(axis=0))
                score[label]=np.mean(sampled,axis=0)
            row=dict(island_id=island,scope=scope,
               regime=g.loc[island,"analysis_regime"],
               spatial_block=g.loc[island,"spatial_block"],
               source_signature=src[island][0],
               n_same_carrier_island=len(local_names),
               n_same_carrier_mainland=len(source_names),
               paired_delta_SC=float(score["island"][0]-score["source"][0]),
               paired_delta_selfing=float(score["island"][1]-score["source"][1]))
            row["paired_delta_selfing_minus_SC"]=(
                row["paired_delta_selfing"]-row["paired_delta_SC"])
            matched_rows.append(row)
    matched_private=pd.DataFrame(matched_rows)
    matched_summary=[]
    if not matched_private.empty:
        for (scope,regime),part in matched_private.groupby(["scope","regime"]):
            if len(part)<10 or part.spatial_block.nunique()<5:
                continue
            entry=dict(scope=scope,regime=regime,
                n_islands_with_20_same_species_each_population=len(part),
                n_spatial_blocks=part.spatial_block.nunique(),
                n_nearest_mainland_sources=part.source_signature.nunique(),
                same_carrier_both_traits_and_same_source=True)
            for j,key in enumerate(("paired_delta_SC","paired_delta_selfing",
                                    "paired_delta_selfing_minus_SC")):
                mean,lo,hi=_block_ci(part,key,draws,seed+700+j)
                entry["block_mean_"+key]=mean
                entry["block_ci025_"+key]=lo
                entry["block_ci975_"+key]=hi
                if part.source_signature.nunique()>=5:
                    source_means=part.groupby("source_signature")[key].mean().to_numpy(float)
                    rng=np.random.default_rng(seed+900+j)
                    boot=source_means[rng.integers(len(source_means),
                                    size=(draws,len(source_means)))].mean(axis=1)
                    entry["source_mean_"+key]=float(source_means.mean())
                    entry["source_ci025_"+key]=float(np.quantile(boot,.025))
                    entry["source_ci975_"+key]=float(np.quantile(boot,.975))
            matched_summary.append(entry)
    return {
     "status":"DIRECT_REPRODUCTIVE_FUNCTIONAL_PROXIES_NOT_POLLINATION_SERVICE",
     "n_original_direct_triple_flora_islands":928,
     "n_original_focal_Bombus_nondetection_bird_detection_candidates":21,
     "n_species_with_direct_assurance_by_trait":{
        t:int(binary[t].notna().sum()) for t in REPRO},
     "reproductive_assurance_functional_assemblage_by_region":out,
     "candidate_potential_pollinator_shift_assurance_support":cand,
     "same_island_native_vs_all_assurance_sensitivity":paired,
     "same_plant_carrier_SC_and_mating_system_source_control":matched_summary,
     "carrier_specific_support_min":sample_size,
     "rarefactions":repeats,
     "same_carrier_morphology_and_assurance_scoring":True,
     "unknown_reproductive_trait_is_not_state_zero":True,
     "source_geography_is_not_an_ancestor":True,
     "no_actual_pollinator_loss_or_shift_measured":True,
     "no_species_island_or_source_level_data_exported":True,
    }
