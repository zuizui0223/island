"""Whole-botanical-family-heldout EXACT SPECIES GloBI cited visitor composition.

This is NOT island-local pollination and NOT evolutionary switching.
Only source-v0.9 original bibliography identities define evidence. Compare
an IDENTICAL species/family sample before/after removing the archive-FIRST
label that dominates the outcome-independent 928-island Direct species panel.
Project five nonnegative visitor-reference shares to the simplex.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import numpy as np
import pandas as pd

GUILDS=("bombus","non_bombus_bees","lepidoptera",
        "flower_visiting_birds","diptera")
TRAITS=("flower_primary_color","floral_form","floral_symmetry")


def _simplex(q):
    u=np.sort(q,axis=1)[:,::-1]
    cs=np.cumsum(u,axis=1)-1.
    k=np.arange(1,q.shape[1]+1)
    rho=np.sum(u-cs/k>0,axis=1)
    return np.maximum(q-cs[np.arange(len(q)),rho-1,None]/rho[:,None],0.)


def _ridge_predict_families(X,Y,groups,alpha):
    pred=np.empty_like(Y)
    for family in np.unique(groups):
        test=np.flatnonzero(groups==family)
        train=np.flatnonzero(groups!=family)
        means=X[train].mean(axis=0)
        s=X[train].std(axis=0)
        s[s<1e-12]=1
        a=(X[train]-means)/s
        b=(X[test]-means)/s
        intercept=Y[train].mean(axis=0)
        beta=np.linalg.solve(
            a.T@a+alpha*np.eye(a.shape[1]),a.T@(Y[train]-intercept))
        pred[test]=intercept+b@beta
    return _simplex(pred)


def species_family_documentation_syndrome_transfer(
    original_documented_units:set[tuple],
    providers_per_bibliography:dict,
    source_direct_ledger:pd.DataFrame,
    source_gbif_taxonomy:pd.DataFrame,
    source_flora_928_species:set[str], *,
    min_distinct_references:int=3,
    min_family_groups:int=15,
    min_species:int=50,
    alpha:float=100.,
    bootstrap:int=2999,
    seed:int=20261112,
)->dict:
    """Return only released aggregates, original plant/reference IDs private."""
    if (min_distinct_references<3 or min_family_groups<15 or
        min_species<50 or bootstrap<199 or alpha<=0):
        raise ValueError("Unsafe original source/citation/family evidence threshold")
    expected={"accepted_species","trait_name","normalized_value",
              "evidence_scope","evidence_quality"}
    if not expected.issubset(source_direct_ledger):
        raise ValueError("Original Direct species floral source missing")
    if not {"genus","family"}.issubset(source_gbif_taxonomy):
        raise ValueError("Source GBIF family taxonomy missing")
    x=source_direct_ledger.loc[
        source_direct_ledger.trait_name.isin(TRAITS) &
        source_direct_ledger.evidence_scope.eq("species_direct") &
        source_direct_ledger.evidence_quality.isin(("high","medium")) &
        source_direct_ledger.normalized_value.notna() &
        ~source_direct_ledger.normalized_value.astype(str).str.contains("|",regex=False)]
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original Direct species floral evidence conflict")
    x=x.pivot(index="accepted_species",columns="trait_name",
              values="normalized_value").dropna(subset=list(TRAITS))
    x=x.loc[x.index.intersection(source_flora_928_species)].copy()
    if x.empty or x.index.has_duplicates:
        raise ValueError("No source Direct fixed-flora species")
    x["genus"]=x.index.astype(str).str.split().str[0]
    taxa=source_gbif_taxonomy[["genus","family"]].dropna().copy()
    taxonomy=taxa.groupby("genus").family.agg(lambda t:tuple(sorted(set(t))))
    taxonomy=taxonomy[taxonomy.map(len).eq(1)].map(lambda t:t[0])
    x["family"]=x.genus.map(taxonomy)
    x=x.dropna(subset=["family"]).sort_index()
    if x.empty:
        raise ValueError("No uniquely resolved source botanical family")
    duplicate={ref for ref,providers in providers_per_bibliography.items()
               if len(providers)>1}
    units={
        (plant,guild,ref,provider)
        for plant,guild,ref,provider,evidence in original_documented_units
        if plant in x.index and ref not in duplicate and guild in GUILDS
        and provider and ref
    }
    if not units:
        return {"status":"NO_CITED_SOURCE_SPECIES", "source_location_visitor_effectiveness_tested":False}
    source_count=Counter(p for _,_,_,p in units)
    removed=sorted(source_count,key=lambda p:(-source_count[p],p))[0]

    def make_frame(v):
        rec=defaultdict(lambda:defaultdict(set))
        archives=defaultdict(set)
        for plant,guild,ref,provider in v:
            rec[plant][guild].add(ref)
            archives[plant].add(provider)
        frame=pd.DataFrame(index=x.index)
        for guild in GUILDS:
            frame[guild]=[len(rec[p][guild]) for p in x.index]
        frame["n_cited_ref_keys"]=[
            len(set().union(*rec[p].values())) if rec[p] else 0
            for p in x.index]
        frame["n_original_archive_labels"]=[len(archives[p]) for p in x.index]
        return frame

    pre=make_frame(units)
    post=make_frame({u for u in units if u[3]!=removed})
    elig=(pre.n_cited_ref_keys.ge(min_distinct_references) &
          post.n_cited_ref_keys.ge(min_distinct_references))
    eligible=x.loc[elig]
    frame_pre=pre.loc[elig]
    frame_post=post.loc[elig]
    fam=eligible.family.to_numpy(str)
    n_family=len(np.unique(fam))
    common=dict(
        n_original_Direct_triple_complete_source_flora_species=int(len(x)),
        n_cited_exact_plant_species_before_omission=int(pre.n_cited_ref_keys.ge(1).sum()),
        n_cited_exact_plant_species_after_omission=int(post.n_cited_ref_keys.ge(1).sum()),
        n_strict_common_species_with_3plus_documentation=int(len(eligible)),
        n_strict_common_botanical_families=int(n_family),
        n_archive_first_source_labels=int(len(source_count)),
        largest_archive_first_label_fraction_of_source_units=float(
            source_count[removed]/len(units)),
        n_crossarchive_bibliographic_key_collisions_removed=int(len(duplicate)),
        same_species_same_family_before_after_sample=True,
        original_archive_first_label_not_independent_study=True,
        source_location_visitor_effectiveness_tested=False,
        evolutionary_pollinator_shift_confirmed=False,
    )
    if len(eligible)<min_species or n_family<min_family_groups:
        return {"status":"NOT_IDENTIFIABLE_AFTER_ORIGINAL_SOURCE_OMISSION",**common}

    form=eligible.floral_form.eq("open_radial").to_numpy(float)
    sym=eligible.floral_symmetry.eq("actinomorphic").to_numpy(float)
    col=eligible.flower_primary_color.eq("green_brown_inconspicuous").to_numpy(float)
    traits=np.column_stack((form,sym,col))
    forms={
        "source_documentation_only":0,
        "plus_open_flower":1,
        "plus_open_and_symmetry":2,
        "plus_three_independent_flower_states":3,
        "plus_three_states_and_pairwise_coordination":6,
    }
    interaction=np.column_stack((form*sym,form*col,sym*col))
    rng=np.random.default_rng(seed)
    pooled={}
    for stage,frame in (("before",frame_pre),("after",frame_post)):
        recorded=frame[list(GUILDS)].to_numpy(float)
        totals=recorded.sum(axis=1)
        if (totals<min_distinct_references).any():
            raise ValueError("Invalid matched original bibliography shares")
        y=recorded/totals[:,None]
        effort=np.column_stack((
            np.log1p(frame.n_cited_ref_keys.to_numpy(float)),
            np.log1p(frame.n_original_archive_labels.to_numpy(float))))
        model_features={
            "source_documentation_only":effort,
            "plus_open_flower":np.column_stack((effort,traits[:,:1])),
            "plus_open_and_symmetry":np.column_stack((effort,traits[:,:2])),
            "plus_three_independent_flower_states":np.column_stack((effort,traits)),
            "plus_three_states_and_pairwise_coordination":
                np.column_stack((effort,traits,interaction)),
        }
        errs={}
        for name,features in model_features.items():
            pred=_ridge_predict_families(features,y,fam,alpha)
            errs[name]=((pred-y)**2).mean(axis=1)
        group_ix=[np.flatnonzero(fam==g) for g in np.unique(fam)]
        baseline=errs["source_documentation_only"]
        results={}
        for name in model_features:
            gain=baseline-errs[name]
            fm=np.array([gain[ids].mean() for ids in group_ix])
            draws=fm[rng.integers(n_family,size=(bootstrap,n_family))].mean(axis=1)
            results[name]=dict(
                n_added_flower_features=forms[name],
                genus_species_equal_MSE=float(errs[name].mean()),
                species_equal_gain_over_source_effort=float(gain.mean()),
                family_equal_gain_over_source_effort=float(fm.mean()),
                family_equal_bootstrap_025=float(np.quantile(draws,.025)),
                family_equal_bootstrap_975=float(np.quantile(draws,.975)),
                n_families_positive=int((fm>0).sum()))
        pooled[stage]=results
    return {
        "status":"COMPLETED_EXPLORATORY_EXACT_SPECIES_FAMILY_LOFO_DOCUMENTARY_SHARES",
        **common, "models_before_source_label_omission":pooled["before"],
        "models_after_source_label_omission":pooled["after"],
        "floral_coordination_is_predictive_not_evolutionary_adaptation":True,
        "original_plant_family_archive_or_citation_identifiers_released":False,
    }
