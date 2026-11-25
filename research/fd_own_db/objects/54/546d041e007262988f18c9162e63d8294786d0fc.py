"""Source-disjoint genus prediction for strict binary island floral/reproductive traits.

One heldout source-component at a time. For the strict all_g0 mode, exclude
ALL original evaluated island species names from ALL mainland GIFT donors and
exclude the five-source mainland regions used by the focal component. Tests
use names as units, not repeated island occurrences. Release only aggregated
metrics; never original taxa, genera, islands or mainland source IDs.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from .source_signature_dependency import _component_labels

STATES={
 "self_compatibility":("self_incompatibility",{"SC":1,"SI":0}),
 "floral_symmetry":("floral_symmetry",{"actinomorphic":1,"zygomorphic":0}),
}
TAUS=(.05,1.,8.)
BASE=(
 "region","source_component_size_rank","n_islands","ledger","scope",
 "holdout_policy","trait","evaluation_set","n_test_species",
 "n_test_island_species_rows","n_train_source_disjoint_trait_species",
 "n_test_species_with_genus_donors","test_species_genus_donor_coverage",
 "brier_global",
)
COLUMNS=(*BASE,*[f"brier_genus_tau_{v:g}" for v in TAUS],
         *[f"brier_gain_tau_{v:g}" for v in TAUS])


def holdout_genus_transfer(
    regimes:pd.DataFrame, source_assignments:pd.DataFrame,
    island_flora:pd.DataFrame, native_mainland:pd.DataFrame,
    ledger:pd.DataFrame, *, ledger_label:str, scope_label:str,
    min_islands:int=40, min_test_species:int=50,
    holdout_policy:str="all_g0",
) -> pd.DataFrame:
    """Compare global native-mainland prior with out-of-source genus prediction.

    holdout_policy='component' is a WEAKER sensitivity: names evaluated in
    OTHER components can still enter donors. 'all_g0' is the strict main test.
    Component connectedness is defined within broad regions and does NOT
    imply independent species histories across components or regions.
    """
    if min_islands<40 or min_test_species<50:
        raise ValueError("Privacy/support minimum 40 islands, 50 evaluated names")
    if holdout_policy not in {"all_g0","component"}:
        raise ValueError("Unknown name holdout policy")
    required=[
        (regimes,{"island_id","analysis_regime"}),
        (source_assignments,{"island_id","source_mode","source_rank","entity_ID"}),
        (island_flora,{"island_id","accepted_species","origin_status"}),
        (native_mainland,{"entity_ID","work_species"}),
        (ledger,{"accepted_species","trait_name","normalized_value"}),
    ]
    for frame,fields in required:
        if not fields.issubset(frame.columns):
            raise ValueError("Missing input columns: "+str(sorted(fields-set(frame.columns))))
    g=regimes[["island_id","analysis_regime"]].copy()
    if g.isna().any().any() or g.island_id.duplicated().any():
        raise ValueError("Invalid unique island regimes")
    g=g.astype(str).reset_index(drop=True)
    s=source_assignments[["island_id","source_mode","source_rank","entity_ID"]].copy()
    s.source_rank=pd.to_numeric(s.source_rank,errors="coerce")
    if s.source_rank.isna().any():
        raise ValueError("Invalid source ranks")
    s=s[s.source_mode.eq("geo_k5") & s.source_rank.between(1,5)].copy()
    if (s.duplicated(["island_id","source_rank"]).any()
        or s.duplicated(["island_id","entity_ID"]).any()):
        raise ValueError("Duplicate GIFT donor region or rank")
    signatures=s.groupby("island_id").entity_ID.agg(lambda v:tuple(sorted(set(v))))
    if not signatures.map(len).eq(5).all():
        raise ValueError("Incomplete five-source signature")
    g=g.join(signatures.rename("signature"),on="island_id")
    if g.signature.isna().any():
        raise ValueError("Missing source labels for focal island")
    for _,sub in g.groupby("analysis_regime"):
        siglist=list(sub.signature.unique())
        _,labels=_component_labels(siglist,threshold=1)
        g.loc[sub.index,"component"]=sub.signature.map(dict(zip(siglist,labels)))
    g.component=g.component.astype(int)
    flora=island_flora[["island_id","accepted_species","origin_status"]].copy()
    if flora.isna().any().any():
        raise ValueError("Missing floral name/native status")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate island x taxon record")
    flora=flora[flora.island_id.isin(g.island_id)].merge(
        g[["island_id","analysis_regime","component"]],
        on="island_id",validate="many_to_one")
    all_names=set(flora.accepted_species)
    mainland=native_mainland[["entity_ID","work_species"]].dropna().drop_duplicates()
    trait_data=ledger[["accepted_species","trait_name","normalized_value"]].copy()
    if trait_data.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Duplicate taxon trait evidence")
    states={}
    for trait,(trait_name,allowed) in STATES.items():
        d=trait_data[trait_data.trait_name.eq(trait_name)&
                     trait_data.normalized_value.isin(allowed)]
        states[trait]=pd.Series(
            d.normalized_value.map(allowed).to_numpy(dtype=float),
            index=d.accepted_species,
        )
    out=[]
    for region,sub in g.groupby("analysis_regime",sort=True):
        for rank,(component,n_islands) in enumerate(sub.component.value_counts().items(),start=1):
            if n_islands<min_islands:
                continue
            islands=set(sub.loc[sub.component.eq(component),"island_id"])
            focal=flora[flora.island_id.isin(islands)]
            excluded_names=(all_names if holdout_policy=="all_g0"
                            else set(focal.accepted_species))
            held_regions=set().union(*(set(x) for x in sub.loc[
                sub.component.eq(component),"signature"]))
            donors=mainland.loc[
                ~mainland.entity_ID.isin(held_regions)
                & ~mainland.work_species.isin(excluded_names)]
            donor_names=pd.Index(donors.work_species.unique())
            for trait,known in states.items():
                for eval_set in ("all_recorded","source_native"):
                    test=focal if eval_set=="all_recorded" else focal.loc[
                        focal.origin_status.eq("native")]
                    test=test.loc[test.accepted_species.isin(known.index)]
                    target_names=pd.Index(test.accepted_species.unique())
                    if len(target_names)<min_test_species:
                        continue
                    y=known.reindex(target_names).to_numpy(dtype=float)
                    if not np.isfinite(y).all():
                        raise ValueError("Invalid heldout target state")
                    training_names=donor_names.intersection(known.index)
                    if len(training_names)<min_test_species:
                        continue
                    if not set(training_names).isdisjoint(target_names):
                        raise AssertionError("Test species leaked into GIFT donor training")
                    training_y=known.reindex(training_names).to_numpy(dtype=float)
                    baseline=float(training_y.mean())
                    gdf=pd.DataFrame({
                        "genus":training_names.str.split().str[0],
                        "state":training_y,
                    }).groupby("genus").state.agg(["count","sum"])
                    target_genera=target_names.str.split().str[0]
                    n=pd.Index(target_genera).map(gdf["count"]).fillna(0).to_numpy(float)
                    k=pd.Index(target_genera).map(gdf["sum"]).fillna(0).to_numpy(float)
                    base_brier=float(np.mean((y-baseline)**2))
                    row=dict(
                        region=region,source_component_size_rank=int(rank),
                        n_islands=int(n_islands),ledger=ledger_label,scope=scope_label,
                        holdout_policy=holdout_policy,trait=trait,
                        evaluation_set=eval_set,n_test_species=int(len(target_names)),
                        n_test_island_species_rows=int(len(test)),
                        n_train_source_disjoint_trait_species=int(len(training_names)),
                        n_test_species_with_genus_donors=int((n>0).sum()),
                        test_species_genus_donor_coverage=float((n>0).mean()),
                        brier_global=base_brier,
                    )
                    for tau in TAUS:
                        pred=(k+tau*baseline)/(n+tau)
                        err=float(np.mean((y-pred)**2))
                        row[f"brier_genus_tau_{tau:g}"]=err
                        row[f"brier_gain_tau_{tau:g}"]=base_brier-err
                    out.append(row)
    return pd.DataFrame(out,columns=COLUMNS)
