"""Original 928-island WCVP/Direct native-status sensitivity for functional diversity.

All metrics are EXACT finite-population expectations of n=30 species-direct
flower morphology rarefactions, not Gaussian KDE. Evaluate marginal Rao Q,
joint Rao Q, and colour/architecture pairing departure versus an exact
fixed-colour-margin shuffle, plus two visitor-access shape proxies.

An original WCVP native TDWG-L3 RANGE is not local island nativeness.
Only geographic-regime, spatial-block, and GIFT-source-signature aggregates
may leave protected original-source CI. No pollinator guild loss is measured.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.regional_native_flora_fd import TRAITS, METRICS, exact_rarefied_flora_fd
from fd.regional_native_provenance_volume import (
    regional_native_candidates,
    frozen_seed,cluster_equal_bootstrap
)

SCOPES=("all_recorded","strict_source_native","regional_WCVP_compatible")


def source_direct_triples(direct:pd.DataFrame)->pd.DataFrame:
    needed={"accepted_species","trait_name","normalized_value",
            "evidence_scope","evidence_quality"}
    if not needed.issubset(direct):
        raise ValueError("Original source Direct flower trait schema changed")
    t=direct.loc[
        direct.trait_name.isin(TRAITS)&
        direct.evidence_scope.eq("species_direct")&
        direct.evidence_quality.isin(("high","medium"))&
        direct.normalized_value.notna()&
        ~direct.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if t.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Source trait evidence conflicting across original accepted species")
    z=t.pivot(index="accepted_species",columns="trait_name",
        values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(z)!=2234:
        raise ValueError("Original 2,234 Direct floral species cohort no longer frozen")
    return z


def _audit_group(frame:pd.DataFrame,column:str,cluster:str,*,
                 seed:int):
    return cluster_equal_bootstrap(
        frame,column,cluster,draws=1999,seed=seed)


def run(original_root:Path,wcvp_status:Path,island_tdwg:Path, *,
        seed:int=20261114)->dict:
    ledger=pd.read_csv(
        original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                            "evidence_scope","evidence_quality"])
    flower=source_direct_triples(ledger)
    flora=pd.read_csv(
        original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status"])
    eligible,gate=strict_direct_928_flora(ledger,flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original source 928 island IDs no longer stable")
    if eligible.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island plant accepted name duplicated")
    geographic=pd.read_csv(
        original_root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"])
    geographic=geographic.loc[geographic.island_id.isin(
        eligible.island_id)].drop_duplicates("island_id")
    if len(geographic)!=928 or geographic.spatial_block.isna().any():
        raise ValueError("Original spatial block/region evidence is incomplete")
    island_meta=geographic.set_index("island_id")
    ranges=pd.read_csv(wcvp_status,dtype=str,
        usecols=["accepted_species","native_l3_codes"]).drop_duplicates()
    mapping=pd.read_csv(island_tdwg,dtype=str,
        usecols=["island_id","tdwg_l3_code","tdwg_match_status"])
    if len(mapping)!=8265 or len(ranges.accepted_species.unique())!=55137:
        raise ValueError("Original WCVP/TDWG August 2026 published source drifted")
    classified=regional_native_candidates(eligible,mapping,ranges)
    if classified.island_id.nunique()!=928:
        raise ValueError("WCVP source propagation lost original 928 island flora")
    reported=classified.groupby("island_id")[
        ["strict_source_native","regional_native_compatible"]].sum()
    if int((reported.strict_source_native>=30).sum())!=87:
        raise ValueError("Source-original 87 strict native islands drifted")
    if int((reported.regional_native_compatible>=30).sum())!=357:
        raise ValueError("WCVP source 357 candidate-native-eligible island cohort drifted")

    src=pd.read_csv(
        original_root/"fixed/source/result/island_source_assignments.csv.gz",
        usecols=["island_id","source_mode","source_rank","entity_ID"],dtype=str)
    src=src.loc[
        src.island_id.isin(eligible.island_id)&
        src.source_mode.eq("geo_k5")].copy()
    src["source_rank"]=pd.to_numeric(src["source_rank"],errors="raise")
    src=src.loc[src.source_rank.between(1,5)]
    if src.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT nearest native source position repeated")
    near=src.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(near.get(i,()))!=5 or len(set(near[i]))!=5
           for i in classified.island_id.unique()):
        raise ValueError("Original GIFT five nearest candidate mainland slots incomplete")
    donor=pd.read_csv(
        original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    donor=donor.loc[donor.work_species.isin(flower.index)].drop_duplicates()
    donor_by_source=donor.groupby("entity_ID").work_species.agg(set).to_dict()
    taxa=set(flower.index)
    byisland={i:g for i,g in classified.groupby("island_id")}
    rows=[]
    for island in sorted(byisland):
        g=byisland[island]
        alltaxa=set(g.accepted_species)
        source=(set().union(*(donor_by_source.get(j,set()) for j in near[island]))&
                 taxa)-alltaxa
        source_score=exact_rarefied_flora_fd(
            flower.loc[sorted(source),list(TRAITS)].to_numpy(dtype=object),
            n=30) if len(source)>=30 else None
        tiers=(
           ("all_recorded",alltaxa),
           ("strict_source_native",set(g.loc[
                g.strict_source_native,"accepted_species"])),
           ("regional_WCVP_compatible",set(g.loc[
                g.regional_native_compatible,"accepted_species"])),
        )
        for scope,names in tiers:
            if len(names)<30:continue
            score=exact_rarefied_flora_fd(
                flower.loc[sorted(names),list(TRAITS)].to_numpy(dtype=object),
                n=30)
            row=dict(
                island_id=island,
                scope=scope,
                regime=island_meta.at[island,"analysis_regime"],
                spatial_block=island_meta.at[island,"spatial_block"],
                source_signature=near[island][0],
                original_direct_plant_records_in_scope=len(names),
                n_candidate_source_distinct_species=len(source))
            for m in METRICS:
                row[m]=score[m]
                row["source_disjoint_delta_"+m]=(
                    score[m]-source_score[m] if source_score is not None else float("nan"))
            rows.append(row)
    data=pd.DataFrame.from_records(rows)
    if len(data.loc[data.scope.eq("all_recorded")])!=928:
        raise ValueError("Original all-flora 928 cohort lost during evaluation")
    regional=[]
    for (region,scope),g in data.groupby(["regime","scope"]):
        d=g.loc[g.n_candidate_source_distinct_species>=30]
        entry=dict(region=region,scope=scope,n_islands=int(len(g)),
            n_source_comparable_islands=int(len(d)),
            n_original_geographic_blocks=int(d.spatial_block.nunique()),
            n_first_GIFT_candidate_source_signatures=int(d.source_signature.nunique()),
            median_original_direct_plant_species_in_scope=float(
                g.original_direct_plant_records_in_scope.median()))
        if len(d)>=10 and d.spatial_block.nunique()>=5 and d.source_signature.nunique()>=5:
            for m in METRICS:
                target="source_disjoint_delta_"+m
                for label,grp in (("block","spatial_block"),
                                  ("source","source_signature")):
                    entry[label+"_"+target]=_audit_group(
                        d,target,grp,seed=frozen_seed(seed,"regional",
                                                     region,scope,m,label))
        regional.append(entry)
    paired=[]
    for left,right in (
        ("regional_WCVP_compatible","all_recorded"),
        ("strict_source_native","all_recorded"),
        ("regional_WCVP_compatible","strict_source_native")
    ):
        a=data.loc[data.scope.eq(left)].set_index("island_id")
        b=data.loc[data.scope.eq(right)].set_index("island_id")
        for region,geo in geographic.groupby("analysis_regime"):
            ids=a.index.intersection(b.index).intersection(geo.island_id)
            if len(ids)<10:continue
            q=a.loc[ids].copy()
            q=q.loc[q["source_disjoint_delta_marginal_rao"].notna()].copy()
            ids=q.index
            if len(ids)<10 or q.spatial_block.nunique()<5 or q.source_signature.nunique()<5:
                continue
            for m in METRICS:
                q["paired_"+m]=(
                    a.loc[ids,"source_disjoint_delta_"+m]-
                    b.loc[ids,"source_disjoint_delta_"+m]).to_numpy(float)
            rec=dict(region=region,contrast=left+"_minus_"+right,
                n_identical_physical_islands=len(ids),
                n_original_spatial_blocks=int(q.spatial_block.nunique()),
                n_first_GIFT_source_signatures=int(q.source_signature.nunique()))
            for m in METRICS:
                for label,grp in (("block","spatial_block"),
                                  ("source","source_signature")):
                    rec[label+"_paired_"+m]=_audit_group(
                        q,"paired_"+m,grp,seed=frozen_seed(
                            seed,"paired",region,left,right,m,label))
            paired.append(rec)
    return dict(
       status="ORIGINAL_DIRECT_WCVP_NATIVE_FLORA_FINITE_RAO_AND_PAIRING",
       n_original_direct_plant_species_global=2234,
       n_actual_original_three_trait_profiles=100,
       n_original_eligible_island_floras=928,
       n_original_strict_native_eligible_islands=87,
       n_WCVP_region_compatible_candidate_eligible_islands=357,
       n_original_island_TDWG_L3_records=8265,
       n_original_WCVP_range_source_species=55137,
       n_direct_flora_source_species_with_available_regional_range=int(
           len(set(classified.accepted_species)&set(ranges.accepted_species))),
       n_original_direct_wcvp_eligible_rows=int(len(classified)),
       n_original_direct_source_native_rows=int(classified.strict_source_native.sum()),
       n_unresolved_WCVP_compatible_upgrade_rows=int(
           classified.regional_WCVP_upgrade_only.sum()),
       numerical_seed=seed,
       method=dict(
          equal_original_Direct_species_represented_by_analytic_rarefaction=30,
          original_species_direct_all_morphology_complete_only=True,
          exact_expected_rao_under_random_without_replacement=True,
          no_species_sampling_monte_carlo_noise=True,
          fixed_colour_marginal_pairing_departure_analytic=True,
          intact_original_form_and_symmetry_pairs=True,
          both_original_same_island_and_original_source_compared=True,
          original_same_island_species_excluded_from_GIFT_mainland_donors=True,
          original_source_GIFT_introduced_taxa_never_upgraded=True,
          WCVP_TDWG_L3_is_regionally_compatible_not_exact_island_native=True,
          two_original_cluster_resampling_draws_per_region=1999,
          original_inferred_absence_of_pollinators_never_permitted=True,
          no_within_plant_original_island_trait_evolution_estimated=True,
       ),
       region_anonymous_native_flora_functional_diversity=regional,
       same_exact_physical_island_provenance_sensitivity=paired,
       original_source_plant_island_region_taxon_records_not_published=True,
       no_actual_pollinator_guild_change_or_pollen_service_measured=True,
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--wcvp-status",type=Path,required=True)
    p.add_argument("--island-tdwg",type=Path,required=True)
    p.add_argument("--seed",type=int,default=20261114)
    p.add_argument("--output",type=Path,required=True)
    z=p.parse_args()
    data=run(z.original_root,z.wcvp_status,z.island_tdwg,
             seed=z.seed)
    text=json.dumps(data,indent=2,sort_keys=True,allow_nan=False)
    if any('"'+k+'"' in text for k in (
        "island_id","accepted_species","entity_ID","work_species",
        "source_signature","spatial_block","longitude","latitude",
        "native_l3_codes","tdwg_l3_code","genus","family")):
        raise ValueError("Rights-protected source taxon or island ID leaked")
    z.output.parent.mkdir(parents=True,exist_ok=True)
    z.output.write_text(text+"\n")
    print(json.dumps({
        "status":data["status"],
        "n_original_928_islands":data["n_original_eligible_island_floras"],
        "n_original_source_native_30plus":data["n_original_strict_native_eligible_islands"],
        "n_regional_native_candidate_30plus":data["n_WCVP_region_compatible_candidate_eligible_islands"],
        "regions":[{
          "region":r["region"],"scope":r["scope"],
          "n_original":r["n_islands"],"n_GIFT_source_comparable":r["n_source_comparable_islands"],
          "n_blocks":r["n_original_geographic_blocks"],
          "n_GIFT_sources":r["n_first_GIFT_candidate_source_signatures"],
        } for r in data["region_anonymous_native_flora_functional_diversity"]],
    },sort_keys=True))


if __name__=="__main__":main()
