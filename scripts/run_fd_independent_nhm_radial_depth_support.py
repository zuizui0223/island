"""Independent NHM text trait × Direct depth source-only feasibility audit.

For 4,505 original GIFT botanical island floras, can independently extracted
NHM symmetry or form plus original (NOT externally remeasured) tube-depth
trait in the SAME different-name mainland candidate species support a regional
actual-source comparison of island A×deep morphology?

Strictly an identification/data opportunity audit, NOT an independent two-axis
replication, pollinator syndrome validation or actual visitor shift.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd_nhm_original_provider_audit import (
    PROVIDERS,BLOB_PINS,read_provider_file,consensus_provider_profiles,
)
from fd_nhm_external_flora_preflight import nhm_git_blob_sha,nhm_external_categories
from fd.real_joint_syndrome_access import make_joint_species
from fd.source_genus_same_island import genus

FLORAL_PAIRS={
    "actinomorphic_deep":("symmetry",frozenset({"actinomorphic"})),
    "zygomorphic_deep":("symmetry",frozenset({"zygomorphic"})),
    "tubular_bell_deep":("form",frozenset({"tubular","bell_campanulate"})),
}
TIERS=("one_or_more_consistent_NHM_corpora",
       "two_distinct_consistent_NHM_corpora")
MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
HASHES={
"direct_only_trait_ledger.csv.gz":
 "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
"chapter1_status_flora.csv.gz":
 "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}


def digest(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def external_hybrid_donors(external_axis:dict[str,str],
                           original_depth:dict[str,str],
                           accepted_gift_names:set[str],
                           every_original_island_plant_name:set[str],
                           selected_positive:frozenset[str]):
    """Only NHM-external A and original-Direct B on SAME external plant name.

    NO donor may occur in ANY of 4,505 original island floras, ever.
    A name with missing depth is not a negative paired state.
    """
    if not isinstance(selected_positive,frozenset):
        raise ValueError("explicit external phenotype state required")
    candidates=(set(external_axis)&set(original_depth)&accepted_gift_names) - every_original_island_plant_name
    if candidates & every_original_island_plant_name:
        raise AssertionError("original island plant leaked into external mainland donor")
    scored={name:int(external_axis[name] in selected_positive and
                     original_depth[name]=="deep") for name in candidates}
    return scored


def original_feasibility(
    ledger:pd.DataFrame,flora:pd.DataFrame,geo:pd.DataFrame,
    mainland:pd.DataFrame,assign:pd.DataFrame,
    NHM:dict[str,pd.DataFrame],
)->dict:
    full=status_groups(flora)
    if len(full)!=4505:raise ValueError("original GIFT 4,505 islands contracted")
    if geo.island_id.duplicated().any() or geo.island_id.nunique()<8200:
        raise ValueError("GSHHG original geographic source invalid")
    georef=geo.set_index("island_id").reindex(sorted(full))
    if int(georef.analysis_regime.notna().sum())!=4504:
        raise ValueError("original GSHHG island reconciliation drift")
    if mainland.duplicated(["entity_ID","work_species"]).any() or (
        mainland.entity_ID.nunique()!=748):
        raise ValueError("original GIFT native mainland source invalid")
    GIFT=set(mainland.work_species)
    original_islands=set().union(*(v["full"] for v in full.values()))
    traits=strict_traits(ledger)
    paired=make_joint_species(traits)
    original_depth=traits["tube_depth_class"]

    by_provider={}
    for key,frame in NHM.items():
        by_provider[key]=nhm_external_categories(frame)[0]
    ext={}
    external_audit={}
    for axis in ("symmetry","form"):
        consensus,confirmed,conflict,_=consensus_provider_profiles(
            by_provider,axis=axis)
        ext[axis]={"one_or_more_consistent_NHM_corpora":consensus,
                   "two_distinct_consistent_NHM_corpora":confirmed}
        direct=traits["floral_symmetry" if axis=="symmetry" else "floral_form"]
        compared=set(consensus)&set(direct)
        external_audit[axis]=dict(
            n_original_consensus_external_names=len(consensus),
            n_original_two_corpus_agreement_names=len(confirmed),
            n_original_conflicted_external_species=conflict,
            n_original_external_and_Direct_axis_same_name=len(compared),
            n_agreement_on_comparable_accepted_names=sum(
                consensus[s]==direct[s] for s in compared
            ),
            n_external_GIFT_source_only_not_on_ANY_original_island=len(
                (set(consensus)&GIFT)-original_islands
            ),
            n_external_GIFT_source_only_depth_Direct_same_name=len(
                (set(consensus)&set(original_depth)&GIFT)-original_islands
            ),
        )

    pool={}
    original_donor_records=[]
    for name,(axis,positive) in FLORAL_PAIRS.items():
        for tier in TIERS:
            result=external_hybrid_donors(
                ext[axis][tier],original_depth,GIFT,original_islands,positive
            )
            pool[(name,tier)]=result
            original_donor_records.append(dict(
                sector=name,external_source_provider_tier=tier,
                n_original_same_name_external_A_and_original_Direct_depth_B=len(result),
                n_independent_external_Direct_depth_axis=0,
                n_original_hybrid_source_positive_joint=sum(result.values()),
                n_original_hybrid_source_botanical_genera=len({
                    genus(x) for x in result
                }),
                n_external_same_name_gift_accepted_sources_not_local_evolution=True
            ))

    A=assign.loc[assign.source_mode.isin(MODES)].copy()
    A.source_rank=pd.to_numeric(A.source_rank,errors="raise").astype(int)
    if A.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("original candidate source region rank duplicated")
    if {mode:int(A.source_mode.eq(mode).sum()) for mode in MODES}!={
        "geo_k5":41325,"geo50_climate10":82190
    }:
        raise ValueError("original frozen source assignment opportunities changed")
    for mode,(_,rank) in MODES.items():
        if not A.loc[A.source_mode.eq(mode),"source_rank"].between(1,rank).all():
            raise ValueError("invalid original source regional rank")
    byisland=A.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True).entity_ID.agg(tuple).to_dict()
    used=set().union(*(set(p) for p in pool.values()))
    external_mainland=mainland.loc[mainland.work_species.isin(used)]
    byregion=external_mainland.groupby("entity_ID").work_species.agg(set).to_dict()

    counters=defaultdict(Counter)
    original_block_support=defaultdict(set)
    original_source_support=defaultdict(set)
    for island,q in full.items():
        reg=georef.at[island,"analysis_regime"]
        reg=str(reg) if pd.notna(reg) else "geography_unresolved"
        if reg not in REGIMES:raise ValueError("original geographical regime unknown")
        block=georef.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        sources={}
        for mode,(label,ranks) in MODES.items():
            regional=byisland.get((island,mode),())
            if regional and len(regional)!=ranks:
                raise ValueError("partial source regional donor list")
            sources[label]=set().union(*(byregion.get(x,set()) for x in regional))
        for sector in FLORAL_PAIRS:
            target_known=paired[sector]["known"]
            for scope,roster in (
                ("all_recorded",q["full"]),("strict_native",q["native"]),
                ("strict_native_nonendemic",q["native_nonendemic"]),
            ):
                target=roster&target_known
                n=len(target)
                for tier in TIERS:
                    key=(reg,scope,sector,tier)
                    c=counters[key]
                    c["n_full_region_islands"]+=1
                    if n>=5:c["n_at_least_five_real_island_pair_carriers"]+=1
                    qsrc={mode:(source&set(pool[(sector,tier)]))
                          for mode,source in sources.items()}
                    for mode in ("geo5","climate10"):
                        nd=len(qsrc[mode])
                        for limit in (1,5,10):
                            if nd>=limit:c[f"n_source_{mode}_hybrid_donors_ge_{limit}"]+=1
                    dgeo={genus(s) for s in qsrc["geo5"]}
                    dcli={genus(s) for s in qsrc["climate10"]}
                    shared=dgeo&dcli
                    eligible={s for s in target if genus(s) in shared}
                    if n>=5 and len(qsrc["geo5"])>=1 and len(qsrc["climate10"])>=1:
                        c["n_islands_with_5plus_island_true_pair_and_both_hybrid_source_donor"]+=1
                    if n>=5 and len(qsrc["geo5"])>=10 and len(qsrc["climate10"])>=10:
                        c["n_islands_with_5plus_island_true_pair_and_both_source10"]+=1
                    if len(eligible)>=5 and len({genus(s) for s in eligible})>=2 and (
                        len(qsrc["geo5"])>=10 and len(qsrc["climate10"])>=10):
                        c["n_islands_with_5plus_same_island_target_congeners_and_2plus_genera_both_source10"]+=1
                        if block is not None:original_block_support[key].add(block)
                        near=byisland.get((island,"geo_k5"),())
                        if near:original_source_support[key].add(near[0])

    rows=[]
    for region in REGIMES:
        for scope in SCOPES:
            for sector in FLORAL_PAIRS:
                for tier in TIERS:
                    key=(region,scope,sector,tier)
                    c=counters[key]
                    nqual=c[
                        "n_islands_with_5plus_same_island_target_congeners_and_2plus_genera_both_source10"
                    ]
                    nblock=len(original_block_support[key])
                    nsource=len(original_source_support[key])
                    rows.append(dict(
                        region=region,scope=scope,structural_pair=sector,
                        source_provider_tier=tier,
                        **{v:int(c[v]) for v in (
                            "n_full_region_islands",
                            "n_at_least_five_real_island_pair_carriers",
                            "n_source_geo5_hybrid_donors_ge_1",
                            "n_source_geo5_hybrid_donors_ge_5",
                            "n_source_geo5_hybrid_donors_ge_10",
                            "n_source_climate10_hybrid_donors_ge_1",
                            "n_source_climate10_hybrid_donors_ge_5",
                            "n_source_climate10_hybrid_donors_ge_10",
                            "n_islands_with_5plus_island_true_pair_and_both_hybrid_source_donor",
                            "n_islands_with_5plus_island_true_pair_and_both_source10",
                            "n_islands_with_5plus_same_island_target_congeners_and_2plus_genera_both_source10",
                        )},
                        n_independent_original_10degree_blocks_on_final_support=nblock,
                        n_original_first_GIFT_sources_on_final_support=nsource,
                        status=("SUPPORT_ONLY_HYBRID_SOURCE_COMPARISON_POSSIBLE"
                                if nqual>=30 and nblock>=10 and nsource>=5
                                else "NOT_IDENTIFIABLE_EXTERNAL_FULL_SYNDROME_VALIDATION"),
                    ))
    if len(rows)!=len(REGIMES)*len(SCOPES)*len(FLORAL_PAIRS)*len(TIERS):
        raise ValueError("original region×native×source hybrid support collapsed")
    return dict(
        status="ORIGINAL_GIFT_4505_EXTERNAL_NHM_SYMMETRY_FORM_PLUS_DIRECT_DEPTH_COVERAGE_AUDIT",
        n_original_GIFT_island_floras=4505,
        n_original_native_mainland_source_regions=748,
        n_original_single_species_real_symmetry_depth_pairs={
            k:len(paired[k]["known"]) for k in FLORAL_PAIRS
        },
        original_external_NHM_category_quality=external_audit,
        hybrid_external_source_only_donor_counts=original_donor_records,
        regional_native_source_overlap=rows,
        n_regional_status_sector_provider_panels=len(rows),
        n_hybrid_source_panels_meeting_actual_comparison_gates=sum(
            r["status"]=="SUPPORT_ONLY_HYBRID_SOURCE_COMPARISON_POSSIBLE"
            for r in rows
        ),
        limitations=dict(
            NHM_does_not_supply_independently_recorded_tube_depth=True,
            mixed_provider_symmetry_form_and_original_Direct_tube_depth_not_independent_pair=True,
            original_4505_islands_retained_as_denominator=True,
            no_missing_gift_plant_trait_or_absence_imputed=True,
            source_only_accepted_species_absent_from_ANY_island_flora=True,
            outside_source_provider_is_independent_release_not_independent_field_study=True,
            actual_pollinator_guild_or_effective_pollen_not_measured=True,
            score_or_effect_never_fitted_below_identification_gates=True,
        ),
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--nhm-directory",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    ledger_file=a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_file=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for f in (ledger_file,flora_file):
        if digest(f)!=HASHES[f.name]:
            raise ValueError("source-protected original botanical source hash changed")
    ledger=pd.read_csv(ledger_file,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    flora=pd.read_csv(flora_file,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    geo=pd.read_csv(a.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assign=pd.read_csv(
        a.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(
        a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    NHM={}
    for label,name in PROVIDERS.items():
        f=a.nhm_directory/name
        if not f.is_file() or nhm_git_blob_sha(f)!=BLOB_PINS[label]:
            raise ValueError("NHM source independent published blob identity changed")
        NHM[label]=read_provider_file(f,label)
    result=original_feasibility(
        ledger,flora,geo,mainland,assign,NHM
    )
    result["original_Direct_species_ledger_sha256"]=digest(ledger_file)
    result["original_GIFT_flora_sha256"]=digest(flora_file)
    body=json.dumps(result,indent=2,ensure_ascii=False,
                    allow_nan=False,sort_keys=True)
    for private in (
        "island_id","accepted_species","entity_ID","work_species",
        "source_signature","spatial_block","source_genus","source_family",
        "latitude","longitude","GBIF_id","photo_id","taxon_id",
    ):
        if '"'+private+'"' in body:
            raise ValueError("source botanist/taxon/island protected identifier leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=result["status"],original_islands=4505,
        provider_comparison_panels=result[
            "n_regional_status_sector_provider_panels"],
        n_identifiable_hybrid_external_source_panels=result[
            "n_hybrid_source_panels_meeting_actual_comparison_gates"],
        independently_remeasured_depth=False,
    ),sort_keys=True))


if __name__=="__main__":
    main()
