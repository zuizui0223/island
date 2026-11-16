"""Protected original 4505 GIFT island NHM numeric tube genus concentration audit.

Original Direct island deep labels are used only to SCORE predictions, NEVER
to calibrate their own botanical genus. The source donor score is trained from
different accepted mainland species' NHM numerical cm. Exports aggregate
region/status/corpus summaries only, not island, taxon or genus identities.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import pandas as pd
import numpy as np

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd.nhm_numeric_source_congener_predict import (
    make_genus_heldout_length_calibrator,calibrated_prob,botanical_genus,
)
from fd.nhm_numeric_genus_concentration import (
    SCOPES,TIERS,REGIMES,MODES,audit_genus_concentration,
)
from run_fd_independent_nhm_numeric_tube_length import (
    parse_nhm_numeric_file,combine_external_provider_lengths,
)
from fd_nhm_original_provider_audit import PROVIDERS,BLOB_PINS
from fd_nhm_external_flora_preflight import nhm_git_blob_sha

PIN={
    "direct_only_trait_ledger.csv.gz":
    "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
    "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
SOURCE_MODES={"geo_k5":("geo5",5),"geo50_climate10":("climate10",10)}


def file_sha(p:Path)->str:
    h=sha256()
    with p.open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):
            h.update(z)
    return h.hexdigest()


def original_genus_private_rows(
    ledger:pd.DataFrame,flora:pd.DataFrame,geography:pd.DataFrame,
    assignments:pd.DataFrame,mainland:pd.DataFrame,
    original_numeric:dict[str,dict],
)->tuple[pd.DataFrame,pd.DataFrame,dict]:
    botanical=status_groups(flora)
    if len(botanical)!=4505:
        raise ValueError("original 4505 GIFT accepted island flora universe unavailable")
    if geography.island_id.duplicated().any() or geography.island_id.nunique()<8200:
        raise ValueError("original GSHHG crosswalk incomplete")
    G=geography.set_index("island_id").reindex(sorted(botanical))
    if int(G.analysis_regime.notna().sum())!=4504:
        raise ValueError("original geographically unknown island incorrectly removed")
    if mainland.duplicated(["entity_ID","work_species"]).any() or (
        mainland.entity_ID.nunique()!=748):
        raise ValueError("original 748 GIFT native source entities changed")
    gift=set(mainland.work_species)
    traits=strict_traits(ledger)
    depth=traits["tube_depth_class"]
    lengths_all,lengths_two,prov=combine_external_provider_lengths(original_numeric)
    cal,cal_meta=make_genus_heldout_length_calibrator(
        lengths_all,depth,gift
    )
    source=assignments.loc[assignments.source_mode.isin(SOURCE_MODES)].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise").astype(int)
    if source.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("duplicate original GIFT source region ranking")
    if {m:int(source.source_mode.eq(m).sum()) for m in SOURCE_MODES}!={
        "geo_k5":41325,"geo50_climate10":82190,
    }:
        raise ValueError("original geographical/climate source assignment drift")
    for mode,(_,k) in SOURCE_MODES.items():
        if not source.loc[source.source_mode.eq(mode),"source_rank"].between(1,k).all():
            raise ValueError("out-of-range original mainland source ranking")
    choices=source.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    bysource=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    @lru_cache(maxsize=16)
    def pooled(ids:tuple[str,...])->frozenset[str]:
        return frozenset().union(*(bysource.get(k,set()) for k in ids))

    private=[]
    rosters=[]
    n_climate=0
    for island in sorted(botanical):
        item=botanical[island]
        full=item["full"]
        reg=G.at[island,"analysis_regime"]
        region=str(reg) if pd.notna(reg) else "geography_unresolved"
        if region not in REGIMES:
            raise ValueError("invalid original GSHHG region")
        block=G.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        geo5_sig=choices.get((island,"geo_k5"),())
        first=geo5_sig[0] if len(geo5_sig)==5 else None
        rosters.append(dict(island_id=island,regime=region))
        sources={}
        for mode,(name,expected) in SOURCE_MODES.items():
            ids=choices.get((island,mode),())
            if ids and len(ids)!=expected:
                raise ValueError("partial original native mainland candidate source")
            sources[name]=(set(pooled(ids))-full) if ids else set()
            if sources[name]&full:
                raise AssertionError("focal island accepted names in independent source")
        if choices.get((island,"geo50_climate10")):
            n_climate+=1
        scopes={
            "all_recorded":full,
            "strict_native":item["native"],
            "strict_native_nonendemic":item["native_nonendemic"],
        }
        original_known_genera={
            botanical_genus(sp) for sp in full
            if sp in depth and depth[sp] in ("shallow","intermediate","deep")
        }
        if not original_known_genera:
            continue
        for tier,lengths in zip(TIERS,(lengths_all,lengths_two)):
            numeric_by_mode={}
            all_count_by_mode={}
            for mode in MODES:
                by_genus=defaultdict(list)
                counts=Counter()
                for species in sources[mode]:
                    g=botanical_genus(species)
                    if g not in original_known_genera:continue
                    counts[g]+=1
                    if species in lengths:
                        fit=cal(g)
                        if not fit["heldout_target_genus_never_in_calibration"]:
                            raise AssertionError("source NHM calibration used target-genus trait")
                        by_genus[g].append(calibrated_prob(
                            float(lengths[species]),fit
                        ))
                numeric_by_mode[mode]=by_genus
                all_count_by_mode[mode]=counts
            supported=set(numeric_by_mode["geo5"])&set(
                numeric_by_mode["climate10"]
            )
            if not supported:continue
            for scope,names in scopes.items():
                target=defaultdict(list)
                orig_count=Counter(botanical_genus(sp) for sp in names)
                for sp in names:
                    if sp in depth and depth[sp] in (
                        "shallow","intermediate","deep"
                    ):
                        g=botanical_genus(sp)
                        if g in supported:
                            target[g].append(int(depth[sp]=="deep"))
                for g,labels in target.items():
                    row=dict(
                        island_id=island,regime=region,
                        scope=scope,tier=tier,
                        genus=g,spatial_block=block,first_source=first,
                        n_target=len(labels),n_target_deep=sum(labels),
                        n_target_original_genus=orig_count[g],
                    )
                    for mode in MODES:
                        p=numeric_by_mode[mode][g]
                        row[mode+"__n_numeric"]=len(p)
                        row[mode+"__n_original_source_genus"]=all_count_by_mode[mode][g]
                        row[mode+"__genus_pdeep"]=float(np.mean(p))
                    private.append(row)
    if n_climate!=4485:
        raise ValueError("original climate10 source opportunity no longer 4485")
    if len(rosters)!=4505:
        raise AssertionError("full 4505 accepted GIFT flora denominator missing")
    if not private:
        raise ValueError("no original source-congener numeric cohort")
    return pd.DataFrame(private),pd.DataFrame(rosters),dict(
        **cal_meta,original_NHM_numeric_provider_summary=prov,
        n_NHM_numeric_accepted_species=len(lengths_all),
        n_NHM_multi_corpus_consistent_numeric_accepted_species=len(lengths_two),
    )


def compare_original_anonymous_reference(
    new:dict,previous:dict,*,max_abs=1e-9
)->dict:
    """Validate newest source genus weighted species-equal score matches
    previously successful original independent run on EXACT same cohorts.
    Never import private IDs from old report (it is anonymous).
    """
    old={(p["region"],p["scope"],p["tier"]):p for p in previous["panels"]}
    diffs=[]
    for p in new["panels"]:
        key=(p["region"],p["scope"],p["tier"])
        q=old[key]
        a=p["all_original_eligible__spatial_block"]
        b=q["all__spatial_block"]
        if a["n_original_floras"]!=q["n_same_island_both_source_mode_qualified"]:
            raise ValueError("source-validated original island sample changed")
        if a["n_original_floras"]!=b["n_islands"]:
            raise ValueError("original grouped geographic cohort inconsistent")
        if a["n_geography_or_source_groups"]!=b["n_groups"]:
            raise ValueError("original original geographic blocks not identical")
        for mode in MODES:
            if a["status"]=="DESCRIPTIVE_GENUS_CONCENTRATION_SENSITIVITY":
                x=a[mode+"__brier_gain_species_equal__block_mean"]
                y=b[mode+"__brier_gain__equal_group_mean"]
                diff=abs(x-y)
                if diff>max_abs:
                    raise ValueError(
                        f"new genus concentration does not reproduce frozen original score: {key} {mode} Δ={diff}"
                    )
                diffs.append(diff)
    return dict(
        original_report_status=previous["status"],
        n_original_paired_region_status_tier_panels_compared=len(new["panels"]),
        n_original_numeric_gains_machine_precision_matched=len(diffs),
        max_original_score_absolute_difference=float(max(diffs)) if diffs else None,
        matching_same_source_cohort_and_original_geographic_blocks=True,
    )


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-root",type=Path,required=True)
    ap.add_argument("--geology-crosswalk",type=Path,required=True)
    ap.add_argument("--nhm-directory",type=Path,required=True)
    ap.add_argument("--original-anonymous-benchmark",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    lp=a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    fp=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for f in (lp,fp):
        if file_sha(f)!=PIN[f.name]:
            raise ValueError("original SHA-verified accepted GIFT flora or traits changed")
    ledger=pd.read_csv(lp,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(fp,dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"
    ])
    geography=pd.read_csv(a.geology_crosswalk,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"
    ])
    assignments=pd.read_csv(
        a.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"]
    )
    mainland=pd.read_csv(
        a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"]
    )
    numeric={}
    for source,filename in PROVIDERS.items():
        path=a.nhm_directory/filename
        if not path.exists() or nhm_git_blob_sha(path)!=BLOB_PINS[source]:
            raise ValueError("NHM five provider numeric source Git blob changed")
        numeric[source],_=parse_nhm_numeric_file(path,source)
    full,roster,source_meta=original_genus_private_rows(
        ledger,flora,geography,assignments,mainland,numeric
    )
    result=audit_genus_concentration(
        full,roster,population=4505,bootstrap=999
    )
    old=json.loads(a.original_anonymous_benchmark.read_text())
    if (old.get("status")!="ORIGINAL_4505_NHM_NUMERIC_LENGTH_SOURCE_CONGENER_DEEP_PREDICTION"
        or old.get("n_original_botanical_flora_islands")!=4505
        or old.get("n_original_GIFT_native_mainland_source_entities")!=748):
        raise ValueError("original successful source benchmark does not match frozen original GIFT flora")
    result["exact_original_release_equivalence"]=compare_original_anonymous_reference(
        result,old
    )
    result.update(
        original_source_numeric_provenance=source_meta,
        n_original_GIFT_mainland_source_entities=748,
        original_GIFT_plant_traits_sha256=file_sha(lp),
        original_botanical_GIFT_flora_sha256=file_sha(fp),
        original_NHM_corpus_blob_guards=BLOB_PINS,
        source_all_plant_and_island_identifiers_not_released=True,
    )
    raw=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,
                   allow_nan=False)
    forbidden=("island_id","accepted_species","work_species","entity_ID",
               "source_signature","spatial_block","source_genus","source_family",
               "latitude","longitude","GBIF_id","photo_id","taxon_id","genus")
    if any('"'+key+'"' in raw for key in forbidden):
        raise ValueError("private original botanical source/plant/island identities leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=result["status"],original_islands=4505,
        n_region_status_tier_panels=len(result["panels"]),
        n_paired_original_gains_validated=result[
            "exact_original_release_equivalence"
        ]["n_original_numeric_gains_machine_precision_matched"],
        n_tropical_source_supported=next(
            z["n_original_islands_source_numeric_paired_support"]
            for z in result["panels"]
            if z["region"]=="tropical" and z["scope"]=="all_recorded"
            and z["tier"]==TIERS[0]
        ),
    ),sort_keys=True))


if __name__=="__main__":
    main()
