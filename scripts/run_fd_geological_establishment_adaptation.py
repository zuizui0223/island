"""Original oceanic-versus-continental floral source assembly × endemic proxy.

GIFT/Direct, WCVP native-range, and independent Roeble/Zizka geology are
kept inside source-protected original CI. Source documentation is not
observed arrivals; endemic congenital floral differences are not adaptation.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.regional_native_provenance_volume import regional_native_candidates
from fd.geological_establishment_adaptation import geology_origin_flora_pilot
from fd.regional_native_flora_fd import TRAITS


def run(original:Path,wcvp_status:Path,island_tdwg:Path,crosswalk:Path,*,
        seed:int=20261114)->dict:
    ledger=pd.read_csv(
        original/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        dtype=str,usecols=["accepted_species","trait_name","normalized_value",
                           "evidence_scope","evidence_quality"])
    raw=pd.read_csv(
        original/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        dtype=str,usecols=["island_id","accepted_species","origin_status"])
    eligible,gate=strict_direct_928_flora(ledger,raw)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Original 928 Direct floral islands altered")
    data=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ]
    if data.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Source Direct accepted plant has conflicting trait states")
    direct=data.pivot(index="accepted_species",columns="trait_name",
          values="normalized_value").dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(direct)!=2234:
        raise ValueError("Original source only 2234 direct species allowed")
    status=pd.read_csv(wcvp_status,dtype=str,
        usecols=["accepted_species","native_l3_codes"]).drop_duplicates()
    tdwg=pd.read_csv(island_tdwg,dtype=str,
        usecols=["island_id","tdwg_l3_code","tdwg_match_status"])
    if len(tdwg)!=8265 or status.accepted_species.nunique()!=55137:
        raise ValueError("Archived WCVP source crosswalk changed")
    classified=regional_native_candidates(eligible,tdwg,status)
    n=classified.groupby("island_id")[
        ["strict_source_native","regional_native_compatible"]].sum()
    if int((n.strict_source_native>=30).sum())!=87 or (
        int((n.regional_native_compatible>=30).sum())!=357):
        raise ValueError("Original 87/357 source native flora cohorts changed")
    geo=pd.read_csv(
        original/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"],dtype=str)
    geo=geo.loc[geo.island_id.isin(eligible.island_id)].drop_duplicates("island_id")
    if len(geo)!=928:
        raise ValueError("Original geographic blocks unavailable")
    geology=pd.read_csv(crosswalk,dtype=str)
    if geology.island_id.nunique()<8000:
        raise ValueError("Independent Roeble/Zizka GSHHG geology coverage invalid")
    assignments=pd.read_csv(
        original/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"])
    assignments=assignments.loc[
        assignments.island_id.isin(eligible.island_id)&
        assignments.source_mode.eq("geo_k5")].copy()
    assignments["source_rank"]=pd.to_numeric(
        assignments.source_rank,errors="raise")
    assignments=assignments.loc[assignments.source_rank.between(1,5)]
    if assignments.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original mainland source candidates repeat")
    five=assignments.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(five.get(k,()))!=5 for k in classified.island_id.unique()):
        raise ValueError("No original 5-nearest source opportunity for some islands")
    donor=pd.read_csv(
        original/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"])
    donor=donor.loc[donor.work_species.isin(direct.index)].drop_duplicates()
    by=donor.groupby("entity_ID").work_species.agg(set).to_dict()
    out=geology_origin_flora_pilot(
        direct,classified,geo,geology,five,by,seed=seed,boots=1999)
    if out["n_original_all_928_Direct_island_floras"]!=928:
        raise ValueError("Original source 928 island support was lost")
    assert any("median_source_distinct_genera" in q for q in out["groups"])
    assert any(q["response"]=="expected_10genus_share" for q in out["geological_origin_associations"])
    assert len(out["geology_matched_geography_and_archipelago_diagnostics"])==6
    assert all(q["n_spatial_and_verified_archipelago_disjoint_pairs"]>=0
               for q in out["geology_matched_geography_and_archipelago_diagnostics"])
    assert all("n_islands_with_any_Direct_endemic" in q for q in out["groups"] if q["scope"]=="all_recorded")
    out["original_seed"]=seed
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--wcvp-status",type=Path,required=True)
    p.add_argument("--island-tdwg",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--seed",type=int,default=20261114)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=run(a.original_root,a.wcvp_status,a.island_tdwg,
          a.geology_crosswalk,seed=a.seed)
    text=json.dumps(r,sort_keys=True,indent=2,allow_nan=False)
    confidential=(
        "accepted_species","island_id","entity_ID","work_species",
        "spatial_block","source_signature","native_l3_codes",
        "tdwg_l3_code","geological_record_id","weigelt_record_id",
        "island_latitude","island_longitude","geological_reference",
        "source_genus","source_family","family_map")
    if any('"'+col+'"' in text for col in confidential):
        raise ValueError("Original taxonomic or geology/source IDs escaped")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n")
    print(json.dumps({
        "status":r["status"],
        "n_geology_labelled":r["n_geology_labelled_source_comparable_all_flora"],
        "n_provenance_origin_groups":len(r["groups"]),
        "n_source_relative_oceanic_minus_continental_tests":sum(
            x["origin_contrast_supported_sampling"]
            for x in r["geological_origin_associations"]),
        "endemic_with_source_congener_pairs":r[
            "n_endemic_original_direct_species_with_some_source_congener"],
        "geology_type":"Roeble/Zizka, NOT historical LGM mainland connectivity",
    },sort_keys=True))


if __name__=="__main__":main()
