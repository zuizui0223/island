"""NHM per-corpus provenance for outside-ALL-4505-island GIFT mainland flora.

The NHM combined.csv does not expose a row-level original provider identifier.
Independently SHA-pin the FIVE published constituent sources and rebuild the
joint flower colour and form evidence without treating duplicate documents
as independent studies or pretending text traits are mainland field measures.

Outputs are aggregate only; no original plant, island, genus, or GIFT source ID.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from collections import defaultdict
import pandas as pd
from fd_island_functional_beta import strict_traits
from fd_nhm_external_flora_preflight import (
    nhm_git_blob_sha,nhm_external_categories,
    single_description_color,accepted_exact_binomial,
    EXTERNAL_COLUMNS,NHM_COMMIT,NHM_LICENSE,
)

PROVIDERS={
    "BHL_text":"bhl.csv",
    "Ecoflora_text":"ecoflora.csv",
    "FloraChina_text":"efloras.flora_of_china.csv",
    "FloraNorthAmerica_text":"efloras.flora_of_north_america.csv",
    "FloraPakistan_text":"efloras.flora_of_pakistan.csv",
}
BLOB_PINS={
    "BHL_text":"f46e15b95e92b1daad5c0ebf2145df73173505be",
    "Ecoflora_text":"dbdc6b4469583743b8b35e26912d21d05b550f32",
    "FloraChina_text":"ec87a99eee2cc1bf9a8868e1d5425654e7974897",
    "FloraNorthAmerica_text":"3710b08d7b041205a6cb5e7b811c08de317e006d",
    "FloraPakistan_text":"7f506417a3ae31fc4daa959a27976ed62c7b1449",
}
EXTERNAL_COLS=("taxon","flower colour","flower shape","flower symmetry",
               "corolla colour","petal colour")


def read_provider_file(path:Path,label:str):
    head=set(pd.read_csv(path,nrows=0).columns)
    if "taxon" not in head:raise ValueError("NHM source must contain taxon")
    use=sorted(set(EXTERNAL_COLS)&head)
    f=pd.read_csv(path,dtype=str,keep_default_na=False,low_memory=False,usecols=use)
    for col in EXTERNAL_COLS:
        if col not in f:f[col]=""
    f["source"]=label
    f["wfo"]=""
    return f[list(EXTERNAL_COLUMNS)]


def consensus_provider_profiles(by_source:dict,*,axis="primary"):
    """Keep only names with observed provider-given category AGREEMENT."""
    report=defaultdict(list)
    for provider,axes in by_source.items():
        for name,state in axes[axis].items():
            report[name].append((provider,state))
    consensus={}
    two_plus={}
    conflicted=0
    for taxon,entries in report.items():
        values={state for _,state in entries}
        if len(values)>1:
            conflicted+=1
            continue
        value=next(iter(values))
        consensus[taxon]=value
        if len({provider for provider,_ in entries})>=2:
            two_plus[taxon]=value
    return consensus,two_plus,conflicted,report


def validate_external_labels(by_provider:dict,original:dict,
                             original_islands:set,
                             gift:dict):
    donor=set(gift)-original_islands
    summary={}
    axis_report={}
    for axis in ("primary","expanded","form","symmetry"):
        same,twoplus,conflicted,_=consensus_provider_profiles(by_provider,axis=axis)
        matched=set(same)&set(original.get(axis,{}))
        agreement=sum(same[n]==original[axis][n] for n in matched)
        summary[axis]=dict(
            n_independent_NHM_provider_unequivocal_names=len(same),
            n_conflicting_between_provider_taxon_names=conflicted,
            n_consistent_two_or_more_source_corpus_names=len(twoplus),
            n_in_GIFT_mainland_not_in_original_4505_islands=len(set(same)&donor),
            n_two_or_more_corpora_in_GIFT_outside_original_4505=len(set(twoplus)&donor),
            n_exact_original_Direct_species_names_comparable=len(matched),
            n_agree_original_direct_category=agreement,
            agreement_if_comparable=agreement/len(matched) if matched else None,
        )
    return summary


def analyze(raw:dict,trait_original:dict,
            flora:pd.DataFrame,mainland:pd.DataFrame):
    if not {"accepted_species","island_id"}<=set(flora) or not (
        {"entity_ID","work_species"}<=set(mainland)):
        raise ValueError("Original protected botanical flora/mainland incomplete")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island species duplicates")
    n_islands=flora.island_id.nunique()
    if n_islands!=4505:raise ValueError("Original 4505 island data missing")
    all_island=set(flora.accepted_species)
    gift=set(mainland.work_species)
    by_provider={}
    provider_info=[]
    for provider,frame in raw.items():
        states,meta=nhm_external_categories(frame)
        by_provider[provider]=states
        raw_name=set(states["primary"])
        provider_info.append(dict(
            source_corpus=provider,
            n_original_description_rows=len(frame),
            n_exact_source_botanical_species=meta["n_exact_taxon_binomial_names"],
            n_strict_primary_flower_colour_species=len(raw_name),
            n_strict_organ_colour_expanded_species=len(states["expanded"]),
            n_shape_species=len(states["form"]),
            n_symmetry_species=len(states["symmetry"]),
            n_primary_colour_GIFT_source_only_names=len((raw_name&gift)-all_island),
            n_expanded_colour_GIFT_source_only_names=len(
                (set(states["expanded"])&gift)-all_island),
            n_singleton_source_nonmissing_colour_conflicts=meta[
                "n_species_with_documentation_conflict_by_axis"]["primary"],
        ))
    original={
        "primary":trait_original["flower_primary_color"],
        "expanded":trait_original["flower_primary_color"],
        "form":trait_original["floral_form"],
        "symmetry":trait_original["floral_symmetry"],
    }
    joined=validate_external_labels(by_provider,original,all_island,
                                    {name:1 for name in gift})
    return dict(
        status="NHM_PROVIDER_ORIGINAL_TEXT_SOURCE_INDEPENDENCE_AUDIT",
        n_original_gift_island_flora=n_islands,
        n_original_distinct_island_accepted_species=len(all_island),
        n_original_mainland_gift_regions=mainland.entity_ID.nunique(),
        n_original_gift_mainland_unique_species=len(gift),
        provider_source_counts=sorted(provider_info,key=lambda x:x["source_corpus"]),
        provider_consensus_source_only_flower_trait_opportunity=joined,
        claims=dict(
            distinct_corpora_not_independent_botanical_field_studies=True,
            source_provider_filename_is_known_original_public_repository_file=True,
            all_4505_original_island_plant_names_excluded_from_external_GIFT_donors=True,
            no_NHM_categorical_source_disagreement_forced_to_single_colour=True,
            primary_colour_vs_organ_sensitivity_separate=True,
            original_all_4505_macro_universe_retained=True,
            not_mainland_local_population_trait_or_pollination_success=True,
            no_original_taxon_island_GIFT_entity_or_geographic_identifiers_exported=True,
        ))


def main():
    p=argparse.ArgumentParser()
    for name in ("original-root","nhm-directory","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    source=a.original_root
    original_traits=source/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    original_flora=source/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    original_gift=source/"fixed/source/source/gift_native_mainland_flora.csv.gz"
    original_sha={
        original_traits:"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        original_flora:"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        original_gift:"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
    }
    import hashlib
    def sha(path):
        h=hashlib.sha256()
        with path.open("rb") as fp:
            for chunk in iter(lambda:fp.read(1<<20),b""):h.update(chunk)
        return h.hexdigest()
    for path,expected in original_sha.items():
        if sha(path)!=expected:
            raise ValueError("Original protected biological dataset SHA drift")
    providers={}
    for label,filename in PROVIDERS.items():
        path=a.nhm_directory/filename
        if not path.is_file() or nhm_git_blob_sha(path)!=BLOB_PINS[label]:
            raise ValueError("Original NHM provider blob identity not matched: "+label)
        providers[label]=read_provider_file(path,label)
    traits=strict_traits(pd.read_csv(original_traits,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original_flora,dtype=str,usecols=["island_id","accepted_species"])
    gift=pd.read_csv(original_gift,dtype=str,usecols=["entity_ID","work_species"])
    out=analyze(providers,traits,flora,gift)
    out["source_metadata"]=dict(
        NHM_original_repository="https://github.com/NaturalHistoryMuseum/plant-traits",
        NHM_pinned_git_commit=NHM_COMMIT,
        NHM_license=NHM_LICENSE,
        NHM_raw_source_git_blobs=BLOB_PINS,
    )
    payload=json.dumps(out,sort_keys=True,indent=2,allow_nan=False)
    for s in ('"island_id"','"accepted_species"','"work_species"',
              '"entity_ID"','"spatial_block"','"source_signature"'):
        if s in payload:raise ValueError("Private protected botanical identity in output")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(payload+"\n",encoding="utf-8")
    print(json.dumps(dict(status=out["status"],
                          n_provider_corpora=len(PROVIDERS),
                          original_global_islands=4505,
                          no_adaptive_causal_flower_claim=True)))


if __name__=="__main__":main()
