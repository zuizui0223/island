"""Independent publicly released NHM botanical-description flower traits × original GIFT islands.

NHM/plant-traits (CC BY 4.0) is an independently RELEASED source description
catalogue, NOT a measured mainland plant population or proof of pollinator
effectiveness. This original rights-protected runner exports ONLY anonymous
aggregate counts and conditional agreement, and NEVER private source names.
No required triple complete cases, no interpolated missing flower categories.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_source_name_heldout_macro_preflight import (
    short_distance_block_disjoint_endpoint_union,
)

NHM_COMMIT="4843c22bb2fc91018d9828dcddebf067c83b1097"
NHM_BLOB="e1cfd2d5c518c2061ec91d1a9e8b64fc7ca69389"
NHM_PATH="output/raw/ocr-traits/combined.csv"
NHM_LICENSE="CC BY 4.0"
EXTERNAL_COLUMNS=(
    "taxon","flower colour","flower shape","flower symmetry",
    "corolla colour","petal colour","source","wfo",
)
MAIN_COLOR="flower colour"
OPTIONAL_COLOR=("corolla colour","petal colour")
COLOR_TOKENS={
    "white":r"white|whitish|cream|ivory|snowy",
    "red_pink":r"red|reddish|rose|rosy|pink|pinkish|scarlet|crimson|magenta|ruby|carmine",
    "yellow_orange":r"yellow|yellowish|orange|golden|gold",
    "blue_purple":r"blue|bluish|purple|purplish|violet|violetish|lilac|lavender|mauve",
    "green_brown_inconspicuous":r"green|greenish|brown|brownish|dull|inconspicuous",
}
FORM_MAP={
    "campanulate":"bell_campanulate",
    "bell-shaped":"bell_campanulate",
    "bell shaped":"bell_campanulate",
    "tubular":"tubular",
    "tube":"tubular",
    "rotate":"open_radial",
    "funnel":"funnel_trumpet",
    "funnelform":"funnel_trumpet",
    "funnel-shaped":"funnel_trumpet",
    "salverform":"salverform",
    "ligulate":"ligulate",
}
SYMMETRY={
    "actinomorphic":"actinomorphic",
    "zygomorphic":"zygomorphic",
    "asymmetric":"asymmetric",
}
EXACT_BINOMIAL=re.compile(r"^[A-Z][a-z-]+ [a-z][a-z-]+$")
COLOR_RE={category:re.compile(r"\b(?:"+pattern+r")\b",re.I)
          for category,pattern in COLOR_TOKENS.items()}


def accepted_exact_binomial(s:str)->str|None:
    if not isinstance(s,str):return None
    value=" ".join(s.strip().split())
    return value if EXACT_BINOMIAL.fullmatch(value) else None


def single_description_color(text:str)->str|None:
    """Refuse mixed/unknown flower colours, unlike substring-first coding."""
    if not isinstance(text,str) or not text.strip():return None
    cleaned=text.lower().strip()
    if any(x in cleaned for x in ("variable","various","many-coloured",
                                   "multicol","unknown","not recorded")):
        return None
    hits=[k for k,regex in COLOR_RE.items() if regex.search(cleaned)]
    if len(hits)!=1:return None
    # Short, descriptive phrases only: reject numbers/long OCR prose.
    if len(cleaned)>85 or re.search(r"\d",cleaned):return None
    return hits[0]


def single_description_form(text:str)->str|None:
    if not isinstance(text,str):return None
    normalized=" ".join(text.lower().replace("_"," ").strip().split())
    return FORM_MAP.get(normalized)


def single_description_symmetry(text:str)->str|None:
    if not isinstance(text,str):return None
    return SYMMETRY.get(text.lower().strip())


def nhm_git_blob_sha(path:Path)->str:
    """Pin the actual raw UTF-8 NHM CSV git BLOB, not a moving main branch."""
    with path.open("rb") as f:b=f.read()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()


def nhm_external_categories(external:pd.DataFrame)->tuple[dict,dict]:
    if not set(EXTERNAL_COLUMNS)<=set(external):
        raise ValueError("NHM original source schema changed")
    frame=external[list(EXTERNAL_COLUMNS)].fillna("").astype(str)
    data=defaultdict(lambda:defaultdict(list))
    names=set()
    for row in frame.itertuples(index=False,name=None):
        value=dict(zip(EXTERNAL_COLUMNS,row))
        name=accepted_exact_binomial(value["taxon"])
        if name is None:continue
        names.add(name)
        primary=single_description_color(value[MAIN_COLOR])
        colours=[single_description_color(value[k]) for k in OPTIONAL_COLOR]
        colours=[v for v in colours if v]
        form=single_description_form(value["flower shape"])
        sym=single_description_symmetry(value["flower symmetry"])
        if primary:data[name]["primary"].append(primary)
        # Only include organ-colours if unambiguous and nonconflicting with
        # directly described whole-flower colour; reported separately.
        if colours and len(set(colours))==1:
            data[name]["organ"].append(colours[0])
        if form:data[name]["form"].append(form)
        if sym:data[name]["symmetry"].append(sym)
    out={x:{} for x in ("primary","expanded","form","symmetry")}
    multi={x:0 for x in out}
    for name,record in data.items():
        for axis in ("primary","form","symmetry"):
            vals=set(record.get(axis,[]))
            if len(vals)==1:out[axis][name]=next(iter(vals))
            elif len(vals)>1:multi[axis]+=1
        p=set(record.get("primary",[]))
        organ=set(record.get("organ",[]))
        combined=p|organ
        if len(combined)==1:out["expanded"][name]=next(iter(combined))
        elif len(combined)>1:multi["expanded"]+=1
    summary={
        "n_external_rows":len(frame),
        "n_exact_taxon_binomial_names":len(names),
        "n_unique_exact_binomial_species_by_axis":{
            axis:len(v) for axis,v in out.items()},
        "n_species_with_documentation_conflict_by_axis":multi,
        "source_documentation_rows_by_name":frame.source.value_counts().head(12).to_dict(),
    }
    return out,summary


def _coverage(original:dict, external:dict)->dict:
    overlap=set(original)&set(external)
    equal=sum(original[k]==external[k] for k in overlap)
    return dict(
        n_direct_names=len(original),
        n_external_names=len(external),
        n_exact_accepted_name_comparable=len(overlap),
        n_equal_categories=equal,
        n_different_categories=len(overlap)-equal,
        fraction_exact_name_agreement=(equal/len(overlap) if overlap else None),
    )


def _external_gift_source_coverage(all_flora:dict,
                                   candidates:list[str],
                                   measured:dict,
                                   source_by_entity:dict,
                                   nearest:dict)->dict:
    total=0
    nearest_genus=0
    n_eligible_floras=0
    n_floras_any_genus=0
    n_floras_source_missing=0
    n_islands_with_local_external_species=0
    source_local_names=set(measured)
    for island in candidates:
        target=all_flora.get(island,set())&set(measured.get("__target_original",{}))
        # Control target names passed through "__target_original" only.
        if not target:continue
        n_eligible_floras+=1
        local=set()
        entities=nearest.get(island,())
        if not entities:n_floras_source_missing+=1
        for e in entities:
            local.update(source_by_entity.get(e,set()))
        local&=source_local_names
        if local:n_islands_with_local_external_species+=1
        local_genera={x.split()[0] for x in local}
        supported=sum(sp.split()[0] in local_genera for sp in target)
        total+=len(target)
        nearest_genus+=supported
        if supported:n_floras_any_genus+=1
    return dict(
        n_target_original_direct_eligible_islands=n_eligible_floras,
        n_direct_original_target_carrier_occurrences=total,
        n_original_floras_with_any_external_colour_source_species=n_islands_with_local_external_species,
        n_islands_with_any_exact_congeneric_external_source_donor=n_floras_any_genus,
        n_original_floras_missing_any_near_GIFT_source=n_floras_source_missing,
        n_original_target_carriers_with_near_source_external_same_genus=nearest_genus,
        fraction_original_target_carriers_with_near_source_external_same_genus=(
            nearest_genus/total if total else None),
    )


def analyze(original_traits:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            mainland:pd.DataFrame,assignments:pd.DataFrame,
            external:pd.DataFrame):
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Missing original GIFT island flora")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicated original island×species")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original botanical status invalid")
    if not {"entity_ID","work_species"}<=set(mainland) or (
            mainland.duplicated(["entity_ID","work_species"]).any()):
        raise ValueError("Native mainland GIFT source identity invalid")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assignments):
        raise ValueError("Original rank-mode geographic source assignments invalid")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("4,505 original botanical flora required")
    all_plant_names=set().union(*full.values())
    external_states,external_summary=nhm_external_categories(external)
    direct_col=original_traits["flower_primary_color"]
    direct_shape=original_traits["floral_form"]
    direct_sym=original_traits["floral_symmetry"]
    labelled=set(direct_col)&set(direct_shape)
    exact_gift=set(mainland.work_species)
    donor_eligibility={}
    for axis,labels in external_states.items():
        all_external=set(labels)
        candidate=all_external&exact_gift
        outside=all_external-all_plant_names
        eligible=candidate-all_plant_names
        donor_eligibility[axis]=dict(
            n_external_directly_described_names=len(all_external),
            n_external_names_in_original_GIFT_native_mainland_lists=len(candidate),
            n_external_names_outside_ALL_4505_island_flora_names=len(outside),
            n_EXTERNAL_MEASURED_SOURCE_NAMES_in_GIFT_not_in_ANY_4505_island=len(eligible),
            n_external_mainland_GIFT_genera_outside_all_islands=len({
                name.split()[0] for name in eligible}),
        )
    crosscheck={
        "original_colour_vs_independent_NHM_primary":_coverage(direct_col,external_states["primary"]),
        "original_colour_vs_independent_NHM_organ_expanded":_coverage(direct_col,external_states["expanded"]),
        "original_form_vs_independent_NHM_shape":_coverage(direct_shape,external_states["form"]),
        "original_symmetry_vs_independent_NHM_symmetry":_coverage(direct_sym,external_states["symmetry"]),
    }
    # No future island test names in independent mainland source donor pool.
    source_names=set(external_states["primary"])&exact_gift-all_plant_names
    ranked=assignments.loc[assignments.source_mode.eq("geo_k5")].copy()
    ranked["source_rank"]=pd.to_numeric(ranked.source_rank,errors="raise")
    ranked=ranked.loc[ranked.source_rank.between(1,5)]
    if ranked.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT near-source IDs duplicated")
    near=ranked.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    relevant=mainland.loc[mainland.work_species.isin(source_names)]
    by_entity=relevant.groupby("entity_ID").work_species.agg(set).to_dict()
    geo=geo.set_index("island_id")
    n5=sorted(i for i in full if len(full[i]&set(direct_col))>=5)
    northmid_10=sorted(i for i in full
        if i in geo.index and geo.at[i,"analysis_regime"]=="northern_midlatitude"
        and len(full[i]&labelled)>=10 and
        len(full[i]&labelled)/len(full[i])>=.10)
    focal215,pairs,blocks=short_distance_block_disjoint_endpoint_union(
        northmid_10,full,geo.reset_index())
    if (len(focal215),pairs,blocks)!=(215,456,29):
        raise ValueError("Original baseline 215/456/29 selected pair drift")
    # Inline source mapping, not reported identities.
    measured={"__target_original":direct_col,**{
        key:value for key,value in external_states["primary"].items()
        if key in source_names}}
    pairs_of_interest={
        "all_direct_flower_colour_n5_original_floras":n5,
        "northern_midlatitude_456_pair_endpoint_215_floras":focal215,
    }
    source_coverage={name:_external_gift_source_coverage(full,ids,measured,by_entity,near)
        for name,ids in pairs_of_interest.items()}
    return dict(
        status="FD_EXTERNAL_NHM_BOTANICAL_DESCRIPTION_SOURCE_PREFLIGHT",
        n_original_botanical_islands=len(full),
        n_original_distinct_island_flora_species=len(all_plant_names),
        n_original_highmedium_colour_direct_species=len(direct_col),
        n_original_highmedium_form_direct_species=len(direct_shape),
        n_original_colour_x_form_direct_species=len(labelled),
        n_GIFT_native_mainland_source_rows=len(mainland),
        n_GIFT_native_mainland_source_regions=mainland.entity_ID.nunique(),
        n_GIFT_native_mainland_distinct_species=len(exact_gift),
        n_eligible_original_flower_colour_n5_islands=len(n5),
        n_original_posthoc_456_source_candidate_islands=len(focal215),
        external_NHM_source=external_summary,
        external_source_direct_trait_opportunity=donor_eligibility,
        independent_text_vs_original_direct_trait_agreement=crosscheck,
        external_source_only_mainland_colour_congeneric_opportunity=source_coverage,
        claims=dict(
            external_NHM_CCBY4_descriptive_text_extraction_not_population_field_measurement=True,
            external_source_names_absent_from_ALL_4505_islands_for_independent_donors=True,
            no_HNM_nlp_mixed_colour_forced_to_single_trait=True,
            no_new_observed_field_pollination_or_ancestry=True,
            original_all_4505_not_replaced_with_215_sensitivity=True,
            no_independent_morphology_source_claim_without_agreement_and_coverage=True,
            no_original_plant_island_source_id_exported=True,
        ))


def sha256(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    for name in ("original-root","geology","nhm-combined","output"):
        ap.add_argument("--"+name,type=Path,required=True)
    args=ap.parse_args()
    original={
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assignment":args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
        "mainland":args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
    }
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "assignment":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
        "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
    }
    for name,path in original.items():
        if not path.is_file() or sha256(path)!=pins[name]:
            raise ValueError("Original archived GIFT/Direct SHA pin changed: "+name)
    blob=nhm_git_blob_sha(args.nhm_combined)
    if blob!=NHM_BLOB:raise ValueError("Independent NHM CSV Git blob SHA drift")
    frame=pd.read_csv(args.nhm_combined,dtype=str,keep_default_na=False,
                      low_memory=False,usecols=list(EXTERNAL_COLUMNS))
    if not 50000<=len(frame)<=120000:
        raise ValueError("Unexpected NHM combined source row universe")
    traits=strict_traits(pd.read_csv(original["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "island_latitude","island_longitude"])
    assignments=pd.read_csv(original["assignment"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    native=pd.read_csv(original["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    output=analyze(traits,flora,geo,native,assignments,frame)
    if (output["n_original_botanical_islands"]!=4505 or
        output["n_GIFT_native_mainland_source_regions"]!=748 or
        output["n_eligible_original_flower_colour_n5_islands"]!=3132):
        raise ValueError("Original full-population Direct botanical baseline changed")
    output["source_metadata"]=dict(
        NHM_public_repository="https://github.com/NaturalHistoryMuseum/plant-traits",
        NHM_source_commit=NHM_COMMIT,
        NHM_combined_source_path=NHM_PATH,
        NHM_combined_git_blob_sha=blob,
        NHM_license=NHM_LICENSE,
        NHM_original_csv_sha256=sha256(args.nhm_combined),
        original_GIFT_archive_source_sha256=pins,
    )
    payload=json.dumps(output,sort_keys=True,indent=2,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"work_species"',
                    '"entity_ID"','"spatial_block"','"island_latitude"',
                    '"island_longitude"','"source_signature"'):
        if private in payload:
            raise ValueError("Original botanical/GIFT private identity exported")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(payload+"\n",encoding="utf-8")
    print(json.dumps(dict(status=output["status"],
                          n_original_island_floras=4505,
                          external_flower_colour_source_species=output[
                            "external_NHM_source"]["n_unique_exact_binomial_species_by_axis"]["primary"],
                          no_causal_pollinator_or_adaptation_inference=True)))


if __name__=="__main__":main()
