"""Source-pinned GloBI v0.9 provider/citation diagnostic for fd.

Works only on SUPPORTS interactions.tsv.gz, never refuted interactions. Source
provider IDs/citations/taxon names stay in memory. Emit aggregate receipt only.
Global documented reference composition != local pollination effectiveness.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import gzip
import hashlib
from pathlib import Path
import pandas as pd
import numpy as np

PIN="095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17"
VERSION="10.5281/zenodo.20546682"
RO={
 "http://purl.obolibrary.org/obo/RO_0002455":("source","target","pollination"),
 "http://purl.obolibrary.org/obo/RO_0002456":("target","source","pollination"),
 "http://purl.obolibrary.org/obo/RO_0002622":("source","target","flower_visit"),
 "http://purl.obolibrary.org/obo/RO_0002623":("target","source","flower_visit"),
}
CHANNELS=("bombus","non_bombus_bees","lepidoptera","flower_visiting_birds","diptera")
BEE_FAMILIES=frozenset(("Andrenidae","Apidae","Colletidae","Halictidae","Megachilidae","Melittidae","Stenotritidae"))
RAW=(
 "sourceTaxonGenusName","sourceTaxonSpeciesName","sourceTaxonFamilyName",
 "sourceTaxonOrderName","sourceTaxonClassName","sourceTaxonKingdomName",
 "targetTaxonGenusName","targetTaxonSpeciesName","targetTaxonFamilyName",
 "targetTaxonOrderName","targetTaxonClassName","targetTaxonKingdomName",
 "interactionTypeId","referenceDoi","referenceCitation","referenceUrl",
 "sourceDOI","sourceArchiveURI","sourceNamespace","sourceCitation",
)
COLS=("plant_genus","visitor_channel","provider_id","reference_id","evidence_class",
      "source_namespace","source_doi","source_citation")

def val(x):
    if x is None or pd.isna(x):return ""
    return " ".join(str(x).strip().split())

def reference(r):
    doi=val(r.get("referenceDoi"))
    if doi:
        doi=doi.lower().removeprefix("https://doi.org/").removeprefix("http://doi.org/").removeprefix("doi:")
        return "doi:"+doi
    if val(r.get("referenceCitation")):return "citation:"+val(r["referenceCitation"]).lower()
    if val(r.get("referenceUrl")):return "url:"+val(r["referenceUrl"])
    # Critical: never count a source dataset identifier as a fake study.
    return ""

def provider(r):
    if val(r.get("sourceArchiveURI")):return "archive:"+val(r["sourceArchiveURI"])
    if val(r.get("sourceDOI")):return "doi:"+val(r["sourceDOI"]).lower()
    if val(r.get("sourceNamespace")) and val(r.get("sourceCitation")):
        return "namespace_citation:"+val(r["sourceNamespace"])+"|"+val(r["sourceCitation"])
    return ""

def channel(r, prefix):
    genus=val(r.get(prefix+"TaxonGenusName"))
    family=val(r.get(prefix+"TaxonFamilyName"))
    order=val(r.get(prefix+"TaxonOrderName"))
    cls=val(r.get(prefix+"TaxonClassName"))
    if genus=="Bombus":return "bombus"
    if family in BEE_FAMILIES and genus!="Bombus":return "non_bombus_bees"
    if order=="Lepidoptera":return "lepidoptera"
    if cls=="Aves":return "flower_visiting_birds"
    if order=="Diptera":return "diptera"
    return ""

def canonicalize_chunk(frame:pd.DataFrame,focus:set[str]|None=None):
    missing=set(RAW)-set(frame)
    if missing:raise ValueError("Missing original GloBI source fields: "+str(sorted(missing)))
    subset=frame[frame.interactionTypeId.fillna("").astype(str).isin(RO)].copy()
    ctr=Counter({"raw_input_rows":len(frame),"flower_relation_rows":len(subset)})
    rows=[]
    for r in subset.fillna("").to_dict("records"):
        visitor,plant,strength=RO[val(r["interactionTypeId"])]
        if val(r.get(plant+"TaxonKingdomName"))!="Plantae":
            ctr["nonplant_or_unknown_kingdom"]+=1;continue
        genus=val(r.get(plant+"TaxonGenusName"))
        if not genus:
            ctr["unresolved_plant_genus"]+=1;continue
        if focus is not None and genus not in focus:
            ctr["out_of_focus_genus"]+=1;continue
        if not val(r.get(visitor+"TaxonSpeciesName")):
            ctr["unresolved_visitor_species"]+=1;continue
        ch=channel(r,visitor)
        if not ch:
            ctr["unresolved_visitor_channel"]+=1;continue
        prov=provider(r)
        if not prov:
            ctr["missing_source_provider"]+=1;continue
        ref=reference(r)
        if not ref:
            ctr["missing_bibliographic_reference"]+=1;continue
        rows.append((genus,ch,prov,ref,strength,val(r.get("sourceNamespace")),
                     val(r.get("sourceDOI")),val(r.get("sourceCitation"))))
    ctr["cited_flower_interaction_rows"]=len(rows)
    return pd.DataFrame.from_records(rows,columns=COLS).drop_duplicates(),dict(ctr)

class ProviderReferenceEvidence:
    def __init__(self):
        self.rows=set()
        self.tallies=Counter()
        # Structural source-provider hierarchy audit: sourceArchiveURI may group
        # many independent sourceNamespace, DOI or sourceCitation records.
        # These IDs never leave the runner; only grouped cardinalities do.
        self.provider_namespaces=defaultdict(set)
        self.provider_dois=defaultdict(set)
        self.provider_citations=defaultdict(set)

    def add(self,canonical:pd.DataFrame,tally:dict|None=None):
        if set(COLS)-set(canonical):raise ValueError("Missing canonical reference fields")
        if tally:self.tallies.update(tally)
        for r in canonical.itertuples(index=False):
            self.rows.add((r.plant_genus,r.visitor_channel,r.provider_id,r.reference_id))
            if r.source_namespace:
                self.provider_namespaces[r.provider_id].add(r.source_namespace)
            if r.source_doi:
                self.provider_dois[r.provider_id].add(r.source_doi)
            if r.source_citation:
                self.provider_citations[r.provider_id].add(r.source_citation)

    def summarize(self,min_refs=3):
        if min_refs<3:raise ValueError("Expected >=3 distinct references")
        if not self.rows:raise ValueError("No GloBI citation-backed study units")
        by_reference=defaultdict(set)
        for _,_,prov,ref in self.rows:by_reference[ref].add(prov)
        collision={ref for ref,providers in by_reference.items() if len(providers)>1}
        unique=[r for r in self.rows if r[3] not in collision]
        if not unique:raise ValueError("No uncollided citation-backed units")
        df=pd.DataFrame.from_records(unique,columns=("plant_genus","visitor_channel","provider_id","reference_id"))
        provider_n=df.provider_id.value_counts()
        excluded=provider_n.index[0]
        original=df.groupby(["plant_genus","visitor_channel"]).reference_id.nunique().unstack(fill_value=0).reindex(columns=CHANNELS,fill_value=0)
        after=df[df.provider_id.ne(excluded)].groupby(["plant_genus","visitor_channel"]).reference_id.nunique().unstack(fill_value=0).reindex(index=original.index,columns=CHANNELS,fill_value=0)
        orig_ok=original.sum(axis=1)>=min_refs
        after_ok=after.sum(axis=1)>=min_refs
        common=orig_ok&after_ok
        x=original.loc[common].div(original.loc[common].sum(axis=1),axis=0)
        y=after.loc[common].div(after.loc[common].sum(axis=1),axis=0)
        tv=(x-y).abs().sum(axis=1)/2
        receipt={
            "source_version_doi":VERSION,"source_sha256":PIN,
            "status":"PROVIDER_COVERAGE_ONLY_NOT_FAMILY_PREDICTION",
            "n_original_interaction_rows":int(self.tallies["raw_input_rows"]),
            "n_matched_citation_source_rows":int(self.tallies["cited_flower_interaction_rows"]),
            "n_matching_genera":int(len(original)),
            "n_distinct_source_providers":int(len(provider_n)),
            "n_crossprovider_repeated_references_excluded":int(len(collision)),
            "n_unique_plant_guild_provider_reference_units":int(len(df)),
            "top_provider_fraction_of_reference_units":float(provider_n.iloc[0]/len(df)),
            "largest_provider_distinct_subnamespace_count":int(len(self.provider_namespaces[excluded])),
            "largest_provider_distinct_source_doi_count":int(len(self.provider_dois[excluded])),
            "largest_provider_distinct_source_citation_count":int(len(self.provider_citations[excluded])),
            "provider_grouping_rule":"sourceArchiveURI_first_else_sourceDOI_else_namespace_plus_sourceCitation",
            "provider_aggregator_granularity_verified":False,
            "n_genera_with_3plus_refs_before":int(orig_ok.sum()),
            "n_genera_with_3plus_refs_after":int(after_ok.sum()),
            "n_common_genera_3plus_refs":int(common.sum()),
            "median_genus_documentation_total_variation_after_omit":float(tv.median()) if len(tv) else None,
            "provider_family_prediction_rerun_completed":False,
            "source_and_taxon_identifiers_released":False,
            "pollination_effectiveness_observed":False,
        }
        return original,after,receipt

def sha256_file(path:Path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()

def analyze_v09_gzip(path:Path,focus:set[str],chunksize=200_000,
                     extra_hook=None,extra_columns=()):
    """Single-source pass with optional PRIVATE exact-species evidence hook.

    A second reader must not silently use a newer snapshot. The original
    SHA256 and original SOURCE-based genus eligibility are verified first.
    """
    if not focus:raise ValueError("Explicit outcome-independent focal genus list required")
    if sha256_file(path)!=PIN:raise ValueError("GloBI v0.9 original SHA mismatch")
    acc=ProviderReferenceEvidence()
    # Original v0.9 may omit event-time and some geography columns that
    # appear in documentation for other GloBI versions. Only original
    # taxonomy/interaction/reference columns are mandatory. Never fabricate
    # an absent optional field or silently substitute another GloBI release.
    with gzip.open(path,"rt",encoding="utf-8") as source:
        available=set(source.readline().rstrip("\n\r").split("\t"))
    missing=set(RAW)-available
    if missing:
        raise ValueError("Pinned v0.9 source lacks required interaction fields: "+
                         str(sorted(missing)))
    present_optional=tuple(c for c in extra_columns if c in available)
    fields=tuple(dict.fromkeys((*RAW,*present_optional)))
    for chunk in pd.read_csv(path,sep="\t",compression="gzip",dtype=str,
                             usecols=fields,chunksize=chunksize,low_memory=False):
        if extra_hook is not None:
            extra_hook(chunk)
        row,tally=canonicalize_chunk(chunk,focus)
        acc.add(row,tally)
    return acc.summarize()
