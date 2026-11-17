"""Mangal original archived network node/positive-interaction island preflight.

Study uses immutable original 4,505 GIFT island source population and
927 actually source-verified GSHHG island polygons, but Mangal sites are
modern georeferenced metadata; source site GPS precision remains unverified.
All original plant/island, bibliography, exact site IDs stay private.
Mangal is a public network AGGREGATOR: overlapping original GloBI or Web of
Life studies have not been de-duplicated, so NEVER call source rows independent
effective pollination observations or guild causal effects.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
from collections import Counter

import pandas as pd
import requests

from fd_island_functional_beta import strict_traits
from fd.independent_pollination_networks import (
    mangal_location_metadata,source_bounded_site_plant_audit
)
from fd.independent_mangal_pollination_links import (
    strict_mangal_pollination_candidate_metadata,
    extract_mangal_positive_plant_links,
)
from fd.gbif_mangal_visitor_taxonomy import (
    target_original_morphology_partner_keys,
    summarize_gbif_resolved_source,
)
from run_fd_independent_network_site_preflight import (
    ORIGINAL_SHA,strict_source_geometry,same_original_island_network_points,
    get_all_mangal,get_bytes,MANGAL_BASE,
)

MAX_NETWORKS=120
MAX_PAGES_PER_NETWORK=4
MAX_TARGET_GBIF_TAXON_LOOKUPS=160


def fetch_scoped(session,section,network_id):
    result=[];source_sha=[]
    for page in range(MAX_PAGES_PER_NETWORK):
        url=f"{MANGAL_BASE}/{section}?network_id={network_id}&count=1000&page={page}"
        blob=get_bytes(session,url,max_bytes=18000000)
        source_sha.append(sha256(blob).hexdigest())
        document=json.loads(blob)
        if isinstance(document,list):
            subset=document
        elif isinstance(document,dict) and isinstance(document.get("data"),list):
            subset=document["data"]
        elif isinstance(document,dict) and isinstance(document.get("results"),list):
            subset=document["results"]
        else:
            raise ValueError("Mangal scoped original API response is not a record array")
        if any(str(z.get("network_id"))!=str(network_id) for z in subset):
            raise ValueError("Mangal did not respect requested original source network ID")
        result.extend(subset)
        if len(subset)<1000:
            return result,dict(
                source_status="COMPLETE_SCOPED_ORIGINAL_NETWORK",
                n_original_rows=len(result),
                source_page_sha256=source_sha,
                page_cap=MAX_PAGES_PER_NETWORK,
            )
    raise ValueError("Mangal original network exceeds closed pagination cap")


def run(original_root:Path,geometry:Path,receipt:Path)->dict:
    original_ledger=original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    original_flora=original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for p in (original_ledger,original_flora):
        if not p.is_file() or sha256(p.read_bytes()).hexdigest()!=ORIGINAL_SHA[p.name]:
            raise ValueError("original 4505 GIFT accepted flora or true Direct trait SHA mismatch")
    direct=pd.read_csv(original_ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality",
    ])
    flora=pd.read_csv(original_flora,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"
    ])
    polygon,full=strict_source_geometry(geometry,receipt,direct,flora)
    flora=flora.loc[flora.island_id.isin(set(polygon.island_id))].copy()
    flora_by_site={
        island:set(g.accepted_species)
        for island,g in flora.groupby("island_id",sort=True)
    }
    traits=strict_traits(direct)
    session=requests.Session()
    source_error=None
    try:
        datasets,ds_meta=get_all_mangal(session,"dataset")
        nets,net_meta=get_all_mangal(session,"network")
        metadata,base=mangal_location_metadata(nets,datasets)
        original_sites=same_original_island_network_points(metadata,polygon)
    except (requests.RequestException,KeyError,ValueError) as exc:
        datasets=[];nets=[];base={}
        ds_meta=net_meta={}
        original_sites=pd.DataFrame(columns=["network_id","island_id"])
        source_error=type(exc).__name__+":"+str(exc)[:160]

    network_by_id={str(x.get("id")):x for x in nets}
    dataset_by_id={str(x.get("id")):x for x in datasets}
    strict_candidates=[]
    skipped_text=0
    for item in original_sites.itertuples(index=False):
        nw=network_by_id.get(item.network_id)
        if nw is None:
            raise AssertionError("original Mangal original site has no network metadata")
        ds=dataset_by_id.get(str(nw.get("dataset_id")),{})
        if not strict_mangal_pollination_candidate_metadata(nw,ds):
            skipped_text+=1
            continue
        strict_candidates.append(item)
    # This is source-level study title/classification, NOT verified visit.
    status=Counter()
    positive=[]
    original_floral_target_partners=[]
    n_interactions=0
    n_source_node=0
    bytype=Counter()
    provider_studies=set()
    errors=[]
    for item in strict_candidates[:MAX_NETWORKS]:
        network_id=item.network_id
        try:
            nodes,nmeta=fetch_scoped(session,"node",network_id)
            relations,imeta=fetch_scoped(session,"interaction",network_id)
            documented,quality=extract_mangal_positive_plant_links(
                nodes,relations,flora_by_site.get(item.island_id,set()),
                network_id=network_id
            )
            status["n_source_original_networks_fetched"]+=1
            n_interactions+=len(relations)
            n_source_node+=len(nodes)
            status["n_original_networks_with_positive_exact_plant_link"]+=int(len(documented)>0)
            status["n_original_accepted_plant_links_positive"]+=len(documented)
            status["n_original_links_taxonomically_classified_guild"]+=quality[
                "n_partner_taxonomy_classified_to_family_or_order"]
            provider_studies.add(str(network_by_id[network_id].get("dataset_id")))
            bytype.update(documented.interaction_type.value_counts().to_dict())
            if len(documented):
                documented["island_id"]=item.island_id
                # Only the source-positive partner taxa of plants with
                # BOTH original Direct floral axes are GBIF-taxonomy targets.
                checked=target_original_morphology_partner_keys(
                    documented,traits
                )
                if len(checked):
                    original_floral_target_partners.append(checked)
                positive.append(documented[["network_id","island_id","plant","guild"]])
        except (requests.RequestException,ValueError,KeyError) as exc:
            status["n_original_networks_source_unavailable_or_unsupported"]+=1
            errors.append(type(exc).__name__+":"+str(exc)[:145])

    rows=pd.concat(positive,ignore_index=True) if positive else pd.DataFrame(
        columns=["network_id","island_id","plant","guild"])
    existing=source_bounded_site_plant_audit(
        rows,flora,traits,global_original_islands=4505
    )
    if original_floral_target_partners:
        private_morph=pd.concat(
            original_floral_target_partners,ignore_index=True
        )
    else:
        private_morph=pd.DataFrame(columns=[
            "network_id","island_id","plant","partner_gbif_key",
            "source_interaction_id","actinomorphic_deep_true",
            "zygomorphic_deep_true","tube_bell_deep_true",
        ])
    original_status=flora.rename(
        columns={"accepted_species":"plant"}
    )[["island_id","plant","origin_status","endemic_status"]]
    private_morph=private_morph.merge(
        original_status,on=["island_id","plant"],
        how="left",validate="many_to_one"
    )
    if len(private_morph) and private_morph.origin_status.isna().any():
        raise ValueError("Mangal positive floral visitor link lacks true original flora native status")
    all_keys=sorted(
        {int(x) for x in private_morph.partner_gbif_key.dropna()}
    )
    selected_keys=all_keys[:MAX_TARGET_GBIF_TAXON_LOOKUPS]
    gbif_taxonomy={}
    gbif_errors=Counter()
    for original_key in selected_keys:
        try:
            blob=get_bytes(
                session,f"https://api.gbif.org/v1/species/{original_key}",
                max_bytes=300000
            )
            datum=json.loads(blob)
            if (not isinstance(datum,dict) or
                int(datum.get("key",0))!=original_key):
                raise ValueError("GBIF taxon source key did not match original Mangal ID")
            gbif_taxonomy[original_key]=datum
        except (requests.RequestException,ValueError,KeyError) as exc:
            gbif_errors[type(exc).__name__]+=1
    gbif_audit=summarize_gbif_resolved_source(
        private_morph,gbif_taxonomy,total_original_islands=4505
    )
    gbif_audit["source_lookup_coverage"]={
        "n_original_true_flower_pair_links_without_source_gbif_key":int(
            private_morph.partner_gbif_key.isna().sum()
        ),
        "n_original_true_flower_pair_unique_gbif_keys":len(all_keys),
        "n_requested_taxon_id_lookups":len(selected_keys),
        "n_original_source_taxa_over_explicit_lookup_cap":max(
            0,len(all_keys)-MAX_TARGET_GBIF_TAXON_LOOKUPS
        ),
        "n_GBIF_taxon_lookups_failed":sum(gbif_errors.values()),
        "n_GBIF_lookup_error_classes":dict(gbif_errors),
        "n_valid_GBIF_taxon_lookups":len(gbif_taxonomy),
        "unmatched_or_missing_original_taxon_keys_not_negative_visitors":True,
    }
    return dict(
        status="ORIGINAL_GIFT_4505_INDEPENDENT_MANGAL_NETWORK_INTERACTION_PREFLIGHT",
        n_original_botanical_GIFT_island_population=4505,
        n_original_geometry_authenticated_island_subset=927,
        source_mangal_metadata_availability=(
            "SOURCE_RETRIEVED" if source_error is None else
            "SOURCE_UNAVAILABLE_DO_NOT_TREAT_AS_NEGATIVE"
        ),
        source_mangal_metadata_error_class=source_error,
        original_mangal_source_metadata=dict(
            **base,
            source_dataset_page_provenance=ds_meta,
            source_network_page_provenance=net_meta,
            n_original_site_keyword_candidates_matching_GIFT_polygons=int(
                len(original_sites)),
            n_original_original_islands_with_keyword_candidates=int(
                original_sites.island_id.nunique()),
            n_original_keyword_candidates_failing_strict_flower_descriptor=skipped_text,
            n_strict_original_publication_described_flower_networks=int(
                len(strict_candidates)),
            n_networks_skipped_explicit_fetch_cap=int(
                max(0,len(strict_candidates)-MAX_NETWORKS)),
        ),
        original_mangal_positive_network_evidence=dict(
            **dict(status),
            n_original_nodes_loaded_from_scoped_public_source=n_source_node,
            n_original_positive_or_zero_type_unfiltered_interactions_loaded=n_interactions,
            n_original_mangal_datasets_with_fetched_candidate_networks=len(provider_studies),
            n_original_positive_island_plant_guild_source_units_before_flora_dedup=len(rows),
            allowed_positive_source_interaction_type_counts=dict(bytype),
            network_errors_count=len(errors),
            network_errors_sample_class_only=errors[:6],
            source_partner_family_order_unknown_not_imputed=True,
            Mangal_source_dataset_not_independent_original_study=True,
        ),
        exact_island_plant_guild_and_real_flower_traits=existing,
        original_Mangal_flower_partner_GBIF_classification=gbif_audit,
        original_direct_trait_ledger_SHA256=sha256(original_ledger.read_bytes()).hexdigest(),
        original_GIFT_island_flora_SHA256=sha256(original_flora.read_bytes()).hexdigest(),
        epistemic_limits=dict(
            Mangal_source_keyword_site_not_proof_of_pollination=True,
            network_node_taxonomy_without_family_does_not_identify_guild=True,
            only_positive_original_interactions_counted=True,
            zero_source_edge_strength_never_pseudoabsent_visit=True,
            entire_4505_original_GIFT_flora_frame_retained=True,
            original_927_georeferenced_subset_not_the_full_botanical_population=True,
            original_source_networks_may_be_in_GloBI_or_Web_of_Life=True,
            no_measured_effective_pollination_seed_set_or_historical_shift=True,
            source_api_snapshot_not_immutable_original_publication_data=True,
        )
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--original-gshhg-geometry",type=Path,required=True)
    p.add_argument("--original-geometry-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    out=run(a.original_root,a.original_gshhg_geometry,
            a.original_geometry_receipt)
    text=json.dumps(out,indent=2,ensure_ascii=False,
                    sort_keys=True,allow_nan=False)
    for key in ("island_id","accepted_species","reference_id","provider_id",
                "latitude","longitude","entity_ID","work_species",
                "source_signature","spatial_block","source_genus",
                "source_family","visitor_species"):
        if '"'+key+'"' in text:
            raise ValueError("protected original botanical GIFT source identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(text+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=out["status"],
        n_original_GIFT_islands=4505,
        n_mangal_strict_source_networks=out["original_mangal_source_metadata"][
            "n_strict_original_publication_described_flower_networks"],
        n_original_mangal_positive_original_island_plant_network_links=out[
            "exact_island_plant_guild_and_real_flower_traits"][
                "n_original_species_network_plant_guild_units"],
        n_original_real_flower_pair_visitors_with_GBIF_taxon_keys=out[
            "original_Mangal_flower_partner_GBIF_classification"][
                "source_lookup_coverage"]["n_original_true_flower_pair_unique_gbif_keys"],
        no_unobserved_visitor_absence_or_shift_inferred=True,
    ),sort_keys=True))


if __name__=="__main__":
    main()
