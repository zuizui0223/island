"""Original GIFT 4505 × SHA-fixed 2026 GloBI positive documentary exact-site audit.

Unlike the original 928 complete-colour/form/symmetry source pilot, this
uses EVERY genuine SAME-SPECIES high/medium Direct TWO-AXIS flower mask
from the full original 4505 GIFT flora while independently REPLAYING the
historical 928 / 927 / 552-source-plant baseline. Unknown animal guild
or unobserved source link is NEVER ecological absence or measured seed set.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import argparse
import json
import pandas as pd
import geopandas as gpd

from fd_island_functional_beta import strict_traits
from island_v2.fd_globi_species_shift_v09 import (
    ExactSpeciesGloBICoverage, GEO_FIELDS,
)
from island_v2.fd_globi_provider_v09 import (
    RAW, PIN, sha256_file,
)
from island_v2.fd_globi_macro_928_locality_v09 import (
    strict_direct_928_flora,macro_exact_928_locality_gate,
)
from fd.fullframe_wol_network_bridge import full_botanical_exact_geometry_gate
from fd.fullframe_globi_local_support import (
    genuine_direct_motifs,strict_source_geography,
    original_source_guild_morphology_availability,
)

PLANTS_SHA="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"
DIRECT_SHA="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"

def _largest_original_documented_archive(source):
    excluded={ref for ref,providers in source.providers_by_reference.items()
              if len(providers)>1}
    valid=[row for row in source.units if row[2] not in excluded]
    counter=Counter(row[3] for row in valid)
    if not counter:
        raise ValueError("Original source has no unambiguous bibliography/provider")
    return sorted(counter,key=lambda k:(-counter[k],k))[0],excluded

def original_full_macro(root:Path,globi:Path,
                        polygon_path:Path,receipt_path:Path,
                        geography_path:Path,*,chunksize=200000):
    ledger=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    plants=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    if sha256_file(ledger)!=DIRECT_SHA or sha256_file(plants)!=PLANTS_SHA:
        raise ValueError("Full GIFT 4505 source species and Direct traits have changed")
    if sha256_file(globi)!=PIN:
        raise ValueError("Immutable GloBI v0.9 original relationship stream SHA changed")
    direct=pd.read_csv(ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(plants,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    if flora.island_id.nunique()!=4505:
        raise ValueError("4505 source botanical islands not preserved")
    previous,_=strict_direct_928_flora(direct,flora)
    old_names=set(previous.accepted_species)
    if previous.island_id.nunique()!=928:
        raise ValueError("Historical 928 morphology-rich island pilot changed")
    states=strict_traits(direct)
    original_names=set(flora.accepted_species)
    motifs=genuine_direct_motifs(states)
    genuine_two_axis_names=(set().union(*(set(mapping) for mapping in motifs.values()))
                            & original_names)
    # The historical 928 source pilot uses its own frozen strict-three
    # documentation ontology. Query all those names solely to reproduce its
    # source baseline, WITHOUT assigning any missing new two-axis phenotype.
    source_names=genuine_two_axis_names | old_names
    n_historical_only_names=len(old_names-genuine_two_axis_names)
    # Original full GIFT source 4505 island and incomplete accepted flora
    # records remain in population, although the GloBI species query has a
    # strict MORPHOLOGY-PARTIALLY-MEASURED source-name opportunity.
    full_source=ExactSpeciesGloBICoverage(source_names)
    old_source=ExactSpeciesGloBICoverage(old_names)
    cols=set(RAW)|set(GEO_FIELDS)
    chunks=0
    for frame in pd.read_csv(globi,sep="\t",compression="gzip",
         dtype=str,usecols=lambda name:name in cols,
         low_memory=False,chunksize=chunksize):
        full_source.add_chunk(frame)
        old_source.add_chunk(frame)
        chunks+=1
    if (full_source.raw["n_original_GloBI_rows_read"]!=24577183
        or old_source.raw["n_original_GloBI_rows_read"]!=24577183):
        raise ValueError("Complete frozen 24.6M GloBI source original row count mismatch")

    receipt=json.loads(receipt_path.read_text())
    if sha256_file(polygon_path)!=receipt.get("geometry_file_SHA256"):
        raise ValueError("Original GSHHG hash-verified polygons altered")
    gpkg=gpd.read_file(polygon_path,layer="islands")
    if gpkg.crs is None or gpkg.crs.to_epsg()!=4326:
        raise ValueError("GSHHG source CRS drift")
    polygons,geometry_audit=full_botanical_exact_geometry_gate(
        gpkg,flora,set(previous.island_id),receipt)
    # The old original-source morphology-conditional exact locality output
    # is independently replayed on the same 927 verified original polygons.
    old_largest,old_bad=_largest_original_documented_archive(old_source)
    original_928_audit=macro_exact_928_locality_gate(
        old_source.georeferenced_documented_relations,
        old_source.providers_by_reference,previous,
        polygon_path,receipt_path,
        largest_global_archive_first_label=old_largest,
    )
    if (original_928_audit["n_islands_matching_original_GSHHG_geometric_id"]!=927
        or original_928_audit["n_distinct_islands_with_species_locality_match"]!=109
        or original_928_audit["n_distinct_island_species_pairs_with_record"]!=552
        or original_928_audit[
            "n_strict_native_flora_island_species_pairs_with_visitor_relation"]!=157):
        raise ValueError("Historical 928 GloBI exact source bibliography site baseline drift")

    biggest,bad_refs=_largest_original_documented_archive(full_source)
    raw_relations=[row for row in full_source.georeferenced_documented_relations
                   if row[2] not in bad_refs]
    linked,link_audit=strict_source_geography(
        raw_relations,polygons[["island_id","geometry"]],flora)
    regions=pd.read_csv(geography_path,dtype=str,
        usecols=["island_id","analysis_regime"])
    if regions.island_id.duplicated().any():
        raise ValueError("Published source GIFT geographical regime duplicated")
    grouped=original_source_guild_morphology_availability(
        linked,flora,states,biggest,
        regions.set_index("island_id").analysis_regime.to_dict())
    if grouped["n_source_positive_original_island_plant_pairs"]<552:
        raise AssertionError("Full GIFT true Direct two-axis subset lost the historical original cited island×plant pairs")

    return dict(
        status="ORIGINAL_FULL_4505_GIFT_TRUE_PAIR_GLOBI_CITED_LOCAL_GUILD_OPPORTUNITY",
        n_original_parent_GIFT_botanical_island_floras=4505,
        n_original_source_highmedium_direct_two_axis_candidate_species=len(genuine_two_axis_names),
        n_historical_replay_only_names_not_used_as_new_direct_motifs=n_historical_only_names,
        n_original_source_query_names_including_replay=len(source_names),
        n_original_historical_triple_direct_source_names=len(old_names),
        original_source_GIFT_flora_SHA256=sha256_file(plants),
        original_source_Direct_traits_SHA256=sha256_file(ledger),
        original_source_GloBI_v09_SHA256=sha256_file(globi),
        n_original_GloBI_rows_checked=int(full_source.raw["n_original_GloBI_rows_read"]),
        n_original_source_chunks=chunks,
        geometry_original_biological_source=geometry_audit,
        old_928_locality_species_source_replay=dict(
            n_old_original_verified_direct_islands=928,
            n_old_exact_polygons=927,
            n_old_original_islands_with_matching_local_documentation=109,
            n_old_original_cited_plant_island_pairs=552,
            n_old_original_documented_literal_native_plant_island_pairs=157,
            replay_passed=True,
        ),
        n_global_excluded_bibliography_labels_multi_archive=len(bad_refs),
        original_source_full_GloBI_reference_level=dict(
            n_original_exact_species_predicate_rows=int(
                full_source.raw["n_exact_species_predicate_rows"]),
            n_original_georeferenced_candidate_relation_records=len(
                full_source.georeferenced_documented_relations),
            n_original_distinct_eligible_species_source_name_labels=len(full_source.names),
            n_original_source_archive_first_labels=len({
                row[3] for row in full_source.units}),
        ),
        geographical_original_GloBI_documentary_join=link_audit,
        original_site_plant_real_morphological_guild_support=grouped,
        biodiversity_and_pollinator_limitation=dict(
            original_global_4505_parent_retained=True,
            GloBI_citation_coordinate_not_direct_field_survey_verified=True,
            same_species_direct_morphology_not_measured_at_local_visitor_record=True,
            original_GIFT_source_plant_native_status_not_reconstructed=True,
            original_GloBI_largest_source_omission_not_independent_literature_replication=True,
            no_visitor_nonobservations_taken_as_animal_absence=True,
            no_sphingidae_versus_butterfly_split_from_coarse_GloBI_channel=True,
            no_pollination_efficiency_seedset_or_historical_guild_shift=True,
        ),
    )

def main():
    ap=argparse.ArgumentParser()
    for key in ("original-root","globi-original-tsv-gz",
                "original-gshhg-geometry","original-geometry-receipt",
                "original-geography","output"):
        ap.add_argument("--"+key,type=Path,required=True)
    args=ap.parse_args()
    result=original_full_macro(
        args.original_root,args.globi_original_tsv_gz,
        args.original_gshhg_geometry,args.original_geometry_receipt,
        args.original_geography)
    body=json.dumps(result,ensure_ascii=False,sort_keys=True,
                    indent=2,allow_nan=False)
    for term in ("accepted_species","island_id","reference_id",
                 "provider_id","latitude","longitude","visitor_species",
                 "network_id","source_signature","source_genus",
                 "source_region_id","original_site_id"):
        if '"'+term+'"' in body:
            raise ValueError("Original exact island, accepted taxon or citation leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    q=result["original_site_plant_real_morphological_guild_support"]
    for mode in ("all_cited","largest_archive_excluded"):
        focus=[{
           "source_mode":c["source_mode"],"region":c["region"],
           "botanical_status":c["botanical_status"],
           "real_species_direct_pair":c["real_species_direct_pair"],
           "original_source_visitor_guild":c["original_source_visitor_guild"],
           "n_cited":c["n_positive_cited_relation_units"],
           "n_islands":c["n_documented_islands"],
           "n_distinct_plants":c["n_documented_accepted_plant_species"],
           "n_source_labels":c["n_documented_archive_first_labels"],
           "n_reference_labels":c["n_documented_reference_labels"],
           "n_pollination_class":c["n_explicit_pollination_relations"],
        } for c in q["cells"] if (
             c["source_mode"]==mode and c["region"] in ("GLOBAL_ALL","tropical",
                  "northern_midlatitude") and
             c["botanical_status"] in ("all_recorded","original_native") and
             c["n_positive_cited_relation_units"]>0 and
             c["original_source_visitor_guild"] in (
               "flower_visiting_birds","lepidoptera","bombus",
               "non_bombus_bees"))]
        print("ANONYMOUS_GLOBI_FULL4505_GUILD="+json.dumps(focus,sort_keys=True))
    print("ANONYMOUS_GLOBI_FULL4505_SUMMARY="+json.dumps(dict(
        n_original_botanical_floras=4505,
        n_verified_original_polygons=result[
            "geometry_original_biological_source"][
                "n_original_island_floras_with_exact_hash_rebuilt_polygon"],
        n_original_direct_pair_species=result[
            "n_original_source_highmedium_direct_two_axis_candidate_species"],
        n_islands_with_cited_source=q[
            "n_source_original_islands_with_some_positive_relation"],
        n_island_plant_pairs_with_cited_source=q[
            "n_source_positive_original_island_plant_pairs"],
        n_original_unambiguous_bibliographic_labels=q[
            "n_original_source_bibliographic_keys"],
        original_old_928_552_paired_source_replay_passed=True,
        no_pollen_or_guild_absence_inference=True),sort_keys=True))
if __name__=="__main__":
    main()
