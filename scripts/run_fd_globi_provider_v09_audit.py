"""Run fd source-provider coverage audit from original pinned GloBI v0.9 TSV.

Input target genera are read from the PRE-EXISTING outcome-independent source
genus breadth artifact. Only a privacy-safe JSON receipt is published.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd

from island_v2.fd_globi_provider_v09 import analyze_v09_gzip, VERSION, PIN
from island_v2.fd_globi_provider_family_transfer import family_transfer
from island_v2.fd_globi_species_shift_v09 import (
    select_native_source_names, ExactSpeciesGloBICoverage, GEO_FIELDS
)
from island_v2.fd_globi_species_shift_spatial_v09 import (
    original_candidate_island_locality_gate
)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--interactions",required=True,type=Path)
    p.add_argument("--frozen-source-genus-csv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--original-direct-ledger",type=Path)
    p.add_argument("--original-island-taxonomy",type=Path)
    p.add_argument("--original-progressive-root",type=Path)
    p.add_argument("--qualified-island-channel-states",type=Path)
    p.add_argument("--original-candidate-island-GSHHG-gpkg",type=Path)
    p.add_argument("--original-candidate-island-GSHHG-receipt",type=Path)
    p.add_argument("--chunksize",type=int,default=200_000)
    a=p.parse_args()
    d=pd.read_csv(a.frozen_source_genus_csv)
    needed={"genus","n_reference_channel_units","n_independent_references"}
    if not needed.issubset(d):raise ValueError("Missing frozen source-genus breadth contract")
    if d.genus.isna().any() or d.genus.duplicated().any():
        raise ValueError("Unresolved or duplicate focal plant genus")
    original_ledger=None
    if a.original_direct_ledger is not None:
        original_ledger=pd.read_csv(a.original_direct_ledger,
            usecols=["accepted_species","trait_name","normalized_value",
                     "evidence_scope","evidence_quality"])
    species_evidence=None
    species_selection=None
    private_native_species_island_pairs=None
    if bool(a.original_progressive_root) != bool(a.qualified_island_channel_states):
        raise ValueError("Both original progressive source and qualified island channel cohort required")
    if a.original_progressive_root is not None:
        if original_ledger is None:
            raise ValueError("Species shift evidence requires original species-direct traits")
        root=a.original_progressive_root
        flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
                           usecols=["island_id","accepted_species","origin_status"])
        states=pd.read_csv(a.qualified_island_channel_states,
                           usecols=["island_id","channel_id","channel_state","source_state"])
        sources=pd.read_csv(root/"fixed/source/result/island_source_assignments.csv.gz",
                           usecols=["island_id","source_mode","source_rank","entity_ID"])
        mainland=pd.read_csv(root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
                            usecols=["entity_ID","work_species"])
        names,species_selection,private_native_species_island_pairs=select_native_source_names(
            original_ledger,flora,states,sources,mainland,
            return_private_pairs=True)
        # Cohort is fixed and all identity-level records remain private;
        # any unexpected change means this is NOT the original 54-species test.
        required={"n_islands_with_30plus_species_direct_flower_traits":928,
                  "n_candidate_islands_qualified_Bombus_nondetection_bird_recorded":21,
                  "n_native_same_name_GIFT_source_island_plant_rows":197,
                  "n_distinct_native_same_name_GIFT_source_plant_species":54,
                  "n_islands_with_native_same_name_GIFT_source_plant":13}
        if any(species_selection.get(k)!=v for k,v in required.items()):
            raise ValueError("Original 2026-11-12 source-matched 54-species cohort has drifted")
        # Join ONLY species-direct, unambiguous, immutable flower form
        # measured in the original source ledger (never a guild-derived
        # pollination-syndrome classification).
        d_form=original_ledger.loc[
            original_ledger.trait_name.eq("floral_form") &
            original_ledger.evidence_scope.eq("species_direct") &
            original_ledger.evidence_quality.isin(("high","medium")) &
            original_ledger.accepted_species.isin(names)].copy()
        if d_form.duplicated("accepted_species").any():
            raise ValueError("Multiple incompatible original species-direct floral form records")
        d_form=d_form.dropna(subset=["normalized_value"])
        d_form=d_form.loc[~d_form.normalized_value.astype(str).str.contains("|",regex=False)]
        forms=d_form.set_index("accepted_species").normalized_value.to_dict()
        if set(forms)!=names:
            raise ValueError("54 original candidate plants lack Direct unambiguous flower form")
        species_evidence=ExactSpeciesGloBICoverage(
            names,floral_form_by_species=forms)
        print("Fixed candidate native/source cohort: 54 species, 13 islands; identifiers kept private",flush=True)
    original, after, result=analyze_v09_gzip(
        a.interactions, set(d.genus.astype(str)), chunksize=a.chunksize,
        extra_hook=(species_evidence.add_chunk if species_evidence else None),
        extra_columns=(GEO_FIELDS if species_evidence else ()),
    )
    if bool(a.original_candidate_island_GSHHG_gpkg) != bool(a.original_candidate_island_GSHHG_receipt):
        raise ValueError("Both original candidate exact-GSHHG geometry and receipt required")
    if a.original_candidate_island_GSHHG_gpkg is not None and species_evidence is None:
        raise ValueError("Geospatial visitor audit requires original strict 54 plant species")
    if species_evidence is not None:
        species_payload={
            "selection":species_selection,
            "GloBI_exact_species":species_evidence.summarize(),
            "status":"EXACT_SPECIES_GLOBAL_GLOBI_EVIDENCE_COVERAGE_NOT_ACTUAL_ISLAND_SOURCE_POLLINATOR_SHIFT",
        }
        if a.original_candidate_island_GSHHG_gpkg is not None:
            species_payload["original_island_polygon_locality_overlap"] = (
                original_candidate_island_locality_gate(
                    species_evidence.georeferenced_documented_relations,
                    species_evidence.providers_by_reference,
                    private_native_species_island_pairs,
                    a.original_candidate_island_GSHHG_gpkg,
                    a.original_candidate_island_GSHHG_receipt,
                )
            )
        result["species_direct_mainland_island_shift_evidence_coverage"]=species_payload
    if bool(a.original_direct_ledger) != bool(a.original_island_taxonomy):
        raise ValueError("Need BOTH the original direct ledger and frozen GBIF taxonomy")
    if a.original_direct_ledger:
        for pth in (a.original_direct_ledger, a.original_island_taxonomy):
            if not pth.is_file():
                raise ValueError("Missing original source-backed morphology or taxonomy")
        original_taxonomy=pd.read_csv(a.original_island_taxonomy,
           usecols=["genus","family"])
        result["provider_omitted_family_transfer"]=family_transfer(
            original,after,original_ledger,d,original_taxonomy)

    result.update({
        "n_outcome_blind_frozen_original_genus_list":int(len(d)),
        "v09_full_TSV_source_version":VERSION,
        "v09_full_TSV_SHA256":PIN,
        "focal_island_traits_used":False,
        # Source universe was fixed before looking at morphology. The
        # prediction DOES use Direct morphology in its feature matrix;
        # avoid implying that no phenotype data entered the prediction.
        "source_universe_selected_without_flower_traits":True,
        "direct_flower_traits_used_for_family_prediction":
            bool(a.original_direct_ledger),
        "native_same_species_GIFT_source_evidence_audit_completed":
            bool(species_evidence),
        "original_13_candidate_island_GloBI_point_polygon_preflight_completed":
            bool(species_evidence and a.original_candidate_island_GSHHG_gpkg),
        "provider_family_prediction_rerun_completed":
            bool(result.get("provider_omitted_family_transfer",{}).get(
                "heldout_original_provider_prediction_completed",False)),
        "provider_heldout_effect_on_flower_morphology_NOT_TESTED":
            not result.get("provider_omitted_family_transfer",{}).get(
                "heldout_original_provider_prediction_completed",False),
        "original_archive_source_granularity_unresolved":True,
    })
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
