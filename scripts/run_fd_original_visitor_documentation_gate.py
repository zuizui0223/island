"""Audit ONLY original anonymous source-verified GIFT GloBI/WOL artifacts.

This stage does not download raw taxa, invent extra observations or attempt
independence-corrected ecology. Input originals have already passed protected
2026-11-18 CI, with SHA-locked GIFT 4505 source and 928 historical replay.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from fd.visitor_documentation_support import audit_original_guild_documentation

GIFT_SHA="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"
DIRECT_SHA="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"
GLOBI_SHA="095e57398420d1e0746e61233c796beedcfb8e03bd0c6a732aa17ae62043bd17"
WOL_SHA="f0109a68510bf347006aa97582dd448d94daf29b91b9061e7d1c2a88850856f3"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--globi",type=Path,required=True)
    ap.add_argument("--wol",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    g=json.loads(a.globi.read_text())
    w=json.loads(a.wol.read_text())
    for label,x,y in (
        ("GIFT floras",g["original_source_GIFT_flora_SHA256"],GIFT_SHA),
        ("Direct floral traits",g["original_source_Direct_traits_SHA256"],DIRECT_SHA),
        ("GloBI original v0.9",g["original_source_GloBI_v09_SHA256"],GLOBI_SHA),
        ("Web of Life catalogue",w["source_web_of_life_catalogue_sha256"],WOL_SHA),
    ):
        if x!=y:raise ValueError("Original protected source/receipt changed: "+label)
    if (g["n_original_GloBI_rows_checked"]!=24577183 or
        w["full_original_botanical_flora_island_visitor_evidence"][
            "n_original_species_network_plant_guild_units"]!=2529):
        raise ValueError("Original anonymous animal documentation support changed")
    output=audit_original_guild_documentation(g,wol=w)
    by={(x["region"],x["botanical_status"]):x for x in
        output["source_qualification_panels"]}
    if (by[("GLOBAL_ALL","all_recorded")][
        "n_preliminary_documentation_gates_passing"]!=12 or
        by[("GLOBAL_ALL","all_recorded")][
        "n_provider_removal_documentation_gates_passing"]!=0 or
        by[("tropical","all_recorded")][
        "n_preliminary_documentation_gates_passing"]!=0):
        raise ValueError("Original pre/post provider gate source counts do not reconcile")
    output["pinned_protected_source_run_ids"]=dict(
        GloBI=44301575888,WOL=44299833563)
    output["original_artifact_archive_sha256_reported_by_GitHub"]=dict(
        GloBI="bf680fdbd32d999d83f0bbfed2fc83f34b9a55e7f32a05ad8ad458997896c241",
        WOL="4febcabe8aff19d9ea1a4edd4c9877c1e2e87c4e2fc9516ea1a562d9e2941099")
    raw=json.dumps(output,indent=2,sort_keys=True,allow_nan=False)
    for key in ("island_id","accepted_species","reference_id","provider_id",
                "latitude","longitude","source_genus","source_region_id",
                "network_id","spatial_block","work_species","entity_ID"):
        if '"'+key+'"' in raw:
            raise ValueError("Anonymous source gate contains original private data")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    for q in output["source_qualification_panels"]:
        print("GUILD_DOC_GATE="+json.dumps({
            k:q[k] for k in (
            "region","botanical_status",
            "n_preselected_animal_guild_x_true_morphology_endpoints",
            "n_preliminary_documentation_gates_passing",
            "n_provider_removal_documentation_gates_passing")},
            sort_keys=True))
    print(json.dumps(dict(status=output["status"],floras=4505,
        n_globi_geographically_cited_original_islands=output[
            "n_original_GloBI_positive_documented_islands"],
        no_effective_pollination_inferred=True)))
if __name__=="__main__":
    main()
