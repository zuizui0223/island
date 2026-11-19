"""Run original GloBI v0.9 source+GSHHG polygon exact-site visitor evidence
against genuine SAME accepted species floral symmetry/form × Direct tube depth.

Every original island, plant, source study and local coordinate ID is kept
inside this SHA-frozen GitHub Actions workspace. Only anonymous availability
cells are exported. NOT an actual effective pollination or historical shift.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd_island_functional_beta import strict_traits
from island_v2.fd_globi_species_shift_v09 import (
    ExactSpeciesGloBICoverage,GEO_FIELDS,
)
from island_v2.fd_globi_provider_v09 import (
    RAW,PIN,sha256_file,
)
from island_v2.fd_globi_macro_928_locality_v09 import (
    strict_direct_928_flora,macro_exact_928_locality_gate,
)
from fd.globi_local_real_syndrome_gate import (
    real_same_species_island_visitor_audit,
)

PLANTS_SHA="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"
DIRECT_SHA="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"


def sha(p:Path)->str:
    d=sha256()
    with p.open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):
            d.update(z)
    return d.hexdigest()


def study(original_root:Path,original_GloBI:Path,
          geometry:Path,geometry_receipt:Path,*,chunksize=200000)->dict:
    original_ledger=original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    original_plants=original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    if sha(original_ledger)!=DIRECT_SHA or sha(original_plants)!=PLANTS_SHA:
        raise ValueError("original botanical source Direct and full GIFT flora changed")
    if sha256_file(original_GloBI)!=PIN:
        raise ValueError("original pinned GloBI v0.9 SUPPORTS source SHA drifted")
    traits=pd.read_csv(original_ledger,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(original_plants,dtype=str,usecols=[
        "island_id","accepted_species","origin_status"
    ])
    eligible,original_frame_receipt=strict_direct_928_flora(traits,flora)
    source=ExactSpeciesGloBICoverage(set(eligible.accepted_species))
    fields=set(RAW)|set(GEO_FIELDS)
    chunks=0
    for chunk in pd.read_csv(
        original_GloBI,sep="\t",compression="gzip",dtype=str,
        usecols=lambda x:x in fields,
        low_memory=False,chunksize=chunksize,
    ):
        source.add_chunk(chunk)
        chunks+=1
    if source.raw["n_original_GloBI_rows_read"]!=24577183:
        raise ValueError("original 24,577,183 source GloBI relation rows drifted")
    bad={ref for ref,p in source.providers_by_reference.items()
         if len(p)>1}
    eligible_units=[r for r in source.units if r[2] not in bad]
    counts=Counter(row[3] for row in eligible_units)
    if not counts:
        raise ValueError("original globally selected GloBI source archives missing")
    largest=sorted(counts,key=lambda k:(-counts[k],k))[0]

    original_geo,private=macro_exact_928_locality_gate(
        source.georeferenced_documented_relations,
        source.providers_by_reference,eligible,
        geometry,geometry_receipt,
        largest_global_archive_first_label=largest,
        return_private_linked=True,
    )
    if (original_geo["n_distinct_islands_with_species_locality_match"]!=109
        or original_geo["n_distinct_island_species_pairs_with_record"]!=552
        or original_geo["n_strict_native_flora_island_species_pairs_with_visitor_relation"]!=157
        or original_geo["n_islands_with_native_plant_visitation_reference"]!=40):
        raise ValueError("old original 928 exact flora×GSHHG×GloBI locality counts drifted")
    audited=real_same_species_island_visitor_audit(
        eligible,private,strict_traits(traits),largest,
        original_global_flora_n=4505,
        verified_polygons_n=original_geo["n_islands_matching_original_GSHHG_geometric_id"],
    )
    return {
        "status":"SOURCE_VERIFIED_ORIGINAL_GIFT_GLOBI_REAL_JOINT_LOCAL_VISITOR_GATE",
        "n_original_main_study_botanical_floras":4505,
        "n_morphology_complete_locality_subset_islands":928,
        "source_original_flora_SHA256":sha(original_plants),
        "source_original_species_Direct_ledger_SHA256":sha(original_ledger),
        "source_original_GloBI_v09_SHA256":sha256_file(original_GloBI),
        "source_GloBI_rows_checked":int(source.raw["n_original_GloBI_rows_read"]),
        "source_GloBI_input_chunks":chunks,
        "original_GloBI_macro_locality_coverage":original_geo,
        "same_plant_real_architecture_guild_coverage":audited,
        "identity_and_exact_site_data_never_serialized":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=Path)
    p.add_argument("--globi-original-tsv-gz",required=True,type=Path)
    p.add_argument("--original-geometry",required=True,type=Path)
    p.add_argument("--geometry-receipt",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    result=study(a.original_root,a.globi_original_tsv_gz,
                 a.original_geometry,a.geometry_receipt)
    raw=json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2,
                   allow_nan=False)
    for name in (
        "island_id","accepted_species","reference_id","provider_id",
        "work_species","entity_ID","latitude","longitude",
        "visitor_species","source_genus","source_family","spatial_block"
    ):
        if '"'+name+'"' in raw:
            raise ValueError("protected GloBI/locality/plant/species identifiers leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    q=result["same_plant_real_architecture_guild_coverage"]
    print(json.dumps({
        "status":result["status"],
        "n_all_original_study_floras":4505,
        "n_original_direct_locality_eligible_islands":928,
        "n_geometry_matching_original_islands":927,
        "n_flower_guild_origin_source_support_cells":q["n_source_guild_trait_plant_origin_support_cells"],
    },sort_keys=True))


if __name__=="__main__":
    main()
