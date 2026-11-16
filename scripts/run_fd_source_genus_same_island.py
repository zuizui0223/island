"""Protected original 4,505 GIFT floral congener donor audit across source modes.

Runs on locked source botanical taxa, mainland regional lists and Direct
species phenotypes within island CI. Never exports names, island identifiers,
source GIFT regions, geography or phylogenetic ancestry assertions.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    source_sectors, status_groups,
)
from fd.source_genus_same_island import (
    MORPHOLOGY, SCOPES, REGIMES,
    compare_island_source_genera_scopes, summarize_source_genus_frame,
)

ORIGINAL_HASHES = {
    "direct_only_trait_ledger.csv.gz":
        "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
    "chapter1_status_flora.csv.gz":
        "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
}
MODE_RANKS = {"geo_k5": 5, "geo50_climate10": 10}


def sha(path: Path) -> str:
    q=sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            q.update(block)
    return q.hexdigest()


def build_private_original_source_genus_frame(
    ledger: pd.DataFrame, flora: pd.DataFrame,
    geography: pd.DataFrame, assignments: pd.DataFrame,
    mainland: pd.DataFrame,
) -> pd.DataFrame:
    original = status_groups(flora)
    if len(original) != 4505:
        raise ValueError("original 4505 GIFT botanical floras not present")
    if (geography.island_id.duplicated().any()
            or geography.island_id.nunique() < 8200):
        raise ValueError("independently published 8264-island crosswalk invalid")
    georef=geography.set_index("island_id").reindex(sorted(original))
    if int(georef.analysis_regime.notna().sum()) != 4504:
        raise ValueError("original one geography-unresolved island changed")
    if assignments.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("GIFT original mainland region ranking duplicated")
    selected=assignments.loc[assignments.source_mode.isin(MODE_RANKS)].copy()
    selected.source_rank=pd.to_numeric(selected.source_rank,errors="raise").astype(int)
    for mode,k in MODE_RANKS.items():
        if not selected.loc[selected.source_mode.eq(mode),
                            "source_rank"].between(1,k).all():
            raise ValueError("original ranked mainland candidate scope drift")
    expected_assign={"geo_k5":41325,"geo50_climate10":82190}
    if {
        mode:int(selected.source_mode.eq(mode).sum()) for mode in MODE_RANKS
    }!=expected_assign:
        raise ValueError("original frozen geographical/climate source ranks changed")
    ranks=selected.sort_values("source_rank").groupby(
        ["island_id","source_mode"],sort=True
    ).entity_ID.agg(tuple).to_dict()
    if mainland.duplicated(["entity_ID","work_species"]).any():
        raise ValueError("original GIFT entity×accepted taxon duplication")
    regions=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    if len(regions)!=748:
        raise ValueError("original GIFT 748 native mainland entities drifted")

    @lru_cache(maxsize=15)
    def source_roster(ids: tuple[str,...]) -> frozenset[str]:
        return frozenset().union(*(regions.get(entity,set()) for entity in ids))

    sectors=source_sectors(strict_traits(ledger))
    if set(MORPHOLOGY)-set(sectors):
        raise ValueError("original Direct real joint-flower morphology unavailable")
    rows=[]
    climate_count=0
    for island in sorted(original):
        current=original[island]
        full=current["full"]
        bymode={}
        rawmodes={}
        for mode,k in MODE_RANKS.items():
            rank=ranks.get((island,mode),())
            if rank and len(rank)!=k:
                raise ValueError("partial native mainland source rank for original island")
            if rank:
                rawmodes[mode]=source_roster(rank)
                bymode[mode]=set(rawmodes[mode])-full
            else:
                bymode[mode]=set()
        if "geo50_climate10" in rawmodes:
            climate_count+=1
        scope_plant_sets={
            "all_recorded":full,
            "strict_native":current["native"],
            "strict_native_nonendemic":current["native_nonendemic"],
            # A separate, carefully selected plant group: status-known native
            # nonendemics whose accepted species name is in BOTH independent
            # GIFT source-region rosters. It is NOT the main study population.
            "native_nonendemic_source_listed_both":
                current["native_nonendemic"]
                & rawmodes.get("geo_k5",frozenset())
                & rawmodes.get("geo50_climate10",frozenset()),
        }
        if set(scope_plant_sets)!=set(SCOPES):
            raise AssertionError("four original provenance strata incomplete")
        region=(str(georef.at[island,"analysis_regime"])
                if pd.notna(georef.at[island,"analysis_regime"])
                else "geography_unresolved")
        if region not in REGIMES:
            raise ValueError("unknown source GSHHG region")
        block=georef.at[island,"spatial_block"]
        block=str(block) if pd.notna(block) else None
        nearest=ranks.get((island,"geo_k5"),())
        signature=nearest[0] if len(nearest)==5 else None
        for sector in MORPHOLOGY:
            scored=sectors[sector]
            results=compare_island_source_genera_scopes(
                full_island_names=full,
                scope_names=scope_plant_sets,
                geo5_candidates=bymode["geo_k5"],
                climate10_candidates=bymode["geo50_climate10"],
                species_states=scored,
            )
            for scope,result in results.items():
                row=dict(
                    island_id=island,regime=region,scope=scope,sector=sector,
                    spatial_block=block,source_signature=signature,
                )
                for key,val in result.items():
                    if key.startswith("shared_genus_"):
                        for name,value in val.items():
                            row[key+"__"+name]=value
                    else:
                        row[key]=val
                rows.append(row)
    if climate_count!=4485:
        raise ValueError("climate10 main original 4485/4505 source coverage changed")
    if len(rows)!=4505*len(SCOPES)*len(MORPHOLOGY):
        raise AssertionError("original all-island source population was narrowed")
    return pd.DataFrame(rows)


def run(root:Path,geography:Path)->dict:
    ledger_path=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    flora_path=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for f in (ledger_path,flora_path):
        if sha(f)!=ORIGINAL_HASHES[f.name]:
            raise ValueError("source sha provenance not original")
    ledger=pd.read_csv(ledger_path,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"
    ])
    flora=pd.read_csv(flora_path,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"
    ])
    geo=pd.read_csv(geography,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"
    ])
    assignments=pd.read_csv(
        root/"fixed/source/result/island_source_assignments.csv.gz",
        dtype=str,usecols=["island_id","source_mode","source_rank","entity_ID"]
    )
    mainland=pd.read_csv(
        root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        dtype=str,usecols=["entity_ID","work_species"]
    )
    private=build_private_original_source_genus_frame(
        ledger,flora,geo,assignments,mainland
    )
    summary=summarize_source_genus_frame(
        private,population=4505,bootstrap=999
    )
    if (summary["n_original_botanical_island_floras"]!=4505 or
        summary["n_regime_scope_morphology_panels"]!=100):
        raise AssertionError("original cohort/genus support audit drifted")
    return dict(
        **summary,
        n_GIFT_native_mainland_source_entities=748,
        original_GIFT_source_region_assignments_geo5=41325,
        original_GIFT_source_region_assignments_climate10=82190,
        original_botanical_flora_sha256=sha(flora_path),
        original_Direct_species_ledger_sha256=sha(ledger_path),
        actual_morphology_no_imputed_joint_traits=True,
        original_4_505_flora_and_taxa_never_exported=True,
    )


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--geology-crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    obj=run(a.original_root,a.geology_crosswalk)
    body=json.dumps(obj,sort_keys=True,indent=2,
                    allow_nan=False,ensure_ascii=False)
    forbidden=(
        "island_id","accepted_species","entity_ID","work_species",
        "source_signature","spatial_block","source_region_id",
        "latitude","longitude","GBIF_id","taxon_id","observer_id",
        "source_genus","source_family",
    )
    if any('"'+key+'"' in body for key in forbidden):
        raise ValueError("source taxon/island provenance leaked to public artifact")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    positive=sum(
        r["one_donor__full_matched_sample__spatial_block"]["status"] ==
        "ORIGINAL_CONGENERIC_SOURCE_DESCRIPTIVE_NOT_CAUSAL"
        for r in obj["panels"]
    )
    print(json.dumps(dict(
        status=obj["status"],
        original_islands=obj["n_original_botanical_island_floras"],
        source_conditional_panels=obj["n_regime_scope_morphology_panels"],
        panels_with_minimum_block_support=positive,
        nonidentification_of_evolution_and_pollinator_effect=True,
    ),sort_keys=True))


if __name__=="__main__":
    main()
