"""Original 4,505 GIFT botanical flora x immutable European DoPI v2 visitors.

This evaluates guild-conditional flower morphology in independently published
positive European visitor records, NOT visits observed on the GIFT islands.
The original accepted plant and publication identifiers remain protected.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups
from fd.dopi_real_motif_multiresolution import make_real_multiresolution_motifs
from fd.dopi_study_genus_equal import genus_equal_dopi_audit
from fd.dopi_site_period_guild import site_time_positive_index,site_period_audit
from fd.dopi_positive_coaccess import summarize_positive_coaccess
from fd.dopi_guild_portfolio_effort_null import (
    positive_visitor_effort_index,guild_richness_motif_audit
)
from fd.dopi_degree_fixed_visitor_null import audit_original_degree_fixed_networks
from fd.dopi_gift_guild_motif_bridge import (
    screen_positive_dopi_rows,summarize_study_guild_motifs,
)
from audit_dopi_v2_positive_visitors import identify_visitor_group

SOURCE_HASHES={
  "direct_only_trait_ledger.csv.gz":
    "9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
  "chapter1_status_flora.csv.gz":
    "6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
  "DoPI_dataset_2.0.csv":
    "b95d3084399ed2c0978af2485801108a27d14e843b9fed2611934ea582b745c8",
}
ORIGINAL_DOPI_CSV_SIZE=384896230


def digest(path:Path)->str:
    h=sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def original_island_flora_dopi_guild_motifs(
    *,direct:pd.DataFrame,flora:pd.DataFrame,
    public_dopi_csv:Path,draws:int=999,
)->dict:
    plants=status_groups(flora)
    if len(plants)!=4505:
        raise ValueError("original 4505 GIFT accepted botanical island population changed")
    all_accepted=set().union(*(z["full"] for z in plants.values()))
    if not all_accepted:
        raise ValueError("original island plant species rosters unavailable")
    motifs=make_real_multiresolution_motifs(strict_traits(direct))
    if len(motifs)!=27:
        raise ValueError("original GIFT full 4 single, 11 paired, 12 three/four-axis floral motifs missing")
    if digest(public_dopi_csv)!=SOURCE_HASHES["DoPI_dataset_2.0.csv"]:
        raise ValueError("DoPI original v2 public source SHA256 drifted")
    if public_dopi_csv.stat().st_size!=ORIGINAL_DOPI_CSV_SIZE:
        raise ValueError("DoPI original v2 exact source byte length changed")

    with public_dopi_csv.open("r",encoding="utf-8-sig",
                              errors="replace",newline="") as f:
        reader=csv.DictReader(f)
        index=screen_positive_dopi_rows(
            reader,columns=set(reader.fieldnames or []),
            island_plant_names=all_accepted,motifs=motifs,
            visitor_classifier=identify_visitor_group,
        )
    report=summarize_study_guild_motifs(index,motifs=motifs,draws=draws)
    # Repeat the exact public positive source stream with genuine publication,
    # published-site and observation-year/month constraints; never infer a
    # location or guild non-detection from absent metadata.
    with public_dopi_csv.open("r",encoding="utf-8-sig",
                              errors="replace",newline="") as handle:
        reader2=csv.DictReader(handle)
        geo_index=site_time_positive_index(
            reader2,fields=set(reader2.fieldnames or []),
            island_names=all_accepted,motifs=motifs,
            classifier=identify_visitor_group,
        )
    report["same_published_site_period"] = site_period_audit(
        geo_index,motifs=motifs,draws=draws,
    )
    # Positive visits by BOTH real guilds on the SAME accepted plant species,
    # inside an exact original publication+published site+year/month.
    # Never turn a missing guild record into ecological absence.
    report["positive_multiguild_flower_coaccess"] = (
        summarize_positive_coaccess(
            geo_index,motifs=motifs,draws=draws,
        )
    )
    # Fixed observation-exposure null: each plant keeps its distinct
    # recorded visitor taxon count, each original published site-month keeps
    # its empirically documented guild distribution. No unreported visitor
    # or a source non-visit is inferred.
    with public_dopi_csv.open("r",encoding="utf-8-sig",
                              errors="replace",newline="") as handle:
        reader3=csv.DictReader(handle)
        effort_index=positive_visitor_effort_index(
            reader3,fields=set(reader3.fieldnames or []),
            GIFT_accepted_names=all_accepted,motifs=motifs,
            visitor_classifier=identify_visitor_group,
        )
    if effort_index["source_row_counts"]["original_DoPI_v2_source_rows"]!=676900:
        raise ValueError("source DoPI raw positive visitor effort rows changed")
    report["visitor_guild_portfolio_effort_null"] = (
        guild_richness_motif_audit(effort_index,motifs=motifs,draws=draws)
    )
    # STRONGER falsification: preserve degree of EACH actual documented
    # positive visitor species, not merely its guild total, and exact number
    # of different recorded visitors per plant. This is NOT a true field
    # preference/absence/pollen efficacy counterfactual.
    report["visitor_species_degree_fixed_network_null"] = (
        audit_original_degree_fixed_networks(
            effort_index,motifs=motifs,draws=draws,site_draws=80,
        )
    )
    if report["same_published_site_period"]["n_original_DOPI_source_rows"]!=676900:
        raise ValueError("original DoPI 2.0 source row count changed in site/time audit")

    # New source-study/plant-genus equal-weight validation: same private
    # documented positive units; no source plant/study identifiers exported.
    report["within_publication_equal_genus"]=(
        genus_equal_dopi_audit(index,motifs=motifs,draws=draws)
    )
    if report["n_original_DOPI_source_rows"]!=676900:
        raise ValueError("published DoPI source-original row count changed")
    report.update(original_GIFT_flora_sha256=SOURCE_HASHES["chapter1_status_flora.csv.gz"],
                  original_GIFT_Direct_phenotype_sha256=SOURCE_HASHES["direct_only_trait_ledger.csv.gz"],
                  original_DoPI_v2_csv_sha256=SOURCE_HASHES["DoPI_dataset_2.0.csv"],
                  original_DoPI_v2_release_doi="10.5281/zenodo.22876365",
                  original_GIFT_botanical_accepted_names_total=len(all_accepted))
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--original-root",required=True,type=Path)
    parser.add_argument("--dopi-csv",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    a=parser.parse_args()
    root=a.original_root
    ld=root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz"
    fl=root/"fixed/canonical/input/chapter1_status_flora.csv.gz"
    for path in (ld,fl):
        if digest(path)!=SOURCE_HASHES[path.name]:
            raise ValueError("original frozen GIFT plant/trait archive hash changed")
    direct=pd.read_csv(ld,dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality",
    ])
    flora=pd.read_csv(fl,dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status",
    ])
    result=original_island_flora_dopi_guild_motifs(
        direct=direct,flora=flora,public_dopi_csv=a.dopi_csv,draws=999,
    )
    raw=json.dumps(result,sort_keys=True,ensure_ascii=False,
                   indent=2,allow_nan=False)
    forbidden=("island_id","accepted_species","plant_species",
               "pollinator_species","work_species","source_signature",
               "spatial_block","doi","publication","entity_ID",
               "latitude","longitude","plant_name","taxon_id",
               "original_study_key","original_species")
    if any('"'+s+'"' in raw for s in forbidden):
        raise ValueError("source-protected botanical/DoPI exact identifier leak")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(raw+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=result["status"],
        n_original_botanical_islands=result["n_original_GIFT_botanical_islands"],
        n_publications=result["n_unique_original_publication_keys"],
        n_matched_island_flora_accepted_plants=result[
            "n_GIFT_accepted_island_plant_species_matching_DoPI"],
        n_study_motif_panels=result["n_qualification_panels"],
        supported_motif_guild_study_tests=sum(
            p["all_study"]["status"]==
            "DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION"
            for p in result["panels"]
        ),
        n_equal_genus_supported=sum(
            q["genus_equal_all"]["status"]=="DESCRIPTIVE_POSITIVE_VISITOR_STUDY_EQUAL_GENUS_FIXED"
            for q in result["within_publication_equal_genus"]["panels"]
        ),
        n_variable_genus_supported=sum(
            q["genus_equal_informative"]["status"]=="DESCRIPTIVE_POSITIVE_VISITOR_STUDY_EQUAL_GENUS_FIXED"
            for q in result["within_publication_equal_genus"]["panels"]
        ),
        n_same_published_site_year_supported=sum(
            p["same_site_positive_visitor"]["status"]==
            "DESCRIPTIVE_POSITIVE_GUILD_COMMON_SITE_PERIOD"
            for p in result["same_published_site_period"]["panels"]
            if p["resolution"]=="same_study_site_year"
        ),
        n_positive_coaccess_supported=sum(
            p["co_documented_vs_positive_union"]["status"]==
            "DESCRIPTIVE_REAL_GUILD_CO_DOCUMENTATION"
            for p in result["positive_multiguild_flower_coaccess"]["panels"]
        ),
        n_genus_fixed_coaccess_supported=sum(
            p["genus_fixed_co_documentation_excess"]["status"]==
            "DESCRIPTIVE_GENUS_CONDITIONED_CO_DOCUMENTATION_EXCESS"
            for p in result["positive_multiguild_flower_coaccess"]["panels"]
        ),
        n_guild_richness_fixed_exposure_supported=sum(
            p["study_site"]["status"]==
            "DESCRIPTIVE_EFFORT_CONDITIONED_POSITIVE_VISITOR_GUILD_PORTFOLIO"
            for p in result["visitor_guild_portfolio_effort_null"]["panels"]
        ),
        n_guild_richness_genus_fixed_supported=sum(
            p["same_genus_study_site"]["status"]==
            "DESCRIPTIVE_EFFORT_CONDITIONED_POSITIVE_VISITOR_GUILD_PORTFOLIO"
            for p in result["visitor_guild_portfolio_effort_null"]["panels"]
        ),
        n_degree_fixed_network_motifs_supported=sum(
            p["status"]==
            "DESCRIPTIVE_SOURCE_VISITOR_DEGREE_FIXED_FLORAL_GUILD_ASSOCIATION"
            for p in result["visitor_species_degree_fixed_network_null"]["panels"]
        ),
        only_positive_study_guild_comparisons=True
    ),sort_keys=True))


if __name__=="__main__":
    main()
