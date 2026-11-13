"""Protected original 928-island true dated phylogeny placement/opportunity audit.

Published original V.PhyloMaker2 GBOTB TPL/LCVP/WP three phylogenies,
pinned commit and unmodified data. Separately source-pin CoL XR exact
unreviewed source-family candidates and Original FD Direct species. Neither
species-name/tip map nor original island/lineage identities are uploaded.

A backbone tip is evidence of a phylogenetic placement; a genus/family
*insertion opportunity* is a HYPOTHESIS, not a measured sequence tip.
"""
from __future__ import annotations
from collections import Counter
import argparse
import json
from pathlib import Path
import subprocess
import pandas as pd

from island_v2.fd_globi_macro_928_locality_v09 import strict_direct_928_flora
from fd.colxr_family_opportunity import strict_candidate_family_map
from fd.discrete_functional_volume import TRAITS

PINNED_ORIGINAL_DIRECT_SPECIES=2234
TREE_COMMIT="7af3fb5152f691af2e4ec9d5e2e467d1b50505e9"
TREE_RDA_SHA256={
    "TPL":"9ac19e04d6016c9e7ee8758524996f40a089271b2fd34371a5dbd607b5dbf2a3",
    "LCVP":"30fd7e9e67bdfff427b336d3b35983f63e284ef118e91f76a7f51a884c9b4e5b",
    "WP":"7fc97970fcd9e69f6506b0c860b1296f4a797635193a8ce4bcb052a83fd1f965",
}


def direct_triplets(ledger:pd.DataFrame):
    needed={"accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"}
    if not needed.issubset(ledger):
        raise ValueError("Original trait source columns missing")
    s=ledger.loc[
       ledger.trait_name.isin(TRAITS)&
       ledger.evidence_scope.eq("species_direct")&
       ledger.evidence_quality.isin(("high","medium"))&
       ledger.normalized_value.notna()&
       ~ledger.normalized_value.astype(str).str.contains("|",regex=False),
       list(needed)].copy()
    if s.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original Direct accepted species/trait conflict")
    a=s.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    a=a.dropna(subset=list(TRAITS))[list(TRAITS)]
    if len(a)!=PINNED_ORIGINAL_DIRECT_SPECIES:
        raise ValueError("Frozen original global Direct 2234 species drift")
    return a


def trial(root:Path,backbone:Path,tax:Path,reader:Path,tmp:Path)->dict:
    orig=pd.read_csv(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        usecols=["accepted_species","trait_name","normalized_value",
                 "evidence_scope","evidence_quality"],dtype=str)
    triples=direct_triplets(orig)
    original_flora=pd.read_csv(root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        usecols=["island_id","accepted_species","origin_status"],dtype=str)
    islg, gate=strict_direct_928_flora(orig,original_flora)
    if gate["n_original_eligible_islands"]!=928:
        raise ValueError("Frozen original 928 original flora ID universe changed")
    geo=pd.read_csv(root/"fixed/isolation/results/purpose_shortest_island_data.csv",
        usecols=["island_id","analysis_regime","spatial_block"],dtype=str)
    if geo.island_id.duplicated().any():
        raise ValueError("Original island source geography duplicate")
    full=triples.copy()
    source_family=pd.read_csv(tax,dtype=str)
    family_map, original_mapping_audit=strict_candidate_family_map(
        source_family,set(triples.index))
    if len(family_map)!=2190:
        raise ValueError("Frozen original 2190 independently mapped family candidates drift")
    original_phenotypes=list(map(tuple,triples.to_numpy(dtype=str)))
    carrier=Counter(original_phenotypes)
    tiers={}
    low=0
    high=0
    for species,profile in zip(triples.index,original_phenotypes):
        n=carrier[profile]
        tier=("rare" if n<=2 else "common" if n>=26 else "middle")
        tiers[species]=tier
        low+=tier=="rare"
        high+=tier=="common"
    if low!=40 or high!=1734:
        raise ValueError("Original rare 40 vs common 1734 source species cohort drift")
    private=tmp/"internal_only"
    private.mkdir(parents=True,exist_ok=True)
    input_csv=private/"original_source_plant_species_and_family.tsv"
    output_csv=private/"original_tree_placement_by_species.tsv"
    row=pd.DataFrame({
        "accepted_species":list(triples.index),
        "source_family":[family_map.get(name,"") for name in triples.index],
        "source_genus":[name.split()[0] for name in triples.index],
        "rarity_tier":[tiers[name] for name in triples.index],
    })
    row.to_csv(input_csv,sep="\t",index=False)
    # It is strictly forbidden to upload either privacy-sensitive internal TSV.
    subprocess.run(["Rscript",str(reader),str(backbone/"data"),
                    str(input_csv),str(output_csv)],
                   check=True,timeout=180)
    matches=pd.read_csv(output_csv,sep="\t",dtype=str)
    if len(matches)!=2234*3 or matches.duplicated(
       ["accepted_species","tree_version"]).any():
        raise ValueError("One original plant per source tree taxonomy match lost")
    if set(matches.accepted_species)!=set(triples.index):
        raise ValueError("Original phylogeny source species differ from Direct input")
    expected={
      "exact_backbone_tip","genus_backbone_insertion_only",
      "family_backbone_insertion_only","unresolved_no_valid_source_family"
    }
    if not set(matches.supported_species_placement).issubset(expected):
        raise ValueError("Unexpected undocumented source phylogenetic graft")
    reporting=[]
    for version,frame in matches.groupby("tree_version"):
        states=frame.supported_species_placement.value_counts().to_dict()
        grouped=dict(tree_version=version,
           n_original_full_megatree_tips={"TPL":74529,"LCVP":73420,"WP":72570}[version],
           n_original_global_direct=2234,
           n_source_backbone_original_species_names_matched_exactly=int(
             frame.backbone_tip_and_info_names_exact.eq("TRUE").sum() if
             frame.backbone_tip_and_info_names_exact.isin(["TRUE","FALSE"]).any()
             else frame.backbone_tip_and_info_names_exact.eq("True").sum()),
           n_source_backbone_exact_tip_and_agreed_family_genus=int(
             states.get("exact_backbone_tip",0)),
           n_missing_exact_with_valid_source_genus_in_backbone=int(
             states.get("genus_backbone_insertion_only",0)),
           n_missing_genus_with_valid_family_in_backbone=int(
             states.get("family_backbone_insertion_only",0)),
           n_unplaceable_or_taxonomically_unverified=int(
             states.get("unresolved_no_valid_source_family",0)),
           placement_of_unseen_species_is_not_empirical_sequence_phylogeny=True)
        for tier,sub in frame.groupby("rarity_tier"):
            reporting.append(dict(type="original_source_species_tier",
                tier=tier,**grouped,
                n_tier_species=int(len(sub)),
                n_tier_exact_phylogenetic_tip=int(
                    sub.supported_species_placement.eq("exact_backbone_tip").sum()),
                n_tier_genus_insertion_hypothesis_only=int(
                    sub.supported_species_placement.eq("genus_backbone_insertion_only").sum()),
                n_tier_family_insertion_hypothesis_only=int(
                    sub.supported_species_placement.eq("family_backbone_insertion_only").sum())))
    islg=islg.merge(geo,on="island_id",how="left",validate="many_to_one")
    if islg.analysis_regime.isna().any():
        raise ValueError("Some original eligible 928 islands lack geographic region")
    coverage=[]
    for version,frame in matches.groupby("tree_version"):
        member=islg.merge(frame[["accepted_species","supported_species_placement",
                                "rarity_tier"]],
            on="accepted_species",how="left",validate="many_to_one")
        for regime,g in member.groupby("analysis_regime"):
            n=g.groupby("island_id").accepted_species.nunique()
            exact=g.loc[g.supported_species_placement.eq(
                "exact_backbone_tip")].groupby("island_id").accepted_species.nunique().reindex(
                    n.index,fill_value=0)
            genus=g.loc[g.supported_species_placement.isin(
                ("exact_backbone_tip","genus_backbone_insertion_only"))].groupby(
                    "island_id").accepted_species.nunique().reindex(n.index,fill_value=0)
            native=g.loc[g.origin_status.eq("native")]
            strict_native=native.groupby("island_id").accepted_species.nunique()
            rare=g.loc[g.rarity_tier.eq("rare")]
            rare_n=rare.groupby("island_id").accepted_species.nunique().reindex(n.index,fill_value=0)
            rare_exact=rare.loc[rare.supported_species_placement.eq(
                "exact_backbone_tip")].groupby("island_id").accepted_species.nunique().reindex(n.index,fill_value=0)
            coverage.append(dict(
                tree_version=version,regime=regime,n_original_direct_30plus_islands=len(n),
                n_spatial_blocks=int(g.spatial_block.nunique()),
                n_islands_exact_tree_tip_5plus=int(exact.ge(5).sum()),
                n_islands_exact_tree_tip_10plus=int(exact.ge(10).sum()),
                n_islands_exact_tree_tip_20plus=int(exact.ge(20).sum()),
                n_islands_exact_or_genus_hypothesis_20plus=int(genus.ge(20).sum()),
                n_islands_exact_or_genus_hypothesis_30plus=int(genus.ge(30).sum()),
                median_share_direct_island_species_exact_tip=float(
                    (exact/n).median()),
                median_share_direct_island_species_exact_or_genus=float(
                    (genus/n).median()),
                n_islands_strict_native_direct_30plus=int(strict_native.ge(30).sum()),
                n_islands_rare_direct_species_present=int(rare_n.gt(0).sum()),
                n_islands_rare_direct_species_and_exact_tree_tip=int(rare_exact.gt(0).sum()),
                exact_tip_only_macro_group_not_representative_of_all_islands=True
            ))
    return {
      "status":"PINNED_GBOTB2_ORIGINAL_DIRECT_928_PHYLOGENY_COVERAGE_ONLY",
      "backbone_commit":TREE_COMMIT,
      "backbone_rda_sha256":TREE_RDA_SHA256,
      "original_CoLXR_exact_candidate_family_species":len(family_map),
      "original_CoLXR_candidate_family_not_promoted":True,
      "original_direct":2234,"original_eligible_islands":928,
      "n_original_rare_profile_carrier_species":low,
      "n_original_common_profile_carrier_species":high,
      "by_original_species_rarity_tier_and_tree":reporting,
      "by_island_regime_and_tree":coverage,
      "source_tree_branch_lengths_present_but_imputed_tips_have_no_experimental_dates":True,
      "within_species_trait_changes_not_measurable":True,
      "original_phylogenetic_branch_length_effect_not_yet_fit":True,
      "phylogeny_tree_taxonomy_triple_original_sites_never_exported":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",type=Path,required=True)
    p.add_argument("--colxr-crosswalk",type=Path,required=True)
    p.add_argument("--backbone-checkout",type=Path,required=True)
    p.add_argument("--r-matcher",type=Path,required=True)
    p.add_argument("--workdir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=trial(a.original_root,a.backbone_checkout,a.colxr_crosswalk,
                 a.r_matcher,a.workdir)
    output=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for prohibited in ('"accepted_species"','"source_genus"','"source_family"',
      '"island_id"','"family"','"genus"','"species"','"tip_label"','"source_signature"'):
        if prohibited in output:
            raise ValueError("Private original plant/tree mapping leaked into public output")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(output+"\n")
    print(json.dumps({
       "status":result["status"],"n_islands":result["original_eligible_islands"],
       "n_original_direct":result["original_direct"],
       "global_phylo_source_tiers":[{
         "tree":t["tree_version"],"tier":t["tier"],
         "tier_species":t["n_tier_species"],
         "exact_backbone_tip_species":t["n_tier_exact_phylogenetic_tip"],
         "potential_genus_only":t["n_tier_genus_insertion_hypothesis_only"],
       }for t in result["by_original_species_rarity_tier_and_tree"]],
       "n_region_coverage_panels":len(result["by_island_regime_and_tree"]),
    },sort_keys=True))
if __name__=="__main__":
    main()
