"""Inspect source-verified anonymous GIFT vs DoPI floral motif direction matching.

The original 4505-island GIFT floral composition and the European DoPI
butterfly-vs-Bombus positive visitor study contrasts are DIFFERENT estimands.
A shared sign does not establish effective pollination or historical island
pollinator shifts. Plant taxonomy / source studies / islands overlap.
"""
from __future__ import annotations
from collections import Counter
from itertools import permutations
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
from scipy.stats import spearmanr,pearsonr

STRUCTURAL=("actinomorphic_deep","zygomorphic_deep","tubular_bell_deep",
            "actinomorphic_tubular_bell","zygomorphic_tubular_bell")
COLOUR_FORM=("red_pink_tubular_bell","yellow_orange_tubular_bell",
             "white_tubular_bell")
MOTIFS=STRUCTURAL+COLOUR_FORM
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
MODES=("geo5","climate10")
DIRECT_SHA="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745"
FLORA_SHA="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"
DOPI_SHA="b95d3084399ed2c0978af2485801108a27d14e843b9fed2611934ea582b745c8"


def read_anonymous(path:Path)->dict:
    if path.suffix.lower()==".zip":
        with ZipFile(path) as z:
            names=[n for n in z.namelist() if n.endswith(".json")]
            if len(names)!=1:raise ValueError("anonymous original artifact requires one JSON")
            return json.loads(z.read(names[0]))
    return json.loads(path.read_text())


def _source_receipts(gift:dict,dopi:dict)->None:
    if not (gift.get("n_full_original_botanical_islands")==4505 and
            gift.get("original_Direct_species_trait_ledger_sha256")==DIRECT_SHA and
            gift.get("original_botanical_flora_sha256")==FLORA_SHA and
            gift.get("n_real_two_trait_pairs")==11):
        raise ValueError("source-original GIFT 4505 or Direct ledger mismatch")
    if not (dopi.get("n_original_GIFT_botanical_islands")==4505 and
            dopi.get("original_DoPI_v2_csv_sha256")==DOPI_SHA and
            dopi.get("n_original_DOPI_source_rows")==676900 and
            dopi.get("n_motif_profiles")==27):
        raise ValueError("DoPI v2 source receipt / original GIFT population mismatch")


def exact_profile_label_benchmark(a,b)->dict:
    """EXACT 8! label shuffle; NOT a valid ecological p-value.

    The true flower predicates overlap and are NONexchangeable botanical
    characters. This is a diagnostic of cherry-picked 3/3 signs only.
    """
    x=np.asarray(a,float);y=np.asarray(b,float)
    if len(x)!=8 or len(y)!=8 or (not np.isfinite(x+y).all()
                                  or (x==0).any() or (y==0).any()):
        raise ValueError("exact source-qualified eight nonzero motifs required")
    nx=np.argsort(np.argsort(x)).astype(float)
    ny=np.argsort(np.argsort(y)).astype(float)
    nx-=nx.mean();ny-=ny.mean()
    denom=float(np.sqrt((nx@nx)*(ny@ny)))
    base_rho=float(nx@ny/denom)
    base_match=int(np.sum(np.sign(x)==np.sign(y)))
    n=0;s=0;r=0
    for k in permutations(range(8)):
        z=np.array(k,int)
        n+=1
        s+=int(np.sum(np.sign(x)==np.sign(y[z]))>=base_match)
        r+=int(abs(float(nx@ny[z]/denom))>=abs(base_rho)-1e-12)
    if n!=40320:raise AssertionError("8! permutations incomplete")
    return dict(matched_signs=base_match,
                n_exact_label_permutations=n,
                sign_tail_fraction=s/n,
                rank_correlation=base_rho,
                abs_rank_tail_fraction=r/n,
                overlapping_profiles_nonexchangeable_not_ecological_p=True)


def audit(gift:dict,dopi:dict)->dict:
    _source_receipts(gift,dopi)
    G={}
    for q in gift["panels"]:
        key=(q["region"],q["plant_scope"],q["documented_joint_sector"])
        if key in G:raise ValueError("duplicated original GIFT panel")
        G[key]=q
    D={}
    for q in dopi["panels"]:
        if (q["guild_a"],q["guild_b"],q["evidence_scope"])!=(
            "butterfly_family","Bombus","positive_visitor_or_pollen_123"
        ):continue
        if q["motif"] in D:raise ValueError("duplicate source DoPI panel")
        D[q["motif"]]=q
    records=[];summaries=[]
    for region in ("tropical","northern_midlatitude",
                   "southern_extratropical","northern_high_latitude"):
      for scope in SCOPES:
       for mode in MODES:
        cells=[]
        for motif in MOTIFS:
            g=G.get((region,scope,motif))
            d=D.get(motif)
            if not g or not d:continue
            p=g["all__spatial_block"];v=d["all_study"]
            if (p["status"]!="DESCRIPTIVE_REAL_SAME_SPECIES_PAIR_ONLY" or
                v["status"]!="DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION"):
                continue
            name="delta_"+mode+"_pair_prevalence"
            x=float(p[name+"_block_equal"])
            y=float(v["mean_guildA_minus_guildB"])
            ci_g=p[name+"_bootstrap95"];ci_d=v["bootstrap95"]
            if not np.isfinite(x+y):raise ValueError("supported source lacks effect")
            shared=d["shared_genus_study"]
            z=dict(region=region,scope=scope,mode=mode,motif=motif,
                 class_name="structural" if motif in STRUCTURAL else "colour_form",
                 gift_relative_difference=x,
                 gift_bootstrap95=ci_g,
                 n_gift_islands=p["n_islands"],
                 dopi_butterfly_minus_bombus=y,
                 dopi_bootstrap95=ci_d,
                 n_dopi_studies=d["n_within_publication_trait_scored_contrasts"],
                 n_dopi_studies_with_positive_motif=d["n_studies_with_positive_motif_either"],
                 matched_genus_dopi_status=shared["status"],
                 matched_genus_dopi_difference=(
                   shared["mean_guildA_minus_guildB"] if
                   shared["status"]=="DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION"
                   else None),
                 sign_match=bool(np.sign(x)==np.sign(y)),
                 both_bootstrap95_exclude_zero=bool(
                    ci_g[0]*ci_g[1]>0 and ci_d[0]*ci_d[1]>0))
            records.append(z);cells.append(z)
        st=[z for z in cells if z["class_name"]=="structural"]
        co=[z for z in cells if z["class_name"]=="colour_form"]
        summary=dict(region=region,scope=scope,mode=mode,
          n_jointly_supported_motifs=len(cells),
          structural_sign_matches=sum(z["sign_match"] for z in st),
          structural_total=len(st),
          colour_form_sign_matches=sum(z["sign_match"] for z in co),
          colour_form_total=len(co),
          total_sign_matches=sum(z["sign_match"] for z in cells),
          n_both_95pct_intervals_excluding_zero=sum(z["both_bootstrap95_exclude_zero"] for z in cells))
        if len(cells)==8:
            a=[z["gift_relative_difference"] for z in cells]
            b=[z["dopi_butterfly_minus_bombus"] for z in cells]
            summary["descriptive_profile_label_shuffle"]=exact_profile_label_benchmark(a,b)
            summary["Pearson_across_nonindependent_motifs"]=float(pearsonr(a,b).statistic)
            summary["Spearman_across_nonindependent_motifs"]=float(spearmanr(a,b).statistic)
        summaries.append(summary)
    if not records:raise ValueError("no original anonymous supported source motif pairs")
    return dict(status="GIFT4505_VS_DOPI2026_EIGHT_MOTIF_SOURCE_CONCORDANCE_AUDIT",
      n_global_GIFT_botanical_island_floras=4505,
      n_compared_anonymous_motif_scenario_rows=len(records),
      motif_rows=records,summaries=summaries,
      limitations=dict(European_DoPI_not_local_tropical_island_guild=True,
        source_and_visitor_frequency_differences_not_same_estimand=True,
        accepted_species_direct_flower_traits_not_local_evolution=True,
        overlapping_motifs_not_independent=True,
        positive_visitor_links_not_effective_pollination=True,
        genus_and_trait_documentation_confounding=True,
        selected_3of3_directional_alignment_not_confirmatory=True,
        no_pollinator_shift_established=True))


def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--gift-aggregate",type=Path,required=True)
    p.add_argument("--dopi-aggregate",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=audit(read_anonymous(args.gift_aggregate),
                 read_anonymous(args.dopi_aggregate))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    q=next(z for z in result["summaries"] if z["region"]=="tropical"
        and z["scope"]=="all_recorded" and z["mode"]=="geo5")
    print(json.dumps({"status":result["status"],"tropical":q},sort_keys=True))

if __name__=="__main__":main()
