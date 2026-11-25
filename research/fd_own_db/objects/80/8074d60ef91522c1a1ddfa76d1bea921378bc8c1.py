"""Rebuild four fd manuscript figures from SHA-pinned, anonymized ORIGINAL island Actions artifacts.

No original GIFT accepted species, island IDs, source entities, locations or
pollinator service observations enter or leave this standalone plotter.
All numbers are observational, not a causal island-syndrome test.
"""
from __future__ import annotations
import argparse,hashlib,json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT={
"composition":("anonymous_island_functional_beta.json","c9b5cceb6930eb667557444edc6f21966ac30564086953ef3d4ff7d161c78494","ISLAND_FUNCTIONAL_COMPOSITION_BETA_DESCRIPTIVE"),
"rao":("anonymous_island_floral_beta.json","7d321333eca23e0374f22d3468785ac2a184e60c1810345c8fd2a53911e5dee2","GLOBAL_ISLAND_FLORAL_COMPOSITION_BETA_NOT_ADAPTATION"),
"pairs":("anonymous_original_true_joint_syndrome_access.json","478e89e0dbdde12c967c1611c70b4b6fc18f500eb1fc5b8832e53da0f7cbaf54","ORIGINAL_4505_GIFT_REAL_PAIRED_SYNDROME_ACCESS_SOURCE"),
"genus":("anonymous_joint_syndrome_genus_4505.json","56b01a4fae6264a313eef8fb728aa6cb4efb5dd9bb22f6084f02ccc304bd8fc8","ORIGINAL_4505_REAL_SYNDROME_PAIR_SOURCE_GENUS_PARTITION")}
REGION=["northern_midlatitude","tropical","southern_extratropical","northern_high_latitude"]
SECTOR=["actinomorphic_deep","zygomorphic_deep","tubular_bell_deep"]
WORDS={"flower_primary_color":"Colour","floral_form":"Form","floral_symmetry":"Symmetry",
       "tube_depth_class":"Tube depth","flower_primary_color__floral_form":"Colour × form",
       "flower_primary_color__floral_symmetry":"Colour × symmetry",
       "floral_form__floral_symmetry":"Form × symmetry",
       "actinomorphic_deep":"Radial × deep","zygomorphic_deep":"Bilateral × deep",
       "tubular_bell_deep":"Tubular/bell × deep"}
RWORDS=["N. midlatitude","Tropical","S. extratropical","N. high latitude"]
NAVY="#24486b";TEAL="#168881";ORANGE="#cf7743";PALE="#75aab0"

def one(rows,**kw):
    q=[v for v in rows if all(v.get(k)==value for k,value in kw.items())]
    if len(q)!=1:raise ValueError("Original anonymous panel was not uniquely identifiable: "+str(kw))
    return q[0]

def load(root):
    out={}; manifest={}
    for name,(fn,expected,status) in INPUT.items():
        raw=(root/fn).read_bytes()
        h=hashlib.sha256(raw).hexdigest()
        if h!=expected:raise ValueError("Original source JSON SHA changed for "+name)
        obj=json.loads(raw)
        if obj.get("status")!=status:raise ValueError("Original source type changed")
        n=(obj.get("n_original_islands") or obj.get("n_original_flora_islands") or
           obj.get("n_full_original_botanical_islands") or
           obj.get("n_original_GIFT_botanical_islands"))
        if n!=4505:raise ValueError("Original full GIFT 4505-island study population changed")
        out[name]=obj;manifest[name]={"sha256":h,"n_original_flora_islands":n}
    if len(out["composition"]["endpoint_coverage"])!=26 or len(out["rao"]["endpoints"])!=26:
        raise ValueError("Original 13×2 finite trait support axes changed")
    if out["pairs"]["n_real_two_trait_pairs"]!=11 or len(out["genus"]["panels"])!=45:
        raise ValueError("Original real plant pair/genus source schema changed")
    return out,manifest

def tidy(ax):
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="x",alpha=.16)
    ax.set_axisbelow(True)
    ax.tick_params(length=0,labelsize=9)

def save(fig,folder,name):
    fig.savefig(folder/(name+".png"),dpi=175,bbox_inches="tight",facecolor="white")
    fig.savefig(folder/(name+".svg"),bbox_inches="tight",facecolor="white")
    plt.close(fig)

def coverage_and_rao(data,out):
    traits=["flower_primary_color","floral_form","floral_symmetry","tube_depth_class",
            "flower_primary_color__floral_form","flower_primary_color__floral_symmetry",
            "floral_form__floral_symmetry"]
    count=[];q=[];records=[]
    for s in traits:
        aa=[];bb=[]
        for scope in ("all_recorded","original_native"):
            x=one(data["composition"]["endpoint_coverage"],scope=scope,endpoint=s)
            y=one(data["rao"]["endpoints"],scope=scope,endpoint=s)
            if x["n_n5_islands"]!=y["n_eligible_islands"]:
                raise ValueError("Original independent coverage and Rao Q not cohort equivalent")
            aa.append(x["n_n5_islands"]);bb.append(y["mean_raw_finite5_rao"])
            records.append({"endpoint":s,"scope":scope,"n_islands":aa[-1],"finite5_Rao_Q":bb[-1]})
        count.append(aa);q.append(bb)
    C=np.array(count);Q=np.array(q);y=np.arange(len(traits))
    fig,(a,b)=plt.subplots(1,2,figsize=(13,5.5),gridspec_kw={"width_ratios":[1.35,1]})
    a.barh(y+.2,C[:,0],height=.36,color=NAVY,label="All recorded")
    a.barh(y-.2,C[:,1],height=.36,color=TEAL,label="Source-known native")
    for i in y:
        a.text(C[i,0]+25,i+.2,str(C[i,0]),va="center",fontsize=8)
        a.text(C[i,1]+25,i-.2,str(C[i,1]),va="center",fontsize=8,color=TEAL)
    a.set_yticks(y,[WORDS[v] for v in traits]);a.invert_yaxis()
    a.set_xlim(0,4505);a.set_xlabel("Original islands with ≥5 directly documented species")
    a.set_title("A. Evidence by trait and provenance",loc="left",fontweight="bold")
    a.legend(frameon=False,fontsize=8,loc="lower right")
    b.plot(Q[:,0],y,"o-",color=NAVY,lw=1.5,label="All recorded")
    b.plot(Q[:,1],y,"s-",color=TEAL,lw=1.5,label="Source-known native")
    b.set_yticks(y,[""]*len(y));b.invert_yaxis();b.set_xlim(0,1)
    b.set_title("B. Rao Q is NOT a hypervolume",loc="left",fontweight="bold")
    b.set_xlabel("Observed finite-five-species categorical Rao Q")
    fig.suptitle("Figure 1  |  The study frame is 4,505 original GIFT botanical islands",
                 x=.01,y=1.04,ha="left",fontweight="bold",fontsize=13)
    for ax in (a,b):tidy(ax)
    fig.tight_layout()
    save(fig,out,"fig1_trait_coverage_and_RaoQ")
    return records

def source_pair_region(data,out):
    rows=[]; fig,axes=plt.subplots(1,3,figsize=(14.5,5),sharey=True)
    for j,sector in enumerate(SECTOR):
        ax=axes[j]; y=np.arange(4)
        for mode,clr,marker,delta in (("geo5",NAVY,"o",-.1),
                                      ("climate10",ORANGE,"s",.1)):
            means=[];errs=[]
            for region in REGION:
                p=one(data["pairs"]["panels"],region=region,plant_scope="all_recorded",
                      documented_joint_sector=sector)
                g=one(data["genus"]["panels"],region=region,scope="all_recorded",
                      sector=sector)
                if p["n_both_source_modes_and_joint_carriers"]!=g["n_same_original_islands_both_mainland_sources"]:
                    raise ValueError("Original source genus and true paired plant samples differ")
                x=p["all__spatial_block"]
                if x["status"]!="DESCRIPTIVE_REAL_SAME_SPECIES_PAIR_ONLY":
                    raise ValueError("Unidentifiable original source pair must not be plotted")
                m=x[f"delta_{mode}_pair_prevalence_block_equal"]*100
                low,high=[a*100 for a in x[f"delta_{mode}_pair_prevalence_bootstrap95"]]
                n=p["n_both_source_modes_and_joint_carriers"]
                if abs(m-100*g["all__spatial_block"][mode+"__observed_difference__equal_group_mean"])>1e-9:
                    raise ValueError("Different source original calculations disagree")
                rows.append(dict(region=region,sector=sector,scope="all_recorded",
                                 original_source=mode,n_islands=n,mean_delta_pp=m,
                                 lower95_pp=low,upper95_pp=high,n_original_blocks=x["n_clusters"]))
                means.append(m);errs.append((m-low,high-m))
            err=np.array(errs).T
            ax.errorbar(means,y+delta,xerr=err,fmt=marker,color=clr,capsize=3,
                        label=("Geographic 5" if mode=="geo5" else "Climate 10"),lw=1.5)
        ax.axvline(0,color="#aab6c1",ls="--")
        ax.set_yticks(y,RWORDS);ax.invert_yaxis()
        ax.set_xlim(-13,22)
        ax.set_title("ABC"[j]+". "+WORDS[sector],loc="left",fontweight="bold")
        ax.set_xlabel("Original island − candidate source (pp)")
        tidy(ax)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.51,1.03),
               frameon=False,ncol=2,fontsize=8.7)
    fig.suptitle("Figure 2  |  Real floral architecture differs across regions",
                 x=.01,y=1.11,ha="left",fontweight="bold",fontsize=13)
    fig.text(.01,-.03,"95% 10° geographical-block resampling; two source definitions share original focal island plants. Exploratory, not a guild-effect test.",fontsize=8,color="#555")
    fig.tight_layout()
    save(fig,out,"fig2_regional_true_flower_pair_directions")
    return rows

def genus_signed(data,out):
    terms=["shared_genus_composition","shared_genus_within_species_states",
           "native_exclusive_genus","source_exclusive_genus_negative"]
    display=["Shared genus relative weights","Shared genus different-species states",
             "Island-only documented genera","Source-only candidate genera"]
    rows=[];fig,axes=plt.subplots(1,3,figsize=(14.7,5.5),sharey=True)
    for j,sector in enumerate(SECTOR):
        p=one(data["genus"]["panels"],region="tropical",scope="all_recorded",sector=sector)
        cell=p["all__spatial_block"]
        if cell["status"]!="EXACT_REAL_PAIRED_FLORA_GENUS_ACCOUNTING_NOT_CAUSAL":
            raise ValueError("Missing original signed genus accounting")
        values=[]
        for term in terms:
            v=100*cell["geo5__"+term+"__equal_group_mean"];values.append(v)
            rows.append(dict(sector=sector,term=term,source="geo5",signed_pp=v,
                             n_islands=p["n_same_original_islands_both_mainland_sources"]))
        total=100*cell["geo5__observed_difference__equal_group_mean"]
        lo,hi=[100*x for x in cell["geo5__observed_difference__bootstrap95"]]
        if abs(sum(values)-total)>1e-8:
            raise ValueError("Four original signed genus terms do not sum exactly")
        rows.append(dict(sector=sector,term="observed_total",source="geo5",
                         signed_pp=total,n_islands=p["n_same_original_islands_both_mainland_sources"]))
        ax=axes[j];ax.barh(range(4),values,color=[NAVY,TEAL,PALE,ORANGE],height=.6)
        ax.axvline(0,color="#bbc4cc",ls="--")
        ax.set_yticks(range(4),display);ax.invert_yaxis();ax.set_xlim(-29,29)
        for y,v in enumerate(values):
            ax.text(v+.5 if v>=0 else v-.5,y,f"{v:+.2f}",
                    ha="left" if v>=0 else "right",va="center",fontsize=8.7)
        ax.set_xlabel("Exact signed island − source contribution (pp)")
        ax.set_title("ABC"[j]+". "+WORDS[sector]+
                     f"\nTotal {total:+.2f} pp  [95% {lo:+.1f}, {hi:+.1f}]",
                     loc="left",fontsize=9.7,fontweight="bold")
        tidy(ax)
    fig.suptitle("Figure 3  |  Tropical observed genus-composition accounting",
                 x=.01,y=1.09,ha="left",fontweight="bold",fontsize=13)
    fig.text(.01,-.03,"Source: closest geographic 5. Four signed values cancel and are NOT fractions of evolution, selection or causal variance.",fontsize=8,color="#555")
    fig.tight_layout()
    save(fig,out,"fig3_tropical_signed_genus_accounting")
    return rows

def beta_profiles(data,out):
    endpoint="flower_primary_color__floral_form"
    summaries=data["composition"]["pair_summaries"]
    a=one(summaries,scope="all_recorded",endpoint=endpoint,
          region="GLOBAL_ALL",geology_pair="all",subset="all")
    b=one(summaries,scope="all_recorded",endpoint=endpoint,
          region="GLOBAL_ALL",geology_pair="all",subset="no_shared_full_flora_species")
    if (a["n_pairs"],b["n_pairs"])!=(2699326,842058):
        raise ValueError("Original 2.7m / 842k real island-pair cohorts changed")
    fig,axes=plt.subplots(1,2,figsize=(13.4,5.5),gridspec_kw={"width_ratios":[1,1.22]})
    x=np.arange(2)
    func=[a["mean_functional_overlap"],b["mean_functional_overlap"]]
    tax=[a["mean_taxonomic_similarity_full"],b["mean_taxonomic_similarity_full"]]
    axes[0].barh(x+.18,func,height=.33,color=NAVY,label="Category-frequency overlap")
    axes[0].barh(x-.18,tax,height=.33,color=ORANGE,label="Species Sørensen–Dice")
    axes[0].set_yticks(x,["All island pairs","No shared accepted plant species"])
    axes[0].invert_yaxis();axes[0].set_xlim(0,.57)
    axes[0].set_title("A. Trait categories can recur across species",loc="left",fontweight="bold")
    axes[0].set_xlabel("Pairwise overlap: distinct category and species metrics")
    axes[0].legend(loc="upper left",bbox_to_anchor=(0,-.18),frameon=False,fontsize=8)
    selected=[("white","tubular"),("white","funnel_trumpet"),
              ("white","open_radial"),("yellow_orange","open_radial")]
    colors=[]
    for region in ("tropical","northern_midlatitude"):
        p=one(data["composition"]["regional_recorded_trait_composition"],
              region=region,scope="all_recorded",endpoint=endpoint,geological_origin="all")
        f={tuple(i["state"]):i["island_equal_mean_recorded_species_fraction"] for i in p["categories"]}
        colors.append([100*f.get(k,0.) for k in selected])
    y=np.arange(len(selected))
    axes[1].barh(y-.18,colors[0],height=.33,color=TEAL,label="Tropical; n=586")
    axes[1].barh(y+.18,colors[1],height=.33,color=NAVY,label="N. midlatitude; n=1,464")
    axes[1].set_yticks(y,["White × tubular","White × funnel/trumpet",
                          "White × open/radial","Yellow/orange × open/radial"])
    axes[1].invert_yaxis()
    axes[1].set_title("B. Frequencies differ among regions",loc="left",fontweight="bold")
    axes[1].set_xlabel("Recorded same-plant colour×form species share (%)")
    axes[1].legend(loc="upper left",bbox_to_anchor=(0,-.18),ncol=2,frameon=False,fontsize=8)
    for ax in axes:tidy(ax)
    fig.suptitle("Figure 4  |  Functional sharing coexists with regional structure",
                 x=.01,y=1.06,ha="left",fontweight="bold",fontsize=13)
    fig.text(.01,-.04,"Species similarity and functional frequency overlap are different metrics. No named pollinator, dispersal or adaptation estimated.",fontsize=8,color="#555")
    fig.tight_layout()
    save(fig,out,"fig4_functional_vs_taxonomic_and_regional_profiles")
    return [{"cohort":n,"n_original_pairs":v["n_pairs"],
             "mean_functional_overlap":v["mean_functional_overlap"],
             "mean_species_sorensen_dice":v["mean_taxonomic_similarity_full"],
             "recorded_profile_turnover":v["mean_category_turnover"],
             "recorded_profile_nestedness":v["mean_category_nestedness"]}
            for n,v in (("all",a),("no_shared_full_flora_species",b))]

def write_table(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    data,manifest=load(a.input)
    figs=[
        ("figure1_trait_coverage_rao.csv",coverage_and_rao),
        ("figure2_regional_pair_source.csv",source_pair_region),
        ("figure3_tropical_signed_genus.csv",genus_signed),
        ("figure4_similarity.csv",beta_profiles)
    ]
    for name,fn in figs:write_table(a.output/name,fn(data,a.output))
    receipt=dict(status="ORIGINAL_GIFT_4505_SOURCE_PINNED_MANUSCRIPT_FIGURES",
                 n_original_islands=4505,n_original_figures=4,
                 original_source_sha256=manifest,
                 limits=dict(native_status_separately_recorded=True,
                             no_species_imputation=True,
                             no_17D_KDE_as_functional_diversity=True,
                             no_pollinator_guild_causality=True,
                             signed_genus_terms_not_causal_variance=True,
                             original_928_complete_case_not_main_population=True))
    (a.output/"anonymous_figure_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":receipt["status"],"n_original_figures":4},sort_keys=True))

if __name__=="__main__": main()
