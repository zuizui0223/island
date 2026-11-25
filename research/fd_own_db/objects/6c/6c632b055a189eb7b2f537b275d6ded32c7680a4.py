"""Real SAME-species floral trait conjunctions: nonexclusive access/signal proxies.

A conjunction may be consistent with more than one visitor functional guild;
none is a pollinator classification, species-specific effective interaction,
historical pollinator shift, local island flower evolution or fitness estimate.
"""
from __future__ import annotations

from hashlib import sha256
import numpy as np
import pandas as pd

# Outcome-blind, explicitly semantic exact species-Direct pair construction.
# No multiplication of plant-trait marginal frequencies across distinct taxa.
JOINT_SPECS = {
    "actinomorphic_deep": (
        "floral_symmetry", frozenset({"actinomorphic"}),
        "tube_depth_class", frozenset({"deep"})),
    "zygomorphic_deep": (
        "floral_symmetry", frozenset({"zygomorphic"}),
        "tube_depth_class", frozenset({"deep"})),
    "tubular_bell_deep": (
        "floral_form", frozenset({"tubular", "bell_campanulate"}),
        "tube_depth_class", frozenset({"deep"})),
    "actinomorphic_tubular_bell": (
        "floral_symmetry", frozenset({"actinomorphic"}),
        "floral_form", frozenset({"tubular", "bell_campanulate"})),
    "zygomorphic_tubular_bell": (
        "floral_symmetry", frozenset({"zygomorphic"}),
        "floral_form", frozenset({"tubular", "bell_campanulate"})),
    # The original accepted Direct colour ledger does NOT have standalone
    # red/orange/cream labels. It has red_pink, yellow_orange and white.
    # This factual source-ontology repair follows an initial pilot in which
    # literal red/orange led to all-zero targets. Colour analyses below are
    # exploratory and cannot discriminate bird from butterfly/moth guilds.
    "red_pink_deep": (
        "flower_primary_color", frozenset({"red_pink"}),
        "tube_depth_class", frozenset({"deep"})),
    "yellow_orange_deep": (
        "flower_primary_color", frozenset({"yellow_orange"}),
        "tube_depth_class", frozenset({"deep"})),
    "white_deep": (
        "flower_primary_color", frozenset({"white"}),
        "tube_depth_class", frozenset({"deep"})),
    "red_pink_tubular_bell": (
        "flower_primary_color", frozenset({"red_pink"}),
        "floral_form", frozenset({"tubular", "bell_campanulate"})),
    "yellow_orange_tubular_bell": (
        "flower_primary_color", frozenset({"yellow_orange"}),
        "floral_form", frozenset({"tubular", "bell_campanulate"})),
    "white_tubular_bell": (
        "flower_primary_color", frozenset({"white"}),
        "floral_form", frozenset({"tubular", "bell_campanulate"})),
}
PRIORITIZED_MORPHOLOGY = (
    "actinomorphic_deep", "zygomorphic_deep", "tubular_bell_deep",
    "actinomorphic_tubular_bell", "zygomorphic_tubular_bell"
)
MODES=("geo_k5","geo50_climate10")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
MIN_ISLAND=5
MIN_SOURCE=10
MIN_REGION_ISLANDS=30
MIN_REGION_BLOCKS=10


def make_joint_species(states: dict[str, dict[str,str]]) -> dict:
    if not {x for spec in JOINT_SPECS.values()
            for x in (spec[0],spec[2])} <= set(states):
        raise ValueError("original Direct pair axis omitted")
    output={}
    for name,(ta,positive_a,tb,positive_b) in JOINT_SPECS.items():
        # Intersection is REQUIRED: a plant with only one axis is NOT coded
        # negative for the missing other trait.
        species=set(states[ta]) & set(states[tb])
        a=frozenset(sp for sp in species if states[ta][sp] in positive_a)
        b=frozenset(sp for sp in species if states[tb][sp] in positive_b)
        # Discrete inclusion means actual same accepted species has BOTH states.
        output[name]=dict(
            known=frozenset(species), A=a,B=b, joint=a&b,
            trait_axes=(ta,tb), positive_a=tuple(sorted(positive_a)),
            positive_b=tuple(sorted(positive_b)),
        )
    return output


def summarize_one_roster(roster: set[str], profile: dict) -> dict:
    known=roster & profile["known"]
    n=len(known)
    if not n:
        return dict(n=0,p_joint=np.nan, p_A=np.nan,p_B=np.nan,
                    coupling=np.nan,n_joint=0,n_A=0,n_B=0)
    a=len(known & profile["A"])
    b=len(known & profile["B"])
    both=len(known & profile["joint"])
    pa=a/n
    pb=b/n
    frac=both/n
    return dict(n=n,p_joint=frac,p_A=pa,p_B=pb,
                coupling=float(frac-pa*pb),
                n_joint=both,n_A=a,n_B=b)


def finite_joint_null_counts(n: int, a: int, b: int) -> tuple[float,float]:
    """Exact fixed margin random-pair expectation and finite variance of k."""
    if not all(isinstance(z,(int,np.integer)) for z in (n,a,b)) or (
        n < 0 or min(a,b)<0 or max(a,b)>n
    ):
        raise ValueError("invalid finite binary flower pair margins")
    if n==0:
        return 0.,0.
    exp=a*b/n
    var=(a*(b/n)*(1-b/n)*((n-a)/(n-1))) if n>1 else 0.
    return float(exp),float(var)


def same_island_pair_source_frame(rows: pd.DataFrame,
                                  *, population: int=4505,
                                  bootstrap: int=999)->dict:
    need={
        "island_id","regime","scope","sector","spatial_block","first_source",
        "island_total","island_direct","geo5_total","geo5_direct",
        "climate10_total","climate10_direct",
        "island_pair_prevalence","geo5_pair_prevalence",
        "climate10_pair_prevalence",
        "island_pair_coupling","geo5_pair_coupling",
        "climate10_pair_coupling",
    }
    if not need <= set(rows):
        raise ValueError(f"missing real pair-source input {need-set(rows)}")
    if rows.duplicated(["island_id","scope","sector"]).any():
        raise ValueError("duplicate island×scope×real-joint-sector")
    if rows.island_id.nunique()!=population:
        raise ValueError("global original botanical island population changed")
    if len(rows)!=population*len(SCOPES)*len(JOINT_SPECS):
        raise ValueError("real joint phenotype/scopes were reduced")
    if not rows.scope.isin(SCOPES).all() or not rows.sector.isin(JOINT_SPECS).all():
        raise ValueError("unregistered source scope or paired phenotype")
    if not rows.regime.isin(REGIMES).all():
        raise ValueError("unregistered real GSHHG geographical regime")
    for z in ("island_total","island_direct","geo5_total","geo5_direct",
              "climate10_total","climate10_direct"):
        q=pd.to_numeric(rows[z],errors="coerce")
        if q.isna().any() or (q<0).any() or (q%1).ne(0).any():
            raise ValueError("invalid original Direct pair carrier counts")
    for prefix in ("island","geo5","climate10"):
        if (rows[f"{prefix}_direct"]>rows[f"{prefix}_total"]).any():
            raise ValueError("Direct pair carrier exceeds original source members")
        for metric in ("pair_prevalence","pair_coupling"):
            q=rows[f"{prefix}_{metric}"]
            if (q.isna()!=rows[f"{prefix}_direct"].eq(0)).any():
                raise ValueError("real pair missingness was imputed")
        if ((rows[f"{prefix}_pair_prevalence"].dropna()<0)|
            (rows[f"{prefix}_pair_prevalence"].dropna()>1)).any():
            raise ValueError("observed paired trait proportion out of range")
    if bootstrap<199:
        raise ValueError("source block uncertainty must have >=199 draws")

    report=[]
    for region in REGIMES:
        for scope in SCOPES:
            for sector in JOINT_SPECS:
                pool=rows.loc[
                    rows.regime.eq(region) & rows.scope.eq(scope)
                    & rows.sector.eq(sector)
                ].copy()
                q=pool.loc[
                    pool.island_direct.ge(MIN_ISLAND)
                    & pool.geo5_direct.ge(MIN_SOURCE)
                    & pool.climate10_direct.ge(MIN_SOURCE)
                    & pool.island_pair_prevalence.notna()
                    & pool.geo5_pair_prevalence.notna()
                    & pool.climate10_pair_prevalence.notna()
                ].copy()
                for metric in ("pair_prevalence","pair_coupling"):
                    q["delta_geo5_"+metric]=(
                        q["island_"+metric]-q["geo5_"+metric])
                    q["delta_climate10_"+metric]=(
                        q["island_"+metric]-q["climate10_"+metric])
                if len(q):
                    frac_island=q.island_direct/q.island_total.replace(0,np.nan)
                    frac_geo=q.geo5_direct/q.geo5_total.replace(0,np.nan)
                    frac_clim=q.climate10_direct/q.climate10_total.replace(0,np.nan)
                    balanced=(frac_island.sub(frac_geo).abs().le(.10)
                              & frac_island.sub(frac_clim).abs().le(.10))
                else:
                    balanced=pd.Series(False,index=q.index)
                result=dict(
                    region=region, plant_scope=scope,
                    documented_joint_sector=sector,
                    prioritized_morphological_conjunction_exploratory=sector in PRIORITIZED_MORPHOLOGY,
                    n_original_botanical_islands=int(len(pool)),
                    n_both_source_modes_and_joint_carriers=int(len(q)),
                    n_original_documentation_balanced_le10pp=int(balanced.sum()),
                    original_joint_axis_names=list(JOINT_SPECS[sector][::2][:2]),
                )
                for tag, subset in (("all",q),
                                    ("doc_gap_le10pp",q.loc[balanced])):
                    for group in ("spatial_block","first_source"):
                        z=subset.dropna(subset=[group]).groupby(group,sort=True)
                        G=z.ngroups
                        obj=dict(n_islands=int(len(subset)),n_clusters=int(G))
                        if len(subset)<MIN_REGION_ISLANDS or G<MIN_REGION_BLOCKS:
                            obj["status"]="NOT_IDENTIFIABLE_PAIR_SOURCE_SUPPORT"
                        else:
                            columns=[
                                f"delta_{mode}_{metric}"
                                for mode in ("geo5","climate10")
                                for metric in ("pair_prevalence","pair_coupling")
                            ]
                            vals=z[columns].mean().to_numpy(float)
                            mean=vals.mean(axis=0)
                            seed=int.from_bytes(sha256(
                                f"{region}|{scope}|{sector}|{tag}|{group}".encode()
                            ).digest()[:8],"big")
                            rng=np.random.default_rng(seed)
                            boot=vals[rng.integers(G,size=(bootstrap,G))].mean(1)
                            for i,col in enumerate(columns):
                                obj[col+"_block_equal"]=float(mean[i])
                                obj[col+"_bootstrap95"]=[
                                    float(x) for x in np.quantile(
                                        boot[:,i],[.025,.975])
                                ]
                            obj["status"]="DESCRIPTIVE_REAL_SAME_SPECIES_PAIR_ONLY"
                        result[tag+"__"+group]=obj
                report.append(result)
    return {
        "status":"ORIGINAL_4505_GIFT_REAL_PAIRED_SYNDROME_ACCESS_SOURCE",
        "n_full_original_botanical_islands":population,
        "n_real_two_trait_pairs":len(JOINT_SPECS),
        "n_original_island_scope_joint_sector_rows":int(len(rows)),
        "n_source_mode_provenance_region_panels":len(report),
        "panels":report,
        "limits":{
            "marginal_flower_traits_never_multiplied_across_different_plants":True,
            "coupling_uses_exact_same_species_fixed_marginal_reference":True,
            "no_guild_identity_predicted_from_morphology":True,
            "not_known_historical_pollinator_shift_or_fitness":True,
            "strict_native_is_selected_nonrandom_status_subset":True,
            "entire_4505_flora_frame_retained":True,
            "morphological_conjunction_not_pollinator_syndrome_classification":True,
            "color_bins_originally_red_pink_yellow_orange_white_not_separate_red_orange":True,
            "original_first_pilot_zero_color_was_invalid_ontology_not_null_evidence":True,
            "color_ontology_mapping_corrected_after_pilot_and_results_exploratory":True,
        },
    }
