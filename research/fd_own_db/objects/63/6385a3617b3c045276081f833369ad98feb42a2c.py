"""Quantify genuine, species-Direct multi-axis floral signal/access motifs.

Every motif is an observed same-accepted-species category conjunction, NEVER a
named effective pollinator syndrome. Colour bins in the original Direct ledger
include red_pink and yellow_orange, not independent red, pink or orange labels.
Models compare potential native GIFT mainland source-only plant lists with
documented island flora; they do not observe arrival, visitors, pollen transfer
or island population adaptation.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256

import numpy as np
import pandas as pd

REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
MODES=("geo5","climate10")
# Anatomical/visual combinations, not discrete assignments to bird, bee or moth.
MOTIFS={
 "red_pink_radial_deep":(
  ("flower_primary_color",("red_pink",)),
  ("floral_symmetry",("actinomorphic",)),
  ("tube_depth_class",("deep",))),
 "yellow_orange_radial_deep":(
  ("flower_primary_color",("yellow_orange",)),
  ("floral_symmetry",("actinomorphic",)),
  ("tube_depth_class",("deep",))),
 "white_radial_deep":(
  ("flower_primary_color",("white",)),
  ("floral_symmetry",("actinomorphic",)),
  ("tube_depth_class",("deep",))),
 "red_pink_tubular_deep":(
  ("flower_primary_color",("red_pink",)),
  ("floral_form",("tubular","bell_campanulate")),
  ("tube_depth_class",("deep",))),
 "yellow_orange_tubular_deep":(
  ("flower_primary_color",("yellow_orange",)),
  ("floral_form",("tubular","bell_campanulate")),
  ("tube_depth_class",("deep",))),
 "white_tubular_deep":(
  ("flower_primary_color",("white",)),
  ("floral_form",("tubular","bell_campanulate")),
  ("tube_depth_class",("deep",))),
 "blue_purple_zygomorphic_deep":(
  ("flower_primary_color",("blue_purple",)),
  ("floral_symmetry",("zygomorphic",)),
  ("tube_depth_class",("deep",))),
 "yellow_orange_zygomorphic_deep":(
  ("flower_primary_color",("yellow_orange",)),
  ("floral_symmetry",("zygomorphic",)),
  ("tube_depth_class",("deep",))),
 "open_radial_shallow":(
  ("floral_form",("open_radial",)),
  ("floral_symmetry",("actinomorphic",)),
  ("tube_depth_class",("shallow",))),
 "red_pink_radial_tubular_deep":(
  ("flower_primary_color",("red_pink",)),
  ("floral_symmetry",("actinomorphic",)),
  ("floral_form",("tubular","bell_campanulate")),
  ("tube_depth_class",("deep",))),
 "yellow_orange_radial_tubular_deep":(
  ("flower_primary_color",("yellow_orange",)),
  ("floral_symmetry",("actinomorphic",)),
  ("floral_form",("tubular","bell_campanulate")),
  ("tube_depth_class",("deep",))),
 "white_radial_tubular_deep":(
  ("flower_primary_color",("white",)),
  ("floral_symmetry",("actinomorphic",)),
  ("floral_form",("tubular","bell_campanulate")),
  ("tube_depth_class",("deep",))),
}
PARTS=("shared_genus_source_weight","shared_genus_different_species_state",
       "island_only_recorded_genus","source_only_recorded_genus_negative")
MIN_ISLAND=5
MIN_SOURCE=10
MIN_REGION_ISLANDS=30
MIN_BLOCKS=10


def _genus(species:str)->str:
    bits=str(species).split()
    return bits[0] if len(bits)>1 else "unparsed_name:"+str(species)


def make_motifs(states:dict[str,dict[str,str]])->dict:
    if not set(k for motif in MOTIFS.values() for k,_ in motif)<=set(states):
        raise ValueError("required real species-Direct floral axes absent")
    out={}
    for name,spec in MOTIFS.items():
        known=set.intersection(*(set(states[axis]) for axis,_ in spec))
        positives=tuple(
            frozenset(sp for sp in known if states[axis][sp] in values)
            for axis,values in spec
        )
        joint=frozenset(set.intersection(*map(set,positives)))
        out[name]=dict(known=frozenset(known), joint=joint,
                       axis_positive=positives, n_axes=len(spec),
                       axes=tuple(axis for axis,_ in spec))
    return out


def _stat(roster:set[str], motif:dict)->dict:
    known=roster & motif["known"]
    n=len(known)
    if not n:
        return dict(n=0,positive=0,freq=np.nan,independent=np.nan,
                    coupling=np.nan)
    k=len(known & motif["joint"])
    product=1.
    for s in motif["axis_positive"]:
        product*=len(known & s)/n
    observed=k/n
    return dict(n=n,positive=k,freq=float(observed),
                independent=float(product),
                coupling=float(observed-product))


def _donor_cache(mainland:set[str],motif:dict)->dict:
    known=mainland & motif["known"]
    return dict(n=len(known),
                counts=Counter(map(_genus,known)),
                positives=Counter(map(_genus,known & motif["joint"])),
                total_positive=len(known & motif["joint"]))


def fourpart(island:set[str],mainland:set[str],motif:dict,
             *,source_cache:dict|None=None)->tuple[float,...]:
    """Exact signed four-component source genus identity, NOT causal shares."""
    I=island & motif["known"]
    source_cache=source_cache if source_cache is not None else _donor_cache(mainland,motif)
    if not I or not source_cache["n"]: return (np.nan,)*4
    count_i=Counter(map(_genus,I))
    count_s=source_cache["counts"]
    pos_i=Counter(map(_genus,I & motif["joint"]))
    pos_s=source_cache["positives"]
    ni=len(I);ns=source_cache["n"]
    shared=set(count_i)&set(count_s)
    a=sum((count_i[g]/ni-count_s[g]/ns)*(pos_s[g]/count_s[g])
          for g in shared)
    b=sum((count_i[g]/ni)*(pos_i[g]/count_i[g]-pos_s[g]/count_s[g])
          for g in shared)
    c=sum(pos_i[g]/ni for g in set(count_i)-set(count_s))
    d=-sum(pos_s[g]/ns for g in set(count_s)-set(count_i))
    observed=len(I & motif["joint"])/ni-source_cache["total_positive"]/ns
    if abs(a+b+c+d-observed)>1e-10:
        raise AssertionError("genuine 3/4-axis genus identity failed")
    return tuple(map(float,(a,b,c,d)))


def make_island_motif_rows(*,island_id:str,region:str,spatial_block:str|None,
                           first_source:str|None,
                           scope_rosters:dict[str,set[str]],
                           sources:dict[str,set[str]],
                           motifs:dict)->list[dict]:
    if region not in REGIMES or not set(SCOPES)<=set(scope_rosters) or not set(MODES)<=set(sources):
        raise ValueError("GIFT source region/scope/mode missing")
    full=scope_rosters["all_recorded"]
    if any(not scope_rosters[key]<=full for key in SCOPES):
        raise ValueError("island native status outside botanical flora")
    if any(set(sources[k])&full for k in MODES):
        raise ValueError("source donor plant is on original focal island flora")
    source_genus_full={mode:Counter(map(_genus,sources[mode])) for mode in MODES}
    scope_genus_full={scope:Counter(map(_genus,scope_rosters[scope])) for scope in SCOPES}
    rows=[]
    for name,motif in motifs.items():
        src={mode:_stat(sources[mode],motif) for mode in MODES}
        source_cache={mode:_donor_cache(sources[mode],motif) for mode in MODES}
        for scope in SCOPES:
            roster=scope_rosters[scope]
            island=_stat(roster,motif)
            island_direct_genera=Counter(map(_genus,roster & motif["known"]))
            item=dict(island_id=island_id,region=region,scope=scope,motif=name,
                      spatial_block=spatial_block,first_source=first_source,
                      n_axes=motif["n_axes"],
                      n_island_total=len(roster),
                      n_island_direct=island["n"],
                      n_island_positive=island["positive"],
                      island_frequency=island["freq"],
                      island_coupling=island["coupling"])
            for mode in MODES:
                s=src[mode]
                shared_annotated_genera=(
                    set(island_direct_genera) & set(source_cache[mode]["counts"])
                )
                num_island=sum(island_direct_genera[g] for g in shared_annotated_genera)
                num_source=sum(source_cache[mode]["counts"][g] for g in shared_annotated_genera)
                den_island=sum(scope_genus_full[scope][g] for g in shared_annotated_genera)
                den_source=sum(source_genus_full[mode][g] for g in shared_annotated_genera)
                doc_gap=(abs(num_island/den_island-num_source/den_source)
                         if den_island and den_source else np.nan)
                item.update({mode+"__n_source_total":len(sources[mode]),
                             mode+"__n_source_direct":s["n"],
                             mode+"__source_frequency":s["freq"],
                             mode+"__source_coupling":s["coupling"],
                             mode+"__n_shared_direct_genera":len(shared_annotated_genera),
                             mode+"__genus_conditioned_documentation_gap":doc_gap})
                if island["n"] and s["n"]:
                    item[mode+"__difference"]=island["freq"]-s["freq"]
                    item[mode+"__coupling_difference"]=(
                        island["coupling"]-s["coupling"])
                    for k,v in zip(PARTS,fourpart(roster,sources[mode],motif,
                                                   source_cache=source_cache[mode])):
                        item[mode+"__"+k]=v
                else:
                    item[mode+"__difference"]=np.nan
                    item[mode+"__coupling_difference"]=np.nan
                    for k in PARTS:item[mode+"__"+k]=np.nan
            rows.append(item)
    if len(rows)!=len(SCOPES)*len(MOTIFS):
        raise AssertionError("full island coverage contract broken")
    return rows


def aggregate_motifs(frame:pd.DataFrame,*,population:int=4505,
                     draws:int=999)->dict:
    need={"island_id","region","scope","motif","spatial_block","first_source",
          "n_island_total","n_island_direct","n_island_positive",
          "island_frequency","island_coupling"}
    fields=[]
    for mode in MODES:
        need|={mode+"__n_source_total",mode+"__n_source_direct",
               mode+"__difference",mode+"__coupling_difference",
               mode+"__source_frequency",
               mode+"__n_shared_direct_genera",
               mode+"__genus_conditioned_documentation_gap"}
        need|={mode+"__"+k for k in PARTS}
        fields.extend([mode+"__difference",mode+"__coupling_difference",
                       *(mode+"__"+k for k in PARTS)])
    if missing:=need-set(frame):
        raise ValueError("missing original GIFT motif records "+str(sorted(missing)))
    if (frame.island_id.nunique()!=population or
        len(frame)!=population*len(SCOPES)*len(MOTIFS) or
        frame.duplicated(["island_id","scope","motif"]).any()):
        raise ValueError("full original botanical island denominator lost")
    if not frame.region.isin(REGIMES).all() or not frame.scope.isin(SCOPES).all() or not frame.motif.isin(MOTIFS).all():
        raise ValueError("unknown GIFT real motif or status/region")
    if draws<199:raise ValueError("too few conditional independent block draws")
    for mode in MODES:
        matrix=frame[[mode+"__"+k for k in PARTS]]
        complete=matrix.notna().all(axis=1)
        if (matrix.notna().any(axis=1)!=complete).any():
            raise ValueError("incomplete signed botanical genus accounting")
        if (frame[mode+"__difference"].notna()!=complete).any():
            raise ValueError("source pair missingness mismatch")
        if complete.any() and not np.allclose(matrix.loc[complete].sum(axis=1),
                                              frame.loc[complete,mode+"__difference"],
                                              atol=1e-10,rtol=0):
            raise AssertionError("true multi-axis signed terms drift")
    result=[]
    for region in REGIMES:
        for scope in SCOPES:
            for motif in MOTIFS:
                p=frame.loc[frame.region.eq(region)&
                            frame.scope.eq(scope)&frame.motif.eq(motif)].copy()
                good=p.n_island_direct.ge(MIN_ISLAND)
                for mode in MODES:
                    good&=p[mode+"__n_source_direct"].ge(MIN_SOURCE)&p[mode+"__difference"].notna()
                q=p.loc[good].copy()
                # Compare trait recording WITHIN genera with actual independently
                # Direct-scored plants on both sides. Full-flora fractions
                # can both be near zero and falsely pass a 10pp gap check.
                balanced=pd.Series(True,index=q.index)
                for mode in MODES:
                    balanced&=(
                        q[mode+"__n_shared_direct_genera"].ge(2)
                        &q[mode+"__genus_conditioned_documentation_gap"].le(.10)
                    )
                panel=dict(region=region,scope=scope,motif=motif,
                           n_original_region_islands=len(p),
                           n_islands_two_source_supported=len(q),
                           n_documentation_balanced=int(balanced.sum()),
                           n_islands_with_at_least_two_shared_annotated_genera_both_modes=int((
                               q.geo5__n_shared_direct_genera.ge(2)&
                               q.climate10__n_shared_direct_genera.ge(2)).sum()),
                           n_positive_island_trait_occurrences=int(q.n_island_positive.sum()),
                           median_known_target_species=float(q.n_island_direct.median()) if len(q) else None)
                for sample,sub in (("all",q),("doc_gap_le10pp",q.loc[balanced])):
                    for grouping in ("spatial_block","first_source"):
                        groups=sub.dropna(subset=[grouping]).groupby(grouping,sort=True)
                        out=dict(n_islands=len(sub),n_groups=groups.ngroups)
                        if len(sub)<MIN_REGION_ISLANDS or groups.ngroups<MIN_BLOCKS:
                            out["status"]="NOT_IDENTIFIABLE_DIRECT_SOURCE_OR_BLOCK_SUPPORT"
                        else:
                            g=groups[fields].mean().to_numpy(float)
                            ave=g.mean(axis=0)
                            seed=int.from_bytes(sha256(
                                f"{region}|{scope}|{motif}|{sample}|{grouping}".encode()
                            ).digest()[:8],"big")
                            rng=np.random.default_rng(seed)
                            bs=g[rng.integers(len(g),size=(draws,len(g)))].mean(1)
                            for j,k in enumerate(fields):
                                out[k+"__equal_group_mean"]=float(ave[j])
                                out[k+"__bootstrap95"]=[
                                    float(z) for z in np.quantile(bs[:,j],[.025,.975])]
                            out["status"]="EXPLORATORY_REAL_SAME_PLANT_MOTIF_SOURCE_COMPOSITION"
                        panel[sample+"__"+grouping]=out
                result.append(panel)
    return dict(status="ORIGINAL_4505_REAL_THREE_FOUR_AXIS_SYNDROME_COMPATIBLE_MOTIFS",
                n_original_GIFT_botanical_islands=population,
                n_original_island_scope_motif_rows=len(frame),
                n_motifs=len(MOTIFS),
                n_reporting_panels=len(result),
                panels=result,
                limits=dict(
                    exact_accepted_species_all_axes_direct=True,
                    no_imputation_missing_flower_traits=True,
                    four_signed_genus_terms_not_causal=True,
                    colour_ontology_red_pink_not_separate_red=True,
                    yellow_orange_not_separate_orange=True,
                    no_specific_bird_butterfly_hawkmoth_bee_pollinator_labels=True,
                    visitor_effectiveness_not_measured=True,
                    candidate_source_not_true_seed_arrivals=True,
                    no_global_complete_trait_requirement=True,
                    genus_conditioned_direct_documentation_balance_both_sources=True,
                    requires_two_shared_annotated_genera_for_strict_documentation=True,
                    original_4505_islands_preserved=True))
