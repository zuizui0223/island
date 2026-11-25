"""Genus-EQUAL, study-matched external visitor/floral morphology diagnostic.

The earlier DoPI 'shared genus' analysis dropped absent genera but pooled all
plants, so *different weights within the shared genera remained*. This module
instead gives EVERY shared genus within a source publication equal weight,
compares the exact same publication cohort to its original pooled/full
contrast, and runs a separate within-genus-variable opportunity gate.

These are source-positive flower visits, NOT true absence/preference,
independent sites, effective pollination, island-local observed visitors
or causal pollen-mediated floral evolution.
"""
from __future__ import annotations
from collections import defaultdict
from hashlib import sha256
import numpy as np

from .dopi_gift_guild_motif_bridge import (
    PAIRS, SCOPES, MIN_DIRECT_PLANTS_EACH_GUILD, MIN_STUDIES,
    botanical_genus,
)

MIN_MOTIF_POSITIVE_STUDIES=10
MIN_SHARED_GENERA_PER_STUDY=2
MIN_INFORMATIVE_GENERA_PER_STUDY=1


def _summarize(studies:list[dict], *, label:str,seed:str,draws:int) -> dict:
    """Study-cluster confidence intervals, not p-values."""
    columns=("raw_study_diff", "shared_genus_pooled_diff",
             "genus_equal_diff", "raw_minus_genus_equal",
             "shared_pooled_minus_genus_equal")
    npos=sum(item["any_positive"] for item in studies)
    out=dict(
        n_studies=len(studies),n_studies_with_motif_positive=npos,
        n_genus_study_positive_units=sum(q["n_common_genera"] for q in studies),
        n_studies_with_both_guilds_same_accepted_plant=sum(
            q["has_same_species"] for q in studies),
        n_studies_with_variable_genus=sum(
            q["n_informative_genera"]>0 for q in studies),
    )
    if len(studies)<MIN_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_STUDIES"
        return out
    if npos<MIN_MOTIF_POSITIVE_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_POSITIVE_MOTIF_STUDY_SUPPORT"
        return out
    values=np.asarray([[q[col] for col in columns] for q in studies],float)
    if not np.isfinite(values).all():
        raise ValueError("invalid matched study+genus flower proportion")
    rng=np.random.default_rng(int.from_bytes(
        sha256(seed.encode()).digest()[:8],"big"))
    samples=values[rng.integers(len(values),size=(draws,len(values)))].mean(axis=1)
    for j,k in enumerate(columns):
        out[k]=float(values[:,j].mean())
        out[k+"_ci95"]=[float(z) for z in np.quantile(samples[:,j],[.025,.975])]
    out["status"]="DESCRIPTIVE_POSITIVE_VISITOR_STUDY_EQUAL_GENUS_FIXED"
    return out


def genus_equal_dopi_audit(index:dict,*,motifs:dict,draws:int=999)->dict:
    """Original DoPI positive source papers; 4 evidence scopes x 5 pairs x 27 motifs.

    No names or bibliography IDs are returned. Each plant can occur in
    multiple guild-positive lists; this is NOT a pairwise non-visit test.
    """
    if draws<199:
        raise ValueError("insufficient whole-publication bootstrap")
    if not motifs or not index.get("no_true_absence_inferred"):
        raise ValueError("positive source-only guild index missing")
    units=index.get("indexed_units")
    if not units:
        raise ValueError("no source-linked GIFT botanical visitor-guild units")
    studies=sorted({k[0] for k in units})
    panels=[]
    for scope in SCOPES:
        for guild_a,guild_b in PAIRS:
            pool=[s for s in studies if
                  (s,guild_a,scope) in units and (s,guild_b,scope) in units]
            for name,motif in motifs.items():
                known=set(motif["known"])
                positive=set(motif["joint"])
                matched=[]
                variable=[]
                for study in pool:
                    a=units[(study,guild_a,scope)]&known
                    b=units[(study,guild_b,scope)]&known
                    if min(len(a),len(b))<MIN_DIRECT_PLANTS_EACH_GUILD:
                        continue
                    gen_a=defaultdict(set);gen_b=defaultdict(set)
                    for sp in a:gen_a[botanical_genus(sp)].add(sp)
                    for sp in b:gen_b[botanical_genus(sp)].add(sp)
                    common=set(gen_a)&set(gen_b)
                    if len(common)<MIN_SHARED_GENERA_PER_STUDY:
                        continue
                    ag=set().union(*(gen_a[g] for g in common))
                    bg=set().union(*(gen_b[g] for g in common))
                    raw_diff=len(a&positive)/len(a)-len(b&positive)/len(b)
                    pooled_diff=len(ag&positive)/len(ag)-len(bg&positive)/len(bg)
                    diffs=[]
                    informative=[]
                    positive_any=False
                    for g in common:
                        aa=gen_a[g];bb=gen_b[g]
                        local_diff=len(aa&positive)/len(aa)-len(bb&positive)/len(bb)
                        diffs.append(local_diff)
                        if (aa|bb)&positive:
                            positive_any=True
                        # Actual genus provides both positive and negative
                        # flower states, with at least two distinct taxa.
                        if 0<len((aa|bb)&positive)<len(aa|bb):
                            informative.append((g,local_diff))
                    equal=float(np.mean(diffs))
                    base=dict(
                        raw_study_diff=float(raw_diff),
                        shared_genus_pooled_diff=float(pooled_diff),
                        genus_equal_diff=equal,
                        raw_minus_genus_equal=float(raw_diff-equal),
                        shared_pooled_minus_genus_equal=float(pooled_diff-equal),
                        n_common_genera=len(common),
                        n_informative_genera=len(informative),
                        any_positive=bool(positive_any),
                        has_same_species=bool(a&b),
                    )
                    matched.append(base)
                    if len(informative)>=MIN_INFORMATIVE_GENERA_PER_STUDY:
                        only_equal=float(np.mean([v for _,v in informative]))
                        # Compare on precisely the same informative genera
                        # (otherwise taxonomy/phenotype sample masks move).
                        va=set().union(*(gen_a[g] for g,_ in informative))
                        vb=set().union(*(gen_b[g] for g,_ in informative))
                        pa=len(va&positive)/len(va)-len(vb&positive)/len(vb)
                        # Shared-genus POOLED reference on same genus subset;
                        # the raw reference equals it in this cohort.
                        variable.append(dict(
                            raw_study_diff=float(pa),
                            shared_genus_pooled_diff=float(pa),
                            genus_equal_diff=only_equal,
                            raw_minus_genus_equal=float(pa-only_equal),
                            shared_pooled_minus_genus_equal=float(pa-only_equal),
                            n_common_genera=len(informative),
                            n_informative_genera=len(informative),
                            any_positive=bool((va|vb)&positive),
                            has_same_species=bool(va&vb),
                        ))
                idmeta=dict(
                    guild_a=guild_a,guild_b=guild_b,
                    floral_motif=name,
                    evidence_scope=scope,
                    n_recorded_flower_axes=motif.get("axis_count",
                        motif.get("n_axes",len(motif.get("trait_axes",())))),
                    n_source_publications_both_positive_guilds=len(pool),
                    n_eligible_publications_with_two_shared_genera=len(matched),
                    n_publications_with_variable_motif_shared_genus=len(variable),
                )
                idmeta["genus_equal_all"]=_summarize(
                    matched,label="all",seed=f"{scope}|{guild_a}|{guild_b}|{name}|all",draws=draws)
                idmeta["genus_equal_informative"]=_summarize(
                    variable,label="informative",
                    seed=f"{scope}|{guild_a}|{guild_b}|{name}|variable",draws=draws)
                panels.append(idmeta)
    return dict(
        status="DOPI_2026_GIFT_TRUE_FLORAL_WITHIN_STUDY_GENUS_EQUAL_AUDIT",
        n_original_GIFT_botanical_islands=4505,
        n_original_external_publication_keys=index["n_distinct_publication_keys"],
        n_original_DoPI_source_rows=index["original_source_row_counts"].get("source_rows",0),
        n_original_GIFT_island_accepted_plants_matching_DoPI=index[
            "n_original_island_accepted_plants_matching_DoPI"],
        n_original_guild_pairs=len(PAIRS),
        n_original_flower_motifs=len(motifs),
        n_evidence_tiers=len(SCOPES),
        n_panels=len(panels),
        panels=panels,
        interpretation=dict(
            European_positive_visitor_source_not_island_local=True,
            does_not_infer_undocumented_guild_absence=True,
            compares_same_study_and_equal_weights_each_shared_genus=True,
            genus_structure_control_not_species_ancestry_or_phylogeny_test=True,
            genotype_or_local_phenotype_not_observed=True,
            genus_variation_gate_descriptive_not_pollinator_fitness=True,
            no_true_effective_guild_shift_inference=True,
            source_study_bootstrap_not_independent_species_events=True,
            no_original_accepted_plant_or_study_identifiers_exported=True,
        ),
    )
