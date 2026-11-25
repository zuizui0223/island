"""Flower motifs vs DoPI observed positive visitor-guild PORTFOLIO RICHNESS.

Controls for source observation effort by EXACT fixed source study + published
site + year + month, each plant's DISTINCT observed visitor taxon count and
the site's overall positive visitor guild counts. The finite-sampling
reference is random reassignment of recorded visitor GUILD labels among
already observed positive plant-visitors, NOT a source-justified statement
that unreported guilds were absent. This is an observation-conditioned null,
not a physical floral access or historical pollinator-shift null.

Source GIFT original botanical 4505 islands remain the parent population;
their plant species, site, publication and observer identities are not
serialized outside the protected original island workflow.
"""
from __future__ import annotations

from collections import Counter,defaultdict
from hashlib import sha256
import math
import re
import numpy as np

from .dopi_site_period_guild import _site, _time
from .dopi_gift_guild_motif_bridge import (
    _s,_study,GUILDS,botanical_genus,
)

MIN_VISITOR_TAXA_PER_PLANT=2
MIN_POSITIVE_AND_NEGATIVE_PER_SITE=2
MIN_STUDIES=10
MIN_INFORMATIVE_MULTIGUILD_STUDIES=10
_BINOMIAL=re.compile(r"^([A-Z][A-Za-z-]+)\s+([a-z][a-z-]+)(?=\s|[.,;()]|$)")
_UNRESOLVED_EPITHETS=frozenset({"sp","spp","cf","aff","indet","unknown","nr"})


def accepted_visitor_binomial(raw:str)->str|None:
    """Keep only real binomial visitor IDs, never inflate 'Bombus sp' breadth."""
    q=str(raw or "").strip().replace("_"," ")
    m=_BINOMIAL.match(q)
    if not m or m.group(2) in _UNRESOLVED_EPITHETS:
        return None
    return (m.group(1)+" "+m.group(2)).casefold()


def expected_positive_guild_richness(n:int,counts:dict[str,int])->float:
    """Exact expectation under fixed visitor-guild marginals (hypergeometric).

    The sample space is the study-site/year/month's actually documented
    GIFT-listed plant×visitor-species positive edges, not all insect taxa.
    """
    if not isinstance(n,int) or not 0<=n<=sum(counts.values()):
        raise ValueError("invalid directly documented plant visitor opportunity")
    if any(not isinstance(v,int) or v<0 for v in counts.values()):
        raise ValueError("invalid original positive guild source counts")
    N=sum(counts.values())
    if not N or not n:return 0.
    result=0.
    for K in counts.values():
        if K<=0:continue
        if n>N-K:
            p_zero=0.
        else:
            # Stable product of ratios instead of gigantic binomial numbers.
            p_zero=1.
            for i in range(n):
                p_zero*=(N-K-i)/(N-i)
        result+=1.-p_zero
    return float(result)


def positive_visitor_effort_index(rows,*,fields:set[str],
                                  GIFT_accepted_names:set[str],
                                  motifs:dict, visitor_classifier)->dict:
    """Source-private plant×original visitor species, site/month positive edges."""
    required={"plant_species","pollinator_species","pollinator_family",
        "pollinator_order","doi","publication","pollination",
        "country","grid_letter","grid_code","latitude","longitude","year","month"}
    if missing:=required-fields:
        raise ValueError(f"missing original DoPI study/site/time/visitor fields {sorted(missing)}")
    if not motifs or not GIFT_accepted_names:
        raise ValueError("original GIFT accepted island botanical flora missing")
    known=set().union(*(set(z["known"]) for z in motifs.values()))
    edges=defaultdict(set)
    counters=Counter()
    for row in rows:
        counters["original_DoPI_v2_source_rows"]+=1
        plant=_s(row.get("plant_species"))
        if plant not in GIFT_accepted_names or plant not in known:
            continue
        counters["original_GIFT_direct_flower_trait_accepted_source_rows"]+=1
        study=_study(row)
        site=_site(row)
        when=_time(row)
        if not study or site is None or when is None or when[1] is None:
            counters["original_study_site_year_month_unsupported"]+=1
            continue
        grade=_s(row.get("pollination"))
        if grade not in ("1","2","3"):
            continue
        animal=_s(row.get("pollinator_species"))
        accepted=accepted_visitor_binomial(animal)
        if accepted is None:
            counters["original_visitor_taxon_not_resolved_to_species"]+=1
            continue
        guild=visitor_classifier(row.get("pollinator_order"),
                                 row.get("pollinator_family"),animal)
        if guild not in GUILDS:
            continue
        study_key=(study,site,when[0],when[1])
        edges[study_key].add((plant,accepted,guild))
        counters["original_evidence_qualified_source_rows"]+=1
    return dict(private_positive_edges=edges,
                n_original_source_studies_qualified=len({k[0] for k in edges}),
                n_original_site_month_units=len(edges),
                source_row_counts=dict(counters),
                documented_guild_absence_not_inferred=True)


def _bootstrap(studies:dict[str,list[tuple[float,float,float]]],
               *,motif_studies:set[str],seed:str,draws:int)->dict:
    n=len(studies)
    out=dict(n_original_source_publications=n,
             n_studies_with_positive_true_multiguild_documentation=len(motif_studies),
             n_original_site_month_comparisons=sum(map(len,studies.values())))
    if n<MIN_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_INDEPENDENT_PUBLICATIONS"
        return out
    if len(motif_studies)<MIN_INFORMATIVE_MULTIGUILD_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_REAL_MULTIGUILD_DOCUMENTATION"
        return out
    arr=np.asarray([np.mean(studies[k],axis=0) for k in sorted(studies)],float)
    if arr.shape!=(n,3) or not np.isfinite(arr).all():
        raise ValueError("nonfinite source visitor richness expectation")
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    bs=arr[rng.integers(n,size=(draws,n))].mean(axis=1)
    for j,label in enumerate(("observed_positive_guild_richness_difference",
                              "expected_difference_fixed_visitor_opportunity",
                              "excess_over_fixed_source_effort_null")):
        out[label]=float(arr[:,j].mean())
        out[label+"_ci95"]=[float(z) for z in np.quantile(bs[:,j],[.025,.975])]
    out["status"]="DESCRIPTIVE_EFFORT_CONDITIONED_POSITIVE_VISITOR_GUILD_PORTFOLIO"
    return out


def guild_richness_motif_audit(index:dict,*,motifs:dict,draws:int=999)->dict:
    """Directly positive source visitor-guild richness vs flower states.

    Cluster at publication level and report separate genus fixed sensitivity.
    Lacking morphological states/visitor taxa/papers is non-identifiable,
    NEVER a zero ecological guild repertoire.
    """
    if draws<199:raise ValueError("insufficient source study bootstrap draws")
    if not index.get("documented_guild_absence_not_inferred"):
        raise ValueError("source positive visitor only input required")
    private=index.get("private_positive_edges")
    if not private or not motifs:
        raise ValueError("no original study-site positive visitor evidence")
    plant_profiles={}
    per_site={}
    for key,edges in private.items():
        per_plant=defaultdict(set)
        source_guild_counts=Counter()
        for plant,visitor,guild in edges:
            per_plant[plant].add((visitor,guild))
            source_guild_counts[guild]+=1
        # Only species with >=2 recorded *different visitor species* can
        # demonstrate a visitor-guild breadth beyond one.
        rows={}
        for plant,visitors in per_plant.items():
            # Same visitor species should not be typed to different guilds.
            if len({v for v,g in visitors})!=len(visitors):
                raise ValueError("source original visitor species has conflicting guild taxonomy")
            n=len(visitors)
            if n<MIN_VISITOR_TAXA_PER_PLANT:
                continue
            actual=len({g for v,g in visitors})
            expected=expected_positive_guild_richness(n,dict(source_guild_counts))
            rows[plant]=(actual,expected,actual-expected)
        if len(rows)>=MIN_POSITIVE_AND_NEGATIVE_PER_SITE*2:
            per_site[key]=rows
    panels=[]
    for motif_name,motif in motifs.items():
        known=set(motif["known"])
        positive=set(motif["joint"])
        all_studies=defaultdict(list)
        genus_studies=defaultdict(list)
        yes_multi=set()
        genus_multi=set()
        n_sites=0;n_genus_sites=0
        n_trait_scored_plant_site_units=0
        n_multi_guild_positive_plant_site_units=0
        n_variation_genus_units=0
        seen_species=set()
        for key,site_data in per_site.items():
            real={p:z for p,z in site_data.items() if p in known}
            positives=[p for p in real if p in positive]
            negatives=[p for p in real if p not in positive]
            if min(len(positives),len(negatives))<MIN_POSITIVE_AND_NEGATIVE_PER_SITE:
                continue
            n_sites+=1
            n_trait_scored_plant_site_units+=len(real)
            seen_species.update(real)
            n_multi_guild_positive_plant_site_units+=sum(real[p][0]>=2 for p in real)
            paper=key[0]
            y=np.mean([real[p] for p in positives],axis=0)
            n=np.mean([real[p] for p in negatives],axis=0)
            all_studies[paper].append(tuple((y-n).tolist()))
            if any(real[p][0]>=2 for p in real):
                yes_multi.add(paper)
            genera=defaultdict(list)
            for plant in real:genera[botanical_genus(plant)].append(plant)
            diffs=[]
            informative_multi=False
            for g,plants in genera.items():
                ap=[p for p in plants if p in positive]
                bp=[p for p in plants if p not in positive]
                if not ap or not bp:
                    continue
                n_variation_genus_units+=1
                v1=np.mean([real[p] for p in ap],axis=0)
                v0=np.mean([real[p] for p in bp],axis=0)
                diffs.append(v1-v0)
                if any(real[p][0]>=2 for p in plants):
                    informative_multi=True
            if diffs:
                n_genus_sites+=1
                genus_studies[paper].append(tuple(
                    np.mean(diffs,axis=0).tolist()))
                if informative_multi:
                    genus_multi.add(paper)
        meta=dict(
            floral_motif=motif_name,
            n_recorded_direct_flower_axes=motif.get("axis_count",
                motif.get("n_axes",len(motif.get("trait_axes",())))),
            n_original_site_months_with_both_motif_states_and_two_visitor_taxa=n_sites,
            n_original_site_months_with_same_genus_both_flower_states=n_genus_sites,
            n_original_variable_flower_genus_site_months=n_variation_genus_units,
            n_original_trait_scored_plant_site_month_units=n_trait_scored_plant_site_units,
            n_real_positive_multiguild_plant_site_month_units=n_multi_guild_positive_plant_site_units,
            n_distinct_trait_scored_species_in_any_eligible_site=len(seen_species))
        meta["study_site"]=_bootstrap(all_studies,
            motif_studies=yes_multi,seed=f"{motif_name}|source",
            draws=draws)
        meta["same_genus_study_site"]=_bootstrap(genus_studies,
            motif_studies=genus_multi,seed=f"{motif_name}|genus",
            draws=draws)
        panels.append(meta)
    return dict(
        status="DOPI_GIFT_POSITIVE_GUILD_RICHNESS_FIXED_EFFORT_SITE_NULL",
        original_GIFT_botanical_island_population=4505,
        n_original_DOPI_source_rows=index["source_row_counts"].get(
            "original_DoPI_v2_source_rows",0),
        n_original_source_studies_qualified=index[
            "n_original_source_studies_qualified"],
        n_original_published_site_month_positive_units=index[
            "n_original_site_month_units"],
        n_scored_source_sites_with_two_or_more_visitor_taxa=int(len(per_site)),
        source_row_counts=index["source_row_counts"],
        n_direct_flower_profiles=len(motifs),
        n_reporting_panels=len(panels),
        panels=panels,
        limitations=dict(
            fixed_plant_positive_visitor_species_count_per_site=True,
            fixed_site_positive_guild_margins=True,
            only_documented_visitor_species_positive_edges=True,
            guild_nondetection_not_biological_absence=True,
            expected_guild_breadth_not_true_ecological_absence=True,
            same_study_published_site_year_month_required=True,
            original_source_plant_abundance_and_survey_effort_unknown=True,
            source_event_counts_not_equal_visitation_intensity=True,
            no_original_site_study_or_accepted_plant_identifiers_exported=True,
            original_GIFT_4505_not_source_site_visitations=True,
            no_island_local_adaptive_syndrome_inference=True))
