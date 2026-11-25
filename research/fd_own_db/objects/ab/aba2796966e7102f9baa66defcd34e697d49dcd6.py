"""Evidence-qualified positive multi-guild co-documentation on true floral traits.

For DoPI v2, an accepted plant may have BOTH guilds positively documented
inside the same original publication + published coordinate/grid + year/month.
That is observed co-documentation of positive visits, NOT an observation that
another guild is absent, not simultaneous visitation of the same flower
individual, not efficient pollination and not island-local evolution.

Compare floral motif frequencies among plants with both guilds documented
against the positively documented union at that exact source location/period.
A separate exact expected positive shared-plant morphology proportion
conditions on observed botanical GENUS counts, motif frequencies and
positive co-documentation counts inside every published site.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256

import numpy as np

from .dopi_gift_guild_motif_bridge import PAIRS, botanical_genus
from .dopi_site_period_guild import RESOLUTIONS, EVIDENCE

MIN_STUDIES=10
MIN_STUDIES_MOTIF_UNION=10
MIN_PLANTS_POSITIVE_UNION=3
MIN_STUDIES_GENUS_INFORMATIVE=10


def _report_bootstrap(study_values:dict[str,list[dict]],
                      *,kind:str,seed:str,draws:int) -> dict:
    n=len(study_values)
    npositive=sum(any(x["motif_in_union"] for x in group)
                  for group in study_values.values())
    nco=sum(any(x["motif_in_both"] for x in group)
            for group in study_values.values())
    nsites=sum(len(group) for group in study_values.values())
    out=dict(n_source_publications=n,
             n_source_publications_with_observed_motif_in_positive_union=npositive,
             n_source_publications_with_observed_motif_in_both_guilds=nco,
             n_source_site_periods_with_positive_both_guild_plant=nsites)
    if n<MIN_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_SOURCE_PUBLICATIONS"
        return out
    if npositive<MIN_STUDIES_MOTIF_UNION:
        out["status"]="NOT_IDENTIFIABLE_MOTIF_POSITIVE_DOCUMENTATION"
        return out
    values=np.asarray([
        np.mean([q["delta"] for q in study_values[paper]])
        for paper in sorted(study_values)
    ],float)
    if not np.isfinite(values).all():
        raise ValueError("nonfinite positive co-visit or botanical genus expected frequency")
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    bs=values[rng.integers(len(values),size=(draws,len(values)))].mean(axis=1)
    out.update(status=kind, study_equal_difference=float(values.mean()),
               study_cluster_bootstrap95=[float(z) for z in np.quantile(bs,[.025,.975])])
    return out


def summarize_positive_coaccess(index:dict,*,motifs:dict,
                                draws:int=999)->dict:
    """Original source-positive co-visitation in same published site and period.

    Each site is a source publication×as-reported geographic key×year
    or a source publication×geographic key×year×month. The
    public output only contains independent publication counts.
    """
    if draws<199:
        raise ValueError("insufficient original DoPI study-cluster bootstrap")
    if not index.get("true_guild_absence_not_inferred"):
        raise ValueError("original positive guild index required, no pseudoabsence")
    units=index.get("private_units")
    if not units:
        raise ValueError("no original source-positive DoPI site observations")
    if not motifs:
        raise ValueError("source GIFT floral morphology Direct states missing")
    # Index source site/year/month availability, never make fake nonvisitation
    # species labels when one of the guilds has no source record at that site.
    both_keys={}
    for r in RESOLUTIONS:
        for e in EVIDENCE:
            per_guild=defaultdict(set)
            for res,tier,key,guild in units:
                if res==r and tier==e:
                    per_guild[guild].add(key)
            for ga,gb in PAIRS:
                both_keys[(r,e,ga,gb)]=sorted(
                    per_guild[ga]&per_guild[gb])
    reports=[]
    for resolution in RESOLUTIONS:
      for tier in EVIDENCE:
       for ga,gb in PAIRS:
        keys=both_keys[(resolution,tier,ga,gb)]
        for name,motif in motifs.items():
            known=set(motif["known"])
            positive=set(motif["joint"])
            raw_by_pub=defaultdict(list)
            within_genus_by_pub=defaultdict(list)
            n_site_both=0
            n_site_with_actual_same_plant=0
            n_site_with_motif_co_documented=0
            n_site_with_informative_genus=0
            n_site_genus_groups_variable=0
            n_distinct_co_visited=set()
            for key in keys:
                a=units[(resolution,tier,key,ga)]&known
                b=units[(resolution,tier,key,gb)]&known
                union=a|b
                shared=a&b
                if len(union)<MIN_PLANTS_POSITIVE_UNION or not shared:
                    continue
                n_site_both+=1
                n_site_with_actual_same_plant+=1
                n_distinct_co_visited.update(shared)
                paper=key[0]
                k_union=len(union&positive)
                k_shared=len(shared&positive)
                if k_shared:n_site_with_motif_co_documented+=1
                # Baseline: the exact same site's union of any positively
                # documented guild plant names. A missing second guild record
                # is "not co-documented", NOT absence of its visitation.
                raw_by_pub[paper].append(dict(
                    delta=k_shared/len(shared)-k_union/len(union),
                    motif_in_union=bool(k_union),
                    motif_in_both=bool(k_shared),
                ))
                by_genus=defaultdict(set)
                for plant in union:by_genus[botanical_genus(plant)].add(plant)
                obs_minus_expected=0.
                n_co_var=0
                var_genera=0
                k_motif_union=0
                k_motif_both=0
                for g,pool in by_genus.items():
                    n=len(pool)
                    k=len(pool&positive)
                    c=len(pool&shared)
                    # Morphology must actually vary within genus AND
                    # documentation must vary (co vs one guild only).
                    if n<=1 or not (0<k<n) or not (0<c<n):
                        continue
                    o=len(pool&shared&positive)
                    obs_minus_expected+=o-c*k/n
                    n_co_var+=c
                    var_genera+=1
                    k_motif_union+=k
                    k_motif_both+=o
                if var_genera:
                    n_site_with_informative_genus+=1
                    n_site_genus_groups_variable+=var_genera
                    within_genus_by_pub[paper].append(dict(
                        delta=obs_minus_expected/n_co_var,
                        motif_in_union=bool(k_motif_union),
                        motif_in_both=bool(k_motif_both),
                    ))
            cell=dict(
                published_location_scope=resolution,
                visitor_evidence=tier,
                guild_a=ga,guild_b=gb,
                floral_motif=name,
                n_original_direct_flower_axes=motif.get("axis_count",
                    motif.get("n_axes",len(motif.get("trait_axes",())))),
                n_source_site_periods_with_both_guilds_positive=len(keys),
                n_source_site_periods_with_at_least_one_same_plant_positive_both=n_site_both,
                n_source_site_periods_with_actual_motif_positive_both=n_site_with_motif_co_documented,
                n_distinct_source_positive_co_documented_plant_names=len(n_distinct_co_visited),
                n_site_periods_with_informative_within_genus_morphology=n_site_with_informative_genus,
                n_informative_source_site_genus_groups=n_site_genus_groups_variable)
            cell["co_documented_vs_positive_union"]=_report_bootstrap(
                raw_by_pub,kind="DESCRIPTIVE_REAL_GUILD_CO_DOCUMENTATION",
                seed=f"{resolution}|{tier}|{ga}|{gb}|{name}|observed",
                draws=draws)
            cell["genus_fixed_co_documentation_excess"]=_report_bootstrap(
                within_genus_by_pub,
                kind="DESCRIPTIVE_GENUS_CONDITIONED_CO_DOCUMENTATION_EXCESS",
                seed=f"{resolution}|{tier}|{ga}|{gb}|{name}|genus_fixed",
                draws=draws)
            reports.append(cell)
    return dict(
        status="DOPI_V2_GIFT_POSITIVE_MULTI_GUILD_COMMON_SITE_MORPHOLOGY",
        original_GIFT_botanical_flora_population=4505,
        n_original_DoPI_records=index["source_scan_counts"].get(
            "n_original_source_rows",0),
        n_original_positive_site_guild_units=index[
            "n_original_site_guild_positive_units"],
        n_original_source_publications_with_site_and_time=index[
            "n_original_source_publication_keys"],
        n_original_motifs=len(motifs),
        n_guild_pairs=len(PAIRS),
        n_panels=len(reports),
        panels=reports,
        limitations=dict(
            same_accepted_plant_positively_documented_by_both_real_guilds=True,
            guild_not_reported_on_same_plant_is_not_true_absence=True,
            exact_published_location_and_year_month_not_same_flower_session=True,
            source_positive_co_documentation_not_pollination_effectiveness=True,
            source_sample_observation_intensity_and_plant_abundance_uncontrolled=True,
            genotype_and_island_local_morphology_not_measured=True,
            genus_fixed_null_conditions_on_actual_positive_site_composition=True,
            all_4505_GIFT_islands_remain_botanical_parent_population=True,
            European_DoPI_not_actual_tropical_island_guild=True,
            no_private_original_site_study_or_taxon_identifier_serialized=True))
