"""Global 8k-island direct-trait functional diversity WITHOUT complete-case gate.

The unit of the main analysis is one island × ONE observed trait. The
three independent flower/other reproductive Direct traits are NOT joined
before computing their own diversity. For two-trait coupling ONLY the
specified pair must be directly observed on the same accepted species.
Full colour+form+symmetry completion is reserved as a historical sensitivity.

Source of original status: protected GIFT island flora with independent
native/introduced/unresolved and endemic labels. Source plant trait states
are species-wide high/medium Direct evidence, NOT verified island-specific
petal colour or reproductive fitness.

No imputation, no invented zeros, no inference of flower pollinator success
or geologically induced adaptive evolution. Source/geography/spatial IDs
remain private. The returned report contains only aggregate counts and
regional 10-degree-block / GIFT-nearest-source cluster-bootstrap statistics.
"""
from __future__ import annotations

from collections import Counter,defaultdict
from hashlib import sha256

import numpy as np
import pandas as pd

from .geological_endemic_trait_opportunity import TRAITS
from .regional_native_provenance_volume import cluster_equal_bootstrap


PAIR_TRAITS=(("flower_primary_color","floral_form"),
             ("flower_primary_color","floral_symmetry"),
             ("floral_form","floral_symmetry"))
SCOPE=("all_recorded","original_native")
N_RAO=5
OBSERVATION_GATE=(2,5,10,30)
GEOGRAPHIC_ORIGINS=("oceanic","continental")


def _seed(*x:object)->int:
    return int.from_bytes(sha256("|".join(map(str,x)).encode()).digest()[:8],"big")


def strict_single_trait_original(ledger:pd.DataFrame)->dict[str,dict[str,str]]:
    expected={"accepted_species","trait_name","normalized_value",
              "evidence_scope","evidence_quality"}
    if not expected.issubset(ledger):
        raise ValueError("Source original Direct botanical evidence columns incomplete")
    q=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()
    ].copy()
    q["normalized_value"]=q.normalized_value.astype(str).str.strip()
    q=q.loc[q.normalized_value.ne("")&
            ~q.normalized_value.str.contains("|",regex=False)&
            ~q.normalized_value.str.lower().isin(("none","nan","null","unknown","unresolved"))]
    if q.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Original Direct source has multiple conflicting species × trait rows")
    if q.accepted_species.isna().any() or q.trait_name.isna().any():
        raise ValueError("Original accepted plant/trait ID missing")
    return {trait:q.loc[q.trait_name.eq(trait)].set_index(
        "accepted_species").normalized_value.to_dict()
        for trait in TRAITS}


def _same_distinct_probability(categories:list[str])->float:
    n=len(categories)
    if n<2:
        raise ValueError("Distinct plant states need 2+ source species")
    counts=Counter(categories)
    return float(sum(v*(v-1) for v in counts.values())/(n*(n-1)))


def single_trait_fd(states:dict[str,str],eligible_names:set[str])->dict:
    """One trait per species, original distinct accepted plant names.

    The unbiased pair mismatch of the entire observed island flora is
    multiplied by 4/5 to equal the EXPECTED five-species finite-sample
    categorical Rao Q under sampling WITHOUT replacement.
    n>=5 is a trait-SPECIFIC reporting threshold, not a complete-case gate.
    """
    names=sorted(set(states)&set(eligible_names))
    n=len(names)
    raw=[states[x] for x in names]
    genus={x:x.split()[0] for x in names if len(x.split())>=2}
    by=defaultdict(Counter)
    for x in names:
        if x in genus:
            by[genus[x]][states[x]]+=1
    out=dict(n_direct_species=n,n_distinct_recorded_genera=len(by),
             n_distinct_original_trait_categories=len(set(raw)))
    if n>=N_RAO:
        out["finite_5_species_rao"]=float((N_RAO-1)/N_RAO*
            (1-_same_distinct_probability(raw)))
    else:
        out["finite_5_species_rao"]=float("nan")
    # Compare one independently observed source species per genus
    # integrated analytically over every possible choice. This does not
    # invent a genus/category, nor control botanical lineage identity.
    if len(by)>=N_RAO:
        diversity_g=len(by)
        weights={g:{s:n/n_all for s,n in counter.items()
                  for n_all in [sum(counter.values())]}
                 for g,counter in by.items()}
        mass=Counter()
        mass_sq=Counter()
        for states_g in weights.values():
            for category,p in states_g.items():
                mass[category]+=p
                mass_sq[category]+=p*p
        p_same_between_genera=sum(
            mass[c]*mass[c]-mass_sq[c] for c in mass
        )/(diversity_g*(diversity_g-1))
        out["finite_5_genus_equal_rao"]=float(
            (N_RAO-1)/N_RAO*(1-p_same_between_genera))
    else:
        out["finite_5_genus_equal_rao"]=float("nan")
    return out


def pair_trait_fd(left:dict[str,str],right:dict[str,str],
                  eligible_names:set[str])->dict:
    """Pairwise colour/form/symmetry coupling; only TWO traits required."""
    names=sorted(set(left)&set(right)&set(eligible_names))
    n=len(names)
    out=dict(n_two_trait_original_species=n)
    if n<N_RAO:
        return dict(**out,finite_5_joint_rao=float("nan"),
                    finite_5_pairing_departure=float("nan"))
    a=[left[x] for x in names]
    b=[right[x] for x in names]
    d=list(zip(a,b))
    finite=(N_RAO-1)/N_RAO
    # Permuting all colour states across original real species cannot
    # change either marginal: expected independently shuffled pair SAME
    # probability = P(same left)*P(same right).
    pa=_same_distinct_probability(a)
    pb=_same_distinct_probability(b)
    pab=_same_distinct_probability(d)
    return dict(**out,
                finite_5_joint_rao=float(finite*(1-pab)),
                finite_5_pairing_departure=float(finite*(pa*pb-pab)))


def original_global_partial_traits(
    ledger:pd.DataFrame,flora:pd.DataFrame,geology:pd.DataFrame,
    assignments:pd.DataFrame,mainland:pd.DataFrame,*,
    seed:int=20261114,bootstrap_draws:int=999,
    return_private_island_rows:bool=False,
)->dict|tuple[dict,pd.DataFrame]:
    """Full original GIFT island frame plus independent geological subset.

    Unlike former >=30-triple flora, ALL islands and ALL per-trait
    available source-Direct measurements enter the coverage denominator.
    Each island is represented ONCE in each trait/scope outcome.
    """
    if bootstrap_draws<999:
        raise ValueError("Original exploratory block bootstrap requires 999+ draws")
    needed={"island_id","accepted_species","origin_status",
            "endemic_status","floristic_status"}
    if not needed.issubset(flora) or flora.duplicated(
        ["island_id","accepted_species"]).any():
        raise ValueError("Original global island botanical flora identity not unique")
    if flora[["island_id","accepted_species"]].isna().any().any():
        raise ValueError("Original island/taxon membership blank")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original flora status not native/introduced/unresolved")
    if ((flora.floristic_status.eq("endemic") &
         (~flora.origin_status.eq("native") |
          ~flora.endemic_status.eq("endemic"))).any()):
        raise ValueError("Original island endemic and native status contradict")
    states=strict_single_trait_original(ledger)
    geokeys={"island_id","geological_origin","analysis_regime","spatial_block"}
    if not geokeys.issubset(geology) or geology.island_id.duplicated().any():
        raise ValueError("Original independent Roeble/Zizka GSHHG geography invalid")
    # The independent GSHHG crosswalk intentionally excludes a continental
    # split/Eurasian component, while the original plant membership may still
    # contain that source key. NEVER drop the island/plant records just
    # because geological classification or geography is unavailable.
    missing_geometry=set(flora.island_id)-set(geology.island_id)
    if len(missing_geometry)>max(30,.01*flora.island_id.nunique()):
        raise ValueError("Too many original flora islands lack independent GSHHG geography")
    if not {"island_id","source_mode","source_rank","entity_ID"}.issubset(assignments):
        raise ValueError("Original five-source GIFT native candidate IDs absent")
    nearby=assignments.loc[assignments.source_mode.eq("geo_k5")].copy()
    nearby["source_rank"]=pd.to_numeric(nearby.source_rank,errors="raise")
    nearby=nearby.loc[nearby.source_rank.between(1,5)]
    if nearby.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Source original island rank and region duplicate")
    near=nearby.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if not {"entity_ID","work_species"}.issubset(mainland):
        raise ValueError("GIFT native mainland source botanical checklist missing")
    native=mainland[["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    geo=geology.set_index("island_id")
    n_original=len(set(flora.island_id))
    if n_original<4000:
        raise ValueError("This is supposed to be global mainland-flora universe; not 928")
    if sum(geo.geological_origin.isin(GEOGRAPHIC_ORIGINS))<400:
        raise ValueError("Independent origin-class source sample too small")
    all_traits={t:set(q) for t,q in states.items()}
    gplots=[]
    completeness=[]
    all_source_metadata=[]
    for island,original in flora.groupby("island_id",sort=True):
        island=str(island)
        own=set(original.accepted_species)
        if island in geo.index:
            g=geo.loc[island]
            tropical_region=str(g["analysis_regime"])
            block=str(g["spatial_block"])
            geology_origin=str(g["geological_origin"])
        else:
            # A missing geomorphic source is NOT botanical trait absence.
            tropical_region="geography_unresolved"
            block="unclassified_GSHHG_source"
            geology_origin="unresolved"
        g5=near.get(island,())
        five_source=(len(g5)==5 and len(set(g5))==5)
        # This candidate source GIFT region is NOT the historical ancestor.
        source=(set().union(*(native.get(x,set()) for x in g5))-own
                if five_source and geology_origin in GEOGRAPHIC_ORIGINS else set())
        if source & own:
            raise AssertionError("Island accepted source plant leaked into GIFT candidate mainland")
        all_source_metadata.append(dict(
            island_id=island,analysis_regime=tropical_region,origin=geology_origin,
            source_5_available=bool(five_source),n_total_recorded_flora=len(own)))
        complete=len(own&set(states["flower_primary_color"])&
                     set(states["floral_form"])&set(states["floral_symmetry"]))
        completeness.append(dict(island_id=island,
           n_complete_direct_three_traits=int(complete)))
        for scope,local in (
            ("all_recorded",own),
            ("original_native",set(original.loc[
                original.origin_status.eq("native"),"accepted_species"]))):
            n_scope=len(local)
            for trait in TRAITS:
                own_fd=single_trait_fd(states[trait],local)
                row=dict(island_id=island,regime=tropical_region,
                    geological_origin=geology_origin,original_scope=scope,
                    spatial_block=block,
                    source_signature=g5[0] if five_source else None,
                    outcome=trait,outcome_type="single",
                    n_original_flora_species=n_scope,
                    observed_direct_coverage=float(own_fd["n_direct_species"]/n_scope)
                        if n_scope else float("nan"),
                    n_recorded_direct_species=own_fd["n_direct_species"],
                    n_recorded_genera=own_fd["n_distinct_recorded_genera"],
                    source_available=five_source,
                    observed_finite5=own_fd["finite_5_species_rao"],
                    observed_genus_equal_finite5=own_fd["finite_5_genus_equal_rao"],
                    observed_two_trait_pairing=float("nan"),
                    n_source_direct_species=0,
                    source_relative_finite5=float("nan"),
                    source_relative_genus_equal_finite5=float("nan"),
                    source_relative_pairing=float("nan"),
                )
                if source and n_scope>=N_RAO:
                    src=single_trait_fd(states[trait],source)
                    row["n_source_direct_species"]=src["n_direct_species"]
                    if (own_fd["n_direct_species"]>=N_RAO and
                        src["n_direct_species"]>=N_RAO):
                        row["source_relative_finite5"]=float(
                            own_fd["finite_5_species_rao"]-
                            src["finite_5_species_rao"])
                    if (np.isfinite(own_fd["finite_5_genus_equal_rao"]) and
                        np.isfinite(src["finite_5_genus_equal_rao"])):
                        row["source_relative_genus_equal_finite5"]=float(
                            own_fd["finite_5_genus_equal_rao"]-
                            src["finite_5_genus_equal_rao"])
                gplots.append(row)
            for first,second in PAIR_TRAITS:
                local_fd=pair_trait_fd(states[first],states[second],local)
                tag=first+"__"+second
                row=dict(island_id=island,regime=tropical_region,
                    geological_origin=geology_origin,original_scope=scope,
                    spatial_block=block,
                    source_signature=g5[0] if five_source else None,
                    outcome=tag,outcome_type="pair",
                    n_original_flora_species=n_scope,
                    observed_direct_coverage=float(
                        local_fd["n_two_trait_original_species"]/n_scope)
                        if n_scope else float("nan"),
                    n_recorded_direct_species=local_fd["n_two_trait_original_species"],
                    n_recorded_genera=None,source_available=five_source,
                    observed_finite5=local_fd["finite_5_joint_rao"],
                    observed_genus_equal_finite5=float("nan"),
                    observed_two_trait_pairing=local_fd["finite_5_pairing_departure"],
                    n_source_direct_species=0,
                    source_relative_finite5=float("nan"),
                    source_relative_genus_equal_finite5=float("nan"),
                    source_relative_pairing=float("nan"))
                if source and n_scope>=N_RAO:
                    src=pair_trait_fd(states[first],states[second],source)
                    row["n_source_direct_species"]=src["n_two_trait_original_species"]
                    if (local_fd["n_two_trait_original_species"]>=N_RAO and
                        src["n_two_trait_original_species"]>=N_RAO):
                        row["source_relative_finite5"]=float(
                            local_fd["finite_5_joint_rao"]-src["finite_5_joint_rao"])
                        row["source_relative_pairing"]=float(
                            local_fd["finite_5_pairing_departure"]-
                            src["finite_5_pairing_departure"])
                gplots.append(row)
    data=pd.DataFrame(gplots)
    comp=pd.DataFrame(completeness)
    if len(data)!=n_original*len(SCOPE)*(len(TRAITS)+len(PAIR_TRAITS)):
        raise ValueError("Original all-island trait-specific records missing")
    n_legacy=int((comp.n_complete_direct_three_traits>=30).sum())
    if n_legacy!=928:
        raise ValueError("Original historic >=30 species three-trait 928-island cohort not preserved")
    # Summary group denominators are ALL island units per region/scope,
    # even if no source direct trait observation. Zero is never used for
    # missing botanical evidence or unestimable diversity.
    grouped=[]
    for (regime,origin,scope,kind,trait),q in data.groupby(
        ["regime","geological_origin","original_scope",
         "outcome_type","outcome"],sort=True):
        good=q.loc[q.n_recorded_direct_species.ge(N_RAO)]
        src=q.loc[q.source_relative_finite5.notna()]
        row=dict(region=regime,geological_origin=origin,
            source_plant_status_scope=scope,
            endpoint_type=kind,trait_name=trait,
            n_all_original_flora_islands=int(len(q)),
            n_islands_with_any_original_direct_trait=int(
                q.n_recorded_direct_species.gt(0).sum()),
            n_islands_with_at_least_two_direct_species=int(
                q.n_recorded_direct_species.ge(2).sum()),
            n_islands_with_at_least_five_direct_species=int(len(good)),
            n_islands_with_at_least_ten_direct_species=int(
                q.n_recorded_direct_species.ge(10).sum()),
            n_islands_with_at_least_thirty_direct_species=int(
                q.n_recorded_direct_species.ge(30).sum()),
            n_source_GIFT_native_regions_five_available=int(
                q.source_available.sum()),
            n_original_source_relative_islands_estimate=int(len(src)),
            n_blocks_source_relative=int(src.spatial_block.nunique()),
            n_first_GIFT_candidate_source_groups=int(
                src.source_signature.nunique()),
            total_observed_island_times_species_Direct=int(
                q.n_recorded_direct_species.sum()),
            mean_individual_island_trait_direct_coverage_fraction=float(
                (q.n_recorded_direct_species/
                 q.n_original_flora_species.replace(0,np.nan)).mean())
                if q.n_original_flora_species.gt(0).any() else None,
            median_observed_direct_species_per_original_island=float(
                q.n_recorded_direct_species.median()),
            median_direct_species_among_estimable_islands=float(
                good.n_recorded_direct_species.median()) if len(good) else None,
        )
        # In general there is no physical/source independence gate in
        # global source-absolute trait observations, so avoid causal CI.
        if len(good):
            row["raw_observed_island_equal_finite5_mean"]=float(
                good.observed_finite5.mean())
            row["n_distinct_blocks_usable_raw"]=int(
                good.spatial_block.nunique())
        if len(src):
            row["source_relative_observed_island_equal_finite5_mean"]=float(
                src.source_relative_finite5.mean())
        grouped.append(row)
    conditional=[]
    geology_data=data.loc[data.geological_origin.isin(GEOGRAPHIC_ORIGINS)].copy()
    for (regime,scope,kind,trait),q in geology_data.groupby(
        ["regime","original_scope","outcome_type","outcome"],sort=True):
        a=q.loc[q.geological_origin.eq("oceanic")]
        b=q.loc[q.geological_origin.eq("continental")]
        n_a=int(a.island_id.nunique())
        n_b=int(b.island_id.nunique())
        for metric in ("observed_direct_coverage","observed_finite5","observed_genus_equal_finite5",
                       "observed_two_trait_pairing","source_relative_finite5",
                       "source_relative_genus_equal_finite5","source_relative_pairing"):
            x=a.loc[a[metric].notna()]
            y=b.loc[b[metric].notna()]
            blocks=[x.spatial_block.nunique(),y.spatial_block.nunique()]
            sources=[x.source_signature.nunique(),y.source_signature.nunique()]
            # Raw metric also requires physical and GIFT-source cluster
            # opportunity for conservative comparable-block claims.
            report=(min(len(x),len(y))>=10 and min(blocks)>=5
                    and min(sources)>=5)
            result=dict(
                region=regime,source_plant_status_scope=scope,
                endpoint_type=kind,trait_name=trait,
                metric=metric,
                n_original_oceanic_island_floras=n_a,
                n_original_continental_island_floras=n_b,
                n_qualified_oceanic_islands=int(len(x)),
                n_qualified_continental_islands=int(len(y)),
                n_qualified_spatial_blocks_oceanic=int(blocks[0]),
                n_qualified_spatial_blocks_continental=int(blocks[1]),
                n_qualified_GIFT_source_groups_oceanic=int(sources[0]),
                n_qualified_GIFT_source_groups_continental=int(sources[1]),
                conditional_sampling_support=bool(report),
                primary_coverage_or_source_gate="per-trait_n5_10_island_5block_5source",
            )
            if report:
                for label,col in (
                    ("physical_block","spatial_block"),
                    ("GIFT_first_source","source_signature")):
                    xa=x.groupby(col)[metric].mean().dropna().to_numpy(float)
                    yb=y.groupby(col)[metric].mean().dropna().to_numpy(float)
                    if min(len(xa),len(yb))<5:
                        continue
                    rng=np.random.default_rng(_seed(
                        seed,"origin",regime,scope,trait,metric,label))
                    bs_a=xa[rng.integers(len(xa),size=(
                        bootstrap_draws,len(xa)))].mean(axis=1)
                    bs_b=yb[rng.integers(len(yb),size=(
                        bootstrap_draws,len(yb)))].mean(axis=1)
                    val=bs_a-bs_b
                    result[label+"_oceanic_minus_continental"]=dict(
                        mean=float(xa.mean()-yb.mean()),
                        ci95=[float(w) for w in np.quantile(val,[.025,.975])],
                        n_independent_origin_clusters=[len(xa),len(yb)],
                    )
            conditional.append(result)
    counts={}
    for trait,x in states.items():
        counts[trait]=int(len(x))
    original_islands_count=(
        data[["island_id","regime","geological_origin"]].drop_duplicates())
    out=dict(
        status="GLOBAL_GIFT_FLORA_PARTIAL_DIRECT_TRAITS_NO_COMPLETE_CASE",
        n_original_global_distinct_island_floras=n_original,
        n_original_flora_islands_without_independent_GSHHG_geography=len(missing_geometry),
        n_original_distinct_islands_with_geological_origin_typed=int(
            original_islands_count.geological_origin.isin(
                GEOGRAPHIC_ORIGINS).sum()),
        n_independently_geologically_classified_GSHHG_islands=int(
            geo.geological_origin.isin(GEOGRAPHIC_ORIGINS).sum()),
        n_original_three_trait_complete30_islands_retained_for_SENSITIVITY=n_legacy,
        n_original_direct_plant_species_per_trait=counts,
        n_original_partial_island_scope_endpoint_rows=int(len(data)),
        n_original_complete_case_three_trait_sensitivity_islands=928,
        marginal_trait_and_pairwise_coverage=grouped,
        oceanic_minus_continental_conditional_traits=conditional,
        source_design=dict(
            original_island_complete_case_gate_removed_for_main_analysis=True,
            full_original_island_flora_not_just_928=True,
            original_native_and_all_recorded_outcomes_always_separate=True,
            all_10_original_Direct_traits_enter_as_available=True,
            original_pairwise_joint_Direct_evidence_needs_exactly_two_traits=True,
            no_three_trait_required_for_pairwise_trait_association=True,
            no_imputation_or_single_trait_missing_as_absence=True,
            direct_source_species_unique_high_medium_multistate_excluded=True,
            finite_n5_rao_same_source_sample_size_without_repeated_rarefaction=True,
            original_GIFT_nearest_five_native_mainland_not_true_ancestry=True,
            original_full_flora_taxa_removed_from_candidate_mainland_sources=True,
            genus_equal_weighting_does_not_control_lineage_history=True,
            missingness_could_be_taxonomically_and_geographically_nonrandom=True,
            original_direct_trait_observation_fraction_by_island_and_scope=True,
            native_status_scope_without_original_native_species_is_NA_not_zero=True,
            independent_block_and_nearest_GIFT_source_group_CIs_unadjusted_for_many_tests=True,
            historical_928_three_trait_sample_is_sensitivity_ONLY=True,
            no_pollinator_effect_or_adaptive_evolution_observed=True,
            source_species_island_geography_IDs_never_public=True,
        )
    )
    # Available ONLY to rights-protected original-source local CI workflows.
    # This DataFrame contains island IDs and never belongs in JSON artifacts.
    if return_private_island_rows:
        return out,data
    return out
