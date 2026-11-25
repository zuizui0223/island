"""Oceanic vs continental: observed floral establishment filtering vs divergence proxy.

The source floristic 'availability' is a checklist, NOT transported seed rain
or a true dispersal/establishment experiment. Island endemic recorded floral
distance from potential source congeners is NOT proven local adaptation.

Geological labels are independently audited Roeble 2024/Zizka 2022 labels,
not inferred from absence of LGM mainland connection. All outcomes are
island/species-direct category evidence; no hypothetical flower profiles.
"""
from __future__ import annotations
from collections import defaultdict
from hashlib import sha256

import numpy as np
import pandas as pd

from .regional_native_flora_fd import TRAITS, METRICS, exact_rarefied_flora_fd
from .regional_native_provenance_volume import cluster_equal_bootstrap
from .native_taxon_matched_pairing import species_genus
from .geological_disjoint_matched import matched_oceanic_continental_source_audit

ORIGINS=("continental","oceanic")
SCOPES=("all_recorded","regional_WCVP_compatible","strict_source_native")
RESULT_METRICS=tuple(METRICS)+(
    "shared_source_genus_fraction",
    "source_genus_explains_local_fraction",
    "source_genus_jaccard",
    "source_genus_sorensen",
    "expected_10genus_share",
    "source_plant_genus_pair_match_probability",
)
GROUP_MIN=(10,5,5)


def _seed(*parts:object)->int:
    return int.from_bytes(sha256(":".join(map(str,parts)).encode()).digest()[:8],"big")


def _dist(a:tuple[str,str,str],b:tuple[str,str,str])->int:
    return sum(x!=y for x,y in zip(a,b))


def candidate_endemic_divergence(
    island_rows:pd.DataFrame,local:set[str],source:set[str],
    original_profiles:dict[str,tuple[str,str,str]],
)->dict:
    """Within-source-genus floral distance of sampled island taxon names.

    Separate origin='endemic' and origin='native' labels. Neither source
    congeners nor native nonendemic plants are reconstructed ancestors.
    """
    donor_by_genus=defaultdict(list)
    for x in source:
        donor_by_genus[species_genus(x)].append(x)
    required={"accepted_species","origin_status","endemic_status","floristic_status"}
    if not required.issubset(island_rows):
        raise ValueError("Original separate GIFT endemic-status evidence missing")
    original=island_rows[list(required)].copy()
    if original.endemic_status.isna().any() or not original.endemic_status.isin(
        ("endemic","nonendemic","unresolved")).all():
        raise ValueError("Original endemicity unknown or invalid source category")
    if ((original.floristic_status.eq("endemic") &
         (~original.origin_status.eq("native") |
          ~original.endemic_status.eq("endemic"))).any() or
        (original.floristic_status.eq("native_nonendemic") &
         (~original.origin_status.eq("native") |
          ~original.endemic_status.eq("nonendemic"))).any()):
        raise ValueError("Separate original origin and endemism evidence conflicts")
    eligible=island_rows.loc[island_rows.accepted_species.isin(local)]
    recorded={"endemic":[],"native":[]}
    totals={"endemic":0,"native":0}
    for taxon,status in zip(eligible.accepted_species,eligible.floristic_status):
        # Endemic is a separate source layer; origin='native' alone is NOT
        # evidence of either endemic or non-endemic status.
        category={"endemic":"endemic","native_nonendemic":"native"}.get(status)
        if category is None:
            continue
        totals[category]+=1
        candidates=donor_by_genus.get(species_genus(taxon),())
        if candidates:
            d=min(_dist(original_profiles[taxon],original_profiles[y])
                  for y in candidates)
            recorded[category].append(d)
    out={}
    for k,ds in recorded.items():
        out["n_"+k+"_original_direct"]=totals[k]
        out["n_"+k+"_source_genus_comparable"]=len(ds)
        out[k+"_source_congener_category_mismatch_mean"]=(
            float(np.mean(ds)) if ds else np.nan)
    if min(map(len,recorded.values()))>=3:
        out["endemic_minus_native_source_congener_mismatch"]=(
            out["endemic_source_congener_category_mismatch_mean"]-
            out["native_source_congener_category_mismatch_mean"])
    else:
        out["endemic_minus_native_source_congener_mismatch"]=np.nan
    return out


def geology_origin_flora_pilot(
    direct:pd.DataFrame, classified:pd.DataFrame, geography:pd.DataFrame,
    independent_geology:pd.DataFrame, nearby:dict[str,tuple[str,...]],
    donors:dict[str,set[str]],*, seed:int=20261114, boots:int=1999
)->dict:
    """Calculates only non-identifying original island-type regional aggregates.

    No interpretation of source genera as seed arrivals, and no inferred
    evolutionary selection from endemic flower-category distances.
    """
    if len(direct)!=2234 or list(direct.columns)!=list(TRAITS):
        raise ValueError("Original 2234 exact source Direct species changed")
    if classified.island_id.nunique()!=928:
        raise ValueError("Original 928 complete island floras changed")
    if not {"island_id","accepted_species","origin_status",
            "strict_source_native","regional_native_compatible"}.issubset(classified):
        raise ValueError("Original botanical origin information incomplete")
    if not {"island_id","geological_origin","weigelt_archipelago",
            "geological_class_source","LGM_connectivity"}.issubset(independent_geology):
        raise ValueError("Independent geological label source not verified")
    if not {"island_id","analysis_regime","spatial_block"}.issubset(geography):
        raise ValueError("Original geographic strata incomplete")
    if independent_geology.island_id.duplicated().any():
        raise ValueError("Repeated geological source island crosswalk")
    if geography.island_id.duplicated().any() or len(geography)!=928:
        raise ValueError("Original source island geography collision")
    if classified.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original source island trait rows duplicated")
    if len(nearby)!=928:
        raise ValueError("GIFT five-source assignment incomplete")
    validate=independent_geology.geological_origin.value_counts().to_dict()
    if min(validate.get(v,0) for v in ORIGINS)<150:
        raise ValueError("Independent Roeble/Zizka geological sample smaller than archived minimum")
    # Join geology ONLY by original island_id; not region centroid approximation here.
    g=geography.merge(
        independent_geology[["island_id","geological_origin","LGM_connectivity",
              "weigelt_archipelago","geological_class_source",
              "area_km2","distance_to_continent_km",
              "climate_pc1","climate_pc2","climate_pc3","climate_pc4"]],
        on="island_id",validate="one_to_one",how="left")
    g=g.set_index("island_id")
    profiles={str(k):tuple(map(str,row))
              for k,row in zip(direct.index,direct.to_numpy(dtype=str))}
    records=[]
    for island,flora in classified.groupby("island_id",sort=True):
        island=str(island)
        if g.at[island,"geological_origin"] not in ORIGINS:
            continue
        if (len(nearby[island])!=5 or len(set(nearby[island]))!=5):
            raise ValueError("Original five GIFT candidate mainland regions missing")
        roster=set(flora.accepted_species)
        source=(set().union(*(donors.get(t,set()) for t in nearby[island]))&
                set(profiles))-roster
        if len(source)<30:
            continue
        source_states=np.asarray([profiles[t] for t in sorted(source)],dtype=object)
        source_score=exact_rarefied_flora_fd(source_states,n=30)
        sg={species_genus(t) for t in source}
        source_count={}
        for t in source:
            genus=species_genus(t)
            source_count[genus]=source_count.get(genus,0)+1
        endemic=candidate_endemic_divergence(flora,roster,source,profiles)
        for scope,names in (
            ("all_recorded",roster),
            ("regional_WCVP_compatible",set(flora.loc[
                flora.regional_native_compatible,"accepted_species"])),
            ("strict_source_native",set(flora.loc[
                flora.strict_source_native,"accepted_species"])),
        ):
            if len(names)<30:continue
            score=exact_rarefied_flora_fd(
                np.asarray([profiles[t] for t in sorted(names)],dtype=object),n=30)
            gn={species_genus(t) for t in names}
            local_count={}
            for t in names:
                genus=species_genus(t)
                local_count[genus]=local_count.get(genus,0)+1
            shared=len(gn & sg)
            # Exactly ten distinct genera drawn independently without
            # replacement from each island and source genus roster:
            # E[#shared]/10 = 10*|intersection|/(|local|*|source|).
            expected_ten=(10*shared/(len(gn)*len(sg))
                          if min(len(gn),len(sg))>=10 else float("nan"))
            plant_genus_pair_probability=sum(
                (local_count[gx]/len(names))*(source_count[gx]/len(source))
                for gx in gn & sg)
            row=dict(
                island_id=island,
                origin=g.at[island,"geological_origin"],
                regime=g.at[island,"analysis_regime"],
                spatial_block=g.at[island,"spatial_block"],
                source_signature=nearby[island][0],
                archipelago=g.at[island,"weigelt_archipelago"],
                geology_source=g.at[island,"geological_class_source"],
                LGM_connectivity=g.at[island,"LGM_connectivity"],
                scope=scope,n_original_direct_species=len(names),
                n_donor_direct_species=len(source),
                n_local_distinct_genera=len(gn),n_source_genera=len(sg),
                n_shared_genera=shared,
                area_km2=g.at[island,"area_km2"],
                distance_to_continent_km=g.at[island,"distance_to_continent_km"],
                climate_pc1=g.at[island,"climate_pc1"],
                climate_pc2=g.at[island,"climate_pc2"],
                climate_pc3=g.at[island,"climate_pc3"],
                climate_pc4=g.at[island,"climate_pc4"],
                shared_source_genus_fraction=shared/len(sg),
                source_genus_explains_local_fraction=shared/len(gn),
                source_genus_jaccard=shared/(len(sg|gn)),
                source_genus_sorensen=2*shared/(len(sg)+len(gn)),
                expected_10genus_share=expected_ten,
                source_plant_genus_pair_match_probability=plant_genus_pair_probability,
            )
            for metric in METRICS:
                row["source_relative_"+metric]=score[metric]-source_score[metric]
            if scope=="all_recorded":
                row.update(endemic)
            records.append(row)
    if not records:
        raise ValueError("No geological-class original 3-trait intersections")
    data=pd.DataFrame(records)
    aggregate=[]
    for (regime,scope,origin),q in data.groupby(
        ["regime","scope","origin"],sort=True
    ):
        r=dict(regime=regime,scope=scope,geological_origin=origin,
               n_geology_source_comparable_islands=len(q),
               n_original_physical_spatial_blocks=q.spatial_block.nunique(),
               n_first_nearest_GIFT_source_signatures=q.source_signature.nunique(),
               n_distinct_labelled_archipelagos=q.archipelago.dropna().replace("",np.nan).nunique(),
               n_Roeble_geological_source_islands=int(q.geology_source.eq("Roeble_2024_510_island_geology").sum()),
               n_Zizka_supplement_geological_source_islands=int(q.geology_source.eq("Zizka_2022_missing_geology_additions").sum()),
               median_direct_flora_species=float(q.n_original_direct_species.median()),
               median_source_direct_species=float(q.n_donor_direct_species.median()),
               median_local_distinct_genera=float(q.n_local_distinct_genera.median()),
               median_source_distinct_genera=float(q.n_source_genera.median()),
               median_shared_genera=float(q.n_shared_genera.median()))
        enabled=(len(q)>=GROUP_MIN[0] and
                q.spatial_block.nunique()>=GROUP_MIN[1] and
                q.source_signature.nunique()>=GROUP_MIN[2])
        r["within_geology_descriptive_gate"]=bool(enabled)
        if enabled:
            for metric in RESULT_METRICS:
                key="source_relative_"+metric if metric in METRICS else metric
                for label,cluster in (("block","spatial_block"),("source","source_signature")):
                    r[label+"_"+metric]=cluster_equal_bootstrap(
                        q,key,cluster,draws=boots,
                        seed=_seed(seed,"group",regime,scope,origin,metric,label))
        if scope=="all_recorded":
            r["n_source_endemic_direct_species"]=int(
                q.n_endemic_original_direct.sum())
            r["n_source_native_direct_species"]=int(
                q.n_native_original_direct.sum())
            r["n_islands_with_any_Direct_endemic"]=int(
                q.n_endemic_original_direct.gt(0).sum())
            r["n_source_congener_comparable_endemic_species"]=int(
                q.n_endemic_source_genus_comparable.sum())
            r["n_source_congener_comparable_native_species"]=int(
                q.n_native_source_genus_comparable.sum())
            r["n_islands_with_3_endemic_and_3_native_congeneric_source_comparisons"]=int(
                q.endemic_minus_native_source_congener_mismatch.notna().sum())
            qq=q.loc[q.endemic_minus_native_source_congener_mismatch.notna()]
            r["n_native_endemic_paired_islands_within_geology"]=len(qq)
            if len(qq)>=10 and qq.spatial_block.nunique()>=5 and qq.source_signature.nunique()>=5:
                for weight,key in (("block","spatial_block"),("source","source_signature")):
                    r[weight+"_endemic_minus_native_original_congener_mismatch"]=cluster_equal_bootstrap(
                        qq,"endemic_minus_native_source_congener_mismatch",key,
                        draws=boots,seed=_seed(seed,"endemic",regime,origin,weight))
        aggregate.append(r)
    between=[]
    for (regime,scope),q in data.groupby(["regime","scope"],sort=True):
        counts=q.origin.value_counts().to_dict()
        c=q.loc[q.origin.eq("continental")]
        o=q.loc[q.origin.eq("oceanic")]
        valid=min(counts.get("oceanic",0),counts.get("continental",0))>=10
        for metric in RESULT_METRICS:
            k="source_relative_"+metric if metric in METRICS else metric
            row=dict(regime=regime,scope=scope,response=metric,
                     n_oceanic_islands=len(o),n_continental_islands=len(c),
                     n_oceanic_blocks=o.spatial_block.nunique(),
                     n_continental_blocks=c.spatial_block.nunique(),
                     n_oceanic_source_signatures=o.source_signature.nunique(),
                     n_continental_source_signatures=c.source_signature.nunique())
            row["origin_contrast_supported_sampling"]=bool(valid and
                min(row["n_oceanic_blocks"],row["n_continental_blocks"])>=5 and
                min(row["n_oceanic_source_signatures"],
                    row["n_continental_source_signatures"])>=5)
            if row["origin_contrast_supported_sampling"]:
                # Both directions summarized as group means; bootstrap within
                # each origin and source/block cluster independently.
                for label,cluster in (("block","spatial_block"),
                                      ("source","source_signature")):
                    a=c.groupby(cluster)[k].mean().dropna().to_numpy(float)
                    b=o.groupby(cluster)[k].mean().dropna().to_numpy(float)
                    if min(len(a),len(b))<5:
                        row[label+"_oceanic_minus_continental"]=None
                        continue
                    rng=np.random.default_rng(_seed(seed,"contrast",regime,scope,metric,label))
                    contrast=b.mean()-a.mean()
                    boots_o=b[rng.integers(len(b),size=(boots,len(b)))].mean(axis=1)
                    boots_c=a[rng.integers(len(a),size=(boots,len(a)))].mean(axis=1)
                    bs=boots_o-boots_c
                    row[label+"_oceanic_minus_continental"]=dict(
                        mean=float(contrast),
                        ci95=[float(v) for v in np.quantile(bs,[.025,.975])],
                        n_geographic_or_source_clusters_by_origin=[len(b),len(a)])
            between.append(row)
    return dict(
        status="GEOLOGICAL_OCEANIC_CONTINENTAL_OBSERVED_FLORAL_ASSEMBLY_AND_ENDEMIC_PROXY",
        n_original_all_928_Direct_island_floras=928,
        original_global_direct_species_complete=2234,
        independent_archived_geology_label_counts=validate,
        n_geology_labelled_source_comparable_all_flora=data.loc[
            data.scope.eq("all_recorded"),"island_id"].nunique(),
        n_endemic_original_direct_species_with_some_source_congener=int(
            data.loc[data.scope.eq("all_recorded"),"n_endemic_source_genus_comparable"].sum()
        ),
        analysis=dict(
            geography_ancillary_Roeble_vs_Zizka_geology_not_LGM=True,
            candidate_GIFT_native_mainland_checklist_not_actual_propagule_source=True,
            establishment_survival_not_observable_from_species_checklists=True,
            genus_composition_and_trait_Direct_coverage_confounding_remain=True,
            WCVP_native_range_not_island_local_native_status=True,
            endemic_source_congener_comparison_is_not_adaptation_or_sister_phylogeny=True,
            original_source_taxon_names_removed_from_all_donors=True,
            entire_original_10degree_block_and_first_source_weighting_separate=True,
            correlated_nested_posthoc_estimands_not_independent_replicates=True,
            no_historical_guild_loss_or_plant_reproductive_fitness_inferred=True,
        ),
        groups=aggregate,
        geological_origin_associations=between,
        geology_matched_geography_and_archipelago_diagnostics=(
            matched_oceanic_continental_source_audit(data)),
        no_raw_island_species_or_geological_row_exported=True,
    )
