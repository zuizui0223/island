"""Global GIFT botanical flora × five GBIF guild co-occurrence source-power audit.

Retains all 4,505 original island plant floras, never requires three
floral states jointly observed. GBIF 'retained' means a qualified animal
occurrence, not a documented visit to the same plant. 'disrupted'
means qualified source-aware NONDETECTION, not historical local extinction.
If no channel row, mark unresolved—not a zero detection.

This is a source-coverage and candidate-guild coexistence power audit.
Neither a pollination shift, replacement, functional compensation, nor
seed-set effectiveness is identified by animal occurrence states.
No private accepted plants, guild taxon names, island GSHHG identifiers,
GIFT source geography or original animal observation rows are exported.
"""
from __future__ import annotations

import pandas as pd

from .global_partial_pollinator_guild_bridge import (
    CHANNELS, normalize_guild_evidence
)
from .global_partial_trait_flora_fd import (
    TRAITS, PAIR_TRAITS, strict_single_trait_original
)

TARGETS=(
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "autonomous_selfing_capacity",
    "flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
)


def guild_cooccurrence_global(
    original_flora:pd.DataFrame,ledger:pd.DataFrame,
    independent_geo:pd.DataFrame,
    original_channels:pd.DataFrame
)->dict:
    fields={"island_id","accepted_species","origin_status"}
    if not fields.issubset(original_flora) or original_flora.duplicated(
        ["island_id","accepted_species"]).any():
        raise ValueError("Frozen global GIFT native/all original island roster corrupted")
    if not original_flora.origin_status.isin(
        ("native","introduced","unresolved")).all():
        raise ValueError("Original GIFT native origin should not be inferred")
    meta_fields={"island_id","analysis_regime","spatial_block","geological_origin"}
    if not meta_fields.issubset(independent_geo) or (
        independent_geo.island_id.duplicated().any()):
        raise ValueError("Published independent geology and original island geo invalid")
    n=int(original_flora.island_id.nunique())
    if n<4000:
        raise ValueError("Global trait denominator collapsed into old 928 complete case")
    states=strict_single_trait_original(ledger)
    # Construct source opportunity sets ONCE, not 4505 times per trait.
    names_by_trait={key:set(src) for key,src in states.items()}
    names_by_pair={a+"__"+b:names_by_trait[a]&names_by_trait[b]
                   for a,b in PAIR_TRAITS}
    keys=original_flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=original_flora.loc[
        original_flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    meta=pd.DataFrame({"island_id":sorted(keys)})
    geo=meta.merge(
        independent_geo[list(meta_fields)],how="left",on="island_id",
        validate="one_to_one")
    geo["analysis_regime"]=geo.analysis_regime.fillna("geography_unresolved")
    geo["spatial_block"]=geo.spatial_block.fillna("unknown_geography")
    geo["geological_origin"]=geo.geological_origin.fillna("unresolved")
    channels=normalize_guild_evidence(original_channels,set(keys))
    wide=channels.pivot(index="island_id",columns="channel_id",values="channel_state")
    wide=wide.loc[sorted(keys),list(CHANNELS)]
    if len(wide)!=n or wide.isna().any().any():
        raise AssertionError("Original island channel record absent not a qualified nondetection")
    cnt=[]
    for island in sorted(keys):
        full=keys[island]
        is_native=native.get(island,set())
        row=dict(island_id=island)
        for t in TRAITS:
            roster=names_by_trait[t]
            row["n_all_direct__"+t]=len(full&roster)
            row["n_native_direct__"+t]=len(is_native&roster)
        for first,second in PAIR_TRAITS:
            tag=first+"__"+second
            roster=names_by_pair[tag]
            row["n_all_direct__"+tag]=len(full&roster)
            row["n_native_direct__"+tag]=len(is_native&roster)
        cnt.append(row)
    evidence=geo.merge(pd.DataFrame(cnt),on="island_id",validate="one_to_one")
    evidence=evidence.set_index("island_id").join(wide)
    if len(evidence)!=n:
        raise AssertionError("Global original source flora lost")
    output=[]
    for region in ["GLOBAL_ALL",*sorted(
        str(r) for r in evidence.analysis_regime.unique()
        if r not in ("","geography_unresolved"))]:
        region_subset=(evidence if region=="GLOBAL_ALL" else
                       evidence.loc[evidence.analysis_regime.eq(region)])
        if region_subset.empty:
            continue
        for geological in ("all","oceanic","continental"):
            q=(region_subset if geological=="all" else
               region_subset.loc[region_subset.geological_origin.eq(geological)])
            if q.empty:continue
            detected=q[list(CHANNELS)].eq("retained")
            non=q[list(CHANNELS)].eq("disrupted")
            unresolved=q[list(CHANNELS)].eq("unresolved")
            record=dict(
                region=region,geological_origin=geological,
                n_ALL_islands_in_original_botanical_flora=int(len(q)),
                n_at_least_one_recorded_guild=int(detected.any(axis=1).sum()),
                n_at_least_two_recorded_guilds=int(detected.sum(axis=1).ge(2).sum()),
                n_all_five_guilds_qualified_detected_or_nondetected=int(
                    (detected|non).all(axis=1).sum()),
                n_any_source_qualified_guild_nondetection=int(non.any(axis=1).sum()),
                n_any_guild_unresolved=int(unresolved.any(axis=1).sum()),
                no_animal_occurrence_deemed_a_proven_pollinator_visit=True,
            )
            by=[]
            for absent in CHANNELS:
                for present in CHANNELS:
                    if absent==present:continue
                    mask=non[absent]&detected[present]
                    controls=detected[absent]&detected[present]
                    outcome=dict(
                        originally_qualified_NONDETECTION_guild=absent,
                        other_guild_GBIF_RECORDED=present,
                        n_both_source_qualified_guild_opportunity=int(mask.sum()),
                        n_islands_with_both_guilds_GBIF_detected=int(controls.sum()),
                        n_10degree_blocks_with_candidate=int(
                            q.loc[mask,"spatial_block"].nunique()),
                    )
                    for tag in TARGETS:
                        for sc,prefix in (("all_recorded","n_all_direct__"),
                                          ("original_native","n_native_direct__")):
                            qualified=(mask&q[prefix+tag].ge(5))
                            outcome[
                                "n_candidate_with_5_direct__"+sc+"__"+tag
                            ]=int(qualified.sum())
                    by.append(outcome)
            record["source_guild_nondetection_plus_other_recorded_pairs"]=by
            output.append(record)
    global_record=next(x for x in output if
        x["region"]=="GLOBAL_ALL" and x["geological_origin"]=="all")
    return dict(
        status="GLOBAL_4505_PARTIAL_FLORA_GBIF_GUILD_CO_OCCURRENCE_NOT_SHIFT",
        n_original_global_plant_flora_islands=n,
        n_original_guilds=len(CHANNELS),
        original_GBIF_evaluable_all_five_guilds=global_record[
            "n_all_five_guilds_qualified_detected_or_nondetected"],
        original_global_guild_co_observation_by_region_and_geology=output,
        original_analyzed_single_trait_and_pairwise_target_evidence=list(TARGETS),
        source_scope=dict(
            all_original_GIFT_islands_stay_in_denom_even_if_guild_unknown=True,
            independent_pollinator_GBIF_nondetection_requires_source_and_effort=True,
            source_animal_occurrence_does_not_establish_actual_flower_visitor=True,
            source_qualified_nondetection_NOT_historical_guild_loss=True,
            guild_co_recorded_not_effective_replacement=True,
            flower_color_pair_requires_only_two_measured_traits=True,
            plant_native_subset_independently_original_reviewed=True,
            observed_plant_species_with_fewer_than_five_trait_states_still_in_data=True,
            no_original_island_plant_taxon_GIFT_source_or_GBIF_animal_ids_exported=True,
        )
    )
