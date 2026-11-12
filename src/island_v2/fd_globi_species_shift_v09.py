"""Private original GloBI v0.9 EXACT plant-species flower interaction coverage.

The original source/island cohort is selected without outcome-based guild
matching. Global GloBI species-linked records and coordinates do NOT verify
a local island visit, successful pollination or historical pollinator shift.
Release only aggregate sample coverage; species/island/reference IDs stay local.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import gzip
from pathlib import Path
import numpy as np
import pandas as pd

from island_v2.fd_globi_provider_v09 import (
    RO, RAW, PIN, VERSION, val, channel, provider, reference, sha256_file,
)

TRAITS=("flower_primary_color","floral_form","floral_symmetry")
GUILDS=("bombus","non_bombus_bees","lepidoptera",
        "flower_visiting_birds","diptera")
GEO_FIELDS=("decimalLatitude","decimalLongitude","localityName",
            "eventDateUnixEpoch")


def select_native_source_names(
    direct_ledger:pd.DataFrame, flora:pd.DataFrame,
    channel_evidence:pd.DataFrame, source_assignments:pd.DataFrame,
    mainland_flora:pd.DataFrame
) -> tuple[set[str],dict]:
    """Select only source-GIFT-name matched native island Direct flowers.

    FIVE plausible GIFT mainland source regions are not historical ancestors.
    Bombus nondetection / bird occurrence is a candidate ISLAND filter, not a
    species-specific bird-for-bee pollination-shift observation.
    """
    required=[
      (direct_ledger,{"accepted_species","trait_name","normalized_value",
                      "evidence_scope","evidence_quality"}),
      (flora,{"island_id","accepted_species","origin_status"}),
      (channel_evidence,{"island_id","channel_id","channel_state","source_state"}),
      (source_assignments,{"island_id","source_mode","source_rank","entity_ID"}),
      (mainland_flora,{"entity_ID","work_species"}),
    ]
    for frame,columns in required:
        if columns-set(frame):
            raise ValueError("Original island/trait/source fields missing")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate island plant identity")
    d=direct_ledger.loc[
        direct_ledger.trait_name.isin(TRAITS) &
        direct_ledger.evidence_scope.eq("species_direct") &
        direct_ledger.evidence_quality.isin(("high","medium"))].copy()
    if d.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Conflicting species-direct original trait")
    d=d.loc[d.normalized_value.notna() &
            ~d.normalized_value.astype(str).str.contains("|",regex=False)]
    wide=d.pivot(index="accepted_species",columns="trait_name",
                 values="normalized_value")
    if not set(TRAITS).issubset(wide):
        raise ValueError("All 3 species-direct flower traits required")
    known=set(wide.dropna(subset=list(TRAITS)).index)
    df=flora.loc[flora.accepted_species.isin(known)]
    counts=df.groupby("island_id").accepted_species.nunique()
    eligible=set(counts[counts>=30].index)
    e=channel_evidence.loc[channel_evidence.island_id.isin(eligible)]
    if (e.duplicated(["island_id","channel_id"]).any() or
        e[["island_id","channel_id","channel_state","source_state"]]
          .isna().any().any()):
        raise ValueError("Missing or duplicate original pollinator states")
    if (e.channel_state.eq("disrupted") &
        e.source_state.ne("available")).any():
        raise ValueError("Unqualified adequate nondetection")
    states=e.pivot(index="island_id",columns="channel_id",
                   values="channel_state")
    if not {"bombus","flower_visiting_birds"}.issubset(states):
        raise ValueError("Required guild evidence missing")
    candidates=set(states.index[
        states.bombus.eq("disrupted") &
        states.flower_visiting_birds.eq("retained")])
    s=source_assignments.loc[
        source_assignments.island_id.isin(candidates) &
        source_assignments.source_mode.eq("geo_k5") &
        source_assignments.source_rank.between(1,5)]
    if (s.duplicated(["island_id","source_rank"]).any() or
        s.duplicated(["island_id","entity_ID"]).any()):
        raise ValueError("Duplicate candidate GIFT source rank")
    sr=s.groupby("island_id").entity_ID.agg(list).to_dict()
    if any(len(sr.get(i,[]))!=5 for i in candidates):
        raise ValueError("Five candidate mainland source regions required")
    gift=mainland_flora.loc[
        mainland_flora.entity_ID.isin(set(s.entity_ID))].dropna().drop_duplicates()
    by_gift=gift.groupby("entity_ID").work_species.agg(set).to_dict()
    pools={i:set().union(*(by_gift.get(e,set()) for e in es))
           for i,es in sr.items()}
    native=df.loc[df.island_id.isin(candidates) &
                  df.origin_status.eq("native")]
    matched=native.loc[
        [name in pools[island]
         for island,name in native[["island_id","accepted_species"]]
                           .itertuples(index=False,name=None)]]
    species=set(matched.accepted_species)
    return species,{
      "n_islands_with_30plus_species_direct_flower_traits":len(eligible),
      "n_candidate_islands_qualified_Bombus_nondetection_bird_recorded":len(candidates),
      "n_candidate_native_island_plant_rows":len(native),
      "n_native_same_name_GIFT_source_island_plant_rows":len(matched),
      "n_distinct_native_same_name_GIFT_source_plant_species":len(species),
      "n_distinct_native_same_name_GIFT_source_plant_genera":
          len({z.split()[0] for z in species}),
      "n_islands_with_native_same_name_GIFT_source_plant":
          int(matched.island_id.nunique()),
      "focal_species_names_exported":False,
    }


class ExactSpeciesGloBICoverage:
    """Incremental citation-backed species × visitor-guild evidence auditor.

    Coordinates can establish that SOME record was georeferenced, not that the
    visit was on the focal island (requires original island polygon overlay).
    """
    def __init__(self,names:set[str]):
        if not names:raise ValueError("No fixed GIFT same-species focal names")
        self.names=set(names)
        self.raw=Counter()
        self.units=set()
        self.providers_by_reference=defaultdict(set)
        self.geocoded_species=set()
        self.dated_species=set()
        self.available_geography_fields=set()
        self.available_event_time_fields=set()
        self.spatial_bins=defaultdict(set)

    def add_chunk(self,df:pd.DataFrame):
        if set(RAW)-set(df):
            raise ValueError("Pinned GloBI required original source fields missing")
        self.available_geography_fields.update(
            set(GEO_FIELDS[:3]).intersection(df.columns))
        self.available_event_time_fields.update(
            set(GEO_FIELDS[3:]).intersection(df.columns))
        self.raw["n_original_GloBI_rows_read"]+=len(df)
        sub=df.loc[df.interactionTypeId.fillna("").isin(RO)]
        fwd=sub.interactionTypeId.isin([
            k for k,v in RO.items() if v[1]=="target"
        ]) & sub.targetTaxonSpeciesName.isin(self.names)
        rev=sub.interactionTypeId.isin([
            k for k,v in RO.items() if v[1]=="source"
        ]) & sub.sourceTaxonSpeciesName.isin(self.names)
        subset=sub.loc[fwd|rev]
        self.raw["n_exact_species_predicate_rows"]+=len(subset)
        for row in subset.fillna("").to_dict("records"):
            v,p,etype=RO[val(row["interactionTypeId"])]
            if val(row.get(p+"TaxonKingdomName"))!="Plantae":
                self.raw["n_plant_kingdom_mismatch"]+=1
                continue
            plant=val(row.get(p+"TaxonSpeciesName"))
            if plant not in self.names:
                continue
            if not val(row.get(v+"TaxonSpeciesName")):
                self.raw["n_unresolved_visitor_species"]+=1
                continue
            guild=channel(row,v)
            if not guild:
                self.raw["n_visitor_outside_predeclared_guilds"]+=1
                continue
            self.raw["n_eligible_species_flower_relation_rows"]+=1
            lat=pd.to_numeric(row.get("decimalLatitude"),errors="coerce")
            lon=pd.to_numeric(row.get("decimalLongitude"),errors="coerce")
            if (np.isfinite(lat) and np.isfinite(lon) and
                -90<=lat<=90 and -180<=lon<=180):
                self.geocoded_species.add(plant)
                self.spatial_bins[plant].add(
                    (int(np.floor(lat*10)),int(np.floor(lon*10))))
                self.raw["n_georeferenced_species_relations"]+=1
            try:
                epoch=float(row.get("eventDateUnixEpoch",""))
                if np.isfinite(epoch):
                    self.dated_species.add(plant)
                    self.raw["n_species_relation_rows_with_event_date"]+=1
            except (TypeError,ValueError):
                pass
            prov=provider(row)
            ref=reference(row)
            if not prov:
                self.raw["n_relation_rows_without_provider"]+=1
            if not ref:
                self.raw["n_relation_rows_without_bibliography"]+=1
            if not prov or not ref:
                continue
            self.units.add((plant,guild,ref,prov,etype))
            self.providers_by_reference[ref].add(prov)

    def summarize(self)->dict:
        duplicate_refs={ref for ref,providers in self.providers_by_reference.items()
                        if len(providers)>1}
        filtered=[r for r in self.units if r[2] not in duplicate_refs]
        by_plant=defaultdict(set)
        by_guild=defaultdict(set)
        by_etype=defaultdict(set)
        by_provider=defaultdict(set)
        for plant,guild,ref,prov,etype in filtered:
            by_plant[plant].add(ref)
            by_guild[(plant,guild)].add(ref)
            by_etype[(plant,etype)].add(ref)
            by_provider[plant].add(prov)
        # A single archive-first source label may subsume many studies.
        # Quantify exactly how much source-backed species-level support
        # survives omitting the largest original provider LABEL; never
        # present this as independent-study holdout.
        provider_units=Counter(row[3] for row in filtered)
        largest_label=provider_units.most_common(1)[0][0] if provider_units else None
        after_units=[row for row in filtered if row[3]!=largest_label]
        after_by_plant=defaultdict(set)
        after_by_guild=defaultdict(set)
        after_by_type=defaultdict(set)
        for plant,guild,ref,prov,etype in after_units:
            after_by_plant[plant].add(ref)
            after_by_guild[(plant,guild)].add(ref)
            after_by_type[(plant,etype)].add(ref)
        by_guild_summary=[]
        for g in GUILDS:
            species={p:references for (p,gg),references
                     in by_guild.items() if gg==g}
            surviving={p:refs for (p,gg),refs in after_by_guild.items()
                       if gg==g}
            by_guild_summary.append({
                "guild":g,
                "n_candidate_species_with_cited_interaction":
                    sum(len(v)>=1 for v in species.values()),
                "n_candidate_species_with_3plus_distinct_bibliography":
                    sum(len(v)>=3 for v in species.values()),
                "n_candidate_species_with_cited_interaction_after_largest_archive_removed":
                    sum(len(v)>=1 for v in surviving.values()),
                "n_candidate_species_with_3plus_distinct_bibliography_after_largest_archive_removed":
                    sum(len(v)>=3 for v in surviving.values()),
            })
        return {
          "source_version":VERSION,"GloBI_v09_SHA256":PIN,
          "n_fixed_candidate_species":len(self.names),
          "n_species_with_at_least_one_exact_cited_flower_interaction":
              sum(len(v)>=1 for v in by_plant.values()),
          "n_species_with_3plus_exact_distinct_bibliography":
              sum(len(v)>=3 for v in by_plant.values()),
          "n_species_with_exact_cited_pollination_class":
              sum(len(v)>=1 for (p,t),v in by_etype.items()
                  if t=="pollination"),
          "n_species_with_exact_cited_flower_visit_class":
              sum(len(v)>=1 for (p,t),v in by_etype.items()
                  if t=="flower_visit"),
          "n_species_with_2plus_distinct_documented_guilds":
              len({p for p in by_plant
                   if sum((p,g) in by_guild for g in GUILDS)>=2}),
          "n_species_with_2plus_provider_archive_labels":
              sum(len(v)>=2 for v in by_provider.values()),
          "n_distinct_archive_first_source_labels":
              len(provider_units),
          "largest_archive_label_fraction_of_cited_species_guild_reference_units":
              (float(provider_units[largest_label]/len(filtered))
               if filtered else None),
          "n_species_with_exact_cited_interactions_after_largest_archive_removed":
              sum(len(v)>=1 for v in after_by_plant.values()),
          "n_species_with_3plus_distinct_bibliography_after_largest_archive_removed":
              sum(len(v)>=3 for v in after_by_plant.values()),
          "n_species_with_exact_pollination_class_after_largest_archive_removed":
              sum(len(v)>=1 for (p,t),v in after_by_type.items()
                  if t=="pollination"),
          "n_species_with_2plus_guilds_after_largest_archive_removed":
              len({p for p in after_by_plant
                   if sum((p,g) in after_by_guild for g in GUILDS)>=2}),
          "archive_level_provider_is_not_proven_independent_study":True,
          "n_species_with_any_georeferenced_relation":
              len(self.geocoded_species),
          "n_species_with_2plus_distinct_coarse_01deg_coordinate_bins":
              sum(len(v)>=2 for v in self.spatial_bins.values()),
          "n_species_with_any_event_dated_relation":
              len(self.dated_species) if self.available_event_time_fields else None,
          "source_geographical_coordinate_fields_available":
              ("decimalLatitude" in self.available_geography_fields and
               "decimalLongitude" in self.available_geography_fields),
          "source_event_time_field_available":
              bool(self.available_event_time_fields),
          "GloBI_event_time_unassessable_when_field_missing":True,
          "n_original_GloBI_rows_read":int(self.raw["n_original_GloBI_rows_read"]),
          "n_exact_species_predicate_rows":int(self.raw["n_exact_species_predicate_rows"]),
          "n_eligible_species_flower_relation_rows":
              int(self.raw["n_eligible_species_flower_relation_rows"]),
          "n_georeferenced_species_relation_rows":
              int(self.raw["n_georeferenced_species_relations"]),
          "n_event_dated_species_relation_rows":
              int(self.raw["n_species_relation_rows_with_event_date"]),
          "n_relations_missing_provider":
              int(self.raw["n_relation_rows_without_provider"]),
          "n_relations_missing_bibliography":
              int(self.raw["n_relation_rows_without_bibliography"]),
          "n_crossprovider_shared_reference_keys_excluded":
              len(duplicate_refs),
          "n_exact_species_guild_provider_bibliography_units":
              len(filtered),
          "per_visitor_guild":by_guild_summary,
          "n_plant_populations_with_verified_island_mainland_flower_visits":0,
          "n_effective_pollen_delivery_observations_in_input":0,
          "n_confirmed_evolutionary_pollinator_shifts":0,
          "coordinates_are_not_island_polygon_membership":True,
          "original_cited_guild_is_not_a_seedset_or_pollination_success_observation":True,
          "no_exact_species_or_island_or_provider_id_released":True,
        }
