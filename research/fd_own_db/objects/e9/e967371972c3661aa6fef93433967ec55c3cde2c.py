"""Conservative survey CHOICE-SET FEASIBILITY, not estimated visitor preference.

A positive recorded visitor and same-species flower count do not establish
that *all* flowers available at that site/date were measured. Check whether
the original EuPPollNet Flower_counts inventory includes other accepted
GIFT flora plants in the SAME source study, network, Site_id, day and
published physical Units, anchored by at least one EXACT positive source
visit-to-same-source-flower-plant identity.

Require the original study×network to have a single WGS84 GPS location
ACROSS ALL original public interaction records; otherwise source Site_id
cannot safely be mapped to one GIFT island. No island, plant, network or
source study identifiers leave source CI. A flower with NO visitor entry
is only an undocumented relation, never a true negative observation.
"""
from __future__ import annotations

import pandas as pd

from .euppollnet_gift_exact_site_bridge import strict_original_accepted_plant
from .fullframe_globi_local_support import genuine_direct_motifs

MOTIFS=("actinomorphic_deep","zygomorphic_deep","tube_bell_deep",
        "white_tube_bell","red_pink_tube_bell")
ANCHOR=("Study_id","Network_id","Site_id","Date","Units","island_id")
LINK=("Study_id","Network_id","Site_id","Date","Units")
KEY="Flower_data_merger"


def _clean_string(q:pd.Series)->pd.Series:
    text=q.astype("string").str.strip()
    return text.mask(text.str.lower().isin(
        ("","none","nan","na","null","unknown")).fillna(True))


def _qualify_all_original_study_network_locations(positive:pd.DataFrame):
    required={"Study_id","Network_id","Latitude","Longitude"}
    if not required<=set(positive):
        raise ValueError("Original full source coordinates unavailable")
    raw=positive[list(required)].copy()
    raw["Latitude"]=pd.to_numeric(raw.Latitude,errors="coerce")
    raw["Longitude"]=pd.to_numeric(raw.Longitude,errors="coerce")
    raw=raw.loc[raw.Latitude.between(-90,90)&raw.Longitude.between(-180,180)]
    if raw.empty:
        return set(),0
    points=raw[["Study_id","Network_id","Latitude","Longitude"]].drop_duplicates()
    counts=points.groupby(["Study_id","Network_id"],dropna=False).size()
    qualified=set(counts.index[counts.eq(1)])
    return qualified,int(counts.gt(1).sum())


def audited_nonvisited_flower_census_opportunity(
    exact_same_species_positive_source:pd.DataFrame,
    original_all_flower_counts:pd.DataFrame,
    original_full_source_positive_visits:pd.DataFrame,
    original_GIFT_full_flora:pd.DataFrame,
    direct_species_traits:dict,*,parent_botanical_floras:int=4505,
)->dict:
    """Return only anonymous candidate choice-set coverage, NOT animal choice.

    exact_same_species_positive_source is the private protected output of
    qualified_rows_and_identity_receipt() and already has validated source
    flower species Plant_species, positive count and study/network/day.
    """
    if parent_botanical_floras!=4505 or original_GIFT_full_flora.island_id.nunique()!=4505:
        raise ValueError("Original GIFT 4505 botanical population dropped")
    req={
        "Study_id","Network_id","Site_id","Date","Units","island_id",
        "plant","Flower_data_merger"
    }
    flower_req={
        "Study_id","Network_id","Site_id","Day","Month","Year",
        "Units","Plant_species","Flower_count","Flower_data_merger"
    }
    if (not req<=set(exact_same_species_positive_source) or
            not flower_req<=set(original_all_flower_counts) or
            not {"island_id","accepted_species"}<=set(original_GIFT_full_flora)):
        raise ValueError("Original exact source site/date GIFT visitor/plant census schema missing")
    if original_GIFT_full_flora.duplicated(
            ["island_id","accepted_species"]).any():
        raise ValueError("Original accepted GIFT botanical list duplicates")

    allowed_networks,multi_geo_count=_qualify_all_original_study_network_locations(
        original_full_source_positive_visits)
    v=exact_same_species_positive_source.copy()
    n_original=len(v)
    v["Site_id"]=_clean_string(v.Site_id)
    v["Units"]=_clean_string(v.Units)
    v=v.loc[v.Site_id.notna()&v.Units.notna()].copy()
    n_source_site_ids=len(v)
    loc_ok=[(study,net) in allowed_networks for study,net in zip(
        v.Study_id,v.Network_id)]
    v=v.loc[loc_ok].copy()
    n_single_coordinate=len(v)
    if len(v)==0:
        return dict(
            status="GIFT_4505_EUPPOLLNET_COMPLETE_SITE_FLOWER_CHOICESET_UNIDENTIFIABLE",
            original_global_flora_parent=4505,
            source_waterfall=dict(
                n_source_positive_same_plant_linked_rows=n_original,
                n_positive_rows_with_source_site_and_unit=n_source_site_ids,
                n_positive_rows_on_single_GPS_source_study_network=n_single_coordinate,
                n_source_study_network_keys_with_multiple_original_GPS=multi_geo_count),
            absent_source_choice_set_is_NOT_no_animals=True,
            no_visitation_preference_or_pollen_delivery_estimated=True)

    anchors=v[list(ANCHOR)].drop_duplicates()
    # The same public flower survey location key cannot map to two original
    # different GIFT island IDs. Fail closed rather than selecting a donor.
    ambiguous=anchors.groupby(list(LINK),dropna=False).island_id.nunique()
    if ambiguous.gt(1).any():
        raise ValueError("Same EuPPollNet source site/day/unit matches multiple GIFT original islands")

    f=original_all_flower_counts.copy()
    f["Site_id"]=_clean_string(f.Site_id)
    f["Units"]=_clean_string(f.Units)
    f[KEY]=_clean_string(f[KEY])
    f=f.loc[f.Site_id.notna()&f.Units.notna()&f[KEY].notna()].copy()
    frequency=f[KEY].value_counts()
    clean_unique=set(frequency.index[frequency.eq(1)])
    f=f.loc[f[KEY].isin(clean_unique)].copy()
    for part in ("Day","Month","Year"):
        f[part]=pd.to_numeric(f[part],errors="coerce")
    f["Date"]=pd.to_datetime(
        dict(year=f.Year,month=f.Month,day=f.Day),
        errors="coerce").dt.strftime("%Y-%m-%d")
    f=f.loc[f.Date.notna()].copy()
    f["Flower_count"]=pd.to_numeric(f.Flower_count,errors="coerce")
    f=f.loc[f.Flower_count.notna()&f.Flower_count.gt(0)].copy()
    names=set(f.Plant_species.dropna().astype(str))
    normalized={x:strict_original_accepted_plant(x) for x in names}
    f["_accepted_flower_plant"]=f.Plant_species.map(normalized)
    f=f.loc[f._accepted_flower_plant.notna()].copy()

    candidates=anchors.merge(f,on=list(LINK),
        how="inner",validate="one_to_many")
    candidate_before_flora=len(candidates)
    candidates=candidates.merge(
        original_GIFT_full_flora[["island_id","accepted_species"]].rename(
            columns={"accepted_species":"_accepted_flower_plant"}),
        on=["island_id","_accepted_flower_plant"],
        how="inner",validate="many_to_one")
    candidate_after_flora=len(candidates)

    # Group site/day/unit at the original survey grain, NOT one
    # network/year pooled across sites, and count availability all plants
    # even if they were not in a positive documented visitor edge.
    group_keys=list(ANCHOR)
    visited_by={key:set(g.plant) for key,g in v.groupby(group_keys,dropna=False)}
    candidate_contexts={}
    for key,g in candidates.groupby(group_keys,dropna=False):
        site_plants=set(g._accepted_flower_plant)
        visitor_plants=visited_by.get(key,set())
        candidate_contexts[key]=dict(
            source_plants=site_plants,
            positive_documented_plants=visitor_plants,
        )
    original_motifs=genuine_direct_motifs(direct_species_traits)
    cells=[]
    for motif in MOTIFS:
        states=original_motifs[motif]
        opportunity=[]
        for key,row in candidate_contexts.items():
            scored=row["source_plants"]&set(states)
            positives={name for name in scored if states[name]}
            negatives=scored-positives
            # "Negatives" here are only NOT in botanical morphology class,
            # never source pollinator nondetections.
            if positives and negatives:
                opportunity.append(key)
        cells.append(dict(
            true_same_species_direct_motif=motif,
            n_original_site_day_unit_source_flora_contexts_with_both_scored_morphology_states=len(opportunity),
            n_original_source_studies_with_both_scored_flower_states=len({
                k[0] for k in opportunity}),
            n_original_GIFT_islands_with_both_scored_flower_states=len({
                k[5] for k in opportunity}),
            only_botanical_morphology_true_false_not_animal_nonvisit=True,
        ))

    source_no_recorded_positive=sum(
        len(row["source_plants"]-row["positive_documented_plants"])
        for row in candidate_contexts.values())
    return dict(
        status="GIFT_4505_EUPPOLLNET_FULL_SITE_DAY_FLOWER_CENSUS_CHOICESET_PREFLIGHT",
        original_global_flora_parent=4505,
        source_waterfall=dict(
            n_source_positive_same_plant_linked_rows=n_original,
            n_positive_rows_with_source_site_and_unit=n_source_site_ids,
            n_positive_rows_on_single_GPS_source_study_network=n_single_coordinate,
            n_source_study_network_keys_with_multiple_original_GPS=multi_geo_count,
            n_original_distinct_positive_anchor_site_day_unit_contexts=len(anchors),
            n_original_all_source_flower_rows_on_anchored_study_site_day_unit_before_GIFT_check=candidate_before_flora,
            n_source_census_flower_rows_with_original_same_island_accepted_species=candidate_after_flora,
            n_distinct_site_day_unit_contexts_with_source_flower_census=len(candidate_contexts),
            n_original_source_studies_with_flower_census_on_GIFT_island=len({
                k[0] for k in candidate_contexts}),
            n_original_GIFT_islands_with_flower_census=len({
                k[5] for k in candidate_contexts}),
            n_source_census_original_flower_species_contexts_without_positive_reference=source_no_recorded_positive,
            n_original_unique_flower_source_keys_in_census=int(candidates[KEY].nunique()) if len(candidates) else 0,
        ),
        true_same_species_botanical_morphology_source_census_opportunities=cells,
        limits=dict(
            original_european_flower_census_published_does_not_establish_complete_all_flower_choice_set=True,
            surveyed_flower_without_recorded_animal_link_not_biological_nondetection=True,
            site_id_is_from_flower_source_not_geocoordinate_independently_validated=True,
            exact_source_network_single_GPS_essential=True,
            no_poisson_offset_or_choice_effect_without_source_survey_effort_and_sampling_completeness=True,
            no_original_private_plant_island_network_study_key_or_coordinates_serialized=True,
        ),
    )
