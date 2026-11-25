"""Pinned-public EuPPollNet flower-availability denominator audit.

This is SOURCE QUALITY/FEASIBILITY, never inferred visitor preference.
The original Flower_counts and Interaction_data join on an explicit
published Flower_data_merger key. Repeated key matches are NOT treated
as independent flower observations or automatically summed/merged.
No original GIFT island data/accepted names are loaded here.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import numpy as np
import pandas as pd
from .euppollnet_public_preflight import (
    SOURCE_COMMIT,SOURCE_FILES,qualify_original_blob,exact_git_blob_sha,
)

FLOWER_SOURCE_FILE="Flower_counts.csv"
FLOWER_SOURCE_SIZE=20315455
FLOWER_GIT_BLOB="8543e446a53f3b7cffece76e1699489e7169f7fd"
JOIN="Flower_data_merger"


def source_original_flower_receipt(data:bytes)->dict:
    if len(data)!=FLOWER_SOURCE_SIZE or exact_git_blob_sha(data)!=FLOWER_GIT_BLOB:
        raise ValueError("Original public EuPPollNet Flower_counts source changed")
    return dict(original_source_file=FLOWER_SOURCE_FILE,
        original_git_commit=SOURCE_COMMIT,
        original_file_bytes=len(data),
        original_git_blob_sha1=FLOWER_GIT_BLOB,
        original_file_sha256=sha256(data).hexdigest())


def _nonmissing(q:pd.Series)->pd.Series:
    v=q.dropna().astype(str).str.strip()
    return v[v.ne("")&~v.str.lower().isin(
        ("na","nan","null","none","unknown"))]


def audit_public_flower_source(flower_blob:bytes, interaction_blob:bytes,
                              *, allow_synthetic_source:bool=False)->dict:
    """Source-original joint-key presence, no abundance-model invention."""
    if not allow_synthetic_source:
        receipt=source_original_flower_receipt(flower_blob)
        visitor=qualify_original_blob("Interaction_data.csv.gz",interaction_blob)
    else:
        receipt=dict(original_source_file="SYNTHETIC",
                     original_file_bytes=len(flower_blob))
        visitor=dict(original_public_source="SYNTHETIC")
    try:
        f=pd.read_csv(BytesIO(flower_blob),dtype=str,low_memory=False,
            encoding="utf-8",encoding_errors="replace")
        i=pd.read_csv(BytesIO(interaction_blob),compression="gzip",
            dtype=str,low_memory=False,encoding="utf-8",
            encoding_errors="replace",
            usecols=lambda c:c==JOIN)
    except (ValueError,UnicodeError,EOFError) as exc:
        raise ValueError("Original public source CSV impossible to parse") from exc
    if JOIN not in f or JOIN not in i:
        raise ValueError("Published Flower_data_merger key missing")
    fkey=_nonmissing(f[JOIN]);ikey=_nonmissing(i[JOIN])
    set_f,set_i=set(fkey),set(ikey)
    common=set_f&set_i
    freq_f=fkey.value_counts()
    freq_i=ikey.value_counts()
    # Any source join key with multiple flower-count rows requires survey
    # grain validation BEFORE use as botanical exposure/availability.
    ambiguous=set(freq_f.index[freq_f.gt(1)])
    potential={x for x in common if x not in ambiguous}
    # Headers only; never publish original per-plant/site/taxon join key values.
    numeric=[c for c in f.columns if c!=JOIN and any(x in c.lower()
             for x in ("flower_count","flowers_count","flower_n",
                       "n_flowers","flower_abundance","flower_number"))]
    availability={}
    for col in numeric:
        q=pd.to_numeric(f[col],errors="coerce")
        good=np.isfinite(q)
        availability[col]=dict(n_nonnegative_numeric=int((q[good]>=0).sum()),
            n_negative_numeric=int((q[good]<0).sum()),
            n_missing_or_nonnumeric=int((~good).sum()))
    return dict(
        status="PUBLIC_EUPPOLLNET_ORIGINAL_FLOWER_EXPOSURE_KEY_FEASIBILITY",
        source_original_flower_receipt=receipt,
        source_original_visitor_receipt=visitor,
        n_original_flowers_raw_rows=int(len(f)),
        n_original_interaction_raw_rows=int(len(i)),
        n_nonempty_flower_join_rows=int(len(fkey)),
        n_nonempty_interaction_join_rows=int(len(ikey)),
        n_distinct_flower_join_keys=int(len(set_f)),
        n_distinct_interaction_join_keys=int(len(set_i)),
        n_original_distinct_join_keys_shared=int(len(common)),
        n_distinct_flower_join_keys_with_multiple_source_rows=int(len(ambiguous)),
        n_original_join_keys_shared_and_unique_flower_count_rows=int(len(potential)),
        n_original_interaction_rows_on_any_shared_join_key=int(
            ikey.isin(common).sum()),
        n_original_interaction_rows_with_unique_flower_key=int(
            ikey.isin(potential).sum()),
        flower_source_column_names=list(f.columns),
        original_quantitative_flower_abundance_column_candidates=availability,
        conclusions=dict(
            source_species_taxonomy_and_original_GIFT_flora_not_yet_joined=True,
            original_plant_flower_abundance_not_assumed_available_for_every_visit=True,
            flower_source_duplicate_join_keys_not_auto_aggregated=True,
            original_source_codebook_needed_to_interpret_count_units=True,
            valid_source_flower_count_is_not_measured_pollen_deposition=True,
            original_European_study_not_tropical_island_replication=True,
            no_effect_coefficient_or_pvalue_estimated=True,
        ))


def qualified_gift_flower_supply_on_positive_visits(
    source_local_visits:pd.DataFrame, original_flower_counts:pd.DataFrame,
    direct_traits:dict, *, original_GIFT_botanical_floras:int=4505,
)->dict:
    """PRIVATE original exact-GIFT-site source dataframe -> anonymous support.

    The input source_local_visits MUST already be joined to independent GIFT
    original exact GSHHG polygon and SAME island accepted botanical species.
    Counts are not visitor choice or availability estimates. Require ONE
    published flower-count row per source key and identical study/network/
    original date to avoid false exposure links. Original source names,
    sites, study IDs and flower keys are NEVER serialized.
    """
    from .fullframe_globi_local_support import genuine_direct_motifs
    from .euppollnet_gift_exact_site_bridge import original_guild
    if original_GIFT_botanical_floras!=4505:
        raise ValueError("Original full 4505 botanical island source population lost")
    required_visits={
        "Study_id","Network_id","Flower_data_merger","Date",
        "island_id","plant","origin_status","endemic_status",
        "Pollinator_order","Pollinator_family","Pollinator_accepted_name",
    }
    required_flower={
        "Flower_data_merger","Study_id","Network_id","Day","Month","Year",
        "Flower_count","Units",
    }
    if not required_visits<=set(source_local_visits) or not required_flower<=set(original_flower_counts):
        raise ValueError("Original field source/GIFT positive visitors or original supply missing")
    visit=source_local_visits.copy()
    fl=original_flower_counts.copy()
    if len(visit)==0:
        return dict(status="ORIGINAL_GIFT_SITE_EUPPOLLNET_NO_POSITIVE_VISITS",
            n_original_qualified_visitors=0,
            original_source_unsampled_not_animal_absence=True)
    if visit[["Study_id","Network_id","island_id","plant"]].isna().any().any():
        raise ValueError("Source field visitor ID or botanical flora site missing")
    # The original source key is a published survey linker, NOT a new species
    # observation. Empty/null source keys are unqualified.
    visit[JOIN]=visit[JOIN].astype("string").str.strip()
    fl[JOIN]=fl[JOIN].astype("string").str.strip()
    visit=visit.loc[visit[JOIN].notna()&visit[JOIN].ne("")].copy()
    fl=fl.loc[fl[JOIN].notna()&fl[JOIN].ne("")].copy()
    frequencies=fl[JOIN].value_counts()
    unique=set(frequencies.index[frequencies.eq(1)])
    n_source_original_flower_keys_ambiguous=int(frequencies.gt(1).sum())
    exact=fl.loc[fl[JOIN].isin(unique)]
    # Never auto-sum multiple positive source visits into extra independent
    # flower-count denominators.
    linked=visit.merge(exact,on=JOIN,how="inner",
        suffixes=("","_flower"),validate="many_to_one")
    n_key_linked=len(linked)
    original_study_network=(linked.Study_id.eq(linked.Study_id_flower)&
            linked.Network_id.eq(linked.Network_id_flower))
    linked=linked.loc[original_study_network].copy()
    n_same_study_network=len(linked)
    visit_date=pd.to_datetime(linked.Date,errors="coerce")
    year=pd.to_numeric(linked.Year,errors="coerce")
    month=pd.to_numeric(linked.Month,errors="coerce")
    day=pd.to_numeric(linked.Day,errors="coerce")
    same=(visit_date.notna()&year.eq(visit_date.dt.year)&
          month.eq(visit_date.dt.month)&day.eq(visit_date.dt.day))
    linked=linked.loc[same].copy()
    n_same_original_calendar_day=len(linked)
    numeric=pd.to_numeric(linked.Flower_count,errors="coerce")
    positive=numeric.notna() & numeric.gt(0)
    units=linked.Units.astype(str).str.strip()
    valid_unit=units.ne("")&~units.str.lower().isin(
        ("nan","none","unknown","na","null"))
    linked=linked.loc[positive&valid_unit].copy()
    n_quantitative_exposure_rows=len(linked)
    # Exact true-same-species floral availability labels only: no genus
    # imputation and no trait-state imputation from different accepted plants.
    motifs=genuine_direct_motifs(direct_traits)
    verified=linked.copy()
    verified["visit_guild"]=[original_guild(o,f,a) for o,f,a in zip(
        verified.Pollinator_order,verified.Pollinator_family,
        verified.Pollinator_accepted_name)]
    cells=[]
    for scope in ("all_recorded","native_label","native_endemic"):
        if scope=="all_recorded":subset=verified
        elif scope=="native_label":subset=verified[
            verified.origin_status.eq("native")]
        else:subset=verified[
            verified.origin_status.eq("native")&
            verified.endemic_status.eq("endemic")]
        for motif,by_species in motifs.items():
            scored={s for s,v in by_species.items() if v}
            yes=subset.loc[subset.plant.isin(scored)]
            for guild in ("bombus","other_bees","butterfly_family",
                          "sphingidae_moth","flies_order","other_lepidoptera"):
                q=yes.loc[yes.visit_guild.eq(guild)]
                cells.append(dict(botanical_scope=scope,
                    genuine_same_species_flower_pair=motif,
                    original_animal_taxonomic_guild=guild,
                    n_original_source_gift_island_qualified_flower_supply_sites=int(
                        q.island_id.nunique()),
                    n_original_source_studies_with_qualified_flower_supply=int(
                        q.Study_id.nunique()),
                    n_original_accepted_species_with_qualified_flower_supply=int(
                        q.plant.nunique()),
                    n_original_qualified_flower_count_linked_visitor_rows=int(len(q))))
    source_units=(linked.Units.fillna("MISSING").astype(str).str.strip().
                  value_counts().to_dict())
    return dict(
        status="ORIGINAL_GIFT_EUPPOLLNET_EXACT_POSITIVE_FLOWER_SUPPLY_LINK_AUDIT",
        parent_original_GIFT_botanical_floras=4505,
        n_original_source_positive_same_island_botanical_visits=len(source_local_visits),
        n_original_positive_source_with_nonmissing_join_key=len(visit),
        n_original_source_flower_keys_ambiguous_not_automatically_summed=n_source_original_flower_keys_ambiguous,
        n_original_visitor_rows_with_unique_original_flower_key=n_key_linked,
        n_original_visitors_with_same_published_study_network_flower_source=n_same_study_network,
        n_original_visitors_with_same_published_study_network_calendar_day=n_same_original_calendar_day,
        n_original_positive_visitors_with_nonzero_numeric_flower_count_and_unit=n_quantitative_exposure_rows,
        n_original_GIFT_islands_with_qualified_flower_exposure=int(linked.island_id.nunique()),
        n_original_source_studies_with_qualified_flower_exposure=int(linked.Study_id.nunique()),
        n_original_accepted_botanical_plants_with_qualified_flower_exposure=int(linked.plant.nunique()),
        original_qualified_source_flower_unit_labels_counts=source_units,
        true_source_same_species_flower_x_animal_positive_exposure_cells=cells,
        limitations=dict(
            original_source_StudyNetworkDate_plant_supply_matches_are_not_independent_guild_preferences=True,
            flower_unit_heterogeneity_has_not_been_normalized=True,
            repeated_visits_do_not_add_independent_flower_denominators=True,
            original_animal_nondetections_and_plant_effort_not_verified=True,
            source_flower_traits_not_measured_on_same_european_field_site=True,
            no_calibrated_pollen_deposition_or_seed_set_estimated=True,
            original_sources_overlap_GloBI_DoPI_WOL_unresolved=True,
            no_private_original_study_network_island_plant_or_flower_source_keys_exported=True))
