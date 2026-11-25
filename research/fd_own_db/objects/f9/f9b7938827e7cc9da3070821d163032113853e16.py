"""Public EuPPollNet source-only site/time/animal taxonomy preflight.

Designed as an ALTERNATIVE ORIGINAL FIELD-SITE DATASET to the aggregated
GloBI positive source corpus; do not assume study independence from DoPI.
This preflight DOES NOT access GIFT island accepted taxon IDs or match
unverified locations; it only inventories publicly available metadata and
freezes the original Git blob byte signature.
"""
from __future__ import annotations
from hashlib import sha1,sha256
from io import BytesIO
import gzip
import pandas as pd

SOURCE_REPO="joseblanuza/EuPPollNet"
SOURCE_COMMIT="31716d21d5aad8c27c1d91dd0620149ef42e8836"
SOURCE_FILES={
  "Interaction_data.csv.gz":(2373239,"876db2c6fc2aec465be0b69f735922800f202bb2"),
  "Metadata.csv":(28483,"f9708e2f2b142ab2437a8be9d16e33aaefdef0ed")
}
CANONICAL=("study_id","network_id","plant_accepted_name",
           "pollinator_accepted_name")


def exact_git_blob_sha(data:bytes)->str:
    # GitHub content API SHA1 is a Git BLOB identity, not plain SHA1(bytes).
    return sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()


def qualify_original_blob(name:str,data:bytes)->dict:
    if name not in SOURCE_FILES or not isinstance(data,bytes):
        raise ValueError("Original frozen EuPPollNet public dataset file unknown")
    expect_n,expect_git=SOURCE_FILES[name]
    if len(data)!=expect_n or exact_git_blob_sha(data)!=expect_git:
        raise ValueError("Original EuPPollNet GitHub source blob identity changed")
    return dict(original_public_source=name,n_source_bytes=len(data),
        github_commit_sha=SOURCE_COMMIT,
        git_blob_sha1=expect_git,
        source_byte_sha256=sha256(data).hexdigest(),
        raw_european_file_not_independent_GIFT_animal_observation=True)


def available_columns(cols)->dict:
    """Use explicit field labels only; unknown taxonomic/site columns missing."""
    col={str(x).strip().lower():str(x).strip() for x in cols}
    def choose(*names):
        return next((col[x] for x in names if x in col),None)
    return dict(
      study_id=choose("study_id","studyid","study"),
      network_id=choose("network_id","networkid","network"),
      plant_name=choose("plant_accepted_name","plant_name","plant_taxon"),
      pollinator_name=choose("pollinator_accepted_name","pollinator_name",
                            "visitor_accepted_name"),
      latitude=choose("latitude","lat","site_latitude","latitude_dd","lat_wgs84"),
      longitude=choose("longitude","lon","long","site_longitude",
                       "longitude_dd","lon_wgs84"),
      country=choose("country","country_name"),
      site=choose("site_id","site","site_code","location","site_number"),
      date=choose("date","observation_date","sample_date"),
      flower_count_key=choose("flower_data_merger"),
      animal_order=choose("pollinator_order","pollinator_accepted_order",
                          "pollinator_taxon_order"),
      animal_family=choose("pollinator_family","pollinator_accepted_family",
                           "pollinator_taxon_family"),
    )


def summarize_public_euppollnet(interaction_data:bytes,
                               metadata:bytes)->dict:
    """Aggregate public source opportunity; never claim original visits on GIFT.

    Full distinct study, network/site/date counts are source metadata only.
    """
    receipt=[qualify_original_blob("Interaction_data.csv.gz",interaction_data),
             qualify_original_blob("Metadata.csv",metadata)]
    try:
        first=pd.read_csv(BytesIO(interaction_data),compression="gzip",
                          nrows=30,low_memory=False,encoding="utf-8",
                          encoding_errors="replace")
        meta=pd.read_csv(BytesIO(metadata),low_memory=False)
    except (ValueError,UnicodeError,EOFError) as exc:
        raise ValueError("Frozen EuPPollNet source CSV not readable") from exc
    columns=available_columns(first.columns)
    network_id=columns["network_id"]
    study_id=columns["study_id"]
    if not study_id or not network_id:
        raise ValueError("Original study and sampling-network IDs absent")
    names=set();networks=set(); sites=set();st=set();animals=set()
    species=set(); years=set();countries=set()
    nrows=0;nvalid_geo=0;nunknown_sites=0;valid_plant=0
    n_unusable_plant_encoding=0
    n_unusable_pollinator_encoding=0
    latc=columns["latitude"];lonc=columns["longitude"]
    for chunk in pd.read_csv(BytesIO(interaction_data),
            compression="gzip",dtype=str,chunksize=150000,
            low_memory=False,encoding="utf-8",encoding_errors="replace"):
        nrows+=len(chunk)
        if chunk[[study_id,network_id]].isna().any().any():
            raise ValueError("Public source has original study/network ID missing")
        if chunk[[study_id,network_id]].astype(str).apply(
                lambda c:c.str.contains("\ufffd",regex=False)).any().any():
            raise ValueError("Source study/network ID contains unrecoverable text encoding")
        st.update(chunk[study_id].unique())
        networks.update(zip(chunk[study_id],chunk[network_id]))
        if columns["plant_name"]:
            candidates=chunk[columns["plant_name"]].dropna()
            invalid=candidates.str.contains("\ufffd",regex=False)
            n_unusable_plant_encoding+=int(invalid.sum())
            good=candidates[~invalid]
            valid_plant+=len(good)
            species.update(good.unique())
        if columns["pollinator_name"]:
            visitor=chunk[columns["pollinator_name"]].dropna()
            invalid=visitor.str.contains("\ufffd",regex=False)
            n_unusable_pollinator_encoding+=int(invalid.sum())
            animals.update(visitor[~invalid].unique())
        if columns["country"]:
            countries.update(chunk[columns["country"]].dropna().unique())
        if columns["date"]:
            years.update(chunk[columns["date"]].dropna().str[:4].unique())
        if latc and lonc:
            lat=pd.to_numeric(chunk[latc],errors="coerce")
            lon=pd.to_numeric(chunk[lonc],errors="coerce")
            ok=lat.between(-90,90)&lon.between(-180,180)&lat.notna()&lon.notna()
            nvalid_geo+=int(ok.sum())
            sites.update(zip(lat[ok].round(6).astype(str),
                             lon[ok].round(6).astype(str)))
            nunknown_sites+=int((~ok).sum())
    if not nrows:
        raise ValueError("No actual European source flower visitor records")
    return dict(
        status="PUBLIC_EUPPOLLNET_ORIGINAL_SITE_EVIDENCE_PREFLIGHT",
        original_public_source_receipts=receipt,
        n_original_source_interaction_rows=nrows,
        n_distinct_original_study_ids=len(st),
        n_distinct_original_study_network_combinations=len(networks),
        n_original_source_distinct_named_plants=len(species),
        n_original_source_distinct_named_pollinators=len(animals),
        n_original_source_distinct_countries=len(countries),
        n_original_source_year_labels=len(years),
        n_original_georeferenced_rows=nvalid_geo if latc else None,
        n_distinct_geographic_coordinate_pairs=len(sites) if latc else None,
        n_ungeoreferenced_source_rows=nunknown_sites if latc else None,
        n_original_rows_with_plant_name=valid_plant,
        n_original_plant_name_rows_quarantined_invalid_UTF8=n_unusable_plant_encoding,
        n_original_pollinator_name_rows_quarantined_invalid_UTF8=n_unusable_pollinator_encoding,
        raw_csv_UTF8_invalid_bytes_not_silently_cast_to_taxonomy=True,
        source_interaction_headers=list(map(str,first.columns)),
        source_metadata_headers=list(map(str,meta.columns)),
        original_metadata_rows=len(meta),
        original_field_contract=columns,
        limits=dict(
          original_european_study_not_global_tropical_data=True,
          source_original_study_overlap_with_GloBI_DoPI_WOL_unknown=True,
          no_original_GIFT_source_island_taxon_or_polygon_join_done=True,
          original_site_coordinates_need_exact_GIFT_GSHHG_polygon_intersection=True,
          published_source_sampling_absence_not_plant_pollination_inefficacy=True,
          no_adaptive_evolution_or_pollen_deposition_test=True))
