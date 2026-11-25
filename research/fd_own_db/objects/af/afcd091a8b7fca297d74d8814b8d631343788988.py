"""Public original study/site metadata completeness, never GIFT source data."""
import gzip
from hashlib import sha256
import pandas as pd
import pytest
from fd.euppollnet_public_preflight import (
    SOURCE_FILES,available_columns,exact_git_blob_sha,
    qualify_original_blob,summarize_public_euppollnet,
)


def test_metadata_requires_original_site_and_study_ids():
    assert available_columns([
      "Study_id","Network_id","Plant_accepted_name",
      "Pollinator_accepted_name","Latitude","Longitude","Date",
      "Pollinator_family"])==dict(
      study_id="Study_id",network_id="Network_id",
      plant_name="Plant_accepted_name",
      pollinator_name="Pollinator_accepted_name",
      latitude="Latitude",longitude="Longitude",country=None,
      site=None,date="Date",flower_count_key=None,
      animal_order=None,animal_family="Pollinator_family")
    assert available_columns(["continent"])["latitude"] is None


def test_blob_hash_must_use_git_object_header_and_original_signed_identity():
    blob=b"example"
    expect=exact_git_blob_sha(blob)
    assert len(expect)==40 and expect!=sha256(blob).hexdigest()[:40]
    with pytest.raises(ValueError):
        qualify_original_blob("Unknown file",blob)
    with pytest.raises(ValueError,match="identity changed"):
        qualify_original_blob("Metadata.csv",blob)


def test_original_public_dataset_counts_and_geography_missingness(monkeypatch):
    source=pd.DataFrame([
        dict(Study_id="A",Network_id="N1",Plant_accepted_name="Alpha alba",
             Pollinator_accepted_name="Bombus terrestris",
             Latitude="44.0",Longitude="8.0",Country="Italy",Date="2020-06-01"),
        dict(Study_id="A",Network_id="N1",Plant_accepted_name="Beta alba",
             Pollinator_accepted_name="Syrphus ribesii",
             Latitude="44.0",Longitude="8.0",Country="Italy",Date="2020-06-02"),
        dict(Study_id="B",Network_id="N1",Plant_accepted_name="Alpha alba",
             Pollinator_accepted_name="Bombus terrestris",
             Latitude="unknown",Longitude="8.0",Country="Italy",Date="2019-05-03"),
    ])
    interactions=gzip.compress(source.to_csv(index=False).encode())
    metadata=b"Study_id,Year\\nA,2020\\nB,2019\\n".replace(b"\\n",b"\n")
    monkeypatch.setitem(SOURCE_FILES,"Interaction_data.csv.gz",
                        (len(interactions),exact_git_blob_sha(interactions)))
    monkeypatch.setitem(SOURCE_FILES,"Metadata.csv",
                        (len(metadata),exact_git_blob_sha(metadata)))
    out=summarize_public_euppollnet(interactions,metadata)
    assert out["n_original_source_interaction_rows"]==3
    assert out["n_distinct_original_study_ids"]==2
    assert out["n_distinct_original_study_network_combinations"]==2
    assert out["n_original_source_distinct_named_plants"]==2
    assert out["n_original_source_distinct_named_pollinators"]==2
    assert out["n_original_georeferenced_rows"]==2
    assert out["n_distinct_geographic_coordinate_pairs"]==1
    assert out["n_ungeoreferenced_source_rows"]==1
    assert out["original_metadata_rows"]==2
    assert out["limits"]["source_original_study_overlap_with_GloBI_DoPI_WOL_unknown"]


def test_invalid_utf8_plant_name_is_quarantined_not_used_as_taxonomic_match(monkeypatch):
    raw=(b"Study_id,Network_id,Plant_accepted_name,Pollinator_accepted_name,Latitude,Longitude\n"
         b"A,N1,Alph\x8c alba,Bombus terrestris,44.0,8.0\n"
         b"A,N1,Beta alba,Syrphus ribesii,44.0,8.0\n")
    interactions=gzip.compress(raw)
    metadata=b"Study_id,Year\nA,2020\n"
    monkeypatch.setitem(SOURCE_FILES,"Interaction_data.csv.gz",
                        (len(interactions),exact_git_blob_sha(interactions)))
    monkeypatch.setitem(SOURCE_FILES,"Metadata.csv",
                        (len(metadata),exact_git_blob_sha(metadata)))
    q=summarize_public_euppollnet(interactions,metadata)
    assert q["n_original_source_interaction_rows"]==2
    assert q["n_original_plant_name_rows_quarantined_invalid_UTF8"]==1
    assert q["n_original_source_distinct_named_plants"]==1
    assert q["raw_csv_UTF8_invalid_bytes_not_silently_cast_to_taxonomy"]
