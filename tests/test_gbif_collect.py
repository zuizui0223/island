import io
import zipfile

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from island_v2.gbif_collect import (
    assign_occurrences_to_islands,
    island_species_table,
    read_block_occurrences,
    summarize_observation_effort,
)


def _islands():
    # Two small, disjoint unit squares far apart.
    return gpd.GeoDataFrame(
        {"island_id": ["isl_a", "isl_b"]},
        geometry=[
            Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]),
            Polygon([(10, 10), (10, 11), (11, 11), (11, 10), (10, 10)]),
        ],
        crs=4326,
    )


def _occurrences():
    # inside A, inside B, and one in neither (buffer-only).
    return pd.DataFrame(
        {
            "gbif_id": ["1", "2", "3", "4"],
            "dataset_key": ["d1", "d1", "d2", "d3"],
            "species": ["Aaa aaa", "Aaa aaa", "Bbb bbb", "Ccc ccc"],
            "scientific_name": ["Aaa aaa L.", "Aaa aaa L.", "Bbb bbb L.", "Ccc ccc L."],
            "decimal_longitude": [0.5, 0.6, 10.5, 5.0],
            "decimal_latitude": [0.5, 0.4, 10.5, 5.0],
            "year": ["2001", "2010", "1999", "2020"],
            "basis_of_record": [
                "PRESERVED_SPECIMEN",
                "HUMAN_OBSERVATION",
                "PRESERVED_SPECIMEN",
                "HUMAN_OBSERVATION",
            ],
        }
    )


def test_assignment_maps_points_to_exact_islands_and_drops_buffer_only():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())

    by_id = dict(zip(assigned["gbif_id"], assigned["island_id"], strict=True))
    assert by_id["1"] == "isl_a"
    assert by_id["2"] == "isl_a"
    assert by_id["3"] == "isl_b"
    # point 4 is inside the buffered catchment but on no real island
    assert pd.isna(by_id["4"])
    # no record is ever duplicated by the spatial join
    assert len(assigned) == 4


def test_island_species_table_is_deduplicated_with_provenance():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())
    table = island_species_table(assigned)

    aaa = table[(table["island_id"] == "isl_a") & (table["species"] == "Aaa aaa")].iloc[0]
    assert aaa["n_records"] == 2
    assert set(aaa["basis_of_record_set"].split("|")) == {"HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"}
    # the buffer-only Ccc ccc record is not attributed to any island
    assert "Ccc ccc" not in set(table["species"])


def test_observation_effort_counts_specimens_and_years():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())
    effort = summarize_observation_effort(assigned).set_index("island_id")

    assert effort.loc["isl_a", "n_records"] == 2
    assert effort.loc["isl_a", "n_species"] == 1
    assert effort.loc["isl_a", "n_datasets"] == 1
    assert effort.loc["isl_a", "n_preserved_specimen"] == 1
    assert effort.loc["isl_a", "year_min"] == 2001
    assert effort.loc["isl_a", "year_max"] == 2010


def test_read_block_occurrences_parses_simple_csv_and_drops_bad_coords(tmp_path):
    header = (
        "gbifID\tdatasetKey\tspecies\tscientificName\ttaxonRank\ttaxonKey\tspeciesKey\t"
        "decimalLatitude\tdecimalLongitude\tcoordinateUncertaintyInMeters\tyear\t"
        "basisOfRecord\toccurrenceStatus\testablishmentMeans"
    )
    good = "1\td1\tAaa aaa\tAaa aaa L.\tSPECIES\t100\t100\t0.5\t0.5\t10\t2001\tPRESERVED_SPECIMEN\tPRESENT\t"
    no_coord = "2\td1\tBbb bbb\tBbb bbb L.\tSPECIES\t200\t200\t\t\t\t2002\tHUMAN_OBSERVATION\tPRESENT\t"
    tsv = "\n".join([header, good, no_coord]) + "\n"

    archive_path = tmp_path / "0000-key.zip"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zf:
        zf.writestr("0000-key.csv", tsv)
    archive_path.write_bytes(buffer.getvalue())

    frame = read_block_occurrences(archive_path)

    assert list(frame["gbif_id"]) == ["1"]  # row without coordinates dropped
    assert frame.iloc[0]["species"] == "Aaa aaa"
    assert frame.iloc[0]["decimal_longitude"] == 0.5
