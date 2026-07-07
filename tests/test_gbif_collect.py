import io
import zipfile

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from island_v2.gbif_collect import (
    assign_occurrences_to_islands,
    deduplicate_occurrences,
    island_species_table,
    island_taxa_table,
    iter_block_occurrences,
    read_block_occurrences,
    summarize_observation_effort,
)

# This is the stable interface produced by island_taxa_table. It must not depend
# on a retired extraction implementation merely because a downstream workflow
# once consumed the same columns.
REQUIRED_TAXON_COLUMNS = {"accepted_species", "genus", "family", "n_islands", "n_records"}


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
            "family": ["Fabaceae", "Fabaceae", "Rosaceae", "Poaceae"],
            "genus": ["Aaa", "Aaa", "Bbb", "Ccc"],
            "species": ["Aaa aaa", "Aaa aaa", "Bbb bbb", "Ccc ccc"],
            "scientific_name": ["Aaa aaa L.", "Aaa aaa L.", "Bbb bbb L.", "Ccc ccc L."],
            "decimal_longitude": [0.5, 0.6, 10.5, 5.0],
            "decimal_latitude": [0.5, 0.4, 10.5, 5.0],
            "coordinate_uncertainty_m": [10, 2000, 50, 100],
            "year": ["2001", "2010", "1999", "2020"],
            "basis_of_record": [
                "PRESERVED_SPECIMEN",
                "HUMAN_OBSERVATION",
                "PRESERVED_SPECIMEN",
                "HUMAN_OBSERVATION",
            ],
            "establishment_means": ["NATIVE", "INTRODUCED", "", "CULTIVATED"],
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


def test_assignment_retains_exact_boundary_points():
    boundary = pd.DataFrame(
        {
            "gbif_id": ["boundary"],
            "decimal_longitude": [0.0],
            "decimal_latitude": [0.5],
        }
    )
    assigned = assign_occurrences_to_islands(boundary, _islands())
    assert assigned.loc[0, "island_id"] == "isl_a"


def test_deduplicate_occurrences_prefers_exact_island_copy_and_audits_overlap():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())
    # Same GBIF record appears in two overlapping regional block archives. The
    # exact-island copy must win over a buffer-only duplicate.
    exact = assigned.iloc[[0]].copy()
    exact["block_id"] = "block_a"
    buffer_only = assigned.iloc[[0]].copy()
    buffer_only["island_id"] = pd.NA
    buffer_only["coordinate_uncertainty_m"] = 1
    buffer_only["block_id"] = "block_b"
    unique = assigned.iloc[[2]].copy()
    unique["block_id"] = "block_b"

    deduplicated, audit = deduplicate_occurrences(pd.concat([buffer_only, exact, unique], ignore_index=True))

    assert len(deduplicated) == 2
    kept = deduplicated.set_index("gbif_id")
    assert kept.loc["1", "island_id"] == "isl_a"
    assert bool(kept.loc["1", "is_cross_block_duplicate"])
    assert audit["n_duplicate_gbif_ids"] == 1
    assert audit["n_duplicate_rows_removed"] == 1
    assert audit["n_cross_block_duplicate_gbif_ids"] == 1


def test_island_species_table_is_deduplicated_with_provenance():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())
    table = island_species_table(assigned)

    aaa = table[(table["island_id"] == "isl_a") & (table["species"] == "Aaa aaa")].iloc[0]
    assert aaa["n_records"] == 2
    assert aaa["n_unique_gbif_ids"] == 2
    assert set(aaa["basis_of_record_set"].split("|")) == {"HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"}
    # the buffer-only Ccc ccc record is not attributed to any island
    assert "Ccc ccc" not in set(table["species"])


def test_island_taxa_table_has_stable_downstream_contract():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())
    taxa = island_taxa_table(assigned)

    assert REQUIRED_TAXON_COLUMNS.issubset(set(taxa.columns))

    by_species = taxa.set_index("accepted_species")
    assert by_species.loc["Aaa aaa", "genus"] == "Aaa"
    assert by_species.loc["Aaa aaa", "family"] == "Fabaceae"
    assert by_species.loc["Aaa aaa", "n_records"] == 2
    # buffer-only species is excluded from the taxa list too
    assert "Ccc ccc" not in by_species.index


def test_observation_effort_counts_basis_quality_and_establishment():
    assigned = assign_occurrences_to_islands(_occurrences(), _islands())
    effort = summarize_observation_effort(assigned).set_index("island_id")

    assert effort.loc["isl_a", "n_records"] == 2
    assert effort.loc["isl_a", "n_unique_gbif_ids"] == 2
    assert effort.loc["isl_a", "n_species"] == 1
    assert effort.loc["isl_a", "n_datasets"] == 1
    assert effort.loc["isl_a", "n_preserved_specimen"] == 1
    assert effort.loc["isl_a", "n_human_observation"] == 1
    assert effort.loc["isl_a", "n_coordinate_uncertainty_gt_1000m"] == 1
    assert effort.loc["isl_a", "n_establishment_native"] == 1
    assert effort.loc["isl_a", "n_establishment_introduced"] == 1
    assert effort.loc["isl_a", "year_min"] == 2001
    assert effort.loc["isl_a", "year_max"] == 2010


def test_read_block_occurrences_parses_simple_csv_and_drops_bad_coords(tmp_path):
    header = (
        "gbifID\tdatasetKey\tfamily\tgenus\tspecies\tscientificName\ttaxonRank\ttaxonKey\tspeciesKey\t"
        "decimalLatitude\tdecimalLongitude\tcoordinateUncertaintyInMeters\tyear\t"
        "basisOfRecord\toccurrenceStatus\testablishmentMeans"
    )
    good = "1\td1\tFabaceae\tAaa\tAaa aaa\tAaa aaa L.\tSPECIES\t100\t100\t0.5\t0.5\t10\t2001\tPRESERVED_SPECIMEN\tPRESENT\t"
    no_coord = "2\td1\tRosaceae\tBbb\tBbb bbb\tBbb bbb L.\tSPECIES\t200\t200\t\t\t\t2002\tHUMAN_OBSERVATION\tPRESENT\t"
    tsv = "\n".join([header, good, no_coord]) + "\n"

    archive_path = tmp_path / "0000-key.zip"
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("0000-key.csv", tsv)
    archive_path.write_bytes(buffer.getvalue())

    frame = read_block_occurrences(archive_path)
    chunks = list(iter_block_occurrences(archive_path, chunksize=1))

    assert list(frame["gbif_id"]) == ["1"]  # row without coordinates dropped
    assert frame.iloc[0]["species"] == "Aaa aaa"
    assert frame.iloc[0]["genus"] == "Aaa"
    assert frame.iloc[0]["family"] == "Fabaceae"
    assert frame.iloc[0]["decimal_longitude"] == 0.5
    assert sum(chunk.attrs["n_source_rows"] for chunk in chunks) == 2
