import pandas as pd

from island_v2.gbif_collect_bounded import _summarize_exact_chunks


def _chunk(rows):
    return pd.DataFrame(rows)


def test_sqlite_spill_summarizes_chunks_and_removes_duplicate_gbif_ids():
    first = _chunk(
        [
            {
                "gbif_id": "1",
                "dataset_key": "dataset_a",
                "family": "Fabaceae",
                "genus": "Aaa",
                "species": "Aaa aaa",
                "coordinate_uncertainty_m": 1000,
                "year": "2001",
                "basis_of_record": "HUMAN_OBSERVATION",
                "establishment_means": "NATIVE",
                "island_id": "island_a",
            },
            {
                "gbif_id": "2",
                "dataset_key": "dataset_a",
                "family": "Fabaceae",
                "genus": "Aaa",
                "species": "Aaa aaa",
                "coordinate_uncertainty_m": 20,
                "year": "2002",
                "basis_of_record": "HUMAN_OBSERVATION",
                "establishment_means": "INTRODUCED",
                "island_id": "island_a",
            },
        ]
    )
    later_lower_uncertainty_copy = _chunk(
        [
            {
                "gbif_id": "1",
                "dataset_key": "dataset_a",
                "family": "Fabaceae",
                "genus": "Aaa",
                "species": "Aaa aaa",
                "coordinate_uncertainty_m": 10,
                "year": "2001",
                "basis_of_record": "PRESERVED_SPECIMEN",
                "establishment_means": "NATIVE",
                "island_id": "island_a",
            }
        ]
    )

    species, effort, taxa, audit = _summarize_exact_chunks(
        [first, later_lower_uncertainty_copy], "block_a"
    )

    row = species.iloc[0]
    assert row["island_id"] == "island_a"
    assert row["species"] == "Aaa aaa"
    assert row["n_records"] == 2
    assert row["n_unique_gbif_ids"] == 2
    assert set(row["basis_of_record_set"].split("|")) == {
        "HUMAN_OBSERVATION",
        "PRESERVED_SPECIMEN",
    }

    effort_row = effort.set_index("island_id").loc["island_a"]
    assert effort_row["n_records"] == 2
    assert effort_row["n_unique_gbif_ids"] == 2
    assert effort_row["median_coordinate_uncertainty_m"] == 15
    assert effort_row["p90_coordinate_uncertainty_m"] == 19
    assert effort_row["n_preserved_specimen"] == 1
    assert effort_row["n_human_observation"] == 1

    taxa_row = taxa.set_index("accepted_species").loc["Aaa aaa"]
    assert taxa_row["genus"] == "Aaa"
    assert taxa_row["family"] == "Fabaceae"
    assert taxa_row["n_islands"] == 1
    assert taxa_row["n_records"] == 2

    assert audit["n_input_rows"] == 3
    assert audit["n_duplicate_rows_removed"] == 1
