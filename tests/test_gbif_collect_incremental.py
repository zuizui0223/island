import importlib

import pandas as pd

collector = importlib.import_module("island_v2.gbif_collect_incremental")


def _campaign() -> dict[str, object]:
    return {
        "ledger": [
            {"block_id": "large", "request_status": "succeeded", "download_key": "k_large"},
            {"block_id": "small", "request_status": "succeeded", "download_key": "k_small"},
            {"block_id": "running", "request_status": "running", "download_key": "k_running"},
        ]
    }


def test_select_batch_skips_completed_blocks_and_starts_with_smaller_member_group():
    manifest = pd.DataFrame(
        {
            "block_id": ["small"],
            "collection_status": ["collected"],
        }
    )
    selected = collector.select_batch(
        _campaign(), manifest, {"small": 1, "large": 100}, max_blocks=1
    )

    assert [entry["block_id"] for entry in selected] == ["large"]


def test_select_batch_prefers_smallest_unprocessed_block():
    selected = collector.select_batch(
        _campaign(), pd.DataFrame(), {"small": 1, "large": 100}, max_blocks=1
    )

    assert [entry["block_id"] for entry in selected] == ["small"]


def test_merge_species_and_taxa_accumulate_across_nonoverlapping_blocks():
    old_species = pd.DataFrame(
        {
            "island_id": ["a"],
            "species": ["Plant one"],
            "n_records": [2],
            "n_unique_gbif_ids": [2],
            "basis_of_record_set": ["PRESERVED_SPECIMEN"],
            "review_status": ["unresolved_taxonomy"],
        }
    )
    new_species = pd.DataFrame(
        {
            "island_id": ["b"],
            "species": ["Plant one"],
            "n_records": [3],
            "n_unique_gbif_ids": [3],
            "basis_of_record_set": ["HUMAN_OBSERVATION"],
            "review_status": ["unresolved_taxonomy"],
        }
    )
    species = collector._merge_species(old_species, new_species)

    assert len(species) == 2
    old_taxa = pd.DataFrame(
        {
            "accepted_species": ["Plant one"],
            "genus": ["Plant"],
            "family": ["Family"],
            "n_islands": [1],
            "n_records": [2],
        }
    )
    new_taxa = old_taxa.assign(n_islands=1, n_records=3)
    taxa = collector._merge_taxa(old_taxa, new_taxa).iloc[0]

    assert taxa["n_islands"] == 2
    assert taxa["n_records"] == 5


def test_merge_effort_rejects_overlapping_islands():
    effort = pd.DataFrame({"island_id": ["same"], "n_records": [1]})

    try:
        collector._merge_effort(effort, effort)
    except Exception as exc:
        assert "overlap" in str(exc)
    else:
        raise AssertionError("Expected overlapping island block products to be rejected")


def test_gzip_species_snapshot_round_trip_and_precedence(tmp_path):
    compressed = tmp_path / collector.CUMULATIVE_SPECIES_SNAPSHOT
    legacy = tmp_path / collector.LEGACY_SPECIES_SNAPSHOT
    expected = pd.DataFrame(
        {
            "island_id": ["island-a"],
            "species": ["Plant one"],
            "n_records": [2],
        }
    )
    collector._write_csv(expected, compressed)
    legacy.write_text("island_id,species,n_records\nlegacy,Legacy plant,1\n", encoding="utf-8")

    restored = collector._read_cumulative_species(tmp_path)

    assert restored.to_dict("records") == [
        {"island_id": "island-a", "species": "Plant one", "n_records": "2"}
    ]
