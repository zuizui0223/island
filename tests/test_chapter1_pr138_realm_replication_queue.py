import pandas as pd

from island_v2.chapter1_pr138_realm_replication_queue import (
    CORE_SYNDROMES,
    build_realm_replication_queue,
)


def test_queue_is_outcome_blind_and_excludes_already_resolved_islands():
    status = pd.DataFrame(
        {
            "island_id": ["resolved", "resolved", "candidate", "candidate", "weak"],
            "accepted_species": ["a", "b", "a", "b", "c"],
            "origin_status": ["native", "native", "unresolved", "unresolved", "unresolved"],
            "status_resolved": [True, True, False, False, False],
        }
    )
    realm = pd.DataFrame(
        {
            "island_id": ["resolved", "candidate", "weak"],
            "biogeographic_realm": ["Nearctic", "Nearctic", "Nearctic"],
            "island_latitude": [40.0, 41.0, 42.0],
            "island_longitude": [-70.0, -71.0, -72.0],
        }
    )
    syndrome_rows = []
    for species in ["a", "b"]:
        for syndrome in CORE_SYNDROMES:
            syndrome_rows.append(
                {
                    "accepted_species": species,
                    "syndrome": syndrome,
                    "syndrome_concordance": 999.0,
                    "soft_membership": 999.0,
                }
            )
    syndrome = pd.DataFrame(syndrome_rows)
    observation = pd.DataFrame(
        {
            "island_id": ["resolved", "candidate", "weak"],
            "flora_recorded": [True, True, True],
            "n_flora_species_recorded": [100, 80, 10],
            "distance_to_continent_km": [20.0, 40.0, 60.0],
            "area_km2": [10.0, 20.0, 30.0],
            "spatial_block": ["a", "b", "c"],
        }
    )

    queue = build_realm_replication_queue(
        status,
        realm,
        syndrome,
        observation,
        target_realm="Nearctic",
        min_flora_species=20,
        max_flora_species=400,
        min_species_per_syndrome=1,
        min_syndromes_meeting_support=5,
    )

    assert list(queue["island_id"]) == ["candidate"]
    assert queue.loc[0, "n_syndromes_meeting_support"] == 5
    assert queue.loc[0, "min_syndrome_species"] == 2
    assert "syndrome_concordance" not in queue.columns
    assert "soft_membership" not in queue.columns
    assert not any("effect" in column for column in queue.columns)
    assert not any("p_value" in column for column in queue.columns)
