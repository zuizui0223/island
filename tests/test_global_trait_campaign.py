from pathlib import Path

import pandas as pd

from island_v2 import global_trait_campaign as campaign


def _config() -> dict:
    return {
        "version": "test",
        "default_batch_size": 3,
        "max_attempts": 2,
        "task_order": [
            "reproductive_wikimedia",
            "reproductive_openalex",
            "floral_access_wikimedia",
            "alternative_guild_openalex",
        ],
        "tasks": {
            "reproductive_wikimedia": {
                "phase": "reproductive",
                "source_kind": "wikimedia_reported_ecology",
                "depends_on": [],
                "eligibility": "all",
                "target_traits": ["pollen_vector_mode", "self_incompatibility"],
            },
            "reproductive_openalex": {
                "phase": "reproductive",
                "source_kind": "openalex_reported_ecology",
                "depends_on": ["reproductive_wikimedia"],
                "eligibility": "all",
                "target_traits": ["pollen_vector_mode", "self_incompatibility"],
            },
            "floral_access_wikimedia": {
                "phase": "floral_access",
                "source_kind": "wikimedia_web_reported",
                "depends_on": ["reproductive_wikimedia", "reproductive_openalex"],
                "eligibility": "machine_biotic_candidate",
                "target_traits": ["floral_symmetry", "floral_form"],
            },
            "alternative_guild_openalex": {
                "phase": "alternative_guild",
                "source_kind": "openalex_reported_ecology",
                "depends_on": ["floral_access_wikimedia"],
                "eligibility": "machine_biotic_candidate",
                "target_traits": ["pollination_functional_guild"],
            },
        },
        "excluded_families": ["Pinaceae"],
    }


def _master() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "genus": "Alpha",
                "family": "FamA",
                "n_islands": 5,
                "n_records": 10,
            },
            {
                "accepted_species": "Alpha two",
                "genus": "Alpha",
                "family": "FamA",
                "n_islands": 4,
                "n_records": 9,
            },
            {
                "accepted_species": "Beta one",
                "genus": "Beta",
                "family": "FamB",
                "n_islands": 3,
                "n_records": 8,
            },
            {
                "accepted_species": "Gamma one",
                "genus": "Gamma",
                "family": "FamC",
                "n_islands": 2,
                "n_records": 7,
            },
        ]
    )


def test_family_balanced_batch_spreads_first_wave_across_families():
    ledger = campaign.reconcile_ledger(_master(), None, _config())
    task = campaign.choose_active_task(ledger, _config())
    batch = campaign.family_balanced_batch(ledger, task, 3)

    assert task == "reproductive_wikimedia"
    assert set(batch["family"]) == {"FamA", "FamB", "FamC"}


def test_phase_barrier_and_biotic_gate_are_fail_closed():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    for task in ("reproductive_wikimedia", "reproductive_openalex"):
        ledger[f"{task}_status"] = "processed"
    ledger.loc[ledger["accepted_species"].eq("Alpha one"), "machine_biotic_candidate"] = True

    prepared = campaign.prepare_dependent_statuses(ledger, config)

    assert campaign.choose_active_task(prepared, config) == "floral_access_wikimedia"
    non_biotic = prepared.loc[~prepared["machine_biotic_candidate"]]
    assert set(non_biotic["floral_access_wikimedia_status"]) == {"not_applicable"}
    batch = campaign.family_balanced_batch(prepared, "floral_access_wikimedia", 10)
    assert batch["accepted_species"].tolist() == ["Alpha one"]


def test_apply_result_marks_zero_hit_as_processed_and_records_direct_biotic():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    batch = ledger.head(2)[["accepted_species", "genus", "family", "n_islands", "n_records"]]
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": batch.iloc[0]["accepted_species"],
                "trait_name": "pollen_vector_mode",
                "candidate_value": "biotic",
                "source_url": "https://example.test",
                "source_citation": "Example",
                "source_excerpt": "Pollinated by bees.",
                "source_type": "paper",
                "evidence_scope": "species_direct",
                "matched_term": "pollinated by bees",
                "confidence": "rule_match",
                "campaign_task": "reproductive_wikimedia",
                "campaign_phase": "reproductive",
                "target_for_task": True,
            }
        ],
        columns=campaign.COMMON_CANDIDATE_COLUMNS,
    )

    updated = campaign.apply_task_result(
        ledger,
        batch,
        candidates,
        errors=[],
        source_species={batch.iloc[0]["accepted_species"]},
        task="reproductive_wikimedia",
        wave_id="wave_00001_reproductive_wikimedia",
        config=config,
    )

    attempted = updated.loc[updated["accepted_species"].isin(batch["accepted_species"])]
    assert set(attempted["reproductive_wikimedia_status"]) == {"processed"}
    assert attempted["reproductive_wikimedia_candidate_count"].sum() == 1
    assert updated.loc[
        updated["accepted_species"].eq(batch.iloc[0]["accepted_species"]),
        "machine_biotic_candidate",
    ].iloc[0]


def test_reconcile_adds_new_species_without_resetting_existing_state():
    config = _config()
    ledger = campaign.reconcile_ledger(_master().head(1), None, config)
    ledger.loc[0, "reproductive_wikimedia_status"] = "processed"
    expanded = campaign.reconcile_ledger(_master().head(2), ledger, config)

    old = expanded.loc[expanded["accepted_species"].eq("Alpha one")].iloc[0]
    new = expanded.loc[expanded["accepted_species"].eq("Alpha two")].iloc[0]
    assert old["reproductive_wikimedia_status"] == "processed"
    assert new["reproductive_wikimedia_status"] == "pending"


def test_load_species_master_excludes_obvious_non_angiosperms(tmp_path: Path):
    master = _master()
    master.loc[len(master)] = {
        "accepted_species": "Pinus example",
        "genus": "Pinus",
        "family": "Pinaceae",
        "n_islands": 1,
        "n_records": 1,
    }
    path = tmp_path / "master.csv"
    master.to_csv(path, index=False)

    loaded = campaign.load_species_master(path, _config())

    assert "Pinus example" not in set(loaded["accepted_species"])
