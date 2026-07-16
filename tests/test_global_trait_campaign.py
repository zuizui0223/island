import json
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
            "alternative_guild_wikimedia",
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
                "required_for_primary_completion": False,
                "depends_on": ["reproductive_wikimedia"],
                "eligibility": "all",
                "target_traits": ["pollen_vector_mode", "self_incompatibility"],
            },
            "floral_access_wikimedia": {
                "phase": "floral_access",
                "source_kind": "wikimedia_web_reported",
                "depends_on": ["reproductive_wikimedia"],
                "eligibility": "machine_biotic_candidate",
                "target_traits": ["floral_symmetry", "floral_form"],
            },
            "alternative_guild_wikimedia": {
                "phase": "alternative_guild",
                "source_kind": "wikimedia_reported_ecology",
                "depends_on": ["floral_access_wikimedia"],
                "eligibility": "machine_biotic_candidate",
                "target_traits": ["pollination_functional_guild"],
            },
            "alternative_guild_openalex": {
                "phase": "alternative_guild",
                "source_kind": "openalex_reported_ecology",
                "required_for_primary_completion": False,
                "depends_on": ["alternative_guild_wikimedia"],
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


def test_production_config_prioritizes_intrinsic_floral_and_reproductive_traits():
    config = campaign.load_config(Path("config/global_trait_campaign.yml"))

    assert config["scope_config_path"] == "config/angiosperm_scope.yml"
    assert config["expected_species"] == 106295
    assert config["task_order"][0] == "floral_access_wikimedia"
    assert config["tasks"]["floral_access_wikimedia"]["eligibility"] == "all"
    required_tasks = [
        config["tasks"][name]
        for name in config["task_order"]
        if config["tasks"][name].get("required_for_primary_completion", True)
    ]
    required_targets = {
        trait
        for task in required_tasks
        for trait in task.get("target_traits", [])
    }
    assert "flower_primary_color" in required_targets
    assert "floral_form" in required_targets
    assert "floral_symmetry" in required_targets
    assert "mating_system" in required_targets
    assert "self_incompatibility" in required_targets
    assert "pollination_functional_guild" not in required_targets
    assert "pollen_vector_mode" in required_targets
    assert "inflorescence_display" in required_targets
    assert "reward_type" in required_targets
    assert "herkogamy" in required_targets
    assert "alternative_guild_wikimedia" not in config["task_order"]


def test_reproductive_rules_cover_new_priority_traits_without_guild_targeting():
    ontology = campaign.ecology.load_ontology(Path("config/trait_ontology.yml"))
    sources = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "source_text": (
                    "Alpha one has pronounced herkogamy, is mainly selfing, "
                    "and is ambophilous with variable dichogamy."
                ),
                "source_url": "https://example.test/alpha",
                "source_citation": "Example flora",
                "source_type": "flora",
                "evidence_scope": "species_direct",
            }
        ]
    )

    candidates, _ = campaign.ecology.extract_candidates_from_text_sources(
        sources,
        ontology,
        extraction_model="test",
        prompt_id="test",
    )
    observed = set(candidates[["trait_name", "candidate_value"]].itertuples(index=False, name=None))

    assert ("herkogamy", "pronounced") in observed
    assert ("mating_system", "predominantly_selfing") in observed
    assert ("pollen_vector_mode", "mixed") in observed
    assert ("dichogamy", "variable") in observed


def test_production_master_uses_exact_angiosperm_scope_denominator() -> None:
    config = campaign.load_config(Path("config/global_trait_campaign.yml"))
    master = campaign.load_species_master(
        Path("data/v2/staging/gbif/collected/island_taxa.csv"), config
    )

    assert len(master) == 106295
    assert master["accepted_species"].is_unique


def test_family_balanced_batch_spreads_first_wave_across_families():
    ledger = campaign.reconcile_ledger(_master(), None, _config())
    task = campaign.choose_active_task(ledger, _config())
    batch = campaign.family_balanced_batch(ledger, task, 3)

    assert task == "reproductive_wikimedia"
    assert set(batch["family"]) == {"FamA", "FamB", "FamC"}


def test_family_balance_never_crosses_a_higher_priority_geographic_stratum():
    ledger = campaign.reconcile_ledger(_master(), None, _config())
    task = campaign.choose_active_task(ledger, _config())
    ledger["scheduling_region_category"] = [
        "reviewed_nh_temperate",
        "southern_hemisphere",
        "reviewed_nh_temperate",
        "reviewed_nh_temperate",
    ]
    ledger["any_reviewed_NH_temperate"] = ["true", "false", "true", "true"]
    ledger["any_NH_temperate_unreviewed"] = [
        "false",
        "false",
        "false",
        "false",
    ]

    batch = campaign.family_balanced_batch(ledger, task, 3)

    assert batch["accepted_species"].tolist() == ["Alpha one", "Beta one", "Gamma one"]
    assert set(batch["scheduling_region_category"]) == {"reviewed_nh_temperate"}


def test_family_balance_never_crosses_a_higher_priority_trait_stratum():
    ledger = campaign.reconcile_ledger(_master(), None, _config())
    task = campaign.choose_active_task(ledger, _config())
    ledger["scheduling_region_category"] = "nh_temperate_unreviewed"
    ledger["any_reviewed_NH_temperate"] = "false"
    ledger["any_NH_temperate_unreviewed"] = "true"
    ledger["highest_priority_unresolved_rank"] = [2, 1, 2, 1]
    ledger["weighted_unresolved_trait_score"] = [10, 1, 10, 1]

    batch = campaign.family_balanced_batch(ledger, task, 2)

    assert batch["highest_priority_unresolved_rank"].tolist() == [1, 1]


def test_reconcile_reopens_historical_not_applicable_after_task_becomes_all_species():
    config = _config()
    config["tasks"]["floral_access_wikimedia"]["eligibility"] = "all"
    config["tasks"]["floral_access_wikimedia"]["depends_on"] = []
    config["reopen_not_applicable_tasks"] = ["floral_access_wikimedia"]
    old = campaign.reconcile_ledger(_master(), None, _config())
    old["floral_access_wikimedia_status"] = "not_applicable"
    old["floral_access_wikimedia_attempts"] = 2

    reopened = campaign.reconcile_ledger(_master(), old, config)
    assert set(reopened["floral_access_wikimedia_status"]) == {"pending"}
    assert set(reopened["floral_access_wikimedia_attempts"]) == {0}


def test_phase_barrier_and_biotic_gate_are_fail_closed():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    ledger["reproductive_wikimedia_status"] = "processed"
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


def test_campaign_summary_is_json_serializable():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    summary = campaign.campaign_summary(ledger, config)

    encoded = json.dumps(summary)

    assert "reproductive_wikimedia" in encoded


def test_per_species_streaming_dependencies_allow_early_downstream():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    first = ledger["accepted_species"].eq("Alpha one")
    second = ledger["accepted_species"].eq("Alpha two")
    ledger.loc[first, "reproductive_wikimedia_status"] = "processed"
    ledger.loc[first, "machine_biotic_candidate"] = True
    ledger.loc[second, "reproductive_wikimedia_status"] = "processed"
    ledger.loc[second, "machine_biotic_candidate"] = True

    prepared = campaign.prepare_dependent_statuses(ledger, config)
    eligible = campaign.task_eligible_mask(prepared, "floral_access_wikimedia", config)
    batch = campaign.family_balanced_batch(prepared, "floral_access_wikimedia", 10, config)

    assert prepared.loc[first, "floral_access_wikimedia_status"].iloc[0] == "pending"
    assert prepared.loc[second, "floral_access_wikimedia_status"].iloc[0] == "pending"
    assert prepared.loc[eligible, "accepted_species"].tolist() == ["Alpha one", "Alpha two"]
    assert batch["accepted_species"].tolist() == ["Alpha one", "Alpha two"]


def test_nonbiotic_gate_closes_only_after_species_dependencies_finish():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    first = ledger["accepted_species"].eq("Alpha one")
    second = ledger["accepted_species"].eq("Alpha two")
    ledger.loc[first, "reproductive_wikimedia_status"] = "processed"
    ledger.loc[second, "reproductive_wikimedia_status"] = "processed"

    prepared = campaign.prepare_dependent_statuses(ledger, config)

    assert prepared.loc[first, "floral_access_wikimedia_status"].iloc[0] == "not_applicable"
    assert prepared.loc[second, "floral_access_wikimedia_status"].iloc[0] == "not_applicable"


def test_optional_openalex_does_not_block_primary_completion():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    ledger["reproductive_wikimedia_status"] = "processed"
    ledger["floral_access_wikimedia_status"] = "not_applicable"
    ledger["alternative_guild_wikimedia_status"] = "not_applicable"

    assert set(ledger["reproductive_openalex_status"]) == {"pending"}
    assert campaign.choose_active_task(ledger, config) == "complete"


def test_wikimedia_ecology_fallback_runs_only_for_zero_target_species(monkeypatch):
    config = _config()
    config["wikipedia_languages"] = ["en", "fr"]
    ledger = campaign.reconcile_ledger(_master(), None, config)
    batch = ledger.head(2)[["accepted_species", "genus", "family", "n_islands", "n_records"]]
    calls: list[tuple[tuple[str, ...], tuple[str, ...]]] = []

    def fake_getter(**_kwargs):
        return object()

    def fake_wikimedia_text_sources(species_df, _getter, max_taxa, wikipedia_languages=None):
        languages = tuple(wikipedia_languages or ["en"])
        species = tuple(species_df["accepted_species"].astype(str).head(max_taxa))
        calls.append((languages, species))
        if languages == ("en",):
            return (
                pd.DataFrame(
                    [
                        {
                            "accepted_species": species[0],
                            "source_text": f"{species[0]} is pollinated by bees.",
                            "source_url": "https://en.wikipedia.org/wiki/Alpha_one",
                            "source_citation": "English source",
                            "source_type": "wikipedia_extract",
                            "evidence_scope": "species_direct",
                        },
                        {
                            "accepted_species": species[1],
                            "source_text": f"{species[1]} is a shrub.",
                            "source_url": "https://en.wikipedia.org/wiki/Alpha_two",
                            "source_citation": "English source",
                            "source_type": "wikipedia_extract",
                            "evidence_scope": "species_direct",
                        },
                    ]
                ),
                [],
            )
        return (
            pd.DataFrame(
                [
                    {
                        "accepted_species": species[0],
                        "source_text": f"{species[0]} est autocompatible.",
                        "source_url": "https://fr.wikipedia.org/wiki/Alpha_two",
                        "source_citation": "French source",
                        "source_type": "wikipedia_extract",
                        "evidence_scope": "species_direct",
                    }
                ]
            ),
            [],
        )

    monkeypatch.setattr(campaign.web_reported, "_httpx_getter", fake_getter)
    monkeypatch.setattr(campaign.ecology, "wikimedia_text_sources", fake_wikimedia_text_sources)

    candidates, holdouts, errors, source_species = campaign.fetch_wikimedia_ecology_candidates(
        batch,
        "reproductive_wikimedia",
        config,
        campaign.ecology.load_ontology(Path("config/trait_ontology.yml")),
    )

    assert calls == [
        (("en",), ("Alpha one", "Alpha two")),
        (("fr",), ("Alpha two",)),
    ]
    assert errors == []
    assert holdouts.empty
    assert source_species == {"Alpha one", "Alpha two"}
    values = set(candidates[["trait_name", "candidate_value"]].itertuples(index=False, name=None))
    assert ("pollen_vector_mode", "biotic") in values
    assert ("self_incompatibility", "SC") in values
