import pandas as pd
import pytest
from typer.testing import CliRunner

from island_v2.local_reproductive_assurance_collection import (
    _explicit_mapping_supported,
    PROVIDERS,
    RA_TRAITS,
    app,
    augment_cleistogamy_cache_candidates,
    bootstrap_provider_checkpoint,
    build_reviewed_search_batch,
    materialize_reviewed_cache,
    refresh_coverage_summary,
)


def legacy_rows(species: list[str]) -> pd.DataFrame:
    rows = []
    for name in species:
        for provider in PROVIDERS:
            enabled = provider not in {"openalex", "world_flora"}
            status = "pending" if not enabled else "no_hit"
            attempts = 0 if not enabled else 1
            if provider == "europe_pmc" and name == "Alpha beta":
                status = "hit"
            rows.append(
                {
                    "species": name,
                    "provider": provider,
                    "policy_version": f"policy:{provider}",
                    "enabled": str(enabled).lower(),
                    "status": status,
                    "attempts": attempts,
                    "last_packet_id": "",
                    "last_error": "",
                    "updated_at": "2026-07-17T00:00:00+00:00",
                }
            )
    return pd.DataFrame(rows)


def axis_rows(species: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": name,
                "genus": name.split()[0],
                "family": "Exampleaceae",
                "colour_vividness_evidence_traits": "flower_primary_color",
                "colour_vividness_direct": True,
                "floral_complexity_evidence_traits": "floral_form",
                "floral_complexity_direct": True,
                "reproductive_assurance_evidence_traits": "",
                "reproductive_assurance_direct": False,
                "all_three_axes_direct": False,
                "missing_direct_axes": "reproductive_assurance",
                "reproductive_assurance_priority": True,
                "provider_state_available": False,
                "checkpoint_scope": "materialized_axis_baseline_only",
            }
            for name in species
        ]
    )


def candidate_row(
    candidate_id: str,
    species: str,
    *,
    trait: str = "self_incompatibility",
    value: str = "SI",
    raw: str = "Alpha beta is self-incompatible.",
) -> dict[str, str]:
    return {
        "candidate_id": candidate_id,
        "accepted_species": species,
        "trait_name": trait,
        "provisional_candidate_value": value,
        "matched_term": "self-incompatible",
        "candidate_class": "reported",
        "provider": "europe_pmc",
        "source_type": "europe_pmc_title_abstract",
        "source_url": "https://europepmc.org/article/MED/1",
        "source_citation": (
            "A Author. (2024). Breeding system of Alpha beta. Example Journal. "
            "DOI:10.1000/example; PMID:1"
        ),
        "raw_description": raw,
        "source_text": f"Breeding system of {species}. {raw}",
        "evidence_scope": "species_direct",
        "wild_or_cultivated": "wild",
        "review_status": "unreviewed",
        "provider_status": "hit",
        "provider_attempts": "1",
        "provider_updated_at": "2026-07-17T00:00:00Z",
        "provider_policy_version": "policy",
        "local_file_path": "/tmp/cache.csv",
        "local_file_hash": "a" * 64,
        "retrieval_date": "2026-07-17",
        "retrieval_date_basis": "provider_updated_at",
        "source_run_id": "123",
        "audit_status": "candidate",
        "audit_reason": "explicit term",
    }


def test_bootstrap_has_exact_species_trait_provider_key():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta", "Gamma delta"]})
    checkpoint = bootstrap_provider_checkpoint(queue, legacy_rows(queue.accepted_species.tolist()))
    assert len(checkpoint) == 2 * len(RA_TRAITS) * len(PROVIDERS)
    assert not checkpoint.duplicated(["species", "trait", "provider"]).any()


def test_europe_pmc_query_completion_expands_without_fabricating_trait_hit():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidates = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "trait_name": "self_incompatibility",
                "provider": "europe_pmc",
            }
        ]
    )
    checkpoint = bootstrap_provider_checkpoint(
        queue,
        legacy_rows(["Alpha beta"]),
        candidates=candidates,
    )
    europe = checkpoint.loc[checkpoint["provider"].eq("europe_pmc")].set_index("trait")
    assert europe.loc["self_incompatibility", "status"] == (
        "provider_pass_complete_with_candidates"
    )
    assert europe.loc["autonomous_selfing_capacity", "status"] == "provider_pass_complete"
    assert europe.loc["cleistogamy", "status"] == "provider_pass_complete"
    assert europe["terminal"].astype(bool).all()


def test_provider_scope_does_not_expand_web_description_to_autonomous_traits():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    checkpoint = bootstrap_provider_checkpoint(queue, legacy_rows(["Alpha beta"]))
    row = checkpoint.loc[
        checkpoint["provider"].eq("web_descriptions")
        & checkpoint["trait"].eq("autonomous_selfing_capacity")
    ].iloc[0]
    assert row["status"] == "not_applicable"
    assert row["migration_basis"] == "provider_does_not_query_trait"


def test_openalex_global_status_preserves_processed_and_transient_retry():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta", "Gamma delta"]})
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "reproductive_openalex_status": "processed",
                "reproductive_openalex_attempts": "1",
                "reproductive_openalex_wave_id": "wave_1",
                "reproductive_openalex_last_error": "",
            },
            {
                "accepted_species": "Gamma delta",
                "reproductive_openalex_status": "retry",
                "reproductive_openalex_attempts": "1",
                "reproductive_openalex_wave_id": "wave_2",
                "reproductive_openalex_last_error": "429",
            },
        ]
    )
    checkpoint = bootstrap_provider_checkpoint(
        queue,
        legacy_rows(queue.accepted_species.tolist()),
        global_ledger=ledger,
    )
    openalex = checkpoint.loc[checkpoint["provider"].eq("openalex")]
    assert set(openalex.loc[openalex["species"].eq("Alpha beta"), "status"]) == {
        "provider_pass_complete"
    }
    assert set(openalex.loc[openalex["species"].eq("Gamma delta"), "status"]) == {"retry"}


def test_one_token_taxon_is_terminal_and_never_queried():
    queue = pd.DataFrame({"accepted_species": ["Agave"]})
    checkpoint = bootstrap_provider_checkpoint(queue, legacy_rows(["Agave"]))
    applicable = checkpoint.loc[checkpoint["provider"].isin({"europe_pmc", "openalex"})]
    assert set(applicable["status"]) == {"invalid_species_name"}
    assert applicable["terminal"].astype(bool).all()


def test_one_token_taxon_preserves_prior_attempt_provenance():
    queue = pd.DataFrame({"accepted_species": ["Agave"]})
    legacy = legacy_rows(["Agave"])
    mask = legacy["provider"].eq("europe_pmc")
    legacy.loc[mask, "attempts"] = 3
    legacy.loc[mask, "last_error"] = "legacy terminal failure"
    checkpoint = bootstrap_provider_checkpoint(queue, legacy)
    row = checkpoint.loc[
        checkpoint["provider"].eq("europe_pmc") & checkpoint["trait"].eq("self_incompatibility")
    ].iloc[0]
    assert row["attempts"] == 3
    assert "legacy terminal failure" in row["last_error"]
    assert row["updated_at"] == "2026-07-17T00:00:00+00:00"


def test_legacy_checkpoint_must_cover_every_provider():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    legacy = legacy_rows(["Alpha beta"]).iloc[:-1]
    with pytest.raises(ValueError, match="exact queue and provider set"):
        bootstrap_provider_checkpoint(queue, legacy)


def test_bootstrap_is_an_explicit_cli_subcommand():
    result = CliRunner().invoke(app, ["bootstrap", "--help"])
    assert result.exit_code == 0
    assert "Create the local trait-scoped ledger" in result.output


def test_blank_attempts_are_treated_as_zero():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    legacy = legacy_rows(["Alpha beta"])
    legacy["attempts"] = legacy["attempts"].astype(object)
    legacy.loc[legacy["provider"].eq("europe_pmc"), "attempts"] = ""
    checkpoint = bootstrap_provider_checkpoint(queue, legacy)
    assert set(checkpoint.loc[checkpoint["provider"].eq("europe_pmc"), "attempts"]) == {0}


def test_openalex_wave_id_is_not_written_as_a_timestamp():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "reproductive_openalex_status": "processed",
                "reproductive_openalex_attempts": "1",
                "reproductive_openalex_wave_id": "wave_1",
                "reproductive_openalex_last_error": "",
            }
        ]
    )
    checkpoint = bootstrap_provider_checkpoint(
        queue,
        legacy_rows(["Alpha beta"]),
        global_ledger=ledger,
    )
    openalex = checkpoint.loc[checkpoint["provider"].eq("openalex")]
    assert set(openalex["updated_at"]) == {""}
    assert set(openalex["last_wave_id"]) == {"wave_1"}


def test_reviewed_cache_materialization_updates_all_requested_outputs():
    species = ["Alpha beta", "Gamma delta"]
    queue = pd.DataFrame({"accepted_species": species})
    candidates = pd.DataFrame(
        [candidate_row("accepted", "Alpha beta"), candidate_row("rejected", "Gamma delta")]
    )
    providers = bootstrap_provider_checkpoint(
        queue,
        legacy_rows(species),
        candidates=candidates,
    )
    # Persisted checkpoints are restored as strings; materialization must
    # normalize numeric and Boolean columns before updating them.
    providers = providers.astype("str")
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=axis_rows(species).astype("str"),
        provider_checkpoint=providers,
        candidates=candidates,
        review_payload={
            "accepted": [
                {
                    "candidate_id": "accepted",
                    "accepted_value": "SI",
                    "review_reason": "Exact species is explicitly self-incompatible.",
                }
            ]
        },
    )
    assert result["records"][["species", "trait", "value"]].to_dict("records") == [
        {"species": "Alpha beta", "trait": "self_incompatibility", "value": "SI"}
    ]
    assert result["records"].iloc[0]["DOI"] == "10.1000/example"
    assert result["records"].iloc[0]["year"] == "2024"
    assert set(result["next_queue"]["accepted_species"]) == {"Gamma delta"}
    assert result["summary"]["newly_completed_species"] == 1
    assert result["summary"]["cumulative_newly_completed_species"] == 1
    assert result["summary"]["latest_batch_newly_completed_species"] == 1
    assert result["summary"]["total_three_axis_completion"] == 1
    assert result["summary"]["cache_accepted_candidate_count"] == 1
    assert result["summary"]["cache_accepted_trait_record_count"] == 1
    assert result["summary"]["cache_rejected_candidate_count"] == 1
    assert result["summary"]["new_record_count"] == 1
    assert result["summary"]["total_record_count"] == 1
    accepted_provider = (
        result["provider_checkpoint"]
        .loc[
            result["provider_checkpoint"]["species"].eq("Alpha beta")
            & result["provider_checkpoint"]["trait"].eq("self_incompatibility")
            & result["provider_checkpoint"]["provider"].eq("europe_pmc")
        ]
        .iloc[0]
    )
    assert accepted_provider["status"] == "accepted_evidence"
    assert accepted_provider["accepted_record_count"] == 1
    alpha_provider_rows = result["provider_checkpoint"].loc[
        result["provider_checkpoint"]["species"].eq("Alpha beta")
    ]
    assert alpha_provider_rows["terminal"].all()
    skipped = alpha_provider_rows.loc[alpha_provider_rows["status"].eq("provider_seed_covered")]
    assert not skipped.empty
    assert (
        skipped["migration_basis"]
        .str.contains(
            "reproductive_assurance_completed_by_another_provider",
            regex=False,
        )
        .all()
    )
    gamma_provider_rows = result["provider_checkpoint"].loc[
        result["provider_checkpoint"]["species"].eq("Gamma delta")
    ]
    assert (~gamma_provider_rows["terminal"]).any()
    assert result["review_audit"]["review_decision"].value_counts().to_dict() == {
        "accepted": 1,
        "rejected": 1,
    }


def test_completion_summary_distinguishes_cumulative_and_latest_batch():
    species = ["Alpha beta", "Gamma delta"]
    queue = pd.DataFrame({"accepted_species": species})
    baseline = axis_rows(species)
    providers = bootstrap_provider_checkpoint(queue, legacy_rows(species))

    alpha_candidate = pd.DataFrame([candidate_row("alpha", "Alpha beta")])
    first = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=baseline,
        provider_checkpoint=providers,
        candidates=alpha_candidate,
        review_payload={
            "accepted": [
                {
                    "candidate_id": "alpha",
                    "accepted_value": "SI",
                    "review_reason": "Exact species is explicitly self-incompatible.",
                }
            ]
        },
    )

    gamma_candidate = pd.DataFrame(
        [
            candidate_row(
                "gamma",
                "Gamma delta",
                value="SC",
                raw="Gamma delta is self-compatible.",
            )
        ]
    )
    second = materialize_reviewed_cache(
        queue=first["next_queue"],
        baseline_checkpoint=first["three_axis_checkpoint"],
        provider_checkpoint=first["provider_checkpoint"],
        candidates=gamma_candidate,
        review_payload={
            "accepted": [
                {
                    "candidate_id": "gamma",
                    "accepted_value": "SC",
                    "review_reason": "Exact species is explicitly self-compatible.",
                }
            ]
        },
        prior_records=first["records"],
    )
    summary = second["summary"]
    assert summary["newly_completed_species"] == 2
    assert summary["newly_completed_species_names"] == species
    assert summary["cumulative_newly_completed_species"] == 2
    assert summary["latest_batch_newly_completed_species"] == 1
    assert summary["latest_batch_newly_completed_species_names"] == ["Gamma delta"]

    refreshed = refresh_coverage_summary(
        prior_summary=summary,
        axes=second["three_axis_checkpoint"],
        queue=second["next_queue"],
        providers=second["provider_checkpoint"],
        records=second["records"],
        acquisition_status_paths=[],
        output_paths=[],
    )
    assert refreshed["newly_completed_species"] == 2
    assert refreshed["cumulative_newly_completed_species"] == 2
    assert refreshed["latest_batch_newly_completed_species"] == 1


def test_zero_yield_materialization_does_not_claim_an_entire_empty_pass():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidates = pd.DataFrame([candidate_row("held", "Alpha beta")])
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=axis_rows(["Alpha beta"]),
        provider_checkpoint=bootstrap_provider_checkpoint(
            queue,
            legacy_rows(["Alpha beta"]),
        ),
        candidates=candidates,
        review_payload={"accepted": []},
    )
    assert result["summary"]["new_record_count"] == 0
    assert not result["summary"]["no_new_direct_evidence_in_entire_pass"]


def test_explicit_negative_autonomous_statement_can_correct_machine_candidate():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidate = candidate_row(
        "negative",
        "Alpha beta",
        trait="autonomous_selfing_capacity",
        value="autonomous",
        raw="Automatic self-pollination did not occur in Alpha beta.",
    )
    candidates = pd.DataFrame([candidate])
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=axis_rows(["Alpha beta"]),
        provider_checkpoint=bootstrap_provider_checkpoint(
            queue,
            legacy_rows(["Alpha beta"]),
            candidates=candidates,
        ),
        candidates=candidates,
        review_payload={
            "accepted": [
                {
                    "candidate_id": "negative",
                    "accepted_value": "absent",
                    "review_reason": "The experiment explicitly states that it did not occur.",
                }
            ]
        },
    )
    assert result["records"].iloc[0]["value"] == "absent"


def test_nonexplicit_value_override_is_rejected():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidates = pd.DataFrame([candidate_row("candidate", "Alpha beta")])
    with pytest.raises(ValueError, match="may override a machine mapping only"):
        materialize_reviewed_cache(
            queue=queue,
            baseline_checkpoint=axis_rows(["Alpha beta"]),
            provider_checkpoint=bootstrap_provider_checkpoint(
                queue,
                legacy_rows(["Alpha beta"]),
                candidates=candidates,
            ),
            candidates=candidates,
            review_payload={
                "accepted": [
                    {
                        "candidate_id": "candidate",
                        "accepted_value": "SC",
                        "review_reason": "unsupported override",
                    }
                ]
            },
        )


def test_undated_source_is_rejected_as_incomplete_provenance():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidate = candidate_row("undated", "Alpha beta")
    candidate["source_citation"] = "Alpha beta species account"
    candidates = pd.DataFrame([candidate])
    with pytest.raises(ValueError, match="lacks a dated publication"):
        materialize_reviewed_cache(
            queue=queue,
            baseline_checkpoint=axis_rows(["Alpha beta"]),
            provider_checkpoint=bootstrap_provider_checkpoint(
                queue,
                legacy_rows(["Alpha beta"]),
                candidates=candidates,
            ),
            candidates=candidates,
            review_payload={
                "accepted": [
                    {
                        "candidate_id": "undated",
                        "accepted_value": "SI",
                        "review_reason": "Explicit but undated species-level statement.",
                    }
                ]
            },
        )


def test_method_only_self_pollination_does_not_establish_autonomy():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidate = candidate_row(
        "method-only",
        "Alpha beta",
        trait="autonomous_selfing_capacity",
        value="autonomous",
        raw="Researchers self-pollinated flowers by hand before bagging them.",
    )
    candidates = pd.DataFrame([candidate])
    with pytest.raises(ValueError, match="lacks an explicit source statement"):
        materialize_reviewed_cache(
            queue=queue,
            baseline_checkpoint=axis_rows(["Alpha beta"]),
            provider_checkpoint=bootstrap_provider_checkpoint(
                queue,
                legacy_rows(["Alpha beta"]),
                candidates=candidates,
            ),
            candidates=candidates,
            review_payload={
                "accepted": [
                    {
                        "candidate_id": "method-only",
                        "accepted_value": "autonomous",
                        "review_reason": "Method language alone must not pass.",
                    }
                ]
            },
        )


def test_seed_production_through_self_pollination_supports_compatibility_only():
    text = (
        "Three hundred and forty-six seeds were obtained through self-pollination from 17 clones."
    )
    assert _explicit_mapping_supported("self_incompatibility", "SC", text)
    assert not _explicit_mapping_supported("autonomous_selfing_capacity", "autonomous", text)


@pytest.mark.parametrize(
    "text",
    [
        (
            "Five inbred lines of Cucumis sativus were generated by "
            "self-pollination across six generations."
        ),
        (
            'The inbred line "112-2" was purified over 8 consecutive '
            "generations of self-pollination."
        ),
        ("Populations were obtained from lines with one to four generations of self-pollination."),
    ],
)
def test_repeated_successful_self_pollination_supports_compatibility_only(text):
    assert _explicit_mapping_supported("self_incompatibility", "SC", text)
    assert not _explicit_mapping_supported("autonomous_selfing_capacity", "autonomous", text)


def test_spontaneous_self_pollination_with_fruit_set_supports_sc_and_autonomy():
    text = "Of the flowers submitted to the spontaneous self pollination test, 39% set fruits."
    assert _explicit_mapping_supported("self_incompatibility", "SC", text)
    assert _explicit_mapping_supported("autonomous_selfing_capacity", "autonomous", text)


def test_explicit_self_fertilization_and_cross_pollination_support_mixed_mating():
    text = (
        "Eugenia uniflora reproduction is based on both self-fertilization "
        "as well as cross-pollination."
    )
    assert _explicit_mapping_supported("mating_system", "mixed_mating", text)


def test_facultative_xenogamous_species_supports_predominant_outcrossing():
    text = "Martynia annua is a facultative xenogamous species."
    assert _explicit_mapping_supported("mating_system", "predominantly_outcrossing", text)


def test_one_source_may_support_multiple_explicit_traits():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidate = candidate_row(
        "multi",
        "Alpha beta",
        trait="autonomous_selfing_capacity",
        value="autonomous",
        raw="Alpha beta is autogamous and predominantly outcrossing.",
    )
    candidates = pd.DataFrame([candidate])
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=axis_rows(["Alpha beta"]),
        provider_checkpoint=bootstrap_provider_checkpoint(
            queue,
            legacy_rows(["Alpha beta"]),
            candidates=candidates,
        ),
        candidates=candidates,
        review_payload={
            "accepted": [
                {
                    "candidate_id": "multi",
                    "accepted_trait": "autonomous_selfing_capacity",
                    "accepted_value": "autonomous",
                    "review_reason": "Explicit autogamy statement.",
                },
                {
                    "candidate_id": "multi",
                    "accepted_trait": "mating_system",
                    "accepted_value": "predominantly_outcrossing",
                    "review_reason": "Explicit predominantly-outcrossing statement.",
                },
            ]
        },
    )
    assert set(result["records"]["trait"]) == {
        "autonomous_selfing_capacity",
        "mating_system",
    }
    assert result["review_audit"].iloc[0]["accepted_mappings"] == (
        "autonomous_selfing_capacity=autonomous|mating_system=predominantly_outcrossing"
    )
    assert result["summary"]["cache_accepted_candidate_count"] == 1
    assert result["summary"]["cache_accepted_trait_record_count"] == 2
    assert result["summary"]["cache_rejected_candidate_count"] == 0


def test_exact_autonomous_absence_and_outcrossing_phrases_can_correct_mapping():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    candidate = candidate_row(
        "multi-correction",
        "Alpha beta",
        trait="self_incompatibility",
        value="SC",
        raw=(
            "Plants were self-compatible but did not self-pollinate autonomously; "
            "the study describes their mating system as outcrossing."
        ),
    )
    candidates = pd.DataFrame([candidate])
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=axis_rows(["Alpha beta"]),
        provider_checkpoint=bootstrap_provider_checkpoint(
            queue,
            legacy_rows(["Alpha beta"]),
            candidates=candidates,
        ),
        candidates=candidates,
        review_payload={
            "accepted": [
                {
                    "candidate_id": "multi-correction",
                    "accepted_trait": "autonomous_selfing_capacity",
                    "accepted_value": "absent",
                    "review_reason": "Explicit autonomous self-pollination absence.",
                },
                {
                    "candidate_id": "multi-correction",
                    "accepted_trait": "mating_system",
                    "accepted_value": "predominantly_outcrossing",
                    "review_reason": "Explicit outcrossing mating-system statement.",
                },
            ]
        },
    )
    assert set(result["records"][["trait", "value"]].itertuples(index=False, name=None)) == {
        ("autonomous_selfing_capacity", "absent"),
        ("mating_system", "predominantly_outcrossing"),
    }


def test_cleistogamy_augmentation_surfaces_review_candidates_without_accepting_them():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    search = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "source_text": ("Alpha beta produces both cleistogamous flowers and open flowers."),
                "source_url": "https://example.org/alpha",
                "source_citation": "A Author. (2024). Alpha beta flowers. Journal.",
                "source_type": "flora_or_monograph",
                "evidence_scope": "species_direct",
            }
        ]
    )
    legacy = legacy_rows(["Alpha beta"])
    candidates = pd.DataFrame([candidate_row("existing", "Alpha beta")])
    augmented = augment_cleistogamy_cache_candidates(
        queue=queue,
        search_evidence=search,
        legacy_provider=legacy,
        existing_candidates=candidates,
        local_file_path="/tmp/search.csv",
        local_file_hash="b" * 64,
    )
    cleistogamy = augmented.loc[augmented["trait_name"].eq("cleistogamy")].iloc[0]
    assert cleistogamy["provisional_candidate_value"] == "facultative"
    assert cleistogamy["provider"] == "web_descriptions"
    assert cleistogamy["review_status"] == "unreviewed"
    assert len(augmented) == 2


def test_reviewed_search_batch_reconstructs_provenance_and_decision():
    queue = pd.DataFrame({"accepted_species": ["Alpha beta"]})
    source_text = "Alpha beta is a primarily outcrossing species."
    search = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "source_text": source_text,
                "source_url": "https://europepmc.org/article/MED/1",
                "source_citation": (
                    "A Author. (2024). Mating system of Alpha beta. Journal. DOI:10.1000/example"
                ),
                "source_type": "europe_pmc_title_abstract",
                "evidence_scope": "species_direct",
            }
        ]
    )
    providers = bootstrap_provider_checkpoint(queue, legacy_rows(["Alpha beta"]))
    candidates, decisions = build_reviewed_search_batch(
        queue=queue,
        search_evidence=search,
        provider_checkpoint=providers,
        selection_payload={
            "accepted": [
                {
                    "source_row_index": 0,
                    "accepted_species": "Alpha beta",
                    "source_url": "https://europepmc.org/article/MED/1",
                    "accepted_trait": "mating_system",
                    "accepted_value": "predominantly_outcrossing",
                    "evidence_excerpt": source_text,
                    "review_reason": "Exact species-level outcrossing statement.",
                }
            ]
        },
        local_file_path="/tmp/search.csv",
        local_file_hash="b" * 64,
    )
    assert len(candidates) == 1
    assert candidates.iloc[0]["local_file_hash"] == "b" * 64
    result = materialize_reviewed_cache(
        queue=queue,
        baseline_checkpoint=axis_rows(["Alpha beta"]),
        provider_checkpoint=providers,
        candidates=candidates,
        review_payload=decisions,
    )
    assert result["records"][["species", "trait", "value"]].to_dict("records") == [
        {
            "species": "Alpha beta",
            "trait": "mating_system",
            "value": "predominantly_outcrossing",
        }
    ]
