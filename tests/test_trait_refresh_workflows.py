from pathlib import Path


WORKFLOW_DIR = Path(".github/workflows")


def _workflow(name: str) -> str:
    return (WORKFLOW_DIR / name).read_text(encoding="utf-8")


def test_refresh_chain_has_one_workflow_run_path_per_source() -> None:
    all_master = _workflow("build-all-master-trait-ledger.yml")
    direct_sources = {
        "ingest-eol-traitbank-bulk.yml": "ingest_eol_traitbank_bulk",
        "ingest-austraits-all-floral-bulk.yml": (
            "Ingest all AusTraits target floral traits"
        ),
        "ingest-gift-direct-traits.yml": "Acquire GIFT direct floral traits",
        "ingest-usda-plants-direct-traits.yml": "Acquire USDA PLANTS direct floral traits",
    }
    for workflow_name, workflow_title in direct_sources.items():
        source = _workflow(workflow_name)
        assert workflow_title in all_master
        assert "gh workflow run build-all-master-trait-ledger.yml" not in source

    efloras = _workflow("fetch-efloras-all.yml")
    efloras_promotion = _workflow("promote-efloras-accounted-run.yml")
    assert 'workflows: ["Fetch eFloras species treatments"]' in efloras_promotion
    assert "gh workflow run promote-efloras-accounted-run.yml" not in efloras
    assert "gh workflow run build-all-master-trait-ledger.yml" not in efloras_promotion

    public_web = _workflow("run-public-web-trait-shards.yml")
    public_web_promotion = _workflow("promote-public-web-accounted-run.yml")
    assert 'workflows: ["run_public_web_trait_shards"]' in public_web_promotion
    assert "gh workflow run promote-public-web-accounted-run.yml" not in public_web
    assert "gh workflow run build-all-master-trait-ledger.yml" not in public_web_promotion
    assert "github.event_name == 'workflow_run'" in public_web_promotion


def test_all_master_pin_is_exact_and_still_requires_production_contract() -> None:
    workflow = _workflow("build-all-master-trait-ledger.yml")
    assert "pinned_source_variable:" in workflow
    assert 'if [[ "$pinned_source" == "$variable" ]]' in workflow
    assert "exact verified production source" in workflow
    assert "lacks production marker" in workflow
    assert 'test "$(jq -r .path' in workflow
    assert "eol-traitbank-production-" in workflow
    assert (
        "resolve_run EOL_RUN_ID ingest-eol-traitbank-bulk.yml" in workflow
    )


def test_public_web_chain_preserves_300_per_shard_contract_and_checkpoints() -> None:
    promotion = _workflow("promote-public-web-accounted-run.yml")
    assert "Continue the full campaign from this immutable source run" in promotion
    assert "public_web_continuation_v1" in promotion
    assert "legacy source run has no continuation contract" in promotion.lower()
    assert 'echo "batch_size=300"' in promotion
    assert '-f batch_size="$BATCH_SIZE"' in promotion
    assert "batch_size=500" not in promotion

    acquisition = _workflow("run-public-web-trait-shards.yml")
    assert 'default: "300"' in acquisition
    assert '"batch_size_scope": "per_shard"' in acquisition
    assert '"maximum_species_advanced_across_selected_shards"' in acquisition
    assert "SEED_CACHE_PREFIX" in acquisition
    assert "needs.plan.outputs.seed_run_id" in acquisition
    assert "batch_size=500" not in acquisition
    assert "Fail when the shard checkpoint contract is missing" in acquisition
    assert "steps.status.outputs.ready != 'true'" in acquisition
    assert '"campaign_status.json",' in acquisition


def test_grounded_validation_requires_exact_v4_fifteen_trait_artifact() -> None:
    materialize = _workflow("materialize-trait-fill-cascade.yml")
    assert "trait_fill_cascade_artifact_checkpoint_v2" in materialize
    assert "trait_fill_cascade_shards_v4" in materialize
    assert '"n_traits": len(target_traits)' in materialize
    assert '"global_fallback_rows"' in materialize

    validation = _workflow("run-trait-inference-validation.yml")
    assert "trait_fill_cascade_artifact_lock_v2" in validation
    assert "trait_fill_cascade_artifact_checkpoint_v2" in validation
    assert "trait_fill_cascade_shards_v4" in validation
    assert "set(manifest['trait_direct_counts']) != expected_traits" in validation
    assert "prediction table does not cover all 15 configured traits" in validation
    assert "grid['n'] + grid['n_unresolved_predictions']" in validation
