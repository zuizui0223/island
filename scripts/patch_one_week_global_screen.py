from pathlib import Path

campaign_path = Path("src/island_v2/global_trait_campaign.py")
text = campaign_path.read_text(encoding="utf-8")
old = '''def choose_active_task(ledger: pd.DataFrame, config: dict[str, Any]) -> str:
    """Choose the first task with at least one currently eligible species."""
    for task_value in config["task_order"]:
        task = str(task_value)
        if task_eligible_mask(ledger, task, config).any():
            return task
    return "complete"
'''
new = '''def choose_active_task(ledger: pd.DataFrame, config: dict[str, Any]) -> str:
    """Choose the first required primary-screen task with eligible species."""
    for task_value in config["task_order"]:
        task = str(task_value)
        if not bool(config["tasks"][task].get("required_for_primary_completion", True)):
            continue
        if task_eligible_mask(ledger, task, config).any():
            return task
    return "complete"
'''
if old not in text:
    raise SystemExit("choose_active_task marker not found")
text = text.replace(old, new, 1)
old = '''            "n_eligible_pending": int(task_eligible_mask(ledger, task, config).sum()),
            "terminal": task_is_terminal(ledger, task),
'''
new = '''            "n_eligible_pending": int(task_eligible_mask(ledger, task, config).sum()),
            "required_for_primary_completion": bool(
                config["tasks"][task].get("required_for_primary_completion", True)
            ),
            "terminal": task_is_terminal(ledger, task),
'''
if old not in text:
    raise SystemExit("campaign summary marker not found")
text = text.replace(old, new, 1)
campaign_path.write_text(text, encoding="utf-8")

config_path = Path("config/global_trait_campaign.yml")
config = config_path.read_text(encoding="utf-8")
config = config.replace(
    "  - floral_access_wikimedia\n  - alternative_guild_openalex\n",
    "  - floral_access_wikimedia\n  - alternative_guild_wikimedia\n  - alternative_guild_openalex\n",
    1,
)
config = config.replace(
    "  reproductive_openalex:\n    phase: reproductive_and_pollen_vector\n",
    "  reproductive_openalex:\n    phase: reproductive_and_pollen_vector\n    required_for_primary_completion: false\n",
    1,
)
config = config.replace(
    "    depends_on:\n      - reproductive_wikimedia\n      - reproductive_openalex\n",
    "    depends_on:\n      - reproductive_wikimedia\n",
    1,
)
marker = '''  alternative_guild_openalex:
    phase: alternative_pollinator_function
'''
replacement = '''  alternative_guild_wikimedia:
    phase: alternative_pollinator_function
    source_kind: wikimedia_reported_ecology
    eligibility: machine_biotic_candidate
    depends_on:
      - floral_access_wikimedia
    pause_seconds: 0.05
    max_retries: 4
    target_traits:
      - pollination_functional_guild

  alternative_guild_openalex:
    phase: alternative_pollinator_function
    required_for_primary_completion: false
'''
if marker not in config:
    raise SystemExit("alternative guild config marker not found")
config = config.replace(marker, replacement, 1)
config = config.replace(
    "    depends_on:\n      - floral_access_wikimedia\n    pause_seconds: 0.25\n",
    "    depends_on:\n      - alternative_guild_wikimedia\n    pause_seconds: 0.25\n",
    1,
)
config_path.write_text(config, encoding="utf-8")

workflow_path = Path(".github/workflows/drive-global-trait-campaign.yml")
workflow = workflow_path.read_text(encoding="utf-8")
workflow = workflow.replace(
    '''      # Two scheduled runs per hour yield at most 720 OpenAlex searches/day:
      # (10 + 5) * 48. This leaves headroom below the free-key 1,000-search budget.
      OPENALEX_BATCH_SIZE: "10"
      FLORAL_BATCH_SIZE: "150"
      GUILD_BATCH_SIZE: "5"
      OPENALEX_API_KEY: ${{ secrets.OPENALEX_API_KEY }}
''',
    '''      # 1,000 species every 30 minutes gives a theoretical 48,000-species/day
      # first-pass capacity; the 105,612-species master therefore fits within a week
      # even with substantial scheduling and network slack.
      FLORAL_BATCH_SIZE: "150"
      GUILD_BATCH_SIZE: "150"
''',
    1,
)
credential_step = '''      - name: Report OpenAlex credential state
        run: |
          if test -n "${OPENALEX_API_KEY}"; then
            echo "OpenAlex free-key lane enabled." >> "$GITHUB_STEP_SUMMARY"
          else
            echo "OpenAlex tasks will remain retryable without consuming attempts because OPENALEX_API_KEY is not configured." >> "$GITHUB_STEP_SUMMARY"
          fi
'''
workflow = workflow.replace(credential_step, "", 1)
workflow = workflow.replace(
    '''          run_task reproductive_wikimedia "${BATCH_SIZE}"
          run_task reproductive_openalex "${OPENALEX_BATCH_SIZE}"
          run_task floral_access_wikimedia "${FLORAL_BATCH_SIZE}"
          run_task alternative_guild_openalex "${GUILD_BATCH_SIZE}"
''',
    '''          run_task reproductive_wikimedia "${BATCH_SIZE}"
          run_task floral_access_wikimedia "${FLORAL_BATCH_SIZE}"
          run_task alternative_guild_wikimedia "${GUILD_BATCH_SIZE}"
''',
    1,
)
workflow = workflow.replace(
    '''                  "reproductive_wikimedia": int(os.environ["BATCH_SIZE"]),
                  "reproductive_openalex": int(os.environ["OPENALEX_BATCH_SIZE"]),
                  "floral_access_wikimedia": int(os.environ["FLORAL_BATCH_SIZE"]),
                  "alternative_guild_openalex": int(os.environ["GUILD_BATCH_SIZE"]),
''',
    '''                  "reproductive_wikimedia": int(os.environ["BATCH_SIZE"]),
                  "floral_access_wikimedia": int(os.environ["FLORAL_BATCH_SIZE"]),
                  "alternative_guild_wikimedia": int(os.environ["GUILD_BATCH_SIZE"]),
''',
    1,
)
workflow = workflow.replace(
    '''                  "OpenAlex credential, quota, 429, 5xx, and transport failures do not consume "
                  "the biological lookup-attempt budget or imply trait absence."
''',
    '''                  "Wikimedia lookup failures do not imply trait absence; optional scholarly "
                  "enrichment is outside the primary one-week completion path."
''',
    1,
)
workflow_path.write_text(workflow, encoding="utf-8")

index_path = Path("src/island_v2/global_candidate_index.py")
index_text = index_path.read_text(encoding="utf-8")
index_text = index_text.replace(
    '''        "floral_access_wikimedia_status",
        "alternative_guild_openalex_status",
''',
    '''        "floral_access_wikimedia_status",
        "alternative_guild_wikimedia_status",
''',
    1,
)
index_text = index_text.replace(
    '''    base = base[
        [
            "accepted_species",
            "family",
            "n_islands",
            "n_records",
            "reproductive_wikimedia_status",
            "reproductive_openalex_status",
            "floral_access_wikimedia_status",
            "alternative_guild_openalex_status",
        ]
    ].drop_duplicates("accepted_species")
''',
    '''    status_columns = [
        "accepted_species",
        "family",
        "n_islands",
        "n_records",
        "reproductive_wikimedia_status",
        "floral_access_wikimedia_status",
        "alternative_guild_wikimedia_status",
    ]
    for optional in ("reproductive_openalex_status", "alternative_guild_openalex_status"):
        if optional in base.columns:
            status_columns.append(optional)
    base = base[status_columns].drop_duplicates("accepted_species")
''',
    1,
)
index_path.write_text(index_text, encoding="utf-8")

test_path = Path("tests/test_global_trait_campaign.py")
tests = test_path.read_text(encoding="utf-8")
tests = tests.replace(
    '''            "floral_access_wikimedia",
            "alternative_guild_openalex",
''',
    '''            "floral_access_wikimedia",
            "alternative_guild_wikimedia",
            "alternative_guild_openalex",
''',
    1,
)
tests = tests.replace(
    '''                "source_kind": "openalex_reported_ecology",
                "depends_on": ["reproductive_wikimedia"],
''',
    '''                "source_kind": "openalex_reported_ecology",
                "required_for_primary_completion": False,
                "depends_on": ["reproductive_wikimedia"],
''',
    1,
)
tests = tests.replace(
    '''                "depends_on": ["reproductive_wikimedia", "reproductive_openalex"],
''',
    '''                "depends_on": ["reproductive_wikimedia"],
''',
    1,
)
tests = tests.replace(
    '''            "alternative_guild_openalex": {
                "phase": "alternative_guild",
                "source_kind": "openalex_reported_ecology",
                "depends_on": ["floral_access_wikimedia"],
''',
    '''            "alternative_guild_wikimedia": {
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
''',
    1,
)
tests = tests.replace(
    '''    for task in ("reproductive_wikimedia", "reproductive_openalex"):
        ledger[f"{task}_status"] = "processed"
''',
    '''    ledger["reproductive_wikimedia_status"] = "processed"
''',
    1,
)
tests = tests.replace(
    '''    ledger.loc[first, "reproductive_openalex_status"] = "processed"
''',
    "",
    1,
)
tests = tests.replace(
    '''    assert prepared.loc[second, "floral_access_wikimedia_status"].iloc[0] == "pending"
    assert prepared.loc[eligible, "accepted_species"].tolist() == ["Alpha one"]
    assert batch["accepted_species"].tolist() == ["Alpha one"]
''',
    '''    assert prepared.loc[second, "floral_access_wikimedia_status"].iloc[0] == "pending"
    assert prepared.loc[eligible, "accepted_species"].tolist() == ["Alpha one", "Alpha two"]
    assert batch["accepted_species"].tolist() == ["Alpha one", "Alpha two"]
''',
    1,
)
tests = tests.replace(
    '''    for task in ("reproductive_wikimedia", "reproductive_openalex"):
        ledger.loc[first, f"{task}_status"] = "processed"
''',
    '''    ledger.loc[first, "reproductive_wikimedia_status"] = "processed"
''',
    1,
)
tests += '''


def test_optional_openalex_does_not_block_primary_completion():
    config = _config()
    ledger = campaign.reconcile_ledger(_master(), None, config)
    ledger["reproductive_wikimedia_status"] = "processed"
    ledger["floral_access_wikimedia_status"] = "not_applicable"
    ledger["alternative_guild_wikimedia_status"] = "not_applicable"

    assert set(ledger["reproductive_openalex_status"]) == {"pending"}
    assert campaign.choose_active_task(ledger, config) == "complete"
'''
test_path.write_text(tests, encoding="utf-8")

index_test_path = Path("tests/test_global_candidate_index.py")
index_tests = index_test_path.read_text(encoding="utf-8")
index_tests = index_tests.replace(
    '''                "alternative_guild_openalex_status": "processed",
''',
    '''                "alternative_guild_wikimedia_status": "processed",
                "alternative_guild_openalex_status": "pending",
''',
    1,
)
index_tests = index_tests.replace(
    '''                "alternative_guild_openalex_status": "pending",
''',
    '''                "alternative_guild_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
''',
    1,
)
index_tests = index_tests.replace(
    '''                "alternative_guild_openalex_status": "pending",
''',
    '''                "alternative_guild_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
''',
    1,
)
index_tests = index_tests.replace(
    '''    assert beta["alternative_guild_openalex_status"] == "pending"
''',
    '''    assert beta["alternative_guild_wikimedia_status"] == "pending"
''',
    1,
)
index_test_path.write_text(index_tests, encoding="utf-8")

docs_path = Path("docs/global_trait_campaign.md")
docs = docs_path.read_text(encoding="utf-8")
docs += '''

## One-week global primary screen

The scheduled primary path is intentionally database-light: 1,000 globally
family-balanced species are screened through Wikimedia every 30 minutes, for a
theoretical capacity of 48,000 species per day. Floral-access and alternative
guild Wikimedia extraction run only for direct machine biotic-vector candidates.
OpenAlex remains an optional scholarly-enrichment lane and is not required for
primary completion. The one-week target refers to complete machine first-pass
coverage, not manual acceptance of every candidate trait.
'''
docs_path.write_text(docs, encoding="utf-8")
