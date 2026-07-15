from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

import island_v2.search_enabled_llm_campaign as campaign_module
from island_v2.search_enabled_llm_campaign import (
    OUTPUT_COLUMNS,
    _stable_shard,
    app,
    canonical_result_sha,
)


def test_stable_shard_is_deterministic() -> None:
    value = _stable_shard("Campanula punctata", 128)
    assert value == _stable_shard("Campanula punctata", 128)
    assert 0 <= value < 128


def test_canonical_result_sha_uses_normalized_output_column_json() -> None:
    row = {
        "species": "  Plantus   alba ",
        "flower_color": "white",
        "flower_shape": "open radial",
        "pollination_guild": "bees",
        "pollination_notes": " bee   visits ",
        "mating_system": "mixed_mating",
        "self_incompatibility": "likely_SC",
        "evidence_type": "inference",
        "confidence": "medium",
        "ignored": "not hashed",
    }
    normalized = {column: " ".join(str(row[column]).strip().split()) for column in OUTPUT_COLUMNS}
    material = json.dumps(
        normalized,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=False,
    )
    assert canonical_result_sha(row) == hashlib.sha256(material.encode()).hexdigest()


def test_prepare_and_ingest_roundtrip(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        [
            {"accepted_species": "Plantus alba", "genus": "Plantus", "family": "Plantaceae"},
            {"accepted_species": "Plantus rubra", "genus": "Plantus", "family": "Plantaceae"},
            {"accepted_species": "Floribunda testii", "genus": "Floribunda", "family": "Floraceae"},
        ]
    ).to_csv(master, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()

    result = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--shard-count",
            "1",
            "--shard-index",
            "0",
            "--batch-size",
            "2",
        ],
    )
    assert result.exit_code == 0, result.output
    jobs = [
        json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()
    ]
    assert len(jobs) == 2
    assert len({job["job_id"] for job in jobs}) == 2

    payloads = []
    for job in jobs:
        row = ",".join(
            [
                job["species"],
                "white",
                "tubular",
                "bees",
                "bee pollination likely",
                "mixed_mating",
                "likely_SC",
                "inference",
                "medium",
            ]
        )
        payloads.append(
            json.dumps(
                {
                    "job_id": job["job_id"],
                    "prompt_sha256": job["prompt_sha256"],
                    "result": row,
                }
            )
        )
    result_file = tmp_path / "results.jsonl"
    result_file.write_text("\n".join(payloads) + "\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign),
            "--results-jsonl",
            str(result_file),
            "--model-provider",
            "test",
            "--model-id",
            "test-model",
        ],
    )
    assert result.exit_code == 0, result.output
    output = pd.read_csv(campaign / "validated_trait_results.csv")
    assert len(output) == 2
    assert set(output["model_provider"]) == {"test"}
    assert set(output["model_id"]) == {"test-model"}
    assert not (campaign / "trait_results.csv").exists()
    ledger = pd.read_csv(campaign / "campaign_ledger.csv")
    assert int((ledger["status"] == "completed").sum()) == 2

    duplicate = runner.invoke(
        app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign),
            "--results-jsonl",
            str(result_file),
            "--model-provider",
            "test",
            "--model-id",
            "test-model",
        ],
    )
    assert duplicate.exit_code == 0, duplicate.output
    duplicate_report = json.loads(duplicate.output)
    assert duplicate_report["n_accepted"] == 0
    assert duplicate_report["n_rejected"] == 2
    assert duplicate_report["n_validated_total"] == 2
    assert len(pd.read_csv(campaign / "validated_trait_results.csv")) == 2
    ledger = pd.read_csv(campaign / "campaign_ledger.csv")
    assert int((ledger["status"] == "completed").sum()) == 2

    result = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--shard-count",
            "1",
            "--shard-index",
            "0",
            "--batch-size",
            "2",
        ],
    )
    assert result.exit_code == 0, result.output
    jobs2 = [
        json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()
    ]
    assert len(jobs2) == 1
    assert jobs2[0]["species"] == "Plantus rubra"


def test_prepare_can_resume_queued_jobs_without_completed_species(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        [
            {"accepted_species": "Plantus alba"},
            {"accepted_species": "Plantus rubra"},
            {"accepted_species": "Plantus viridis"},
        ]
    ).to_csv(master, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()
    common = [
        "prepare",
        "--master-csv",
        str(master),
        "--campaign-dir",
        str(campaign),
        "--shard-count",
        "1",
        "--shard-index",
        "0",
        "--batch-size",
        "2",
    ]

    first = runner.invoke(app, common)
    assert first.exit_code == 0, first.output
    first_jobs = [
        json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()
    ]

    resumed = runner.invoke(app, [*common, "--resume-queued"])
    assert resumed.exit_code == 0, resumed.output
    resumed_jobs = [
        json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()
    ]
    assert [job["job_id"] for job in resumed_jobs] == [job["job_id"] for job in first_jobs]
    ledger = pd.read_csv(campaign / "campaign_ledger.csv")
    resumed_species = {job["species"] for job in resumed_jobs}
    assert set(ledger.loc[ledger["attempts"] == 2, "species"]) == resumed_species

    result_rows = []
    for job in resumed_jobs:
        csv_row = ",".join(
            [
                job["species"],
                "unknown",
                "unknown",
                "unknown",
                "unknown",
                "unknown",
                "unknown",
                "inference",
                "low",
            ]
        )
        result_rows.append(
            json.dumps(
                {
                    "job_id": job["job_id"],
                    "prompt_sha256": job["prompt_sha256"],
                    "result": csv_row,
                }
            )
        )
    result_file = tmp_path / "resumed_results.jsonl"
    result_file.write_text("\n".join(result_rows) + "\n", encoding="utf-8")
    ingested = runner.invoke(
        app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign),
            "--results-jsonl",
            str(result_file),
            "--model-provider",
            "test",
            "--model-id",
            "test-model",
        ],
    )
    assert ingested.exit_code == 0, ingested.output

    after_completion = runner.invoke(app, [*common, "--resume-queued"])
    assert after_completion.exit_code == 0, after_completion.output
    remaining_jobs = [
        json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()
    ]
    assert len(remaining_jobs) == 1
    assert remaining_jobs[0]["species"] == "Plantus viridis"
    assert remaining_jobs[0]["species"] not in resumed_species


def test_ingest_rejects_missing_or_stale_prompt_hash(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame([{"accepted_species": "Plantus alba"}]).to_csv(master, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()

    prepared = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--shard-count",
            "1",
        ],
    )
    assert prepared.exit_code == 0, prepared.output
    job = json.loads((campaign / "jobs_shard_0000.jsonl").read_text().strip())
    csv_row = ",".join(
        [
            job["species"],
            "white",
            "tubular",
            "bees",
            "bee pollination likely",
            "mixed_mating",
            "likely_SC",
            "inference",
            "medium",
        ]
    )

    for name, payload, expected_error in (
        (
            "missing",
            {"job_id": job["job_id"], "result": csv_row},
            "missing prompt_sha256",
        ),
        (
            "stale",
            {
                "job_id": job["job_id"],
                "prompt_sha256": "0" * 64,
                "result": csv_row,
            },
            "prompt_sha256 mismatch",
        ),
    ):
        if name == "stale":
            retried = runner.invoke(
                app,
                [
                    "prepare",
                    "--master-csv",
                    str(master),
                    "--campaign-dir",
                    str(campaign),
                    "--shard-count",
                    "1",
                    "--retry-invalid",
                ],
            )
            assert retried.exit_code == 0, retried.output
        result_file = tmp_path / f"{name}.jsonl"
        result_file.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        ingested = runner.invoke(
            app,
            [
                "ingest",
                "--campaign-dir",
                str(campaign),
                "--results-jsonl",
                str(result_file),
                "--model-provider",
                "test",
                "--model-id",
                "test-model",
            ],
        )
        assert ingested.exit_code == 0, ingested.output
        report = json.loads(ingested.output)
        assert report["n_accepted"] == 0
        assert report["n_rejected"] == 1
        rejected = pd.read_csv(campaign / "rejected_results.csv")
        assert expected_error in rejected.iloc[-1]["error"]
        validated = pd.read_csv(campaign / "validated_trait_results.csv")
        assert validated.empty

    ledger = pd.read_csv(campaign / "campaign_ledger.csv")
    assert ledger.loc[0, "status"] == "invalid"


def test_prepare_scopes_to_applicable_species_and_only_unknown_fields(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        [
            {"accepted_species": "Alpha one", "family": "Aaceae"},
            {"accepted_species": "Beta two", "family": "Baceae"},
            {"accepted_species": "Gamma three", "family": "Gaceae"},
        ]
    ).to_csv(master, index=False)
    applicability = tmp_path / "applicability.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two", "Gamma three"],
            "angiosperm_analysis_eligible": [True, True, False],
        }
    ).to_csv(applicability, index=False)
    current = tmp_path / "current.csv"
    pd.DataFrame(
        [
            {
                "species": "Alpha one",
                "flower_color": "white",
                "flower_shape": "unknown",
                "pollination_guild": "unknown",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
            },
            {
                "species": "Beta two",
                "flower_color": "white",
                "flower_shape": "tubular",
                "pollination_guild": "bees",
                "mating_system": "mixed_mating",
                "self_incompatibility": "SC",
            },
            {
                "species": "Gamma three",
                "flower_color": "unknown",
                "flower_shape": "unknown",
                "pollination_guild": "unknown",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
            },
        ]
    ).to_csv(current, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()
    prepared = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--applicability-csv",
            str(applicability),
            "--current-traits-csv",
            str(current),
            "--shard-count",
            "1",
        ],
    )
    assert prepared.exit_code == 0, prepared.output
    report = json.loads(prepared.output)
    assert report["n_global_species"] == 3
    assert report["n_trait_applicable_species"] == 2
    assert report["n_trait_not_applicable_species"] == 1
    assert report["n_species_all_fields_known"] == 1
    jobs = [
        json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()
    ]
    assert [job["species"] for job in jobs] == ["Alpha one"]
    assert jobs[0]["target_fields"] == (
        "flower_shape|pollination_guild|mating_system|self_incompatibility"
    )

    # A provider response cannot overwrite the already-known, non-target colour.
    bad_result = tmp_path / "bad.jsonl"
    bad_result.write_text(
        json.dumps(
            {
                "job_id": jobs[0]["job_id"],
                "prompt_sha256": jobs[0]["prompt_sha256"],
                "result": (
                    "Alpha one,red,tubular,bees,bee inference,mixed_mating,likely_SC,inference,low"
                ),
            }
        )
        + "\n",
        encoding="utf-8",
    )
    ingested = runner.invoke(
        app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign),
            "--results-jsonl",
            str(bad_result),
            "--model-provider",
            "test",
            "--model-id",
            "test-model",
        ],
    )
    assert ingested.exit_code == 0, ingested.output
    assert json.loads(ingested.output)["n_rejected"] == 1
    rejected = pd.read_csv(campaign / "rejected_results.csv")
    assert "non-target field flower_color must be unknown" in rejected.loc[0, "error"]


def test_ingest_current_valid_result_wins_over_stale_line_for_same_job(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame([{"accepted_species": "Plantus alba"}]).to_csv(master, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()
    prepared = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--shard-count",
            "1",
        ],
    )
    assert prepared.exit_code == 0, prepared.output
    job = json.loads((campaign / "jobs_shard_0000.jsonl").read_text().strip())
    result_row = (
        "Plantus alba,white,tubular,bees,bee inference,mixed_mating,likely_SC,inference,medium"
    )
    result_file = tmp_path / "mixed_hashes.jsonl"
    result_file.write_text(
        "\n".join(
            [
                json.dumps(
                    {
                        "job_id": job["job_id"],
                        "prompt_sha256": "0" * 64,
                        "result": result_row,
                    }
                ),
                json.dumps(
                    {
                        "job_id": job["job_id"],
                        "job_sha256": job["job_sha256"],
                        "prompt_sha256": job["prompt_sha256"],
                        "result": result_row,
                    }
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    ingested = runner.invoke(
        app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign),
            "--results-jsonl",
            str(result_file),
            "--model-provider",
            "test",
            "--model-id",
            "test-model",
        ],
    )
    assert ingested.exit_code == 0, ingested.output
    report = json.loads(ingested.output)
    assert report["n_accepted"] == 1
    assert report["n_rejected"] == 1
    rejected = pd.read_csv(campaign / "rejected_results.csv", dtype=str).fillna("")
    assert "prompt_sha256 mismatch" in rejected.iloc[-1]["error"]
    validated = pd.read_csv(campaign / "validated_trait_results.csv", dtype=str).fillna("")
    assert len(validated) == 1
    assert validated.iloc[0]["result_sha256"] == canonical_result_sha(validated.iloc[0])
    ledger = pd.read_csv(campaign / "campaign_ledger.csv", dtype=str).fillna("")
    history = pd.read_csv(campaign / "job_history.csv", dtype=str).fillna("")
    assert ledger.iloc[0]["status"] == "completed"
    assert history.iloc[0]["status"] == "completed"
    assert history.iloc[0]["result_sha256"] == validated.iloc[0]["result_sha256"]


def test_multiwave_jobs_preserve_each_validated_field_and_job_provenance(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame([{"accepted_species": "Plantus alba", "family": "Plantaceae"}]).to_csv(
        master, index=False
    )
    current = tmp_path / "current.csv"
    pd.DataFrame(
        [
            {
                "species": "Plantus alba",
                "flower_color": "unknown",
                "flower_shape": "open",
                "pollination_guild": "bees",
                "mating_system": "mixed_mating",
                "self_incompatibility": "SC",
            }
        ]
    ).to_csv(current, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()

    def prepare_wave() -> dict[str, str]:
        result = runner.invoke(
            app,
            [
                "prepare",
                "--master-csv",
                str(master),
                "--campaign-dir",
                str(campaign),
                "--current-traits-csv",
                str(current),
                "--shard-count",
                "1",
            ],
        )
        assert result.exit_code == 0, result.output
        return json.loads((campaign / "jobs_shard_0000.jsonl").read_text().strip())

    def ingest_wave(job: dict[str, str], csv_row: str, filename: str) -> None:
        result_path = tmp_path / filename
        result_path.write_text(
            json.dumps(
                {
                    "job_id": job["job_id"],
                    "prompt_sha256": job["prompt_sha256"],
                    "result": csv_row,
                }
            )
            + "\n",
            encoding="utf-8",
        )
        result = runner.invoke(
            app,
            [
                "ingest",
                "--campaign-dir",
                str(campaign),
                "--results-jsonl",
                str(result_path),
                "--model-provider",
                "test",
                "--model-id",
                "test-model",
            ],
        )
        assert result.exit_code == 0, result.output
        assert json.loads(result.output)["n_accepted"] == 1

    color_job = prepare_wave()
    assert color_job["target_fields"] == "flower_color"
    ingest_wave(
        color_job,
        "Plantus alba,red,unknown,unknown,color evidence,unknown,unknown,flora,high",
        "color.jsonl",
    )

    pd.DataFrame(
        [
            {
                "species": "Plantus alba",
                "flower_color": "red",
                "flower_shape": "unknown",
                "pollination_guild": "bees",
                "mating_system": "mixed_mating",
                "self_incompatibility": "SC",
            }
        ]
    ).to_csv(current, index=False)
    shape_job = prepare_wave()
    assert shape_job["target_fields"] == "flower_shape"
    assert shape_job["job_id"] != color_job["job_id"]
    ingest_wave(
        shape_job,
        "Plantus alba,unknown,tubular,unknown,shape evidence,unknown,unknown,flora,high",
        "shape.jsonl",
    )

    validated = pd.read_csv(campaign / "validated_trait_results.csv", dtype=str).fillna("")
    assert len(validated) == 2
    assert set(validated["job_id"]) == {color_job["job_id"], shape_job["job_id"]}
    assert set(validated.loc[validated["flower_color"].ne("unknown"), "flower_color"]) == {"red"}
    assert set(validated.loc[validated["flower_shape"].ne("unknown"), "flower_shape"]) == {
        "tubular"
    }
    assert all(
        row["result_sha256"] == canonical_result_sha(row) for row in validated.to_dict("records")
    )
    history = pd.read_csv(campaign / "job_history.csv", dtype=str).fillna("")
    assert len(history) == 2
    assert set(history["status"]) == {"completed"}
    assert set(history["job_sha256"]) == set(validated["job_sha256"])


def test_prepare_purges_results_when_actual_prompt_policy_becomes_stale(
    tmp_path: Path,
    monkeypatch,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame([{"accepted_species": "Plantus alba"}]).to_csv(master, index=False)
    campaign = tmp_path / "campaign"
    runner = CliRunner()
    first = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--shard-count",
            "1",
        ],
    )
    assert first.exit_code == 0, first.output
    first_job = json.loads((campaign / "jobs_shard_0000.jsonl").read_text().strip())
    result_file = tmp_path / "result.jsonl"
    result_file.write_text(
        json.dumps(
            {
                "job_id": first_job["job_id"],
                "prompt_sha256": first_job["prompt_sha256"],
                "result": (
                    "Plantus alba,white,tubular,bees,bee evidence,mixed_mating,"
                    "likely_SC,inference,medium"
                ),
            }
        )
        + "\n",
        encoding="utf-8",
    )
    ingested = runner.invoke(
        app,
        [
            "ingest",
            "--campaign-dir",
            str(campaign),
            "--results-jsonl",
            str(result_file),
            "--model-provider",
            "test",
            "--model-id",
            "test-model",
        ],
    )
    assert ingested.exit_code == 0, ingested.output

    monkeypatch.setattr(
        campaign_module,
        "CONFIDENCE",
        {*campaign_module.CONFIDENCE, "very_high"},
    )
    refreshed = runner.invoke(
        app,
        [
            "prepare",
            "--master-csv",
            str(master),
            "--campaign-dir",
            str(campaign),
            "--shard-count",
            "1",
        ],
    )
    assert refreshed.exit_code == 0, refreshed.output
    report = json.loads(refreshed.output)
    assert report["n_stale_validated_results_removed"] == 1
    assert pd.read_csv(campaign / "validated_trait_results.csv").empty
    new_job = json.loads((campaign / "jobs_shard_0000.jsonl").read_text().strip())
    assert new_job["job_id"] != first_job["job_id"]
    assert new_job["prompt_sha256"] == first_job["prompt_sha256"]
    assert new_job["prompt_policy_sha256"] != first_job["prompt_policy_sha256"]
    history = pd.read_csv(campaign / "job_history.csv", dtype=str).fillna("")
    assert set(history["status"]) == {"stale", "queued"}
    input_manifest = json.loads(
        (campaign / "campaign_input_manifest.json").read_text(encoding="utf-8")
    )
    assert input_manifest["master_csv_sha256"] == hashlib.sha256(master.read_bytes()).hexdigest()
    assert input_manifest["prompt_policy_sha256"] == report["prompt_policy_sha256"]
