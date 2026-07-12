from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from island_v2.search_enabled_llm_campaign import _stable_shard, app


def test_stable_shard_is_deterministic() -> None:
    value = _stable_shard("Campanula punctata", 128)
    assert value == _stable_shard("Campanula punctata", 128)
    assert 0 <= value < 128


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
    jobs = [json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()]
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
        payloads.append(json.dumps({"job_id": job["job_id"], "result": row}))
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
    output = pd.read_csv(campaign / "trait_results.csv")
    assert len(output) == 2
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
    jobs2 = [json.loads(line) for line in (campaign / "jobs_shard_0000.jsonl").read_text().splitlines()]
    assert len(jobs2) == 1
    assert jobs2[0]["species"] == "Plantus rubra"
