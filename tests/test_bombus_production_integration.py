from __future__ import annotations

from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from island_v2.bombus_absence_evidence import app as absence_app
from island_v2.bombus_channel_score import app as channel_app
from island_v2.bombus_niche_hypervolume import app as niche_app


def test_raw_environment_and_occurrence_inputs_produce_channel_scores(
    tmp_path: Path,
) -> None:
    runner = CliRunner()
    occurrences_path = tmp_path / "island_pollinator_occurrences.csv"
    pd.DataFrame(
        [
            {
                "island_id": "island_1",
                "family": "Apidae",
                "genus": "Bombus",
                "dataset_key": "dataset_1",
                "year": 2025,
                "decimal_latitude": 34.0,
                "decimal_longitude": 139.0,
                "gbif_id": "1",
            }
        ]
    ).to_csv(occurrences_path, index=False)

    final_dir = tmp_path / "final"
    absence = runner.invoke(
        absence_app,
        [
            "run",
            "--occurrences-csv",
            str(occurrences_path),
            "--output-dir",
            str(final_dir),
        ],
    )
    assert absence.exit_code == 0, absence.output

    training_path = tmp_path / "bombus_occurrence_environment.csv"
    training_rows = []
    for index in range(20):
        training_rows.append(
            {
                "bombus_species": "Bombus alpha",
                "bio1": 10.0 + index * 0.1,
                "bio12": 900.0 + index * 5.0 + (index % 3),
            }
        )
    pd.DataFrame(training_rows).to_csv(training_path, index=False)

    targets_path = tmp_path / "island_source_pool_environment.csv"
    pd.DataFrame(
        [
            {
                "island_id": "island_1",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 951.0,
            }
        ]
    ).to_csv(targets_path, index=False)

    niche_dir = tmp_path / "niche"
    niche = runner.invoke(
        niche_app,
        [
            "run",
            "--occurrence-environment-csv",
            str(training_path),
            "--island-source-pool-environment-csv",
            str(targets_path),
            "--environment-columns",
            "bio1,bio12",
            "--output-dir",
            str(niche_dir),
        ],
    )
    assert niche.exit_code == 0, niche.output

    channel_dir = tmp_path / "channel"
    channel = runner.invoke(
        channel_app,
        [
            "run",
            "--diagnostics-csv",
            str(final_dir / "island_bombus_occurrence_evidence.csv"),
            "--environmental-scores-csv",
            str(niche_dir / "bombus_species_environmental_compatibility.csv"),
            "--output-dir",
            str(channel_dir),
        ],
    )
    assert channel.exit_code == 0, channel.output

    result = pd.read_csv(channel_dir / "bombus_channel_continuous_scores.csv").iloc[0]
    assert result["observation_availability"] > 0.5
    assert pd.notna(result["environmental_compatibility"])
    assert pd.notna(result["bombus_channel_availability"])
