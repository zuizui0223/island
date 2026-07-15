from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.combine_public_web_artifacts import combine_public_web_run
from island_v2.v1_category_traits import OUTPUT_COLUMNS


def _write_shard(root: Path, index: int, species: str) -> None:
    shard = root / f"public-web-trait-shard-{index}-123"
    cumulative = shard / "cumulative"
    cumulative.mkdir(parents=True)
    (shard / "campaign_manifest.json").write_text(
        json.dumps(
            {
                "contract_version": "public_web_9col_shards_v4",
                "master_fingerprint": "fingerprint",
                "n_global_species": 2,
                "shard_index": index,
                "shard_count": 2,
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "species": species,
                "genus": species.split()[0],
                "family": "Exampleaceae",
                "shard_index": index,
                "status": "completed",
                "attempts": 1,
                "last_packet_id": "batch_000001",
                "last_error": "",
                "result_sha256": f"result-{index}",
                "updated_at": "2026-07-12T00:00:00Z",
            }
        ]
    ).to_csv(shard / "checkpoint.csv", index=False)
    pd.DataFrame(
        [
            {
                "species": species,
                "flower_color": "red" if index == 0 else "unknown",
                "flower_shape": "tubular" if index == 0 else "unknown",
                "pollination_guild": "bees" if index == 0 else "unknown",
                "pollination_notes": "direct example" if index == 0 else "unknown",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
                "evidence_type": "flora" if index == 0 else "inference",
                "confidence": "medium" if index == 0 else "low",
            }
        ],
        columns=OUTPUT_COLUMNS,
    ).to_csv(cumulative / "trait_results.csv", index=False)
    pd.DataFrame(
        [
            {
                "species": species,
                "field": "flower_color",
                "value": "red",
                "matched_term": "red",
                "source_type": "flora_or_monograph",
                "source_url": f"https://flora.test/{index}",
                "source_citation": "Example flora",
                "source_excerpt": "Flowers red.",
                "evidence_type": "flora",
                "confidence": "medium",
                "source_tier": "A_flora",
                "rule_id": "color_red",
            }
        ]
    ).to_csv(cumulative / "trait_evidence.csv", index=False)


def test_combine_validates_population_and_preserves_exact_rows(tmp_path: Path) -> None:
    root = tmp_path / "artifacts"
    _write_shard(root, 0, "Plantus alpha")
    _write_shard(root, 1, "Plantus beta")
    metadata = tmp_path / "metadata.jsonl"
    metadata.write_text(
        "\n".join(
            json.dumps(
                {
                    "id": 100 + index,
                    "name": f"public-web-trait-shard-{index}-123",
                    "digest": f"sha256:digest{index}",
                    "expired": False,
                    "workflow_run": {"head_sha": "head"},
                }
            )
            for index in range(2)
        )
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "combined"
    report = combine_public_web_run(
        artifact_root=root,
        output_dir=output,
        source_run_id=123,
        expected_shards=2,
        expected_species=2,
        artifact_metadata_jsonl=metadata,
    )

    assert report["checkpoint_status_counts"] == {"completed": 2}
    assert report["n_result_species"] == 2
    assert report["n_species_any_requested_value"] == 1
    assert report["n_trait_evidence_rows"] == 2
    assert report["source_artifact_count"] == 2
    result = pd.read_csv(output / "combined_trait_results.csv", dtype=str)
    assert list(result.columns) == OUTPUT_COLUMNS
    assert set(result["species"]) == {"Plantus alpha", "Plantus beta"}
    checkpoint = pd.read_csv(output / "combined_checkpoint.csv.gz", dtype=str)
    assert checkpoint["species"].nunique() == 2
