from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_p1c_matched_genus_aggregate_v2 import (
    P1CAggregateError,
    aggregate_existing_shards,
)


def _config(path: Path, *, shards: int = 2, per_shard: int = 2) -> None:
    cfg = {
        "contract": "chapter1_p1_assembly_depth_defense_v1",
        "pinned_primary_artifact": {"artifact_digest": "sha256:abc"},
        "p1c_matched_complexity_genus_null": {
            "permutations": shards * per_shard,
            "execution_shards": shards,
            "permutations_per_shard": per_shard,
            "minimum_valid_permutations": shards * per_shard,
            "random_seed": 100,
            "permutation_seed_rule": "master_seed_plus_zero_based_permutation_id",
            "primary_pass_rule": {"max_randomization_p_value": 0.05},
        },
    }
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")


def _shard(root: Path, shard: int, observed: float, values: list[float]) -> None:
    start = shard * len(values)
    rows = []
    for offset, value in enumerate(values):
        pid = start + offset
        rows.append(
            {
                "permutation_id": pid,
                "seed": 100 + pid,
                "valid": True,
                "primary_median_conditional_attenuation": value,
            }
        )
    pd.DataFrame(rows).to_csv(root / f"p1c_null_shard_{shard:02d}.csv.gz", index=False)
    manifest = {
        "contract": "chapter1_p1c_matched_genus_null_shard_v1",
        "shard_id": shard,
        "n_rows": len(values),
        "n_valid": len(values),
        "observed_primary_statistic": observed,
        "source_artifact_digest": "sha256:abc",
    }
    (root / f"p1c_null_shard_{shard:02d}_manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )


def test_aggregate_accepts_numerical_equivalence_without_regenerating(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg.yml"
    _config(cfg)
    _shard(tmp_path, 0, 0.7200000000000000, [0.1, 0.2])
    _shard(tmp_path, 1, 0.7200000000000001, [0.3, 0.4])
    out = tmp_path / "out"
    result = aggregate_existing_shards(
        shard_dir=tmp_path,
        p1_config_path=cfg,
        output_dir=out,
        source_run_id=123,
    )
    assert result["permutations_regenerated"] is False
    assert result["valid_permutations"] == 4
    assert result["observed_statistic_spread_across_shards"] < 1e-12
    assert result["source_permutation_run_id"] == 123


def test_aggregate_rejects_material_observed_disagreement(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg.yml"
    _config(cfg)
    _shard(tmp_path, 0, 0.72, [0.1, 0.2])
    _shard(tmp_path, 1, 0.720001, [0.3, 0.4])
    with pytest.raises(P1CAggregateError, match="exceeds tolerance"):
        aggregate_existing_shards(
            shard_dir=tmp_path,
            p1_config_path=cfg,
            output_dir=tmp_path / "out",
        )
